"""Remove the hypothetical formylhopane reductase reaction (RXN-13531) and retrace every query.

The test removes the three MetaCyc hopanoid reactions examined in Table S8; only
RXN-13531 occurs in the final networks (as a broader-criteria reaction). Writes the
Table S1 native-control files.

    python analyses/03_upstream_tracing/formylhopane_reductase_test.py \\
        --networks results/networks --out results/tracing
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.controls import alias_roots, controlled_discovery  # noqa: E402
from biosynthesis.graph import (  # noqa: E402
    CRITERIA,
    CRITERIA_NAMES,
    VERSIONS,
    identity_queries,
    load_network,
    node_label,
    without_reactions,
)
from biosynthesis.panel import load_identities  # noqa: E402
from biosynthesis.tables import write_csv  # noqa: E402

REMOVED_REACTIONS = {
    "RXN-13529": "HpnG hydrolysis equation; later characterization supports phosphorolysis",
    "RXN-13530": "putative spontaneous interconversion of ribosylhopane and formylhopane",
    "RXN-13531": "described in MetaCyc as a hypothetical reductase reaction",
}
QUERY_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "biomarker_number",
    "name",
    "identity_id",
    "identity",
    "in_network",
    "query_node",
    "node",
    "producing_reactions_with_rxn13531",
    "producing_reactions_without_rxn13531",
    "traced_network_changed",
    "compounds_removed",
    "reactions_removed",
)
GROUPS = (
    ("accepted", "accepted"),
    ("candidate", "candidate_not_accepted"),
    ("no_candidate", "no_candidate"),
)
CONDITIONS = (("with_rxn13531", "before"), ("without_rxn13531", "after"))


def query_changes(network, targets: list[str], declarations: dict) -> dict[str, dict]:
    """Trace each query with and without the removed reactions and record what changes."""
    controlled, removed = without_reactions(network, REMOVED_REACTIONS)
    roots = alias_roots(network, declarations, set(targets))
    changes = {}
    for target in sorted(set(targets)):
        evaluated = roots.get(target, target)
        before = controlled_discovery(network, evaluated)
        affected = bool(set(removed) & set(before["edge_ids"]))
        after = controlled_discovery(controlled, evaluated) if affected else before
        changes[target] = {
            "query_node": node_label(network.nodes[target]),
            "node": node_label(network.nodes[evaluated]),
            "before": len(network.producers.get(evaluated, ())),
            "after": len(controlled.producers.get(evaluated, ())),
            "traced_network_changed": before != after,
            "compounds_removed": len(
                set(before["required_nodes"]) - set(after["required_nodes"])
            ),
            "stored_removed": len(
                set(before["stored_nodes"]) - set(after["stored_nodes"])
            ),
            "reactions_removed": len(set(before["edge_ids"]) - set(after["edge_ids"])),
        }
    return changes


def query_rows(queries: list[dict], changes: dict) -> list[dict]:
    """One row per criteria, network version, identity and query state."""
    rows = []
    for criteria in CRITERIA:
        for version in VERSIONS:
            for query in (row for row in queries if row["criteria"] == criteria):
                identity = query["identity"]
                for target in query["targets"] or [None]:
                    row = {
                        "network_version": version,
                        "reaction_criteria": CRITERIA_NAMES[criteria],
                        "biomarker_number": identity["biomarker_number"],
                        "name": identity["name"],
                        "identity_id": identity["identity_id"],
                        "identity": identity["identity"],
                        "in_network": target is not None,
                    }
                    if target is not None:
                        change = changes[version, criteria][target]
                        row.update(
                            query_node=change["query_node"],
                            node=change["node"],
                            producing_reactions_with_rxn13531=change["before"],
                            producing_reactions_without_rxn13531=change["after"],
                            traced_network_changed=change["traced_network_changed"],
                            compounds_removed=change["compounds_removed"],
                            reactions_removed=change["reactions_removed"],
                        )
                    rows.append(row)
    return rows


def biomarker_rows(rows: list[dict]) -> list[dict]:
    """One row per network version, criteria and biomarker: identities, traced nodes, producers."""
    order = {name: criteria for criteria, name in CRITERIA_NAMES.items()}
    groups = sorted(
        {
            (
                row["network_version"],
                order[row["reaction_criteria"]],
                row["biomarker_number"],
            )
            for row in rows
        },
        key=lambda key: (
            {"primary": "primary", "with_conditional": "sensitivity"}[key[0]],
            key[1],
            key[2],
        ),
    )
    output = []
    for version, criteria, number in groups:
        subset = [
            row
            for row in rows
            if (
                row["network_version"],
                order[row["reaction_criteria"]],
                row["biomarker_number"],
            )
            == (version, criteria, number)
        ]
        panel = {
            "network_version": version,
            "reaction_criteria": CRITERIA_NAMES[criteria],
            "biomarker_number": number,
            "name": subset[0]["name"],
        }
        for label, identity in GROUPS:
            selected = [row for row in subset if row["identity"] == identity]
            panel[f"{label}_identities"] = len({row["identity_id"] for row in selected})
            panel[f"{label}_identity_nodes_in_network"] = len(
                {row["node"] for row in selected if row["in_network"]}
            )
            if identity == "no_candidate":
                continue
            for suffix, _ in CONDITIONS:
                panel[f"{label}_identity_nodes_with_producer_{suffix}"] = len(
                    {
                        row["node"]
                        for row in selected
                        if row["in_network"]
                        and row[f"producing_reactions_{suffix}"] > 0
                    }
                )
        for suffix, _ in CONDITIONS:
            panel[f"accepted_identity_has_producer_{suffix}"] = (
                panel[f"accepted_identity_nodes_with_producer_{suffix}"] > 0
            )
            panel[f"candidate_identity_has_producer_{suffix}"] = (
                panel[f"candidate_identity_nodes_with_producer_{suffix}"] > 0
            )
            panel[f"has_producer_{suffix}"] = any(
                row["in_network"] and row[f"producing_reactions_{suffix}"] > 0
                for row in subset
            )
        output.append(panel)
    return output


BIOMARKER_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "biomarker_number",
    "name",
    "accepted_identities",
    "accepted_identity_nodes_in_network",
    "accepted_identity_nodes_with_producer_with_rxn13531",
    "accepted_identity_nodes_with_producer_without_rxn13531",
    "candidate_identities",
    "candidate_identity_nodes_in_network",
    "candidate_identity_nodes_with_producer_with_rxn13531",
    "candidate_identity_nodes_with_producer_without_rxn13531",
    "no_candidate_identities",
    "no_candidate_identity_nodes_in_network",
    "accepted_identity_has_producer_with_rxn13531",
    "candidate_identity_has_producer_with_rxn13531",
    "has_producer_with_rxn13531",
    "accepted_identity_has_producer_without_rxn13531",
    "candidate_identity_has_producer_without_rxn13531",
    "has_producer_without_rxn13531",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--networks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    queries = identity_queries(args.networks, load_identities())
    declarations = json.loads((args.networks / "corrections.json").read_text())[
        "participant_aliases"
    ]
    changes = {}
    for criteria in CRITERIA:
        for version in VERSIONS:
            network = load_network(args.networks, criteria, version)
            targets = [
                t
                for row in queries
                if row["criteria"] == criteria
                for t in row["targets"]
            ]
            changes[version, criteria] = query_changes(network, targets, declarations)
            del network
            print(
                f"{CRITERIA_NAMES[criteria]} {version}: {time.monotonic() - started:.0f} s",
                flush=True,
            )
    rows = query_rows(queries, changes)
    write_csv(args.out / "TableS1_tracing_native_control.csv", rows, QUERY_COLUMNS)
    write_csv(
        args.out / "TableS1_producing_reactions_native_control.csv",
        biomarker_rows(rows),
        BIOMARKER_COLUMNS,
    )
    stored_differs = sum(
        change["stored_removed"] != change["compounds_removed"]
        for values in changes.values()
        for change in values.values()
    )
    print(
        f"Wrote {len(rows)} query rows in {time.monotonic() - started:.0f} s; "
        f"queries where removed stored compounds differ from removed traced compounds: "
        f"{stored_differs}",
        flush=True,
    )


if __name__ == "__main__":
    main()
