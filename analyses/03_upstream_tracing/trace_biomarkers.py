"""Trace every biomarker identity upstream and enumerate its inclusion-minimal reaction sets.

Writes Table S1 (tracing and producing reactions), Table S6 (representative
targets) and one record per traced query, used by the supplementary figures.

    python analyses/03_upstream_tracing/trace_biomarkers.py --networks results/networks \\
        --pairs results/pairs --out results/tracing
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gc
import json
import multiprocessing
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.enumeration import MAX_SEARCH_STATES, enumerate_sets, reaction_graph  # noqa: E402
from biosynthesis.graph import (  # noqa: E402
    CRITERIA,
    CRITERIA_NAMES,
    VERSIONS,
    identity_queries,
    load_network,
    node_label,
    traced_roots,
)
from biosynthesis.panel import DATA, load_identities  # noqa: E402
from biosynthesis.tables import read_csv, write_csv, write_jsonl  # noqa: E402
from biosynthesis.tracing import (
    NO_PRODUCER,
    PRODUCER_EXCLUDED,
    classify_frontier,
    discover,
)  # noqa: E402

TRACING_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "biomarker_number",
    "name",
    "compound_id",
    "identity_id",
    "identity",
    "source_compound_ids",
    "node",
    "equivalent_compound_ids",
    "status",
    "producing_reactions",
    "compounds_traced",
    "reactions_traced",
    "compounds_without_producer",
    "cycles",
    "unidentified_participants",
    "reactions_with_unidentified_participants",
    "reaction_sets",
    "reaction_sets_found",
    "search_states_examined",
    "search_states_left",
    "sets_consuming_target",
    "coefficients_complete",
)
REPRESENTATIVE_COLUMNS = tuple(
    c for c in TRACING_COLUMNS if c != "equivalent_compound_ids"
)
PRODUCING_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "biomarker_number",
    "name",
    "identity",
    "identity_variants",
    "target_nodes_in_network",
    "has_producer",
    "accepted_identity_has_producer",
    "candidate_identity_has_producer",
    "metacyc_pathways",
)
NOT_IN_NETWORK = "not_in_network"

_NETWORK = None


def evaluate(root: str) -> dict:
    """Trace one query state and enumerate its reaction sets in the network of this process."""
    network = _NETWORK
    discovery = discover(network, [root])
    frontier = classify_frontier(network, discovery)
    enumeration = enumerate_sets(
        reaction_graph(network, root), max_states=MAX_SEARCH_STATES
    )
    statuses = [row["status"] for row in frontier["endpoints"]]
    return {
        "node_id": root,
        "complete": discovery["complete"],
        "compounds_traced": len(discovery["required_nodes"]),
        "reactions_traced": len(discovery["edge_ids"]),
        "compounds_without_producer": len(frontier["endpoints"]),
        "endpoints_no_producer": statuses.count(NO_PRODUCER),
        "endpoints_producer_excluded": statuses.count(PRODUCER_EXCLUDED),
        "cycles": len(frontier["cyclic_components"]),
        "unidentified_participants": len(frontier["unresolved_nodes"]),
        "reactions_with_unidentified_participants": len(
            frontier["unresolved_edge_requirements"]
        ),
        **{key: value for key, value in enumeration.items() if key != "reaction_sets"},
    }


def run_queries(network, roots: list[str], workers: int) -> dict[str, dict]:
    """Evaluate each distinct root once, in parallel worker processes when asked."""
    global _NETWORK
    _NETWORK = network
    gc.freeze()  # the network is read only; keep the garbage collector off its objects
    try:
        if workers <= 1:
            return {root: evaluate(root) for root in roots}
        context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            return dict(zip(roots, pool.map(evaluate, roots, chunksize=1)))
    finally:
        _NETWORK = None
        gc.unfreeze()


def representative_compounds(path: Path) -> set[str]:
    """The six target compounds of the representative precursor classes."""
    rows = read_csv(path)
    for column in ("compound_id", "target_local_id", "target_compound_id"):
        if column in rows[0]:
            return {row[column] for row in rows}
    raise ValueError("No compound column in " + str(path))


def pathway_counts(pairs: Path) -> dict[tuple[str, int], int]:
    """MetaCyc pathways containing a producing reaction of each biomarker, by criteria."""
    names = {name: criteria for criteria, name in CRITERIA_NAMES.items()}
    return {
        (names[row["reaction_criteria"]], int(row["biomarker_number"])): int(
            row["pathway_count"]
        )
        for row in read_csv(Path(pairs) / "biomarker_enzyme_sets.csv")
    }


def tracing_rows(queries, results, roots, labels, producers, aliases) -> list[dict]:
    """One Table S1 row per network version, criteria, identity and query state."""
    output = []
    for version in VERSIONS:
        for criteria in CRITERIA:
            for query in (row for row in queries if row["criteria"] == criteria):
                identity = query["identity"]
                for target in query["targets"] or [None]:
                    row = {
                        "network_version": version,
                        "reaction_criteria": CRITERIA_NAMES[criteria],
                        "biomarker_number": identity["biomarker_number"],
                        "name": identity["name"],
                        "compound_id": identity["compound_id"],
                        "identity_id": identity["identity_id"],
                        "identity": identity["identity"],
                        "source_compound_ids": identity["metacyc_compounds"],
                        "producing_reactions": 0,
                    }
                    if target is None:
                        row["status"] = NOT_IN_NETWORK
                    else:
                        root = roots[version, criteria][target]
                        result = results[version, criteria][root]
                        row.update(
                            node=labels[version, criteria][root],
                            equivalent_compound_ids=aliases[version, criteria].get(
                                target
                            ),
                            status=result["status"],
                            producing_reactions=producers[version, criteria][root],
                            reaction_sets=result["sets"],
                            reaction_sets_found=result["sets_found"],
                            **{
                                key: result[key]
                                for key in TRACING_COLUMNS
                                if key in result and key not in ("status",)
                            },
                        )
                    output.append(row)
    return output


def producing_rows(tracing: list[dict], pathways: dict) -> list[dict]:
    """One row per network version, criteria and biomarker: does any identity have a producer."""
    output = []
    for version in VERSIONS:
        for criteria in CRITERIA:
            name = CRITERIA_NAMES[criteria]
            for number in sorted({row["biomarker_number"] for row in tracing}):
                chosen = [
                    row
                    for row in tracing
                    if row["network_version"] == version
                    and row["reaction_criteria"] == name
                    and row["biomarker_number"] == number
                ]
                states = {row["identity"] for row in chosen}
                accepted = {state == "accepted" for state in states}
                if len(accepted) != 1:
                    raise ValueError("Accepted and other identities in one biomarker")
                accepted = accepted.pop()
                producer = any(row["producing_reactions"] > 0 for row in chosen)
                output.append(
                    {
                        "network_version": version,
                        "reaction_criteria": name,
                        "biomarker_number": number,
                        "name": chosen[0]["name"],
                        "identity": "/".join(sorted(states)),
                        "identity_variants": len(
                            {row["identity_id"] for row in chosen}
                        ),
                        "target_nodes_in_network": len(
                            {row["node"] for row in chosen if row.get("node")}
                        ),
                        "has_producer": producer,
                        "accepted_identity_has_producer": accepted and producer,
                        "candidate_identity_has_producer": not accepted and producer,
                        "metacyc_pathways": pathways[criteria, number],
                    }
                )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--networks", type=Path, required=True)
    parser.add_argument(
        "--pairs",
        type=Path,
        required=True,
        help="folder with biomarker_enzyme_sets.csv",
    )
    parser.add_argument(
        "--representatives", type=Path, default=DATA / "representative_targets.csv"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    started = time.monotonic()
    queries = identity_queries(args.networks, load_identities())
    declarations = json.loads((args.networks / "corrections.json").read_text())[
        "participant_aliases"
    ]
    results, roots, labels, producers, aliases, records = {}, {}, {}, {}, {}, []
    for criteria in CRITERIA:
        for version in VERSIONS:
            network = load_network(args.networks, criteria, version)
            targets = [
                t
                for row in queries
                if row["criteria"] == criteria
                for t in row["targets"]
            ]
            mapping = traced_roots(network, declarations, targets)
            unique = sorted(set(mapping.values()))
            key = version, criteria
            results[key] = run_queries(network, unique, args.workers)
            roots[key] = mapping
            labels[key] = {root: node_label(network.nodes[root]) for root in unique}
            producers[key] = {
                root: len(network.producers.get(root, ())) for root in unique
            }
            aliases[key] = {
                target: [
                    network.nodes[target]["identity"]["local_id"],
                    network.nodes[root]["identity"]["local_id"],
                ]
                for target, root in mapping.items()
                if root != target
            }
            for root in unique:
                records.append(
                    {
                        "network_version": version,
                        "reaction_criteria": CRITERIA_NAMES[criteria],
                        "node": labels[key][root],
                        **results[key][root],
                    }
                )
            print(
                f"{CRITERIA_NAMES[criteria]} {version}: {len(unique)} queries, "
                f"{time.monotonic() - started:.0f} s",
                flush=True,
            )
            del network
    for criteria in CRITERIA:
        named = {}
        for version in VERSIONS:
            for state, label in labels[version, criteria].items():
                if named.setdefault(label, state) != state:
                    raise ValueError(
                        "Two traced states share a readable label: " + label
                    )
    tracing = tracing_rows(queries, results, roots, labels, producers, aliases)
    write_csv(args.out / "TableS1_tracing.csv", tracing, TRACING_COLUMNS)
    targets = representative_compounds(args.representatives)
    representatives = [row for row in tracing if row["compound_id"] in targets]
    if {row["compound_id"] for row in representatives} != targets:
        raise ValueError("A representative target is not a biomarker compound")
    write_csv(
        args.out / "TableS6_representative_targets.csv",
        representatives,
        REPRESENTATIVE_COLUMNS,
    )
    write_csv(
        args.out / "TableS1_producing_reactions.csv",
        producing_rows(tracing, pathway_counts(args.pairs)),
        PRODUCING_COLUMNS,
    )
    write_jsonl(args.out / "query_results.jsonl", records)
    print(
        f"Wrote {len(tracing)} tracing rows and {len(records)} query records "
        f"in {time.monotonic() - started:.0f} s",
        flush=True,
    )


if __name__ == "__main__":
    main()
