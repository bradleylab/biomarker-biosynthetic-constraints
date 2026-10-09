"""Matched controls for upstream tracing (Table S9).

Each control changes one tracing rule: stop at the designated starting compounds,
follow only reactions placed in a MetaCyc pathway (with and without stopping),
follow a single pathway, or, under broader criteria, accept only reactions that
meet the strict criteria. A further check supplies the starting-compound and
pathway annotations without using them. Pathway placements and the starting
compounds are taken from the MetaCyc networks before literature records and
corrections are added.

    python analyses/07_matched_controls/matched_controls.py --metacyc $METACYC \\
        --networks results/networks --out results/controls
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

from biosynthesis.controls import (  # noqa: E402
    CONTROL_NAMES,
    UNCONTROLLED,
    alias_roots,
    annotations_change_result,
    controls,
    pathway_ids,
    pathway_masks,
    starting_identifiers,
    stopped_states,
    strict_occurrences,
)
from biosynthesis.graph import (  # noqa: E402
    CRITERIA,
    CRITERIA_NAMES,
    STRICT,
    VERSIONS,
    identity_queries,
    load_base_network,
    load_network,
    node_label,
)
from biosynthesis.panel import DATA, load_identities  # noqa: E402
from biosynthesis.tables import read_csv, write_csv  # noqa: E402

QUERY_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "query_node",
    "node",
    "control",
    "compounds_traced",
    "compounds_stored",
    "reactions_traced",
    "reactions_removed_by_control",
    "membership_changed_by_control",
)
CONTEXT_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "query_node",
    "node",
    "pathway",
    "pathway_reaction_count",
    "pathway_and_stop_reaction_count",
    "reactions_removed_by_pathway_restriction",
    "additional_reactions_removed_by_stopping",
)
ANNOTATION_COLUMNS = (
    "network_version",
    "reaction_criteria",
    "query_node",
    "node",
    "annotation_changes_result",
)

_STATE = None


def evaluate(root: str) -> tuple[dict, bool]:
    """Run all controls for one traced state and test whether unused annotations change it."""
    network, stops, masks, strict_accepted = _STATE
    result = controls(network, root, stops, masks, strict_accepted)
    changed = annotations_change_result(network, root, result[UNCONTROLLED], stops)
    return result, changed


def run_controls(state: tuple, roots: list[str], workers: int) -> dict:
    """Evaluate each distinct traced state once, in parallel worker processes when asked."""
    global _STATE
    _STATE = state
    gc.freeze()  # the network is read only; keep the garbage collector off its objects
    try:
        if workers <= 1:
            return {root: evaluate(root) for root in roots}
        context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            return dict(zip(roots, pool.map(evaluate, roots, chunksize=1)))
    finally:
        _STATE = None
        gc.unfreeze()


def rows_for(
    version: str, criteria: str, network, targets: list[str], roots: dict, results: dict
) -> tuple[list, list, list]:
    """Table rows of one network version and criteria, in the order of the query identifiers."""
    queries, contexts, annotations = [], [], []
    for target in sorted(targets):
        root = roots.get(target, target)
        result, changed = results[root]
        meta = {
            "network_version": version,
            "reaction_criteria": CRITERIA_NAMES[criteria],
            "query_node": node_label(network.nodes[target]),
            "node": node_label(network.nodes[root]),
        }
        uncut = result[UNCONTROLLED]
        for name in CONTROL_NAMES:
            if name not in result:
                continue
            controlled = result[name]
            queries.append(
                {
                    **meta,
                    "control": name,
                    "compounds_traced": len(controlled["required_nodes"]),
                    "compounds_stored": len(controlled["stored_nodes"]),
                    "reactions_traced": len(controlled["edge_ids"]),
                    "reactions_removed_by_control": len(uncut["edge_ids"])
                    - len(controlled["edge_ids"]),
                    "membership_changed_by_control": controlled != uncut,
                }
            )
        for pathway, context in result["pathway_contexts"].items():
            only, stopped = context["pathway_only"], context["pathway_and_stop"]
            contexts.append(
                {
                    **meta,
                    "pathway": pathway,
                    "pathway_reaction_count": len(only["edge_ids"]),
                    "pathway_and_stop_reaction_count": len(stopped["edge_ids"]),
                    "reactions_removed_by_pathway_restriction": len(uncut["edge_ids"])
                    - len(only["edge_ids"]),
                    "additional_reactions_removed_by_stopping": len(only["edge_ids"])
                    - len(stopped["edge_ids"]),
                }
            )
        annotations.append({**meta, "annotation_changes_result": changed})
    return queries, contexts, annotations


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metacyc", type=Path, required=True)
    parser.add_argument("--networks", type=Path, required=True)
    parser.add_argument(
        "--starting-compounds", type=Path, default=DATA / "starting_compounds.csv"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    started = time.monotonic()
    queries = identity_queries(args.networks, load_identities())
    declarations = json.loads((args.networks / "corrections.json").read_text())[
        "participant_aliases"
    ]
    pathways = pathway_ids(args.metacyc / "pathways.dat")
    starting = [row["compound_id"] for row in read_csv(args.starting_compounds)]
    base_strict = load_base_network(args.networks, STRICT)
    identifiers = starting_identifiers(base_strict, starting)
    strict_accepted = strict_occurrences(base_strict)
    del base_strict
    tables = ([], [], [])
    blocks = {}
    for criteria in CRITERIA:
        base = load_base_network(args.networks, criteria)
        stops = stopped_states(base, identifiers)
        base_masks = pathway_masks(base, pathways)
        del base
        targets = sorted(
            {
                t
                for row in queries
                if row["criteria"] == criteria
                for t in row["targets"]
            }
        )
        for version in VERSIONS:
            network = load_network(args.networks, criteria, version)
            masks = {
                label: set(values) & network.edges.keys()
                for label, values in base_masks.items()
            }
            roots = alias_roots(network, declarations, set(targets))
            unique = sorted({roots.get(target, target) for target in targets})
            results = run_controls(
                (network, stops, masks, strict_accepted), unique, args.workers
            )
            blocks[version, criteria] = rows_for(
                version, criteria, network, targets, roots, results
            )
            del network, results
            print(
                f"{CRITERIA_NAMES[criteria]} {version}: {len(unique)} traced states, "
                f"{len(stops)} stopped states, {time.monotonic() - started:.0f} s",
                flush=True,
            )
    for criteria in CRITERIA:
        for version in VERSIONS:
            for table, rows in zip(tables, blocks[version, criteria]):
                table.extend(rows)
    write_csv(args.out / "TableS9_per_query_controls.csv", tables[0], QUERY_COLUMNS)
    write_csv(args.out / "TableS9_per_context_controls.csv", tables[1], CONTEXT_COLUMNS)
    write_csv(
        args.out / "TableS9_annotation_invariance.csv", tables[2], ANNOTATION_COLUMNS
    )
    print(
        f"Wrote {len(tables[2])} queries in {time.monotonic() - started:.0f} s",
        flush=True,
    )


if __name__ == "__main__":
    main()
