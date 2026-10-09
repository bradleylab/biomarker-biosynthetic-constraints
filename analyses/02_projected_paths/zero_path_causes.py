"""Why compounds have no projected path, and what the two hub exemptions change (Table S2).

Reads the primary simplified graph written by projected_paths.py, recounts paths
with beta-carotene and 2,3-oxidosqualene treated like any other potential hub, and
records why each compound without a path in the primary graph has none. Writes:

  TableS2_zero_path_causes.csv
  TableS2_exemption_sensitivity.csv
  zero_path_summary.json   the counts quoted in the text, including the split of
                           the queried compounds that are not graph nodes

  python analyses/02_projected_paths/zero_path_causes.py \\
      --paths results/projected_paths --out results/projected_paths
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis import paths, simplified_graph  # noqa: E402
from biosynthesis.tables import read_csv, write_csv, write_json  # noqa: E402

PRIMARY = "primary"
NOT_A_PRODUCT = "not_a_product"
NO_RECORDS = "no_producing_records"
ALL_REACTANTS_EXCLUDED = "every_producing_record_lacks_an_eligible_reactant"
SELF_LINKS_ONLY = "only_self_links_remain"
EXEMPTION_COLUMNS = [
    "compound_id",
    "queried_compound_id",
    "num_paths_primary_graph",
    "num_paths_without_exemptions",
]
CAUSE_COLUMNS = [
    "compound_id",
    "queried_compound_id",
    "cause",
    "producing_records",
    "excluded_reactants",
]


def zero_path_cause(
    query: str, products: dict, network: dict, precursors: dict
) -> tuple[str, list[str]]:
    """Why a compound has no projected path, with the reactants its records lost."""
    if query not in products:
        return NOT_A_PRODUCT, []
    if not products[query]:
        return NO_RECORDS, []
    if query not in network:
        lost = sorted({r for row in products[query] for r in row.get("reactant", [])})
        return ALL_REACTANTS_EXCLUDED, lost
    if not precursors.get(query):
        return SELF_LINKS_ONLY, []
    raise ValueError(f"Compound has a selected precursor but no path: {query}")


def cause_rows(labels: list[dict], products: dict, network: dict) -> list[dict]:
    """One row per compound of the comparison set without a path."""
    precursors = simplified_graph.precursor_links(network)
    rows = []
    for row in labels:
        if int(row["num_paths"]) != 0:
            continue
        query = row["queried_compound_id"]
        cause, lost = zero_path_cause(query, products, network, precursors)
        rows.append(
            {
                "compound_id": row["compound_id"],
                "queried_compound_id": query,
                "cause": cause,
                "producing_records": len(products.get(query, [])),
                "excluded_reactants": ";".join(lost),
            }
        )
    return rows


def summarize(
    labels: list[dict],
    with_exemptions: dict[str, int],
    without: dict[str, int],
    panel: list[dict],
    causes: list[dict],
    products: dict,
    network: dict,
) -> dict:
    """Counts over distinct queries, comparison-set compounds and biomarkers."""
    changed = {q for q in with_exemptions if with_exemptions[q] != without[q]}
    panel_changed = [row for row in panel if row["queried_compound_id"] in changed]
    nodes = set(network) | {
        value for rows in network.values() for row in rows for value in row["reactant"]
    }
    precursors = simplified_graph.precursor_links(network)
    absent = sorted(set(with_exemptions) - nodes)
    absent_causes = Counter(
        zero_path_cause(query, products, network, precursors)[0] for query in absent
    )
    return {
        "queried_compounds": len(with_exemptions),
        "queried_compounds_changed": len(changed),
        "queried_compounds_fewer_paths": sum(
            without[q] < with_exemptions[q] for q in changed
        ),
        "queried_compounds_more_paths": sum(
            without[q] > with_exemptions[q] for q in changed
        ),
        "queried_compounds_losing_every_path": sum(
            with_exemptions[q] > 0 and without[q] == 0 for q in changed
        ),
        "comparison_compounds": len(labels),
        "comparison_compounds_fewer_paths": sum(
            without[row["queried_compound_id"]]
            < with_exemptions[row["queried_compound_id"]]
            for row in labels
        ),
        "comparison_compounds_more_paths": sum(
            without[row["queried_compound_id"]]
            > with_exemptions[row["queried_compound_id"]]
            for row in labels
        ),
        "comparison_compounds_losing_every_path": sum(
            with_exemptions[row["queried_compound_id"]] > 0
            and without[row["queried_compound_id"]] == 0
            for row in labels
        ),
        "biomarkers": len(panel),
        "biomarkers_fewer_paths": len(panel_changed),
        "biomarkers_fewer_paths_by_class": dict(
            sorted(Counter(row["class"] for row in panel_changed).items())
        ),
        "biomarker_class_sizes": dict(
            sorted(Counter(row["class"] for row in panel).items())
        ),
        "biomarkers_losing_every_path": sorted(
            row["name"]
            for row in panel_changed
            if without[row["queried_compound_id"]] == 0
        ),
        "biomarkers_losing_every_path_counts": sorted(
            {
                with_exemptions[row["queried_compound_id"]]
                for row in panel_changed
                if without[row["queried_compound_id"]] == 0
            }
        ),
        "zero_path_compounds": len(causes),
        "zero_path_compounds_by_cause": dict(
            sorted(Counter(r["cause"] for r in causes).items())
        ),
        "queried_compounds_not_in_graph": len(absent),
        "queried_compounds_not_in_graph_by_cause": dict(sorted(absent_causes.items())),
    }


def main(projection: Path, out: Path, workers: int) -> None:
    graph_dir = projection / "simplified_graph"
    with open(graph_dir / f"products_{PRIMARY}.json", encoding="utf-8") as handle:
        products = json.load(handle)
    with open(graph_dir / "molecular_weights.json", encoding="utf-8") as handle:
        masses = json.load(handle)
    network, _ = simplified_graph.select_precursors(products, masses)
    without_exemptions, _ = simplified_graph.select_precursors(
        products, masses, hub_exemptions=frozenset()
    )
    labels = [
        row
        for row in read_csv(projection / "TableS2_projection.csv")
        if row["network_version"] == PRIMARY
    ]
    panel = [
        row
        for row in read_csv(projection / "TableS1_projection.csv")
        if row["network_version"] == PRIMARY
    ]
    with_exemptions = {
        row["queried_compound_id"]: int(row["num_paths"]) for row in labels
    }
    queries = sorted(with_exemptions)
    counted = paths.measure_targets({"without": without_exemptions}, queries, workers)
    without = {query: int(counted["without"][query]["num_paths"]) for query in queries}
    exemptions = [
        {
            "compound_id": row["compound_id"],
            "queried_compound_id": row["queried_compound_id"],
            "num_paths_primary_graph": with_exemptions[row["queried_compound_id"]],
            "num_paths_without_exemptions": without[row["queried_compound_id"]],
        }
        for row in labels
    ]
    causes = cause_rows(labels, products, network)
    write_csv(out / "TableS2_exemption_sensitivity.csv", exemptions, EXEMPTION_COLUMNS)
    write_csv(out / "TableS2_zero_path_causes.csv", causes, CAUSE_COLUMNS)
    summary = summarize(
        labels, with_exemptions, without, panel, causes, products, network
    )
    write_json(out / "zero_path_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--paths", type=Path, required=True, help="output of projected_paths.py"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    main(args.paths, args.out, args.workers)
