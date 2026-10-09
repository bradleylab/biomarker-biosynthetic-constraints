"""Projected paths of the comparison set and the biomarkers (Tables 1, S1 and S2).

Builds the simplified compound graph from the MetaCyc flat files and the
reaction network (both versions: primary, and with conditional literature
connections), measures every compound of the comparison set, and writes:

  TableS1_projection.csv               the 86 biomarkers
  TableS2_projection.csv               the 10,473 compounds of the comparison set
  TableS2_population_summaries.csv     counts, medians and correlations
  Table1_projection_class_means.csv    class means and standard deviations
  Projection_graph_and_comparison_scope.csv   graph and comparison-set counts
  simplified_graph/                    product records and molecular weights read
                                       by the later scripts (derived from MetaCyc;
                                       not for distribution)

Usage:

  python analyses/02_projected_paths/projected_paths.py \\
      --metacyc METACYC_DIR --networks results/networks --out results/projected_paths
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import pandas as pd  # noqa: E402

from biosynthesis import paths, simplified_graph, statistics  # noqa: E402
from biosynthesis.panel import load_identities  # noqa: E402
from biosynthesis.tables import read_csv, write_csv, write_json  # noqa: E402

DATA = REPO / "data"
CLASS_LABELS = {"Triterpenoid": "Other triterpenoid"}
PATH_STATES = {0: "none", 1: "one"}
MEASUREMENT_COLUMNS = [
    "psi_bits",
    "cumulative_branching_score",
    "num_paths",
    "lambda",
    "max_depth",
    "sigma",
    "paths",
]
TABLE_S1_COLUMNS = [
    "network_version",
    "compound_id",
    "class",
    "name",
    "biomarker_number",
    "identity",
    *MEASUREMENT_COLUMNS,
    "queried_compound_id",
]
TABLE_S2_COLUMNS = [
    "network_version",
    "compound_id",
    *MEASUREMENT_COLUMNS,
    "queried_compound_id",
]
SUMMARY_COLUMNS = [
    "network_version",
    "population",
    "positive_path_n",
    "one_path_n",
    "multi_path_n",
    "multi_path_psi_median",
    "multi_path_psi_mean",
    "multi_path_psi_max",
    "multi_path_lambda_median",
    "multi_path_sigma_median",
    "few_paths_far_branching_n",
    "few_paths_far_branching_percent",
    "pearson_psi_lambda",
    "pearson_psi_sigma",
    "pearson_lambda_sigma",
]
CLASS_COLUMNS = [
    "network_version",
    "class",
    "n",
    "n_with_path",
    "mean_psi",
    "sd_psi",
    "mean_lambda",
    "sd_lambda",
    "mean_sigma",
    "sd_sigma",
]
SCOPE_COLUMNS = [
    "network_version",
    "comparison_compounds_n",
    "queried_compounds_n",
    "biomarkers_n",
    "descriptor_cohort_n",
    "products_including_empty_n",
    "products_with_records_n",
    "products_with_precursor_n",
    "graph_nodes_n",
    "graph_nodes_with_link_n",
    "graph_nodes_outside_comparison_set_n",
    "comparison_compounds_not_in_graph_n",
    "queried_compounds_not_in_graph_n",
]


def comparison_set() -> list[dict]:
    """The comparison set in its fixed order, with query and ordering columns."""
    return read_csv(DATA / "comparison_set.csv")


def biomarker_rows(comparison: list[dict]) -> list[dict]:
    """The 86 biomarker compounds in Table S1 order, with number, name, class and identity."""
    by_compound: dict[str, dict] = {}
    for record in load_identities():
        entry = by_compound.setdefault(
            record["compound_id"],
            {
                "compound_id": record["compound_id"],
                "biomarker_number": record["biomarker_number"],
                "name": record["name"],
                "class": CLASS_LABELS.get(record["class"], record["class"]),
                "identity": record["identity"],
            },
        )
        if entry["identity"] != record["identity"]:
            raise ValueError(f"Mixed identity states for {record['compound_id']}")
    ordered = sorted(
        (row for row in comparison if row["table_s1_row"]),
        key=lambda row: int(row["table_s1_row"]),
    )
    return [
        {
            **by_compound[row["compound_id"]],
            "queried_compound_id": row["queried_compound_id"],
        }
        for row in ordered
    ]


def measurement_cells(raw: dict) -> dict:
    """Table cells of one compound; compounds without a path have blank measurements."""
    count = raw["num_paths"]
    has_path = count > 0
    return {
        "psi_bits": paths.psi_bits(count),
        "cumulative_branching_score": raw["total_complexity"] if has_path else None,
        "num_paths": int(count),
        "lambda": raw["mean_branching_depth_normalized"] if has_path else None,
        "max_depth": int(raw["max_depth"]) if has_path else None,
        "sigma": raw["fraction_shared_80"] if has_path else None,
        "paths": PATH_STATES.get(int(count), "several"),
    }


def graph_scope(
    version: str,
    products: dict,
    network: dict,
    comparison: list[dict],
    biomarkers: int,
    descriptor_cohort: int,
) -> dict:
    """Counts of products, graph nodes and comparison compounds outside the graph."""
    links = simplified_graph.precursor_links(network)
    selected_products = {key for key, rows in network.items() if rows}
    nodes = selected_products | {
        value for rows in network.values() for row in rows for value in row["reactant"]
    }
    compounds = {row["compound_id"] for row in comparison}
    queried = {row["queried_compound_id"] for row in comparison}
    return {
        "network_version": version,
        "comparison_compounds_n": len(compounds),
        "queried_compounds_n": len(queried),
        "biomarkers_n": biomarkers,
        "descriptor_cohort_n": descriptor_cohort,
        "products_including_empty_n": len(products),
        "products_with_records_n": sum(bool(rows) for rows in products.values()),
        "products_with_precursor_n": len(selected_products),
        "graph_nodes_n": len(nodes),
        "graph_nodes_with_link_n": len(
            set(links) | {value for rows in links.values() for value in rows}
        ),
        "graph_nodes_outside_comparison_set_n": len(nodes - compounds),
        "comparison_compounds_not_in_graph_n": len(compounds - nodes),
        "queried_compounds_not_in_graph_n": len(queried - nodes),
    }


def descriptor_cohort_size() -> int:
    from biosynthesis.descriptors import select_cohort

    records = pd.read_csv(DATA / "descriptors.csv", dtype={"compound_id": str})
    return int(select_cohort(records).selected.sum())


def main(metacyc: Path, networks: Path, out: Path, workers: int) -> None:
    started = time.monotonic()
    comparison = comparison_set()
    biomarkers = biomarker_rows(comparison)
    products, masses, summary = simplified_graph.build_product_records(
        metacyc, networks
    )
    networks_built, selection = {}, {}
    for version in simplified_graph.NETWORK_VERSIONS:
        networks_built[version], selection[version] = (
            simplified_graph.select_precursors(products[version], masses)
        )
    graph_dir = out / "simplified_graph"
    for version in simplified_graph.NETWORK_VERSIONS:
        write_json(graph_dir / f"products_{version}.json", products[version])
    write_json(graph_dir / "molecular_weights.json", masses)
    write_json(
        graph_dir / "graph_summary.json", {**summary, "precursor_selection": selection}
    )
    print(f"Graphs built in {time.monotonic() - started:.0f} s", flush=True)

    queries = sorted({row["queried_compound_id"] for row in comparison})
    measured = paths.measure_targets(networks_built, queries, workers)
    print(
        f"{len(queries)} queries measured in {time.monotonic() - started:.0f} s",
        flush=True,
    )

    table_s1, table_s2, summaries, classes, scopes = [], [], [], [], []
    cohort_size = descriptor_cohort_size()
    by_extraction = sorted(comparison, key=lambda row: int(row["extraction_order"]))
    for version in simplified_graph.NETWORK_VERSIONS:
        raw = measured[version]
        for row in comparison:
            table_s2.append(
                {
                    "network_version": version,
                    "compound_id": row["compound_id"],
                    **measurement_cells(raw[row["queried_compound_id"]]),
                    "queried_compound_id": row["queried_compound_id"],
                }
            )
        for row in biomarkers:
            table_s1.append(
                {
                    "network_version": version,
                    **row,
                    **measurement_cells(raw[row["queried_compound_id"]]),
                }
            )
        # Population statistics follow extraction_order of data/comparison_set.csv.
        frame = paths.metric_frame(
            [
                {"compound_id": row["compound_id"], **raw[row["queried_compound_id"]]}
                for row in by_extraction
            ]
        )
        for population, subset in (
            ("all", frame),
            (
                "score_at_most_1e6",
                frame.loc[frame.cumulative_branching_score.le(statistics.SCORE_LIMIT)],
            ),
        ):
            summaries.append(
                {
                    "network_version": version,
                    "population": population,
                    **statistics.population_summary(subset),
                }
            )
        panel = pd.DataFrame(
            [
                {"class": row["class"], "compound_id": row["compound_id"]}
                for row in biomarkers
            ]
        ).merge(
            frame[["compound_id", "psi_bits", "num_paths", "lambda", "sigma"]],
            on="compound_id",
            how="left",
            validate="1:1",
        )
        classes.extend(
            {"network_version": version, **row}
            for row in statistics.class_summaries(panel)
        )
        scopes.append(
            graph_scope(
                version,
                products[version],
                networks_built[version],
                comparison,
                len(biomarkers),
                cohort_size,
            )
        )
    write_csv(out / "TableS1_projection.csv", table_s1, TABLE_S1_COLUMNS)
    write_csv(out / "TableS2_projection.csv", table_s2, TABLE_S2_COLUMNS)
    write_csv(out / "TableS2_population_summaries.csv", summaries, SUMMARY_COLUMNS)
    write_csv(out / "Table1_projection_class_means.csv", classes, CLASS_COLUMNS)
    write_csv(out / "Projection_graph_and_comparison_scope.csv", scopes, SCOPE_COLUMNS)
    print(f"Done in {time.monotonic() - started:.0f} s")


def read_projection(path: Path) -> dict[str, list[dict]]:
    """Rows of a projection table by network version, measurements as numbers."""
    result: dict[str, list[dict]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for column in ("psi_bits", "cumulative_branching_score", "lambda", "sigma"):
                row[column] = float(row[column]) if row[column] else float("nan")
            row["num_paths"] = int(row["num_paths"])
            result.setdefault(row["network_version"], []).append(row)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--metacyc", type=Path, required=True, help="MetaCyc 29.5 flat files"
    )
    parser.add_argument(
        "--networks", type=Path, required=True, help="reaction network folder"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    main(args.metacyc, args.networks, args.out, args.workers)
