"""The descriptor cohort and its selection (Table S4).

The cohort is selected from the 7,380 compounds with descriptor records by the
path counts and cumulative branching scores of the reference version of the
simplified graph (the reference graph; values in data/descriptors.csv). Writes:

  TableS4_descriptor_selection.csv   every compound with a descriptor record
  TableS4_descriptor_cohort.csv      descriptor values, Böttcher status and the
                                     projected-path measurements of both graph versions

  python analyses/02_projected_paths/descriptor_cohort.py \\
      --paths results/projected_paths --out results/descriptors
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import pandas as pd  # noqa: E402

from biosynthesis.descriptors import select_cohort  # noqa: E402
from biosynthesis.simplified_graph import NETWORK_VERSIONS  # noqa: E402
from biosynthesis.statistics import FEATURES  # noqa: E402
from biosynthesis.tables import read_csv, write_csv  # noqa: E402

DATA = REPO / "data"
SELECTION_COLUMNS = [
    "compound_id",
    "num_paths_reference_graph",
    "cumulative_score_reference_graph",
    "has_path_reference_graph",
    "selected",
]
COHORT_COLUMNS = [
    "network_version",
    "compound_id",
    *FEATURES,
    "psi_bits",
    "bottchscore_unscreened",
    "bottcher_status",
    "bottcher_warnings",
    "lambda",
    "sigma",
    "num_paths",
    "queried_compound_id",
]


def descriptor_records() -> list[dict]:
    """Rows of data/descriptors.csv with every value kept as written."""
    return read_csv(DATA / "descriptors.csv")


def selection_table(records: list[dict]) -> pd.DataFrame:
    """The selection criteria of every compound with a descriptor record."""
    frame = pd.DataFrame(
        {
            "compound_id": [row["compound_id"] for row in records],
            "num_paths_reference_graph": pd.to_numeric(
                [row["num_paths_reference_graph"] for row in records]
            ),
            "cumulative_score_reference_graph": [
                float(row["cumulative_score_reference_graph"]) for row in records
            ],
        }
    )
    return select_cohort(frame)


def cohort_rows(
    records: list[dict], selected: set[str], projection: Path
) -> list[dict]:
    """Cohort rows of both graph versions, in cohort order."""
    responses: dict[str, dict[str, dict]] = {}
    for row in read_csv(projection / "TableS2_projection.csv"):
        responses.setdefault(row["network_version"], {})[row["compound_id"]] = row
    rows = []
    for version in NETWORK_VERSIONS:
        for record in records:
            if record["compound_id"] not in selected:
                continue
            measured = responses[version][record["compound_id"]]
            rows.append(
                {
                    "network_version": version,
                    "compound_id": record["compound_id"],
                    **{feature: record[feature] for feature in FEATURES},
                    "psi_bits": measured["psi_bits"],
                    "bottchscore_unscreened": record["bottchscore_unscreened"],
                    "bottcher_status": record["bottcher_status"],
                    "bottcher_warnings": record["bottcher_warnings"],
                    "lambda": measured["lambda"],
                    "sigma": measured["sigma"],
                    "num_paths": measured["num_paths"],
                    "queried_compound_id": measured["queried_compound_id"],
                }
            )
    return rows


def main(projection: Path, out: Path) -> None:
    records = descriptor_records()
    selection = selection_table(records)
    write_csv(
        out / "TableS4_descriptor_selection.csv",
        (
            {
                "compound_id": row.compound_id,
                "num_paths_reference_graph": int(row.num_paths_reference_graph),
                "cumulative_score_reference_graph": float(
                    row.cumulative_score_reference_graph
                ),
                "has_path_reference_graph": bool(row.has_path_reference_graph),
                "selected": bool(row.selected),
            }
            for row in selection.itertuples()
        ),
        SELECTION_COLUMNS,
    )
    selected = set(selection.loc[selection.selected, "compound_id"])
    rows = cohort_rows(records, selected, projection)
    write_csv(out / "TableS4_descriptor_cohort.csv", rows, COHORT_COLUMNS)
    print(f"{len(selection)} descriptor records, {len(selected)} selected")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--paths", type=Path, required=True, help="output of projected_paths.py"
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.paths, args.out)
