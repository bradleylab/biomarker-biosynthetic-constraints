"""Table S5: the values of the Brocks et al. (2023) workbooks, one row per sample record and fossil class.

Usage:
  python analyses/05_source_compatibility/read_workbooks.py --workbooks DIR --out results/sources

--workbooks is the folder holding brocks2023_SI1_MOESM4.xlsx and
brocks2023_SI2_MOESM5.xlsx (Supplementary Information of Brocks et al., 2023).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.tables import write_csv  # noqa: E402
from biosynthesis.workbooks import (  # noqa: E402
    OBSERVATION_COLUMNS,
    assign_samples,
    load_layout,
    read_crosswalk,
    read_workbooks,
)

TABLE_S5 = "TableS5_source_observations.csv"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--workbooks", type=Path, required=True, help="Folder with the two workbooks."
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder.")
    parser.add_argument(
        "--data", type=Path, default=REPO / "data", help="The repository data folder."
    )
    args = parser.parse_args()

    layout = load_layout(args.data / "geology" / "workbook_layout.yaml")
    rows = assign_samples(
        read_workbooks(layout, args.workbooks),
        layout,
        read_crosswalk(args.data / "geology" / "sample_crosswalk.csv"),
    )
    count = write_csv(args.out / TABLE_S5, rows, OBSERVATION_COLUMNS)
    samples = {row["canonical_sample_key"] for row in rows if row["included_sample"]}
    print(f"{TABLE_S5}: {count} rows; {len(samples)} samples.")


if __name__ == "__main__":
    main()
