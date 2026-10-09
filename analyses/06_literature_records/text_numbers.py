"""Print the numbers of the text that come from the Table S8 files on MetaCyc reactions.

python analyses/06_literature_records/text_numbers.py --results results
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.tables import read_csv  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--results", type=Path, required=True, help="results folder, e.g. results"
    )
    args = parser.parse_args()
    corrections = read_csv(
        args.results / "literature" / "TableS8_metacyc_corrections.csv"
    )
    test = read_csv(args.results / "literature" / "TableS8_native_control.csv")
    changed = sum(row["action"] in ("Replaced", "Excluded") for row in corrections)
    retained = sum(row["action"].startswith("Retained") for row in corrections)
    print(
        f"MetaCyc reactions replaced or removed (main text, Results and Methods): {changed}"
    )
    print(
        f"MetaCyc reactions retained with unresolved requirements recorded (Methods): {retained}"
    )
    print(
        f"MetaCyc hopanoid reactions in the single-reaction test (Supplementary Data legend): {len(test)}"
    )


if __name__ == "__main__":
    main()
