"""Table 1: projected-path measurements by biomarker class (primary graph), as Markdown.

Reads Table1_projection_class_means.csv and writes Table1.md with the number of
biomarkers and the number with a path in each class, and the mean (SD) of psi,
lambda and sigma. A dash
marks a class with fewer than two evaluable records.

  python analyses/09_figures/table1.py --results results --out results/projected_paths
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.tables import read_csv  # noqa: E402

CAPTION = (
    "**Table 1. Projected-path measurements by biomarker class.** Total and evaluable "
    "records are shown with means and standard deviations (SD) for path multiplicity (ψ), "
    "normalized branch depth (λ) and fraction shared (σ). Statistics include only "
    "biomarkers with at least one projected path; single-path biomarkers enter at ψ = 0, "
    "λ = 0 and σ = 1. A dash indicates fewer than two evaluable records."
)


def format_cell(mean: str, sd: str) -> str:
    return f"{float(mean):.2f} ({'—' if sd == '' else f'{float(sd):.2f}'})"


def table1(rows: list[dict]) -> list[str]:
    lines = [
        (
            "| Biomarker class | Total records | Evaluable records | ψ (bits), mean (SD) | "
            "λ, mean (SD) | σ, mean (SD) |"
        ),
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['class']} | {row['n']} | {row['n_with_path']} | "
            f"{format_cell(row['mean_psi'], row['sd_psi'])} | "
            f"{format_cell(row['mean_lambda'], row['sd_lambda'])} | "
            f"{format_cell(row['mean_sigma'], row['sd_sigma'])} |"
        )
    return lines


def main(results: Path, out: Path) -> None:
    rows = [
        row
        for row in read_csv(
            results / "projected_paths/Table1_projection_class_means.csv"
        )
        if row["network_version"] == "primary"
    ]
    lines = table1(rows)
    out.mkdir(parents=True, exist_ok=True)
    (out / "Table1.md").write_text(
        CAPTION + "\n\n" + "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("\n".join(lines[2:]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.results, args.out)
