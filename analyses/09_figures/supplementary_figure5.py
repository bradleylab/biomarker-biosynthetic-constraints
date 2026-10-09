"""Supplementary Figure 5: psi, lambda and sigma of the compounds with several paths.

Uses the compounds with several paths and a cumulative branching score of at most
10^6. Scatter panels give the Pearson correlations of Table S2; histogram dashed
lines mark medians. Writes SupplementaryFigure5.{png,pdf} (primary graph),
SupplementaryFigure5_with_conditional.{png,pdf} and SupplementaryFigure5.json.

  python analyses/09_figures/supplementary_figure5.py --results results --out results/figures
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from biosynthesis.simplified_graph import NETWORK_VERSIONS  # noqa: E402
from biosynthesis.statistics import SCORE_LIMIT  # noqa: E402
from biosynthesis.tables import read_csv, write_json  # noqa: E402

DATA = REPO / "data"
BLUE = "#3274A1"
NAME = "SupplementaryFigure5"
LABELS = {"psi_bits": "ψ (bits)", "lambda": "λ", "sigma": "σ"}
PAIRS = [("psi_bits", "lambda"), ("psi_bits", "sigma"), ("lambda", "sigma")]
CORRELATION_COLUMNS = {
    ("psi_bits", "lambda"): "pearson_psi_lambda",
    ("psi_bits", "sigma"): "pearson_psi_sigma",
    ("lambda", "sigma"): "pearson_lambda_sigma",
}


def file_name(version: str) -> str:
    return NAME if version == "primary" else f"{NAME}_{version}"


def save(fig, out: Path, name: str) -> None:
    for extension in ("png", "pdf"):
        metadata = {"Title": name}
        if extension == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        fig.savefig(
            out / f"{name}.{extension}", dpi=220, bbox_inches="tight", metadata=metadata
        )
    plt.close(fig)


def plotted_compounds(projection: Path, version: str) -> pd.DataFrame:
    """Compounds with several paths and a score of at most 10^6, in extraction_order of data/comparison_set.csv."""
    order = {
        row["compound_id"]: int(row["extraction_order"])
        for row in read_csv(DATA / "comparison_set.csv")
    }
    rows = [
        row
        for row in read_csv(projection / "TableS2_projection.csv")
        if row["network_version"] == version
    ]
    rows.sort(key=lambda row: order[row["compound_id"]])
    table = pd.DataFrame(
        {
            "compound_id": [row["compound_id"] for row in rows],
            "num_paths": [int(row["num_paths"]) for row in rows],
            **{
                column: [
                    float(row[column]) if row[column] else float("nan") for row in rows
                ]
                for column in (
                    "psi_bits",
                    "lambda",
                    "sigma",
                    "cumulative_branching_score",
                )
            },
        }
    )
    return table[
        (table.num_paths > 1) & (table.cumulative_branching_score <= SCORE_LIMIT)
    ]


def draw(projection: Path, out: Path, version: str) -> dict:
    selected = plotted_compounds(projection, version)
    summary = next(
        row
        for row in read_csv(projection / "TableS2_population_summaries.csv")
        if row["network_version"] == version
        and row["population"] == "score_at_most_1e6"
    )
    if len(selected) != int(summary["multi_path_n"]):
        raise ValueError("Plotted compounds differ from the population summary")
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    values = {"n": len(selected), "pearson_r": {}, "median": {}}
    for ax, (x, y) in zip(axes[0], PAIRS):
        ax.scatter(
            selected[x], selected[y], s=4, alpha=0.2, color=BLUE, rasterized=True
        )
        r = float(summary[CORRELATION_COLUMNS[(x, y)]])
        ax.set(xlabel=LABELS[x], ylabel=LABELS[y], title=f"Pearson r = {r:.3f}")
        values["pearson_r"][f"{x}_{y}"] = r
    for ax, field in zip(axes[1], LABELS):
        median = float(selected[field].median())
        ax.hist(selected[field], bins="auto", color=BLUE)
        ax.axvline(median, color="black", ls="--")
        ax.set(xlabel=LABELS[field], ylabel="Compounds", title=f"Median = {median:.3f}")
        values["median"][field] = median
    fig.tight_layout()
    save(fig, out, file_name(version))
    values["compounds"] = selected[
        ["compound_id", "psi_bits", "lambda", "sigma"]
    ].to_dict("list")
    return values


def main(results: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    projection = results / "projected_paths"
    values = {version: draw(projection, out, version) for version in NETWORK_VERSIONS}
    write_json(out / f"{NAME}.json", values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.results, args.out)
