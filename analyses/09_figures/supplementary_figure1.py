"""Supplementary Figure 1: the 13 molecular descriptors against path multiplicity psi.

Each panel shows every compound of the descriptor cohort with both values, with
the Spearman correlation and its denominator from Table S4 (main analysis, all
cohort compounds). Writes SupplementaryFigure1.{png,pdf} (primary graph),
SupplementaryFigure1_with_conditional.{png,pdf} and SupplementaryFigure1.json.

  python analyses/09_figures/supplementary_figure1.py --results results --out results/figures
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
from biosynthesis.statistics import FEATURES, numeric  # noqa: E402
from biosynthesis.tables import read_csv, write_json  # noqa: E402

BLUE = "#245f8e"
NAME = "SupplementaryFigure1"
LABELS = dict(
    zip(
        FEATURES,
        [
            "Böttcher score",
            "Synthetic-accessibility score",
            "Bertz CT",
            "Stereocenters (count)",
            "Fraction sp³ carbon",
            "Rings (count)",
            "Rotatable bonds (count)",
            "H-bond acceptors (count)",
            "H-bond donors (count)",
            "Heavy atoms (count)",
            "TPSA (Å²)",
            "logP",
            "Molecular weight (Da)",
        ],
    )
)


def file_name(version: str) -> str:
    return NAME if version == "primary" else f"{NAME}_{version}"


def draw(descriptors: Path, out: Path, version: str) -> list[dict]:
    cohort = pd.read_csv(
        descriptors / "TableS4_descriptor_cohort.csv", dtype=str, keep_default_na=False
    )
    data = numeric(
        cohort.loc[cohort.network_version.eq(version)].reset_index(drop=True),
        FEATURES,
    )
    data["psi_bits"] = [float(v) if v else float("nan") for v in data["psi_bits"]]
    rows = [
        row
        for row in read_csv(descriptors / "TableS4_associations.csv")
        if row["network_version"] == version
        and row["cohort"] == "all"
        and row["target"] == "psi_bits"
    ]
    fig, axes = plt.subplots(4, 4, figsize=(14, 12))
    panels = []
    for ax, row in zip(axes.flat, rows):
        descriptor, n, rho = (
            row["descriptor"],
            int(row["n"]),
            float(row["spearman_rho"]),
        )
        valid = data[descriptor].notna() & data.psi_bits.notna()
        if valid.sum() != n:
            raise ValueError(f"{descriptor}: plotted compounds differ from Table S4")
        ax.scatter(
            data.loc[valid, descriptor],
            data.loc[valid, "psi_bits"],
            s=3,
            alpha=0.12,
            color=BLUE,
            rasterized=True,
        )
        ax.set(
            xlabel=LABELS[descriptor],
            ylabel="ψ (bits)",
            title=f"ρ = {rho:.3f}; n = {n:,}",
        )
        panels.append(
            {
                "descriptor": descriptor,
                "n": n,
                "rho": rho,
                "included_ids": data.loc[valid, "compound_id"].tolist(),
            }
        )
    for ax in list(axes.flat)[len(rows) :]:
        ax.axis("off")
    fig.tight_layout()
    name = file_name(version)
    for extension in ("png", "pdf"):
        metadata = {"Title": name}
        if extension == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        fig.savefig(
            out / f"{name}.{extension}", dpi=220, bbox_inches="tight", metadata=metadata
        )
    plt.close(fig)
    return panels


def main(results: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    descriptors = results / "descriptors"
    write_json(
        out / f"{NAME}.json",
        {version: draw(descriptors, out, version) for version in NETWORK_VERSIONS},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.results, args.out)
