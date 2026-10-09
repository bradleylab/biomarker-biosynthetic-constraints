"""Supplementary Figure 6: psi, lambda and sigma of the biomarkers by class.

Points are biomarkers with at least one path, black bars mark class medians, and
labels give the number with a path and the class size. Writes
SupplementaryFigure6.{png,pdf} (primary graph),
SupplementaryFigure6_with_conditional.{png,pdf} and SupplementaryFigure6.json.

  python analyses/09_figures/supplementary_figure6.py --results results --out results/figures
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
from biosynthesis.tables import read_csv, write_json  # noqa: E402

BLUE = "#3274A1"
NAME = "SupplementaryFigure6"
AXES = [
    ("psi_bits", "Path multiplicity, ψ (bits)"),
    ("lambda", "Normalized branch depth, λ"),
    ("sigma", "Fraction shared, σ"),
]


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


def biomarkers(projection: Path, version: str) -> pd.DataFrame:
    rows = [
        row
        for row in read_csv(projection / "TableS1_projection.csv")
        if row["network_version"] == version
    ]
    return pd.DataFrame(
        {
            "compound_id": [row["compound_id"] for row in rows],
            "class": [row["class"] for row in rows],
            "num_paths": [int(row["num_paths"]) for row in rows],
            **{
                field: [
                    float(row[field]) if row[field] else float("nan") for row in rows
                ]
                for field, _ in AXES
            },
        }
    )


def draw(projection: Path, out: Path, version: str) -> list[dict]:
    table = biomarkers(projection, version)
    classes = list(dict.fromkeys(table["class"]))
    values = []
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
    for ax, (field, label) in zip(axes, AXES):
        for i, group in enumerate(classes):
            source = table[table["class"] == group]
            selected = source.loc[source[field].notna(), field]
            if not source.loc[source[field].isna(), "num_paths"].eq(0).all():
                raise ValueError("A compound with a path lacks a measurement")
            ax.scatter(selected, [i] * len(selected), s=23, alpha=0.55, color=BLUE)
            median = float(selected.median()) if len(selected) else None
            if len(selected):
                ax.plot([median] * 2, [i - 0.22, i + 0.22], color="black", lw=2)
            values.append(
                {
                    "class": group,
                    "measurement": field,
                    "total": len(source),
                    "with_path": len(selected),
                    "median": median,
                    "values": dict(
                        zip(source.loc[selected.index, "compound_id"], selected)
                    ),
                }
            )
        ax.set_xlabel(label)
        ax.grid(axis="x", alpha=0.2)
    axes[0].set_yticks(
        range(len(classes)),
        [
            f"{c} ({table[table['class'] == c].psi_bits.notna().sum()}/{(table['class'] == c).sum()})"
            for c in classes
        ],
    )
    axes[0].invert_yaxis()
    fig.tight_layout()
    save(fig, out, file_name(version))
    return values


def main(results: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    projection = results / "projected_paths"
    write_json(
        out / f"{NAME}.json",
        {version: draw(projection, out, version) for version in NETWORK_VERSIONS},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.results, args.out)
