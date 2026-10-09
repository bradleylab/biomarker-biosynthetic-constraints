"""Draw Supplementary Figure 4: pair matrices of enzymes and pathways.

Each cell is a pair of biomarkers, classified from the set sizes in Table S3: no
record for one or both, no shared enzyme or pathway, overlap that is not identical,
or identical sets. The top row compares enzymes and the bottom row pathways, under
strict and broader criteria. A dot marks a direct-precursor pair in the row of the
precursor and the column of its product; a pair in which each is a precursor of the
other has a dot in both cells. Writes SupplementaryFigure4.pdf, .png and
SupplementaryFigure4_values.json (the plotted values).

Usage:
  python analyses/09_figures/supplementary_figure4.py --results results --out results/figures
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from biosynthesis import enzymes
from biosynthesis.tables import read_csv, write_json

PAIRS = "pairs"
PANEL_SIZE = 86
COLORS = ("#eeeeee", "#808080", "#3274A1", "#CC8B25")
NAMES = (
    "No record for one or both",
    "No shared enzyme or pathway",
    "Overlap, not identical",
    "Identical sets",
)
KEYS = (
    "no_record_for_one_or_both",
    "no_shared_enzyme_or_pathway",
    "overlap_not_identical",
    "identical_sets",
)
TITLES = {"strict": "Strict criteria", "broader": "Broader criteria"}
UNITS = {
    "enzyme": ("enzyme records", "{c}_{k}_enzyme_count"),
    "pathway": ("pathways", "{c}_{k}_pathway_count"),
}
PANELS = (
    ("strict", "enzyme"),
    ("broader", "enzyme"),
    ("strict", "pathway"),
    ("broader", "pathway"),
)
PRECURSOR_LABEL = "Direct precursor (row biomarker to column biomarker)"
AXIS_LABEL = "Biomarker number (Table S1)"


def classify(row: dict, criteria: str, unit: str) -> int:
    """0 no record for one or both, 1 nothing shared, 2 overlap, 3 identical sets."""
    a, b, shared = (
        int(row[UNITS[unit][1].format(c=criteria, k=k)])
        for k in ("left", "right", "shared")
    )
    if not a or not b:
        return 0
    if not shared:
        return 1
    return 3 if a == b == shared else 2


def matrices(rows: list[dict]) -> dict[tuple[str, str], np.ndarray]:
    found = {}
    for criteria, unit in PANELS:
        matrix = np.zeros((PANEL_SIZE, PANEL_SIZE), dtype=int)
        for row in rows:
            i, j = (
                int(row["left_biomarker_number"]) - 1,
                int(row["right_biomarker_number"]) - 1,
            )
            matrix[i, j] = matrix[j, i] = classify(row, criteria, unit)
        found[(criteria, unit)] = matrix
    return found


def precursor_cells(rows: list[dict], criteria: str) -> list[tuple[int, int]]:
    """(row, column) cells: the precursor's row and its product's column."""
    cells = []
    for row in rows:
        mark = row[f"{criteria}_precursor_relation"]
        left, right = (
            int(row["left_biomarker_number"]) - 1,
            int(row["right_biomarker_number"]) - 1,
        )
        if mark in (enzymes.LEFT_PRECURSOR, enzymes.EACH_PRECURSOR):
            cells.append((left, right))
        if mark in (enzymes.RIGHT_PRECURSOR, enzymes.EACH_PRECURSOR):
            cells.append((right, left))
    return cells


def counts(matrix: np.ndarray) -> dict[str, int]:
    upper = matrix[np.triu_indices(PANEL_SIZE, k=1)]
    return {key: int((upper == i).sum()) for i, key in enumerate(KEYS)}


def draw(rows: list[dict], out: Path) -> dict:
    found = matrices(rows)
    cmap = ListedColormap(COLORS)
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    panels = {}
    for ax, (criteria, unit) in zip(axes.flat, PANELS):
        ax.imshow(
            found[(criteria, unit)],
            cmap=cmap,
            norm=BoundaryNorm(np.arange(-0.5, 4.5), 4),
            interpolation="nearest",
        )
        cells = precursor_cells(rows, criteria)
        ax.scatter(
            [c for _, c in cells],
            [r for r, _ in cells],
            s=5,
            c="black",
            edgecolors="white",
            linewidths=0.3,
            zorder=3,
        )
        ticks = list(range(0, PANEL_SIZE, 10))
        ax.set_xticks(ticks, [t + 1 for t in ticks])
        ax.set_yticks(ticks, [t + 1 for t in ticks])
        ax.set(
            xlabel=AXIS_LABEL,
            ylabel=AXIS_LABEL,
            title=f"{TITLES[criteria]}: {UNITS[unit][0]}",
        )
        panels[f"{criteria}_{unit}"] = {
            **counts(found[(criteria, unit)]),
            "direct_precursor_cells": len(cells),
        }
    handles = [Patch(facecolor=cmap(i), label=name) for i, name in enumerate(NAMES)]
    handles.append(
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markerfacecolor="black",
            markeredgecolor="white",
            markeredgewidth=0.3,
            markersize=4,
            label=PRECURSOR_LABEL,
        )
    )
    fig.legend(handles=handles, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.06, 1, 0.98))
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        metadata = {"CreationDate": None, "ModDate": None} if ext == "pdf" else None
        fig.savefig(
            out / f"SupplementaryFigure4.{ext}",
            dpi=300,
            bbox_inches="tight",
            metadata=metadata,
        )
    plt.close(fig)
    return {
        "panels": panels,
        "matrices": {f"{c}_{u}": matrix.tolist() for (c, u), matrix in found.items()},
    }


def main(args: argparse.Namespace) -> None:
    rows = read_csv(args.results / PAIRS / "TableS3_pair_comparisons.csv")
    values = draw(rows, args.out)
    write_json(args.out / "SupplementaryFigure4_values.json", values)
    print(json.dumps(values["panels"], indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True, help="results folder")
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    main(parser.parse_args())
