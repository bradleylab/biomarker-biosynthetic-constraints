"""Supplementary Figure 2 and its detail pages: search status of every biomarker identity.

Reads Table S1 tracing (primary graph rows) and writes SupplementaryFigure2.{png,pdf,svg},
SupplementaryFigure2_DETAIL.pdf and the plotted values.

    python analyses/09_figures/supplementary_figure2.py --results results --out results/figures
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "supplementary-figure-2"
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.tables import read_csv, write_json  # noqa: E402

CRITERIA = ("strict", "broader")
CRITERIA_LABELS = ("Strict criteria", "Broader criteria")
CODES = ("not_in_network", "exhausted_no_producer", "exhausted_producer", "open")
COLORS = ("#dddddd", "#b7c4cd", "#3274A1", "#CC8B25")
SHORT = ("Not in network", "Exhausted: no producer", "Exhausted: producer", "Open")
CANDIDATE_NOTE = "* Identity not accepted"
IDENTITY_LABELS = {
    "accepted": "Accepted",
    "candidate_not_accepted": "Candidate, not accepted",
    "no_candidate": "No candidate",
}


def display_status(row: dict) -> str:
    """Search status shown for one Table S1 row."""
    if row["status"] == "not_in_network":
        return "not_in_network"
    if row["status"] == "open":
        return "open"
    if row["status"] != "exact":
        raise ValueError("Unrecognized search status: " + row["status"])
    return (
        "exhausted_no_producer"
        if int(row["producing_reactions"]) == 0
        else "exhausted_producer"
    )


def save_figure(fig, out: Path, name: str) -> None:
    """Write PNG, PDF and SVG versions of a figure without creation dates."""
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        metadata = (
            {"CreationDate": None, "ModDate": None}
            if ext == "pdf"
            else ({"Date": None} if ext == "svg" else None)
        )
        fig.savefig(
            out / f"{name}.{ext}", dpi=300, bbox_inches="tight", metadata=metadata
        )
    plt.close(fig)


def overview(states: list[dict], out: Path) -> None:
    """One colored column per set of criteria, one row per identity and node."""
    fig, axes = plt.subplots(1, 2, figsize=(17, 25))
    for ax, criteria, title in zip(axes, CRITERIA, CRITERIA_LABELS):
        selected = [row for row in states if row["reaction_criteria"] == criteria]
        ax.imshow(
            [[CODES.index(row["display"])] for row in selected],
            aspect="auto",
            cmap=ListedColormap(COLORS),
            vmin=0,
            vmax=3,
        )
        labels = [
            f"{row['biomarker_number']}: {row['name']}"
            + (" *" if row["identity"] != "accepted" else "")
            for row in selected
        ]
        ax.set_yticks(range(len(labels)), labels, fontsize=7)
        ax.set_xticks([0], ["Search"])
        ax.set_title(title)
        for i, row in enumerate(selected):
            code = CODES.index(row["display"])
            ax.text(
                0,
                i,
                SHORT[code],
                ha="center",
                va="center",
                fontsize=6,
                color="white" if code == 2 else "black",
            )
    fig.text(0.5, 0.005, CANDIDATE_NOTE, ha="center", fontsize=11)
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    save_figure(fig, out, "SupplementaryFigure2")


def detail(states: list[dict], out: Path) -> None:
    """Table pages listing every identity and node with its search status."""
    with PdfPages(
        out / "SupplementaryFigure2_DETAIL.pdf",
        metadata={
            "Title": "Supplementary Figure 2 detail",
            "CreationDate": None,
            "ModDate": None,
        },
    ) as pdf:
        for criteria, title in zip(CRITERIA, CRITERIA_LABELS):
            selected = [row for row in states if row["reaction_criteria"] == criteria]
            for start in range(0, len(selected), 18):
                batch = selected[start : start + 18]
                fig, ax = plt.subplots(figsize=(8.27, 11.69))
                ax.axis("off")
                fig.subplots_adjust(left=0.06, right=0.96, top=0.87, bottom=0.12)
                cells = []
                for row in batch:
                    node = row["node"] or "not in network"
                    cells.append(
                        [
                            textwrap.fill(
                                f"{row['biomarker_number']}: {row['name']}", 37
                            )
                            + f"\nIdentity {row['identity_id']}; node {node}",
                            textwrap.fill(IDENTITY_LABELS[row["identity"]], 20),
                            textwrap.fill(SHORT[CODES.index(row["display"])], 17),
                        ]
                    )
                table = ax.table(
                    cellText=cells,
                    colLabels=("Biomarker and identity", "Identity", "Search"),
                    colWidths=(0.55, 0.24, 0.21),
                    cellLoc="left",
                    bbox=(0, 0, 1, 1),
                )
                table.auto_set_font_size(False)
                table.set_fontsize(7.5)
                for (i, j), cell in table.get_celld().items():
                    cell.set_edgecolor("#b9c1c9")
                    cell.set_linewidth(0.3)
                    if i == 0:
                        cell.set_facecolor("#e9eef3")
                    elif j == 2:
                        code = CODES.index(batch[i - 1]["display"])
                        cell.set_facecolor(COLORS[code])
                        if code == 2:
                            cell.get_text().set_color("white")
                fig.text(
                    0.06,
                    0.95,
                    "Supplementary Figure 2 — Biomarker identities",
                    fontsize=13,
                    weight="bold",
                )
                fig.text(0.06, 0.92, title, fontsize=12)
                fig.text(
                    0.06,
                    0.89,
                    f"Rows {start + 1}–{start + len(batch)} of {len(selected)}",
                    fontsize=9,
                )
                fig.text(
                    0.06,
                    0.09,
                    "Primary graph; rows in Table S1 order within each set of criteria.\n"
                    "Identity and node identifiers as in Table S1.",
                    fontsize=8,
                    va="top",
                )
                pdf.savefig(fig)
                plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = read_csv(args.results / "tracing" / "TableS1_tracing.csv")
    states = [
        {**row, "display": display_status(row)}
        for row in rows
        if row["network_version"] == "primary"
    ]
    overview(states, args.out)
    detail(states, args.out)
    write_json(
        args.out / "SupplementaryFigure2_values.json",
        {
            "rows_by_criteria": {
                criteria: {
                    code: Counter(
                        row["display"]
                        for row in states
                        if row["reaction_criteria"] == criteria
                    )[code]
                    for code in CODES
                }
                for criteria in CRITERIA
            },
            "identity_rows": len(states),
            "candidate_rows": sum(row["identity"] != "accepted" for row in states),
            "rows": [
                {
                    key: row[key]
                    for key in (
                        "reaction_criteria",
                        "biomarker_number",
                        "identity_id",
                        "node",
                        "display",
                    )
                }
                for row in states
            ],
        },
    )
    print(f"Supplementary Figure 2: {len(states)} identity rows", flush=True)


if __name__ == "__main__":
    main()
