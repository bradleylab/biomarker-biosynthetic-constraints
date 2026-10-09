"""Supplementary Figure 3: compounds without an accepted producer reached by tracing.

Counts the endpoints of every distinct traced query in the primary graph, separately
for compounds with no producer in the network and compounds whose producing
reactions were excluded by the criteria. A compound shared by several queries is
counted once per query.

    python analyses/09_figures/supplementary_figure3.py --results results --out results/figures
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "supplementary-figure-3"
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.tables import read_jsonl, write_json  # noqa: E402

CRITERIA = ("strict", "broader")
CRITERIA_LABELS = ("Strict criteria", "Broader criteria")
CATEGORIES = (
    ("endpoints_no_producer", "No producer in network"),
    ("endpoints_producer_excluded", "Producer excluded by criteria"),
)
COLOR = "#3274A1"


def endpoint_counts(records: list[dict]) -> dict:
    """Endpoints summed over the distinct primary-graph queries of each set of criteria."""
    counts = {}
    for criteria in CRITERIA:
        chosen = [
            row
            for row in records
            if row["network_version"] == "primary"
            and row["reaction_criteria"] == criteria
        ]
        if not all(row["complete"] for row in chosen):
            raise ValueError("Tracing did not finish for a query")
        counts[criteria] = {
            "queries": len(chosen),
            **{field: sum(row[field] for row in chosen) for field, _ in CATEGORIES},
        }
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    counts = endpoint_counts(
        read_jsonl(args.results / "tracing" / "query_results.jsonl")
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, criteria, title in zip(axes, CRITERIA, CRITERIA_LABELS):
        values = [counts[criteria][field] for field, _ in CATEGORIES]
        ax.barh([label for _, label in CATEGORIES], values, color=COLOR)
        for i, value in enumerate(values):
            ax.text(value, i, f" {value:,}", va="center")
        ax.set(xlabel="Endpoints, summed over queries", title=title)
        ax.margins(x=0.2)
    fig.tight_layout()
    args.out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        metadata = (
            {"CreationDate": None, "ModDate": None}
            if ext == "pdf"
            else ({"Date": None} if ext == "svg" else None)
        )
        fig.savefig(
            args.out / f"SupplementaryFigure3.{ext}",
            dpi=300,
            bbox_inches="tight",
            metadata=metadata,
        )
    plt.close(fig)
    write_json(
        args.out / "SupplementaryFigure3_values.json",
        {
            "distinct_queries": sum(value["queries"] for value in counts.values()),
            "endpoint_counts": counts,
        },
    )
    print(f"Supplementary Figure 3: {counts}", flush=True)


if __name__ == "__main__":
    main()
