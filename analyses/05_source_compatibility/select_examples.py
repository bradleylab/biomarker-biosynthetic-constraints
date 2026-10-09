"""Select the two geological samples shown in Figure 4.

At a detection threshold of zero with low signals included, each sample has
one set of detected fossil classes. Two detection patterns are considered:
(1) a hopanoid class (aromatic hopanoids or hopanes) and cyclosteranes, with
no conventional sterane (cholestanes, ergostanes, stigmastanes) detected; and
(2) a hopanoid class and a conventional sterane, with no cyclosteranes
detected. For each pattern the most frequent set of detected classes is found
(ties would be broken by the sorted class identifiers), and the example is the
first sample of that set in sorted order of sample label. Writes
figure4_examples.csv and prints the counts.

Usage:
  python analyses/05_source_compatibility/select_examples.py --results results/sources
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.sensitivity import zero_threshold  # noqa: E402
from biosynthesis.sources import STRICT  # noqa: E402
from biosynthesis.tables import read_jsonl, write_csv  # noqa: E402

HOPANOIDS = frozenset({"si1_aromatic_hopanoids", "si2_hopanes"})
CONVENTIONAL_STERANES = frozenset(
    {"si2_cholestanes", "si2_ergostanes", "si2_stigmastanes"}
)
CYCLOSTERANES = frozenset({"si2_cyclosteranes"})
PATTERNS = (
    (
        "hopanoids and cyclosteranes, no conventional steranes",
        HOPANOIDS,
        CYCLOSTERANES,
        CONVENTIONAL_STERANES,
    ),
    (
        "hopanoids and conventional steranes, no cyclosteranes",
        HOPANOIDS,
        CONVENTIONAL_STERANES,
        CYCLOSTERANES,
    ),
)


def detected_at_zero(scenarios) -> dict[str, tuple[str, ...]]:
    """The set of detected classes of each sample at a threshold of zero with low signals included."""
    detected = {}
    for row in scenarios:
        if (
            row["assignment_criteria"] == STRICT
            and row["low_signal_mode"] == "include"
            and zero_threshold(row)
        ):
            detected[row["sample_id"]] = tuple(row["detected_analyte_ids"])
    return detected


def sample_labels(samples) -> dict[str, str]:
    """The workbook label of each sample (from its first workbook row)."""
    return {
        row["sample_id"]: row["source_records"][0]["source_sample_id_text"].strip()
        for row in samples
    }


def select(
    detected: dict[str, tuple[str, ...]], labels: dict[str, str]
) -> list[dict[str, object]]:
    """The most frequent set of detected classes of each pattern and its first sample by label."""
    rows = []
    for name, either, required, absent in PATTERNS:
        counts = Counter(
            classes
            for classes in detected.values()
            if set(classes) & either
            and set(classes) & required
            and not set(classes) & absent
        )
        ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        (classes, count), runner_up = ranked[0], ranked[1][1] if len(ranked) > 1 else 0
        members = sorted(
            (labels[sample_id], sample_id)
            for sample_id, found in detected.items()
            if found == classes
        )
        rows.append(
            {
                "pattern": name,
                "detected_analyte_ids": ";".join(classes),
                "samples_with_this_set": count,
                "next_most_frequent_set": runner_up,
                "sample_label": members[0][0],
                "sample_id": members[0][1],
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--results",
        type=Path,
        required=True,
        help="Folder written by compatible_combinations.py.",
    )
    args = parser.parse_args()
    rows = select(
        detected_at_zero(read_jsonl(args.results / "sample_scenarios.jsonl")),
        sample_labels(read_jsonl(args.results / "samples.jsonl")),
    )
    write_csv(args.results / "figure4_examples.csv", rows, tuple(rows[0]))
    for row in rows:
        print(
            f"Samples with {row['pattern']}: the most frequent set of detected classes "
            f"({row['detected_analyte_ids'].replace(';', ', ')}) occurs in {row['samples_with_this_set']} samples "
            f"(next most frequent: {row['next_most_frequent_set']}); first by label: {row['sample_label']}."
        )


if __name__ == "__main__":
    main()
