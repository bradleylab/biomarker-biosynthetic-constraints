"""Table S6: the assignments of fossil classes to precursor classes and of precursor classes to source categories.

Writes five files from data/geology/: the fossil-to-precursor links
(TableS6_analyte_precursor_map.csv), the precursor-to-source links
(TableS6_source_compatibility.csv), the source profile recorded for context
(TableS6_source_profiles.csv), the cited sources (TableS6_citations.csv) and
the 14 combinations of fossil class, precursor class and source category used
in the analysis (TableS6_compatibility_links.csv). A combination is strict only
when both of its links are strict.

Usage:
  python analyses/05_source_compatibility/source_tables.py --out results/sources
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.sources import BROADER, STRICT  # noqa: E402
from biosynthesis.tables import read_csv, write_csv  # noqa: E402

# Columns of data/geology/precursor_to_source.csv not repeated in Table S6
# (included_in_analysis is true in every row).
COMPATIBILITY_COLUMNS = (
    "compatibility_id",
    "precursor_id",
    "source_class_id",
    "disposition",
    "assignment_criteria",
    "evidence_basis",
    "modern_compatibility_only",
    "overlap_group",
    "citation_ids",
    "evidence_locator",
    "evidence_summary",
    "conflict_status",
    "interpretation_note",
)
LINK_COLUMNS = (
    "analyte_id",
    "precursor_id",
    "source_class_id",
    "assignment_criteria",
    "mapping_record_id",
    "compatibility_record_id",
    "mapping_citation_ids",
    "compatibility_citation_ids",
)


def compatibility_links(mappings, compatibilities) -> list[dict[str, str]]:
    """Join the links used in the analysis through their precursor class."""
    used_mappings = [row for row in mappings if row["included_in_analysis"] == "true"]
    by_precursor: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in compatibilities:
        if row["included_in_analysis"] == "true":
            by_precursor[row["precursor_id"]].append(row)
    links = []
    for mapping in used_mappings:
        for compatibility in by_precursor[mapping["precursor_id"]]:
            strict = (
                mapping["assignment_criteria"]
                == compatibility["assignment_criteria"]
                == STRICT
            )
            links.append(
                {
                    "analyte_id": mapping["analyte_id"],
                    "precursor_id": mapping["precursor_id"],
                    "source_class_id": compatibility["source_class_id"],
                    "assignment_criteria": STRICT if strict else BROADER,
                    "mapping_record_id": mapping["mapping_id"],
                    "compatibility_record_id": compatibility["compatibility_id"],
                    "mapping_citation_ids": mapping["citation_ids"],
                    "compatibility_citation_ids": compatibility["citation_ids"],
                }
            )
    return sorted(
        links,
        key=lambda row: (
            row["analyte_id"],
            row["source_class_id"],
            row["precursor_id"],
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder.")
    parser.add_argument(
        "--data", type=Path, default=REPO / "data", help="The repository data folder."
    )
    args = parser.parse_args()

    folder = args.data / "geology"
    tables = {}
    for name, source in (
        ("TableS6_analyte_precursor_map.csv", "fossil_to_precursor.csv"),
        ("TableS6_source_compatibility.csv", "precursor_to_source.csv"),
        ("TableS6_source_profiles.csv", "source_profiles.csv"),
        ("TableS6_citations.csv", "citations.csv"),
    ):
        tables[name] = read_csv(folder / source)
    columns = {name: tuple(rows[0]) for name, rows in tables.items()}
    columns["TableS6_source_compatibility.csv"] = COMPATIBILITY_COLUMNS
    for name, rows in tables.items():
        write_csv(args.out / name, rows, columns[name])
    links = compatibility_links(
        tables["TableS6_analyte_precursor_map.csv"],
        tables["TableS6_source_compatibility.csv"],
    )
    write_csv(args.out / "TableS6_compatibility_links.csv", links, LINK_COLUMNS)
    used = [
        row
        for row in tables["TableS6_analyte_precursor_map.csv"]
        if row["included_in_analysis"] == "true"
    ]
    print(
        f"{len(used)} fossil-to-precursor links used ({len({r['analyte_id'] for r in used})} fossil classes, "
        f"{len({r['precursor_id'] for r in used})} precursor classes); "
        f"{len(tables['TableS6_source_compatibility.csv'])} precursor-to-source links; "
        f"{len(links)} combinations."
    )


if __name__ == "__main__":
    main()
