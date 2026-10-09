"""Compare every pair of the 86 biomarkers by enzymes and by pathways (Table S3).

For each set of reaction criteria (strict, broader), the producing reactions of a
biomarker are the reaction records that directly produce one of its MetaCyc
compounds and that the primary network keeps; biomarkers with several candidate
identities are pooled over them. The enzyme set of a biomarker holds the ENZYME of
every enzrxns.dat record whose REACTION is a producing reaction. The pathway set
holds the pathways that list a producing reaction and whose recorded direction for
it agrees with production of the biomarker.

Writes, to the folder given by --out:

- direct_production_records.jsonl: the producing records of every biomarker identity.
- pathway_anchors.jsonl: the pathways of those records, with their recorded direction.
- TableS3_pair_comparisons.csv: Table S3.
- biomarker_enzyme_sets.csv: producing reactions, enzymes and pathway counts of
  each biomarker, read by the figure scripts.

Usage:
  python analyses/04_enzyme_comparisons/compare_pairs.py --metacyc DIR \\
    --networks results/networks --out results/pairs
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

from biosynthesis import enzymes, panel
from biosynthesis.tables import write_csv, write_jsonl

TABLE_NAME = "TableS3_pair_comparisons.csv"
SETS_NAME = "biomarker_enzyme_sets.csv"
PRODUCTION_NAME = "direct_production_records.jsonl"
ANCHORS_NAME = "pathway_anchors.jsonl"
REACTION_RECORDS = "reaction_records.jsonl"
PRIMARY_EDGES = {
    "strict": Path("strict_primary/edges.jsonl"),
    "broader": Path("broader_primary/edges.jsonl"),
}
IDENTITY_COLUMNS = [
    f"{side}_{key}"
    for side in ("left", "right")
    for key in ("biomarker_number", "compound_id", "name", "biomarker_class")
]
ENZYME_COLUMNS = (
    "left_enzyme_count",
    "right_enzyme_count",
    "shared_enzyme_count",
    "shared_enzymes",
    "enzyme_left_retained_fraction",
    "enzyme_right_retained_fraction",
    "enzyme_left_information_bits",
    "enzyme_right_information_bits",
    "enzyme_jaccard",
    "enzyme_relation",
)
PATHWAY_COLUMNS = (
    "pathway_relation",
    "left_pathway_count",
    "right_pathway_count",
    "shared_pathway_count",
    "shared_pathways",
    "pathway_left_retained_fraction",
    "pathway_right_retained_fraction",
    "pathway_left_information_bits",
    "pathway_right_information_bits",
    "pathway_jaccard",
)
SETS_COLUMNS = (
    "reaction_criteria",
    "biomarker_number",
    "compound_id",
    "name",
    "literature_defined",
    "producing_reactions",
    "reactions_without_enzyme",
    "enzyme_count",
    "enzymes",
    "pathway_count",
)


def table_columns() -> list[str]:
    """Column order of Table S3."""
    columns = list(IDENTITY_COLUMNS)
    for criteria in enzymes.CRITERIA:
        columns += [f"{criteria}_{name}" for name in ENZYME_COLUMNS]
        columns += [f"{criteria}_{name}" for name in PATHWAY_COLUMNS]
        columns.append(f"{criteria}_precursor_relation")
    return columns


def comparison_cells(result: dict, criteria: str, unit: str) -> dict:
    """Cells of one set comparison, for enzymes or pathways."""
    if unit == "enzyme":
        values = {
            "left_enzyme_count": result["left_count"],
            "right_enzyme_count": result["right_count"],
            "shared_enzyme_count": result["shared_count"],
            "shared_enzymes": result["shared"],
        }
    else:
        values = {
            "left_pathway_count": result["left_count"],
            "right_pathway_count": result["right_count"],
            "shared_pathway_count": result["shared_count"],
            "shared_pathways": result["shared"],
        }
    values |= {
        f"{unit}_left_retained_fraction": result["left_fraction"],
        f"{unit}_right_retained_fraction": result["right_fraction"],
        f"{unit}_left_information_bits": result["left_bits"],
        f"{unit}_right_information_bits": result["right_bits"],
        f"{unit}_jaccard": result["jaccard"],
        f"{unit}_relation": result["relation"],
    }
    return {f"{criteria}_{key}": value for key, value in values.items()}


def pair_rows(biomarkers: list[dict], sets: dict, reactants: dict) -> list[dict]:
    """One row per unordered pair of biomarkers, lower biomarker number on the left."""
    rows = []
    for left, right in itertools.combinations(biomarkers, 2):
        row = {}
        for side, biomarker in (("left", left), ("right", right)):
            row |= {
                f"{side}_biomarker_number": biomarker["biomarker_number"],
                f"{side}_compound_id": biomarker["compound_id"],
                f"{side}_name": biomarker["name"],
                f"{side}_biomarker_class": biomarker["class"],
            }
        pair = (left["biomarker_number"], right["biomarker_number"])
        for criteria in enzymes.CRITERIA:
            for unit, key in (("enzyme", "enzymes"), ("pathway", "pathways")):
                found = sets[criteria][key]
                result = enzymes.compare(found[pair[0]], found[pair[1]])
                row |= comparison_cells(result, criteria, unit)
            row[f"{criteria}_precursor_relation"] = enzymes.precursor_relation(
                (pair[0], left["compound_id"]),
                (pair[1], right["compound_id"]),
                reactants[criteria],
            )
        rows.append(row)
    return rows


def set_rows(
    biomarkers: list[dict], sets: dict, index: enzymes.EnzymeIndex
) -> list[dict]:
    """Producing reactions, enzymes and pathway count of each biomarker."""
    rows = []
    for criteria in enzymes.CRITERIA:
        for biomarker in biomarkers:
            number = biomarker["biomarker_number"]
            reactions = sets[criteria]["reactions"][number]
            rows.append(
                {
                    "reaction_criteria": criteria,
                    "biomarker_number": number,
                    "compound_id": biomarker["compound_id"],
                    "name": biomarker["name"],
                    "literature_defined": enzymes.is_literature(
                        biomarker["compound_id"]
                    ),
                    "producing_reactions": sorted(reactions),
                    "reactions_without_enzyme": sorted(
                        r for r in reactions if not index.enzymes({r})
                    ),
                    "enzyme_count": len(sets[criteria]["enzymes"][number]),
                    "enzymes": sorted(sets[criteria]["enzymes"][number]),
                    "pathway_count": len(sets[criteria]["pathways"][number]),
                }
            )
    return rows


def main(args: argparse.Namespace) -> None:
    identities = panel.load_identities(args.identities)
    biomarkers = panel.load_biomarkers(args.identities)
    numbers = [row["biomarker_number"] for row in biomarkers]
    compounds = {str(c) for row in identities for c in row["metacyc_compounds"]}
    index = enzymes.records_by_compound(args.networks / REACTION_RECORDS, compounds)
    production = enzymes.production_records(identities, index)
    records = {
        str(row["source_occurrence_id"]): row for rows in index.values() for row in rows
    }
    anchors = enzymes.pathway_anchors(
        production, records, enzymes.pathway_ids(args.metacyc / "pathways.dat")
    )
    write_jsonl(args.out / PRODUCTION_NAME, production)
    write_jsonl(args.out / ANCHORS_NAME, anchors)

    enzyme_index = enzymes.enzyme_index(args.metacyc)
    sets, reactants = {}, {}
    for criteria in enzymes.CRITERIA:
        kept = enzymes.kept_occurrences(args.networks / PRIMARY_EDGES[criteria])
        sets[criteria] = enzymes.biomarker_sets(
            criteria, production, anchors, kept, numbers
        )
        sets[criteria]["enzymes"] = {
            n: enzyme_index.enzymes(r) for n, r in sets[criteria]["reactions"].items()
        }
        reactants[criteria] = enzymes.precursors_by_biomarker(
            criteria, production, kept
        )
    count = write_csv(
        args.out / TABLE_NAME, pair_rows(biomarkers, sets, reactants), table_columns()
    )
    write_csv(
        args.out / SETS_NAME, set_rows(biomarkers, sets, enzyme_index), SETS_COLUMNS
    )
    print(
        f"{len(production)} producing records and {len(anchors)} pathway placements; "
        f"{count} pairs written to {args.out / TABLE_NAME}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--metacyc", type=Path, required=True, help="MetaCyc 29.5 data folder"
    )
    parser.add_argument(
        "--networks",
        type=Path,
        required=True,
        help="network folder with reaction_records.jsonl and <criteria>_primary/edges.jsonl",
    )
    parser.add_argument(
        "--identities",
        type=Path,
        default=panel.IDENTITIES,
        help="biomarker identities (default: %(default)s)",
    )
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    main(parser.parse_args())
