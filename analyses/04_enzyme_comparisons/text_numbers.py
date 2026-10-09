"""Print the numbers of the enzyme and pathway comparisons quoted in the paper.

Reads Table S3 and the biomarker sets written by compare_pairs.py, and the MetaCyc
flat files for the names and organisms of the enzymes of the worked examples.

Usage:
  python analyses/04_enzyme_comparisons/text_numbers.py --results results --metacyc DIR
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter
from pathlib import Path

from biosynthesis import enzymes
from biosynthesis.tables import read_csv

PAIRS = "pairs"
CRITERIA = ("strict", "broader")
SHARING = ("overlap", "identical")
PRECURSOR_PAIRS = (
    enzymes.LEFT_PRECURSOR,
    enzymes.RIGHT_PRECURSOR,
    enzymes.EACH_PRECURSOR,
)
# Compounds of the worked examples: beta-amyrin and alpha-amyrin; diploptene,
# tetrahymanol and diplopterol.
BETA_AMYRIN, ALPHA_AMYRIN = "CPD-6948", "CPD-8250"
DIPLOPTENE, TETRAHYMANOL, DIPLOPTEROL = "HOP-2229-ENE", "CPD-11533", "HOPAN-22-OL"
# Taxa of classes.dat used to name the organism groups of the example enzymes
# (classes.dat names TAX-1521262 "Polypodiidae <ferns>").
BACTERIA, FERNS, CILIATES = "TAX-2", "TAX-1521262", "TAX-5878"
GROUPS = {
    BACTERIA: "bacteria",
    FERNS: "ferns",
    CILIATES: "ciliates",
    "TAX-2759": "eukaryotes",
}


def pair_row(rows: list[dict], first: str, second: str) -> dict:
    """The Table S3 row of two biomarkers given by compound identifier."""
    found = [
        row
        for row in rows
        if {row["left_compound_id"], row["right_compound_id"]} == {first, second}
    ]
    if len(found) != 1:
        raise ValueError(f"expected one row for {first} and {second}")
    return found[0]


def side(row: dict, compound: str) -> str:
    return "left" if row["left_compound_id"] == compound else "right"


def criteria_counts(rows: list[dict], criteria: str) -> dict:
    """Pair counts of one set of criteria."""
    enzyme = Counter(row[f"{criteria}_enzyme_relation"] for row in rows)
    pathway = Counter(row[f"{criteria}_pathway_relation"] for row in rows)
    precursor = Counter(row[f"{criteria}_precursor_relation"] for row in rows)
    sharing = [row for row in rows if row[f"{criteria}_enzyme_relation"] in SHARING]
    bits = [
        float(row[f"{criteria}_enzyme_{s}_information_bits"])
        for row in sharing
        for s in ("left", "right")
    ]
    marked = [
        row for row in rows if row[f"{criteria}_precursor_relation"] in PRECURSOR_PAIRS
    ]
    both = enzyme["disjoint"] + enzyme["overlap"] + enzyme["identical"]
    return {
        "both_with_enzymes": both,
        "without_enzymes_one_or_both": len(rows) - both,
        "sharing_enzyme": enzyme["overlap"] + enzyme["identical"],
        "identical_enzymes": enzyme["identical"],
        "disjoint_enzymes": enzyme["disjoint"],
        "left_empty": enzyme["left_empty"],
        "right_empty": enzyme["right_empty"],
        "both_empty": enzyme["both_empty"],
        "median_bits": statistics.median(bits),
        "median_directions": len(bits),
        "overlap_with_positive_reduction": sum(
            any(
                float(row[f"{criteria}_enzyme_{s}_information_bits"]) > 0
                for s in ("left", "right")
            )
            for row in sharing
        ),
        "sharing_pathway": pathway["overlap"] + pathway["identical"],
        "identical_pathways": pathway["identical"],
        "precursor_pairs": len(marked),
        "precursor_pairs_sharing_enzyme": sum(
            row[f"{criteria}_enzyme_relation"] in SHARING for row in marked
        ),
        "precursor_not_tested": precursor[enzymes.NOT_TESTED],
        "precursor_no": precursor[enzymes.NOT_PRECURSOR],
        "identical_enzyme_pairs": [
            (row["left_name"], row["right_name"])
            for row in rows
            if row[f"{criteria}_enzyme_relation"] == "identical"
        ],
    }


def ancestors(taxon: str, classes: dict) -> set[str]:
    """Every TYPES parent of a taxon in classes.dat, at any depth."""
    seen, stack = set(), [taxon]
    while stack:
        current = stack.pop()
        if current not in seen and current in classes:
            seen.add(current)
            stack.extend(classes[current].get("TYPES", []))
    return seen


def protein_name(enzyme: str, proteins: dict, genes: dict) -> str:
    """Protein name of an enzyme, or its gene names when the protein has none."""
    names = proteins[enzyme].get("COMMON-NAME") or proteins[enzyme].get("SYNONYMS")
    if names:
        return enzymes.display_name(names[0])
    gene_names = [
        enzymes.display_name(genes[g]["COMMON-NAME"][0])
        for g in proteins[enzyme].get("GENE", [])
        if genes.get(g, {}).get("COMMON-NAME")
    ]
    return f"no protein name; gene {', '.join(gene_names)}" if gene_names else "no name"


def enzyme_details(
    found: set[str], index: enzymes.EnzymeIndex, genes: dict, reactions: set
) -> list[dict]:
    """Name, organism, organism group and producing reactions of each enzyme."""
    details = []
    for enzyme in sorted(found):
        taxa, _ = index.taxa(enzyme)
        lineage = set().union(*(ancestors(t, index.classes) for t in taxa))
        details.append(
            {
                "enzyme": enzyme,
                "name": protein_name(enzyme, index.proteins, genes),
                "organisms": [
                    enzymes.display_name(index.classes[t]["COMMON-NAME"][0])
                    for t in sorted(taxa)
                ],
                "groups": [
                    label for taxon, label in GROUPS.items() if taxon in lineage
                ],
                "reactions": sorted(
                    r for r in reactions if enzyme in index.enzymes({r})
                ),
            }
        )
    return details


def numbers(results: Path, metacyc: Path) -> dict:
    """Every number of the enzyme and pathway comparisons quoted in the text."""
    rows = read_csv(results / PAIRS / "TableS3_pair_comparisons.csv")
    sets = {
        (row["reaction_criteria"], row["compound_id"]): row
        for row in read_csv(results / PAIRS / "biomarker_enzyme_sets.csv")
    }

    def enzyme_set(criteria: str, compound: str) -> set[str]:
        return set(json.loads(sets[(criteria, compound)]["enzymes"]))

    def reaction_set(criteria: str, compound: str) -> set[str]:
        return set(json.loads(sets[(criteria, compound)]["producing_reactions"]))

    found = {
        "pairs": len(rows),
        "biomarkers": len(
            {row["left_biomarker_number"] for row in rows}
            | {row["right_biomarker_number"] for row in rows}
        ),
    }
    for criteria in CRITERIA:
        found |= {
            f"{criteria}_{k}": v for k, v in criteria_counts(rows, criteria).items()
        }

    amyrin = pair_row(rows, BETA_AMYRIN, ALPHA_AMYRIN)
    beta, alpha = side(amyrin, BETA_AMYRIN), side(amyrin, ALPHA_AMYRIN)
    found |= {
        "beta_amyrin_enzymes": int(amyrin[f"strict_{beta}_enzyme_count"]),
        "alpha_amyrin_enzymes": int(amyrin[f"strict_{alpha}_enzyme_count"]),
        "amyrin_shared_enzymes": int(amyrin["strict_shared_enzyme_count"]),
        "beta_amyrin_bits": float(amyrin[f"strict_enzyme_{beta}_information_bits"]),
        "alpha_amyrin_bits": float(amyrin[f"strict_enzyme_{alpha}_information_bits"]),
        "amyrin_same_under_both_criteria": all(
            enzyme_set("strict", c) == enzyme_set("broader", c)
            and reaction_set("strict", c) == reaction_set("broader", c)
            for c in (BETA_AMYRIN, ALPHA_AMYRIN)
        ),
        "beta_amyrin_pathways": int(sets[("strict", BETA_AMYRIN)]["pathway_count"]),
    }
    hopanoids = pair_row(rows, DIPLOPTENE, TETRAHYMANOL)
    diplopterol = pair_row(rows, DIPLOPTENE, DIPLOPTEROL)
    found |= {
        "diploptene_enzymes": int(
            hopanoids[f"strict_{side(hopanoids, DIPLOPTENE)}_enzyme_count"]
        ),
        "tetrahymanol_enzymes": int(
            hopanoids[f"strict_{side(hopanoids, TETRAHYMANOL)}_enzyme_count"]
        ),
        "diploptene_tetrahymanol_shared_enzymes": max(
            int(hopanoids[f"{c}_shared_enzyme_count"]) for c in CRITERIA
        ),
        "diploptene_tetrahymanol_shared_pathways": json.loads(
            hopanoids["strict_shared_pathways"]
        ),
        "diploptene_diplopterol_enzymes": int(
            diplopterol["strict_shared_enzyme_count"]
        ),
        "diploptene_diplopterol_relations": [
            diplopterol[f"{c}_enzyme_relation"] for c in CRITERIA
        ],
    }
    identical = [
        row
        for row in rows
        for c in CRITERIA
        if row[f"{c}_enzyme_relation"] == "identical"
    ]
    disjoint = [
        row
        for row in rows
        for c in CRITERIA
        if row[f"{c}_enzyme_relation"] == "disjoint"
    ]
    found["identical_rows_f1_i0"] = all(
        row[f"{c}_enzyme_{s}_retained_fraction"] == "1.0"
        and row[f"{c}_enzyme_{s}_information_bits"] == "0.0"
        for row in identical
        for c in CRITERIA
        if row[f"{c}_enzyme_relation"] == "identical"
        for s in ("left", "right")
    )
    found["disjoint_rows_f0_i_blank"] = all(
        row[f"{c}_enzyme_{s}_retained_fraction"] == "0.0"
        and row[f"{c}_enzyme_{s}_information_bits"] == ""
        for row in disjoint
        for c in CRITERIA
        if row[f"{c}_enzyme_relation"] == "disjoint"
        for s in ("left", "right")
    )

    index = enzymes.enzyme_index(metacyc)
    genes = enzymes.read_dat(metacyc / "genes.dat", {"COMMON-NAME"})
    reactions = reaction_set("broader", BETA_AMYRIN) | reaction_set(
        "broader", ALPHA_AMYRIN
    )
    found["alpha_amyrin_enzyme_details"] = enzyme_details(
        enzyme_set("broader", ALPHA_AMYRIN), index, genes, reactions
    )
    hopanoid_reactions = reaction_set("broader", DIPLOPTENE) | reaction_set(
        "broader", TETRAHYMANOL
    )
    found["diploptene_enzyme_details"] = enzyme_details(
        enzyme_set("broader", DIPLOPTENE), index, genes, hopanoid_reactions
    )
    found["tetrahymanol_enzyme_details"] = enzyme_details(
        enzyme_set("broader", TETRAHYMANOL), index, genes, hopanoid_reactions
    )
    found["diploptene_bacterial_enzymes"] = sum(
        "bacteria" in d["groups"] for d in found["diploptene_enzyme_details"]
    )
    found["diploptene_fern_enzymes"] = sum(
        "ferns" in d["groups"] for d in found["diploptene_enzyme_details"]
    )
    pathways = enzymes.read_dat(metacyc / "pathways.dat", {"COMMON-NAME"})
    found["shared_pathway_names"] = {
        p: enzymes.display_name(pathways[p]["COMMON-NAME"][0])
        for p in found["diploptene_tetrahymanol_shared_pathways"]
    }
    reaction_records = enzymes.read_dat(
        metacyc / "reactions.dat", {"LEFT", "RIGHT", "EC-NUMBER"}
    )
    found["example_reactions"] = {
        r: " + ".join(reaction_records[r].get("LEFT", []))
        + " = "
        + " + ".join(reaction_records[r].get("RIGHT", []))
        for r in sorted(reactions | hopanoid_reactions)
    }
    return found


def print_enzymes(title: str, details: list[dict]) -> None:
    print(title)
    for d in details:
        print(
            f"  {d['enzyme']}: {d['name']}; {', '.join(d['organisms'])} "
            f"({', '.join(d['groups'])}); producing reactions {', '.join(d['reactions'])}"
        )


def report(n: dict) -> None:
    """Print the numbers as plain sentences."""
    print(f"Pairs: {n['pairs']:,} unordered pairs of {n['biomarkers']} biomarkers.")
    print(
        f"Beta-amyrin has {n['beta_amyrin_enzymes']} enzymes and alpha-amyrin "
        f"{n['alpha_amyrin_enzymes']}; {n['amyrin_shared_enzymes']} are shared. Assuming a "
        f"shared enzyme narrows beta-amyrin by {n['beta_amyrin_bits']:.3f} bits and "
        f"alpha-amyrin by {n['alpha_amyrin_bits']:.3f} bits. Same sets under strict and "
        f"broader criteria: {n['amyrin_same_under_both_criteria']}."
    )
    print(
        f"Pathways containing a producing reaction of beta-amyrin (strict): {n['beta_amyrin_pathways']}."
    )
    print_enzymes("Enzymes of alpha-amyrin:", n["alpha_amyrin_enzyme_details"])
    print(
        f"Diploptene has {n['diploptene_enzymes']} enzymes ({n['diploptene_bacterial_enzymes']} "
        f"from bacteria, {n['diploptene_fern_enzymes']} from a fern) and tetrahymanol "
        f"{n['tetrahymanol_enzymes']}; shared: {n['diploptene_tetrahymanol_shared_enzymes']}."
    )
    print_enzymes("Enzymes of diploptene:", n["diploptene_enzyme_details"])
    print_enzymes("Enzymes of tetrahymanol:", n["tetrahymanol_enzyme_details"])
    print(
        "Pathways shared by diploptene and tetrahymanol (strict): "
        + "; ".join(f"{p} {name}" for p, name in n["shared_pathway_names"].items())
    )
    print("Producing reactions of the examples, as written in reactions.dat:")
    for reaction, equation in n["example_reactions"].items():
        print(f"  {reaction}: {equation}")
    print(
        f"Diploptene and diplopterol share {n['diploptene_diplopterol_enzymes']} enzymes; "
        f"relation under strict and broader criteria: {n['diploptene_diplopterol_relations']}."
    )
    for c in CRITERIA:
        print(
            f"{c.capitalize()} criteria: {n[f'{c}_both_with_enzymes']:,} pairs with enzymes on "
            f"both sides ({n[f'{c}_without_enzymes_one_or_both']:,} without for one or both); "
            f"{n[f'{c}_sharing_enzyme']} share an enzyme, {n[f'{c}_identical_enzymes']} with "
            f"identical sets ({n[f'{c}_identical_enzyme_pairs']}); "
            f"{n[f'{c}_overlap_with_positive_reduction']} with a positive reduction in at least "
            f"one direction; median reduction {n[f'{c}_median_bits']:.2f} bits over "
            f"{n[f'{c}_median_directions']} directions."
        )
        print(
            f"  {n[f'{c}_disjoint_enzymes']:,} disjoint, {n[f'{c}_left_empty']:,} with only the "
            f"left set empty, {n[f'{c}_right_empty']:,} with only the right set empty, "
            f"{n[f'{c}_both_empty']:,} with both empty."
        )
        print(
            f"  {n[f'{c}_precursor_pairs']} precursor pairs, of which "
            f"{n[f'{c}_precursor_pairs_sharing_enzyme']} share an enzyme; "
            f"{n[f'{c}_precursor_no']:,} tested pairs are not precursor pairs and "
            f"{n[f'{c}_precursor_not_tested']:,} pairs were not tested."
        )
        print(
            f"  {n[f'{c}_sharing_pathway']} pairs share a pathway, "
            f"{n[f'{c}_identical_pathways']} with identical pathway sets."
        )
    print(
        "Identical enzyme sets have retained fractions of 1 and reductions of 0: "
        f"{n['identical_rows_f1_i0']}. Disjoint sets have retained fractions of 0 and no "
        f"reduction: {n['disjoint_rows_f0_i_blank']}."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True, help="results folder")
    parser.add_argument(
        "--metacyc", type=Path, required=True, help="MetaCyc 29.5 data folder"
    )
    args = parser.parse_args()
    report(numbers(args.results, args.metacyc))
