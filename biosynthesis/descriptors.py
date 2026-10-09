"""Molecular descriptors: cohort selection, Böttcher score treatments and calculation.

The cohort and score rules work on data/descriptors.csv. The calculation
functions need RDKit (and the Böttcher scorer with Open Babel); they import them
only when called.
"""

from __future__ import annotations

import os
import sys
from typing import Mapping

import numpy as np
import pandas as pd

# Selection of the descriptor cohort from the reference graph.
REFERENCE_SCORE_LIMIT = 1_000_000
CURRENCY_COMPOUNDS = frozenset(
    {
        "PROTON",
        "WATER",
        "ATP",
        "ADP",
        "AMP",
        "NAD",
        "NADH",
        "NADP",
        "NADPH",
        "CO-A",
        "Acceptor",
        "Donor",
        "FAD",
        "FADH2",
        "PPI",
        "Pi",
        "CO2",
        "O2",
        "Reduced-ferredoxins",
        "Oxidized-ferredoxins",
        "Reduced-flavodoxins",
        "UDP",
        "UTP",
        "GDP",
        "GTP",
        "CTP",
        "CDP",
        "TTP",
        "AMMONIA",
    }
)

# Versions used to compute data/descriptors.csv: RDKit for the 12 descriptors, RDKit for
# structural eligibility, and the Open Babel wheel (library 3.1.0) of the Böttcher scorer.
RDKIT_DESCRIPTOR_VERSION = "2023.03.3"
RDKIT_ELIGIBILITY_VERSION = "2025.09.4"
OPENBABEL_WHEEL_VERSION = "3.1.1.23"
OPENBABEL_LIBRARY_VERSION = "3.1.0"
BOTTCHER_SCORER_REVISION = "6f1368cd349f14a3e719258f50a71fea0b526052"

# Böttcher scoring (forlilab/bottchscore): first pass and the larger symmetry search.
BOTTCHER_SETTINGS = dict(disable_mesomer=True, automorp_memory_maxsize=10000)
LARGER_SYMMETRY_SEARCH = 3_000_000
# Warning labels and the scorer messages that identify them.
WARNING_PATTERNS = {
    "automorphism_memory_limit": "memory limit exceeded",
    "cis_trans_stereochemistry": "Error in cis/trans stereochemistry",
    "aromatic_kekulization": "Failed to kekulize",
    "atomic_score_error": "Error calculating complexity: atom_idx",
}
FAILURE_WARNINGS = ("automorphism_memory_limit", "atomic_score_error")

# Böttcher status of each cohort compound.
NO_WARNING, RERUN, WITHHELD, NO_SCORE = (
    "no_warning",
    "rerun_larger_search",
    "withheld",
    "no_score",
)
# Score treatments compared in Table S4, in table order.
SCORE_TREATMENTS = (
    "unscreened",
    "larger_search_values",
    "withhold_failures",
    "withhold_all_warnings",
    "exclude_flagged_compounds",
    "main_analysis",
)


def select_cohort(records: pd.DataFrame) -> pd.DataFrame:
    """Selection table: a compound is selected when it has a path in the reference graph,
    a cumulative score no greater than the limit, and is not a currency compound.

    `records` has compound_id and the reference num_paths and cumulative score.
    """
    table = records[["compound_id"]].copy()
    table["num_paths_reference_graph"] = records["num_paths_reference_graph"]
    table["cumulative_score_reference_graph"] = records[
        "cumulative_score_reference_graph"
    ]
    matched = table.num_paths_reference_graph.notna()
    table["has_path_reference_graph"] = table.num_paths_reference_graph.gt(0)
    within = table.cumulative_score_reference_graph.le(REFERENCE_SCORE_LIMIT)
    currency = records.compound_id.isin(CURRENCY_COMPOUNDS)
    table["selected"] = matched & table.has_path_reference_graph & within & ~currency
    assert table.compound_id.is_unique
    return table


def treatment_sets(cohort: pd.DataFrame) -> tuple[dict[str, float], set[str], set[str]]:
    """Larger-search replacement scores, failure-flagged and all flagged compounds.

    Failures are warnings of a memory-limited symmetry search or a failed atomic
    term; flagged compounds have any warning.
    """
    rerun = cohort.loc[cohort.bottcher_status.eq(RERUN)]
    replacements = {
        compound: float(value)
        for compound, value in zip(rerun.compound_id, rerun.bottchscore)
    }
    warned = cohort.loc[cohort.bottcher_warnings.ne("")]
    all_flags = set(warned.compound_id)
    failures = set(
        warned.loc[
            warned.bottcher_warnings.str.contains("|".join(FAILURE_WARNINGS)),
            "compound_id",
        ]
    )
    return replacements, failures, all_flags


def classify_warnings(text: str) -> list[str]:
    """Warning labels found in the scorer messages of one compound."""
    labels = [label for label, pattern in WARNING_PATTERNS.items() if pattern in text]
    return labels or ["unclassified"]


def bottcher_status(
    unscreened: str, warnings: str, replacement: str | None
) -> tuple[str, str]:
    """Status and analysis score of one compound from its unscreened score and warnings.

    A compound with warnings keeps a score only when the larger symmetry search
    gave a score with no failed term, no message and a completed search.
    """
    if warnings and replacement is not None:
        return RERUN, replacement
    if warnings:
        return WITHHELD, ""
    if unscreened == "":
        return NO_SCORE, ""
    return NO_WARNING, unscreened


def quiet_rdkit() -> None:
    """Report only RDKit errors, not warnings about the structures it reads."""
    from rdkit import RDLogger

    RDLogger.logger().setLevel(RDLogger.ERROR)


def sa_scorer():
    """The synthetic-accessibility scorer of the RDKit Contrib directory."""
    from rdkit.Chem import RDConfig

    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer

    return sascorer


def compute_structural_descriptors(smiles: str, sascorer=None):
    """The 12 RDKit descriptors of a SMILES string, or None when it does not parse."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, GraphDescriptors, rdMolDescriptors

    sascorer = sascorer or sa_scorer()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        sa_score = sascorer.calculateScore(mol)
    except Exception:
        sa_score = np.nan
    try:
        bertz_ct = GraphDescriptors.BertzCT(mol)
    except Exception:
        bertz_ct = np.nan
    return {
        "sa_score": sa_score,
        "bertz_ct": bertz_ct,
        "stereocenters": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "fsp3": Descriptors.FractionCSP3(mol),
        "ring_count": Descriptors.RingCount(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "hb_acceptors": Descriptors.NumHAcceptors(mol),
        "hb_donors": Descriptors.NumHDonors(mol),
        "heavy_atoms": Descriptors.HeavyAtomCount(mol),
        "tpsa": Descriptors.TPSA(mol),
        "logp": Descriptors.MolLogP(mol),
        "mw": Descriptors.MolWt(mol),
    }


def structural_eligibility(smiles: object) -> Mapping[str, object]:
    """Whether a structure parses, is connected, contains carbon and has two heavy atoms."""
    from rdkit import Chem

    if not isinstance(smiles, str) or not smiles.strip():
        return {"structurally_eligible": False, "eligibility_reason": "missing_smiles"}
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return {
            "structurally_eligible": False,
            "eligibility_reason": "rdkit_parse_failure",
        }
    component_count = len(Chem.GetMolFrags(molecule))
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in molecule.GetAtoms())
    heavy_atom_count = molecule.GetNumHeavyAtoms()
    reasons: list[str] = []
    if component_count != 1:
        reasons.append("disconnected_structure")
    if carbon_count < 1:
        reasons.append("no_carbon_atom")
    if heavy_atom_count < 2:
        reasons.append("fewer_than_two_heavy_atoms")
    return {
        "structurally_eligible": not reasons,
        "eligibility_reason": ";".join(reasons),
    }
