"""Associations of the molecular descriptors with the projected-path measurements (Table S4).

For both graph versions and two cohorts (all compounds, and the structurally
eligible ones), computes Spearman correlations of each descriptor with psi, lambda
and sigma and an OLS regression of psi on the 13 standardized descriptors, under
six treatments of the Böttcher score. Writes:

  TableS4_associations.csv              the main analysis
  TableS4_score_treatment_summary.csv   regression fit and strongest psi correlation per treatment
  TableS4_diagnostic_sensitivity.csv    every statistic under every treatment, with
                                        differences from the unscreened scores

  python analyses/02_projected_paths/descriptor_associations.py \\
      --cohort results/descriptors --out results/descriptors
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biosynthesis.descriptors import treatment_sets  # noqa: E402
from biosynthesis.simplified_graph import NETWORK_VERSIONS  # noqa: E402
from biosynthesis.statistics import (  # noqa: E402
    FEATURES,
    TARGETS,
    descriptor_statistics,
    numeric,
    strongest_psi_correlation,
)
from biosynthesis.tables import write_csv  # noqa: E402

DATA = REPO / "data"
MAIN = "main_analysis"
STATISTIC_COLUMNS = [
    "analysis_type",
    "descriptor",
    "target",
    "n",
    "spearman_rho",
    "p_value",
    "analysis_rows",
    "response_evaluable_n",
    "complete_case_n",
    "feature_count",
    "r_squared",
    "adjusted_r_squared",
]
ASSOCIATION_COLUMNS = ["network_version", "cohort", *STATISTIC_COLUMNS]
SUMMARY_COLUMNS = [
    "network_version",
    "cohort",
    "score_treatment",
    "analysis_rows",
    "regression_complete_case_n",
    "r_squared",
    "adjusted_r_squared",
    "largest_absolute_psi_rho",
    "strongest_psi_descriptor",
    "signed_strongest_psi_rho",
]
COMPARED = ["n", "spearman_rho", "complete_case_n", "r_squared", "adjusted_r_squared"]
SENSITIVITY_COLUMNS = [
    "network_version",
    "score_treatment",
    "cohort",
    *STATISTIC_COLUMNS,
    *["unscreened_" + field for field in COMPARED],
    *["delta_" + field for field in COMPARED],
]
INTEGER_COLUMNS = {
    "n",
    "analysis_rows",
    "response_evaluable_n",
    "complete_case_n",
    "feature_count",
    "unscreened_n",
    "unscreened_complete_case_n",
    "delta_n",
    "delta_complete_case_n",
    "regression_complete_case_n",
}
# Treatments whose differences from the unscreened scores are not reported: the
# whole-record exclusion and the main analysis change the compared populations.
NO_DIFFERENCES = ("exclude_flagged_compounds", MAIN)


def plain(column: str, value: object) -> object:
    """A table cell: missing as None, counts as integers, other numbers as floats."""
    if value is None or value == "":
        return value
    if isinstance(value, (float, np.floating)) and math.isnan(value):
        return None
    if column in INTEGER_COLUMNS:
        return int(value)
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def eligible_compounds() -> set[str]:
    """Compounds whose structures parse, are connected, contain carbon and have two heavy atoms."""
    records = pd.read_csv(DATA / "descriptors.csv", dtype=str, keep_default_na=False)
    return set(records.loc[records.structurally_eligible.eq("true"), "compound_id"])


def responses(cohort: pd.DataFrame) -> pd.DataFrame:
    """psi_bits, lambda and sigma of the cohort as numbers; missing without a path."""
    frame = cohort[["compound_id"]].copy()
    for target in TARGETS:
        frame[target] = [float(value) if value else np.nan for value in cohort[target]]
    return frame


def treatment_scenarios(cohort: pd.DataFrame) -> dict[str, tuple]:
    """(use analysis scores, replacement scores, withheld scores, excluded compounds)."""
    replacements, failures, all_flags = treatment_sets(cohort)
    return {
        "unscreened": (False, {}, set(), set()),
        "larger_search_values": (False, replacements, set(), set()),
        "withhold_failures": (False, {}, failures, set()),
        "withhold_all_warnings": (False, {}, all_flags, set()),
        "exclude_flagged_compounds": (False, {}, set(), all_flags),
        MAIN: (True, {}, set(), set()),
    }


def version_statistics(
    cohort: pd.DataFrame, eligible: set[str], version: str
) -> tuple[list[dict], list[dict]]:
    """Statistics rows and summary rows of one graph version under every treatment."""
    unscreened = cohort.copy()
    unscreened["bottchscore"] = unscreened["bottchscore_unscreened"]
    base = numeric(unscreened, FEATURES)
    accepted = numeric(cohort, FEATURES)
    measured = responses(cohort)
    results, summary = [], []
    for treatment, (use_current, replace, mask, excluded) in treatment_scenarios(
        cohort
    ).items():
        descriptors = (accepted if use_current else base)[
            ["compound_id", *FEATURES]
        ].copy()
        replace_mask = descriptors.compound_id.isin(replace)
        descriptors.loc[replace_mask, "bottchscore"] = descriptors.loc[
            replace_mask, "compound_id"
        ].map(replace)
        descriptors.loc[descriptors.compound_id.isin(mask), "bottchscore"] = np.nan
        data = descriptors.merge(measured, on="compound_id", how="left", validate="1:1")
        data = data.loc[~data.compound_id.isin(excluded)]
        for population, subset in (
            ("all", data),
            ("structurally_eligible", data.loc[data.compound_id.isin(eligible)]),
        ):
            rows, _ = descriptor_statistics(subset)
            results.extend(
                {
                    "network_version": version,
                    "score_treatment": treatment,
                    "cohort": population,
                    **row,
                }
                for row in rows
            )
            fit = next(
                row for row in rows if row["analysis_type"] == "ols_standardized"
            )
            strongest = strongest_psi_correlation(rows)
            summary.append(
                {
                    "network_version": version,
                    "cohort": population,
                    "score_treatment": treatment,
                    "analysis_rows": len(subset),
                    "regression_complete_case_n": fit["complete_case_n"],
                    "r_squared": fit["r_squared"],
                    "adjusted_r_squared": fit["adjusted_r_squared"],
                    "largest_absolute_psi_rho": abs(strongest["spearman_rho"]),
                    "strongest_psi_descriptor": strongest["descriptor"],
                    "signed_strongest_psi_rho": strongest["spearman_rho"],
                }
            )
    return results, summary


def with_differences(results: list[dict]) -> pd.DataFrame:
    """Each statistic with its unscreened value and the difference from it."""
    table = pd.DataFrame(results)
    keys = ["network_version", "cohort", "analysis_type", "descriptor", "target"]
    reference = table.loc[
        table.score_treatment.eq("unscreened"), keys + COMPARED
    ].rename(columns={field: "unscreened_" + field for field in COMPARED})
    table = table.merge(reference, on=keys, validate="m:1")
    for field in COMPARED:
        table["delta_" + field] = table[field] - table["unscreened_" + field]
    table.loc[
        table.score_treatment.isin(NO_DIFFERENCES),
        ["unscreened_" + field for field in COMPARED]
        + ["delta_" + field for field in COMPARED],
    ] = np.nan
    return table


def write_table(path: Path, rows: list[dict], columns: list[str]) -> None:
    write_csv(
        path,
        (
            {column: plain(column, row.get(column)) for column in columns}
            for row in rows
        ),
        columns,
    )


def main(cohort_dir: Path, out: Path) -> None:
    cohort = pd.read_csv(
        cohort_dir / "TableS4_descriptor_cohort.csv", dtype=str, keep_default_na=False
    )
    eligible = eligible_compounds()
    results, summary = [], []
    for version in NETWORK_VERSIONS:
        rows, fits = version_statistics(
            cohort.loc[cohort.network_version.eq(version)].reset_index(drop=True),
            eligible,
            version,
        )
        results.extend(rows)
        summary.extend(fits)
    associations = [row for row in results if row["score_treatment"] == MAIN]
    write_table(out / "TableS4_associations.csv", associations, ASSOCIATION_COLUMNS)
    write_table(out / "TableS4_score_treatment_summary.csv", summary, SUMMARY_COLUMNS)
    sensitivity = with_differences(results).to_dict("records")
    write_table(
        out / "TableS4_diagnostic_sensitivity.csv", sensitivity, SENSITIVITY_COLUMNS
    )
    print(f"{len(associations)} association rows, {len(sensitivity)} treatment rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--cohort", type=Path, required=True, help="output of descriptor_cohort.py"
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.cohort, args.out)
