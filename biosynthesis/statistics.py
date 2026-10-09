"""Population summaries, class means, and descriptor correlations and regressions."""

from __future__ import annotations

from typing import Mapping

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

FEATURES = [
    "bottchscore",
    "sa_score",
    "bertz_ct",
    "stereocenters",
    "fsp3",
    "ring_count",
    "rotatable_bonds",
    "hb_acceptors",
    "hb_donors",
    "heavy_atoms",
    "tpsa",
    "logp",
    "mw",
]
TARGETS = ["psi_bits", "lambda", "sigma"]
# Population summaries are also given after excluding compounds whose cumulative
# branching score exceeds this value.
SCORE_LIMIT = 1_000_000
# Few paths that branch far upstream and overlap heavily: psi <= 2, lambda > 0.7, sigma > 0.5.
SCREEN_PSI_MAX, SCREEN_LAMBDA_MIN, SCREEN_SIGMA_MIN = 2, 0.7, 0.5


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Convert the named columns to numbers; blank cells are the only missing token."""
    result = frame.copy()
    for column in columns:
        result[column] = pd.to_numeric(
            result[column].replace("", np.nan), errors="raise"
        )
        assert not np.isinf(result[column].dropna()).any()
    return result


def descriptor_statistics(frame: pd.DataFrame) -> tuple[list[dict], dict]:
    """Spearman correlations of each descriptor with each measurement, and the OLS fit.

    Each correlation uses the compounds with both values; the regression of psi on
    the 13 standardized descriptors uses the complete cases.
    """
    results = []
    accounting = dict(
        analysis_rows=len(frame),
        missing_ids_by_column={
            c: frame.loc[frame[c].isna(), "compound_id"].tolist()
            for c in FEATURES + TARGETS
        },
    )
    for feature in FEATURES:
        for target in TARGETS:
            valid = frame[feature].notna() & frame[target].notna()
            rho, p = spearmanr(frame.loc[valid, feature], frame.loc[valid, target])
            results.append(
                dict(
                    analysis_type="spearman",
                    descriptor=feature,
                    target=target,
                    n=int(valid.sum()),
                    spearman_rho=float(rho),
                    p_value=float(p),
                )
            )
    valid = frame[FEATURES + ["psi_bits"]].notna().all(axis=1)
    complete = frame.loc[valid]
    design = StandardScaler().fit_transform(complete[FEATURES])
    response = complete.psi_bits.to_numpy()
    model = LinearRegression().fit(design, response)
    r_squared = float(model.score(design, response))
    n, k = design.shape
    results.append(
        dict(
            analysis_type="ols_standardized",
            descriptor="",
            target="",
            analysis_rows=len(frame),
            response_evaluable_n=int(frame.psi_bits.notna().sum()),
            complete_case_n=n,
            feature_count=k,
            r_squared=r_squared,
            adjusted_r_squared=1 - (1 - r_squared) * (n - 1) / (n - k - 1),
        )
    )
    accounting["regression_excluded_ids"] = frame.loc[~valid, "compound_id"].tolist()
    accounting["regression_included_ids"] = complete.compound_id.tolist()
    accounting["standardized_coefficients"] = dict(zip(FEATURES, model.coef_.tolist()))
    accounting["intercept"] = float(model.intercept_)
    return results, accounting


def population_summary(frame: pd.DataFrame) -> dict[str, object]:
    """Counts, medians, mean, maximum and Pearson correlations of a compound population.

    Counts cover compounds with at least one path, one path and several paths; the
    other values use the compounds with several paths.
    """
    positive = frame.loc[frame["num_paths"].gt(0)].copy()
    multi = positive.loc[positive["num_paths"].gt(1)].copy()
    correlations: dict[str, float] = {}
    for left, right in [
        ("psi_bits", "lambda"),
        ("psi_bits", "sigma"),
        ("lambda", "sigma"),
    ]:
        valid = multi[left].notna() & multi[right].notna()
        name = f"pearson_{left.split('_')[0]}_{right.split('_')[0]}"
        correlations[name] = float(
            stats.pearsonr(multi.loc[valid, left], multi.loc[valid, right]).statistic
        )
    screen = (
        multi["psi_bits"].le(SCREEN_PSI_MAX)
        & multi["lambda"].gt(SCREEN_LAMBDA_MIN)
        & multi["sigma"].gt(SCREEN_SIGMA_MIN)
    )
    return {
        "positive_path_n": int(len(positive)),
        "one_path_n": int(positive["num_paths"].eq(1).sum()),
        "multi_path_n": int(len(multi)),
        "multi_path_psi_median": float(multi["psi_bits"].median()),
        "multi_path_psi_mean": float(multi["psi_bits"].mean()),
        "multi_path_psi_max": float(multi["psi_bits"].max()),
        "multi_path_lambda_median": float(multi["lambda"].median()),
        "multi_path_sigma_median": float(multi["sigma"].median()),
        "few_paths_far_branching_n": int(screen.sum()),
        "few_paths_far_branching_percent": float(100 * screen.mean()),
        **correlations,
    }


def class_summaries(panel: pd.DataFrame) -> list[dict]:
    """Size, number with a path, and mean and SD (n - 1) of psi, lambda and sigma per class.

    Classes appear in the order of the panel rows; compounds without a path have
    missing values and do not enter the means.
    """
    rows = []
    for label, group in panel.groupby("class", sort=False):
        row = {
            "class": label,
            "n": len(group),
            "n_with_path": int(group.num_paths.gt(0).sum()),
        }
        for short, column in (
            ("psi", "psi_bits"),
            ("lambda", "lambda"),
            ("sigma", "sigma"),
        ):
            row[f"mean_{short}"] = float(group[column].mean())
            row[f"sd_{short}"] = float(group[column].std(ddof=1))
        rows.append(row)
    return rows


def strongest_psi_correlation(rows: list[Mapping]) -> Mapping:
    """The correlation with psi of largest absolute value (first one on ties)."""
    return max(
        (row for row in rows if row["target"] == "psi_bits"),
        key=lambda row: abs(row["spearman_rho"]),
    )
