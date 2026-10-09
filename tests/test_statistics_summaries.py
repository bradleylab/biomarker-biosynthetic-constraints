"""Population summaries, class means and the descriptor correlations and regression."""

import math

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

from biosynthesis.descriptors import bottcher_status, select_cohort, treatment_sets
from biosynthesis.paths import metric_frame
from biosynthesis.statistics import (
    FEATURES,
    class_summaries,
    descriptor_statistics,
    population_summary,
)


def raw(num_paths, score, lam, sigma):
    return {
        "num_paths": float(num_paths),
        "total_complexity": float(score),
        "mean_branching_depth_normalized": float(lam),
        "fraction_shared_80": float(sigma),
        "max_depth": 1.0,
    }


def test_class_summaries_use_compounds_with_a_path():
    frame = metric_frame(
        [raw(0, 0, 0, 0), raw(2, 2, 0.5, 0.8), raw(4, 8, 0.25, 0.4), raw(1, 0, 0, 1)]
    )
    frame["class"] = ["A", "A", "A", "B"]
    rows = class_summaries(frame)
    assert [row["class"] for row in rows] == ["A", "B"]
    first, second = rows
    assert (first["n"], first["n_with_path"]) == (3, 2)
    assert first["mean_psi"] == 1.5 and first["mean_lambda"] == 0.375
    assert first["sd_psi"] == pytest.approx(math.sqrt(0.5))
    assert (second["mean_psi"], second["mean_sigma"]) == (0.0, 1.0)
    assert math.isnan(second["sd_psi"])


def test_population_summary_counts_and_screen():
    frame = metric_frame(
        [
            raw(0, 0, 0, 0),
            raw(1, 0, 0, 1),
            raw(2, 2, 0.8, 0.9),
            raw(4, 8, 0.2, 0.3),
            raw(8, 24, 0.5, 0.6),
        ]
    )
    summary = population_summary(frame)
    assert (
        summary["positive_path_n"],
        summary["one_path_n"],
        summary["multi_path_n"],
    ) == (4, 1, 3)
    assert summary["multi_path_psi_median"] == 2.0
    assert summary["few_paths_far_branching_n"] == 1  # psi = 1, lambda 0.8, sigma 0.9
    assert summary["few_paths_far_branching_percent"] == pytest.approx(100 / 3)


def synthetic_descriptors(rows=60, seed=7):
    rng = np.random.default_rng(seed)
    frame = pd.DataFrame(rng.normal(size=(rows, len(FEATURES))), columns=FEATURES)
    frame.insert(0, "compound_id", [f"C{i}" for i in range(rows)])
    frame["psi_bits"] = frame[FEATURES].sum(axis=1) + rng.normal(size=rows)
    frame["lambda"] = rng.uniform(size=rows)
    frame["sigma"] = rng.uniform(size=rows)
    frame.loc[3, "bottchscore"] = np.nan
    frame.loc[5, ["psi_bits", "lambda", "sigma"]] = np.nan
    return frame


def test_descriptor_statistics_against_a_separate_calculation():
    frame = synthetic_descriptors()
    rows, _ = descriptor_statistics(frame)
    for row in rows:
        if row["analysis_type"] == "spearman":
            pair = frame[[row["descriptor"], row["target"]]].dropna()
            assert len(pair) == row["n"]
            rho = np.corrcoef(rankdata(pair.iloc[:, 0]), rankdata(pair.iloc[:, 1]))[
                0, 1
            ]
            assert math.isclose(rho, row["spearman_rho"], rel_tol=1e-12, abs_tol=1e-12)
        else:
            complete = frame[FEATURES + ["psi_bits"]].dropna()
            assert len(complete) == row["complete_case_n"] == 58
            x, y = complete[FEATURES].to_numpy(), complete.psi_bits.to_numpy()
            x = (x - x.mean(axis=0)) / x.std(axis=0)
            design = np.column_stack([np.ones(len(x)), x])
            fit = design @ np.linalg.lstsq(design, y, rcond=None)[0]
            r_squared = 1 - np.square(y - fit).sum() / np.square(y - y.mean()).sum()
            assert math.isclose(
                r_squared, row["r_squared"], rel_tol=1e-12, abs_tol=1e-12
            )
            adjusted = 1 - (1 - r_squared) * (len(y) - 1) / (len(y) - x.shape[1] - 1)
            assert math.isclose(
                adjusted, row["adjusted_r_squared"], rel_tol=1e-12, abs_tol=1e-12
            )


def test_cohort_selection_rule():
    records = pd.DataFrame(
        {
            "compound_id": ["A", "B", "C", "ATP"],
            "num_paths_reference_graph": [3, 0, 5, 2],
            "cumulative_score_reference_graph": [4.0, 0.0, 2e6, 1.0],
        }
    )
    assert select_cohort(records).selected.tolist() == [True, False, False, False]


def test_bottcher_status_and_treatments():
    assert bottcher_status("12.5", "", None) == ("no_warning", "12.5")
    assert bottcher_status("", "", None) == ("no_score", "")
    assert bottcher_status("12.5", "cis_trans_stereochemistry", None) == (
        "withheld",
        "",
    )
    assert bottcher_status("12.5", "automorphism_memory_limit", "13.0") == (
        "rerun_larger_search",
        "13.0",
    )
    cohort = pd.DataFrame(
        {
            "compound_id": ["A", "B", "C", "D"],
            "bottchscore": ["1.0", "", "13.0", ""],
            "bottcher_status": [
                "no_warning",
                "withheld",
                "rerun_larger_search",
                "withheld",
            ],
            "bottcher_warnings": [
                "",
                "cis_trans_stereochemistry",
                "automorphism_memory_limit",
                "atomic_score_error",
            ],
        }
    )
    replacements, failures, flagged = treatment_sets(cohort)
    assert replacements == {"C": 13.0}
    assert failures == {"C", "D"}
    assert flagged == {"B", "C", "D"}
