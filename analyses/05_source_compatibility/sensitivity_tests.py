"""Sensitivity tests of the source-compatibility analysis (Table S7, last five files).

Reads sample_scenarios.jsonl and detected_class_sets.jsonl (written by
compatible_combinations.py) and writes:

- TableS7_thresholds_by_sample.csv: for each sample, set of assignments and
  low-signal treatment, how the result changes across detection thresholds,
  with scenarios without detections kept apart;
- TableS7_delegation_by_sample.csv: for each sample under strict assignments
  at a threshold of zero with low signals included, the smallest number of
  detected classes that must come from sources outside the five categories to
  change the minimum number of categories, and the kind of change;
- TableS7_evidence_removal.csv, TableS7_omissions.csv, TableS7_additions.csv:
  the results of removing one evidence record, omitting one detected class and
  adding one undetected class;
- delegation_summary.csv: for each set of detected classes, the smallest
  delegation that changes each conclusion;
- delegated_class_sets.csv: every smallest delegation behind
  TableS7_delegation_by_sample.csv.

A change of evaluability, of the minimum number of categories and of the
compatible combinations are recorded as separate kinds of change, which may
overlap.

Usage:
  python analyses/05_source_compatibility/sensitivity_tests.py --results results/sources --out results/sources
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.sensitivity import outcome, transition_kind, zero_threshold  # noqa: E402
from biosynthesis.sources import BROADER, STRICT  # noqa: E402
from biosynthesis.tables import read_jsonl, write_csv  # noqa: E402

COMMON = ("analysis_id", "assignment_criteria", "detected_analyte_ids")
BASELINE_CONCLUSION = "min_categories"


def joined(values) -> str:
    return ";".join(str(value) for value in values)


def json_text(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def per_analysis_tables(analyses):
    """Evidence removal, omission, addition and delegation-summary rows of every analysis."""
    removals, omissions, additions, summary = [], [], [], []
    for analysis in analyses:
        common = {
            "analysis_id": analysis["analysis_id"],
            "assignment_criteria": analysis["assignment_criteria"],
            "detected_analyte_ids": joined(analysis["detected_analyte_ids"]),
        }
        for row in analysis["evidence_removal_test"]["removals"]:
            removals.append(
                {
                    **common,
                    "changed_conclusions": joined(row["changed_conclusions"]),
                    "evidence_kind": row["evidence_kind"],
                    "evidence_record_id": row["evidence_record_id"],
                    "removed_assignment_options": row["removed_assignment_options"],
                }
            )
        for row in analysis["omission_test"]["omissions"]:
            omissions.append(
                {
                    **common,
                    "analyte_id": row["analyte_id"],
                    "minimal_combinations_no_longer_minimal": row[
                        "minimal_combinations_no_longer_minimal"
                    ],
                    "compatible_combinations_gained": row[
                        "compatible_combinations_gained"
                    ],
                    "new_combinations": row["new_combinations"],
                }
            )
        test = analysis["addition_test"]
        for row in test["additions"] if test else ():
            additions.append(
                {
                    "analysis_id": analysis["analysis_id"],
                    "candidate_analyte_id": row["candidate_analyte_id"],
                    "status": row["status"],
                    "reason_code": row["reason_code"],
                    "baseline_compatible_combinations": row[
                        "baseline_compatible_combinations"
                    ],
                    "compatible_combinations_after_addition": row[
                        "compatible_combinations_after_addition"
                    ],
                    "eliminated_combinations": row["eliminated_combinations"],
                    "reduction_bits": row["reduction_bits"],
                }
            )
        for boundary in analysis["delegation_test"]["boundaries"]:
            summary.append(
                {
                    "analysis_id": analysis["analysis_id"],
                    "conclusion_id": boundary["conclusion_id"],
                    "baseline_value": json_text(boundary["baseline_value"]),
                    "min_delegated_classes": boundary["min_delegated_classes"],
                    "delegations_found": len(boundary["changing_delegations"]),
                }
            )
    additions.sort(key=lambda row: (row["analysis_id"], row["candidate_analyte_id"]))
    summary.sort(key=lambda row: (row["analysis_id"], row["conclusion_id"]))
    return removals, omissions, additions, summary


def delegation_by_sample(scenarios, analyses, samples):
    """Smallest delegation changing the minimum number of categories, per sample (strict, zero, included)."""
    baselines = [
        row
        for row in scenarios
        if row["assignment_criteria"] == STRICT
        and row["low_signal_mode"] == "include"
        and zero_threshold(row)
    ]
    if sorted(row["sample_id"] for row in baselines) != samples:
        raise ValueError("expected one zero-threshold scenario per sample")
    delegation, details = [], []
    for scenario in sorted(baselines, key=lambda row: row["sample_id"]):
        analysis = analyses[scenario["analysis_id"]]
        base = analysis["result"]
        (boundary,) = [
            row
            for row in analysis["delegation_test"]["boundaries"]
            if row["conclusion_id"] == BASELINE_CONCLUSION
        ]
        kinds = set()
        for index, item in enumerate(boundary["changing_delegations"]):
            changed = item["summary"]
            kind = transition_kind(
                base["status"],
                base["min_categories"],
                changed["status"],
                changed["min_categories"],
            )
            kinds.add(kind)
            details.append(
                {
                    "sample_id": scenario["sample_id"],
                    "scenario_id": scenario["scenario_id"],
                    "analysis_id": scenario["analysis_id"],
                    "delegation_index": index,
                    "min_delegated_classes": boundary["min_delegated_classes"],
                    "delegated_analyte_ids": joined(item["delegated_analyte_ids"]),
                    "remaining_analyte_ids": joined(changed["detected_analyte_ids"]),
                    "before_status": base["status"],
                    "before_minimum": base["min_categories"],
                    "after_status": changed["status"],
                    "after_minimum": changed["min_categories"],
                    "transition_kind": kind,
                }
            )
        delegation.append(
            {
                "sample_id": scenario["sample_id"],
                "baseline_scenario_id": scenario["scenario_id"],
                "baseline_analysis_id": scenario["analysis_id"],
                "assignment_criteria": STRICT,
                "low_signal_mode": "include",
                "baseline_status": base["status"],
                "baseline_minimum": base["min_categories"],
                "min_delegated_classes": boundary["min_delegated_classes"],
                "transition_kinds": joined(sorted(kinds)),
                "delegations_found": len(boundary["changing_delegations"]),
            }
        )
    return delegation, details


def thresholds_by_sample(scenarios, analyses, samples):
    """How the result of each sample changes across detection thresholds."""
    rows = []
    for criteria in (STRICT, BROADER):
        for low in ("exclude", "include"):
            for sample in samples:
                group = [
                    row
                    for row in scenarios
                    if row["assignment_criteria"] == criteria
                    and row["low_signal_mode"] == low
                    and row["sample_id"] == sample
                ]
                nonempty = [row for row in group if row["detected_analyte_ids"]]
                empty = [row for row in group if not row["detected_analyte_ids"]]
                results = [analyses[row["analysis_id"]]["result"] for row in nonempty]
                evaluable = [
                    result for result in results if result["status"] == "evaluable"
                ]
                outcomes = {outcome(result) for result in results}
                minima = {result["min_categories"] for result in evaluable}
                rows.append(
                    {
                        "sample_id": sample,
                        "assignment_criteria": criteria,
                        "low_signal_mode": low,
                        "scenario_count": len(group),
                        "nonempty_scenario_count": len(nonempty),
                        "has_nonempty_and_empty": bool(nonempty and empty),
                        "always_empty": not nonempty,
                        "nonempty_outcome_changes": len(outcomes) > 1,
                        "nonempty_evaluability_changes": len(
                            {result["status"] for result in results}
                        )
                        > 1,
                        "evaluable_nonempty_minimum_changes": len(minima) > 1,
                        "evaluable_nonempty_result_changes": len(
                            {outcome(result) for result in evaluable}
                        )
                        > 1,
                        "evaluable_nonempty_minima": joined(sorted(minima)),
                        "scenario_ids": joined(
                            sorted(row["scenario_id"] for row in group)
                        ),
                        "nonempty_analysis_ids": joined(
                            sorted({row["analysis_id"] for row in nonempty})
                        ),
                        "nonempty_outcomes": [
                            json.loads(item) for item in sorted(outcomes)
                        ],
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
    parser.add_argument("--out", type=Path, required=True, help="Output folder.")
    args = parser.parse_args()

    scenarios = read_jsonl(args.results / "sample_scenarios.jsonl")
    analysis_list = read_jsonl(args.results / "detected_class_sets.jsonl")
    analyses = {row["analysis_id"]: row for row in analysis_list}
    samples = sorted({row["sample_id"] for row in scenarios})
    out = args.out

    removals, omissions, additions, summary = per_analysis_tables(analysis_list)
    write_csv(
        out / "TableS7_evidence_removal.csv",
        removals,
        (
            *COMMON,
            "changed_conclusions",
            "evidence_kind",
            "evidence_record_id",
            "removed_assignment_options",
        ),
    )
    write_csv(
        out / "TableS7_omissions.csv",
        omissions,
        (
            *COMMON,
            "analyte_id",
            "minimal_combinations_no_longer_minimal",
            "compatible_combinations_gained",
            "new_combinations",
        ),
    )
    write_csv(
        out / "TableS7_additions.csv",
        additions,
        (
            "analysis_id",
            "candidate_analyte_id",
            "status",
            "reason_code",
            "baseline_compatible_combinations",
            "compatible_combinations_after_addition",
            "eliminated_combinations",
            "reduction_bits",
        ),
    )
    write_csv(
        out / "delegation_summary.csv",
        summary,
        (
            "analysis_id",
            "conclusion_id",
            "baseline_value",
            "min_delegated_classes",
            "delegations_found",
        ),
    )
    delegation, details = delegation_by_sample(scenarios, analyses, samples)
    write_csv(
        out / "TableS7_delegation_by_sample.csv", delegation, tuple(delegation[0])
    )
    write_csv(out / "delegated_class_sets.csv", details, tuple(details[0]))
    thresholds = thresholds_by_sample(scenarios, analyses, samples)
    if sum(row["scenario_count"] for row in thresholds) != len(scenarios):
        raise ValueError("scenario count differs from the sum over samples")
    write_csv(
        out / "TableS7_thresholds_by_sample.csv", thresholds, tuple(thresholds[0])
    )

    kinds = Counter(row["transition_kinds"] for row in delegation)
    print(
        f"{len(removals)} evidence removals, {len(omissions)} omissions, {len(additions)} additions; "
        f"delegation by sample: {', '.join(f'{n} {k}' for k, n in sorted(kinds.items()))}."
    )


if __name__ == "__main__":
    main()
