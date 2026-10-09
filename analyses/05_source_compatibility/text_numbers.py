"""Print the numbers of the geological source-compatibility analysis quoted in the paper.

Reads the outputs of read_workbooks.py, source_tables.py,
compatible_combinations.py, sensitivity_tests.py and select_examples.py.
Unless stated otherwise, results are at a detection threshold of zero with
low-signal values included.

Usage:
  python analyses/05_source_compatibility/text_numbers.py --results results/sources
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biosynthesis.sensitivity import outcome_signature, zero_threshold  # noqa: E402
from biosynthesis.sources import BROADER, STRICT  # noqa: E402
from biosynthesis.tables import read_csv, read_jsonl  # noqa: E402
from select_examples import detected_at_zero, sample_labels, select  # noqa: E402

HOPANOID_CATEGORY = "hopanoid_producing_bacterial_sources"


def baseline(scenarios, criteria: str, low: str = "include") -> dict[str, dict]:
    """The zero-threshold scenario of each sample."""
    return {
        s["sample_id"]: s
        for s in scenarios
        if s["assignment_criteria"] == criteria
        and s["low_signal_mode"] == low
        and zero_threshold(s)
    }


def compute(results: Path) -> list[tuple[str, object, str]]:
    """Every number, as (name, value, sentence)."""
    s5 = read_csv(results / "TableS5_source_observations.csv")
    mappings = read_csv(results / "TableS6_analyte_precursor_map.csv")
    compatibilities = read_csv(results / "TableS6_source_compatibility.csv")
    links = read_csv(results / "TableS6_compatibility_links.csv")
    profiles = read_csv(results / "TableS6_source_profiles.csv")
    removals = read_csv(results / "TableS7_evidence_removal.csv")
    delegation = read_csv(results / "TableS7_delegation_by_sample.csv")
    thresholds = read_csv(results / "TableS7_thresholds_by_sample.csv")
    scenarios = read_jsonl(results / "sample_scenarios.jsonl")
    analyses = {
        a["analysis_id"]: a for a in read_jsonl(results / "detected_class_sets.jsonl")
    }
    samples = read_jsonl(results / "samples.jsonl")
    out: list[tuple[str, object, str]] = []

    def add(name, value, sentence):
        out.append((name, value, sentence.format(value=value)))

    def result(scenario):
        return analyses[scenario["analysis_id"]]["result"]

    # Observations and assignments.
    included = [r for r in s5 if r["included_sample"] == "true"]
    workbooks_by_sample: dict[str, set[str]] = {}
    for r in included:
        workbooks_by_sample.setdefault(r["canonical_sample_key"], set()).add(
            r["workbook_id"]
        )
    add("workbooks", len({r["workbook_id"] for r in s5}), "Workbooks read: {value}.")
    add(
        "samples",
        len(workbooks_by_sample),
        "Samples (named geological samples after joining the workbooks): {value}.",
    )
    add(
        "samples_one_workbook",
        sum(len(v) == 1 for v in workbooks_by_sample.values()),
        "Samples present in only one workbook: {value}.",
    )
    add(
        "low_signal_values",
        sum(r["low_signal"] == "true" for r in s5),
        "Values marked as low signals in Table S5: {value}.",
    )
    used = [r for r in mappings if r["included_in_analysis"] == "true"]
    add(
        "excluded_categories",
        sum(r["assignment_criteria"] == "excluded" for r in mappings),
        "Workbook source-interpretation categories excluded: {value}.",
    )
    add(
        "fossil_classes",
        len({r["analyte_id"] for r in used}),
        "Fossil classes used: {value}.",
    )
    add(
        "measurement_groups",
        len(
            {
                (r["measurement_channel"], r["denominator_group"])
                for r in included
                if r["value_used"] == "true"
            }
        ),
        "Measurement groups of the fossil classes used: {value}.",
    )
    add(
        "precursor_assignments",
        len(used),
        "Fossil-to-precursor assignments used: {value}.",
    )
    add(
        "precursor_classes",
        len({r["precursor_id"] for r in used}),
        "Precursor classes: {value}.",
    )
    cyclosterane = [r for r in used if r["analyte_id"] == "si2_cyclosteranes"]
    add(
        "cyclosterane_strict_assignments",
        sum(r["assignment_criteria"] == STRICT for r in cyclosterane),
        "Strict precursor assignments of cyclosteranes: {value} of "
        + str(len(cyclosterane))
        + ".",
    )
    add("source_links", len(compatibilities), "Precursor-to-source links: {value}.")
    add(
        "combinations_of_links",
        len(links),
        "Combinations of fossil class, precursor class and source category: {value}.",
    )
    categories = sorted({r["source_class_id"] for r in compatibilities})
    add("source_categories", len(categories), "Source categories: {value}.")
    add(
        "nonempty_category_combinations",
        2 ** len(categories) - 1,
        "Nonempty combinations of the categories: {value}.",
    )
    add(
        "source_profiles",
        len(profiles),
        "Source profiles recorded for context: {value}.",
    )
    targets = {
        r["evidence_record_id"] for r in removals if r["evidence_kind"] == "route"
    }
    add(
        "representative_targets",
        len(targets),
        "Representative reaction targets: {value}.",
    )
    for criteria in (BROADER, STRICT):
        records = Counter(
            kind
            for kind, _ in {
                (r["evidence_kind"], r["evidence_record_id"])
                for r in removals
                if r["assignment_criteria"] == criteria
            }
        )
        add(
            f"evidence_records_{criteria}",
            sum(records.values()),
            f"Evidence records tested by removal under {criteria} assignments: {{value}} "
            f"({records['mapping']} precursor assignments, {records['compatibility']} source links, "
            f"{records['route']} representative targets).",
        )
        for kind in ("mapping", "compatibility", "route"):
            out.append((f"evidence_records_{criteria}_{kind}", records[kind], ""))

    # Scenarios.
    statuses = Counter(result(s)["status"] for s in scenarios)
    add("scenarios", len(scenarios), "Sample scenarios: {value}.")
    add(
        "detected_sets",
        len({tuple(s["detected_analyte_ids"]) for s in scenarios}),
        "Distinct sets of detected classes (including none): {value}.",
    )
    for status in ("evaluable", "not_evaluable", "unconstrained"):
        add(
            f"scenarios_{status}",
            statuses[status],
            f"Scenarios {status.replace('_', ' ')}: {{value}}.",
        )
    settings = {(s["assignment_criteria"], s["low_signal_mode"]) for s in scenarios}
    add(
        "settings",
        len(settings),
        "Combinations of assignments and low-signal treatment: {value}.",
    )
    empty_per_setting = Counter(
        (s["sample_id"], s["assignment_criteria"], s["low_signal_mode"])
        for s in scenarios
        if not s["detected_analyte_ids"]
    )
    add(
        "one_empty_scenario_per_sample_and_setting",
        set(empty_per_setting.values()) == {1}
        and len(empty_per_setting) == len(workbooks_by_sample) * len(settings),
        "Each sample has exactly one scenario without detections under each setting: {value}.",
    )
    add(
        "all_five_compatible_when_evaluable",
        all(
            tuple(categories) in {tuple(c) for c in a["result"]["combinations"]}
            for a in analyses.values()
            if a["result"]["status"] == "evaluable"
        ),
        "The combination of all five categories is compatible in every evaluable analysis: {value}.",
    )

    # Baseline results.
    for criteria in (STRICT, BROADER):
        rows = [result(s) for s in baseline(scenarios, criteria).values()]
        add(
            f"{criteria}_samples",
            len(rows),
            f"{criteria.capitalize()} assignments, samples: {{value}}.",
        )
        add(
            f"{criteria}_not_evaluable",
            sum(r["status"] == "not_evaluable" for r in rows),
            f"{criteria.capitalize()}: not evaluable {{value}}.",
        )
        add(
            f"{criteria}_evaluable",
            sum(r["status"] == "evaluable" for r in rows),
            f"{criteria.capitalize()}: evaluable {{value}}.",
        )
        add(
            f"{criteria}_single_category",
            sum(r["min_categories"] == 1 for r in rows),
            f"{criteria.capitalize()}: a single category {{value}}.",
        )
        add(
            f"{criteria}_two_categories",
            sum(r["min_categories"] == 2 for r in rows),
            f"{criteria.capitalize()}: two categories {{value}}.",
        )

    # Examples.
    for row in select(detected_at_zero(scenarios), sample_labels(samples)):
        key = (
            "examples_hopanoids_cyclosteranes"
            if "cyclosteranes, no" in row["pattern"]
            else "examples_hopanoids_conventional"
        )
        add(
            key,
            row["samples_with_this_set"],
            f"Samples with {row['pattern']} sharing the most frequent set: {{value}} (next most frequent: "
            f"{row['next_most_frequent_set']}); example {row['sample_label']}.",
        )
        sample_id, label = row["sample_id"], row["sample_label"]
        records = next(m for m in samples if m["sample_id"] == sample_id)[
            "source_records"
        ]
        (age,) = {r["age_central_ma"] for r in records}
        add(
            f"{label}_age_ma",
            f"{age:,.0f}",
            f"{label}: about {{value}} Ma, "
            + ", ".join(
                sorted(
                    {
                        r["formation_or_basin_text"]
                        for r in records
                        if r["formation_or_basin_text"]
                    }
                    | {
                        r["geological_era_text"]
                        for r in records
                        if r["geological_era_text"]
                    }
                )
            )
            + ".",
        )
        for criteria in (STRICT, BROADER):
            r = result(baseline(scenarios, criteria)[sample_id])
            minimal = r["minimal_combinations"]
            add(
                f"{label}_{criteria}_status",
                r["status"],
                f"{label}, {criteria}: {{value}}.",
            )
            add(
                f"{label}_{criteria}_minimal",
                len(minimal),
                f"{label}, {criteria}: minimal combinations {{value}}.",
            )
            add(
                f"{label}_{criteria}_compatible",
                r["compatible_combinations"],
                f"{label}, {criteria}: compatible combinations {{value}}.",
            )
            add(
                f"{label}_{criteria}_minimum",
                r["min_categories"],
                f"{label}, {criteria}: minimum number of categories {{value}}.",
            )
            add(
                f"{label}_{criteria}_minimal_with_hopanoid_bacteria",
                sum(HOPANOID_CATEGORY in c for c in minimal),
                f"{label}, {criteria}: minimal combinations including hopanoid-producing bacteria {{value}}.",
            )
            add(
                f"{label}_{criteria}_minimal_sizes",
                sorted({len(c) for c in minimal}),
                f"{label}, {criteria}: sizes of the minimal combinations {{value}}.",
            )

    # Sensitivity tests (threshold zero, low signals included).
    for criteria in (STRICT, BROADER):
        base = baseline(scenarios, criteria)
        evaluable = [s for s in base.values() if result(s)["status"] == "evaluable"]
        changed = [
            s
            for s in base.values()
            if any(
                r["changed_conclusions"]
                for r in analyses[s["analysis_id"]]["evidence_removal_test"]["removals"]
            )
        ]
        add(
            f"{criteria}_evaluable_with_essential_record",
            len(changed),
            f"{criteria.capitalize()}: samples with an evidence record whose removal changes the result {{value}} "
            f"(evaluable samples: {len(evaluable)}).",
        )
    kinds = Counter(r["transition_kinds"] for r in delegation)
    add(
        "delegation_became_evaluable",
        kinds["became_evaluable"],
        "Delegation, first change: became evaluable {value}.",
    )
    add(
        "delegation_became_unconstrained",
        kinds["became_empty_unconstrained"],
        "Delegation, first change: no detected class left {value}.",
    )
    add(
        "delegation_minimum_reduced",
        sum(
            r["transition_kinds"] == "numerical_reduction_within_evaluable"
            and r["baseline_minimum"] == "2"
            for r in delegation
        ),
        "Delegation, first change: minimum reduced from two to one {value}.",
    )
    base_strict = baseline(scenarios, STRICT)
    emptied = [
        r for r in delegation if r["transition_kinds"] == "became_empty_unconstrained"
    ]
    add(
        "delegation_unconstrained_classes",
        sorted(
            {len(base_strict[r["sample_id"]]["detected_analyte_ids"]) for r in emptied}
        ),
        "Detected classes of the samples left without a detected class: {value}.",
    )
    add(
        "delegation_unconstrained_were_single_category",
        all(
            r["baseline_minimum"] == "1"
            and int(r["min_delegated_classes"])
            == len(base_strict[r["sample_id"]]["detected_analyte_ids"])
            for r in emptied
        ),
        "Those samples had a single category and all their classes were delegated: {value}.",
    )
    allowance = Counter(r["min_delegated_classes"] for r in delegation)
    add(
        "delegation_one_class",
        allowance["1"],
        "Delegating one class changed the result for {value} samples.",
    )
    add(
        "delegation_two_classes",
        allowance["2"],
        "Two classes were needed for {value} samples.",
    )

    flags = [
        r
        for r in thresholds
        if r["assignment_criteria"] == STRICT and r["low_signal_mode"] == "exclude"
    ]
    for flag, name in (
        ("has_nonempty_and_empty", "threshold_nonempty_and_empty"),
        ("always_empty", "threshold_always_empty"),
        ("nonempty_outcome_changes", "threshold_result_changes"),
        ("nonempty_evaluability_changes", "threshold_evaluability_changes"),
        ("evaluable_nonempty_minimum_changes", "threshold_minimum_changes"),
    ):
        add(
            name,
            sum(r[flag] == "true" for r in flags),
            f"Strict, low signals excluded, across thresholds: {flag.replace('_', ' ')} {{value}}.",
        )

    strict, broader = baseline(scenarios, STRICT), baseline(scenarios, BROADER)
    changed = [
        k
        for k in strict
        if outcome_signature(result(strict[k])) != outcome_signature(result(broader[k]))
    ]
    add(
        "criteria_changed",
        len(changed),
        "Strict to broader assignments changed the result for {value} samples.",
    )
    add(
        "criteria_became_evaluable",
        sum(
            result(strict[k])["status"] != "evaluable"
            and result(broader[k])["status"] == "evaluable"
            for k in changed
        ),
        "Of these, became evaluable: {value}.",
    )
    add(
        "criteria_same_minimum",
        sum(
            result(strict[k])["status"] == result(broader[k])["status"] == "evaluable"
            and result(strict[k])["min_categories"]
            == result(broader[k])["min_categories"]
            for k in changed
        ),
        "Of these, evaluable under both with the same minimum: {value}.",
    )
    for criteria in (STRICT, BROADER):
        included_, excluded_ = (
            baseline(scenarios, criteria, "include"),
            baseline(scenarios, criteria, "exclude"),
        )
        add(
            f"low_signal_changed_{criteria}",
            sum(
                outcome_signature(result(included_[k]))
                != outcome_signature(result(excluded_[k]))
                for k in included_
            ),
            f"Low-signal treatment changed the result ({criteria}) for {{value}} samples.",
        )
    expanded = reduced = 0
    for s in strict.values():
        analysis = analyses[s["analysis_id"]]
        expanded += any(
            (r["compatible_combinations_gained"] or 0) > 0
            for r in analysis["omission_test"]["omissions"]
        )
        test = analysis["addition_test"]
        reduced += bool(test) and any(
            r["status"] == "exact"
            and r["reduction_bits"] is not None
            and r["reduction_bits"] > 0
            for r in test["additions"]
        )
    add(
        "omission_expanded",
        expanded,
        "Strict: omitting one detected class enlarged the compatible combinations for {value} samples.",
    )
    add(
        "addition_reduced",
        reduced,
        "Strict: detecting an undetected class would reduce them for {value} samples.",
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--results",
        type=Path,
        required=True,
        help="Folder of the source-compatibility outputs.",
    )
    args = parser.parse_args()
    for _, _, sentence in compute(args.results):
        if sentence:
            print(sentence)


if __name__ == "__main__":
    main()
