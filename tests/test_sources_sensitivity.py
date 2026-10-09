"""Threshold sweep and sensitivity tests on a small example."""

import math

from biosynthesis.sensitivity import (
    Observation,
    addition_test,
    delegation_test,
    evidence_removal_test,
    global_thresholds,
    omission_test,
    presence_panels,
    transition_kind,
)
from biosynthesis.sources import AssignmentOption

CATEGORIES = ("bacteria", "eukaryotes", "hypothesis")
OPTIONS = (
    AssignmentOption("hopane", "hopanoids", "bacteria", "m_hop", "c_hop_bac", "r_hop"),
    AssignmentOption("sterane", "sterols", "eukaryotes", "m_ste", "c_ste_euk", "r_ste"),
    AssignmentOption("sterane", "sterols", "bacteria", "m_ste", "c_ste_bac", "r_ste"),
    AssignmentOption(
        "cyclosterane", "cyclo", "hypothesis", "m_cyc", "c_cyc_hyp", "r_cyc"
    ),
)


def observation(sample, analyte, value, low=False):
    return Observation(sample, analyte, "channel", "sum", value, low)


def test_thresholds_give_one_scenario_per_set_of_detected_classes():
    rows = [
        observation("s1", "hopane", 0.6),
        observation("s1", "sterane", 0.3, low=True),
        observation("s2", "hopane", 0.1),
    ]
    thresholds = global_thresholds(rows)
    assert thresholds == {"channel|sum": (0.0, 0.1, 0.3, 0.6)}
    included = presence_panels(
        sample_id="s1",
        observations=rows,
        thresholds=thresholds,
        include_low_signal=True,
    )
    assert [p.detected_analyte_ids for p in included] == [
        (),
        ("hopane",),
        ("hopane", "sterane"),
    ]
    first = included[2].threshold_region_sets[0][0]
    assert (first.threshold_lower_inclusive, first.threshold_upper_exclusive) == (
        0.0,
        0.3,
    )
    excluded = presence_panels(
        sample_id="s1",
        observations=rows,
        thresholds=thresholds,
        include_low_signal=False,
    )
    assert [p.detected_analyte_ids for p in excluded] == [(), ("hopane",)]


def arguments(detected):
    return dict(
        detected_analyte_ids=detected,
        assignment_options=OPTIONS,
        source_categories=CATEGORIES,
    )


def test_delegation_finds_the_smallest_change_of_the_minimum():
    test = delegation_test(**arguments(["hopane", "cyclosterane"]))
    (boundary,) = [b for b in test.boundaries if b.conclusion_id == "min_categories"]
    assert boundary.baseline_value == 2
    assert boundary.min_delegated_classes == 1
    assert sorted(d.delegated_analyte_ids for d in boundary.changing_delegations) == [
        ("cyclosterane",),
        ("hopane",),
    ]
    assert len(test.delegations) == 4


def test_removing_a_record_removes_every_assignment_that_uses_it():
    test = evidence_removal_test(**arguments(["sterane", "hopane"]))
    kinds = [(r.evidence_kind, r.evidence_record_id) for r in test.removals]
    assert ("route", "r_ste") in kinds and ("mapping", "m_hop") in kinds
    route = next(r for r in test.removals if r.evidence_record_id == "r_ste")
    assert len(route.removed_assignment_options) == 2
    assert route.summary.status == "not_evaluable"
    # Bacteria alone account for both classes, so losing the eukaryotic link changes nothing.
    eukaryotes = next(r for r in test.removals if r.evidence_record_id == "c_ste_euk")
    assert eukaryotes.changed_conclusions == ()
    bacteria = next(r for r in test.removals if r.evidence_record_id == "c_ste_bac")
    assert bacteria.summary.compatible_combinations == 2
    assert "min_categories" in bacteria.changed_conclusions


def test_omission_can_only_add_and_addition_can_only_remove_combinations():
    omitted = omission_test(**arguments(["sterane", "hopane"]))
    by_class = {r.analyte_id: r for r in omitted.omissions}
    assert by_class["hopane"].compatible_combinations_gained == 2
    added = addition_test(
        candidate_analyte_ids=["cyclosterane"], **arguments(["sterane", "hopane"])
    )
    (result,) = added.additions
    assert (
        result.baseline_compatible_combinations,
        result.compatible_combinations_after_addition,
    ) == (4, 2)
    assert result.reduction_bits == math.log2(4) - math.log2(2)


def test_transition_kinds():
    assert transition_kind("not_evaluable", None, "evaluable", 1) == "became_evaluable"
    assert (
        transition_kind("evaluable", 1, "unconstrained", None)
        == "became_empty_unconstrained"
    )
    assert (
        transition_kind("evaluable", 2, "evaluable", 1)
        == "numerical_reduction_within_evaluable"
    )
