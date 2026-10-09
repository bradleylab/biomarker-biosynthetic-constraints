"""Detection-threshold sweep and the sensitivity tests of the source-compatibility analysis.

The sweep turns each sample's reported abundances into scenarios: one per set
of detected fossil classes, for each low-signal treatment and set of
assignments. The tests recompute the compatible combinations after delegating
detected classes to sources outside the five categories, removing one evidence
record, omitting one detected class, or adding one undetected class.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from itertools import combinations
from typing import Any, Iterable, Literal, Mapping, Sequence, TypeAlias

from biosynthesis.sources import (
    BROADER,
    STRICT,
    AssignmentOption,
    CombinationResult,
    compatible_combinations,
)

ConclusionValue: TypeAlias = (
    str | int | bool | None | tuple[str, ...] | tuple[tuple[str, ...], ...]
)
ESSENTIAL_PREFIX = "in_every_minimal:"


# --- Detection-threshold sweep -----------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Observation:
    """One numeric value of a fossil class in a sample."""

    sample_id: str
    analyte_id: str
    measurement_channel: str
    denominator_group: str
    reported_abundance: float
    low_signal: bool

    @property
    def group_id(self) -> str:
        return f"{self.measurement_channel}|{self.denominator_group}"


@dataclass(frozen=True, slots=True)
class GroupThresholdState:
    """The classes detected in one measurement group over one interval of thresholds."""

    group_id: str
    threshold_lower_inclusive: float
    threshold_upper_exclusive: float | None
    detected_analyte_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PresencePanel:
    """One set of detected classes of a sample and the threshold intervals that give it."""

    detected_analyte_ids: tuple[str, ...]
    threshold_region_sets: tuple[tuple[GroupThresholdState, ...], ...]


@dataclass(frozen=True, slots=True)
class CriteriaAssignment:
    """An assignment option and the strictest criteria under which it is used."""

    option: AssignmentOption
    criteria: str


@dataclass(frozen=True, slots=True)
class RouteState:
    """The representative reaction target of one precursor class and its route state."""

    route_record_id: str
    precursor_id: str
    state: str
    reason_code: str


def _sorting_text(value: object) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def global_thresholds(
    observations: Sequence[Observation],
) -> dict[str, tuple[float, ...]]:
    """Zero and every distinct positive abundance of each measurement group."""
    values: dict[str, set[float]] = defaultdict(lambda: {0.0})
    for row in observations:
        if row.reported_abundance > 0:
            values[row.group_id].add(row.reported_abundance)
    return {
        group: tuple(sorted(thresholds)) for group, thresholds in sorted(values.items())
    }


def _group_threshold_states(
    *,
    sample_rows: Sequence[Observation],
    group_id: str,
    thresholds: Sequence[float],
    include_low_signal: bool,
) -> tuple[GroupThresholdState, ...]:
    group_rows = tuple(row for row in sample_rows if row.group_id == group_id)

    def detected(threshold: float) -> tuple[str, ...]:
        return tuple(
            sorted(
                {
                    row.analyte_id
                    for row in group_rows
                    if row.reported_abundance > threshold
                    and (include_low_signal or not row.low_signal)
                }
            )
        )

    states: list[GroupThresholdState] = []
    start = float(thresholds[0])
    current = detected(start)
    for threshold in thresholds[1:]:
        next_panel = detected(float(threshold))
        if next_panel == current:
            continue
        states.append(
            GroupThresholdState(
                group_id=group_id,
                threshold_lower_inclusive=start,
                threshold_upper_exclusive=float(threshold),
                detected_analyte_ids=current,
            )
        )
        start = float(threshold)
        current = next_panel
    states.append(
        GroupThresholdState(
            group_id=group_id,
            threshold_lower_inclusive=start,
            threshold_upper_exclusive=None,
            detected_analyte_ids=current,
        )
    )
    return tuple(states)


def presence_panels(
    *,
    sample_id: str,
    observations: Sequence[Observation],
    thresholds: Mapping[str, Sequence[float]],
    include_low_signal: bool,
) -> tuple[PresencePanel, ...]:
    """The distinct sets of detected classes of one sample over all combinations of group thresholds."""
    sample_rows = tuple(row for row in observations if row.sample_id == sample_id)
    if not sample_rows:
        return ()
    per_group = tuple(
        _group_threshold_states(
            sample_rows=sample_rows,
            group_id=group_id,
            thresholds=group_thresholds,
            include_low_signal=include_low_signal,
        )
        for group_id, group_thresholds in sorted(thresholds.items())
    )
    regions_by_panel: dict[tuple[str, ...], list[tuple[GroupThresholdState, ...]]] = (
        defaultdict(list)
    )
    for region in itertools.product(*per_group):
        detected = tuple(
            sorted(
                {analyte for state in region for analyte in state.detected_analyte_ids}
            )
        )
        regions_by_panel[detected].append(tuple(region))
    return tuple(
        PresencePanel(
            detected_analyte_ids=detected,
            threshold_region_sets=tuple(
                sorted(
                    regions,
                    key=lambda region: _sorting_text([asdict(item) for item in region]),
                )
            ),
        )
        for detected, regions in sorted(regions_by_panel.items())
    )


def zero_threshold(scenario: Mapping[str, Any]) -> bool:
    """True for the scenario that holds a threshold of zero in every measurement group."""
    return any(
        region and all(state["threshold_lower_inclusive"] == 0 for state in region)
        for region in scenario["threshold_regions"]
    )


# --- Sensitivity tests -------------------------------------------------------------------


@dataclass(frozen=True)
class Delegation:
    """The result after delegating a given subset of detected classes."""

    delegated_count: int
    delegated_analyte_ids: tuple[str, ...]
    summary: CombinationResult


@dataclass(frozen=True)
class ChangingDelegation:
    """A smallest delegation that changes one conclusion, with the changed value."""

    delegated_analyte_ids: tuple[str, ...]
    changed_value: ConclusionValue
    summary: CombinationResult


@dataclass(frozen=True)
class Boundary:
    """The smallest number of delegated classes that changes one conclusion."""

    conclusion_id: str
    baseline_value: ConclusionValue
    min_delegated_classes: int
    changing_delegations: tuple[ChangingDelegation, ...]


@dataclass(frozen=True)
class DelegationTest:
    baseline: CombinationResult
    delegations: tuple[Delegation, ...]
    boundaries: tuple[Boundary, ...]


@dataclass(frozen=True)
class EvidenceRemoval:
    """The result after removing one evidence record and every assignment that uses it."""

    evidence_kind: Literal["mapping", "compatibility", "route"]
    evidence_record_id: str
    removed_assignment_options: tuple[AssignmentOption, ...]
    summary: CombinationResult
    changed_conclusions: tuple[str, ...]


@dataclass(frozen=True)
class EvidenceRemovalTest:
    baseline: CombinationResult
    removals: tuple[EvidenceRemoval, ...]


@dataclass(frozen=True)
class Addition:
    """The result if one currently undetected class were detected."""

    candidate_analyte_id: str
    status: Literal["exact", "not_evaluable"]
    reason_code: str | None
    baseline_compatible_combinations: int | None
    compatible_combinations_after_addition: int | None
    eliminated_combinations: tuple[tuple[str, ...], ...]
    remaining_combinations: tuple[tuple[str, ...], ...]
    reduction_bits: float | None
    summary: CombinationResult


@dataclass(frozen=True)
class AdditionTest:
    baseline: CombinationResult
    additions: tuple[Addition, ...]


@dataclass(frozen=True)
class Omission:
    """The result after omitting one detected class."""

    analyte_id: str
    summary: CombinationResult
    compatible_combinations_gained: int | None
    new_combinations: tuple[tuple[str, ...], ...]
    minimal_combinations_no_longer_minimal: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class OmissionTest:
    baseline: CombinationResult
    omissions: tuple[Omission, ...]


def delegation_test(
    *,
    detected_analyte_ids: Sequence[str],
    assignment_options: Sequence[AssignmentOption],
    source_categories: Sequence[str],
) -> DelegationTest:
    """Delegate every subset of the detected classes and record the smallest change of each conclusion."""
    baseline = compatible_combinations(
        detected_analyte_ids=detected_analyte_ids,
        assignment_options=assignment_options,
        source_categories=source_categories,
    )
    analytes = baseline.detected_analyte_ids
    delegations = tuple(
        Delegation(
            delegated_count=count,
            delegated_analyte_ids=delegated,
            summary=compatible_combinations(
                detected_analyte_ids=_without(analytes, delegated),
                assignment_options=assignment_options,
                source_categories=baseline.source_categories,
            ),
        )
        for count in range(len(analytes) + 1)
        for delegated in combinations(analytes, count)
    )
    baseline_conclusions = conclusions(baseline)
    delegation_conclusions = tuple(
        (item, conclusions(item.summary)) for item in delegations
    )
    boundaries: list[Boundary] = []
    for conclusion_id, baseline_value in sorted(baseline_conclusions.items()):
        changed = [
            (item, values)
            for item, values in delegation_conclusions
            if values[conclusion_id] != baseline_value
        ]
        if not changed:
            continue
        smallest = min(item.delegated_count for item, _ in changed)
        boundaries.append(
            Boundary(
                conclusion_id=conclusion_id,
                baseline_value=baseline_value,
                min_delegated_classes=smallest,
                changing_delegations=tuple(
                    ChangingDelegation(
                        delegated_analyte_ids=item.delegated_analyte_ids,
                        changed_value=values[conclusion_id],
                        summary=item.summary,
                    )
                    for item, values in changed
                    if item.delegated_count == smallest
                ),
            )
        )
    return DelegationTest(
        baseline=baseline, delegations=delegations, boundaries=tuple(boundaries)
    )


def evidence_removal_test(
    *,
    detected_analyte_ids: Sequence[str],
    assignment_options: Sequence[AssignmentOption],
    source_categories: Sequence[str],
) -> EvidenceRemovalTest:
    """Remove each fossil-to-precursor link, precursor-to-source link and representative target in turn."""
    baseline = compatible_combinations(
        detected_analyte_ids=detected_analyte_ids,
        assignment_options=assignment_options,
        source_categories=source_categories,
    )
    options = tuple(sorted(set(assignment_options)))
    relevant = tuple(
        option
        for option in options
        if option.analyte_id in baseline.detected_analyte_ids
    )
    specs: list[tuple[Literal["mapping", "compatibility", "route"], str, str]] = []
    for kind, field in (
        ("mapping", "mapping_record_id"),
        ("compatibility", "compatibility_record_id"),
        ("route", "route_record_id"),
    ):
        specs.extend(
            (kind, record_id, field)
            for record_id in sorted(
                {
                    getattr(option, field)
                    for option in relevant
                    if getattr(option, field)
                }
            )
        )
    baseline_conclusions = conclusions(baseline)
    removals: list[EvidenceRemoval] = []
    for kind, record_id, field in specs:
        removed = tuple(
            option for option in relevant if getattr(option, field) == record_id
        )
        retained = tuple(
            option for option in options if getattr(option, field) != record_id
        )
        summary = compatible_combinations(
            detected_analyte_ids=baseline.detected_analyte_ids,
            assignment_options=retained,
            source_categories=baseline.source_categories,
        )
        values = conclusions(summary)
        removals.append(
            EvidenceRemoval(
                evidence_kind=kind,
                evidence_record_id=record_id,
                removed_assignment_options=removed,
                summary=summary,
                changed_conclusions=tuple(
                    key
                    for key in sorted(baseline_conclusions)
                    if values[key] != baseline_conclusions[key]
                ),
            )
        )
    return EvidenceRemovalTest(baseline=baseline, removals=tuple(removals))


def addition_test(
    *,
    detected_analyte_ids: Sequence[str],
    candidate_analyte_ids: Sequence[str],
    assignment_options: Sequence[AssignmentOption],
    source_categories: Sequence[str],
) -> AdditionTest:
    """Recompute the combinations as if each undetected class were detected."""
    baseline = compatible_combinations(
        detected_analyte_ids=detected_analyte_ids,
        assignment_options=assignment_options,
        source_categories=source_categories,
    )
    candidates = tuple(sorted(set(candidate_analyte_ids)))
    if not all(candidates) or set(candidates) & set(baseline.detected_analyte_ids):
        raise ValueError("added classes must be named and not already detected")
    return AdditionTest(
        baseline=baseline,
        additions=tuple(
            _addition(
                baseline=baseline,
                candidate_analyte_id=candidate,
                assignment_options=assignment_options,
            )
            for candidate in candidates
        ),
    )


def omission_test(
    *,
    detected_analyte_ids: Sequence[str],
    assignment_options: Sequence[AssignmentOption],
    source_categories: Sequence[str],
) -> OmissionTest:
    """Recompute the combinations after omitting each detected class in turn."""
    baseline = compatible_combinations(
        detected_analyte_ids=detected_analyte_ids,
        assignment_options=assignment_options,
        source_categories=source_categories,
    )
    baseline_feasible = set(baseline.combinations)
    baseline_minimal = set(baseline.minimal_combinations)
    omissions: list[Omission] = []
    for analyte_id in baseline.detected_analyte_ids:
        summary = compatible_combinations(
            detected_analyte_ids=_without(baseline.detected_analyte_ids, (analyte_id,)),
            assignment_options=assignment_options,
            source_categories=baseline.source_categories,
        )
        if baseline.status == "evaluable" and summary.status == "evaluable":
            summary_feasible = set(summary.combinations)
            gained = len(summary_feasible) - len(baseline_feasible)
            new = tuple(sorted(summary_feasible - baseline_feasible, key=_cover_key))
            no_longer_minimal = tuple(
                sorted(
                    baseline_minimal - set(summary.minimal_combinations), key=_cover_key
                )
            )
        else:
            gained = None
            new = ()
            no_longer_minimal = ()
        omissions.append(
            Omission(
                analyte_id=analyte_id,
                summary=summary,
                compatible_combinations_gained=gained,
                new_combinations=new,
                minimal_combinations_no_longer_minimal=no_longer_minimal,
            )
        )
    return OmissionTest(baseline=baseline, omissions=tuple(omissions))


def _addition(
    *,
    baseline: CombinationResult,
    candidate_analyte_id: str,
    assignment_options: Sequence[AssignmentOption],
) -> Addition:
    summary = compatible_combinations(
        detected_analyte_ids=(*baseline.detected_analyte_ids, candidate_analyte_id),
        assignment_options=assignment_options,
        source_categories=baseline.source_categories,
    )
    if baseline.status != "evaluable" or summary.status != "evaluable":
        return Addition(
            candidate_analyte_id=candidate_analyte_id,
            status="not_evaluable",
            reason_code="not_evaluable_after_addition",
            baseline_compatible_combinations=baseline.compatible_combinations,
            compatible_combinations_after_addition=summary.compatible_combinations,
            eliminated_combinations=(),
            remaining_combinations=(),
            reduction_bits=None,
            summary=summary,
        )
    baseline_covers = set(baseline.combinations)
    conditioned_covers = set(summary.combinations)
    if not conditioned_covers <= baseline_covers:
        raise ValueError("adding a detected class added compatible combinations")
    return Addition(
        candidate_analyte_id=candidate_analyte_id,
        status="exact",
        reason_code=None,
        baseline_compatible_combinations=baseline.compatible_combinations,
        compatible_combinations_after_addition=summary.compatible_combinations,
        eliminated_combinations=tuple(
            sorted(baseline_covers - conditioned_covers, key=_cover_key)
        ),
        remaining_combinations=tuple(sorted(conditioned_covers, key=_cover_key)),
        reduction_bits=math.log2(len(baseline_covers))
        - math.log2(len(conditioned_covers)),
        summary=summary,
    )


def conclusions(summary: CombinationResult) -> dict[str, ConclusionValue]:
    """The conclusions whose change the sensitivity tests record."""
    evaluable = summary.status == "evaluable"
    values: dict[str, ConclusionValue] = {
        "status": summary.status,
        "compatible_combinations": summary.compatible_combinations,
        "single_category_combinations": summary.single_category_combinations,
        "min_categories": summary.min_categories,
        "minimal_combinations": summary.minimal_combinations,
        "categories_in_every_minimal": (
            summary.categories_in_every_minimal if evaluable else None
        ),
    }
    for category in summary.source_categories:
        values[f"{ESSENTIAL_PREFIX}{category}"] = (
            category in summary.categories_in_every_minimal if evaluable else None
        )
    return values


def _cover_key(cover: tuple[str, ...]) -> tuple[int, tuple[str, ...]]:
    return len(cover), cover


def _without(values: tuple[str, ...], removed: tuple[str, ...]) -> tuple[str, ...]:
    removed_set = set(removed)
    return tuple(value for value in values if value not in removed_set)


# --- Scenarios and analyses --------------------------------------------------------------


def criteria_options(
    assignments: Sequence[CriteriaAssignment], criteria: str
) -> tuple[AssignmentOption, ...]:
    """The assignment options used under strict or broader criteria."""
    return tuple(
        row.option
        for row in assignments
        if row.criteria == STRICT or (criteria == BROADER and row.criteria == BROADER)
    )


def route_context(
    detected_analytes: Sequence[str],
    options: Sequence[AssignmentOption],
    route_states: Sequence[RouteState],
) -> tuple[dict[str, object], ...]:
    """The representative targets of the precursor classes reachable from the detected classes."""
    precursors = {
        option.precursor_id
        for option in options
        if option.analyte_id in set(detected_analytes)
    }
    return tuple(
        asdict(state) for state in route_states if state.precursor_id in precursors
    )


def analyze_detected_classes(
    *,
    detected: tuple[str, ...],
    criteria: str,
    assignments: Sequence[CriteriaAssignment],
    source_categories: Sequence[str],
    all_analyte_ids: Sequence[str],
    route_states: Sequence[RouteState],
) -> dict[str, object]:
    """The combinations and all sensitivity tests for one set of detected classes and criteria."""
    options = criteria_options(assignments, criteria)
    arguments = dict(
        detected_analyte_ids=detected,
        assignment_options=options,
        source_categories=source_categories,
    )
    candidates = tuple(sorted(set(all_analyte_ids) - set(detected)))
    addition = (
        addition_test(candidate_analyte_ids=candidates, **arguments)
        if candidates
        else None
    )
    return {
        "assignment_criteria": criteria,
        "detected_analyte_ids": list(detected),
        "result": asdict(compatible_combinations(**arguments)),
        "delegation_test": asdict(delegation_test(**arguments)),
        "evidence_removal_test": asdict(evidence_removal_test(**arguments)),
        "omission_test": asdict(omission_test(**arguments)),
        "addition_test": asdict(addition) if addition else None,
        "route_context": route_context(detected, options, route_states),
    }


def scenarios_and_analyses(
    *,
    sample_ids: Sequence[str],
    observations: Sequence[Observation],
    assignments: Sequence[CriteriaAssignment],
    source_categories: Sequence[str],
    all_analyte_ids: Sequence[str],
    route_states: Sequence[RouteState],
    low_signal_treatments: Iterable[str],
    criteria_list: Iterable[str],
) -> tuple[
    list[dict[str, object]], dict[tuple[str, tuple[str, ...]], dict[str, object]]
]:
    """Sweep thresholds, low-signal treatments and criteria for every sample; analyze each distinct set once."""
    thresholds = global_thresholds(observations)
    scenarios: list[dict[str, object]] = []
    analyses: dict[tuple[str, tuple[str, ...]], dict[str, object]] = {}
    for sample_id in sample_ids:
        has_rows = any(row.sample_id == sample_id for row in observations)
        for low_signal_mode in sorted(low_signal_treatments):
            for criteria in sorted(criteria_list):
                if not has_rows:
                    scenarios.append(
                        {
                            "sample_id": sample_id,
                            "presence_profile": None,
                            "threshold_regions": [],
                            "low_signal_mode": low_signal_mode,
                            "assignment_criteria": criteria,
                            "detected_analyte_ids": [],
                        }
                    )
                    continue
                for panel in presence_panels(
                    sample_id=sample_id,
                    observations=observations,
                    thresholds=thresholds,
                    include_low_signal=low_signal_mode == "include",
                ):
                    key = (criteria, panel.detected_analyte_ids)
                    if key not in analyses:
                        analyses[key] = analyze_detected_classes(
                            detected=panel.detected_analyte_ids,
                            criteria=criteria,
                            assignments=assignments,
                            source_categories=source_categories,
                            all_analyte_ids=all_analyte_ids,
                            route_states=route_states,
                        )
                    scenarios.append(
                        {
                            "sample_id": sample_id,
                            "presence_profile": panel.detected_analyte_ids,
                            "threshold_regions": [
                                [asdict(item) for item in region]
                                for region in panel.threshold_region_sets
                            ],
                            "low_signal_mode": low_signal_mode,
                            "assignment_criteria": criteria,
                            "detected_analyte_ids": list(panel.detected_analyte_ids),
                        }
                    )
    return scenarios, analyses


def assign_identifiers(
    scenarios: list[dict[str, object]],
    analyses: dict[tuple[str, tuple[str, ...]], dict[str, object]],
    *,
    criteria_order: Sequence[str] = (STRICT, BROADER),
    low_signal_order: Sequence[str] = ("include", "exclude"),
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Number detected-class sets (P01...), analyses (A01...) and scenarios (S0001...).

    Sets of detected classes are numbered by the number of classes and then by
    their sorted identifiers; analyses by criteria (strict first) and then set;
    scenarios by sample, criteria, low-signal treatment (included first) and set.
    """
    profiles = sorted(
        {key[1] for key in analyses}, key=lambda detected: (len(detected), detected)
    )
    profile_number = {
        detected: index for index, detected in enumerate(profiles, start=1)
    }
    profile_ids = {
        detected: f"P{number:02d}" for detected, number in profile_number.items()
    }
    analysis_keys = sorted(
        analyses, key=lambda key: (criteria_order.index(key[0]), profile_number[key[1]])
    )
    analysis_ids = {
        key: f"A{index:02d}" for index, key in enumerate(analysis_keys, start=1)
    }
    analysis_rows = [
        {
            "analysis_id": analysis_ids[key],
            "presence_profile_id": profile_ids[key[1]],
            **analyses[key],
        }
        for key in analysis_keys
    ]

    def scenario_key(row: Mapping[str, object]) -> tuple:
        profile = row["presence_profile"]
        return (
            row["sample_id"],
            criteria_order.index(str(row["assignment_criteria"])),
            low_signal_order.index(str(row["low_signal_mode"])),
            0 if profile is None else profile_number[profile],
        )

    scenario_rows = []
    for index, row in enumerate(sorted(scenarios, key=scenario_key), start=1):
        profile = row["presence_profile"]
        scenario_rows.append(
            {
                "scenario_id": f"S{index:04d}",
                "sample_id": row["sample_id"],
                "low_signal_mode": row["low_signal_mode"],
                "assignment_criteria": row["assignment_criteria"],
                "presence_profile_id": None
                if profile is None
                else profile_ids[profile],
                "threshold_regions": row["threshold_regions"],
                "detected_analyte_ids": row["detected_analyte_ids"],
                "analysis_id": (
                    None
                    if profile is None
                    else analysis_ids[(str(row["assignment_criteria"]), profile)]
                ),
            }
        )
    return scenario_rows, analysis_rows


# --- Results by sample -------------------------------------------------------------------


def transition_kind(before: str, old: int | None, after: str, new: int | None) -> str:
    """Name the change of a result: evaluability, an empty set of classes, or the minimum number of categories."""
    for state, value in ((before, old), (after, new)):
        if (state == "evaluable") != (isinstance(value, int) and value >= 1):
            raise ValueError(f"status and minimum disagree: {state}/{value}")
    if after == "unconstrained" and before != after:
        return "became_empty_unconstrained"
    if before != "evaluable" and after == "evaluable":
        return "became_evaluable"
    if before == "evaluable" and after == "not_evaluable":
        return "became_not_evaluable"
    if before == after == "evaluable" and new != old:
        return (
            "numerical_reduction_within_evaluable"
            if new < old
            else "numerical_increase_within_evaluable"
        )
    return "no_minimum_change"


def outcome(result: Mapping[str, Any]) -> str:
    """Status, minimum, single-category combinations and minimal combinations, as JSON."""
    return json.dumps(
        [
            result["status"],
            result["min_categories"],
            result["single_category_combinations"],
            result["minimal_combinations"],
        ],
        separators=(",", ":"),
    )


def outcome_signature(result: Mapping[str, Any]) -> tuple[object, ...]:
    """The same outcome as a tuple, for comparing two results."""
    return (
        result["status"],
        result["min_categories"],
        tuple(result["single_category_combinations"]),
        tuple(tuple(row) for row in result["minimal_combinations"]),
    )
