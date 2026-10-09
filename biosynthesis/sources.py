"""Compatible combinations of source categories for a set of detected fossil classes.

A combination of source categories is compatible when every detected fossil
class has at least one assignment (fossil class to precursor class to source
category) whose category is in the combination. The number of compatible
combinations counts sets of source categories.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, Literal, Sequence

import yaml

STRICT = "strict"
BROADER = "broader"
USED = "used"


@dataclass(frozen=True)
class Link:
    """One fossil-to-precursor or precursor-to-source link of the assignment tables."""

    record_id: str
    disposition: str
    criteria: str
    precursor_id: str
    included: bool
    analyte_id: str = ""
    source_class_id: str = ""


@dataclass(frozen=True, order=True)
class AssignmentOption:
    """One assignment of a fossil class to a source category through a precursor class."""

    analyte_id: str
    precursor_id: str
    source_class_id: str
    mapping_record_id: str
    compatibility_record_id: str
    route_record_id: str = ""


@dataclass(frozen=True)
class AnalyteAssignment:
    """The assignment chosen for one detected class within a minimal combination."""

    analyte_id: str
    precursor_id: str
    source_class_id: str
    mapping_record_id: str
    compatibility_record_id: str
    route_record_id: str


@dataclass(frozen=True)
class MinimalCombinationAssignment:
    """A minimal combination and one supporting assignment for each detected class."""

    source_class_ids: tuple[str, ...]
    assignments: tuple[AnalyteAssignment, ...]


@dataclass(frozen=True)
class CombinationResult:
    """All compatible, minimal and smallest combinations for one set of detected classes."""

    status: Literal["evaluable", "not_evaluable", "unconstrained"]
    reason_code: str | None
    detected_analyte_ids: tuple[str, ...]
    unassigned_analyte_ids: tuple[str, ...]
    source_categories: tuple[str, ...]
    compatible_combinations: int | None
    combinations: tuple[tuple[str, ...], ...]
    single_category_combinations: tuple[str, ...]
    min_categories: int | None
    minimal_combinations: tuple[tuple[str, ...], ...]
    smallest_combinations: tuple[tuple[str, ...], ...]
    categories_in_every_minimal: tuple[str, ...]
    categories_in_some_minimal: tuple[str, ...]
    categories_in_every_smallest: tuple[str, ...]
    categories_in_some_smallest: tuple[str, ...]
    categories_with_same_support: tuple[tuple[str, ...], ...]
    log2_combinations_status: Literal["exact", "undefined"]
    log2_combinations: float | None
    minimal_combination_assignments: tuple[MinimalCombinationAssignment, ...]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"expected true or false, found {value!r}")
    return normalized == "true"


def load_links(
    path: Path, kind: Literal["mapping", "compatibility"]
) -> tuple[Link, ...]:
    """Read fossil_to_precursor.csv (kind 'mapping') or precursor_to_source.csv (kind 'compatibility')."""
    id_column = "mapping_id" if kind == "mapping" else "compatibility_id"
    return tuple(
        Link(
            record_id=row[id_column].strip(),
            disposition=row["disposition"].strip(),
            criteria=row["assignment_criteria"].strip(),
            precursor_id=row["precursor_id"].strip(),
            included=_bool(row["included_in_analysis"]),
            analyte_id=row.get("analyte_id", "").strip(),
            source_class_id=row.get("source_class_id", "").strip(),
        )
        for row in _read_csv(path)
    )


def load_source_categories(path: Path) -> tuple[str, ...]:
    """Read the identifiers of the five source categories, sorted."""
    value = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return tuple(
        sorted(str(row["source_class_id"]) for row in value["source_categories"])
    )


def assignment_options(
    mappings: Sequence[Link], compatibilities: Sequence[Link], *, criteria: str
) -> tuple[AssignmentOption, ...]:
    """Join the links used under strict or broader criteria into assignment options, sorted."""
    if criteria not in {STRICT, BROADER}:
        raise ValueError("criteria must be strict or broader")

    def included(row_criteria: str) -> bool:
        return row_criteria == STRICT or (
            criteria == BROADER and row_criteria == BROADER
        )

    used_mappings = [
        row
        for row in mappings
        if row.disposition == USED and row.included and included(row.criteria)
    ]
    compatibilities_by_precursor: dict[str, list[Link]] = defaultdict(list)
    for row in compatibilities:
        if row.included and included(row.criteria):
            compatibilities_by_precursor[row.precursor_id].append(row)
    options = {
        AssignmentOption(
            analyte_id=mapping.analyte_id,
            precursor_id=mapping.precursor_id,
            source_class_id=compatibility.source_class_id,
            mapping_record_id=mapping.record_id,
            compatibility_record_id=compatibility.record_id,
        )
        for mapping in used_mappings
        for compatibility in compatibilities_by_precursor[mapping.precursor_id]
    }
    return tuple(sorted(options))


def compatible_combinations(
    *,
    detected_analyte_ids: Sequence[str],
    assignment_options: Sequence[AssignmentOption],
    source_categories: Sequence[str],
) -> CombinationResult:
    """Enumerate the compatible, minimal and smallest combinations of source categories."""
    analytes = _normalize_ids(detected_analyte_ids)
    categories = _normalize_ids(source_categories)
    if not categories:
        raise ValueError("no source categories")
    options = tuple(sorted(set(assignment_options)))
    for option in options:
        if option.source_class_id not in categories:
            raise ValueError(
                f"assignment to an unknown source category: {option.source_class_id}"
            )

    if not analytes:
        return _no_result(
            status="unconstrained",
            reason_code="no_detections",
            analytes=analytes,
            unassigned=(),
            categories=categories,
        )

    options_by_analyte: dict[str, tuple[AssignmentOption, ...]] = {}
    for analyte_id in analytes:
        options_by_analyte[analyte_id] = tuple(
            option for option in options if option.analyte_id == analyte_id
        )
    unassigned = tuple(
        analyte_id
        for analyte_id, analyte_options in options_by_analyte.items()
        if not analyte_options
    )
    if unassigned:
        return _no_result(
            status="not_evaluable",
            reason_code="class_without_assignment",
            analytes=analytes,
            unassigned=unassigned,
            categories=categories,
        )

    feasible = _compatible(categories, options_by_analyte)
    feasible_sets = tuple(frozenset(cover) for cover in feasible)
    minimal = tuple(
        cover
        for index, cover in enumerate(feasible)
        if not any(
            other < feasible_sets[index]
            for other_index, other in enumerate(feasible_sets)
            if other_index != index
        )
    )
    minimum_size = min(len(cover) for cover in feasible)
    smallest = tuple(cover for cover in feasible if len(cover) == minimum_size)
    every_minimal, some_minimal = _contributors(minimal)
    every_smallest, some_smallest = _contributors(smallest)
    supporting = tuple(
        _supporting_assignments(cover, options_by_analyte) for cover in minimal
    )
    count = len(feasible)

    return CombinationResult(
        status="evaluable",
        reason_code=None,
        detected_analyte_ids=analytes,
        unassigned_analyte_ids=(),
        source_categories=categories,
        compatible_combinations=count,
        combinations=feasible,
        single_category_combinations=tuple(
            cover[0] for cover in feasible if len(cover) == 1
        ),
        min_categories=minimum_size,
        minimal_combinations=minimal,
        smallest_combinations=smallest,
        categories_in_every_minimal=every_minimal,
        categories_in_some_minimal=some_minimal,
        categories_in_every_smallest=every_smallest,
        categories_in_some_smallest=some_smallest,
        categories_with_same_support=_same_support(
            analytes, categories, options_by_analyte
        ),
        log2_combinations_status="exact",
        log2_combinations=math.log2(count),
        minimal_combination_assignments=supporting,
    )


def _normalize_ids(values: Iterable[str]) -> tuple[str, ...]:
    normalized = tuple(sorted({str(value).strip() for value in values}))
    if "" in normalized:
        raise ValueError("empty identifier")
    return normalized


def _compatible(
    categories: tuple[str, ...],
    options_by_analyte: dict[str, tuple[AssignmentOption, ...]],
) -> tuple[tuple[str, ...], ...]:
    feasible: list[tuple[str, ...]] = []
    for size in range(1, len(categories) + 1):
        for cover in combinations(categories, size):
            cover_set = set(cover)
            if all(
                any(option.source_class_id in cover_set for option in options)
                for options in options_by_analyte.values()
            ):
                feasible.append(cover)
    return tuple(feasible)


def _contributors(
    covers: Sequence[tuple[str, ...]],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    cover_sets = [set(cover) for cover in covers]
    in_every = set.intersection(*cover_sets)
    in_any = set.union(*cover_sets)
    return tuple(sorted(in_every)), tuple(sorted(in_any - in_every))


def _same_support(
    analytes: tuple[str, ...],
    categories: tuple[str, ...],
    options_by_analyte: dict[str, tuple[AssignmentOption, ...]],
) -> tuple[tuple[str, ...], ...]:
    signatures: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for category in categories:
        signature = tuple(
            analyte
            for analyte in analytes
            if any(
                option.source_class_id == category
                for option in options_by_analyte[analyte]
            )
        )
        if signature:
            signatures[signature].append(category)
    return tuple(
        tuple(group) for _, group in sorted(signatures.items()) if len(group) > 1
    )


def _supporting_assignments(
    cover: tuple[str, ...],
    options_by_analyte: dict[str, tuple[AssignmentOption, ...]],
) -> MinimalCombinationAssignment:
    # The first assignment in sorted order is reported for each detected class.
    cover_set = set(cover)
    assignments = tuple(
        _as_assignment(
            min(
                option
                for option in options_by_analyte[analyte_id]
                if option.source_class_id in cover_set
            )
        )
        for analyte_id in sorted(options_by_analyte)
    )
    return MinimalCombinationAssignment(source_class_ids=cover, assignments=assignments)


def _as_assignment(option: AssignmentOption) -> AnalyteAssignment:
    return AnalyteAssignment(
        analyte_id=option.analyte_id,
        precursor_id=option.precursor_id,
        source_class_id=option.source_class_id,
        mapping_record_id=option.mapping_record_id,
        compatibility_record_id=option.compatibility_record_id,
        route_record_id=option.route_record_id,
    )


def _no_result(
    *,
    status: Literal["not_evaluable", "unconstrained"],
    reason_code: str,
    analytes: tuple[str, ...],
    unassigned: tuple[str, ...],
    categories: tuple[str, ...],
) -> CombinationResult:
    return CombinationResult(
        status=status,
        reason_code=reason_code,
        detected_analyte_ids=analytes,
        unassigned_analyte_ids=unassigned,
        source_categories=categories,
        compatible_combinations=None,
        combinations=(),
        single_category_combinations=(),
        min_categories=None,
        minimal_combinations=(),
        smallest_combinations=(),
        categories_in_every_minimal=(),
        categories_in_some_minimal=(),
        categories_in_every_smallest=(),
        categories_in_some_smallest=(),
        categories_with_same_support=(),
        log2_combinations_status="undefined",
        log2_combinations=None,
        minimal_combination_assignments=(),
    )
