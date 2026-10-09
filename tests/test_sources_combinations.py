"""Compatible, minimal and smallest combinations of source categories on a small example."""

import math

from biosynthesis.sources import (
    AssignmentOption,
    Link,
    assignment_options,
    compatible_combinations,
)

CATEGORIES = ("bacteria", "eukaryotes", "hypothesis")


def option(analyte, category, mapping="m", compatibility=None):
    return AssignmentOption(
        analyte,
        "p_" + analyte,
        category,
        f"{mapping}_{analyte}",
        compatibility or f"c_{category}",
    )


OPTIONS = (
    option("hopane", "bacteria"),
    option("sterane", "bacteria"),
    option("sterane", "eukaryotes"),
    option("cyclosterane", "eukaryotes"),
    option("cyclosterane", "hypothesis"),
)


def test_combinations_minimal_and_smallest():
    result = compatible_combinations(
        detected_analyte_ids=["sterane", "hopane", "cyclosterane"],
        assignment_options=OPTIONS,
        source_categories=CATEGORIES,
    )
    assert result.status == "evaluable"
    assert result.combinations == (
        ("bacteria", "eukaryotes"),
        ("bacteria", "hypothesis"),
        ("bacteria", "eukaryotes", "hypothesis"),
    )
    assert result.compatible_combinations == 3
    assert result.minimal_combinations == (
        ("bacteria", "eukaryotes"),
        ("bacteria", "hypothesis"),
    )
    assert result.min_categories == 2
    assert result.categories_in_every_minimal == ("bacteria",)
    assert result.log2_combinations == math.log2(3)
    first = result.minimal_combination_assignments[0]
    assert [a.analyte_id for a in first.assignments] == [
        "cyclosterane",
        "hopane",
        "sterane",
    ]
    # The sterane can come from either category of the combination; the first in sorted order is reported.
    assert first.assignments[2].source_class_id == "bacteria"


def test_unassigned_class_is_not_evaluable_and_no_detection_is_unconstrained():
    missing = compatible_combinations(
        detected_analyte_ids=["hopane", "unknown"],
        assignment_options=OPTIONS,
        source_categories=CATEGORIES,
    )
    assert (missing.status, missing.reason_code, missing.unassigned_analyte_ids) == (
        "not_evaluable",
        "class_without_assignment",
        ("unknown",),
    )
    empty = compatible_combinations(
        detected_analyte_ids=[],
        assignment_options=OPTIONS,
        source_categories=CATEGORIES,
    )
    assert (empty.status, empty.compatible_combinations) == ("unconstrained", None)


def test_broader_criteria_add_links_to_the_strict_ones():
    mappings = (
        Link("m1", "used", "strict", "p1", True, analyte_id="a"),
        Link("m2", "used", "broader", "p2", True, analyte_id="a"),
        Link(
            "m3",
            "excluded_source_interpretation",
            "excluded",
            "",
            False,
            analyte_id="b",
        ),
    )
    compatibilities = (
        Link("c1", "used", "strict", "p1", True, source_class_id="s1"),
        Link("c2", "used", "broader", "p2", True, source_class_id="s2"),
    )
    strict = assignment_options(mappings, compatibilities, criteria="strict")
    broader = assignment_options(mappings, compatibilities, criteria="broader")
    assert [o.source_class_id for o in strict] == ["s1"]
    assert [o.source_class_id for o in broader] == ["s1", "s2"]
