"""Set comparison and precursor relation of biomarker pairs."""

import math

from biosynthesis import enzymes


def test_overlap_gives_fractions_bits_and_jaccard():
    result = enzymes.compare({"a", "b", "c", "d"}, {"c", "d", "e"})
    assert result["relation"] == "overlap"
    assert result["shared"] == ["c", "d"]
    assert result["left_fraction"] == 0.5
    assert result["left_bits"] == 1.0
    assert math.isclose(result["right_bits"], math.log2(3 / 2))
    assert result["jaccard"] == 2 / 5


def test_identical_sets_give_no_reduction():
    result = enzymes.compare({"a", "b"}, {"a", "b"})
    assert result["relation"] == "identical"
    assert (result["left_fraction"], result["left_bits"]) == (1.0, 0.0)
    assert (result["right_fraction"], result["right_bits"]) == (1.0, 0.0)


def test_disjoint_sets_have_zero_fraction_and_no_bits():
    result = enzymes.compare({"a"}, {"b", "c"})
    assert result["relation"] == "disjoint"
    assert result["left_fraction"] == 0.0
    assert result["left_bits"] is None and result["right_bits"] is None
    assert result["jaccard"] == 0.0


def test_empty_sides_are_reported_separately():
    assert enzymes.compare(set(), set())["relation"] == "both_empty"
    assert enzymes.compare(set(), {"a"})["relation"] == "left_empty"
    assert enzymes.compare({"a"}, set())["relation"] == "right_empty"
    both = enzymes.compare(set(), set())
    assert both["jaccard"] is None and both["left_fraction"] is None
    one = enzymes.compare(set(), {"a"})
    assert one["jaccard"] is None
    assert one["left_fraction"] is None and one["right_fraction"] is None


def test_precursor_relation():
    reactants = {1: {"B"}, 2: {"A"}, 3: set()}
    assert (
        enzymes.precursor_relation((1, "A"), (2, "B"), reactants)
        == "each_precursor_of_other"
    )
    assert (
        enzymes.precursor_relation((1, "A"), (3, "C"), {3: {"A"}})
        == "left_precursor_of_right"
    )
    assert (
        enzymes.precursor_relation((1, "A"), (3, "C"), {1: {"C"}})
        == "right_precursor_of_left"
    )
    assert enzymes.precursor_relation((1, "A"), (3, "C"), reactants) == "no"
    assert (
        enzymes.precursor_relation((1, "A"), (4, "SUPP-CPD-X"), reactants)
        == "not_tested"
    )
