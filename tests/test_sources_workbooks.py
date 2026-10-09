"""Reading workbook cells: ages, value types and sample labels."""

from biosynthesis.workbooks import (
    FossilClassColumn,
    Layout,
    WorkbookLayout,
    assign_samples,
    normalize_sample_label,
    parse_age,
    value_type,
)


def test_ages_keep_their_kind_and_bounds():
    assert parse_age(541).kind == "exact"
    closed = parse_age("635-720 Ma")
    assert (closed.kind, closed.lower_ma, closed.upper_ma, closed.central_ma) == (
        "closed_interval",
        635.0,
        720.0,
        None,
    )
    uncertain = parse_age("1100 ± 50")
    assert (uncertain.lower_ma, uncertain.upper_ma, uncertain.central_ma) == (
        1050.0,
        1150.0,
        1100.0,
    )
    assert parse_age("<750").kind == "upper_bound"
    approximate = parse_age("~750")
    assert (approximate.kind, approximate.central_ma, approximate.status) == (
        "approximate",
        750.0,
        "approximate_no_interval",
    )
    assert parse_age(None).kind == "missing"


def test_value_types_separate_zero_blank_text_and_low_signal():
    assert value_type(None, italic=False) == "blank"
    assert value_type("  ", italic=False) == "blank"
    assert value_type(0, italic=False) == "reported_zero"
    assert value_type(0.2, italic=False) == "reported_positive"
    assert value_type(0.2, italic=True) == "reported_positive_low_signal"
    assert value_type("n.d.", italic=False) == "qualified_text"


def test_sample_labels_are_compared_without_case_or_punctuation():
    assert (
        normalize_sample_label("B04026 - i")
        == normalize_sample_label("B04026-I")
        == "b04026i"
    )


def row(workbook, number, label, analyte, used=True):
    return {
        "workbook_id": workbook,
        "sample_key": f"{workbook}:S:{number}",
        "source_row": number,
        "source_coordinate": f"C{number}",
        "source_sample_id_text": label,
        "analyte_id": analyte,
        "class_used": used,
    }


def test_joined_samples_standards_and_unused_classes():
    layout = Layout(
        workbooks=(
            WorkbookLayout(
                "a",
                "a.xlsx",
                "S",
                1,
                9,
                "A",
                "B",
                frozenset(),
                (FossilClassColumn("C", "x", "g", "d", True),),
            ),
        ),
        standard_labels=frozenset({"Std"}),
        compilation_labels=frozenset({"From elsewhere"}),
    )
    rows = [
        row("a", 2, "11J024 i", "x"),
        row("b", 3, "11J024 I", "y"),
        row("a", 4, "Std", "x"),
        row("a", 5, "B1", "bin", used=False),
    ]
    crosswalk = [
        {
            "normalized_sample_label": "11j024i",
            "si1_sample_key": "a:S:2",
            "si2_sample_key": "b:S:3",
            "joined": "true",
        }
    ]
    out = {r["sample_key"]: r for r in assign_samples(rows, layout, crosswalk)}
    assert (
        out["a:S:2"]["canonical_sample_key"]
        == out["b:S:3"]["canonical_sample_key"]
        == "brocks_physical:11j024i"
    )
    assert out["a:S:2"]["sample_record"] == "both_workbooks"
    assert out["a:S:2"]["source_rows"] == "a:S:2;b:S:3"
    assert out["a:S:4"]["included_sample"] is False
    assert out["a:S:4"]["exclusion_reason"] == "standard_not_geological_sample"
    assert out["a:S:5"]["value_used"] is False
    assert out["a:S:5"]["exclusion_reason"] == "source_interpretation_category"
    assert out["a:S:5"]["sample_record"] == "one_workbook"
