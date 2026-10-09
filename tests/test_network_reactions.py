"""Flat-file reading, compound identity states and reaction records, on small synthetic files."""

from biosynthesis import metacyc, reactions

COMPOUNDS = b"""UNIQUE-ID - A
INCHI - InChI=1S/A
INCHI-KEY - InChIKey=A
CHEMICAL-FORMULA - (C 2)
CHEMICAL-FORMULA - (H 6)
//
UNIQUE-ID - B
CHEMICAL-FORMULA - (C 2)
//
UNIQUE-ID - C
SMILES - CC
//
UNIQUE-ID - D
//
"""
CLASSES = b"""UNIQUE-ID - Small-Molecule-Reactions
//
UNIQUE-ID - Chemical-Reactions
TYPES - Small-Molecule-Reactions
//
"""
PATHWAYS = b"""UNIQUE-ID - PWY-1
REACTION-LAYOUT - (RXN-1 (:LEFT-PRIMARIES A) (:DIRECTION :L2R) (:RIGHT-PRIMARIES B))
//
"""
ENZRXNS = b"""UNIQUE-ID - ENZRXN-1
REACTION - RXN-1
REACTION-DIRECTION - LEFT-TO-RIGHT
//
"""
REACTIONS = b"""UNIQUE-ID - RXN-1
TYPES - Chemical-Reactions
LEFT - A
^COEFFICIENT - 2
RIGHT - B
^COEFFICIENT - 2
REACTION-DIRECTION - LEFT-TO-RIGHT
IN-PATHWAY - PWY-1
ENZYMATIC-REACTION - ENZRXN-1
EC-NUMBER - EC-1.1.1.1
//
UNIQUE-ID - RXN-2
TYPES - Chemical-Reactions
LEFT - A
RIGHT - C
REACTION-DIRECTION - REVERSIBLE
//
UNIQUE-ID - RXN-3
TYPES - Chemical-Reactions
LEFT - A
^COEFFICIENT - n
RIGHT - D
//
"""


def index():
    return reactions.load_source_index(
        {
            "reactions.dat": REACTIONS,
            "compounds.dat": COMPOUNDS,
            "classes.dat": CLASSES,
            "pathways.dat": PATHWAYS,
            "enzrxns.dat": ENZRXNS,
        }
    )


def records():
    return list(metacyc.iter_file_records(REACTIONS, "reactions.dat"))


def test_reader_keeps_fields_annotations_and_continuations():
    source = (
        b"UNIQUE-ID - X\nCOMMENT - first\n/second\nLEFT - A\n^COEFFICIENT - 2\n//\n"
    )
    (record,) = metacyc.iter_file_records(source, "x.dat")
    names = [field.field_name_raw for field in record.field_occurrences]
    assert names == [b"UNIQUE-ID", b"COMMENT", b"LEFT"]
    assert record.field_occurrences[1].value_raw == b"first\n/second"
    (annotation,) = record.annotation_occurrences
    assert (
        annotation.attached_to_field_occurrence_id
        == record.field_occurrences[2].field_occurrence_id
    )
    assert record.record_locator == f"x.dat#bytes=0:{len(source)}"


def test_compound_tiers_follow_structure_fields():
    tiers = {
        row["compound_id"]: row["chemical_identity_completeness_tier"]
        for row in metacyc.iter_compound_records(COMPOUNDS)
    }
    assert tiers == {
        "A": metacyc.TIER_I1,
        "B": metacyc.TIER_I3,
        "C": metacyc.TIER_I2,
        "D": metacyc.TIER_I4,
    }
    states = reactions.compound_identity_states(COMPOUNDS)
    assert states["C"].blocking_code == "compound_identity_incomplete_structure"
    assert states["A"].blocking_code is None


def test_strict_reaction_has_one_oriented_reaction_and_a_key():
    record = reactions.reaction_record(records()[0], index())
    assert record["strict_status"] == reactions.STRICT_STATUS
    assert record["in_pathway"] is True
    (oriented,) = record["directed_orientation_records"]
    assert oriented["identity"]["orientation"] == "left_to_right"
    assert oriented["exact_transformation_key"].startswith("sha256:")
    assert {
        o["provenance_kind"]
        for o in record["source_evidence"]["direction_observations"]
    } == {
        "direct_reaction_direction",
        "pathway_reaction_layout_direction",
        "enzymatic_reaction_direction",
    }


def test_incomplete_structure_and_symbolic_coefficient_are_not_strict():
    second = reactions.reaction_record(records()[1], index())
    assert "compound_identity_incomplete_structure" in second["strict_reason_codes"]
    assert second["source_reaction_v2"] is None
    third = reactions.reaction_record(records()[2], index())
    assert third["symbolic_coefficient"] is True
    assert "symbolic_coefficient" in third["strict_reason_codes"]
    assert "missing_direction_observation" in third["outer_blocker_codes"]


def test_transformation_key_ignores_common_factors():
    def directed(coefficient):
        occurrence = lambda local, locator: {  # noqa: E731
            "occurrence_locator": locator,
            "compound": {
                "namespace": "metacyc",
                "local_id": local,
                "compartment": None,
                "chemical_state": "s",
            },
            "coefficient": {"status": "exact", "value": coefficient},
            "source_annotations": {"status": "empty", "sha256": None},
        }
        source = {
            "left_occurrences": [occurrence("A", "l")],
            "right_occurrences": [occurrence("B", "r")],
        }
        return {"source_reaction": source, "orientation": "left_to_right"}

    assert reactions.exact_transformation_key(
        directed("2"), "x"
    ) == reactions.exact_transformation_key(directed("1"), "x")
    assert reactions.exact_transformation_key(
        directed("1"), "x"
    ) != reactions.exact_transformation_key(directed("1"), "y")


def test_coefficients():
    annotation = lambda value: metacyc.AnnotationOccurrence(  # noqa: E731
        "a", 1, b"COEFFICIENT", value, metacyc.SourceSpan("x", 0, 1, 1, 1), value, "f"
    )
    assert reactions._coefficient([]) == ({"status": "exact", "value": "1"}, None)
    assert reactions._coefficient([annotation(b"4/2")]) == (
        {"status": "exact", "value": "2"},
        None,
    )
    assert (
        reactions._coefficient([annotation(b"0")])[1]
        == "malformed_or_nonpositive_coefficient"
    )
    assert reactions._coefficient([annotation(b"n+1")])[1] == "symbolic_coefficient"
