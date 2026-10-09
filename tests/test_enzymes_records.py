"""Producing records, pathway placements and biomarker sets from a small synthetic reaction."""

import base64

from biosynthesis import enzymes


def b64(text: str) -> str:
    return base64.b64encode(text.encode("ascii")).decode("ascii")


def occurrence(compound: str) -> dict:
    return {
        "compound": {
            "local_id": compound,
            "namespace": "metacyc",
            "chemical_state": "s",
        }
    }


def observation(family: str, record: str, orientations: list[str]) -> dict:
    return {
        "source_family": family,
        "source_record_local_id": record,
        "interpreted_orientations": orientations,
    }


# Reaction RXN-1: S = P + W, recorded left to right by the reaction and by PWY-A,
# right to left by PWY-B; one directed record meets the strict criteria.
RECORD = {
    "reaction_local_id": "RXN-1",
    "source_occurrence_id": "occ-1",
    "structural_reaction": {
        "left_entity_ids": ["S"],
        "right_entity_ids": ["P", "W"],
        "possible_orientations": ["left_to_right", "right_to_left"],
    },
    "source_evidence": {
        "direction_observations": [
            observation("reactions.dat", "RXN-1", ["left_to_right"]),
            observation("pathways.dat", "PWY-A", ["left_to_right"]),
            observation("pathways.dat", "PWY-B", ["right_to_left"]),
        ],
        "field_occurrences": [
            {
                "field_occurrence_id": "f1",
                "field_name_raw_base64": b64("IN-PATHWAY"),
                "value_raw_base64": b64("PWY-A"),
            },
            {
                "field_occurrence_id": "f2",
                "field_name_raw_base64": b64("IN-PATHWAY"),
                "value_raw_base64": b64("PWY-B"),
            },
            {
                "field_occurrence_id": "f3",
                "field_name_raw_base64": b64("IN-PATHWAY"),
                "value_raw_base64": b64("NOT-A-PATHWAY"),
            },
            {
                "field_occurrence_id": "f4",
                "field_name_raw_base64": b64("COMMENT"),
                "value_raw_base64": b64("PWY-A"),
            },
        ],
    },
    "directed_orientation_records": [
        {
            "strict_core_status": enzymes.STRICT_STATUS,
            "directed_occurrence_id": "d1",
            "exact_transformation_key": "x1",
            "directed_reaction_v2": {
                "orientation": "left_to_right",
                "source_reaction": {
                    "left_occurrences": [occurrence("S")],
                    "right_occurrences": [occurrence("P"), occurrence("W")],
                },
            },
        }
    ],
}
IDENTITY = {
    "identity_id": "01",
    "biomarker_number": 1,
    "compound_id": "P",
    "identity": "accepted",
    "metacyc_compounds": ["P"],
}


def test_direct_production_under_both_criteria():
    rows = enzymes.direct_production(IDENTITY, [RECORD])
    assert [(r["reaction_criteria"], r["orientation"]) for r in rows] == [
        ("strict", "left_to_right"),
        ("broader", "left_to_right"),
    ]
    assert rows[0]["precursor_compound_ids"] == ["S"]
    assert rows[0]["exact_transformation_id"] == "x1"
    assert rows[1]["directed_occurrence_id"] == "d1"
    # Broader criteria take every orientation stated by any direction observation,
    # here the right-to-left layout of PWY-B, which produces S.
    reverse = enzymes.direct_production(
        {**IDENTITY, "metacyc_compounds": ["S"]}, [RECORD]
    )
    assert [(r["reaction_criteria"], r["orientation"]) for r in reverse] == [
        ("broader", "right_to_left")
    ]
    assert (
        enzymes.direct_production({**IDENTITY, "metacyc_compounds": ["Q"]}, [RECORD])
        == []
    )


def test_pathway_direction_states():
    observations = RECORD["source_evidence"]["direction_observations"]
    assert (
        enzymes.pathway_direction(observations, "PWY-A", "left_to_right")
        == "consistent"
    )
    assert (
        enzymes.pathway_direction(observations, "PWY-B", "left_to_right") == "opposite"
    )
    assert (
        enzymes.pathway_direction(observations, "PWY-C", "left_to_right")
        == "not_recorded"
    )


def test_anchors_and_sets():
    production = enzymes.direct_production(IDENTITY, [RECORD])
    anchors = enzymes.pathway_anchors(production, {"occ-1": RECORD}, {"PWY-A", "PWY-B"})
    found = {
        (a["reaction_criteria"], a["pathway_id"], a["pathway_direction"])
        for a in anchors
    }
    assert found == {
        ("strict", "PWY-A", "consistent"),
        ("strict", "PWY-B", "opposite"),
        ("broader", "PWY-A", "consistent"),
        ("broader", "PWY-B", "opposite"),
    }
    kept = {("occ-1", "left_to_right")}
    sets = enzymes.biomarker_sets("broader", production, anchors, kept, [1, 2])
    assert sets["reactions"] == {1: {"RXN-1"}, 2: set()}
    assert sets["pathways"] == {1: {"PWY-A"}, 2: set()}
    assert enzymes.biomarker_sets("broader", production, anchors, set(), [1])[
        "reactions"
    ] == {1: set()}
    assert enzymes.precursors_by_biomarker("strict", production, kept) == {1: {"S"}}
