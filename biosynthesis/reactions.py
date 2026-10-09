"""Reaction records built from the MetaCyc flat files.

Each reaction of ``reactions.dat`` becomes one record that keeps the source
fields with their byte ranges, the participants with their coefficients and
compound identity states, the direction observations from the reaction,
enzymatic-reaction and pathway files, and the oriented reactions that follow
from them. A reaction meets the strict criteria when it is a small-molecule
reaction whose every participant is a compound with an exact identity state
and an exact numeric coefficient, and whose source directions do not conflict.
"""

from __future__ import annotations

import base64
import hashlib
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from biosynthesis.metacyc import (
    TIER_I1,
    TIER_I2,
    TIER_I3,
    TIER_I4,
    TIER_IX,
    FlatFileRecord,
    AnnotationOccurrence,
    FieldOccurrence,
    base64_text,
    iter_compound_records,
    iter_file_records,
    read_file,
)
from biosynthesis.records import canonical_line, line_content_id

# Identifier of the MetaCyc 29.5 export: the SHA-256 of its description. Stored
# in the network records; part of the record identifiers.
SOURCE_SNAPSHOT_ID = (
    "sha256:baea792164bb24177240d61148db86b3e4495c53afb1062a2f18f74f5e3589e3"
)

# Stored in the network records; part of the record identifiers.
REACTION_RECORD_SCHEMA_VERSION = "reaction-record/v1"
SOURCE_EVIDENCE_SCHEMA_VERSION = "metacyc-reaction-source-evidence/v1"
IDENTITY_STATE_SCHEMA_VERSION = "metacyc-source-identity-state/v1"
LOCATION_SIGNATURE_SCHEMA_VERSION = "compartment-state-signature/v1"
LOCATION_STATE_SCHEMA_VERSION = "compartment-state/v1"
STRICT_STATUS = "strict"
NOT_STRICT_STATUS = "not_strict"
NOT_SMALL_MOLECULE = "not_small_molecule_reaction"
NO_ORIENTED_REACTION = "missing_exact_directed_v2_projection"
UNPARSED_FIELD_STATUS = "retained_unparsed_not_used_for_filter_or_weight"
DIRECTION_INTERPRETATION = "source_supported_orientation_union_projection"
LOCATION_PARSED = "unresolved_source_location_evidence_not_semantically_parsed"
LOCATION_ABSENT = "unresolved_no_source_location_semantics"

REACTION_FILES = (
    "reactions.dat",
    "compounds.dat",
    "classes.dat",
    "pathways.dat",
    "enzrxns.dat",
)
SMALL_MOLECULE_CLASS = "Small-Molecule-Reactions"

_DIRECT_LTR = frozenset({b"LEFT-TO-RIGHT", b"PHYSIOL-LEFT-TO-RIGHT"})
_DIRECT_RTL = frozenset({b"RIGHT-TO-LEFT", b"PHYSIOL-RIGHT-TO-LEFT"})
_REVERSIBLE = b"REVERSIBLE"
_PATHWAY_REACTION_RE = re.compile(rb"^\(([^\s()]+)")
_PATHWAY_DIRECTION_RE = re.compile(rb"\(:DIRECTION\s+(:L2R|:R2L)\)")
_NUMERIC_LIKE_RE = re.compile(r"^[+-]?(?:[0-9]+(?:\.[0-9]+)?|[0-9]+/[0-9]+)$")
_INTEGER_RE = re.compile(r"^[1-9][0-9]*$")
_FRACTION_RE = re.compile(r"^([1-9][0-9]*)/([1-9][0-9]*)$")
_TRUE_LITERALS = frozenset({b"T", b":TRUE", b"TRUE"})
_PARSED_REACTION_FIELDS = frozenset(
    {
        b"UNIQUE-ID",
        b"TYPES",
        b"LEFT",
        b"RIGHT",
        b"REACTION-DIRECTION",
        b"ENZYMATIC-REACTION",
        b"EC-NUMBER",
        b"SPONTANEOUS?",
        b"IN-PATHWAY",
    }
)
# Problems that keep a reaction out of the strict criteria and out of the
# oriented-reaction projection. Stored in the network records.
_BLOCKING_CODES = frozenset(
    {
        "ambiguous_compound_identity",
        "class_participant_not_v2_projectable",
        "conflicting_direction_observations",
        "compound_identity_conflicting_or_invalid_source_fields",
        "compound_identity_formula_only_not_exact",
        "compound_identity_incomplete_structure",
        "compound_identity_state_missing",
        "dangling_capability_reference",
        "malformed_or_nonpositive_coefficient",
        "missing_direction_observation",
        "missing_literal_side",
        "multiple_coefficient_annotations",
        "source_reference_missing",
        "unqualified_participant_not_v2_projectable",
        "unknown_direction_literal",
    }
)


@dataclass(frozen=True, slots=True)
class CompoundIdentityState:
    """The identity state of one compound, derived from its structure fields."""

    compound_id: str
    identity_state_id: str
    identity_status: str
    completeness_tier: str
    blocking_code: str | None
    compound_record_id: str
    compound_record_bytes: bytes

    @property
    def interpretation(self) -> str:
        return (
            f"{self.identity_status}|tier={self.completeness_tier}|"
            f"state={self.identity_state_id}|control={self.compound_record_id}"
        )


@dataclass(frozen=True, slots=True)
class SourceIndex:
    """What the reaction records need from the compound, class, pathway and enzymatic-reaction files."""

    reactions_sha256: str
    compound_states: Mapping[str, CompoundIdentityState]
    class_ids: frozenset[str]
    small_molecule_class_ids: frozenset[str]
    pathway_ids: frozenset[str]
    enzymatic_reaction_ids: frozenset[str]
    linked_direction_observations: Mapping[str, tuple[Mapping[str, object], ...]]


def compound_identity_states(source: bytes) -> dict[str, CompoundIdentityState]:
    """Assign each compound of ``compounds.dat`` an identity state.

    Compounds with a complete InChI, InChIKey and formula (I1), a formula only
    (I3) or no structure field (I4) have an exact identity; compounds with
    incomplete (I2) or conflicting (IX) structure fields do not.
    """
    states: dict[str, CompoundIdentityState] = {}
    for record in iter_compound_records(source):
        compound_id = str(record["compound_id"])
        tier = str(record["chemical_identity_completeness_tier"])
        if tier == TIER_I1:
            status = "identity_exact"
            basis = "source_declared_complete_inchi_inchikey_valid_formula"
            blocking = None
        elif tier in {TIER_I3, TIER_I4}:
            status = "identity_exact"
            basis = (
                "exact_source_frame_no_structure_formula_retained"
                if tier == TIER_I3
                else "exact_source_frame_no_structure"
            )
            blocking = None
        elif tier == TIER_I2:
            status = "identity_quarantined"
            basis = "incomplete_source_structure_fields"
            blocking = "compound_identity_incomplete_structure"
        elif tier == TIER_IX:
            status = "identity_quarantined"
            basis = "conflicting_or_invalid_source_fields"
            blocking = "compound_identity_conflicting_or_invalid_source_fields"
        else:
            raise ValueError(f"unknown identity tier {tier!r} for {compound_id}")
        state_payload = {
            "schema_version": IDENTITY_STATE_SCHEMA_VERSION,
            "compound_id": compound_id,
            "identity_status": status,
            "completeness_tier": tier,
            "identity_basis": basis,
            "outer_blocker_code": blocking,
            "control_covariate_record": record,
        }
        states[compound_id] = CompoundIdentityState(
            compound_id=compound_id,
            identity_state_id=line_content_id(state_payload),
            identity_status=status,
            completeness_tier=tier,
            blocking_code=blocking,
            compound_record_id=str(record["control_covariate_record_id"]),
            compound_record_bytes=canonical_line(dict(record)),
        )
    return states


def load_source_index(files: Mapping[str, bytes]) -> SourceIndex:
    """Index the compound, class, pathway and enzymatic-reaction files."""
    states = compound_identity_states(files["compounds.dat"])
    class_parents = _class_parents(files["classes.dat"])
    small_molecule = _class_descendants(class_parents, root=SMALL_MOLECULE_CLASS)
    pathway_ids, pathway_directions = _pathway_index(files["pathways.dat"])
    enzrxn_ids, enzrxn_directions = _enzrxn_index(files["enzrxns.dat"])
    linked: dict[str, list[Mapping[str, object]]] = {}
    for reaction_id, observation in (*pathway_directions, *enzrxn_directions):
        linked.setdefault(reaction_id, []).append(observation)
    return SourceIndex(
        reactions_sha256=hashlib.sha256(files["reactions.dat"]).hexdigest(),
        compound_states=states,
        class_ids=frozenset(class_parents),
        small_molecule_class_ids=frozenset(small_molecule),
        pathway_ids=frozenset(pathway_ids),
        enzymatic_reaction_ids=frozenset(enzrxn_ids),
        linked_direction_observations={
            key: tuple(sorted(value, key=lambda item: str(item["observation_id"])))
            for key, value in linked.items()
        },
    )


def build_reaction_records(metacyc: Path) -> list[bytes]:
    """Return one canonical JSON line per reaction, sorted by record identifier."""
    files = {name: read_file(metacyc, name) for name in REACTION_FILES}
    index = load_source_index(files)
    lines: list[tuple[str, bytes]] = []
    for record in iter_file_records(files["reactions.dat"], "reactions.dat"):
        value = reaction_record(record, index)
        lines.append((str(value["source_occurrence_id"]), canonical_line(value)))
    lines.sort(key=lambda item: item[0])
    return [line for _, line in lines]


def reaction_record(record: FlatFileRecord, index: SourceIndex) -> dict[str, object]:
    """Build the record of one reaction: source evidence, participants, directions and orientations."""
    reaction_id = _unique_id(record)
    fields_by_name = _fields_by_name(record)
    annotation_by_id = {
        item.annotation_occurrence_id: item for item in record.annotation_occurrences
    }
    participants: list[dict[str, object]] = []
    left_ids: list[str] = []
    right_ids: list[str] = []
    problems: set[str] = set()
    symbolic = False
    v2_left: list[dict[str, object]] = []
    v2_right: list[dict[str, object]] = []
    participant_ordinal = 0
    for side_raw, side_name, structural_ids, v2_side in (
        (b"LEFT", "left", left_ids, v2_left),
        (b"RIGHT", "right", right_ids, v2_right),
    ):
        for field in fields_by_name.get(side_raw, ()):
            participant_ordinal += 1
            raw_id = _decode_identifier(field.value_raw)
            structural_ids.append(raw_id)
            identity_state = index.compound_states.get(raw_id)
            if identity_state is not None:
                entity_kind = "compound"
                entity_local_id: str | None = raw_id
                interpretation = identity_state.interpretation
                if identity_state.blocking_code is not None:
                    problems.add(identity_state.blocking_code)
            elif raw_id in index.class_ids:
                entity_kind = "class"
                entity_local_id = raw_id
                interpretation = "class_only_not_v2_projectable"
                problems.add("class_participant_not_v2_projectable")
            else:
                entity_kind = "unqualified"
                entity_local_id = None
                interpretation = "unqualified_not_v2_projectable"
                problems.add("unqualified_participant_not_v2_projectable")
            annotations = [
                annotation_by_id[item] for item in field.annotation_occurrence_ids
            ]
            coefficient_annotations = [
                item
                for item in annotations
                if item.annotation_name_raw == b"COEFFICIENT"
            ]
            participant_id = line_content_id(
                {
                    "reaction_local_id": reaction_id,
                    "field_occurrence_id": field.field_occurrence_id,
                    "literal_side": side_name,
                    "ordinal": participant_ordinal,
                    "identity_state_id": (
                        identity_state.identity_state_id
                        if identity_state is not None
                        else None
                    ),
                }
            )
            participants.append(
                {
                    "participant_occurrence_id": participant_id,
                    "field_occurrence_id": field.field_occurrence_id,
                    "literal_side": side_name,
                    "participant_raw_base64": base64_text(field.value_raw),
                    "entity_local_id": entity_local_id,
                    "entity_kind": entity_kind,
                    "coefficient_annotation_occurrence_ids": [
                        item.annotation_occurrence_id
                        for item in coefficient_annotations
                    ],
                    "interpretation_status": interpretation,
                }
            )
            coefficient, coefficient_problem = _coefficient(coefficient_annotations)
            if coefficient_problem == "symbolic_coefficient":
                symbolic = True
            elif coefficient_problem is not None:
                problems.add(coefficient_problem)
            compartment_annotations = [
                item
                for item in annotations
                if item.annotation_name_raw == b"COMPARTMENT"
            ]
            if len(compartment_annotations) > 1:
                problems.add("ambiguous_compound_identity")
                compartment = None
            elif compartment_annotations:
                compartment = _decode_identifier(compartment_annotations[0].value_raw)
            else:
                compartment = None
            if entity_kind == "compound" and coefficient is not None:
                v2_side.append(
                    {
                        "occurrence_locator": field.field_occurrence_id,
                        "compound": {
                            "namespace": "metacyc",
                            "local_id": raw_id,
                            "compartment": compartment,
                            "chemical_state": identity_state.identity_state_id,
                        },
                        "coefficient": coefficient,
                        "source_annotations": _annotation_identity(
                            annotations, identity_state
                        ),
                    }
                )
    if not left_ids or not right_ids:
        problems.add("missing_literal_side")

    directions = [
        *(
            direction_observation(
                source_family="reactions.dat",
                source_record_local_id=reaction_id,
                source_locator=field.field_occurrence_id,
                raw_literal=field.value_raw,
                provenance_kind="direct_reaction_direction",
            )
            for field in fields_by_name.get(b"REACTION-DIRECTION", ())
        ),
        *index.linked_direction_observations.get(reaction_id, ()),
    ]
    known_orientations = sorted(
        {
            o
            for observation in directions
            for o in observation["interpreted_orientations"]
        }
    )
    has_reversible_observation = any(
        set(observation["interpreted_orientations"])
        == {"left_to_right", "right_to_left"}
        for observation in directions
    )
    if not directions:
        problems.add("missing_direction_observation")
    elif any(not observation["interpreted_orientations"] for observation in directions):
        problems.add("unknown_direction_literal")
    if {"left_to_right", "right_to_left"} <= set(
        known_orientations
    ) and not has_reversible_observation:
        problems.add("conflicting_direction_observations")
    possible_orientations = (
        ("left_to_right", "right_to_left")
        if not known_orientations
        or "missing_direction_observation" in problems
        or "unknown_direction_literal" in problems
        else tuple(known_orientations)
    )

    type_ids = {
        _decode_identifier(field.value_raw)
        for field in fields_by_name.get(b"TYPES", ())
    }
    if type_ids - index.class_ids - {SMALL_MOLECULE_CLASS}:
        problems.add("source_reference_missing")
    small_molecule = bool(type_ids & index.small_molecule_class_ids)
    pathway_refs = {
        _decode_identifier(field.value_raw)
        for field in fields_by_name.get(b"IN-PATHWAY", ())
    }
    in_pathway = small_molecule and bool(pathway_refs & index.pathway_ids)
    enzrxn_refs = {
        _decode_identifier(field.value_raw)
        for field in fields_by_name.get(b"ENZYMATIC-REACTION", ())
    }
    if enzrxn_refs - index.enzymatic_reaction_ids:
        problems.add("dangling_capability_reference")

    evidence = normalize_source_evidence(
        {
            "schema_version": SOURCE_EVIDENCE_SCHEMA_VERSION,
            "source_member_path": "reactions.dat",
            "source_member_sha256": index.reactions_sha256,
            "record_start_byte": record.start_byte,
            "record_end_byte": record.end_byte,
            "record_start_line": record.start_line,
            "record_end_line": record.end_line,
            "record_bytes_sha256": hashlib.sha256(record.raw_bytes).hexdigest(),
            "reaction_local_id": reaction_id,
            "field_occurrences": [
                _field_projection(item) for item in record.field_occurrences
            ],
            "annotation_occurrences": [
                _annotation_projection(item) for item in record.annotation_occurrences
            ],
            "participant_occurrences": participants,
            "direction_observations": directions,
            "unrecognized_field_occurrence_ids": [
                item.field_occurrence_id
                for item in record.field_occurrences
                if item.field_name_raw not in _PARSED_REACTION_FIELDS
            ],
            "parse_issue_codes": [],
        }
    )
    evidence_id = str(evidence["evidence_id"])

    blocking = tuple(sorted(problems & _BLOCKING_CODES))
    source_reaction: dict[str, object] | None = None
    oriented: list[dict[str, object]] = []
    if not blocking:
        capability_status, alternatives = _capabilities(
            reaction_id=reaction_id,
            ec_numbers=[
                _decode_identifier(field.value_raw)
                for field in fields_by_name.get(b"EC-NUMBER", ())
            ],
            spontaneous_values=[
                field.value_raw for field in fields_by_name.get(b"SPONTANEOUS?", ())
            ],
        )
        source_reaction = normalize_source_reaction(
            {
                "source_snapshot_id": SOURCE_SNAPSHOT_ID,
                "reaction_local_id": reaction_id,
                "record_locator": record.record_locator,
                "left_occurrences": v2_left,
                "right_occurrences": v2_right,
                "source_annotations": {
                    "status": "identified",
                    "sha256": evidence_id.removeprefix("sha256:"),
                },
                "direction_evidence": {
                    "source_locator": evidence_id,
                    "source_literal": "|".join(known_orientations),
                    "interpretation": DIRECTION_INTERPRETATION,
                },
                "capability_status": capability_status,
                "capability_alternatives": alternatives,
                "unsupported_features": ["symbolic_coefficient"] if symbolic else [],
            }
        )
        oriented = [
            {"source_reaction": source_reaction, "orientation": orientation}
            for orientation in known_orientations
        ]

    reason_codes = set(blocking)
    if not small_molecule:
        reason_codes.add(NOT_SMALL_MOLECULE)
    if symbolic:
        reason_codes.add("symbolic_coefficient")
    if not reason_codes and (source_reaction is None or not oriented):
        reason_codes.add(NO_ORIENTED_REACTION)
    reason_codes = sorted(reason_codes)
    strict = not reason_codes

    facets = unparsed_facets(evidence)
    orientation_records = []
    for directed in oriented:
        identity = {
            "source_snapshot_id": SOURCE_SNAPSHOT_ID,
            "reaction_local_id": reaction_id,
            "orientation": directed["orientation"],
        }
        orientation_record: dict[str, object] = {
            "directed_occurrence_id": line_content_id(
                {
                    "source_occurrence_id": evidence_id,
                    "directed_reaction_identity": identity,
                    "directed_reaction_v2": directed,
                }
            ),
            "identity": identity,
            "directed_reaction_v2": directed,
            "strict_core_status": STRICT_STATUS if strict else NOT_STRICT_STATUS,
            "exact_transformation_key": None,
        }
        if strict:
            state = location_state(evidence, directed, facets)
            orientation_record["exact_transformation_key"] = exact_transformation_key(
                directed, str(state["state_id"])
            )
        orientation_records.append(orientation_record)
    if strict and (source_reaction is None or not orientation_records):
        raise ValueError(f"strict reaction without an oriented reaction: {reaction_id}")

    # Key names and values are stored in the network records and are part of
    # the record identifiers.
    return {
        "schema_version": REACTION_RECORD_SCHEMA_VERSION,
        "source_occurrence_id": evidence_id,
        "source_snapshot_id": SOURCE_SNAPSHOT_ID,
        "reaction_local_id": reaction_id,
        "source_evidence": evidence,
        "structural_reaction": {
            "reaction_local_id": reaction_id,
            "source_evidence_id": evidence_id,
            "left_entity_ids": list(left_ids),
            "right_entity_ids": list(right_ids),
            "possible_orientations": list(possible_orientations),
            "blocker_codes": sorted(
                {*blocking, *(["symbolic_coefficient"] if symbolic else [])}
            ),
        },
        "source_reaction_v2": source_reaction,
        "unparsed_evidence_facets": facets,
        "directed_orientation_records": sorted(
            orientation_records, key=lambda item: str(item["directed_occurrence_id"])
        ),
        "small_molecule_reaction": small_molecule,
        "in_pathway": in_pathway,
        "symbolic_coefficient": symbolic,
        "outer_blocker_codes": list(blocking),
        "strict_status": STRICT_STATUS if strict else NOT_STRICT_STATUS,
        "strict_reason_codes": reason_codes,
    }


def normalize_source_evidence(value: Mapping[str, object]) -> dict[str, object]:
    """Put the source evidence of one reaction in its stored order and identify it.

    Fields and annotations follow their ordinals, participants follow the
    ordinal of their field, direction observations follow their identifiers,
    and the identifier is the SHA-256 of the result.
    """
    fields = [
        {**item, "annotation_occurrence_ids": sorted(item["annotation_occurrence_ids"])}
        for item in value["field_occurrences"]
    ]
    fields.sort(key=lambda item: item["ordinal"])
    annotations = [dict(item) for item in value["annotation_occurrences"]]
    annotations.sort(key=lambda item: item["ordinal"])
    participants = [
        {
            **item,
            "coefficient_annotation_occurrence_ids": sorted(
                item["coefficient_annotation_occurrence_ids"]
            ),
        }
        for item in value["participant_occurrences"]
    ]
    directions = [
        {**item, "interpreted_orientations": sorted(item["interpreted_orientations"])}
        for item in value["direction_observations"]
    ]
    directions.sort(key=lambda item: item["observation_id"])

    field_by_id = {item["field_occurrence_id"]: item for item in fields}
    annotation_by_id = {item["annotation_occurrence_id"]: item for item in annotations}

    def by_ordinal(annotation_id: str) -> int:
        return annotation_by_id.get(annotation_id, {"ordinal": -1})["ordinal"]

    for field in fields:
        field["annotation_occurrence_ids"] = sorted(
            field["annotation_occurrence_ids"], key=by_ordinal
        )
    for participant in participants:
        participant["coefficient_annotation_occurrence_ids"] = sorted(
            participant["coefficient_annotation_occurrence_ids"], key=by_ordinal
        )
    participants.sort(
        key=lambda item: (
            field_by_id[item["field_occurrence_id"]]["ordinal"],
            item["participant_occurrence_id"],
        )
    )
    normalized: dict[str, object] = {
        "schema_version": value["schema_version"],
        "source_member_path": value["source_member_path"],
        "source_member_sha256": value["source_member_sha256"],
        "record_start_byte": value["record_start_byte"],
        "record_end_byte": value["record_end_byte"],
        "record_start_line": value["record_start_line"],
        "record_end_line": value["record_end_line"],
        "record_bytes_sha256": value["record_bytes_sha256"],
        "reaction_local_id": value["reaction_local_id"],
        "field_occurrences": fields,
        "annotation_occurrences": annotations,
        "participant_occurrences": participants,
        "direction_observations": directions,
        "unrecognized_field_occurrence_ids": sorted(
            value["unrecognized_field_occurrence_ids"]
        ),
        "parse_issue_codes": sorted(value["parse_issue_codes"]),
    }
    normalized["evidence_id"] = line_content_id(normalized)
    return normalized


def normalize_source_reaction(value: Mapping[str, object]) -> dict[str, object]:
    """Put a source reaction in its stored order: participants by source locator, capabilities by identifier."""
    return {
        "source_snapshot_id": value["source_snapshot_id"],
        "reaction_local_id": value["reaction_local_id"],
        "record_locator": value["record_locator"],
        "left_occurrences": _sorted_occurrences(value["left_occurrences"]),
        "right_occurrences": _sorted_occurrences(value["right_occurrences"]),
        "source_annotations": dict(value["source_annotations"]),
        "direction_evidence": dict(value["direction_evidence"]),
        "capability_status": value["capability_status"],
        "capability_alternatives": _sorted_alternatives(
            value["capability_alternatives"]
        ),
        "unsupported_features": sorted(value["unsupported_features"]),
    }


def oriented_sides(directed: Mapping[str, object]) -> tuple[list, list]:
    """Return the reactant and product occurrences of an oriented reaction."""
    source = directed["source_reaction"]
    if directed["orientation"] == "left_to_right":
        return source["left_occurrences"], source["right_occurrences"]
    return source["right_occurrences"], source["left_occurrences"]


def exact_transformation_key(
    directed: Mapping[str, object], location_state_id: str
) -> str:
    """Identify the chemical transformation: participants with summed, reduced integer coefficients and location state."""
    reactants, products = oriented_sides(directed)
    left = _aggregate_side(reactants)
    right = _aggregate_side(products)
    left, right = _integer_coefficients(left, right)
    return line_content_id(
        {"reactants": left, "products": right, "compartment_state": location_state_id}
    )


def location_state(
    evidence: Mapping[str, object],
    directed: Mapping[str, object],
    facets: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Record the compartments of the participants and any raw ``RXN-LOCATIONS`` field, without inferring movement between compartments."""
    reactants, products = oriented_sides(directed)
    location_facets = [
        dict(item) for item in facets if item["source_field_name"] == "RXN-LOCATIONS"
    ]
    semantic_status = LOCATION_PARSED if location_facets else LOCATION_ABSENT
    compartments = [
        {
            "directed_side": side,
            "occurrence_locator": occurrence["occurrence_locator"],
            "compound_local_id": occurrence["compound"]["local_id"],
            "compartment": occurrence["compound"]["compartment"],
        }
        for side, occurrences in (("reactant", reactants), ("product", products))
        for occurrence in occurrences
    ]
    signature = {
        "schema_version": LOCATION_SIGNATURE_SCHEMA_VERSION,
        "semantic_status": semantic_status,
        "source_location_evidence_literals": sorted(
            (
                {
                    "field_name_raw_base64": item["field_name_raw_base64"],
                    "value_raw_base64": item["value_raw_base64"],
                }
                for item in location_facets
            ),
            key=lambda item: (
                str(item["field_name_raw_base64"]),
                str(item["value_raw_base64"]),
            ),
        ),
        "participant_compartment_state": sorted(
            (
                {
                    "directed_side": item["directed_side"],
                    "compound_local_id": item["compound_local_id"],
                    "compartment": item["compartment"],
                }
                for item in compartments
            ),
            key=lambda item: (
                str(item["directed_side"]),
                str(item["compound_local_id"]),
                "" if item["compartment"] is None else str(item["compartment"]),
            ),
        ),
        "compartment_change_inferred": False,
    }
    return {
        "schema_version": LOCATION_STATE_SCHEMA_VERSION,
        "state_id": line_content_id(signature),
        "semantic_status": semantic_status,
        "supporting_source_evidence_id": evidence["evidence_id"],
        "source_location_evidence_facets": location_facets,
        "participant_compartments": compartments,
        "compartment_change_inferred": False,
    }


def unparsed_facets(evidence: Mapping[str, object]) -> list[dict[str, object]]:
    """List the source fields the reaction records do not interpret, with their decoded names."""
    unparsed_ids = set(evidence["unrecognized_field_occurrence_ids"])
    facets: list[dict[str, object]] = []
    for field in evidence["field_occurrences"]:
        if field["field_occurrence_id"] not in unparsed_ids:
            continue
        name_base64 = field["field_name_raw_base64"]
        try:
            name = base64.b64decode(name_base64, validate=True).decode("ascii")
            name_status = "decoded_ascii"
        except UnicodeDecodeError:
            name = None
            name_status = "opaque_non_ascii_retained_base64"
        facets.append(
            {
                "field_occurrence_id": field["field_occurrence_id"],
                "source_field_name": name,
                "source_field_name_status": name_status,
                "field_name_raw_base64": name_base64,
                "value_raw_base64": field["value_raw_base64"],
                "semantic_status": UNPARSED_FIELD_STATUS,
            }
        )
    return sorted(facets, key=lambda item: str(item["field_occurrence_id"]))


def direction_observation(
    *,
    source_family: str,
    source_record_local_id: str,
    source_locator: str,
    raw_literal: bytes,
    provenance_kind: str,
) -> dict[str, object]:
    """One direction statement from a reaction, enzymatic-reaction or pathway record, with the orientations it supports."""
    payload: dict[str, object] = {
        "source_family": source_family,
        "source_record_local_id": source_record_local_id,
        "source_locator": source_locator,
        "raw_literal_base64": base64_text(raw_literal),
        "interpreted_orientations": list(
            _interpret_direction(source_family, raw_literal)
        ),
        "provenance_kind": provenance_kind,
    }
    return {"observation_id": line_content_id(payload), **payload}


def _interpret_direction(source_family: str, raw_literal: bytes) -> tuple[str, ...]:
    if source_family == "pathways.dat":
        match = _PATHWAY_DIRECTION_RE.search(raw_literal)
        if match is None:
            return ()
        return ("left_to_right",) if match.group(1) == b":L2R" else ("right_to_left",)
    if raw_literal in _DIRECT_LTR:
        return ("left_to_right",)
    if raw_literal in _DIRECT_RTL:
        return ("right_to_left",)
    if raw_literal == _REVERSIBLE:
        return ("left_to_right", "right_to_left")
    return ()


def _class_parents(source: bytes) -> dict[str, tuple[str, ...]]:
    parents: dict[str, tuple[str, ...]] = {}
    for record in iter_file_records(source, "classes.dat"):
        parents[_unique_id(record)] = tuple(
            _decode_identifier(item.value_raw)
            for item in _fields_by_name(record).get(b"TYPES", ())
        )
    return parents


def _class_descendants(parents: Mapping[str, Sequence[str]], *, root: str) -> set[str]:
    memo: dict[str, bool] = {root: True}

    def reaches_root(node: str, active: frozenset[str]) -> bool:
        if node in memo:
            return memo[node]
        if node in active:
            return False
        result = any(
            parent == root or reaches_root(parent, active | {node})
            for parent in parents.get(node, ())
        )
        memo[node] = result
        return result

    return {node for node in {*parents, root} if reaches_root(node, frozenset())}


def _pathway_index(source: bytes) -> tuple[set[str], list[tuple[str, dict]]]:
    ids: set[str] = set()
    seeds: list[tuple[str, dict]] = []
    for record in iter_file_records(source, "pathways.dat"):
        pathway_id = _unique_id(record)
        ids.add(pathway_id)
        for field in _fields_by_name(record).get(b"REACTION-LAYOUT", ()):
            match = _PATHWAY_REACTION_RE.search(field.value_raw)
            if match is None:
                continue
            seeds.append(
                (
                    _decode_identifier(match.group(1)),
                    direction_observation(
                        source_family="pathways.dat",
                        source_record_local_id=pathway_id,
                        source_locator=field.field_occurrence_id,
                        raw_literal=field.value_raw,
                        provenance_kind="pathway_reaction_layout_direction",
                    ),
                )
            )
    return ids, seeds


def _enzrxn_index(source: bytes) -> tuple[set[str], list[tuple[str, dict]]]:
    ids: set[str] = set()
    seeds: list[tuple[str, dict]] = []
    for record in iter_file_records(source, "enzrxns.dat"):
        enzrxn_id = _unique_id(record)
        ids.add(enzrxn_id)
        fields = _fields_by_name(record)
        for reaction_id in [
            _decode_identifier(item.value_raw) for item in fields.get(b"REACTION", ())
        ]:
            for field in fields.get(b"REACTION-DIRECTION", ()):
                seeds.append(
                    (
                        reaction_id,
                        direction_observation(
                            source_family="enzrxns.dat",
                            source_record_local_id=enzrxn_id,
                            source_locator=field.field_occurrence_id,
                            raw_literal=field.value_raw,
                            provenance_kind="enzymatic_reaction_direction",
                        ),
                    )
                )
    return ids, seeds


def _coefficient(
    annotations: Sequence[AnnotationOccurrence],
) -> tuple[dict | None, str | None]:
    """Read a ``COEFFICIENT`` annotation: 1 when absent, a reduced rational, or a symbolic value."""
    if not annotations:
        return {"status": "exact", "value": "1"}, None
    if len(annotations) > 1:
        return None, "multiple_coefficient_annotations"
    try:
        literal = annotations[0].value_raw.decode("ascii")
    except UnicodeDecodeError:
        return None, "malformed_or_nonpositive_coefficient"
    if _NUMERIC_LIKE_RE.fullmatch(literal):
        try:
            value = Fraction(literal)
        except (ValueError, ZeroDivisionError):
            return None, "malformed_or_nonpositive_coefficient"
        if value <= 0:
            return None, "malformed_or_nonpositive_coefficient"
        text = (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
        return {"status": "exact", "value": text}, None
    if not literal or literal.strip() != literal:
        return None, "malformed_or_nonpositive_coefficient"
    return {"status": "symbolic", "value": literal}, "symbolic_coefficient"


def _capabilities(
    *, reaction_id: str, ec_numbers: Sequence[str], spontaneous_values: Sequence[bytes]
) -> tuple[str, list[dict[str, object]]]:
    """List the EC numbers that can catalyze a reaction; none for spontaneous reactions."""
    if any(value in _TRUE_LITERALS for value in spontaneous_values):
        return "none_required", []
    capability_id = {
        "namespace": "metacyc-reaction-capability",
        "local_id": reaction_id,
    }
    ecs = sorted(set(ec_numbers))
    if not ecs:
        return (
            "unscored",
            [
                {
                    "alternative_id": {
                        "namespace": "metacyc-reaction-capability-alternative",
                        "local_id": reaction_id,
                    },
                    "required_capabilities": [
                        {"capability_id": capability_id, "model_ids": []}
                    ],
                }
            ],
        )
    return (
        "modeled",
        [
            {
                "alternative_id": {
                    "namespace": "metacyc-ec-alternative",
                    "local_id": f"{reaction_id}|{ec}",
                },
                "required_capabilities": [
                    {
                        "capability_id": capability_id,
                        "model_ids": [{"namespace": "ec-number", "local_id": ec}],
                    }
                ],
            }
            for ec in ecs
        ],
    )


def _annotation_identity(
    annotations: Sequence[AnnotationOccurrence], identity_state: CompoundIdentityState
) -> dict[str, object]:
    """Hash the annotations of one participant together with its compound identity state."""
    value = {
        "participant_annotations": [
            {
                "annotation_occurrence_id": item.annotation_occurrence_id,
                "annotation_name_raw_base64": base64_text(item.annotation_name_raw),
                "value_raw_base64": base64_text(item.value_raw),
                "start_byte": item.start_byte,
                "end_byte": item.end_byte,
                "start_line": item.start_line,
                "attached_to_field_occurrence_id": item.attached_to_field_occurrence_id,
            }
            for item in annotations
        ],
        "compound_identity": {
            "identity_state_id": identity_state.identity_state_id,
            "control_covariate_record_id": identity_state.compound_record_id,
            "control_covariate_record_bytes_sha256": hashlib.sha256(
                identity_state.compound_record_bytes
            ).hexdigest(),
        },
    }
    return {
        "status": "identified",
        "sha256": hashlib.sha256(canonical_line(value)).hexdigest(),
    }


def _sorted_occurrences(
    values: Iterable[Mapping[str, object]],
) -> list[dict[str, object]]:
    occurrences = [
        {
            "occurrence_locator": item["occurrence_locator"],
            "compound": dict(item["compound"]),
            "coefficient": dict(item["coefficient"]),
            "source_annotations": dict(item["source_annotations"]),
        }
        for item in values
    ]
    occurrences.sort(
        key=lambda item: (
            item["occurrence_locator"],
            tuple(
                "" if part is None else str(part)
                for part in (
                    item["compound"]["namespace"],
                    item["compound"]["local_id"],
                    item["compound"]["compartment"],
                    item["compound"]["chemical_state"],
                )
            ),
            canonical_line(item),
        )
    )
    return occurrences


def _entity_key(entity: Mapping[str, str]) -> tuple[str, str]:
    return entity["namespace"], entity["local_id"]


def _sorted_alternatives(
    values: Iterable[Mapping[str, object]],
) -> list[dict[str, object]]:
    alternatives = []
    for item in values:
        capabilities = [
            {
                "capability_id": dict(capability["capability_id"]),
                "model_ids": sorted(
                    (dict(model) for model in capability["model_ids"]), key=_entity_key
                ),
            }
            for capability in item["required_capabilities"]
        ]
        capabilities.sort(
            key=lambda capability: _entity_key(capability["capability_id"])
        )
        alternatives.append(
            {
                "alternative_id": dict(item["alternative_id"]),
                "required_capabilities": capabilities,
            }
        )
    alternatives.sort(
        key=lambda alternative: _entity_key(alternative["alternative_id"])
    )
    return alternatives


def _positive_fraction(value: object) -> Fraction:
    if isinstance(value, int) and not isinstance(value, bool) and value > 0:
        return Fraction(value, 1)
    if isinstance(value, str) and _INTEGER_RE.fullmatch(value):
        return Fraction(int(value), 1)
    match = _FRACTION_RE.fullmatch(str(value))
    if match is None:
        raise ValueError(f"coefficient is not an exact positive rational: {value!r}")
    numerator, denominator = (int(token) for token in match.groups())
    return Fraction(numerator, denominator)


def _aggregate_side(
    occurrences: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    totals: dict[bytes, Fraction] = defaultdict(Fraction)
    compounds: dict[bytes, dict[str, object]] = {}
    for occurrence in occurrences:
        coefficient = occurrence["coefficient"]
        if coefficient.get("status") != "exact":
            raise ValueError("a transformation key needs exact coefficients")
        compound = {
            key: occurrence["compound"][key]
            for key in ("namespace", "local_id", "compartment", "chemical_state")
        }
        key = canonical_line(compound)
        compounds[key] = compound
        totals[key] += _positive_fraction(coefficient.get("value"))
    return [
        {"compound": compounds[key], "coefficient": totals[key]}
        for key in sorted(totals)
    ]


def _integer_coefficients(
    left: list[dict], right: list[dict]
) -> tuple[list[dict], list[dict]]:
    """Scale the coefficients of both sides to the smallest integers with the same ratios."""
    fractions = [item["coefficient"] for item in [*left, *right]]
    common_denominator = reduce(math.lcm, (value.denominator for value in fractions), 1)
    integers = [
        value.numerator * (common_denominator // value.denominator)
        for value in fractions
    ]
    divisor = reduce(gcd, integers)
    normalized = [value // divisor for value in integers]
    left_values = [
        {"compound": item["compound"], "coefficient": normalized[cursor]}
        for cursor, item in enumerate(left)
    ]
    right_values = [
        {"compound": item["compound"], "coefficient": normalized[len(left) + cursor]}
        for cursor, item in enumerate(right)
    ]
    return left_values, right_values


def _field_projection(field: FieldOccurrence) -> dict[str, object]:
    return {
        "field_occurrence_id": field.field_occurrence_id,
        "ordinal": field.ordinal,
        "field_name_raw_base64": base64_text(field.field_name_raw),
        "value_raw_base64": base64_text(field.value_raw),
        "start_byte": field.start_byte,
        "end_byte": field.end_byte,
        "start_line": field.start_line,
        "annotation_occurrence_ids": list(field.annotation_occurrence_ids),
    }


def _annotation_projection(annotation: AnnotationOccurrence) -> dict[str, object]:
    return {
        "annotation_occurrence_id": annotation.annotation_occurrence_id,
        "ordinal": annotation.ordinal,
        "annotation_name_raw_base64": base64_text(annotation.annotation_name_raw),
        "value_raw_base64": base64_text(annotation.value_raw),
        "start_byte": annotation.start_byte,
        "end_byte": annotation.end_byte,
        "start_line": annotation.start_line,
        "attached_to_field_occurrence_id": annotation.attached_to_field_occurrence_id,
    }


def _unique_id(record: FlatFileRecord) -> str:
    values = [
        field.value_raw
        for field in record.field_occurrences
        if field.field_name_raw == b"UNIQUE-ID"
    ]
    if len(values) != 1:
        raise ValueError(
            f"record at {record.record_locator} needs exactly one UNIQUE-ID"
        )
    return _decode_identifier(values[0])


def _fields_by_name(record: FlatFileRecord) -> dict[bytes, tuple[FieldOccurrence, ...]]:
    grouped: dict[bytes, list[FieldOccurrence]] = {}
    for field in record.field_occurrences:
        grouped.setdefault(field.field_name_raw, []).append(field)
    return {key: tuple(value) for key, value in grouped.items()}


def _decode_identifier(raw: bytes) -> str:
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")
