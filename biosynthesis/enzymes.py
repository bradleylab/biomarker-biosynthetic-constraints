"""Producing reactions, pathways and enzymes of the biomarkers, and the comparison of pairs.

For each biomarker identity and each set of reaction criteria, the direct
production records are the reaction records of the network that produce one of
the identity's MetaCyc compounds. A pathway is counted for a producing record when
the reaction lists it under IN-PATHWAY and the pathway's recorded direction for the
reaction agrees with production of the biomarker. The enzymes of a biomarker are
the ENZYME values of every enzrxns.dat record whose REACTION is a producing
reaction kept in the primary network.
"""

from __future__ import annotations

import base64
import binascii
import html
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from biosynthesis.reactions import STRICT_STATUS

CRITERIA = ("strict", "broader")
ORIENTATIONS = ("left_to_right", "right_to_left")
PATHWAY_FIELD = "IN-PATHWAY"
PATHWAY_SOURCE = "pathways.dat"
# Direction of a pathway's reaction layout relative to the orientation of a record.
CONSISTENT, MIXED, OPPOSITE, NOT_RECORDED = (
    "consistent",
    "mixed",
    "opposite",
    "not_recorded",
)
LITERATURE_PREFIX = "SUPP-CPD-"

BOTH_EMPTY, LEFT_EMPTY, RIGHT_EMPTY = "both_empty", "left_empty", "right_empty"
DISJOINT, OVERLAP, IDENTICAL = "disjoint", "overlap", "identical"
LEFT_PRECURSOR, RIGHT_PRECURSOR, EACH_PRECURSOR = (
    "left_precursor_of_right",
    "right_precursor_of_left",
    "each_precursor_of_other",
)
NOT_PRECURSOR, NOT_TESTED = "no", "not_tested"


def read_jsonl(path: Path) -> Iterable[dict]:
    """Yield the JSON objects of a JSON-lines file."""
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


# --- Reaction records that produce the biomarkers -----------------------------------


def records_by_compound(path: Path, compound_ids: set[str]) -> dict[str, list[dict]]:
    """Reaction records that list one of the compounds on either side, by compound."""
    index: dict[str, list[dict]] = defaultdict(list)
    for row in read_jsonl(path):
        structural = row.get("structural_reaction") or {}
        participants = {
            str(value)
            for side in ("left_entity_ids", "right_entity_ids")
            for value in structural.get(side, [])
        }
        for compound in sorted(participants & compound_ids):
            index[compound].append(row)
    return dict(index)


def oriented(orientation: str, left: Sequence, right: Sequence) -> tuple[list, list]:
    """Reactants and products of a reaction in the given orientation."""
    if orientation == "left_to_right":
        return list(left), list(right)
    if orientation == "right_to_left":
        return list(right), list(left)
    raise ValueError(f"unknown orientation: {orientation}")


def recorded_orientations(record: Mapping) -> list[str]:
    """Orientations stated by the direction observations of a reaction record."""
    return sorted(
        {
            str(orientation)
            for observation in record.get("source_evidence", {}).get(
                "direction_observations", []
            )
            for orientation in observation.get("interpreted_orientations", [])
            if orientation in ORIENTATIONS
        }
    )


def _local_id(occurrence: Mapping) -> str:
    return str(occurrence["compound"]["local_id"])


def _production_record(identity: Mapping, record: Mapping, **fields) -> dict:
    return {
        "identity_id": identity["identity_id"],
        "biomarker_number": int(identity["biomarker_number"]),
        "compound_id": identity["compound_id"],
        "identity": identity["identity"],
        "reaction_id": record["reaction_local_id"],
        "source_occurrence_id": record["source_occurrence_id"],
        **fields,
    }


def direct_production(identity: Mapping, records: Sequence[Mapping]) -> list[dict]:
    """Producing records of one identity under the strict and broader criteria.

    Under strict criteria each directed record that meets them is a producing record
    when one of the identity's compounds is among its products; under broader
    criteria each orientation stated by the reaction's direction observations is.
    """
    targets = set(map(str, identity.get("metacyc_compounds", [])))
    if not targets:
        return []
    output = []
    for record in records:
        structural = record.get("structural_reaction") or {}
        left_ids = list(map(str, structural.get("left_entity_ids", [])))
        right_ids = list(map(str, structural.get("right_entity_ids", [])))
        strict_by_orientation: dict[str, list[Mapping]] = defaultdict(list)
        for directed in record.get("directed_orientation_records", []):
            if directed.get("strict_core_status") != STRICT_STATUS:
                continue
            reaction = directed["directed_reaction_v2"]
            orientation = str(reaction["orientation"])
            strict_by_orientation[orientation].append(directed)
            precursors, products = oriented(
                orientation,
                reaction["source_reaction"]["left_occurrences"],
                reaction["source_reaction"]["right_occurrences"],
            )
            product_ids = [_local_id(value) for value in products]
            matched = sorted(targets & set(product_ids))
            if not matched:
                continue
            output.append(
                _production_record(
                    identity,
                    record,
                    reaction_criteria="strict",
                    orientation=orientation,
                    target_compound_ids=matched,
                    precursor_compound_ids=[_local_id(value) for value in precursors],
                    product_compound_ids=product_ids,
                    directed_occurrence_id=directed["directed_occurrence_id"],
                    exact_transformation_id=directed["exact_transformation_key"],
                )
            )
        for orientation in recorded_orientations(record):
            precursor_ids, product_ids = oriented(orientation, left_ids, right_ids)
            matched = sorted(targets & set(product_ids))
            if not matched:
                continue
            strict = strict_by_orientation.get(orientation, [])
            exact_ids = sorted(
                {
                    str(value["exact_transformation_key"])
                    for value in strict
                    if value.get("exact_transformation_key")
                }
            )
            directed_ids = sorted(
                {
                    str(value["directed_occurrence_id"])
                    for value in strict
                    if value.get("directed_occurrence_id")
                }
            )
            output.append(
                _production_record(
                    identity,
                    record,
                    reaction_criteria="broader",
                    orientation=orientation,
                    target_compound_ids=matched,
                    precursor_compound_ids=precursor_ids,
                    product_compound_ids=product_ids,
                    directed_occurrence_id=(
                        directed_ids[0] if len(directed_ids) == 1 else None
                    ),
                    exact_transformation_id=exact_ids[0]
                    if len(exact_ids) == 1
                    else None,
                )
            )
    unique = {_canonical(row): row for row in output}
    return [unique[key] for key in sorted(unique, key=lambda key: _order(unique[key]))]


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _order(row: Mapping) -> tuple:
    return (
        CRITERIA.index(row["reaction_criteria"]),
        str(row["reaction_id"]),
        str(row["orientation"]),
        _canonical(row),
    )


def production_records(
    identities: Sequence[Mapping], index: Mapping[str, Sequence[Mapping]]
) -> list[dict]:
    """Producing records of every identity, in biomarker and identity order."""
    found = []
    for identity in identities:
        adjacent = {}
        for compound in identity.get("metacyc_compounds", []):
            for record in index.get(str(compound), []):
                adjacent[str(record["source_occurrence_id"])] = record
        found += direct_production(
            identity, [adjacent[key] for key in sorted(adjacent)]
        )
    return sorted(
        found,
        key=lambda row: (
            row["biomarker_number"],
            str(row["identity_id"]),
            *_order(row),
        ),
    )


# --- Pathways of the producing records -----------------------------------------------


def pathway_ids(path: Path) -> set[str]:
    """UNIQUE-ID values of pathways.dat."""
    found = set()
    with open(path, "rb") as handle:
        for raw in handle:
            line = raw.rstrip(b"\r\n")
            if line.startswith(b"UNIQUE-ID - "):
                found.add(line[len(b"UNIQUE-ID - ") :].decode("ascii"))
    return found


def decode_ascii(value: object) -> str | None:
    """Text of a base64-encoded ASCII value, or None when it is not one."""
    if not isinstance(value, str):
        return None
    try:
        return base64.b64decode(value, validate=True).decode("ascii")
    except (ValueError, binascii.Error, UnicodeDecodeError):
        return None


def pathway_fields(record: Mapping, pathways: set[str]) -> list[dict]:
    """IN-PATHWAY values of a reaction record that name a pathway of pathways.dat."""
    found = []
    for field in record["source_evidence"]["field_occurrences"]:
        if decode_ascii(field.get("field_name_raw_base64")) != PATHWAY_FIELD:
            continue
        pathway = decode_ascii(field.get("value_raw_base64"))
        if pathway is not None and pathway in pathways:
            found.append(
                {
                    "field_occurrence_id": field["field_occurrence_id"],
                    "pathway_id": pathway,
                }
            )
    return found


def pathway_direction(
    observations: Sequence[Mapping], pathway: str, orientation: str
) -> str:
    """Whether a pathway's recorded layout of the reaction agrees with an orientation."""
    opposite = ORIENTATIONS[1] if orientation == ORIENTATIONS[0] else ORIENTATIONS[0]
    local = [
        observation
        for observation in observations
        if observation.get("source_family") == PATHWAY_SOURCE
        and observation.get("source_record_local_id") == pathway
    ]
    supporting = any(
        orientation in o.get("interpreted_orientations", []) for o in local
    )
    opposing = any(opposite in o.get("interpreted_orientations", []) for o in local)
    if supporting and opposing:
        return MIXED
    if supporting:
        return CONSISTENT
    if opposing:
        return OPPOSITE
    return NOT_RECORDED


def directed_orientation(directed: Mapping) -> str | None:
    """Orientation of a directed reaction record."""
    value = directed.get("orientation")
    if value is None:
        value = directed.get("identity", {}).get("orientation")
    if value is None:
        value = directed.get("directed_reaction_v2", {}).get("orientation")
    return str(value) if value is not None else None


def pathway_placements(record: Mapping, pathways: set[str]) -> dict[tuple, list[dict]]:
    """Pathway placements of a reaction record, keyed by orientation and criteria.

    Each placement gives the pathway, its recorded direction, the products of the
    reaction in that orientation and, under strict criteria, the directed record.
    """
    structural = record["structural_reaction"]
    possible = set(map(str, structural.get("possible_orientations", [])))
    strict = {
        str(directed_orientation(item))
        for item in record.get("directed_orientation_records", [])
        if item.get("strict_core_status") == STRICT_STATUS
        and directed_orientation(item) in ORIENTATIONS
    }
    observations = record["source_evidence"].get("direction_observations", [])
    placements: dict[tuple, list[dict]] = defaultdict(list)
    for field in pathway_fields(record, pathways):
        for orientation in sorted((possible | strict) & set(ORIENTATIONS)):
            common = {
                **field,
                "pathway_direction": pathway_direction(
                    observations, field["pathway_id"], orientation
                ),
            }
            _, products = oriented(
                orientation,
                structural["left_entity_ids"],
                structural["right_entity_ids"],
            )
            placements[(orientation, "broader")].append(
                {
                    **common,
                    "transformation_id": None,
                    "product_ids": list(map(str, products)),
                }
            )
            for directed in record.get("directed_orientation_records", []):
                if directed_orientation(directed) != orientation:
                    continue
                if directed.get("strict_core_status") != STRICT_STATUS:
                    continue
                reaction = directed["directed_reaction_v2"]["source_reaction"]
                _, products = oriented(
                    orientation,
                    reaction["left_occurrences"],
                    reaction["right_occurrences"],
                )
                placements[(orientation, "strict")].append(
                    {
                        **common,
                        "transformation_id": directed["exact_transformation_key"],
                        "product_ids": [_local_id(item) for item in products],
                    }
                )
    return placements


def pathway_anchors(
    production: Sequence[Mapping], records: Mapping[str, Mapping], pathways: set[str]
) -> list[dict]:
    """Pathway of each producing record, with the pathway's recorded direction.

    A producing record is placed in a pathway when its reaction lists the pathway
    under IN-PATHWAY and, under strict criteria, the placement concerns the same
    directed record. Records without a pathway get no row.
    """
    placements_by_record: dict[str, dict] = {}
    anchors = []
    for row in production:
        source = str(row["source_occurrence_id"])
        if source not in placements_by_record:
            placements_by_record[source] = pathway_placements(records[source], pathways)
        candidates = placements_by_record[source].get(
            (row["orientation"], row["reaction_criteria"]), []
        )
        if row["reaction_criteria"] == "strict":
            candidates = [
                item
                for item in candidates
                if item["transformation_id"] == row["exact_transformation_id"]
            ]
        targets = set(map(str, row["target_compound_ids"]))
        for placement in candidates:
            matched = sorted({c for c in placement["product_ids"] if c in targets})
            if not matched:
                continue
            anchors.append(
                {
                    key: row[key]
                    for key in (
                        "identity_id",
                        "biomarker_number",
                        "compound_id",
                        "reaction_criteria",
                        "reaction_id",
                        "source_occurrence_id",
                        "orientation",
                    )
                }
                | {
                    "pathway_id": placement["pathway_id"],
                    "field_occurrence_id": placement["field_occurrence_id"],
                    "pathway_direction": placement["pathway_direction"],
                    "target_compound_ids": matched,
                }
            )
    return sorted(anchors, key=lambda row: (row["biomarker_number"], _canonical(row)))


# --- Sets of each biomarker -----------------------------------------------------------


def kept_occurrences(edges_path: Path) -> set[tuple[str, str]]:
    """Reaction record and orientation of every MetaCyc reaction of a network.

    Literature reactions carry no MetaCyc pathway or enzyme record and are skipped.
    """
    return {
        (edge["source_occurrence_id"], edge["orientation"])
        for edge in read_jsonl(edges_path)
        if not edge.get("supplementary_record")
    }


def is_kept(record: Mapping, kept: set[tuple[str, str]]) -> bool:
    """Whether the network keeps the record's reaction in the record's orientation."""
    return (record["source_occurrence_id"], record["orientation"]) in kept


def biomarker_sets(
    criteria: str,
    production: Sequence[Mapping],
    anchors: Sequence[Mapping],
    kept: set[tuple[str, str]],
    numbers: Iterable[int],
) -> dict[str, dict[int, set[str]]]:
    """Producing reactions and pathways of each biomarker under one set of criteria."""
    sets = {
        name: {number: set() for number in numbers}
        for name in ("reactions", "pathways")
    }
    for anchor in anchors:
        if (
            anchor["reaction_criteria"] == criteria
            and anchor["pathway_direction"] == CONSISTENT
            and is_kept(anchor, kept)
        ):
            sets["pathways"][anchor["biomarker_number"]].add(anchor["pathway_id"])
    for row in production:
        if row["reaction_criteria"] == criteria and is_kept(row, kept):
            sets["reactions"][row["biomarker_number"]].add(row["reaction_id"])
    return sets


def precursors_by_biomarker(
    criteria: str, production: Sequence[Mapping], kept: set[tuple[str, str]]
) -> dict[int, set[str]]:
    """Reactants of the kept producing records of each biomarker."""
    found: dict[int, set[str]] = defaultdict(set)
    for row in production:
        if row["reaction_criteria"] == criteria and is_kept(row, kept):
            found[row["biomarker_number"]].update(row["precursor_compound_ids"])
    return dict(found)


def is_literature(compound_id: str) -> bool:
    """Whether a biomarker is defined from the literature rather than by a MetaCyc compound."""
    return compound_id.startswith(LITERATURE_PREFIX)


def precursor_relation(
    left: tuple[int, str], right: tuple[int, str], reactants: Mapping[int, set[str]]
) -> str:
    """Whether either biomarker of a pair is a reactant of a producing reaction of the other."""
    (left_number, left_id), (right_number, right_id) = left, right
    if is_literature(left_id) or is_literature(right_id):
        return NOT_TESTED
    forward = left_id in reactants.get(right_number, set())
    backward = right_id in reactants.get(left_number, set())
    if forward and backward:
        return EACH_PRECURSOR
    return LEFT_PRECURSOR if forward else RIGHT_PRECURSOR if backward else NOT_PRECURSOR


def compare(left: set[str], right: set[str]) -> dict:
    """Relation, sizes, retained fractions, reductions in bits and Jaccard index of two sets."""
    shared, union = left & right, left | right

    def directional(side: set[str]) -> tuple[float | None, float | None]:
        if not side:
            return None, None
        bits = math.log2(len(side) / len(shared)) if shared else None
        return len(shared) / len(side), bits

    if not left and not right:
        relation = BOTH_EMPTY
    elif not left:
        relation = LEFT_EMPTY
    elif not right:
        relation = RIGHT_EMPTY
    elif not shared:
        relation = DISJOINT
    else:
        relation = IDENTICAL if left == right else OVERLAP
    left_fraction, left_bits = directional(left)
    right_fraction, right_bits = directional(right)
    jaccard = len(shared) / len(union) if union else None
    if relation in (BOTH_EMPTY, LEFT_EMPTY, RIGHT_EMPTY):
        # A pair with an empty set has no fractions, reductions or Jaccard index.
        left_fraction = right_fraction = left_bits = right_bits = jaccard = None
    return {
        "relation": relation,
        "left_count": len(left),
        "right_count": len(right),
        "shared_count": len(shared),
        "shared": sorted(shared),
        "left_fraction": left_fraction,
        "right_fraction": right_fraction,
        "left_bits": left_bits,
        "right_bits": right_bits,
        "jaccard": jaccard,
    }


# --- MetaCyc enzyme records -------------------------------------------------------------


def read_dat(path: Path, attributes: set[str]) -> dict[str, dict[str, list[str]]]:
    """Records of a MetaCyc attribute-value file, keyed by UNIQUE-ID, with the given attributes."""
    records, current = {}, defaultdict(list)
    with open(path, encoding="latin-1") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "//":
                records[current["UNIQUE-ID"][0]] = dict(current)
                current = defaultdict(list)
            elif not line.startswith(("#", "/")):
                key, separator, value = line.partition(" - ")
                if separator and (key == "UNIQUE-ID" or key in attributes):
                    current[key].append(value)
    return records


def display_name(text: str) -> str:
    """A MetaCyc name without HTML tags or entities."""
    return html.unescape(re.sub("<[^>]+>", "", text))


@dataclass(frozen=True)
class EnzymeIndex:
    """Enzymes of each MetaCyc reaction and the organisms they are recorded from."""

    by_reaction: dict[str, set[str]]
    proteins: dict[str, dict[str, list[str]]]
    classes: dict[str, dict[str, list[str]]]

    def enzymes(self, reactions: Iterable[str]) -> set[str]:
        """Enzymes recorded for any of the reactions."""
        return {e for r in reactions for e in self.by_reaction.get(r, ())}

    def components(self, enzyme: str) -> set[str]:
        """Every protein that is a component of a complex, at any depth."""
        found, stack = set(), list(self.proteins[enzyme].get("COMPONENTS", []))
        while stack:
            part = stack.pop()
            if part not in found:
                found.add(part)
                stack.extend(self.proteins.get(part, {}).get("COMPONENTS", []))
        return found

    def taxa(self, enzyme: str) -> tuple[set[str], bool]:
        """Recorded SPECIES of an enzyme, or else of its components (flagged True)."""
        recorded = set(self.proteins[enzyme].get("SPECIES", []))
        if recorded:
            return recorded, False
        parts = self.components(enzyme)
        return {
            s for p in parts for s in self.proteins.get(p, {}).get("SPECIES", [])
        }, True


def enzyme_index(metacyc: Path) -> EnzymeIndex:
    """Enzymes of each reaction from enzrxns.dat, with proteins.dat and classes.dat."""
    by_reaction = defaultdict(set)
    for record in read_dat(metacyc / "enzrxns.dat", {"REACTION", "ENZYME"}).values():
        for reaction in record["REACTION"]:
            by_reaction[reaction].update(record["ENZYME"])
    proteins = read_dat(
        metacyc / "proteins.dat",
        {"SPECIES", "COMPONENTS", "COMMON-NAME", "SYNONYMS", "GENE"},
    )
    classes = read_dat(metacyc / "classes.dat", {"TYPES", "COMMON-NAME", "STRAIN-NAME"})
    return EnzymeIndex(dict(by_reaction), proteins, classes)
