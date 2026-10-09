"""Reaction networks: base networks from the reaction records, and final networks.

The base networks keep every reaction record as one edge per supported
orientation, under two criteria: strict (exact compound identity states and
coefficients) and broader (MetaCyc compound identifiers, any source-supported
direction). The final networks add the literature reactions, apply the
decisions on MetaCyc reactions, and come in two versions: supported literature
connections only (primary) and supported plus conditional connections.
"""

from __future__ import annotations

import base64
import csv
import json
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping

from biosynthesis.metacyc import iter_file_records
from biosynthesis.reactions import STRICT_STATUS
from biosynthesis.records import canonical_bytes, content_id

# Criteria names, stored in the network records; part of the record identifiers.
STRICT = "strict_admitted"
BROADER = "source_supported_direction"
CRITERIA = (STRICT, BROADER)
CRITERIA_FILE_NAMES = {STRICT: "strict", BROADER: "broader"}
ORIENTATIONS = ("left_to_right", "right_to_left")
VERSIONS = ("primary", "with_conditional")

# Stored in the network records; part of the record identifiers.
EXACT_STATE = "exact_state"
SOURCE_IDENTIFIER = "source_identifier_not_exact_state"
BASE_LITERATURE_NAMESPACE = "hinkston-supplementary"
LITERATURE_NAMESPACE = "hinkston-source-resolved-intermediate"
COMPARTMENT_ANNOTATION = "Q09NUEFSVE1FTlQ="  # base64 of COMPARTMENT
IN_PATHWAY_FIELD = "SU4tUEFUSFdBWQ=="  # base64 of IN-PATHWAY
ALIAS_IDENTITY_FIELDS = (b"CHEMICAL-FORMULA", b"COMMON-NAME", b"INCHI", b"INCHI-KEY")
# Reason added to the screened-reaction rows of reactions that a decision removes.
REPLACED = "replaced_by_literature_record"

# Text used only for Table S8, not stored in the corrections.
TABLE_ONLY_KEYS = ("table_reason", "single_reaction_test")


# ---------------------------------------------------------------------------
# Network records


def node(identity: dict, status: str) -> dict:
    """A compound node: its identity, identity status and their SHA-256 identifier."""
    return {
        "node_id": content_id({"identity": identity, "identity_status": status}),
        "identity": identity,
        "identity_status": status,
    }


def orient(left: list, right: list, direction: str) -> tuple[list, list]:
    """Return (inputs, outputs) of a reaction in the given orientation."""
    if direction == "left_to_right":
        return left, right
    if direction == "right_to_left":
        return right, left
    raise ValueError(f"Unknown direction: {direction}")


@dataclass
class Network:
    """Nodes, edges and screened-reaction rows under one criteria, with the producing edges of each node."""

    tier: str
    nodes: dict[str, dict]
    edges: dict[str, dict]
    screened: tuple[dict, ...]
    producers: dict[str, tuple[str, ...]]

    @classmethod
    def from_records(
        cls,
        *,
        tier: str,
        nodes: Iterable[dict],
        edges: Iterable[dict],
        screened: Iterable[dict] = (),
    ) -> Network:
        """Index node and edge records; a repeated identical edge is kept once."""
        by_node, by_edge = {}, {}
        producers = defaultdict(list)
        for row in nodes:
            identifier = row["node_id"]
            if identifier in by_node and by_node[identifier] != row:
                raise ValueError("Conflicting node records")
            by_node[identifier] = deepcopy(row)
        for row in edges:
            identifier = row["edge_id"]
            if any(
                value not in by_node
                for side in ("inputs", "outputs")
                for value in row[side]
            ):
                raise ValueError("Edge refers to a missing node")
            if identifier in by_edge:
                if by_edge[identifier] != row:
                    raise ValueError("Conflicting duplicate edge")
                continue
            by_edge[identifier] = deepcopy(row)
            for product in set(row["outputs"]):
                producers[product].append(identifier)
        return cls(
            tier,
            by_node,
            by_edge,
            tuple(deepcopy(list(screened))),
            {key: tuple(sorted(values)) for key, values in producers.items()},
        )


def read_rows(path: Path) -> Iterator[dict]:
    """Read JSON lines."""
    with path.open("rb") as stream:
        for line in stream:
            yield json.loads(line)


def write_rows(path: Path, rows: Iterable[dict]) -> None:
    """Write canonical JSON lines."""
    with path.open("wb") as stream:
        for row in rows:
            stream.write(canonical_bytes(row) + b"\n")


# ---------------------------------------------------------------------------
# Base networks


def _sides(source: dict) -> tuple[dict, list[str]]:
    """Group the recorded participants of a reaction by literal side."""
    raw = source["source_evidence"]["participant_occurrences"]
    sides = {"left": [], "right": []}
    for item in raw:
        sides[item["literal_side"]].append(item)
    return sides, [] if all(sides.values()) else ["recorded_participant_side_empty"]


def _broader_node(source: dict, participant: dict, resolved: dict | None) -> dict:
    """The node of one participant under the broader criteria: MetaCyc identifier and compartment, no identity state."""
    local = participant.get("entity_local_id")
    if not isinstance(local, str) or not local:
        return node(
            {
                "source_occurrence_id": source["source_occurrence_id"],
                "field_occurrence_id": participant["field_occurrence_id"],
                "raw": participant.get("participant_raw_base64"),
            },
            "unresolved_source_entity",
        )
    compound = (resolved or {}).get("compound", {})
    raw_compartments = [
        item.get("value_raw_base64")
        for item in source["source_evidence"].get("annotation_occurrences", [])
        if item.get("attached_to_field_occurrence_id")
        == participant["field_occurrence_id"]
        and item.get("annotation_name_raw_base64") == COMPARTMENT_ANNOTATION
    ]
    compartment = compound.get("compartment")
    try:
        values = {
            base64.b64decode(value, validate=True).decode("ascii")
            for value in raw_compartments
        }
        if (
            "" in values
            or len(values) > 1
            or (values and compartment is not None and compartment not in values)
        ):
            raise ValueError("Conflicting compartment annotations")
        if values:
            compartment = next(iter(values))
    except (ValueError, TypeError, UnicodeError):
        return node(
            {
                "namespace": "metacyc",
                "local_id": local,
                "source_occurrence_id": source["source_occurrence_id"],
                "field_occurrence_id": participant["field_occurrence_id"],
                "raw_compartment_annotations": raw_compartments,
                "parsed_compartment": compartment,
            },
            "unresolved_compartment",
        )
    # A null compartment marks an unqualified MetaCyc identifier.
    identity = {
        "namespace": "metacyc",
        "local_id": local,
        "compartment": compartment,
        "chemical_state": None,
    }
    return node(identity, SOURCE_IDENTIFIER)


def _edge(
    source, tier, direction, left, right, participants, observations, limitations
) -> dict:
    inputs, outputs = orient(left, right, direction)
    identity = {
        "source_occurrence_id": source["source_occurrence_id"],
        "tier": tier,
        "orientation": direction,
    }
    evidence = source["source_evidence"]
    return {
        **identity,
        "edge_id": content_id(identity),
        "reaction_local_id": source["reaction_local_id"],
        "inputs": [item["node_id"] for item in inputs],
        "outputs": [item["node_id"] for item in outputs],
        "participant_records": deepcopy(participants),
        "direction_observations": deepcopy(observations),
        "source_snapshot_id": source["source_snapshot_id"],
        "limitations": sorted(set(limitations)),
        "source_evidence_ref": {
            key: evidence.get(key)
            for key in (
                "evidence_id",
                "source_member_path",
                "source_member_sha256",
                "record_bytes_sha256",
                "record_start_byte",
                "record_end_byte",
            )
        },
        "pathway_fields": [
            deepcopy(field)
            for field in evidence.get("field_occurrences", [])
            if field.get("field_name_raw_base64") == IN_PATHWAY_FIELD
        ],
        "annotation_occurrences": deepcopy(evidence.get("annotation_occurrences", [])),
    }


def reaction_edges(source: dict) -> dict[str, list[dict]]:
    """Turn one reaction record into nodes, edges and screened-reaction rows under both criteria.

    A reaction becomes one edge per source-supported orientation. Under the
    strict criteria only reactions that meet them are kept, with nodes for
    exact compound identity states; under the broader criteria every reaction
    with a supported direction and two non-empty sides is kept, with nodes for
    MetaCyc identifiers. Each orientation not kept gets a row with the reasons.
    """
    raw, participant_errors = _sides(source)
    reaction = source.get("source_reaction_v2") or {}
    parsed = {side: reaction.get(f"{side}_occurrences", []) for side in raw}
    resolved = {
        item["occurrence_locator"]: item
        for values in parsed.values()
        for item in values
    }
    broad = {
        side: [
            _broader_node(source, item, resolved.get(item["field_occurrence_id"]))
            for item in values
        ]
        for side, values in raw.items()
    }
    exact = {
        side: [node(deepcopy(item["compound"]), EXACT_STATE) for item in values]
        for side, values in parsed.items()
    }
    nodes = {
        item["node_id"]: item
        for sides in (broad, exact)
        for values in sides.values()
        for item in values
    }
    observations = source["source_evidence"].get("direction_observations", [])
    directions = sorted(
        {
            direction
            for observation in observations
            for direction in observation.get("interpreted_orientations", [])
            if direction in ORIENTATIONS
        }
    )
    edges, screened = [], []
    blockers = sorted(
        set(
            source.get("strict_reason_codes", [])
            + source.get("outer_blocker_codes", [])
        )
    )
    for tier in CRITERIA:
        candidate_nodes = exact if tier == STRICT else broad
        all_nodes = sorted(
            {item["node_id"] for values in candidate_nodes.values() for item in values}
        )
        for direction in directions or [None]:
            reasons = list(participant_errors)
            if direction is None:
                reasons.append("direction_unresolved")
            matching = [
                item
                for item in source.get("directed_orientation_records", [])
                if item.get("directed_reaction_v2", {}).get("orientation") == direction
                and item.get("strict_core_status") == STRICT_STATUS
            ]
            if tier == STRICT:
                if not matching:
                    reasons.extend(blockers or ["does_not_meet_strict_criteria"])
                if any(
                    Counter(item.get("occurrence_locator") for item in parsed[side])
                    != Counter(item["field_occurrence_id"] for item in raw[side])
                    for side in raw
                ):
                    reasons.append("incomplete_strict_participant_projection")
            kept = []
            if not reasons:
                support = [
                    item
                    for item in observations
                    if direction in item.get("interpreted_orientations", [])
                ]
                edge = _edge(
                    source,
                    tier,
                    direction,
                    candidate_nodes["left"],
                    candidate_nodes["right"],
                    {"raw": raw, "resolved": parsed},
                    support,
                    [] if tier == STRICT else blockers,
                )
                edges.append(edge)
                kept = [edge["edge_id"]]
            possible_products = (
                orient(candidate_nodes["left"], candidate_nodes["right"], direction)[1]
                if direction
                else []
            )
            screened.append(
                {
                    "source_occurrence_id": source["source_occurrence_id"],
                    "reaction_local_id": source["reaction_local_id"],
                    "tier": tier,
                    "orientation": direction,
                    "kept_edge_ids": kept,
                    "reasons": sorted(set(reasons)),
                    "known_product_nodes": sorted(
                        {item["node_id"] for item in possible_products}
                    ),
                    "direction_unresolved_participant_nodes": all_nodes
                    if direction is None
                    else [],
                }
            )
    return {
        "nodes": sorted(nodes.values(), key=lambda row: row["node_id"]),
        "edges": edges,
        "screened": screened,
    }


def base_literature_edges(
    literature: dict, native_ids: set[str]
) -> dict[str, list[dict]]:
    """Add the literature reactions of the base literature file as broader edges, left to right.

    Participants that are MetaCyc compounds keep the MetaCyc namespace; the
    others get the literature namespace. Entries without an identifier are
    section comments and get a screened-reaction row only.
    """
    nodes, edges, screened = {}, [], []
    for ordinal, record in enumerate(literature["reactions"]):
        if "id" not in record:
            screened.append(
                {
                    "source_record_ordinal": ordinal,
                    "tier": BROADER,
                    "reasons": ["non_reaction_comment_record"],
                    "kept_edge_ids": [],
                    "known_product_nodes": [],
                    "direction_unresolved_participant_nodes": [],
                }
            )
            continue
        identifier = record["id"]
        if record.get("reversible") is True:
            raise ValueError(
                "A reversible literature reaction needs a stated direction rule"
            )
        sides = []
        for name in ("substrates", "products"):
            side = []
            for value in record[name]:
                row = node(
                    {
                        "namespace": "metacyc"
                        if value in native_ids
                        else BASE_LITERATURE_NAMESPACE,
                        "local_id": value,
                        "chemical_state": None,
                        "compartment": None,
                    },
                    SOURCE_IDENTIFIER,
                )
                nodes[row["node_id"]] = row
                side.append(row["node_id"])
            sides.append(side)
        occurrence_id = content_id({"reaction_id": identifier})
        identity = {
            "source_occurrence_id": occurrence_id,
            "orientation": "left_to_right",
            "tier": BROADER,
        }
        edge = {
            **identity,
            "edge_id": content_id(identity),
            "reaction_local_id": identifier,
            "inputs": sides[0],
            "outputs": sides[1],
            "supplementary_record": deepcopy(record),
            "source_record_ordinal": ordinal,
            "coefficient_semantics": "participant_multiplicity_only",
            "limitations": [
                "literature_curated_not_strict",
                "structured_primary_source_locator_absent",
            ],
        }
        cofactors = record.get("cofactors", [])
        edge["unresolved_requirements"] = [
            {
                "role": "cofactor",
                "source_literal": value,
                "status": "untyped_source_requirement",
            }
            for value in cofactors
        ]
        if cofactors:
            edge["limitations"].append("auxiliary_cofactor_requirements_untyped")
        edges.append(edge)
        for tier in CRITERIA:
            screened.append(
                {
                    "source_occurrence_id": occurrence_id,
                    "reaction_local_id": identifier,
                    "tier": tier,
                    "orientation": "left_to_right",
                    "kept_edge_ids": [edge["edge_id"]] if tier == BROADER else [],
                    "reasons": []
                    if tier == BROADER
                    else ["literature_record_not_strict"],
                    "known_product_nodes": sides[1] if tier == BROADER else [],
                    "direction_unresolved_participant_nodes": [],
                }
            )
    return {
        "nodes": sorted(nodes.values(), key=lambda row: row["node_id"]),
        "edges": edges,
        "screened": screened,
    }


def compound_ids(compounds: bytes) -> set[str]:
    """The UNIQUE-ID of every compound in ``compounds.dat``."""
    identifiers = set()
    for line in compounds.splitlines(keepends=True):
        if line.startswith(b"UNIQUE-ID - "):
            identifiers.add(
                line[len(b"UNIQUE-ID - ") :].rstrip(b"\r\n").decode("ascii")
            )
    return identifiers


def build_base_networks(
    reaction_records: Path, compounds: bytes, base_literature: bytes
) -> dict[str, Network]:
    """Build the strict and broader base networks from the reaction records and the base literature file."""
    nodes, edges, screened = {}, [], []
    native_ids = compound_ids(compounds)
    for line_number, source in enumerate(read_rows(reaction_records), 1):
        result = reaction_edges(source)
        for row in result["nodes"]:
            if row["node_id"] in nodes and nodes[row["node_id"]] != row:
                raise ValueError("Conflicting node identity")
            nodes[row["node_id"]] = row
            if row["identity"].get("namespace") == "metacyc":
                native_ids.add(row["identity"]["local_id"])
        for row in result["edges"] + result["screened"]:
            row["source_row"] = {
                "input": "source_occurrences",
                "line": line_number,
                "canonical_record_id": content_id(source),
            }
        edges.extend(result["edges"])
        screened.extend(result["screened"])
    literature = base_literature_edges(json.loads(base_literature), native_ids)
    nodes.update({row["node_id"]: row for row in literature["nodes"]})
    edges.extend(literature["edges"])
    screened.extend(literature["screened"])
    return {
        tier: Network.from_records(
            tier=tier,
            nodes=[
                row
                for row in nodes.values()
                if (row["identity_status"] == EXACT_STATE) == (tier == STRICT)
            ],
            edges=[row for row in edges if row["tier"] == tier],
            screened=[row for row in screened if row["tier"] == tier],
        )
        for tier in CRITERIA
    }


def write_base_networks(networks: Mapping[str, Network], folder: Path) -> None:
    """Write nodes and edges sorted by identifier, and the screened-reaction rows in build order."""
    folder.mkdir(parents=True, exist_ok=True)
    for tier, network in networks.items():
        name = CRITERIA_FILE_NAMES[tier]
        write_rows(
            folder / f"{name}.nodes.jsonl",
            (network.nodes[key] for key in sorted(network.nodes)),
        )
        write_rows(
            folder / f"{name}.edges.jsonl",
            (network.edges[key] for key in sorted(network.edges)),
        )
        write_rows(folder / f"{name}.screened_reactions.jsonl", network.screened)


def load_base_network(folder: Path, tier: str) -> Network:
    """Read one base network back from its files."""
    name = CRITERIA_FILE_NAMES[tier]
    return Network.from_records(
        tier=tier,
        nodes=list(read_rows(folder / f"{name}.nodes.jsonl")),
        edges=list(read_rows(folder / f"{name}.edges.jsonl")),
        screened=list(read_rows(folder / f"{name}.screened_reactions.jsonl")),
    )


def local_id_index(networks: Mapping[str, Network]) -> dict[str, dict[str, list[str]]]:
    """For each criteria, the node identifiers of every compound identifier."""
    index = {}
    for tier, network in networks.items():
        by_local = defaultdict(list)
        for key, record in network.nodes.items():
            local = record["identity"].get("local_id")
            if local is not None:
                by_local[local].append(key)
        index[tier] = by_local
    return index


def identity_row(
    identity: dict, tier: str, network: Network, by_local: Mapping, local_ids: set[str]
) -> dict:
    """The nodes that carry an identity's compound identifiers, and the edges producing them."""
    candidates = sorted({key for local in local_ids for key in by_local.get(local, [])})
    return {
        "identity_variant": identity,
        "tier": tier,
        "candidate_node_ids": candidates,
        "direct_producing_edge_ids": sorted(
            {edge for key in candidates for edge in network.producers.get(key, ())}
        ),
    }


def identity_accounting(
    identities: list[dict], networks: Mapping[str, Network]
) -> list[dict]:
    """One row per biomarker identity and criteria.

    Under the broader criteria the panel compound identifier is looked up
    together with the MetaCyc compounds of the identity; candidate states are
    listed, not merged.
    """
    index = local_id_index(networks)
    rows = []
    for identity in identities:
        for tier, network in networks.items():
            local_ids = set(identity["metacyc_compounds"])
            if tier == BROADER:
                local_ids.add(identity["compound_id"])
            rows.append(identity_row(identity, tier, network, index[tier], local_ids))
    return rows


# ---------------------------------------------------------------------------
# Curated corrections


def _stored(value):
    """The stored form of a data-file value: the Table S8 text dropped, recursively."""
    if isinstance(value, list):
        return [_stored(item) for item in value]
    if not isinstance(value, dict):
        return value
    return {
        key: _stored(item) for key, item in value.items() if key not in TABLE_ONLY_KEYS
    }


def _alias_identity_fields(
    compounds: bytes, compound_ids: set[str]
) -> dict[str, dict[str, list[str]]]:
    """Structure fields of the MetaCyc compounds an alias points to, read from ``compounds.dat``."""
    fields: dict[str, dict[str, list[str]]] = {}
    for record in iter_file_records(compounds, "compounds.dat"):
        identifier = record.field_occurrences[0].value_raw.decode("utf-8")
        if identifier not in compound_ids:
            continue
        fields[identifier] = {
            name.decode("ascii"): [
                item.value_raw.decode("utf-8")
                for item in record.field_occurrences
                if item.field_name_raw == name
            ]
            for name in ALIAS_IDENTITY_FIELDS
        }
    return fields


def load_corrections(data: Path, base_literature: bytes, compounds: bytes) -> dict:
    """Assemble the curated corrections from the data files, in the stored form used by the network build."""
    literature = json.loads(
        (data / "literature_reactions.json").read_text(encoding="utf-8")
    )
    metacyc = json.loads(
        (data / "metacyc_corrections.json").read_text(encoding="utf-8")
    )
    with (data / "compound_aliases.csv").open(newline="", encoding="utf-8") as handle:
        aliases = {
            row["literature_id"]: row["metacyc_compound"]
            for row in csv.DictReader(handle)
        }
    records = [_stored(record) for record in literature["records"]]
    states = [
        {"namespace": LITERATURE_NAMESPACE, **_stored(state)}
        for state in literature["compounds"]
    ]
    native = [_stored(decision) for decision in metacyc["reactions"]]
    native_evidence = [
        _stored(row) for row in metacyc["literature_evidence_on_native_reactions"]
    ]
    alias_fields = _alias_identity_fields(
        compounds, {row["to"] for row in metacyc["compound_aliases"]}
    )
    alias_decisions = [
        {**_stored(row), "native_identity_fields": alias_fields[row["to"]]}
        for row in metacyc["compound_aliases"]
    ]
    base = json.loads(base_literature)
    return {
        "expected_supplementary_ids": sorted(
            row["id"] for row in base["reactions"] if "id" in row
        ),
        "records": records,
        "states": states,
        "native_decisions": native,
        "compound_decisions": [_stored(row) for row in metacyc["compounds"]],
        "participant_aliases": aliases,
        "literature_evidence_on_native_reactions": native_evidence,
        "alias_decisions": alias_decisions,
    }


def save_json(path: Path, value) -> None:
    """Write indented JSON with sorted keys."""
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


# ---------------------------------------------------------------------------
# Final networks


def resolve_participants(nodes: dict, aliases: dict) -> dict[str, str]:
    """Map each compound identifier to its single node without a compartment; aliases point to their target's node."""
    choices = {}
    for value in nodes.values():
        identity = value["identity"]
        literal = identity.get("local_id")
        # Literature participants have no compartment and resolve only to the
        # unlocalized node; compartment-specific MetaCyc nodes stay distinct.
        if literal is None or identity.get("compartment") is not None:
            continue
        choices.setdefault(literal, set()).add(value["node_id"])
    resolved = {
        key: next(iter(values)) for key, values in choices.items() if len(values) == 1
    }
    for old, target in aliases.items():
        if old not in choices or target not in resolved:
            raise ValueError("Alias refers to an absent or ambiguous compound: " + old)
        resolved[old] = resolved[target]
    return resolved


def make_edge(
    conversion: dict, record: dict, tier: str, participants: dict[str, str]
) -> dict:
    """One literature connection as an edge, identified by the SHA-256 of its record and connection identifiers."""
    occurrence = content_id(
        {"record_id": record["record_id"], "conversion_id": conversion["id"]}
    )
    edge = {
        "source_occurrence_id": occurrence,
        "orientation": "left_to_right",
        "tier": tier,
    }
    edge["edge_id"] = content_id(edge)
    edge.update(
        reaction_local_id=conversion["id"],
        supplementary_record={
            "id": record["record_id"],
            "substrates": conversion["inputs"],
            "products": conversion["outputs"],
        },
        connection_type=conversion["connection_type"],
        enzyme=conversion.get("enzyme"),
        sources=deepcopy(record["sources"]),
        evidence=deepcopy(conversion.get("evidence", {})),
        unresolved_requirements=deepcopy(conversion.get("unresolved", [])),
        balance_checked=False,
        coefficient_semantics="documented_participant_multiplicity",
    )
    for side in ("inputs", "outputs"):
        try:
            edge[side] = [participants[identifier] for identifier in conversion[side]]
        except KeyError as error:
            raise ValueError(
                "Undeclared or ambiguous participant " + str(error)
            ) from error
    return edge


def build_final_network(
    base: Network, config: dict, version: str
) -> tuple[Network, dict]:
    """Apply the literature records and MetaCyc decisions to one base network.

    Under the broader criteria the base literature edges are removed and the
    literature connections added (conditional ones only in the version with
    conditional connections); MetaCyc reactions marked replace or exclude are
    removed, and annotated ones carry the decision. Returns the network and a
    record of the changes.
    """
    records = config["records"]
    removed = set()
    if base.tier == BROADER:
        removed.update(
            identifier
            for identifier, edge in base.edges.items()
            if edge.get("supplementary_record")
        )
    native_decisions = config.get("native_decisions", [])
    native_ids = {
        row["native_id"]
        for row in native_decisions
        if row["action"] in {"replace", "exclude"}
    }
    removed.update(
        identifier
        for identifier, edge in base.edges.items()
        if edge["reaction_local_id"] in native_ids
    )
    nodes = deepcopy(base.nodes)
    added_nodes = []
    for state in config.get("states", []):
        # A literature compound is identified by its identifier; its description
        # is stored beside the identity.
        value = node(
            {
                "namespace": state["namespace"],
                "local_id": state["local_id"],
                "compartment": state.get("compartment"),
            },
            state["identity_status"],
        )
        value["chemical_state"] = state["chemical_state"]
        if value["node_id"] in nodes and nodes[value["node_id"]] != value:
            raise ValueError("Conflicting literature compound state")
        if value["node_id"] not in nodes:
            nodes[value["node_id"]] = value
            added_nodes.append(value)
    literals = {value["identity"].get("local_id") for value in nodes.values()}
    active_aliases = {
        old: target
        for old, target in config.get("participant_aliases", {}).items()
        if old in literals
    }
    participants = resolve_participants(nodes, active_aliases)
    retained_by_id = {
        key: deepcopy(value) for key, value in base.edges.items() if key not in removed
    }
    changed = {}

    def annotate(native_id: str, decision: dict, unresolved: list) -> list[str]:
        """Attach a decision to the named MetaCyc reaction without changing its chemistry."""
        matches = []
        for identifier, edge in retained_by_id.items():
            if edge["reaction_local_id"] != native_id:
                continue
            changed.setdefault(identifier, {"before": deepcopy(base.edges[identifier])})
            edge.setdefault("literature_annotations", []).append(deepcopy(decision))
            requirements = edge.setdefault("unresolved_requirements", [])
            for requirement in unresolved:
                if requirement not in requirements:
                    requirements.append(deepcopy(requirement))
            matches.append(identifier)
        return matches

    for decision in native_decisions:
        if decision["action"] == "annotate":
            annotate(decision["native_id"], decision, decision.get("unresolved", []))
    native_evidence = {
        (row["record_id"], row["conversion_id"]): row
        for row in config.get("literature_evidence_on_native_reactions", [])
    }
    added, inventory = [], []
    for record in records:
        selected, retained_evidence_ids = [], []
        if base.tier == BROADER:
            for conversion in record["conversions"]:
                if (
                    version == "primary"
                    and conversion["connection_type"] == "conditional"
                ):
                    continue
                link = native_evidence.get((record["record_id"], conversion["id"]))
                if link:
                    evidence = {
                        "record_id": record["record_id"],
                        "conversion": deepcopy(conversion),
                        "sources": deepcopy(record["sources"]),
                        "finding": link["finding"],
                    }
                    matches = annotate(
                        link["native_id"],
                        evidence,
                        link.get("unresolved", conversion.get("unresolved", [])),
                    )
                    if matches:
                        retained_evidence_ids.extend(matches)
                        continue
                    raise ValueError(
                        "A MetaCyc reaction named by a literature record is absent"
                    )
                selected.append(make_edge(conversion, record, base.tier, participants))
        added.extend(selected)
        inventory.append(
            {
                "record_id": record["record_id"],
                "kept_edge_ids": [edge["edge_id"] for edge in selected],
                "retained_native_evidence_edge_ids": sorted(set(retained_evidence_ids)),
                "known_product_nodes": sorted(
                    set(p for edge in selected for p in edge["outputs"])
                ),
                "finding": record["finding"],
                "remaining": record.get("remaining", []),
            }
        )
    for decision in native_decisions:
        for conversion in decision.get("conversions", []):
            # Replacement connections enter only the broader networks.
            if base.tier != BROADER and not conversion.get(
                "meets_strict_criteria", False
            ):
                continue
            if version == "primary" and conversion["connection_type"] == "conditional":
                continue
            added.append(make_edge(conversion, decision, base.tier, participants))
    retained = list(retained_by_id.values())
    screened = deepcopy(list(base.screened))
    for row in screened:
        old_ids = row.get("kept_edge_ids", [])
        if set(old_ids) & removed:
            row["kept_edge_ids"] = [
                identifier for identifier in old_ids if identifier not in removed
            ]
            row.setdefault("reasons", []).append(REPLACED)
    for edge in added:
        screened.append(
            {
                "tier": base.tier,
                "reaction_local_id": edge["reaction_local_id"],
                "source_occurrence_id": edge["source_occurrence_id"],
                "kept_edge_ids": [edge["edge_id"]],
                "known_product_nodes": edge["outputs"],
                "reasons": [],
                "direction_unresolved_participant_nodes": [],
            }
        )
    network = Network.from_records(
        tier=base.tier, nodes=nodes.values(), edges=retained + added, screened=screened
    )
    for key, value in changed.items():
        value["after"] = network.edges[key]
    return network, {
        "removed_edge_ids": sorted(removed),
        "added_nodes": added_nodes,
        "added_edges": added,
        "qualified_native_edges": changed,
        "inventory": inventory,
        # Constant flags: edges not named by a decision and the base nodes are unchanged.
        "unaffected_edges_unchanged": True,
        "base_nodes_unchanged": True,
        "aliases": active_aliases,
        "aliases_not_in_tier": sorted(
            set(config.get("participant_aliases", {})) - set(active_aliases)
        ),
    }


def write_final_network(network: Network, changes: dict, folder: Path) -> None:
    """Write nodes, edges and screened-reaction rows as JSON lines, and the record of changes."""
    folder.mkdir(parents=True, exist_ok=True)
    for part, records in (
        ("nodes", network.nodes.values()),
        ("edges", network.edges.values()),
        ("screened_reactions", network.screened),
    ):
        with (folder / f"{part}.jsonl").open("w") as handle:
            for record in records:
                handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
    save_json(folder / "changes.json", changes)
