"""Matched controls for upstream tracing: stopping at starting compounds, pathway-only reactions,
strict reactions on the broader network, and unused annotations."""

from __future__ import annotations

import base64
import binascii
from collections import defaultdict, deque
from copy import copy
from pathlib import Path
from typing import Iterable

from biosynthesis.graph import BROADER, Network
from biosynthesis.tracing import discover

UNCONTROLLED = "uncontrolled"
STOP_AT_STARTING = "stop_at_starting_compounds"
PATHWAY_ONLY = "pathway_reactions_only"
PATHWAY_AND_STOP = "pathway_reactions_and_stop"
STRICT_REACTIONS = "strict_reactions_broader_compounds"
CONTROL_NAMES = (STOP_AT_STARTING, PATHWAY_ONLY, PATHWAY_AND_STOP, STRICT_REACTIONS)
ORIENTATIONS = ("left_to_right", "right_to_left")


def membership(result: dict) -> dict:
    """The compounds traced, the compounds stored and the reactions of a trace."""
    return {k: result[k] for k in ("required_nodes", "stored_nodes", "edge_ids")}


def controlled_discovery(
    network: Network,
    target: str,
    *,
    stop_nodes: set[str] | None = None,
    allowed_edges: set[str] | None = None,
) -> dict:
    """Trace the target using only the allowed reactions and without expanding stopped compounds.

    A stopped product does not disable a reaction reached through another product.
    None means no restriction; an empty set allows no reaction.
    """
    stops = set() if stop_nodes is None else set(stop_nodes)
    allowed = set(network.edges) if allowed_edges is None else set(allowed_edges)
    if (
        target not in network.nodes
        or not stops <= network.nodes.keys()
        or not allowed <= network.edges.keys()
    ):
        raise ValueError("Unknown target, stopped node, or allowed edge")
    required, stored, used = {target}, {target}, set()
    queue = deque([target])
    while queue:
        compound = queue.popleft()
        if compound in stops:
            continue
        for key in network.producers.get(compound, ()):
            if key not in allowed or key in used:
                continue
            used.add(key)
            edge = network.edges[key]
            stored.update(edge["outputs"])
            for precursor in sorted(set(edge["inputs"])):
                stored.add(precursor)
                if precursor not in required:
                    required.add(precursor)
                    queue.append(precursor)
    return dict(
        required_nodes=sorted(required),
        stored_nodes=sorted(stored),
        edge_ids=sorted(used),
    )


def decode_base64_ascii(value: object) -> tuple[str | None, str]:
    """Decode a base64-encoded source text to ASCII; return None and a reason when decoding fails."""
    if not isinstance(value, str):
        return None, "not_string"
    try:
        raw = base64.b64decode(value, validate=True)
    except (ValueError, binascii.Error):
        return None, "invalid_base64"
    try:
        return raw.decode("ascii"), "decoded_ascii"
    except UnicodeDecodeError:
        return None, "non_ascii"


def pathway_ids(path: Path) -> set[str]:
    """Read the UNIQUE-ID of every record of the MetaCyc pathways file."""
    records: set[str] = set()
    unique = None
    with open(path, "rb") as handle:
        for raw_line in handle:
            line = raw_line.rstrip(b"\r\n")
            if line.startswith(b"UNIQUE-ID - "):
                if unique is not None:
                    raise ValueError("Multiple UNIQUE-ID fields in one pathway record")
                unique = line[len(b"UNIQUE-ID - ") :].decode("ascii")
            if line == b"//":
                if unique is not None:
                    if unique in records:
                        raise ValueError(f"Duplicate pathway UNIQUE-ID: {unique}")
                    records.add(unique)
                unique = None
    if unique is not None:
        raise ValueError("Unterminated final pathway record")
    return records


def pathway_direction(
    observations: list[dict], pathway_id: str, orientation: str
) -> str:
    """Classify the direction a pathway gives a reaction occurrence: consistent, mixed, conflicting or absent."""
    opposite = ORIENTATIONS[1] if orientation == ORIENTATIONS[0] else ORIENTATIONS[0]
    local = [
        observation
        for observation in observations
        if observation.get("source_family") == "pathways.dat"
        and observation.get("source_record_local_id") == pathway_id
    ]
    supporting = [
        o for o in local if orientation in o.get("interpreted_orientations", [])
    ]
    opposing = [o for o in local if opposite in o.get("interpreted_orientations", [])]
    if supporting and opposing:
        return "pathway_local_mixed"
    if supporting:
        return "pathway_local_consistent"
    if opposing:
        return "pathway_local_conflict_only"
    return "pathway_local_absent"


def pathway_masks(network: Network, pathways: set[str]) -> dict[str, set[str]]:
    """For each MetaCyc pathway, the reaction occurrences it contains in a consistent direction."""
    masks = defaultdict(set)
    for key, edge in network.edges.items():
        for field in edge.get("pathway_fields", []):
            label, state = decode_base64_ascii(field["value_raw_base64"])
            if state != "decoded_ascii" or label not in pathways:
                continue
            direction = pathway_direction(
                edge["direction_observations"], label, edge["orientation"]
            )
            if direction == "pathway_local_consistent":
                masks[label].add(key)
    return dict(masks)


def starting_identifiers(
    strict_network: Network, compounds: Iterable[str]
) -> set[tuple]:
    """The starting compounds that have an exact state in the strict network, as (namespace, identifier)."""
    wanted = set(compounds)
    return {
        (node["identity"]["namespace"], node["identity"]["local_id"])
        for node in strict_network.nodes.values()
        if node["identity"].get("namespace") == "metacyc"
        and node["identity"].get("local_id") in wanted
        and node["identity_status"] == "exact_state"
    }


def stopped_states(network: Network, identifiers: set[tuple]) -> set[str]:
    """Every state, in any compartment or chemical state, of the given starting compounds."""
    return {
        key
        for key, node in network.nodes.items()
        if "namespace" in node["identity"]
        and "local_id" in node["identity"]
        and (node["identity"]["namespace"], node["identity"]["local_id"]) in identifiers
    }


def strict_occurrences(strict_network: Network) -> set[tuple[str, str]]:
    """The source reaction occurrences, with orientation, accepted under the strict criteria."""
    return {
        (edge["source_occurrence_id"], edge["orientation"])
        for edge in strict_network.edges.values()
    }


def controls(
    network: Network, target: str, stops: set, masks: dict, strict_accepted: set
) -> dict:
    """Trace the target without restriction and under each control, one rule changed at a time."""
    cache = {}

    def calculate(allowed=None, stopped=None):
        key = (
            None if allowed is None else frozenset(allowed),
            frozenset(stopped or ()),
        )
        if key not in cache:
            cache[key] = controlled_discovery(
                network, target, allowed_edges=allowed, stop_nodes=stopped
            )
        return cache[key]

    consistent = set().union(*masks.values()) if masks else set()
    root_edges = set(network.producers.get(target, ()))
    contexts = {
        label: dict(
            pathway_only=calculate(edges), pathway_and_stop=calculate(edges, stops)
        )
        for label, edges in sorted(masks.items())
        if edges & root_edges
    }
    result = {
        UNCONTROLLED: calculate(),
        STOP_AT_STARTING: calculate(stopped=stops),
        PATHWAY_ONLY: calculate(consistent),
        PATHWAY_AND_STOP: calculate(consistent, stops),
        "pathway_contexts": contexts,
    }
    if network.criteria == BROADER:
        allowed = {
            key
            for key, e in network.edges.items()
            if (e["source_occurrence_id"], e["orientation"]) in strict_accepted
        }
        result[STRICT_REACTIONS] = calculate(allowed)
    return result


def annotations_change_result(
    network: Network, target: str, expected: dict, stops: set
) -> bool:
    """Whether supplying starting-compound or pathway annotations, without using them, changes a trace."""
    unlabeled = copy(network)
    unlabeled.edges = {
        key: {**edge, "pathway_fields": [], "annotation_occurrences": []}
        for key, edge in network.edges.items()
    }
    return (
        membership(discover(network, [target], supplied_annotations=stops)) != expected
        or membership(discover(network, [target], supplied_annotations=network.nodes))
        != expected
        or membership(discover(unlabeled, [target])) != expected
    )


def alias_roots(
    network: Network, declarations: dict, selected: set[str]
) -> dict[str, str]:
    """Map query states whose identifier has a declared alias to the single MetaCyc state of the alias."""
    resolved = {}
    for historical, node in network.nodes.items():
        if historical not in selected:
            continue
        identity = node["identity"]
        local_id = identity.get("local_id")
        if local_id not in declarations:
            continue
        candidates = [
            (key, candidate)
            for key, candidate in network.nodes.items()
            if candidate["identity"].get("namespace") == "metacyc"
            and candidate["identity"].get("local_id") == declarations[local_id]
        ]
        if len(candidates) != 1:
            raise ValueError(f"Declared alias has no unique state: {local_id}")
        target, native = candidates[0]
        if node["identity_status"] != native["identity_status"] or any(
            identity.get(field) != native["identity"].get(field)
            for field in ("chemical_state", "compartment")
        ):
            raise ValueError("Alias would change the chemical state or compartment")
        resolved[historical] = target
    return resolved
