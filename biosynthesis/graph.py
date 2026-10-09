"""Reaction networks read from the network files: compound states, reactions and producers."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable

# Criteria names are stored in the network records; part of the record identifiers.
STRICT = "strict_admitted"
BROADER = "source_supported_direction"
CRITERIA = (STRICT, BROADER)

# Network versions as stored in the folder names of the network files.
VERSIONS = ("primary", "with_conditional")
FOLDERS = {
    (STRICT, "primary"): "strict_primary",
    (STRICT, "with_conditional"): "strict_with_conditional",
    (BROADER, "primary"): "broader_primary",
    (BROADER, "with_conditional"): "broader_with_conditional",
}
CRITERIA_NAMES = {STRICT: "strict", BROADER: "broader"}


@dataclass
class Network:
    """One reaction network: compound states, reactions, rejected occurrences and producers.

    `rejected` holds the rows of screened_reactions.jsonl (one per source reaction
    occurrence and orientation, with the reasons it was not accepted, if any).
    """

    criteria: str
    nodes: dict[str, dict]
    edges: dict[str, dict]
    rejected: tuple[dict, ...]
    producers: dict[str, tuple[str, ...]]


def build_network(
    criteria: str,
    nodes: Iterable[dict],
    edges: Iterable[dict],
    rejected: Iterable[dict] = (),
) -> Network:
    """Index compound states and reactions, and list the reactions that produce each compound."""
    by_node, by_edge = {}, {}
    producers = defaultdict(list)
    for row in nodes:
        by_node[row["node_id"]] = row
    for row in edges:
        identifier = row["edge_id"]
        if identifier in by_edge:
            continue
        by_edge[identifier] = row
        for product in set(row["outputs"]):
            producers[product].append(identifier)
    return Network(
        criteria,
        by_node,
        by_edge,
        tuple(rejected),
        {key: tuple(sorted(values)) for key, values in producers.items()},
    )


def read_jsonl(path: Path) -> Iterable[dict]:
    """Yield the records of a JSON-lines file."""
    with open(path, "rb") as handle:
        for line in handle:
            yield json.loads(line)


def load_network(networks: Path, criteria: str, version: str) -> Network:
    """Read one of the four final networks from the network folder."""
    folder = Path(networks) / FOLDERS[criteria, version]
    return build_network(
        criteria,
        read_jsonl(folder / "nodes.jsonl"),
        read_jsonl(folder / "edges.jsonl"),
        read_jsonl(folder / "screened_reactions.jsonl"),
    )


def load_base_network(networks: Path, criteria: str) -> Network:
    """Read the MetaCyc network of one set of criteria, before literature records and corrections."""
    folder = Path(networks) / "base_graph"
    prefix = CRITERIA_NAMES[criteria]
    return build_network(
        criteria,
        read_jsonl(folder / f"{prefix}.nodes.jsonl"),
        read_jsonl(folder / f"{prefix}.edges.jsonl"),
        read_jsonl(folder / f"{prefix}.screened_reactions.jsonl"),
    )


def identity_queries(networks: Path, identities: list[dict]) -> list[dict]:
    """The query states of every biomarker identity under each set of criteria, in file order.

    Each row gives the criteria, the identity record (from the biomarker identity
    data) and the network states that represent the identity (empty when it has none).
    """
    by_id = {row["identity_id"]: row for row in identities}
    by_key = {
        (
            row["biomarker_number"],
            tuple(sorted(row["metacyc_compounds"])),
            tuple(sorted(row["candidate_inchikeys"])),
        ): row
        for row in identities
    }
    queries = []
    for row in read_jsonl(Path(networks) / "base_graph" / "identity_accounting.jsonl"):
        variant = row["identity_variant"]
        if "identity_id" in variant:
            identity = by_id[variant["identity_id"]]
        else:
            identity = by_key[
                (
                    variant["panel_ordinal"],
                    tuple(sorted(variant["source_compound_ids"])),
                    tuple(sorted(variant["candidate_inchikeys"])),
                )
            ]
        queries.append(
            {
                "criteria": row["tier"],
                "identity": identity,
                "targets": list(row["candidate_node_ids"]),
            }
        )
    return queries


def unlocalized_states(nodes: dict, aliases: dict) -> dict[str, str]:
    """Map each compound identifier to its single state without compartment, adding declared aliases."""
    choices = {}
    for value in nodes.values():
        identity = value["identity"]
        literal = identity.get("local_id")
        if literal is None or identity.get("compartment") is not None:
            continue
        choices.setdefault(literal, set()).add(value["node_id"])
    resolved = {
        key: next(iter(values)) for key, values in choices.items() if len(values) == 1
    }
    for old, target in aliases.items():
        if old not in choices or target not in resolved:
            raise ValueError("Alias refers to an absent or ambiguous state: " + old)
        resolved[old] = resolved[target]
    return resolved


def traced_roots(network: Network, aliases: dict, targets: Iterable[str]) -> dict:
    """Return, for each query state, the state traced: itself, or the state of its declared alias."""
    literals = {value["identity"].get("local_id") for value in network.nodes.values()}
    active = {old: new for old, new in aliases.items() if old in literals}
    canonical = unlocalized_states(network.nodes, active)
    roots = {}
    for target in targets:
        literal = network.nodes[target]["identity"].get("local_id")
        roots[target] = canonical[literal] if literal in active else target
    return roots


def node_label(node: dict) -> str:
    """Readable label of a compound state: its identifier, with the compartment when there is one."""
    identity = node["identity"]
    label = identity.get("local_id") or node["node_id"]
    if identity.get("compartment"):
        label += "[" + identity["compartment"] + "]"
    return label


def without_reactions(
    network: Network, reasons: dict[str, str]
) -> tuple[Network, dict]:
    """Remove every occurrence of the named reactions and record it as not accepted."""
    removed = {
        key: edge
        for key, edge in network.edges.items()
        if edge["reaction_local_id"] in reasons
    }
    # "kept_edge_ids" and the other keys are those stored in the network records.
    rejected = [
        {
            **row,
            "kept_edge_ids": [
                key for key in row.get("kept_edge_ids", []) if key not in removed
            ],
        }
        for row in network.rejected
    ]
    rejected.extend(
        {
            "tier": network.criteria,
            "reaction_local_id": edge["reaction_local_id"],
            "known_product_nodes": edge["outputs"],
            "kept_edge_ids": [],
            "reasons": ["removed_in_test", reasons[edge["reaction_local_id"]]],
            "direction_unresolved_participant_nodes": [],
        }
        for edge in removed.values()
    )
    controlled = build_network(
        network.criteria,
        network.nodes.values(),
        (edge for key, edge in network.edges.items() if key not in removed),
        rejected,
    )
    return controlled, removed
