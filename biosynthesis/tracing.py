"""Upstream tracing: every producing reaction and every reactant, back to closure."""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Iterable

import networkx as nx

from biosynthesis.graph import Network
from biosynthesis.records import content_id

NO_PRODUCER = "no_producer"
PRODUCER_EXCLUDED = "producer_excluded"


def discover(
    network: Network, roots: Iterable[str], *, supplied_annotations: Iterable[str] = ()
) -> dict:
    """Trace all producing reactions of the roots and all their reactants until nothing new is found.

    Annotations given in `supplied_annotations` are recorded but never used to stop tracing.
    """
    if isinstance(roots, (str, bytes)):
        raise ValueError("roots must be a collection of node identifiers")
    selected = sorted(set(roots))
    if not selected or any(key not in network.nodes for key in selected):
        raise ValueError("Every root must be a state of the network")
    annotations = sorted(set(supplied_annotations))
    required, stored, expanded, edge_ids = set(selected), set(selected), set(), set()
    queue = deque(selected)
    while queue:
        compound = queue.popleft()
        expanded.add(compound)
        for edge_id in network.producers.get(compound, ()):
            if edge_id in edge_ids:
                continue
            edge = network.edges[edge_id]
            edge_ids.add(edge_id)
            stored.update(edge["outputs"])
            for precursor in sorted(set(edge["inputs"])):
                stored.add(precursor)
                if precursor not in required:
                    required.add(precursor)
                    queue.append(precursor)
    return {
        "roots": selected,
        "required_nodes": sorted(required),
        "stored_nodes": sorted(stored),
        "expanded_nodes": sorted(expanded),
        "pending_nodes": sorted(queue),
        "edge_ids": sorted(edge_ids),
        "complete": not queue,
        "supplied_annotations": annotations,
    }


def classify_frontier(network: Network, result: dict) -> dict:
    """Find the cycles, the compounds without a producer and the unidentified participants of a trace.

    Cycles are the strongly connected components of the compound-reaction dependency
    graph. A compound without an accepted producer is an endpoint; it is marked as
    having a producer excluded by the criteria when an occurrence producing it was
    not accepted.
    """
    required = set(result["required_nodes"])
    expanded = set(result["expanded_nodes"])
    dependency = nx.DiGraph()
    dependency.add_nodes_from(("compound", key) for key in required)
    unresolved = []
    for key in result["edge_ids"]:
        edge = network.edges[key]
        for product in set(edge["outputs"]) & required:
            dependency.add_edge(("compound", product), ("reaction", key))
        for precursor in set(edge["inputs"]):
            dependency.add_edge(("reaction", key), ("compound", precursor))
        if edge.get("unresolved_requirements"):
            unresolved.append(
                {"edge_id": key, "requirements": edge["unresolved_requirements"]}
            )
    cyclic = []
    for component in nx.strongly_connected_components(dependency):
        if len(component) > 1:
            cyclic.append(
                {
                    "compound_nodes": sorted(
                        key for kind, key in component if kind == "compound"
                    ),
                    "edge_ids": sorted(
                        key for kind, key in component if kind == "reaction"
                    ),
                }
            )
    cyclic.sort(key=content_id)
    excluded = defaultdict(list)
    for row in network.rejected:
        if not row.get("reasons"):
            continue
        for key in set(row.get("known_product_nodes", [])):
            excluded[key].append(row)
    endpoints, unresolved_nodes = [], []
    for key in sorted(required):
        status = network.nodes[key]["identity_status"]
        if status.startswith("unresolved"):
            unresolved_nodes.append({"node_id": key, "status": status})
        if key not in expanded or network.producers.get(key):
            continue
        endpoints.append(
            {
                "node_id": key,
                "status": PRODUCER_EXCLUDED if excluded[key] else NO_PRODUCER,
            }
        )
    return {
        "cyclic_components": cyclic,
        "endpoints": endpoints,
        "unresolved_nodes": unresolved_nodes,
        "unresolved_edge_requirements": unresolved,
    }
