"""Upstream tracing on small synthetic networks."""

from biosynthesis.graph import (
    BROADER,
    build_network,
    node_label,
    traced_roots,
    without_reactions,
)
from biosynthesis.tracing import (
    NO_PRODUCER,
    PRODUCER_EXCLUDED,
    classify_frontier,
    discover,
)


def node(name, compartment=None, status="exact_state", namespace="metacyc"):
    return {
        "node_id": "n:" + name + ("@" + compartment if compartment else ""),
        "identity": {
            "namespace": namespace,
            "local_id": name,
            "compartment": compartment,
            "chemical_state": None,
        },
        "identity_status": status,
    }


def edge(name, inputs, outputs, **extra):
    return {
        "edge_id": "e:" + name,
        "reaction_local_id": name,
        "source_occurrence_id": "o:" + name,
        "orientation": "left_to_right",
        "inputs": ["n:" + x for x in inputs],
        "outputs": ["n:" + x for x in outputs],
        **extra,
    }


def cyclic_network(rejected=()):
    nodes = [node(x) for x in "ABCD"]
    edges = [edge("R1", "BC", "A"), edge("R2", "D", "B"), edge("R3", "A", "D")]
    return build_network(BROADER, nodes, edges, rejected)


def test_tracing_follows_every_reactant_to_closure():
    result = discover(cyclic_network(), ["n:A"])
    assert result["required_nodes"] == ["n:A", "n:B", "n:C", "n:D"]
    assert result["edge_ids"] == ["e:R1", "e:R2", "e:R3"]
    assert result["complete"] and result["pending_nodes"] == []


def test_supplied_annotations_do_not_stop_tracing():
    network = cyclic_network()
    plain = discover(network, ["n:A"])
    annotated = discover(network, ["n:A"], supplied_annotations=["n:B"])
    assert plain["required_nodes"] == annotated["required_nodes"]
    assert plain["edge_ids"] == annotated["edge_ids"]


def test_frontier_reports_cycles_and_endpoints():
    network = cyclic_network()
    frontier = classify_frontier(network, discover(network, ["n:A"]))
    assert len(frontier["cyclic_components"]) == 1
    assert frontier["cyclic_components"][0]["compound_nodes"] == ["n:A", "n:B", "n:D"]
    assert frontier["endpoints"] == [{"node_id": "n:C", "status": NO_PRODUCER}]


def test_endpoint_with_a_rejected_producer_is_marked():
    rejected = [{"known_product_nodes": ["n:C"], "reasons": ["not accepted"]}]
    network = cyclic_network(rejected)
    frontier = classify_frontier(network, discover(network, ["n:A"]))
    assert frontier["endpoints"] == [{"node_id": "n:C", "status": PRODUCER_EXCLUDED}]


def test_unidentified_participants_are_counted():
    nodes = [node("A"), node("X", status="unresolved_identity")]
    edges = [edge("R1", "X", "A", unresolved_requirements=[{"role": "donor"}])]
    network = build_network(BROADER, nodes, edges)
    frontier = classify_frontier(network, discover(network, ["n:A"]))
    assert len(frontier["unresolved_nodes"]) == 1
    assert len(frontier["unresolved_edge_requirements"]) == 1


def test_removed_reaction_leaves_the_network_and_is_recorded():
    network = cyclic_network()
    controlled, removed = without_reactions(network, {"R2": "test"})
    assert list(removed) == ["e:R2"]
    assert "n:B" not in controlled.producers
    assert discover(controlled, ["n:A"])["edge_ids"] == ["e:R1"]
    assert controlled.rejected[-1]["known_product_nodes"] == ["n:B"]


def test_alias_traces_the_state_of_the_declared_compound():
    nodes = [
        node("LIT-X", namespace="literature"),
        node("CPD-X"),
        node("CPD-X", "CCO-IN"),
    ]
    network = build_network(BROADER, nodes, [edge("R1", ["CPD-X@CCO-IN"], ["CPD-X"])])
    roots = traced_roots(network, {"LIT-X": "CPD-X"}, ["n:LIT-X", "n:CPD-X"])
    assert roots == {"n:LIT-X": "n:CPD-X", "n:CPD-X": "n:CPD-X"}
    assert node_label(nodes[2]) == "CPD-X[CCO-IN]"
