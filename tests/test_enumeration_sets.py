"""Bounded enumeration of inclusion-minimal reaction sets on small synthetic networks."""

from biosynthesis.enumeration import (
    MAX_SEARCH_STATES,
    FeasibilityIndex,
    enumerate_sets,
    reaction_graph,
)
from biosynthesis.graph import BROADER, STRICT, build_network


def node(name):
    return {
        "node_id": "n:" + name,
        "identity": {"namespace": "metacyc", "local_id": name, "compartment": None},
        "identity_status": "exact_state",
    }


def literature_edge(name, inputs, outputs, occurrence=None):
    """A literature connection: participants as declared, merged by participants only."""
    return {
        "edge_id": "e:" + (occurrence or name),
        "reaction_local_id": name,
        "source_occurrence_id": "o:" + (occurrence or name),
        "orientation": "left_to_right",
        "supplementary_record": {"id": name},
        "inputs": ["n:" + x for x in inputs],
        "outputs": ["n:" + x for x in outputs],
    }


def network(edges, criteria=BROADER):
    names = {x[2:] for e in edges for side in ("inputs", "outputs") for x in e[side]}
    return build_network(criteria, [node(x) for x in sorted(names)], edges)


def test_search_limit_is_one_hundred_states():
    assert MAX_SEARCH_STATES == 100


def test_two_alternative_producers_give_two_sets():
    net = network([literature_edge("R1", "B", "A"), literature_edge("R2", "C", "A")])
    result = enumerate_sets(reaction_graph(net, "n:A"))
    assert result["status"] == "exact"
    assert result["sets"] == result["sets_found"] == 2
    assert result["search_states_left"] == 0


def test_required_reactant_with_a_producer_must_be_produced():
    edges = [
        literature_edge("R1", "BC", "A"),
        literature_edge("R2", "D", "B"),
        literature_edge("R3", "E", "B"),
    ]
    result = enumerate_sets(reaction_graph(network(edges), "n:A"))
    assert result["status"] == "exact"
    assert (
        result["sets"] == 2
    )  # {R1, R2} and {R1, R3}; C has no producer and stays unresolved


def test_identical_occurrences_merge_into_one_reaction():
    edges = [
        literature_edge("R1", "B", "A", "first"),
        literature_edge("R1", "B", "A", "second"),
    ]
    graph = reaction_graph(network(edges), "n:A")
    assert len(graph.reactions) == 1
    assert enumerate_sets(graph)["sets"] == 1


def test_requirement_descriptions_do_not_enter_the_signature():
    def requirement(text):
        return [{"role": "donor", "status": "unresolved", "source_literal": text}]

    first = literature_edge("R1", "B", "A", "first")
    second = literature_edge("R1", "B", "A", "second")
    first["unresolved_requirements"] = requirement("one description")
    second["unresolved_requirements"] = requirement("another description")
    graph = reaction_graph(network([first, second]), "n:A")
    assert len(graph.reactions) == 1


def test_open_search_reports_a_lower_bound():
    edges = [literature_edge("R1", "BC", "A")]
    edges += [literature_edge(f"B{i}", f"X{i}", "B") for i in range(4)]
    edges += [literature_edge(f"C{i}", f"Y{i}", "C") for i in range(4)]
    result = enumerate_sets(reaction_graph(network(edges), "n:A"), max_states=3)
    assert result["status"] == "open"
    assert result["sets"] is None
    assert result["search_states_examined"] == 3
    assert result["search_states_left"] > 0


def test_cycle_through_the_target_is_counted():
    edges = [literature_edge("R1", "B", "A"), literature_edge("R2", "A", "B")]
    result = enumerate_sets(reaction_graph(network(edges), "n:A"))
    assert result["sets"] == 1
    assert result["sets_consuming_target"] == 1


def test_feasibility_removes_reactions_whose_input_cannot_be_produced():
    edges = [literature_edge("R1", "B", "A"), literature_edge("R2", "C", "B")]
    graph = reaction_graph(network(edges), "n:A")
    index = FeasibilityIndex.compile(graph)
    keys = {graph.reactions[k]["occurrence_edge_ids"][0]: k for k in graph.reactions}
    assert index.viable(frozenset(keys.values()), frozenset(["n:A"]))
    assert not index.viable(frozenset([keys["e:R1"]]), frozenset(["n:A"]))


def test_coefficients_are_complete_only_under_strict_criteria_with_exact_values():
    edges = [literature_edge("R1", "B", "A")]
    assert not enumerate_sets(reaction_graph(network(edges, STRICT), "n:A"))[
        "coefficients_complete"
    ]
    exact = dict(literature_edge("R1", "B", "A"))
    del exact["supplementary_record"]
    exact["participant_records"] = {
        "resolved": {
            "left": [{"coefficient": {"status": "exact", "value": "1"}}],
            "right": [{"coefficient": {"status": "exact", "value": "1"}}],
        }
    }
    assert enumerate_sets(reaction_graph(network([exact], STRICT), "n:A"))[
        "coefficients_complete"
    ]
