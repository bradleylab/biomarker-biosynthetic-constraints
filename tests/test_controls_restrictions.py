"""Matched controls for upstream tracing on small synthetic networks."""

import base64

from biosynthesis.controls import (
    PATHWAY_ONLY,
    STOP_AT_STARTING,
    STRICT_REACTIONS,
    UNCONTROLLED,
    alias_roots,
    annotations_change_result,
    controlled_discovery,
    controls,
    membership,
    pathway_direction,
    pathway_ids,
    pathway_masks,
    starting_identifiers,
    stopped_states,
    strict_occurrences,
)
from biosynthesis.graph import BROADER, STRICT, build_network


def node(name, compartment=None, namespace="metacyc", status="exact_state"):
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


def edge(name, inputs, outputs, pathway=None, direction=("left_to_right",)):
    fields, observations = [], []
    if pathway:
        fields = [{"value_raw_base64": base64.b64encode(pathway.encode()).decode()}]
        observations = [
            {
                "source_family": "pathways.dat",
                "source_record_local_id": pathway,
                "interpreted_orientations": list(direction),
            }
        ]
    return {
        "edge_id": "e:" + name,
        "reaction_local_id": name,
        "source_occurrence_id": "o:" + name,
        "orientation": "left_to_right",
        "pathway_fields": fields,
        "direction_observations": observations,
        "inputs": ["n:" + x for x in inputs],
        "outputs": ["n:" + x for x in outputs],
    }


def chain(criteria=BROADER):
    """A <- B <- C <- D, with ATP required by the first reaction."""
    nodes = [node(x) for x in ("A", "B", "C", "D", "ATP")]
    edges = [
        edge("R1", ["B", "ATP"], ["A"], pathway="PWY-1"),
        edge("R2", ["C"], ["B"], pathway="PWY-1"),
        edge("R3", ["D"], ["C"]),
    ]
    return build_network(criteria, nodes, edges)


def test_stopping_keeps_the_stopped_compound_but_not_its_producers():
    result = controlled_discovery(chain(), "n:A", stop_nodes={"n:B"})
    assert result["edge_ids"] == ["e:R1"]
    assert result["required_nodes"] == ["n:A", "n:ATP", "n:B"]


def test_allowed_reactions_restrict_tracing():
    result = controlled_discovery(chain(), "n:A", allowed_edges={"e:R1", "e:R3"})
    assert result["edge_ids"] == ["e:R1"]


def test_pathway_masks_keep_direction_consistent_placements():
    network = chain()
    assert pathway_masks(network, {"PWY-1"}) == {"PWY-1": {"e:R1", "e:R2"}}
    assert pathway_masks(network, set()) == {}
    assert (
        pathway_direction(
            network.edges["e:R1"]["direction_observations"], "PWY-1", "right_to_left"
        )
        == "pathway_local_conflict_only"
    )


def test_pathway_identifiers_are_read_from_the_pathways_file(tmp_path):
    path = tmp_path / "pathways.dat"
    path.write_bytes(
        b"# comment\nUNIQUE-ID - PWY-1\nTYPES - X\n//\nUNIQUE-ID - PWY-2\n//\n"
    )
    assert pathway_ids(path) == {"PWY-1", "PWY-2"}


def test_controls_change_one_rule_at_a_time():
    network = chain()
    masks = pathway_masks(network, {"PWY-1"})
    stops = stopped_states(network, {("metacyc", "C")})
    strict_accepted = {("o:R1", "left_to_right")}
    result = controls(network, "n:A", stops, masks, strict_accepted)
    assert len(result[UNCONTROLLED]["edge_ids"]) == 3
    assert result[STOP_AT_STARTING]["edge_ids"] == ["e:R1", "e:R2"]
    assert result[PATHWAY_ONLY]["edge_ids"] == ["e:R1", "e:R2"]
    assert result[STRICT_REACTIONS]["edge_ids"] == ["e:R1"]
    assert list(result["pathway_contexts"]) == ["PWY-1"]
    assert STRICT_REACTIONS not in controls(
        chain(STRICT), "n:A", stops, masks, strict_accepted
    )


def test_unused_annotations_do_not_change_tracing():
    network = chain()
    expected = membership(controlled_discovery(network, "n:A"))
    assert not annotations_change_result(network, "n:A", expected, {"n:C"})


def test_starting_compounds_stop_every_state_of_compounds_in_the_strict_network():
    strict = build_network(STRICT, [node("ATP"), node("FAD")], [])
    identifiers = starting_identifiers(strict, ["ATP", "MENAQUINONE"])
    assert identifiers == {("metacyc", "ATP")}
    broader = build_network(
        BROADER, [node("ATP"), node("ATP", "CCO-IN"), node("MENAQUINONE")], []
    )
    assert stopped_states(broader, identifiers) == {"n:ATP", "n:ATP@CCO-IN"}


def test_strict_occurrences_pair_source_and_orientation():
    assert strict_occurrences(chain(STRICT)) == {
        ("o:R1", "left_to_right"),
        ("o:R2", "left_to_right"),
        ("o:R3", "left_to_right"),
    }


def test_alias_resolves_to_the_single_metacyc_state():
    nodes = [node("LIT-A", namespace="literature"), node("CPD-A")]
    network = build_network(BROADER, nodes, [])
    assert alias_roots(network, {"LIT-A": "CPD-A"}, {"n:LIT-A"}) == {
        "n:LIT-A": "n:CPD-A"
    }
