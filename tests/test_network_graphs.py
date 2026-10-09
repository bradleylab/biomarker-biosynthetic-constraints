"""Base and final networks on small synthetic records."""

import pytest

from biosynthesis import network


def compound(local):
    return network.node(
        {
            "namespace": "metacyc",
            "local_id": local,
            "compartment": None,
            "chemical_state": None,
        },
        network.SOURCE_IDENTIFIER,
    )


def base_network():
    nodes = [compound(name) for name in ("A", "B", "C")]
    ids = {row["identity"]["local_id"]: row["node_id"] for row in nodes}
    edges = []
    for reaction, left, right in (("RXN-1", "A", "B"), ("RXN-2", "B", "C")):
        identity = {
            "source_occurrence_id": reaction,
            "tier": network.BROADER,
            "orientation": "left_to_right",
        }
        edges.append(
            {
                **identity,
                "edge_id": network.content_id(identity),
                "reaction_local_id": reaction,
                "inputs": [ids[left]],
                "outputs": [ids[right]],
            }
        )
    screened = [{"kept_edge_ids": [edge["edge_id"]], "reasons": []} for edge in edges]
    return network.Network.from_records(
        tier=network.BROADER, nodes=nodes, edges=edges, screened=screened
    ), ids


def config():
    record = {
        "record_id": "SUPP-RXN-X",
        "finding": "f",
        "sources": [{"doi": "10.1/x", "locator": "Fig. 1"}],
        "conversions": [
            {
                "id": "X-1",
                "connection_type": "supported",
                "inputs": ["C"],
                "outputs": ["SUPP-CPD-X"],
                "unresolved": [],
            },
            {
                "id": "X-2",
                "connection_type": "conditional",
                "inputs": ["SUPP-CPD-X"],
                "outputs": ["Y"],
                "unresolved": [],
            },
        ],
    }
    states = [
        {
            "namespace": network.LITERATURE_NAMESPACE,
            "local_id": name,
            "chemical_state": name,
            "identity_status": "s",
        }
        for name in ("SUPP-CPD-X", "Y")
    ]
    decisions = [
        {"native_id": "RXN-2", "action": "exclude", "finding": "g", "sources": []},
        {
            "native_id": "RXN-1",
            "action": "annotate",
            "finding": "h",
            "sources": [],
            "unresolved": [],
        },
    ]
    return {
        "records": [network._stored(record)],
        "states": states,
        "native_decisions": decisions,
        "participant_aliases": {"Z": "A"},
        "literature_evidence_on_native_reactions": [],
    }


def test_producers_are_indexed_by_output():
    base, ids = base_network()
    assert base.producers[ids["B"]] == (base.screened[0]["kept_edge_ids"][0],)


def test_primary_version_keeps_supported_connections_only():
    base, ids = base_network()
    primary, changes = network.build_final_network(base, config(), "primary")
    conditional, _ = network.build_final_network(base, config(), "with_conditional")
    assert [edge["reaction_local_id"] for edge in changes["added_edges"]] == ["X-1"]
    assert sorted(edge["reaction_local_id"] for edge in conditional.edges.values()) == [
        "RXN-1",
        "X-1",
        "X-2",
    ]
    assert "RXN-2" not in {edge["reaction_local_id"] for edge in primary.edges.values()}
    annotated = next(
        edge for edge in primary.edges.values() if edge["reaction_local_id"] == "RXN-1"
    )
    assert annotated["literature_annotations"][0]["finding"] == "h"
    replaced = [row for row in primary.screened if network.REPLACED in row["reasons"]]
    assert len(replaced) == 1 and replaced[0]["kept_edge_ids"] == []


def test_identifiers_do_not_depend_on_descriptive_text():
    base, _ = base_network()
    first, _ = network.build_final_network(base, config(), "primary")
    changed = config()
    changed["records"][0]["finding"] = "other"
    changed["states"][0]["chemical_state"] = "other description"
    second, _ = network.build_final_network(base, changed, "primary")
    assert set(first.edges) == set(second.edges)
    assert set(first.nodes) == set(second.nodes)
    renamed = config()
    renamed["records"][0]["conversions"][0]["id"] = "X-3"
    third, _ = network.build_final_network(base, renamed, "primary")
    assert set(first.edges) != set(third.edges)


def test_alias_resolves_to_target_node():
    base, ids = base_network()
    nodes = dict(base.nodes)
    extra = compound("Z")
    nodes[extra["node_id"]] = extra
    resolved = network.resolve_participants(nodes, {"Z": "A"})
    assert resolved["Z"] == ids["A"]
    with pytest.raises(ValueError):
        network.resolve_participants(nodes, {"Q": "A"})


def test_stored_form_drops_table_text_only():
    record = {
        "connection_type": "supported",
        "unresolved": [
            {
                "role": "donor",
                "local_id": "UNRESOLVED-X",
                "status": "unresolved",
                "source_literal": "s",
            }
        ],
        "table_reason": "only for the table",
    }
    stored = network._stored(record)
    assert "table_reason" not in stored
    assert stored["unresolved"] == record["unresolved"]
    assert stored["connection_type"] == "supported"
