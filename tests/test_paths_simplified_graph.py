"""Building the simplified graph: product records, corrections and precursor selection."""

from biosynthesis.simplified_graph import (
    HUB_LIMIT,
    add_literature_connections,
    orientations,
    precursor_links,
    product_records,
    select_precursors,
    supporting_reactions,
)


def test_orientations_follow_the_recorded_direction():
    left, right = ["A"], ["B"]
    record = {"LEFT": left, "RIGHT": right}
    assert orientations({**record, "REACTION-DIRECTION": ["LEFT-TO-RIGHT"]}) == [
        (left, right)
    ]
    assert orientations(
        {**record, "REACTION-DIRECTION": ["PHYSIOL-RIGHT-TO-LEFT"]}
    ) == [(right, left)]
    assert orientations({**record, "REACTION-DIRECTION": ["REVERSIBLE"]}) == [
        (left, right),
        (right, left),
    ]
    assert orientations(record) == []
    assert (
        orientations({**record, "REACTION-DIRECTION": ["IRREVERSIBLE-LEFT-TO-RIGHT"]})
        == []
    )


def test_product_records_split_by_product_and_ec_number():
    reactions = {
        "R1": {
            "REACTION-DIRECTION": ["LEFT-TO-RIGHT"],
            "LEFT": ["A", "B"],
            "RIGHT": ["C", "D"],
            "EC-NUMBER": ["EC-1.1.1.1", "EC-2.2.2.2"],
        },
        "R2": {
            "REACTION-DIRECTION": ["LEFT-TO-RIGHT"],
            "LEFT": ["A", "B"],
            "RIGHT": ["C"],
            "EC-NUMBER": ["EC-1.1.1.1"],
        },
        "R3": {"REACTION-DIRECTION": ["REVERSIBLE"], "LEFT": ["C"], "RIGHT": ["C"]},
    }
    pathways = {"P1": {"REACTION-LIST": ["R1", "R2", "R3"]}}
    products = product_records(reactions, pathways, ["P1"])
    assert products["C"] == [
        {"reactant": ["A", "B"], "ec_number": ["EC-1.1.1.1"]},
        {"reactant": ["A", "B"], "ec_number": ["EC-2.2.2.2"]},
    ]
    assert len(products["D"]) == 2


def test_heaviest_reactant_is_selected():
    products = {"P": [{"reactant": ["light", "heavy", "unknown"]}]}
    network, _ = select_precursors(products, {"light": 10.0, "heavy": 300.0})
    assert network["P"][0]["reactant"] == ["heavy"]


def test_ties_and_missing_weights_use_identifier_order():
    products = {"P": [{"reactant": ["b", "a"]}], "Q": [{"reactant": ["z", "y"]}]}
    network, _ = select_precursors(products, {"a": 5.0, "b": 5.0})
    assert network["P"][0]["reactant"] == ["a"]
    assert network["Q"][0]["reactant"] == ["y"]


def test_hubs_carriers_and_exemptions():
    products = {f"P{i}": [{"reactant": ["HUB", f"X{i}"]}] for i in range(HUB_LIMIT + 1)}
    products["Q"] = [{"reactant": ["FE+2"]}]
    masses = {"HUB": 1000.0, **{f"X{i}": 1.0 for i in range(HUB_LIMIT + 1)}}
    network, summary = select_precursors(products, masses)
    assert network["P0"][0]["reactant"] == ["X0"]
    assert "Q" not in network  # its only reactant is a carrier
    assert summary["hub_count"] == 1
    exempted, _ = select_precursors(products, masses, hub_exemptions=frozenset({"HUB"}))
    assert exempted["P0"][0]["reactant"] == ["HUB"]


def test_precursor_links_drop_self_links():
    network = {
        "A": [{"reactant": ["A"]}, {"reactant": ["B"]}],
        "C": [{"reactant": ["C"]}],
    }
    assert precursor_links(network) == {"A": ["B"]}


def test_corrections_and_literature_connections():
    reactions = {
        "R1": {"REACTION-DIRECTION": ["LEFT-TO-RIGHT"], "LEFT": ["A"], "RIGHT": ["B"]},
        "R2": {"REACTION-DIRECTION": ["LEFT-TO-RIGHT"], "LEFT": ["C"], "RIGHT": ["D"]},
        "R3": {"REACTION-DIRECTION": ["LEFT-TO-RIGHT"], "LEFT": ["C"], "RIGHT": ["D"]},
    }
    native = {
        "B": [{"reactant": ["A"], "ec_number": []}],
        "D": [{"reactant": ["C"], "ec_number": []}],
    }
    nodes = {
        "n1": {"identity": {"local_id": "E"}},
        "n2": {"identity": {"local_id": "F"}},
        "n3": {"identity": {"local_id": "A"}},
        "n4": {"identity": {"local_id": "B"}},
    }
    edge = {"reaction_local_id": "L1", "source_occurrence_id": "o", "edge_id": "e"}
    changes = {
        "added_edges": [
            {**edge, "inputs": ["n1"], "outputs": ["n2"]},
            {**edge, "inputs": ["n3"], "outputs": ["n4"]},
        ]
    }
    products = add_literature_connections(
        native, {"R1", "R2"}, supporting_reactions(reactions), changes, nodes
    )
    assert products["B"] == [
        {
            "reactant": ["A"],
            "ec_number": [],
            "_source": "literature_connection",
            "_literature_evidence": [
                {"reaction_id": "L1", "source_occurrence_id": "o", "edge_id": "e"}
            ],
        }
    ]  # R1 excluded, then restored by a literature connection
    assert products["D"] == [
        {"reactant": ["C"], "ec_number": []}
    ]  # R3 still supports it
    assert products["F"][0]["reactant"] == ["E"]
