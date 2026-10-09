"""Path counts and measurements on small graphs, including the worked example of Figure 2B."""

import math

import pytest

from biosynthesis.paths import (
    FIGURE_2B_EXAMPLE,
    compute_dp_metrics,
    node_measurements,
    normalize_result,
    psi_bits,
    target_graph,
)


def network(links: dict[str, list[str]]) -> dict:
    """A simplified graph: each product with one record per selected precursor."""
    return {
        product: [{"reactant": [p]} for p in precursors]
        for product, precursors in links.items()
    }


WORKED_EXAMPLE = FIGURE_2B_EXAMPLE


def test_worked_example_of_figure_2b():
    result = compute_dp_metrics(WORKED_EXAMPLE, "B")
    assert result["num_paths"] == 2
    assert psi_bits(result["num_paths"]) == 1.0
    assert result["mean_branching_depth_normalized"] == 0.5
    assert result["max_depth"] == 4  # endpoint T is four links from the target
    assert result["fraction_shared_80"] == pytest.approx(5 / 5.5)
    assert round(result["fraction_shared_80"], 2) == 0.91
    assert result["total_complexity"] == pytest.approx(
        2.0
    )  # two paths through n3, log2(2)


def test_node_measurements_agree_with_path_counter():
    graph = target_graph(WORKED_EXAMPLE, "B")
    ours = node_measurements(graph, "B")
    method = compute_dp_metrics(WORKED_EXAMPLE, "B")
    assert ours["num_paths"] == method["num_paths"]
    assert ours["lambda"] == pytest.approx(method["mean_branching_depth_normalized"])
    assert ours["sigma"] == pytest.approx(method["fraction_shared_80"])
    shared = {node["id"] for node in ours["nodes"] if node["q"] >= 0.8}
    assert shared == {"B", "n4", "n3", "n1", "T"}


def test_single_path():
    result = compute_dp_metrics(network({"C": ["B"], "B": ["A"]}), "C")
    assert result["num_paths"] == 1
    assert psi_bits(result["num_paths"]) == 0.0
    assert result["mean_branching_depth_normalized"] == 0.0
    assert result["fraction_shared_80"] == 1.0
    assert result["total_complexity"] == 0.0


@pytest.mark.parametrize(
    "links, target",
    [
        ({"B": ["A"]}, "A"),  # a precursor only, never a product
        ({"B": ["A"]}, "Z"),  # absent from the graph
        ({"A": ["A"]}, "A"),  # only a self-link
    ],
)
def test_target_without_a_precursor_has_zero_paths(links, target):
    result = compute_dp_metrics(network(links), target)
    assert result is None
    values = normalize_result(result)
    assert values["num_paths"] == 0.0
    assert psi_bits(values["num_paths"]) is None


def test_cycle_closing_link_is_removed():
    # A <- B, B <- A (closes a cycle), B <- C: one path C -> B -> A remains.
    result = compute_dp_metrics(network({"A": ["B"], "B": ["A", "C"]}), "A")
    assert result["num_paths"] == 1


def test_parallel_branches_multiply():
    # Two independent two-way branches in series give four paths.
    links = {
        "T": ["a", "b"],
        "a": ["m"],
        "b": ["m"],
        "m": ["c", "d"],
        "c": ["s"],
        "d": ["s"],
    }
    result = compute_dp_metrics(network(links), "T")
    assert result["num_paths"] == 4
    assert math.log2(result["num_paths"]) == 2.0
