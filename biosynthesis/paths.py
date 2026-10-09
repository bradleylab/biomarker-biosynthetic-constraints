"""Projected paths in the simplified compound graph and their measurements.

For a target compound, the selected precursor links are followed upstream, links
that close a cycle are removed by a depth-first search that visits precursors in
identifier order, and paths from the endpoints to the target are counted by
dynamic programming. From the counts come path multiplicity (psi, bits), normalized
branch depth (lambda), fraction shared (sigma) and the cumulative branching score.
"""

from __future__ import annotations

import math
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

RAW_METRIC_FIELDS = (
    "total_complexity",
    "num_paths",
    "mean_branching_depth_normalized",
    "max_depth",
    "fraction_shared_50",
    "fraction_shared_80",
    "fraction_shared_90",
    "n_consensus_nodes_80",
)
# Nodes on at least this fraction of paths are shared.
SHARED_FRACTION = 0.80
# The worked example of Figure 2B: one precursor chain with a branch two steps
# below the target B.
FIGURE_2B_EXAMPLE = {
    "B": [{"reactant": ["n4"]}],
    "n4": [{"reactant": ["n3"]}],
    "n3": [{"reactant": ["n1"]}, {"reactant": ["n2"]}],
    "n2": [{"reactant": ["n1"]}],
    "n1": [{"reactant": ["T"]}],
}


def compute_dp_metrics(network: Mapping, compound: str) -> Optional[Dict]:
    """Path count, psi sum, lambda, depth and sigma of one target, or None without a path.

    The network maps each product to its records, each with a one-element
    'reactant' list holding the selected precursor. Depth is the shortest distance
    from the target, and the greatest depth includes the endpoints.
    """
    if compound not in network:
        return None

    # Build main-reactant graph (first reactant per reaction)
    main_reactant_graph = {}
    for product, reactions in network.items():
        main_reactants = set()
        for rxn in reactions:
            reactants = rxn.get("reactant", [])
            if reactants and reactants[0] != product:
                main_reactants.add(reactants[0])
        if main_reactants:
            main_reactant_graph[product] = sorted(main_reactants)

    if compound not in main_reactant_graph:
        return None

    # Find reachable nodes and remove links that close a cycle, giving an acyclic graph
    reachable = set()
    dag_edges = defaultdict(list)
    visited = set()
    rec_stack = set()

    def dfs_dag(node):
        if node in visited:
            return
        visited.add(node)
        rec_stack.add(node)
        reachable.add(node)

        for reactant in main_reactant_graph.get(node, []):
            if reactant in rec_stack:
                # A link back into the current search path closes a cycle: not kept
                pass
            else:
                dag_edges[node].append(reactant)
                if reactant not in visited:
                    dfs_dag(reactant)

        rec_stack.remove(node)

    dfs_dag(compound)

    # Use the acyclic version of the graph
    main_reactant_graph = dict(dag_edges)

    # Endpoints: nodes with no reachable precursor
    terminals = set()
    for node in reachable:
        if node not in main_reactant_graph or not main_reactant_graph[node]:
            terminals.add(node)
        else:
            has_reachable = any(
                m in reachable for m in main_reactant_graph.get(node, [])
            )
            if not has_reachable:
                terminals.add(node)

    if not terminals:
        return None

    # Topological order by depth-first search
    visited = set()
    order = []

    def dfs(node):
        if node in visited:
            return
        visited.add(node)
        for main in main_reactant_graph.get(node, []):
            if main in reachable:
                dfs(main)
        order.append(node)

    dfs(compound)

    # paths_from[node] = paths from node to compound
    paths_from = defaultdict(float)
    paths_from[compound] = 1.0

    for node in reversed(order):
        for main in main_reactant_graph.get(node, []):
            if main in reachable:
                paths_from[main] += paths_from[node]

    # paths_to[node] = paths from the endpoints to node
    paths_to = defaultdict(float)
    for t in terminals:
        paths_to[t] = 1.0

    for node in order:
        for main in main_reactant_graph.get(node, []):
            if main in reachable:
                paths_to[node] += paths_to[main]

    total_paths = paths_to[compound]

    if total_paths == 0:
        return None

    # Cumulative branching score: sum of log2(in-degree) weighted by paths through the node
    psi = 0.0
    for node in reachable:
        if node in main_reactant_graph:
            in_deg = len([m for m in main_reactant_graph[node] if m in reachable])
            if in_deg > 1:
                paths_through = paths_to[node] * paths_from[node]
                psi += paths_through * np.log2(in_deg)

    # Depth from compound; endpoints receive a depth too
    depth = {compound: 0.0}
    for node in reversed(order):
        if node != compound:
            products = [
                p
                for p, mains in main_reactant_graph.items()
                if p in reachable and node in mains
            ]
            if products:
                depth[node] = 1 + min(depth.get(p, 0) for p in products)

    max_depth = max(depth.values()) if depth else 1.0

    # lambda: weighted mean normalized depth of the branch points
    lambda_sum = 0.0
    lambda_weight = 0.0
    for node in reachable:
        if node in main_reactant_graph:
            in_deg = len([m for m in main_reactant_graph[node] if m in reachable])
            if in_deg > 1:
                paths_through = paths_to[node] * paths_from[node]
                norm_depth = depth.get(node, 0) / max_depth if max_depth > 0 else 0
                lambda_sum += paths_through * norm_depth
                lambda_weight += paths_through

    lambda_val = lambda_sum / lambda_weight if lambda_weight > 0 else 0.0

    # sigma at several thresholds
    def compute_sigma(threshold):
        consensus = sorted(
            [
                n
                for n in reachable
                if (paths_to[n] * paths_from[n]) / total_paths >= threshold
            ]
        )
        e_cons = sum((paths_to[n] * paths_from[n]) / total_paths for n in consensus)
        e_len = sum((paths_to[n] * paths_from[n]) / total_paths for n in reachable)
        return e_cons / e_len if e_len > 0 else 0.0, consensus

    sigma_50, cons_50 = compute_sigma(0.50)
    sigma_80, cons_80 = compute_sigma(0.80)
    sigma_90, cons_90 = compute_sigma(0.90)

    return {
        "total_complexity": psi,
        "avg_complexity": psi / total_paths if total_paths > 0 else 0,
        "num_paths": total_paths,
        "mean_branching_depth_normalized": lambda_val,
        "fraction_shared_50": sigma_50,
        "fraction_shared_80": sigma_80,
        "fraction_shared_90": sigma_90,
        "n_consensus_nodes_80": len(cons_80),
        "n_reachable": len(reachable),
        "n_terminals": len(terminals),
        "max_depth": max_depth,
        "consensus_nodes_80": cons_80[:20],
    }


def normalize_result(result: Mapping | None) -> dict[str, float]:
    """The measurements as floats; a target without a path gets zeros."""
    if result is None:
        return {field: 0.0 for field in RAW_METRIC_FIELDS}
    return {field: float(result[field]) for field in RAW_METRIC_FIELDS}


def psi_bits(num_paths: float) -> float | None:
    """Path multiplicity in bits, or None without a path."""
    return math.log2(num_paths) if num_paths > 0 else None


_NETWORKS: dict[str, Mapping] = {}


def _measure(target: str) -> tuple[str, dict]:
    return target, {
        name: normalize_result(compute_dp_metrics(network, target))
        for name, network in _NETWORKS.items()
    }


def measure_targets(
    networks: Mapping[str, Mapping], targets: Sequence[str], workers: int = 1
) -> dict[str, dict[str, dict[str, float]]]:
    """Measure every target in every network: {network name: {target: measurements}}.

    Targets are measured independently, so the number of worker processes does not
    change any value.
    """
    _NETWORKS.clear()
    _NETWORKS.update(networks)
    results: dict[str, dict] = {name: {} for name in networks}
    if workers > 1:
        with mp.get_context("fork").Pool(workers) as pool:
            outcomes: Iterable = pool.imap(_measure, targets, chunksize=8)
            for target, values in outcomes:
                for name, value in values.items():
                    results[name][target] = value
    else:
        for target in targets:
            for name, value in _measure(target)[1].items():
                results[name][target] = value
    _NETWORKS.clear()
    return results


def metric_frame(rows: Sequence[Mapping]) -> pd.DataFrame:
    """Rows of raw measurements as a frame with psi_bits, lambda and sigma added.

    Lambda and sigma are missing for compounds without a path, as is psi_bits.
    """
    frame = pd.DataFrame(rows)
    frame["psi_bits"] = frame["num_paths"].map(psi_bits)
    frame["cumulative_branching_score"] = frame["total_complexity"]
    frame["lambda"] = frame["mean_branching_depth_normalized"]
    frame["sigma"] = frame["fraction_shared_80"]
    frame.loc[frame.num_paths.eq(0), ["lambda", "sigma"]] = np.nan
    return frame


def target_graph(network: Mapping, target: str) -> dict[str, list[str]]:
    """The acyclic precursor graph of one target, as used by compute_dp_metrics."""
    selected = {}
    for product, reactions in network.items():
        mains = {
            r["reactant"][0]
            for r in reactions
            if r.get("reactant") and r["reactant"][0] != product
        }
        if mains:
            selected[product] = sorted(mains)
    dag: dict[str, list[str]] = defaultdict(list)
    visited, stack = set(), set()

    def visit(node: str) -> None:
        visited.add(node)
        stack.add(node)
        for precursor in selected.get(node, []):
            if precursor in stack:
                continue  # a link back into the current path closes a cycle
            dag[node].append(precursor)
            if precursor not in visited:
                visit(precursor)
        stack.remove(node)

    visit(target)
    return {node: dag.get(node, []) for node in visited}


def node_measurements(dag: Mapping[str, list[str]], target: str) -> dict:
    """Paths through each node, depths, and the path count, lambda and sigma of a graph."""
    order, seen = [], set()

    def topo(node: str) -> None:
        seen.add(node)
        for precursor in dag[node]:
            if precursor not in seen:
                topo(precursor)
        order.append(node)

    topo(target)
    endpoints = {n for n, precursors in dag.items() if not precursors}
    to_target = defaultdict(float, {target: 1.0})
    for node in reversed(order):
        for precursor in dag[node]:
            to_target[precursor] += to_target[node]
    from_endpoints = defaultdict(float, {n: 1.0 for n in endpoints})
    for node in order:
        for precursor in dag[node]:
            from_endpoints[node] += from_endpoints[precursor]
    total = from_endpoints[target]
    depth = {target: 0}
    for node in reversed(order):
        for precursor in dag[node]:
            depth[precursor] = min(
                depth.get(precursor, depth[node] + 1), depth[node] + 1
            )
    deepest = max(depth.values())
    through = {n: from_endpoints[n] * to_target[n] for n in dag}
    branches = [n for n in dag if len(dag[n]) > 1]
    weight = sum(through[n] for n in branches)
    lam = (
        sum(through[n] * depth[n] / deepest for n in branches) / weight
        if weight
        else 0.0
    )
    q = {n: through[n] / total for n in dag}
    sigma = sum(v for v in q.values() if v >= SHARED_FRACTION) / sum(q.values())
    return {
        "num_paths": total,
        "lambda": lam,
        "sigma": sigma,
        "nodes": [
            {
                "id": n,
                "depth": depth[n],
                "q": q[n],
                "endpoint": n in endpoints,
                "branch_point": n in branches,
                "target": n == target,
            }
            for n in sorted(dag, key=lambda n: (depth[n], n))
        ],
        "edges": [[n, p] for n in sorted(dag) for p in dag[n]],
    }
