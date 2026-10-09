"""Bounded enumeration of inclusion-minimal reaction sets that account for a target compound.

A reaction set accounts for the target when every selected reaction contributes to
the target and every required compound with an accepted producer is produced by a
selected reaction. A set is inclusion-minimal when no reaction can be removed. The
search visits at most `MAX_SEARCH_STATES` states per query; a search that stops at
the limit reports only a lower bound on the number of sets.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction

import networkx as nx

from biosynthesis.graph import BROADER, STRICT, Network
from biosynthesis.records import content_id
from biosynthesis.tracing import discover

MAX_SEARCH_STATES = 100
# Fields of an unresolved requirement that distinguish reactions; its free-text
# description does not.
REQUIREMENT_FIELDS = ("role", "status", "local_id")

Row = tuple[str, ...]


@dataclass
class ReactionGraph:
    """The traced network of one target, with occurrences of the same reaction merged."""

    target: str
    criteria: str
    nodes: dict
    reactions: dict
    producers: dict


def _signature(network: Network, edge: dict) -> dict:
    """Describe a reaction occurrence by its participants and coefficients, so that repeats merge."""
    signature = {
        "tier": network.criteria,
        "requirements": sorted(
            (
                {key: item[key] for key in REQUIREMENT_FIELDS if key in item}
                for item in edge.get("unresolved_requirements", [])
            ),
            key=content_id,
        ),
    }
    if "supplementary_record" in edge:
        signature["semantics"] = "supplementary_declared_multiplicity"
        for side in ("inputs", "outputs"):
            signature[side] = sorted(Counter(edge[side]).items())
        return signature
    parsed = edge.get("participant_records", {}).get("resolved", {})
    oriented = (
        ("left", "right")
        if edge["orientation"] == "left_to_right"
        else ("right", "left")
    )
    try:
        for side, source_side in zip(("inputs", "outputs"), oriented):
            records = parsed[source_side]
            if network.criteria == BROADER:
                by_locator = {
                    record["occurrence_locator"]: record for record in records
                }
                if len(by_locator) != len(records):
                    raise ValueError("Duplicate coefficient locator")
                records = [
                    by_locator[record["field_occurrence_id"]]
                    for record in edge["participant_records"]["raw"][source_side]
                ]
            if len(records) != len(edge[side]):
                raise ValueError("Incomplete coefficient projection")
            coefficients = defaultdict(Fraction)
            for key, record in zip(edge[side], records):
                if record["coefficient"]["status"] != "exact":
                    raise ValueError("Unresolved coefficient")
                value = Fraction(record["coefficient"]["value"])
                if value <= 0:
                    raise ValueError("Nonpositive coefficient")
                coefficients[key] += value
            signature[side] = sorted(
                (key, str(value)) for key, value in coefficients.items()
            )
        signature["semantics"] = "recorded_exact_coefficients"
    except (KeyError, ValueError, TypeError, ZeroDivisionError):
        signature.update(
            semantics="unresolved_equivalence",
            occurrence_edge_id=edge["edge_id"],
            inputs=sorted(Counter(edge["inputs"]).items()),
            outputs=sorted(Counter(edge["outputs"]).items()),
        )
    return signature


def reaction_graph(network: Network, target: str) -> ReactionGraph:
    """Trace the target and merge reaction occurrences with identical participants and coefficients."""
    discovered = discover(network, [target])
    reactions = {}
    for edge_id in discovered["edge_ids"]:
        edge = network.edges[edge_id]
        signature = _signature(network, edge)
        key = content_id(signature)
        if key not in reactions:
            reactions[key] = {
                "signature": signature,
                "inputs": sorted(set(edge["inputs"])),
                "outputs": sorted(set(edge["outputs"])),
                "occurrence_edge_ids": [],
            }
        reactions[key]["occurrence_edge_ids"].append(edge_id)
    for row in reactions.values():
        row["occurrence_edge_ids"].sort()
    nodes = {key: network.nodes[key] for key in discovered["stored_nodes"]}
    producers = defaultdict(set)
    for key, row in reactions.items():
        for product in row["outputs"]:
            producers[product].add(key)
    return ReactionGraph(
        target,
        network.criteria,
        nodes,
        reactions,
        {key: tuple(sorted(values)) for key, values in producers.items()},
    )


@dataclass(frozen=True)
class FeasibilityIndex:
    """Participants of each reaction, compiled once for repeated tests of reaction subsets.

    A subset is structurally feasible when, after repeatedly disabling reactions
    with an input that has an accepted producer but no surviving selected producer,
    the roots (and the inputs of fixed reactions) still have a selected producer.
    Compounds without any producer stay unresolved and do not disable a reaction.
    """

    nodes: frozenset[str]
    producers: dict[str, frozenset[str]]
    inputs: dict[str, frozenset[str]]
    outputs: dict[str, tuple[str, ...]]

    @classmethod
    def compile(cls, graph: ReactionGraph) -> FeasibilityIndex:
        """Copy the dependency relation of a reaction graph."""
        outputs = {key: [] for key in graph.reactions}
        for compound, producers in graph.producers.items():
            if compound in graph.nodes:
                for key in set(producers):
                    outputs[key].append(compound)
        return cls(
            frozenset(graph.nodes),
            {key: frozenset(values) for key, values in graph.producers.items()},
            {key: frozenset(row["inputs"]) for key, row in graph.reactions.items()},
            {key: tuple(values) for key, values in outputs.items()},
        )

    def viable(
        self,
        allowed: frozenset[str],
        roots: frozenset[str],
        fixed: frozenset[str] = frozenset(),
    ) -> bool:
        """Test whether the roots survive when only the allowed reactions can produce compounds."""
        if not roots <= self.nodes or any(
            not self.producers.get(root, frozenset()) & allowed for root in roots
        ):
            return False
        must_survive = set(roots)
        for key in fixed:
            must_survive.update(self.inputs[key])
        required = set(must_survive)
        counts: dict[str, int] = defaultdict(int)
        consumers: dict[str, list[str]] = defaultdict(list)
        for key in allowed:
            required.update(self.inputs[key])
            for compound in self.inputs[key]:
                consumers[compound].append(key)
            for compound in self.outputs[key]:
                counts[compound] += 1
        invalid = (required - self.nodes) | {
            compound
            for compound in required
            if self.producers.get(compound) and not counts.get(compound)
        }
        if invalid & must_survive:
            return False
        pending = list(invalid)
        live = set(allowed)
        while pending:
            for key in consumers.get(pending.pop(), ()):
                if key not in live:
                    continue
                live.remove(key)
                for product in self.outputs[key]:
                    counts[product] -= 1
                    if not counts[product] and product not in invalid:
                        if product in must_survive:
                            return False
                        invalid.add(product)
                        pending.append(product)
        return True


@dataclass
class SearchState:
    """Depth-first search over reaction subsets: states left, states seen, sets found."""

    processed: int
    pending: list[Row]
    seen: set[Row]
    found: set[Row]

    @classmethod
    def start(cls) -> SearchState:
        """The search begins from the empty reaction set."""
        return cls(0, [()], set(), set())

    def advance(
        self,
        graph: ReactionGraph,
        roots: frozenset[str],
        fixed: frozenset[str],
        budget: int,
        feasibility: FeasibilityIndex,
    ) -> None:
        """Process up to `budget` search states in a fixed deterministic order."""
        for _ in range(budget):
            if not self.pending:
                break
            added = frozenset(self.pending.pop())
            self.processed += 1
            if not any(set(row) <= added for row in self.found):
                selected = fixed | added
                required = set(roots)
                for key in selected:
                    required.update(graph.reactions[key]["inputs"])
                missing = sorted(
                    c
                    for c in required
                    if graph.producers.get(c) and not set(graph.producers[c]) & selected
                )
                if not missing:
                    if feasibility.viable(selected, roots, fixed) and all(
                        not feasibility.viable(selected - {key}, roots, fixed)
                        for key in added
                    ):
                        self.found.add(tuple(sorted(added)))
                else:
                    for key in reversed(graph.producers[missing[0]]):
                        successor = tuple(sorted(added | {key}))
                        if successor not in self.seen:
                            self.seen.add(successor)
                            self.pending.append(successor)


def consumes_target(graph: ReactionGraph, selected: frozenset) -> bool:
    """Whether the target lies on a cycle of the selected reactions (a reaction consumes it)."""
    required = {graph.target} | {
        node for key in selected for node in graph.reactions[key]["inputs"]
    }
    dependency = nx.DiGraph()
    for key in selected:
        row = graph.reactions[key]
        for product in set(row["outputs"]) & required:
            dependency.add_edge(("compound", product), ("reaction", key))
        for reactant in row["inputs"]:
            dependency.add_edge(("reaction", key), ("compound", reactant))
    for component in nx.strongly_connected_components(dependency):
        if len(component) > 1 and ("compound", graph.target) in component:
            return True
    return False


def enumerate_sets(
    graph: ReactionGraph, *, max_states: int = MAX_SEARCH_STATES
) -> dict:
    """Enumerate inclusion-minimal reaction sets for the target, within the search-state limit."""
    if type(max_states) is not int or max_states <= 0:
        raise ValueError("max_states must be a positive integer")
    state = SearchState.start()
    state.advance(
        graph,
        frozenset([graph.target]),
        frozenset(),
        max_states,
        FeasibilityIndex.compile(graph),
    )
    found = sorted(state.found)
    exact = not state.pending
    return {
        "status": "exact" if exact else "open",
        "sets": len(found) if exact else None,
        "sets_found": len(found),
        "sets_consuming_target": sum(
            consumes_target(graph, frozenset(row)) for row in found
        ),
        "coefficients_complete": graph.criteria == STRICT
        and all(
            row["signature"]["semantics"] == "recorded_exact_coefficients"
            for row in graph.reactions.values()
        ),
        "search_states_examined": state.processed,
        "search_states_left": len(state.pending),
        "reaction_sets": [list(row) for row in found],
    }
