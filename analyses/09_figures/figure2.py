"""Figure 2B and 2C: projected graphs of diploptene and cholesterol and the worked example.

Extracts the target-specific acyclic graphs of diploptene (HOP-2229-ENE) and
cholesterol from the primary simplified graph, measures them, checks the values
against the path counter of Table S1, and draws each graph: black target, gray
endpoints, blue nodes on at least 80% of paths and orange nodes on fewer, with
precursors above their products. The worked example of panel B is drawn the
same way. Writes Figure2_diploptene.pdf, Figure2_cholesterol.pdf,
Figure2_worked_example.pdf and Figure2.json (graphs and values).

The layout layers nodes by their longest distance from the target, routes links
that skip layers through invisible waypoints, and orders each layer by
barycenter sweeps followed by adjacent swaps that reduce crossings (Sugiyama et
al., 1981).

  python analyses/09_figures/figure2.py --results results --out results/figures
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from biosynthesis.paths import (  # noqa: E402
    FIGURE_2B_EXAMPLE,
    SHARED_FRACTION,
    compute_dp_metrics,
    node_measurements,
    target_graph,
)
from biosynthesis.simplified_graph import select_precursors  # noqa: E402
from biosynthesis.tables import write_json  # noqa: E402

matplotlib.rcParams["pdf.fonttype"] = 42

TARGETS = {"HOP-2229-ENE": "Diploptene", "CHOLESTEROL": "Cholesterol"}
NODE_COLORS = {  # Okabe and Ito (2008) palette
    "target": "#000000",
    "endpoint": "#999999",
    "shared": "#0072b2",  # on at least 80% of paths
    "other": "#e69f00",  # on fewer than 80% of paths
}
POINTS_PER_INCH = 72
SWEEPS = 12
# Width, height and node diameter (points) of each drawing.
SIZES = {
    "Diploptene": (80, 110, 10),
    "Cholesterol": (175, 230, 7),
    "Worked example": (80, 110, 10),
}
TOLERANCE = 1e-9


def rgb(hex_color: str) -> tuple[float, float, float]:
    value = hex_color.lstrip("#")
    return tuple(int(value[i : i + 2], 16) / 255 for i in (0, 2, 4))


def layers(nodes: list[dict], edges: list[list[str]]) -> dict[str, int]:
    """Longest-path layering from the target, so every arrow points downward."""
    target = next(n["id"] for n in nodes if n["target"])
    layer = {target: 0}
    changed = True
    while changed:
        changed = False
        for product, precursor in edges:
            if product in layer and layer.get(precursor, -1) < layer[product] + 1:
                layer[precursor] = layer[product] + 1
                changed = True
    return layer


def route(edges: list[list[str]], layer: dict[str, int]) -> tuple[dict, dict]:
    """Split arrows that skip layers into one-layer segments through waypoints."""
    chains, level = {}, dict(layer)
    for product, precursor in edges:
        chain = [precursor]
        for step in range(layer[precursor] - 1, layer[product], -1):
            waypoint = f"~{precursor}>{product}@{step}"
            level[waypoint] = step
            chain.append(waypoint)
        chains[(product, precursor)] = chain + [product]
    return chains, level


def crossings(
    upper: list[str], lower: list[str], segments: list[tuple[str, str]]
) -> int:
    up = {n: i for i, n in enumerate(upper)}
    low = {n: i for i, n in enumerate(lower)}
    pairs = [(up[a], low[b]) for a, b in segments if a in up and b in low]
    return sum(
        1
        for i, (a1, b1) in enumerate(pairs)
        for a2, b2 in pairs[i + 1 :]
        if (a1 - a2) * (b1 - b2) < 0
    )


def order(chains: dict, level: dict[str, int]) -> dict[int, list[str]]:
    """Order the nodes of each layer to reduce crossings."""
    rows = defaultdict(list)
    for node, lvl in sorted(level.items()):
        rows[lvl].append(node)
    segments = [(c[i], c[i + 1]) for c in chains.values() for i in range(len(c) - 1)]
    neighbors = defaultdict(list)
    for a, b in segments:
        neighbors[a].append(b)
        neighbors[b].append(a)
    levels = sorted(rows)

    def index() -> dict[str, float]:
        return {n: i - (len(r) - 1) / 2 for r in rows.values() for i, n in enumerate(r)}

    for sweep in range(SWEEPS):
        position = index()
        for lvl in levels if sweep % 2 == 0 else reversed(levels):
            row = rows[lvl]
            key = {
                n: sum(position[m] for m in neighbors[n]) / len(neighbors[n])
                for n in row
            }
            row.sort(key=lambda n: (key[n], n))
            position = index()

    def local(lvl: int) -> int:
        total = 0
        if lvl + 1 in rows:
            total += crossings(rows[lvl + 1], rows[lvl], segments)
        if lvl - 1 in rows:
            total += crossings(rows[lvl], rows[lvl - 1], segments)
        return total

    improved = True
    while improved:
        improved = False
        for lvl in levels:
            row = rows[lvl]
            for i in range(len(row) - 1):
                before = local(lvl)
                row[i], row[i + 1] = row[i + 1], row[i]
                if local(lvl) < before:
                    improved = True
                else:
                    row[i], row[i + 1] = row[i + 1], row[i]
    return rows


def kind(node: dict) -> str:
    if node["target"]:
        return "target"
    if node["endpoint"]:
        return "endpoint"
    return "shared" if node["q"] >= SHARED_FRACTION else "other"


def draw(name: str, graph: dict, out: Path) -> None:
    width, height, diameter = SIZES[name]
    nodes, edges = graph["nodes"], graph["edges"]
    chains, level = route(edges, layers(nodes, edges))
    rows = order(chains, level)
    widest = max(len(r) for r in rows.values())
    levels = max(rows)
    margin = diameter / 2 + 1
    step_x = (width - 2 * margin) / max(widest - 1, 1)

    def point(node: str) -> tuple[float, float]:
        row = rows[level[node]]
        px = width / 2 + (row.index(node) - (len(row) - 1) / 2) * step_x
        py = margin + level[node] / levels * (height - 2 * margin)
        return px, py

    fig = plt.figure(figsize=(width / POINTS_PER_INCH, height / POINTS_PER_INCH))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")
    radius = diameter / 2
    for chain in chains.values():
        points = [point(n) for n in chain]
        for (x0, y0), (x1, y1) in pairwise(points[:-1]):
            ax.plot([x0, x1], [y0, y1], color="black", lw=0.8, solid_capstyle="round")
        start = points[-2]
        first_segment = len(points) == 2
        ax.annotate(
            "",
            xy=points[-1],
            xytext=start,
            arrowprops={
                "arrowstyle": "-|>,head_length=0.35,head_width=0.18",
                "color": "black",
                "lw": 0.8,
                "shrinkA": radius + 0.5 if first_segment else 0,
                "shrinkB": radius + 0.5,
            },
        )
    for node in nodes:
        px, py = point(node["id"])
        ax.add_patch(
            plt.Circle(
                (px, py),
                radius,
                facecolor=rgb(NODE_COLORS[kind(node)]),
                edgecolor="black",
                lw=1.0,
                zorder=3,
            )
        )
    fig.savefig(out, transparent=True, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def check(name: str, ours: dict, method: dict) -> None:
    """The extracted graph gives the same values as the path counter."""
    for key, theirs in (
        ("num_paths", "num_paths"),
        ("lambda", "mean_branching_depth_normalized"),
        ("sigma", "fraction_shared_80"),
    ):
        if abs(ours[key] - method[theirs]) > TOLERANCE:
            raise ValueError(f"{name}: {key} {ours[key]} differs from {method[theirs]}")


def main(results: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    graph_dir = results / "projected_paths/simplified_graph"
    with open(graph_dir / "products_primary.json", encoding="utf-8") as handle:
        products = json.load(handle)
    with open(graph_dir / "molecular_weights.json", encoding="utf-8") as handle:
        masses = json.load(handle)
    network, _ = select_precursors(products, masses)
    graphs = {}
    for compound, label in TARGETS.items():
        ours = node_measurements(target_graph(network, compound), compound)
        check(label, ours, compute_dp_metrics(network, compound))
        graphs[label] = {"compound_id": compound, **ours}
    example = node_measurements(target_graph(FIGURE_2B_EXAMPLE, "B"), "B")
    check("Worked example", example, compute_dp_metrics(FIGURE_2B_EXAMPLE, "B"))
    graphs["Worked example"] = {"compound_id": None, **example}
    for label, graph in graphs.items():
        draw(label, graph, out / f"Figure2_{label.lower().replace(' ', '_')}.pdf")
    write_json(out / "Figure2.json", {"graphs": graphs})
    for label, graph in graphs.items():
        print(
            label,
            graph["num_paths"],
            round(graph["lambda"], 3),
            round(graph["sigma"], 3),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.results, args.out)
