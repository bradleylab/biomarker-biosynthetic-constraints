"""Draw Figure 3: overlap of recorded enzymes between biomarker pairs.

Panel A gives the outcomes for pairs in which both members have at least one
enzyme record, under strict and broader criteria. Panel B draws the two worked
examples (beta-amyrin with alpha-amyrin, diploptene with tetrahymanol) twice: in
the upper row with circle areas proportional to the number of enzyme records, and
in the lower row proportional to the number of projected paths in the simplified
graph with conditional connections. The organism of each enzyme is the COMMON-NAME
in classes.dat of its SPECIES in proteins.dat. Writes Figure3.pdf, Figure3.png and
Figure3_values.json (the plotted values).

Usage:
  python analyses/09_figures/figure3.py --results results --metacyc DIR --out results/figures
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from biosynthesis import enzymes
from biosynthesis.paths import target_graph
from biosynthesis.simplified_graph import select_precursors
from biosynthesis.tables import read_csv, write_json

PAIRS = "pairs"
# The lower row uses the graph with conditional connections, which includes the
# proposed bacterial conversion of diploptene to tetrahymanol (Table S8).
PATH_GRAPH = "projected_paths/simplified_graph"
PATH_VERSION = "with_conditional"
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 8,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    }
)
# Okabe-Ito colors, as in Figures 1 and 2.
BLUE, ORANGE, SKY = "#0072B2", "#E69F00", "#56B4E9"
LIGHT = "#D9D9D9"
TEXT_GRAY = "#555555"
MM = 1 / 25.4
WIDTH = 183 * MM
CRITERIA = (("strict", "Strict criteria"), ("broader", "Broader criteria"))
OUTCOMES = {
    "No shared enzyme": "no_shared_enzyme",
    "Overlap, not identical": "overlap_not_identical",
    "Identical sets": "identical_sets",
}
MISSING = "no_record_for_one_or_both"
EXAMPLES = (("CPD-6948", "CPD-8250"), ("HOP-2229-ENE", "CPD-11533"))
DISPLAY = {
    "CPD-6948": "β-Amyrin",
    "CPD-8250": "α-Amyrin",
    "HOP-2229-ENE": "Diploptene",
    "CPD-11533": "Tetrahymanol",
}
# Shared precursor named in the lower row.
PRECURSOR_DISPLAY = {"EPOXYSQUALENE": "oxidosqualene"}
# Descriptions of the enzymes named in panel B, keyed by the producing reactions
# they catalyze. RXN-7570 and RXN-8434 cyclize 2,3-oxidosqualene to beta- and
# alpha-amyrin, so an enzyme recorded for both is a multifunctional oxidosqualene
# cyclase; 5.4.99.17-RXN is squalene-hopene cyclase (EC 5.4.99.17); RXN-12337 is the
# squalene-tetrahymanol cyclase reaction of Tetrahymena (MetaCyc names the enzyme
# tetrahymanol synthase). Another reaction set raises KeyError.
DESCRIPTIONS = {
    ("RXN-7570", "RXN-8434"): (
        "multifunctional oxidosqualene cyclase",
        "multifunctional oxidosqualene cyclases",
    ),
    ("5.4.99.17-RXN",): ("squalene–hopene cyclase", "squalene–hopene cyclases"),
    ("RXN-12337",): ("squalene–tetrahymanol cyclase", "squalene–tetrahymanol cyclases"),
}


def panel_label(ax, letter: str, title: str, y: float = 1.0) -> None:
    ax.text(
        0, y, letter, transform=ax.transAxes, fontsize=11, weight="bold", va="bottom"
    )
    ax.text(
        0.045, y, title, transform=ax.transAxes, fontsize=9, weight="bold", va="bottom"
    )


def pair_outcomes(rows: list[dict], criteria: str) -> dict[str, int]:
    """Pairs with no record for one or both, no shared enzyme, overlap and identical sets."""
    counts = {name: 0 for name in (MISSING, *OUTCOMES.values())}
    for row in rows:
        a, b, shared = (
            int(row[f"{criteria}_{key}_enzyme_count"])
            for key in ("left", "right", "shared")
        )
        if a == 0 or b == 0:
            counts[MISSING] += 1
        elif shared == 0:
            counts["no_shared_enzyme"] += 1
        elif a == shared and b == shared:
            counts["identical_sets"] += 1
        else:
            counts["overlap_not_identical"] += 1
    return counts


def draw_outcomes(ax, rows: list[dict]) -> dict:
    colors = {
        "No shared enzyme": LIGHT,
        "Overlap, not identical": SKY,
        "Identical sets": BLUE,
    }
    summary = {}
    widest = 0
    for i, (criteria, label) in enumerate(CRITERIA):
        counts = pair_outcomes(rows, criteria)
        summary[criteria] = counts
        y, left = len(CRITERIA) - 1 - i, 0
        for category, color in colors.items():
            n = counts[OUTCOMES[category]]
            ax.barh(
                y, n, left=left, height=0.55, color=color, edgecolor="black", lw=0.6
            )
            left += n
        widest = max(widest, left)
        evaluated = left
        ax.text(
            -30, y, f"{label}\n{evaluated:,} pairs", ha="right", va="center", fontsize=8
        )
        shared = counts["overlap_not_identical"] + counts["identical_sets"]
        ax.text(
            counts["no_shared_enzyme"] / 2,
            y,
            f"{counts['no_shared_enzyme']:,}",
            ha="center",
            va="center",
            fontsize=7.5,
        )
        ax.text(
            evaluated + 40,
            y,
            f"{shared:,} share an enzyme\n({counts['identical_sets']} identical)",
            ha="left",
            va="center",
            fontsize=7.5,
        )
    ax.set_xlim(-10, widest * 1.32)
    ax.set_ylim(-0.6, len(CRITERIA) - 0.25)
    ax.set_yticks([])
    ax.set_xlabel(
        "Biomarker pairs in which both members have an enzyme record", fontsize=8
    )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(labelsize=7)
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black", lw=0.6)
        for c in colors.values()
    ]
    ax.legend(
        handles,
        list(colors),
        loc="upper right",
        frameon=False,
        fontsize=7,
        ncol=1,
        bbox_to_anchor=(1.0, 1.12),
    )
    return summary


def bits_text(n_from: int, n_to: int) -> str:
    value = math.log2(n_from / n_to)
    if value.is_integer():
        return f"{value:.0f} bit" if value == 1 else f"{value:.0f} bits"
    return f"{value:.2f} bits"


def plural(n: int) -> str:
    return "enzyme" if n == 1 else "enzymes"


def italic(text: str) -> str:
    """Italic Arial through mathtext, with hyphens kept outside the italic segments."""
    parts = [part.replace(" ", r"\;") for part in text.split("-")]
    return "-".join(rf"$\mathit{{{part}}}$" for part in parts)


def abbreviated(name: str) -> str:
    """Binomial with the genus abbreviated, in italics."""
    genus, epithet = name.split(" ", 1)
    return italic(f"{genus[0]}. {epithet}")


def enzyme_organisms(found: set[str], index: enzymes.EnzymeIndex) -> dict[str, str]:
    """Organism name of each enzyme: its one SPECIES in proteins.dat, named in classes.dat."""
    names = {}
    for enzyme in sorted(found):
        taxa = index.proteins[enzyme].get("SPECIES", [])
        if len(taxa) != 1:
            raise ValueError(f"expected one recorded SPECIES for {enzyme}")
        names[enzyme] = enzymes.display_name(index.classes[taxa[0]]["COMMON-NAME"][0])
    return names


def described(found: set[str], reactions: set[str], index: enzymes.EnzymeIndex) -> dict:
    """Count, description and organisms of a set of enzymes."""
    singular, many = DESCRIPTIONS[tuple(sorted(reactions))]
    organisms = list(dict.fromkeys(enzyme_organisms(found, index).values()))
    return {
        "count": len(found),
        "description": singular if len(found) == 1 else many,
        "organisms": organisms,
        "organism_text": ", ".join(abbreviated(name) for name in organisms),
    }


def circle_titles(ax, first: str, second: str, a: set, b: set, layout: dict) -> None:
    style = {"ha": "center", "fontsize": 8, "weight": "bold"}
    ax.text(
        layout["first_title_x"],
        1.1,
        f"{DISPLAY[first]}\n{len(a)} {plural(len(a))}",
        va="bottom",
        color=BLUE,
        **style,
    )
    ax.text(
        layout["second_title_x"],
        layout["second_title_y"],
        f"{DISPLAY[second]}\n{len(b)} {plural(len(b))}",
        va="bottom",
        color=ORANGE,
        **style,
    )


def euler(ax, pair: tuple[str, str], sets: dict, reactions: dict, index) -> dict:
    """Two-set overlap diagram with circle areas proportional to enzyme counts."""
    first, second = pair
    a, b = sets[first], sets[second]
    shared = a & b
    r2 = math.sqrt(len(b) / len(a))
    disjoint = not shared
    c2 = 1.0 + 0.25 + r2 if disjoint else 1.15
    ax.add_patch(
        Circle((0, 0), 1.0, facecolor=BLUE, alpha=0.18, edgecolor=BLUE, lw=1.2)
    )
    ax.add_patch(
        Circle((c2, 0), r2, facecolor=ORANGE, alpha=0.22, edgecolor=ORANGE, lw=1.2)
    )
    count = {"ha": "center", "va": "center", "fontsize": 10, "weight": "bold"}
    record = {first: len(a), second: len(b), "shared": sorted(shared)}
    if disjoint:
        ax.text(0, 0, str(len(a)), **count)
        ax.text(c2, 0, str(len(b)), **count)
        circle_titles(
            ax,
            first,
            second,
            a,
            b,
            {"first_title_x": 0, "second_title_x": c2, "second_title_y": r2 + 0.1},
        )
        sides = {n: described(sets[n], reactions[n], index) for n in pair}
        ax.text(
            0.55,
            -1.2,
            "No shared enzyme:\nno recorded enzyme makes both",
            ha="center",
            va="top",
            fontsize=7,
        )
        lines = []
        for name, values in sides.items():
            text = f"{DISPLAY[name]}: {values['count']} {values['description']}"
            organisms = f"({values['organism_text']})"
            lines += (
                [f"{text} {organisms}"]
                if len(values["organisms"]) == 1
                else [text, organisms]
            )
        ax.text(
            0.55,
            -2.25,
            "\n".join(lines),
            ha="center",
            va="top",
            fontsize=7,
            color=TEXT_GRAY,
        )
        record["sides"] = {
            n: {k: v for k, v in s.items() if k != "organism_text"}
            for n, s in sides.items()
        }
    else:
        ax.text(-0.35, 0, str(len(a - b)), **count)
        ax.text(0.78, 0, str(len(shared)), **count)
        ax.text(c2 + r2 * 0.45, 0, str(len(b - a)), **count)
        circle_titles(
            ax,
            first,
            second,
            a,
            b,
            {
                "first_title_x": -0.25,
                "second_title_x": c2 + 0.35,
                "second_title_y": 1.1,
            },
        )
        common = described(shared, reactions[first] | reactions[second], index)
        heading = "Shared enzyme:" if len(shared) == 1 else "Shared enzymes:"
        ax.text(
            0.55,
            -1.2,
            f"{heading}\n{common['count']} {common['description']}\n({common['organism_text']})",
            ha="center",
            va="top",
            fontsize=7,
        )
        reduction = (
            "Assuming a shared enzyme:\n"
            f"{DISPLAY[first]}: from {len(a)} to {len(shared)} ({bits_text(len(a), len(shared))})\n"
            f"{DISPLAY[second]}: from {len(b)} to {len(shared)} ({bits_text(len(b), len(shared))})"
        )
        ax.text(
            0.55, -2.25, reduction, ha="center", va="top", fontsize=7, color=TEXT_GRAY
        )
        record["shared_description"] = {
            k: v for k, v in common.items() if k != "organism_text"
        }
    ax.set_xlim(-1.6, 2.7)
    ax.set_ylim(-3.2, 1.75)
    ax.set_aspect("equal")
    ax.axis("off")
    return record


def path_set(dag: dict[str, list[str]], target: str) -> set[tuple[str, ...]]:
    """Every projected path of a target, listed from the target to an endpoint."""
    found = set()

    def walk(path: tuple[str, ...]) -> None:
        precursors = dag[path[-1]]
        if not precursors:
            found.add(path)
        for precursor in precursors:
            walk(path + (precursor,))

    walk((target,))
    return found


def path_comparison(network: dict, pair: tuple[str, str]) -> dict:
    """Projected paths of each member of a pair and the paths they share.

    If one member lies on the other's paths, the shared paths are the product's
    paths through it, which are the paths left if the product was made from the
    other member. Otherwise the shared paths are those identical upstream of the
    two members.
    """
    paths = {c: path_set(target_graph(network, c), c) for c in pair}
    counts = {c: len(paths[c]) for c in pair}
    for product, precursor in (pair, pair[::-1]):
        through = {p for p in paths[product] if precursor in p[1:]}
        if through:
            precursor_paths = {p[p.index(precursor) :] for p in through}
            return {
                "relation": "precursor_and_product",
                "product": product,
                "precursor": precursor,
                "counts": counts,
                "remaining": {product: len(through), precursor: len(precursor_paths)},
            }
    upstream = [{p[1:] for p in paths[c]} for c in pair]
    shared = upstream[0] & upstream[1]
    if shared != upstream[0] or shared != upstream[1]:
        raise ValueError(f"{pair}: expected identical paths upstream of both members")
    (precursor,) = {p[0] for p in shared}
    return {
        "relation": "same_upstream_paths",
        "precursor": precursor,
        "counts": counts,
        "remaining": dict(counts),
    }


def path_circles(ax, pair: tuple[str, str], comparison: dict) -> None:
    """Two-set diagram of projected paths, with circle areas proportional to path counts."""
    first, second = pair
    colors = {first: (BLUE, 0.18), second: (ORANGE, 0.22)}
    counts, remaining = comparison["counts"], comparison["remaining"]
    count = {"ha": "center", "va": "center", "fontsize": 10, "weight": "bold"}
    title = {"ha": "center", "va": "bottom", "fontsize": 8, "weight": "bold"}

    def circle(name: str, x: float, r: float) -> None:
        color, alpha = colors[name]
        ax.add_patch(
            Circle((x, 0), r, facecolor=color, alpha=alpha, edgecolor=color, lw=1.2)
        )

    def heading(name: str, x: float, y: float) -> None:
        n = counts[name]
        text = f"{DISPLAY[name]}\n{n} {'path' if n == 1 else 'paths'}"
        ax.text(x, y, text, color=colors[name][0], **title)

    if comparison["relation"] == "precursor_and_product":
        product, precursor = comparison["product"], comparison["precursor"]
        r = math.sqrt(counts[precursor] / counts[product])
        inner = 1.0 - r - 0.08
        circle(product, 0, 1.0)
        circle(precursor, inner, r)
        ax.text(-0.6, 0, str(counts[product] - remaining[product]), **count)
        ax.text(inner, 0, str(remaining[product]), **count)
        heading(product, -0.25, 1.1)
        heading(precursor, 1.75, 0.15)
        shared_text = (
            f"Shared paths:\n{remaining[product]} of the {counts[product]} "
            f"{DISPLAY[product].lower()} paths pass through {DISPLAY[precursor].lower()}"
        )
        assumption = f"Assuming {DISPLAY[product].lower()} was made from {DISPLAY[precursor].lower()}:"
        order = (product, precursor)
    else:
        for name in pair:
            circle(name, 0.55, 1.0)
        ax.text(0.55, 0, str(remaining[first]), **count)
        heading(first, -0.05, 1.1)
        heading(second, 1.15, 1.1)
        shared_text = (
            f"Shared paths:\nthe same {remaining[first]} paths to "
            f"{PRECURSOR_DISPLAY[comparison['precursor']]}"
        )
        assumption = "Assuming one organism made both:"
        order = pair
    ax.text(0.55, -1.2, shared_text, ha="center", va="top", fontsize=7)
    lines = [assumption] + [
        f"{DISPLAY[n]}: from {counts[n]} to {remaining[n]} ({bits_text(counts[n], remaining[n])})"
        for n in order
    ]
    ax.text(
        0.55,
        -2.25,
        "\n".join(lines),
        ha="center",
        va="top",
        fontsize=7,
        color=TEXT_GRAY,
    )
    ax.set_xlim(-1.6, 2.7)
    ax.set_ylim(-3.2, 1.75)
    ax.set_aspect("equal")
    ax.axis("off")


def path_network(results: Path) -> dict:
    """The simplified graph with conditional connections, one precursor per record."""
    folder = results / PATH_GRAPH
    with open(folder / f"products_{PATH_VERSION}.json", encoding="utf-8") as handle:
        products = json.load(handle)
    with open(folder / "molecular_weights.json", encoding="utf-8") as handle:
        masses = json.load(handle)
    network, _ = select_precursors(products, masses)
    return network


def broader_sets(path: Path) -> tuple[dict, dict]:
    """Enzymes and producing reactions of each biomarker under broader criteria."""
    found, reactions = {}, {}
    for row in read_csv(path):
        if row["reaction_criteria"] == "broader":
            found[row["compound_id"]] = set(json.loads(row["enzymes"]))
            reactions[row["compound_id"]] = set(json.loads(row["producing_reactions"]))
    return found, reactions


def figure3(
    rows: list[dict], sets: dict, reactions: dict, index, network: dict, out: Path
) -> dict:
    height = 195
    fig = plt.figure(figsize=(WIDTH, height * MM))

    def box(x: float, bottom: float, w: float, h: float) -> tuple:
        return (x, bottom / height, w, h / height)

    ax_a = fig.add_axes(box(0.17, 150, 0.80, 30))
    panel_label(ax_a, "A", "Enzyme overlap between biomarker pairs", y=1.2)
    outcomes = draw_outcomes(ax_a, rows)
    fig.text(0.02, 139 / height, "B", fontsize=11, weight="bold")
    fig.text(0.056, 139 / height, "Worked examples", fontsize=9, weight="bold")
    row_title = {"ha": "center", "fontsize": 8.5}
    fig.text(0.5, 132 / height, "Enzyme records", **row_title)
    fig.text(
        0.5, 64 / height, "Projected paths, with conditional connections", **row_title
    )
    examples, paths = {}, {}
    for i, (first, second) in enumerate(EXAMPLES):
        x = 0.02 + 0.50 * i
        key = f"{first}|{second}"
        examples[key] = euler(
            fig.add_axes(box(x, 68, 0.46, 62)), (first, second), sets, reactions, index
        )
        paths[key] = path_comparison(network, (first, second))
        path_circles(fig.add_axes(box(x, 0, 0.46, 62)), (first, second), paths[key])
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out / "Figure3.pdf", dpi=300, metadata={"CreationDate": None, "ModDate": None}
    )
    fig.savefig(out / "Figure3.png", dpi=300)
    plt.close(fig)
    return {
        "outcomes": outcomes,
        "examples": examples,
        "paths": {"network_version": PATH_VERSION, "examples": paths},
    }


def main(args: argparse.Namespace) -> None:
    rows = read_csv(args.results / PAIRS / "TableS3_pair_comparisons.csv")
    sets, reactions = broader_sets(args.results / PAIRS / "biomarker_enzyme_sets.csv")
    values = figure3(
        rows,
        sets,
        reactions,
        enzymes.enzyme_index(args.metacyc),
        path_network(args.results),
        args.out,
    )
    write_json(args.out / "Figure3_values.json", values)
    print(json.dumps(values["outcomes"], indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True, help="results folder")
    parser.add_argument(
        "--metacyc", type=Path, required=True, help="MetaCyc 29.5 data folder"
    )
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    main(parser.parse_args())
