"""Figure 4: strict and broader assignments side by side for two geological samples.

For each sample, both cells repeat the same fossil-precursor-source diagram and
draw only the links used under that set of assignments; links used only under
broader assignments are styled by the kind of evidence behind them (tentative
precursor assignment, enzyme capacity, source hypothesis). The minimal
combinations of source categories are written beneath each diagram, each with
the broader evidence it needs. The two samples are those chosen by
select_examples.py. Text is set in Arial.

Writes Figure4.pdf, Figure4.png (300 dpi) and Figure4_values.json (the links,
combinations, labels and colors drawn).

Usage:
  python analyses/09_figures/figure4.py --results results --out results/figures

--results is the folder holding sources/, the output folder of the
source-compatibility analyses.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from matplotlib.textpath import TextToPath  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.sensitivity import zero_threshold  # noqa: E402
from biosynthesis.sources import BROADER, STRICT  # noqa: E402
from biosynthesis.tables import read_jsonl, write_json  # noqa: E402

matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 8,
        # Subscripts such as C$_{27}$ use Arial too, not matplotlib's default math font.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    }
)
GRAY = "#999999"
LIGHT = "#D9D9D9"
TEXT_GRAY = "#555555"
MM = 1 / 25.4  # inches per millimeter
PT_TO_MM = 25.4 / 72

FOSSIL = {
    "si1_aromatic_hopanoids": "Aromatic hopanoids",
    "si2_hopanes": "Hopanes",
    "si2_cholestanes": "Cholestanes (C$_{27}$)",
    "si2_ergostanes": "Ergostanes (C$_{28}$)",
    "si2_stigmastanes": "Stigmastanes (C$_{29}$)",
    "si2_cyclosteranes": "Cyclosteranes",
}
CONVENTIONAL_STERANES = ("si2_cholestanes", "si2_ergostanes", "si2_stigmastanes")
PRECURSOR = {
    "hopanoid_biological_precursor_pool": "Biological hopanoids",
    "c27_cholesterol_related_sterol_pool": "C$_{27}$ sterols",
    "c28_ergosterol_related_sterol_pool": "C$_{28}$ sterols",
    "c29_stigmasterol_related_sterol_pool": "C$_{29}$ sterols",
    "cycloartane_skeleton_sterol_pool": "Cycloartane-skeleton sterols",
    "lanosterol_related_sterol_pool": "Lanosterol-related sterols",
}
SOURCE = {
    "hopanoid_producing_bacterial_sources": "Hopanoid-producing\nbacteria",
    "bacterial_conventional_sterol_sources": "Bacteria producing\nconventional sterols",
    "crown_eukaryotic_crown_sterol_sources": "Crown eukaryotes producing\nconventional sterols",
    "bacterial_protosterol_or_ursterol_sources": "Bacteria producing\nprotosterols or ursterols",
    "stem_eukaryotic_protosterol_or_ursterol_hypotheses": "Hypothesized stem eukaryotes\nproducing protosterols or ursterols",
}
SOURCES = list(SOURCE)
VALUE_NOTE = {
    "reported_zero": "reported as zero",
    "blank": "blank in the workbook",
}

# Display names for the three kinds of broader-only evidence, keyed by the
# evidence_basis values of Table S6. Any other broader-only basis stops the script.
STRICT_KIND = "Strict assignment"
BROADER_KIND = {
    "primary_structural_and_diagenetic_interpretation": "Tentative precursor assignment",
    "functional_smt_assay_and_genomic_context": "Enzyme capacity",
    "explicit_fossil_source_hypothesis": "Source hypothesis",
}
KIND_ORDER = (STRICT_KIND, *BROADER_KIND.values())
# Line colors: Okabe-Ito orange and blue plus a mid gray, with distinct dash
# patterns, so no distinction rests on color alone; check_palette() reports the
# smallest CIELAB distance between them.
LINE_STYLE = {
    STRICT_KIND: {"color": "black", "ls": "-", "lw": 0.9},
    "Tentative precursor assignment": {
        "color": "#E69F00",
        "ls": (0, (3.2, 1.6)),
        "lw": 1.2,
    },
    "Enzyme capacity": {"color": "#0072B2", "ls": (0, (4.0, 1.4, 1.0, 1.4)), "lw": 1.2},
    "Source hypothesis": {"color": GRAY, "ls": (0, (1.0, 1.3)), "lw": 1.6},
}
FILL = {"fossil": "#eebf99", "precursor": "#bbe1f6", "source": "#f9f4b3"}
# Dichromacy simulation matrices of Machado, Oliveira and Fernandes (2009), full
# severity, applied in linear RGB.
SIMULATION = {
    "protanopia": [
        [0.152286, 1.052583, -0.204868],
        [0.114503, 0.786281, 0.099216],
        [-0.003882, -0.048116, 1.051998],
    ],
    "deuteranopia": [
        [0.367322, 0.860646, -0.227968],
        [0.280085, 0.672501, 0.047413],
        [-0.011820, 0.042940, 0.968881],
    ],
    "tritanopia": [
        [1.255528, -0.076749, -0.178779],
        [-0.078411, 0.930809, 0.147602],
        [0.004733, 0.691367, 0.303900],
    ],
}
MIN_LINE_COLOR_DISTANCE = (
    30.0  # CIELAB units between any two line colors, every simulation
)
TIER_LABEL = {STRICT: "Strict assignments", BROADER: "Broader assignments"}
COLUMN_HEADERS = ("Detected fossil classes", "Precursor classes", "Source categories")


# --- Data ---------------------------------------------------------------------------


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def cite_text(record: dict[str, str]) -> str:
    surnames = [a.strip().split()[-1] for a in record["authors"].split(";")]
    if len(surnames) == 1:
        lead = surnames[0]
    elif len(surnames) == 2:
        lead = f"{surnames[0]} and {surnames[1]}"
    else:
        lead = f"{surnames[0]} et al."
    return f"{lead}, {record['year']}"


@dataclass
class Evidence:
    links: list[dict[str, str]]
    mapping_kind: dict[tuple[str, str], str]
    compat_kind: dict[tuple[str, str], str]
    citations: dict[str, list[str]]  # evidence kind -> formatted citations


def evidence_kind(row: dict[str, str], citations: dict, kind_citations: dict) -> str:
    basis = row["evidence_basis"]
    if row["assignment_criteria"] == STRICT:
        if basis in BROADER_KIND:
            raise ValueError(f"strict link {row} has a broader-only evidence basis")
        return STRICT_KIND
    if row["assignment_criteria"] != BROADER or basis not in BROADER_KIND:
        raise ValueError(f"unknown broader evidence basis: {basis}")
    kind = BROADER_KIND[basis]
    for cid in row["citation_ids"].split(";"):
        kind_citations.setdefault(kind, set()).add(cite_text(citations[cid]))
    return kind


def load_evidence(sources: Path) -> Evidence:
    citations = {
        r["citation_id"]: r for r in read_rows(sources / "TableS6_citations.csv")
    }
    kind_citations: dict[str, set[str]] = {}
    mapping_kind = {
        (r["analyte_id"], r["precursor_id"]): evidence_kind(
            r, citations, kind_citations
        )
        for r in read_rows(sources / "TableS6_analyte_precursor_map.csv")
        if r["assignment_criteria"] != "excluded"
    }
    compat_kind = {
        (r["precursor_id"], r["source_class_id"]): evidence_kind(
            r, citations, kind_citations
        )
        for r in read_rows(sources / "TableS6_source_compatibility.csv")
    }
    links = read_rows(sources / "TableS6_compatibility_links.csv")
    for r in links:
        kinds = {
            mapping_kind[(r["analyte_id"], r["precursor_id"])],
            compat_kind[(r["precursor_id"], r["source_class_id"])],
        }
        if (r["assignment_criteria"] == STRICT) != (kinds == {STRICT_KIND}):
            raise ValueError(f"link criteria differ from its evidence records: {r}")
    return Evidence(
        links,
        mapping_kind,
        compat_kind,
        {k: sorted(v) for k, v in kind_citations.items()},
    )


@dataclass
class Tier:
    analysis_id: str
    status: str
    compatible_count: int | None
    covers: list[tuple[str, ...]]  # minimal combinations, in the order of Table S7
    essential: tuple[str, ...]  # categories in every minimal combination
    missing: tuple[str, ...]  # detected classes with no assignment
    supporting: dict[tuple[str, ...], list[dict]]
    required: dict[tuple[str, ...], tuple[str, ...]] = field(default_factory=dict)


@dataclass
class Sample:
    name: str
    sample_id: str
    title: str
    absent_note: str
    detected: list[str]
    precursors: list[str]
    links: list[dict[str, str]]
    tiers: dict[str, Tier]


def load_tables(sources: Path) -> dict:
    return {
        "analyses": {
            r["analysis_id"]: r
            for r in read_jsonl(sources / "detected_class_sets.jsonl")
        },
        "scenarios": read_jsonl(sources / "sample_scenarios.jsonl"),
        "samples": {r["sample_id"]: r for r in read_jsonl(sources / "samples.jsonl")},
        "observations": read_rows(sources / "TableS5_source_observations.csv"),
        "covers": read_rows(sources / "TableS7_compatible_combinations.csv"),
        "states": {
            r["analysis_id"]: r
            for r in read_rows(sources / "TableS7_analysis_states.csv")
        },
        "supporting_rows": read_rows(
            sources / "TableS7_minimal_combination_assignments.csv"
        ),
        "examples": read_rows(sources / "figure4_examples.csv"),
    }


def selected_scenarios(sample_id: str, tables: dict) -> list[dict]:
    selected = [
        s
        for s in tables["scenarios"]
        if s["sample_id"] == sample_id
        and s["low_signal_mode"] == "include"
        and zero_threshold(s)
    ]
    if sorted(s["assignment_criteria"] for s in selected) != sorted((STRICT, BROADER)):
        raise ValueError(
            f"expected one zero-threshold scenario per set of assignments for {sample_id}"
        )
    return selected


def sample_title(sample_id: str, samples: dict) -> str:
    """Basin or formation, era and age of a sample, as written above its panel."""
    records = samples[sample_id]["source_records"]
    basin = next(
        (r["formation_or_basin_text"] for r in records if r["formation_or_basin_text"]),
        "",
    )
    era = next(
        (r["geological_era_text"] for r in records if r["geological_era_text"]), ""
    )
    ages = {r["age_central_ma"] for r in records}
    if len(ages) != 1:
        raise ValueError(f"conflicting ages for {sample_id}")
    place = ", ".join(part for part in (basin, era) if part)
    return f"{place}, about {next(iter(ages)):,.0f} Ma"


def plain(label: str, lower: bool = False) -> str:
    text = label.replace("\n", " ").replace("$_{", "").replace("}$", "")
    return text[0].lower() + text[1:] if lower else text


def absent_note(sample_id: str, detected: list[str], tables: dict) -> str:
    """The workbook value of the fossil classes not detected in a sample (zero or blank)."""
    values: dict[str, set[str]] = {}
    for row in tables["observations"]:
        if row["canonical_sample_key"] == sample_id and row["value_used"] == "true":
            values.setdefault(row["analyte_id"], set()).add(row["value_type"])
    absent = [a for a in FOSSIL if a not in detected]
    if absent == list(CONVENTIONAL_STERANES):
        name = "Conventional steranes"
    else:
        name = " and ".join(plain(FOSSIL[a]).split(" (")[0] for a in absent)
    kinds = {kind for a in absent for kind in values.get(a, {"not measured"})}
    if len(kinds) != 1 or next(iter(kinds)) not in VALUE_NOTE:
        raise ValueError(
            f"undetected classes of {sample_id} are neither all zero nor all blank: {kinds}"
        )
    return f"{name}: {VALUE_NOTE[next(iter(kinds))]}"


def load_sample(name: str, sample_id: str, tables: dict, evidence: Evidence) -> Sample:
    selected = selected_scenarios(sample_id, tables)
    detected = sorted(
        {a for s in selected for a in s["detected_analyte_ids"]}, key=list(FOSSIL).index
    )
    links = [r for r in evidence.links if r["analyte_id"] in detected]
    precursors = [p for p in PRECURSOR if any(r["precursor_id"] == p for r in links)]
    tiers = {}
    for scenario in selected:
        result = tables["analyses"][scenario["analysis_id"]]["result"]
        supporting = {
            tuple(w["source_class_ids"]): w["assignments"]
            for w in result["minimal_combination_assignments"]
        }
        tiers[scenario["assignment_criteria"]] = Tier(
            scenario["analysis_id"],
            result["status"],
            result["compatible_combinations"],
            [tuple(c) for c in result["minimal_combinations"]],
            tuple(result["categories_in_every_minimal"]),
            tuple(result["unassigned_analyte_ids"]),
            supporting,
        )
    return Sample(
        name,
        sample_id,
        sample_title(sample_id, tables["samples"]),
        absent_note(sample_id, detected, tables),
        detected,
        precursors,
        links,
        tiers,
    )


# --- Combinations and the evidence they need ----------------------------------------


def accepted(sample: Sample, tier: str) -> list[dict[str, str]]:
    return [
        r for r in sample.links if tier == BROADER or r["assignment_criteria"] == STRICT
    ]


def enumerate_covers(detected: list[str], links: list[dict]) -> tuple[list, list]:
    """Every compatible combination of the five categories, and the minimal ones."""
    reach = {
        a: {r["source_class_id"] for r in links if r["analyte_id"] == a}
        for a in detected
    }
    compatible = [
        frozenset(combo)
        for n in range(1, len(SOURCES) + 1)
        for combo in itertools.combinations(SOURCES, n)
        if all(reach[a] & set(combo) for a in detected)
    ]
    minimal = [c for c in compatible if not any(other < c for other in compatible)]
    return compatible, minimal


def link_kinds(row: dict[str, str], evidence: Evidence) -> set[str]:
    return {
        evidence.mapping_kind[(row["analyte_id"], row["precursor_id"])],
        evidence.compat_kind[(row["precursor_id"], row["source_class_id"])],
    } - {STRICT_KIND}


def required_kinds(
    sample: Sample, tier: str, cover: tuple, evidence: Evidence
) -> tuple:
    """Broader evidence kinds without which some detected class loses every link into the combination."""
    needed = set()
    for a in sample.detected:
        routes = [
            r
            for r in accepted(sample, tier)
            if r["analyte_id"] == a and r["source_class_id"] in cover
        ]
        if not routes:
            raise ValueError(
                f"{cover} does not account for {a} in {sample.name} ({tier})"
            )
        kinds = [link_kinds(r, evidence) for r in routes]
        needed |= {k for k in set().union(*kinds) if all(k in ks for ks in kinds)}
    return tuple(k for k in KIND_ORDER if k in needed)


def check_tier(sample: Sample, tier: str, tables: dict, evidence: Evidence) -> None:
    """Check the drawn links against Table S7 and record the broader evidence each combination needs."""
    t = sample.tiers[tier]
    state = tables["states"][t.analysis_id]
    rows = [r for r in tables["covers"] if r["analysis_id"] == t.analysis_id]
    compatible, minimal = enumerate_covers(sample.detected, accepted(sample, tier))
    if t.status != "evaluable":
        orphans = tuple(
            a
            for a in sample.detected
            if not any(r["analyte_id"] == a for r in accepted(sample, tier))
        )
        if (
            orphans != t.missing
            or compatible
            or rows
            or t.covers
            or state["status"] != "not_evaluable"
            or state["reason_code"] != "class_without_assignment"
        ):
            raise ValueError(f"not-evaluable result differs for {sample.name} ({tier})")
        return
    minimal_covers = {frozenset(c) for c in t.covers}
    table_minimal = {
        frozenset(r["source_class_ids"].split(";"))
        for r in rows
        if r["minimal"] == "true"
    }
    if not (minimal_covers == table_minimal == set(minimal)):
        raise ValueError(f"minimal combinations differ for {sample.name} ({tier})")
    if not (
        t.compatible_count
        == len(rows)
        == len(compatible)
        == int(state["compatible_combinations"])
    ):
        raise ValueError(f"compatible counts differ for {sample.name} ({tier})")
    if set(t.essential) != frozenset.intersection(*minimal_covers):
        raise ValueError(
            f"categories in every minimal combination differ for {sample.name} ({tier})"
        )
    link_keys = {
        (
            r["analyte_id"],
            r["mapping_record_id"],
            r["compatibility_record_id"],
            r["source_class_id"],
        )
        for r in accepted(sample, tier)
    }
    table_supporting = {
        tuple(r["source_class_ids"].split(";"))
        for r in tables["supporting_rows"]
        if r["analysis_id"] == t.analysis_id
    }
    for cover in t.covers:
        assignments = t.supporting[cover]
        if cover not in table_supporting or {
            a["analyte_id"] for a in assignments
        } != set(sample.detected):
            raise ValueError(
                f"supporting assignments missing or incomplete for {cover}"
            )
        for a in assignments:
            key = (
                a["analyte_id"],
                a["mapping_record_id"],
                a["compatibility_record_id"],
                a["source_class_id"],
            )
            if key not in link_keys:
                raise ValueError(f"supporting assignment is not a drawn link: {a}")
        t.required[cover] = required_kinds(sample, tier, cover, evidence)


def cielab(color: str, matrix) -> np.ndarray:
    srgb = np.array(to_rgb(color))
    linear = np.where(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055) ** 2.4)
    linear = np.clip(np.array(matrix) @ linear, 0, 1)
    to_xyz = np.array(
        [[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]]
    )
    xyz = to_xyz @ linear / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def check_palette() -> dict[str, float]:
    """Smallest CIELAB distance between any two line colors, with normal and simulated vision."""
    worst = {}
    for name, matrix in {"normal": np.eye(3).tolist(), **SIMULATION}.items():
        labs = [cielab(style["color"], matrix) for style in LINE_STYLE.values()]
        worst[name] = round(
            float(
                min(np.linalg.norm(a - b) for a, b in itertools.combinations(labs, 2))
            ),
            1,
        )
    if min(worst.values()) < MIN_LINE_COLOR_DISTANCE:
        raise ValueError(f"line colors too close under simulated color vision: {worst}")
    return worst


def numbering(sample: Sample) -> dict[tuple, int]:
    """One number per distinct combination, shared by both sets of assignments (strict first)."""
    numbers: dict[tuple, int] = {}
    for tier in (STRICT, BROADER):
        for cover in sample.tiers[tier].covers:
            numbers.setdefault(cover, len(numbers) + 1)
    return numbers


# --- Text in millimeters -------------------------------------------------------------


def measure(text: str, size: float, weight: str = "normal") -> float:
    width, _, _ = TextToPath().get_text_width_height_descent(
        text,
        FontProperties(family="Arial", size=size, weight=weight),
        ismath="$" in text,
    )
    return width * PT_TO_MM


def wrap_label(label: str, max_mm: float, size: float) -> str:
    """Wrap a box label into the fewest lines that fit, balancing their measured widths."""
    words = label.replace("\n", " ").split()
    for n_lines in range(1, len(words) + 1):
        best = None
        for cuts in itertools.combinations(range(1, len(words)), n_lines - 1):
            bounds = (0, *cuts, len(words))
            lines = [" ".join(words[i:j]) for i, j in zip(bounds, bounds[1:])]
            widest = max(measure(line, size) for line in lines)
            if widest <= max_mm and (best is None or widest < best[0]):
                best = (widest, lines)
        if best:
            return "\n".join(best[1])
    return "\n".join(words)


def label_text(label: str, lower: bool = False) -> str:
    text = label.replace("\n", " ")
    return text[0].lower() + text[1:] if lower else text


@dataclass(frozen=True)
class Run:
    text: str = ""
    color: str = "black"
    weight: str = "normal"
    swatch: str | None = (
        None  # evidence kind drawn as a short line sample before the text
    )
    lead: str = (
        ""  # a word kept on the same line before the line sample ("needs", "and")
    )


def draw_runs(
    ax,
    x: float,
    y: float,
    width: float,
    runs: list[Run],
    size: float,
    pitch: float,
    indent: float = 0.0,
) -> float:
    """Draw mixed-style text with word wrapping; returns the y below the last line."""
    swatch_mm, space = 5.0, measure("a a", size) - measure("aa", size)
    cursor, line_x = x, x
    for run in runs:
        # A line sample and its label stay together on one line; plain text wraps by word.
        pieces = (
            [run]
            if run.swatch
            else [Run(w, run.color, run.weight) for w in run.text.split(" ") if w]
        )
        for piece in pieces:
            lead_w = (
                measure(piece.lead, size, piece.weight) + space if piece.lead else 0.0
            )
            text_w = measure(piece.text, size, piece.weight) if piece.text else 0.0
            w = lead_w + (
                swatch_mm + (space + text_w if piece.text else 0.0)
                if piece.swatch
                else text_w
            )
            if cursor > line_x and cursor + w > x + width:
                y += pitch
                line_x = cursor = x + indent
            text_x = cursor
            if piece.lead:
                ax.text(
                    cursor,
                    y,
                    piece.lead,
                    ha="left",
                    va="top",
                    fontsize=size,
                    color=piece.color,
                    weight=piece.weight,
                )
                text_x = cursor + lead_w
            if piece.swatch:
                ax.plot(
                    [text_x, text_x + swatch_mm],
                    [y + size * PT_TO_MM * 0.55] * 2,
                    solid_capstyle="butt",
                    **LINE_STYLE[piece.swatch],
                )
                text_x += swatch_mm + space
            if piece.text:
                ax.text(
                    text_x,
                    y,
                    piece.text,
                    ha="left",
                    va="top",
                    fontsize=size,
                    color=piece.color,
                    weight=piece.weight,
                )
            cursor += w + space
    return y + pitch


def combination_runs(sample: Sample, tier: str, numbers: dict) -> list[list[Run]]:
    """The written answer for one set of assignments: a header line, then one line per minimal combination."""
    t = sample.tiers[tier]
    if t.status != "evaluable":
        missing = " and ".join(label_text(FOSSIL[a], lower=True) for a in t.missing)
        return [
            [
                Run(
                    f"Not evaluable: no strict assignment links {missing} to a precursor class, "
                    "so no combination of source categories accounts for every detected class.",
                    TEXT_GRAY,
                )
            ]
        ]
    n, m = len(t.covers), t.compatible_count
    common = t.essential if n > 1 else ()
    if n == 1:
        header = f"1 minimal combination ({m} compatible):"
    elif common:
        header = (
            f"{n} minimal combinations ({m} compatible). Each includes "
            f"{' and '.join(label_text(SOURCE[s], lower=True) for s in common)}, plus:"
        )
    else:
        header = f"{n} minimal combinations ({m} compatible):"
    lines = [[Run(header)]]
    for cover in sorted(t.covers, key=numbers.get):
        members = sorted((s for s in cover if s not in common), key=SOURCES.index)
        runs = [
            Run(str(numbers[cover]), weight="bold"),
            Run(" + ".join(label_text(SOURCE[s], lower=True) for s in members)),
        ]
        if tier == BROADER and cover in sample.tiers[STRICT].covers:
            runs.append(Run("(as under strict)", TEXT_GRAY))
        for i, kind in enumerate(t.required[cover]):
            runs.append(
                Run(kind.lower(), TEXT_GRAY, swatch=kind, lead="and" if i else "needs")
            )
        lines.append(runs)
    return lines


def draw_combinations(ax, x, y, width, sample, tier, numbers, size, pitch) -> float:
    for i, runs in enumerate(combination_runs(sample, tier, numbers)):
        if i == 0 or runs[0].weight != "bold":
            y = draw_runs(ax, x, y, width, runs, size, pitch)
            continue
        number, rest = runs[0], runs[1:]
        ax.text(x, y, number.text, ha="left", va="top", fontsize=size, weight="bold")
        y = draw_runs(ax, x + 3.2, y, width - 3.2, rest, size, pitch)
    return y


# --- Drawing ------------------------------------------------------------------------


def axes_mm(fig, x: float, y: float, w: float, h: float):
    """An axes placed in millimeters from the top-left corner, with data units of 1 mm."""
    fig_w, fig_h = fig.get_size_inches() * 25.4
    ax = fig.add_axes((x / fig_w, 1 - (y + h) / fig_h, w / fig_w, h / fig_h))
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return ax


def box(ax, cx, cy, w, h, label, fill, size, active=True) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (cx - w / 2, cy - h / 2),
            w,
            h,
            boxstyle="round,pad=0,rounding_size=1.2",
            facecolor=fill if active else "white",
            edgecolor="black" if active else LIGHT,
            lw=0.9 if active else 0.6,
        )
    )
    ax.text(
        cx,
        cy,
        wrap_label(label, w - 2.0, size),
        ha="center",
        va="center",
        fontsize=size,
        color="black" if active else GRAY,
        linespacing=1.1,
    )


def segment(ax, x0, y0, x1, y1, kind) -> None:
    ax.plot([x0, x1], [y0, y1], zorder=0, solid_capstyle="butt", **LINE_STYLE[kind])


def anchors(pairs: list[tuple], start_y: dict, end_y: dict, span: float = 0.42) -> dict:
    """Spread the ends of lines that meet the same box side, ordered by where they come from."""
    ends = {}
    for side, own, other in (("start", start_y, end_y), ("end", end_y, start_y)):
        grouped = {}
        for pair in pairs:
            key = pair[0] if side == "start" else pair[1]
            grouped.setdefault(key, []).append(pair)
        for key, group in grouped.items():
            group.sort(
                key=lambda pr: -(other[pr[1]] if side == "start" else other[pr[0]])
            )
            n = len(group)
            for i, pair in enumerate(group):
                offset = 0 if n == 1 else span / 2 - i * span / (n - 1)
                ends[(side, pair)] = own[key] + offset
    return ends


@dataclass(frozen=True)
class Geometry:
    columns: tuple[
        tuple[float, float], ...
    ]  # (center x, width) for fossil, precursor, source
    top: float  # y of the first source row
    pitch: float
    box_h: float
    header_y: float
    size: float


def diagram(
    ax, sample: Sample, tier: str, g: Geometry, evidence: Evidence, headers: bool = True
):
    ys = {s: g.top + i * g.pitch for i, s in enumerate(SOURCES)}
    span = (len(SOURCES) - 1) * g.pitch

    def spread(n: int) -> list[float]:
        return (
            [g.top + span / 2]
            if n == 1
            else [g.top + i * span / (n - 1) for i in range(n)]
        )

    yf = dict(zip(sample.detected, spread(len(sample.detected))))
    yp = dict(zip(sample.precursors, spread(len(sample.precursors))))
    links = accepted(sample, tier)
    first = sorted({(r["analyte_id"], r["precursor_id"]) for r in links})
    second = sorted({(r["precursor_id"], r["source_class_id"]) for r in links})
    (xf, wf), (xp, wp), (xs, ws) = g.columns
    ends = anchors(first, yf, yp, span=g.box_h * 0.55)
    for a, p in first:
        segment(
            ax,
            xf + wf / 2,
            ends[("start", (a, p))],
            xp - wp / 2,
            ends[("end", (a, p))],
            evidence.mapping_kind[(a, p)],
        )
    ends = anchors(second, yp, ys, span=g.box_h * 0.55)
    for p, s in second:
        segment(
            ax,
            xp + wp / 2,
            ends[("start", (p, s))],
            xs - ws / 2,
            ends[("end", (p, s))],
            evidence.compat_kind[(p, s)],
        )
    linked_p, linked_s = {p for _, p in first}, {s for _, s in second}
    for a in sample.detected:
        box(ax, xf, yf[a], wf, g.box_h, FOSSIL[a], FILL["fossil"], g.size)
    for p in sample.precursors:
        box(
            ax,
            xp,
            yp[p],
            wp,
            g.box_h,
            PRECURSOR[p],
            FILL["precursor"],
            g.size,
            p in linked_p,
        )
    for s in SOURCES:
        box(
            ax, xs, ys[s], ws, g.box_h, SOURCE[s], FILL["source"], g.size, s in linked_s
        )
    if headers:
        for (x, _), label in zip(g.columns, COLUMN_HEADERS):
            ax.text(
                x,
                g.header_y,
                label,
                ha="center",
                va="center",
                fontsize=g.size,
                weight="bold",
            )
    return ys, g.top + span + g.box_h / 2


def legend(ax, evidence: Evidence, size: float, pitch: float, width: float) -> None:
    gap = Run(" " * 3)  # non-breaking spaces keep a visible gap between legend items
    runs = [Run("Lines are accepted evidence links:", weight="bold"), gap]
    for kind in KIND_ORDER:
        cites = evidence.citations.get(kind)
        runs += [
            Run(f"{kind} ({'; '.join(cites)})" if cites else kind, swatch=kind),
            gap,
        ]
    y = draw_runs(ax, 1, 0.5, width, runs, size, pitch)
    draw_runs(
        ax,
        1,
        y,
        width,
        [
            Run(
                "Line samples beside a combination mark the broader evidence it needs.",
                TEXT_GRAY,
            )
        ],
        size,
        pitch,
    )


def draw(samples: list[Sample], evidence: Evidence, out: Path) -> list[Path]:
    cell_w, gap, size, pitch = 88, 4, 6, 2.5
    # Each panel is as tall as its diagram and list of combinations; panel A lists three.
    cell_h = {"A": 71.5, "B": 62.5}
    legend_h = 9
    fig_w = 183  # double-column width, mm
    fig_h = 6.5 + sum(4.5 + h + 2 for h in cell_h.values()) + legend_h
    # The precursor column is wide enough for "Cycloartane-skeleton" on one line.
    g = Geometry(
        ((10.5, 19), (36.5, 23), (71, 30)),
        top=10,
        pitch=9.2,
        box_h=7.6,
        header_y=4.0,
        size=size,
    )
    fig = plt.figure(figsize=(fig_w * MM, fig_h * MM))
    x_cells = (3, 3 + cell_w + gap)
    header = axes_mm(fig, 0, 0, fig_w, 6)
    for x0, tier in zip(x_cells, (STRICT, BROADER)):
        header.text(
            x0 + cell_w / 2,
            3.2,
            TIER_LABEL[tier],
            ha="center",
            va="center",
            fontsize=8,
            weight="bold",
        )
    y = 6.5
    for letter, sample in zip("AB", samples):
        title = axes_mm(fig, 0, y, fig_w, 4.5)
        title.text(1.5, 2.25, letter, fontsize=10, weight="bold", va="center")
        heading = f"{sample.name}  ·  {sample.title}"
        title.text(8, 2.25, heading, fontsize=8, weight="bold", va="center")
        title.text(
            8 + measure(heading, 8, "bold") + 4,
            2.25,
            sample.absent_note,
            fontsize=6.5,
            color=TEXT_GRAY,
            va="center",
        )
        y += 4.5
        numbers = numbering(sample)
        for x0, tier in zip(x_cells, (STRICT, BROADER)):
            ax = axes_mm(fig, x0, y, cell_w, cell_h[letter])
            _, bottom = diagram(ax, sample, tier, g, evidence)
            draw_combinations(
                ax, 0.5, bottom + 3.0, cell_w - 1.5, sample, tier, numbers, size, pitch
            )
        y += cell_h[letter] + 2
    legend(axes_mm(fig, 0, y, fig_w, fig_h - y), evidence, size, pitch, fig_w - 4)
    out.mkdir(parents=True, exist_ok=True)
    pdf, png = out / "Figure4.pdf", out / "Figure4.png"
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return [pdf, png]


# --- Main ---------------------------------------------------------------------------


def plotted_values(samples: list[Sample], evidence: Evidence) -> dict:
    return {
        "evidence_kind_citations": evidence.citations,
        "samples": {
            s.name: {
                "sample_id": s.sample_id,
                "title": s.title,
                "absent_note": s.absent_note,
                "detected": s.detected,
                "links_drawn": [
                    {
                        k: r[k]
                        for k in (
                            "analyte_id",
                            "precursor_id",
                            "source_class_id",
                            "assignment_criteria",
                        )
                    }
                    for r in s.links
                ],
                "assignments": {
                    tier: {
                        "analysis_id": t.analysis_id,
                        "status": t.status,
                        "compatible_combinations": t.compatible_count,
                        "minimal_combinations": [list(c) for c in t.covers],
                        "in_every_minimal_combination": list(t.essential),
                        "detected_without_assignment": list(t.missing),
                        "broader_evidence_required": {
                            " + ".join(c): list(t.required[c]) for c in t.covers
                        },
                    }
                    for tier, t in s.tiers.items()
                },
            }
            for s in samples
        },
        "line_colors": {kind: style["color"] for kind, style in LINE_STYLE.items()},
        "box_fills": FILL,
        "smallest_cielab_distance_between_line_colors": check_palette(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--results", type=Path, required=True, help="Folder holding sources/."
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder.")
    args = parser.parse_args()
    sources = args.results / "sources"
    check_palette()
    tables, evidence = load_tables(sources), load_evidence(sources)
    samples = [
        load_sample(row["sample_label"], row["sample_id"], tables, evidence)
        for row in tables["examples"]
    ]
    for sample in samples:
        for tier in (STRICT, BROADER):
            check_tier(sample, tier, tables, evidence)
    paths = draw(samples, evidence, args.out)
    write_json(args.out / "Figure4_values.json", plotted_values(samples, evidence))
    for sample in samples:
        for tier, t in sample.tiers.items():
            print(
                sample.name,
                tier,
                t.status,
                t.compatible_count,
                len(t.covers),
                "minimal combinations",
            )
    print(*[p.name for p in paths], sep="\n")


if __name__ == "__main__":
    main()
