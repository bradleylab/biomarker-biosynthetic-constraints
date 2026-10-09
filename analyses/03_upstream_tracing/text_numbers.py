"""Print the numbers on upstream tracing, enumeration and literature records quoted in the paper.

python analyses/03_upstream_tracing/text_numbers.py --results results
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.enumeration import MAX_SEARCH_STATES  # noqa: E402
from biosynthesis.tables import read_csv  # noqa: E402


def select(rows: list[dict], **values) -> list[dict]:
    """Rows whose columns equal the given values."""
    return [
        row for row in rows if all(row[key] == value for key, value in values.items())
    ]


def shared_network(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Split rows at the largest gap in the number of compounds traced: (shared network, others)."""
    sizes = sorted({int(row["compounds_traced"]) for row in rows})
    gaps = [(b - a, b) for a, b in zip(sizes, sizes[1:])]
    threshold = max(gaps)[1] if gaps else sizes[0]
    large = [row for row in rows if int(row["compounds_traced"]) >= threshold]
    return large, [row for row in rows if int(row["compounds_traced"]) < threshold]


def span(values) -> str:
    """A value, or the range of values, with thousands separators."""
    values = list(values)
    low, high = min(values), max(values)
    return f"{low:,}" if low == high else f"{low:,}–{high:,}"


def accepted_with_producer(
    tracing: list[dict], version: str, criteria: str
) -> list[dict]:
    return [
        row
        for row in select(
            tracing,
            network_version=version,
            reaction_criteria=criteria,
            identity="accepted",
        )
        if int(row["producing_reactions"]) > 0
    ]


def biomarkers(rows: list[dict]) -> set[str]:
    return {row["biomarker_number"] for row in rows}


def enzymes_named(connections: str) -> set[str]:
    """The enzyme names of the connections of a Table S8 record (a name or a list of names each)."""
    names = set()
    for item in json.loads(connections):
        value = item.get("enzyme")
        names.update(value if isinstance(value, list) else [value])
    return names


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    folder = args.results / "tracing"
    tracing = read_csv(folder / "TableS1_tracing.csv")
    producing = read_csv(folder / "TableS1_producing_reactions.csv")
    representative = read_csv(folder / "TableS6_representative_targets.csv")
    native = read_csv(folder / "TableS1_tracing_native_control.csv")
    native_biomarkers = read_csv(
        folder / "TableS1_producing_reactions_native_control.csv"
    )
    literature = read_csv(
        args.results / "literature" / "TableS8_literature_reactions.csv"
    )

    print("Biomarker identities")
    first = select(tracing, network_version="primary", reaction_criteria="strict")
    identity_of = {}
    for row in first:
        identity_of.setdefault(row["biomarker_number"], set()).add(row["identity"])
    counts = Counter("/".join(sorted(states)) for states in identity_of.values())
    print(
        f"  {len(identity_of)} biomarkers: {counts['accepted']} assigned a single MetaCyc "
        f"compound, {len(identity_of) - counts['accepted']} defined from the literature, of which "
        f"{counts['no_candidate']} have no candidate identity."
    )
    variants = Counter(
        row["identity"] for row in {r["identity_id"]: r for r in first}.values()
    )
    print(
        f"  {sum(variants.values())} identity variants: {variants['accepted']} accepted, "
        f"{variants['candidate_not_accepted']} candidate identities not accepted, "
        f"{variants['no_candidate']} without a candidate."
    )
    primary = select(tracing, network_version="primary")
    records = {(row["reaction_criteria"], row["identity_id"]) for row in primary}
    doubled = Counter((row["reaction_criteria"], row["identity_id"]) for row in primary)
    print(
        f"  Under the two sets of criteria: {len(records)} identity records and {len(primary)} rows "
        f"per network version ({len(tracing)} rows in Table S1 tracing); "
        f"{sum(v == 2 for v in doubled.values())} identity records have two chemical-state rows."
    )

    print("Producing reactions (accepted identities)")
    for version in ("primary", "with_conditional"):
        for criteria in ("strict", "broader"):
            rows = select(
                producing,
                network_version=version,
                reaction_criteria=criteria,
                identity="accepted",
            )
            print(
                f"  {version}, {criteria}: {sum(r['accepted_identity_has_producer'] == 'true' for r in rows)}"
                f" of {len(rows)} accepted biomarkers have at least one producing reaction."
            )
    gained = [
        r["name"]
        for r in select(
            producing,
            network_version="with_conditional",
            reaction_criteria="broader",
            identity="accepted",
        )
        if r["accepted_identity_has_producer"] == "true"
        and select(
            producing,
            network_version="primary",
            reaction_criteria="broader",
            biomarker_number=r["biomarker_number"],
        )[0]["accepted_identity_has_producer"]
        == "false"
    ]
    print(
        f"  Gained a producer only with conditional connections (broader): {len(gained)} "
        f"({', '.join(gained)})."
    )

    print(
        "Shared networks (primary graph, accepted biomarkers with a producing reaction)"
    )
    for criteria in ("broader", "strict"):
        rows = accepted_with_producer(tracing, "primary", criteria)
        large, small = shared_network(rows)
        print(
            f"  {criteria}: {len(biomarkers(large))} of {len(biomarkers(rows))} biomarkers trace into "
            f"a network of {span(int(r['compounds_traced']) for r in large)} compounds with "
            f"{span(int(r['compounds_without_producer']) for r in large)} compounds without an "
            f"accepted producer; the other {len(biomarkers(small))}: "
            f"{', '.join(sorted({r['name'] for r in small}))} "
            f"({span(int(r['compounds_traced']) for r in small)} compounds traced)."
        )
        exact = [r for r in rows if r["status"] == "exact"]
        print(
            f"  {criteria}: the search finished (exact) for {len(biomarkers(exact))} of "
            f"{len(biomarkers(rows))} biomarkers with a producing reaction."
        )

    print("Search status of all rows")
    print(
        f"  Search limit: {MAX_SEARCH_STATES} search states per query; largest number examined: "
        f"{max(int(r['search_states_examined']) for r in tracing if r['search_states_examined'])}."
    )
    for version, criteria in (
        ("primary", "strict"),
        ("primary", "broader"),
        ("with_conditional", "strict"),
        ("with_conditional", "broader"),
    ):
        rows = select(tracing, network_version=version, reaction_criteria=criteria)
        status = Counter(row["status"] for row in rows)
        print(
            f"  {version}, {criteria}: {status['exact']} exact, {status['open']} open, "
            f"{status['not_in_network']} not in network."
        )
    for criteria in ("strict", "broader"):
        exact = select(
            tracing,
            network_version="primary",
            reaction_criteria=criteria,
            status="exact",
        )
        producer = [r for r in exact if int(r["producing_reactions"]) > 0]
        others = sorted({r["name"] for r in producer if r["identity"] != "accepted"})
        print(
            f"  primary, {criteria}: of {len(exact)} exact rows, "
            f"{len(exact) - len(producer)} have no direct producing reaction and {len(producer)} "
            f"at least one ({sum(r['identity'] == 'accepted' for r in producer)} for accepted "
            f"biomarkers"
            + (f"; the others: {', '.join(others)}" if others else "")
            + ")."
        )
    distinct = {
        (r["network_version"], r["reaction_criteria"], r["node"])
        for r in tracing
        if r["node"]
    }
    for version in ("primary", "with_conditional"):
        print(
            f"  {version}: {sum(v == version for v, _, _ in distinct)} distinct traced queries."
        )

    print("Representative targets (Table S6)")
    for version in ("primary", "with_conditional"):
        for criteria in ("strict", "broader"):
            rows = select(
                representative, network_version=version, reaction_criteria=criteria
            )
            described = "; ".join(
                f"{r['name']}: {r['status']}, {int(r['compounds_traced']):,} compounds, "
                f"{int(r['compounds_without_producer']):,} without an accepted producer, "
                f"{r['producing_reactions']} producing reactions"
                for r in rows
            )
            print(f"  {version}, {criteria}: {described}.")

    print("Formylhopane reductase (RXN-13531) removed")
    for row in native:
        if row["in_network"] == "true" and (
            row["producing_reactions_with_rxn13531"]
            != row["producing_reactions_without_rxn13531"]
        ):
            print(
                f"  {row['network_version']}, {row['reaction_criteria']}: {row['name']}: "
                f"{row['producing_reactions_with_rxn13531']} producing reactions with it, "
                f"{row['producing_reactions_without_rxn13531']} without."
            )
    changed = sum(
        r["has_producer_with_rxn13531"] != r["has_producer_without_rxn13531"]
        for r in native_biomarkers
    )
    print(f"  Biomarker rows whose producer status changes: {changed}.")

    print("Literature reaction records (Table S8)")
    disposition = Counter(row["disposition"] for row in literature)
    print(
        f"  {len(literature)} records: {disposition['supported']} with only supported connections, "
        f"{disposition['supported_and_conditional']} with both, {disposition['conditional']} with "
        f"only conditional connections; {disposition['supported'] + disposition['supported_and_conditional']}"
        f" contain at least one supported connection."
    )
    methyl = [
        row["record_id"]
        for row in literature
        if {"Cgm", "Pgm"} <= enzymes_named(row["substrate_product_connections_json"])
    ]
    print(f"  Records documenting Cgm and Pgm: {len(methyl)} ({', '.join(methyl)}).")


if __name__ == "__main__":
    main()
