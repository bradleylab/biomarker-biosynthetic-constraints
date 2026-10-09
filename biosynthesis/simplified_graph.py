"""The simplified compound graph: one selected precursor per product-specific reaction record.

The graph is built from the reactions of the MetaCyc base pathways outside the
excluded pathway classes, split into one record per product and EC number, with
the literature connections and MetaCyc corrections of the reaction network added.
For each record, the reactant with the greatest recorded molecular weight is kept
as the precursor after hubs and carriers are excluded.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence

EXCLUDED_PATHWAY_CLASSES = frozenset(
    {
        "Energy-Metabolism",
        "Super-Pathways",
        "Amino-Acid-Biosynthesis",
        "Nucleotide-Biosynthesis",
        "Macromolecule-Modification",
    }
)
# Compounds that are reactants in more than this many product-specific records are hubs.
HUB_LIMIT = 20
# Exempt from the hub exclusion: precursors of the aromatic carotenoids and of the
# sterols and cyclic triterpenoids.
HUB_EXEMPTIONS = frozenset({"CPD1F-129", "EPOXYSQUALENE"})
TRANSFER_AND_REDOX_CARRIERS = frozenset(
    {
        "PLASTOQUINONE",
        "PLASTOQUINOL",
        "PLASTOQUINONE-9",
        "PLASTOQUINOL-9",
        "Menaquinones",
        "Ubiquinones",
        "Demethylmenaquinones",
        "MENAQUINONE",
        "MENAQUINOL",
        "UBIQUINONE-8",
        "UBIQUINOL-8",
        "Reduced-2Fe-2S-Ferredoxins",
        "Oxidized-2Fe-2S-Ferredoxins",
        "FE+3",
        "FE+2",
    }
)
FLAT_FILES = ("classes.dat", "pathways.dat", "reactions.dat", "compounds.dat")
# The simplified graph takes its literature connections from the network built
# under the broader reaction criteria.
NETWORK_FOLDERS = {
    "primary": "broader_primary",
    "with_conditional": "broader_with_conditional",
}
NETWORK_VERSIONS = tuple(NETWORK_FOLDERS)
EXCLUDING_ACTIONS = frozenset({"replace", "exclude"})


def read_flat_file(path: Path) -> dict[str, dict[str, list[str]]]:
    """Read a MetaCyc attribute-value file, keeping repeated slots in file order.

    Comment and annotation lines (starting with '/', '#' or '^') are skipped.
    """
    result = {}
    for block in Path(path).read_text(encoding="latin-1").split("\n//\n"):
        fields: dict[str, list[str]] = {}
        for line in block.splitlines():
            if " - " not in line or line.startswith(("/", "#", "^")):
                continue
            key, value = line.split(" - ", 1)
            fields.setdefault(key, []).append(value)
        identifiers = fields.get("UNIQUE-ID", [])
        if identifiers:
            if len(identifiers) != 1 or identifiers[0] in result:
                raise ValueError(f"Ambiguous record identifier: {identifiers}")
            result[identifiers[0]] = fields
    return result


def orientations(record: Mapping[str, list[str]]) -> list[tuple[list[str], list[str]]]:
    """Return the (reactants, products) orientations a reaction contributes.

    Reactions with one left-to-right or right-to-left direction give one
    orientation, reversible reactions give both, and all other reactions none.
    """
    directions = record.get("REACTION-DIRECTION", [])
    if len(directions) != 1:
        return []
    left, right = record.get("LEFT", []), record.get("RIGHT", [])
    if directions[0] in {"LEFT-TO-RIGHT", "PHYSIOL-LEFT-TO-RIGHT"}:
        return [(left, right)]
    if directions[0] in {"RIGHT-TO-LEFT", "PHYSIOL-RIGHT-TO-LEFT"}:
        return [(right, left)]
    if directions[0] == "REVERSIBLE":
        return [(left, right), (right, left)]
    return []


def precursor_links(network: Mapping[str, Sequence[Mapping]]) -> dict[str, list[str]]:
    """Map each product to its sorted selected precursors, without self-links."""
    result = {}
    for product, rows in network.items():
        inputs = sorted(
            {
                row["reactant"][0]
                for row in rows
                if row.get("reactant") and row["reactant"][0] != product
            }
        )
        if inputs:
            result[product] = inputs
    return result


def pathway_domain(classes: Mapping, pathways: Mapping) -> list[str]:
    """Base pathways (no sub-pathways) none of whose classes descends from an excluded class."""

    @lru_cache(None)
    def ancestors(identifier: str) -> frozenset[str]:
        return frozenset(
            {identifier}
            | {
                value
                for parent in classes.get(identifier, {}).get("TYPES", [])
                for value in ancestors(parent)
            }
        )

    return sorted(
        identifier
        for identifier, record in pathways.items()
        if not {
            value for parent in record.get("TYPES", []) for value in ancestors(parent)
        }
        & EXCLUDED_PATHWAY_CLASSES
        and not record.get("SUB-PATHWAYS")
    )


def product_records(
    reactions: Mapping, pathways: Mapping, domain: Sequence[str]
) -> dict:
    """One record per product, reactant list and EC number from the pathway reactions.

    Records that would make a product from itself alone are skipped, and repeated
    records of a product are kept once.
    """
    products: dict[str, list[dict]] = {}
    for pathway in domain:
        for identifier in pathways[pathway].get("REACTION-LIST", []):
            if identifier not in reactions:
                raise ValueError(
                    f"Pathway has non-reaction component: {pathway}/{identifier}"
                )
            record = reactions[identifier]
            codes = [[value.strip("|")] for value in record.get("EC-NUMBER", [])] or [
                []
            ]
            for inputs, outputs in orientations(record):
                inputs = [value.strip().strip("|") for value in inputs]
                for product in outputs:
                    product = product.strip().strip("|")
                    rows = products.setdefault(product, [])
                    if set(inputs) == {product}:
                        continue
                    for code in codes:
                        row = {"reactant": inputs, "ec_number": code}
                        if row not in rows:
                            rows.append(row)
    return products


def molecular_weights(compounds: Mapping) -> dict[str, float | None]:
    """The recorded molecular weight of each compound when it is a single finite value."""
    masses = {}
    for identifier, record in compounds.items():
        literals = record.get("MOLECULAR-WEIGHT", [])
        value = None
        if len(literals) == 1:
            try:
                numeric = float(literals[0])
                value = numeric if math.isfinite(numeric) else None
            except ValueError:
                pass
        masses[identifier] = value
    return masses


def signature(product: str, row: Mapping) -> tuple:
    """Product, sorted reactants and EC number of a record."""
    return (
        product,
        tuple(sorted(value.strip("|") for value in row["reactant"])),
        tuple(value.strip("|") for value in row.get("ec_number", [])),
    )


def supporting_reactions(reactions: Mapping) -> dict[tuple, set[str]]:
    """The MetaCyc reactions that give each record signature, over all reactions."""
    result = defaultdict(set)
    for identifier, record in reactions.items():
        for inputs, outputs in orientations(record):
            codes = [[value] for value in record.get("EC-NUMBER", [])] or [[]]
            for output in outputs:
                for code in codes:
                    if set(inputs) != {output}:
                        result[
                            signature(output, {"reactant": inputs, "ec_number": code})
                        ].add(identifier)
    return result


def add_literature_connections(
    native: Mapping,
    excluded: set[str],
    support: Mapping,
    changes: Mapping,
    nodes: Mapping,
) -> dict:
    """Apply the MetaCyc corrections and add the literature connections of a network.

    A MetaCyc record is removed when every reaction that gives it is excluded or
    replaced by a correction. Each literature connection adds a record per product
    unless the product already has a record with the same set of reactants.
    """
    products: dict[str, list[dict]] = {}
    for product, rows in native.items():
        products[product] = []
        for row in rows:
            identifiers = support.get(signature(product, row), set())
            if not (identifiers & excluded and not identifiers - excluded):
                products[product].append(deepcopy(row))
    for edge in changes["added_edges"]:
        inputs = [nodes[value]["identity"]["local_id"] for value in edge["inputs"]]
        outputs = [nodes[value]["identity"]["local_id"] for value in edge["outputs"]]
        if any(value is None for value in inputs + outputs):
            raise ValueError("A literature compound has no identifier")
        for product in outputs:
            rows = products.setdefault(product, [])
            match = next(
                (row for row in rows if set(row["reactant"]) == set(inputs)), None
            )
            evidence = {
                "reaction_id": edge["reaction_local_id"],
                "source_occurrence_id": edge["source_occurrence_id"],
                "edge_id": edge["edge_id"],
            }
            if match is not None:
                match.setdefault("_literature_evidence", []).append(evidence)
            else:
                rows.append(
                    {
                        "reactant": inputs,
                        "ec_number": [],
                        "_source": "literature_connection",
                        "_literature_evidence": [evidence],
                    }
                )
    return products


def select_precursors(
    product_dict: Mapping[str, Sequence[Mapping[str, Any]]],
    masses: Mapping[str, Any],
    hub_exemptions: frozenset[str] = HUB_EXEMPTIONS,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Keep one precursor per record: the remaining reactant of greatest molecular weight.

    Hubs (reactants of more than HUB_LIMIT records, except the exemptions) and the
    transfer and redox carriers are removed first. Ties and reactants without a
    recorded weight are resolved in identifier order; records with no remaining
    reactant are omitted, and products left without records are dropped.
    """
    out_degree: Counter[str] = Counter()
    for reactions in product_dict.values():
        for reaction in reactions:
            out_degree.update(reaction.get("reactant", []))
    high_degree_hubs = {
        compound
        for compound, degree in out_degree.items()
        if degree > HUB_LIMIT and compound not in hub_exemptions
    }
    excluded = high_degree_hubs | set(TRANSFER_AND_REDOX_CARRIERS)
    selection_counts: Counter[str] = Counter()
    network: dict[str, list[dict[str, Any]]] = {}
    total_reactions = 0
    for product, reactions in product_dict.items():
        filtered_reactions: list[dict[str, Any]] = []
        for reaction in reactions:
            total_reactions += 1
            original = list(reaction.get("reactant", []))
            candidates = sorted({item for item in original if item not in excluded})
            if not candidates:
                selection_counts["omitted_no_remaining_reactant"] += 1
                continue
            known = []
            for candidate in candidates:
                value = masses.get(candidate)
                if value is None:
                    continue
                try:
                    numeric = float(value)
                except (TypeError, ValueError):
                    continue
                if math.isfinite(numeric):
                    known.append((numeric, candidate))
            if known:
                maximum = max(value for value, _ in known)
                tied = sorted(
                    candidate for value, candidate in known if value == maximum
                )
                selected = tied[0]
                selection_counts["greatest_molecular_weight"] += 1
                if len(tied) > 1:
                    selection_counts["weight_tie_identifier_order"] += 1
            else:
                selected = candidates[0]
                selection_counts["no_recorded_weight_identifier_order"] += 1
            if selected == product:
                selection_counts["selected_self_link"] += 1
            new_reaction = dict(reaction)
            new_reaction["reactant"] = [selected]
            filtered_reactions.append(new_reaction)
        if filtered_reactions:
            network[product] = filtered_reactions
    summary = {
        "product_count": len(product_dict),
        "record_count": total_reactions,
        "hub_count": len(high_degree_hubs),
        "hubs": sorted(high_degree_hubs),
        "hub_exemptions": sorted(hub_exemptions),
        "exempted_hubs": sorted(
            compound
            for compound in hub_exemptions
            if out_degree.get(compound, 0) > HUB_LIMIT
        ),
        "carriers_seen_as_reactants": sorted(
            set(TRANSFER_AND_REDOX_CARRIERS) & set(out_degree)
        ),
        "products_with_precursor": len(network),
        "selection_counts": dict(sorted(selection_counts.items())),
    }
    return network, summary


def read_network_changes(networks: Path, version: str) -> tuple[dict, dict]:
    """The node identities and the added literature connections of one network version."""
    folder = Path(networks) / NETWORK_FOLDERS[version]
    with open(folder / "nodes.jsonl", encoding="utf-8") as handle:
        nodes = {
            row["node_id"]: row
            for row in (json.loads(line) for line in handle if line.strip())
        }
    with open(folder / "changes.json", encoding="utf-8") as handle:
        changes = json.load(handle)
    return nodes, changes


def read_corrections(networks: Path) -> dict:
    """The curated corrections file of the reaction network."""
    with open(Path(networks) / "corrections.json", encoding="utf-8") as handle:
        return json.load(handle)


def excluded_native_reactions(corrections: Mapping) -> set[str]:
    """MetaCyc reactions that a correction excludes or replaces."""
    return {
        row["native_id"]
        for row in corrections["native_decisions"]
        if row["action"] in EXCLUDING_ACTIONS
    }


def build_product_records(metacyc: Path, networks: Path) -> tuple[dict, dict, dict]:
    """Product records of both network versions, molecular weights and pathway counts.

    Returns ({version: products}, masses, summary).
    """
    classes, pathways, reactions, compounds = (
        read_flat_file(Path(metacyc) / name) for name in FLAT_FILES
    )
    domain = pathway_domain(classes, pathways)
    native = product_records(reactions, pathways, domain)
    masses = molecular_weights(compounds)
    excluded = excluded_native_reactions(read_corrections(networks))
    support = supporting_reactions(reactions)
    products = {}
    for version in NETWORK_VERSIONS:
        nodes, changes = read_network_changes(networks, version)
        products[version] = add_literature_connections(
            native, excluded, support, changes, nodes
        )
    summary = {
        "base_pathways": len(domain),
        "metacyc_product_keys": len(native),
        "metacyc_records": sum(map(len, native.values())),
        "excluded_metacyc_reactions": sorted(excluded),
    }
    return products, masses, summary
