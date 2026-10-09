"""Table S8: the literature reaction records, their connections and how many enter each network.

    python analyses/06_literature_records/literature_reactions.py --networks results/networks \\
        --out results/literature
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.graph import BROADER, FOLDERS  # noqa: E402
from biosynthesis.tables import write_csv  # noqa: E402

COLUMNS = (
    "record_id",
    "disposition",
    "connections_primary_graph",
    "connections_with_conditional",
    "conditional_connections",
    "substrate_product_connections_json",
    "finding",
    "remaining_scientific_uncertainty",
    "unresolved_requirements_json",
    "source_dois",
    "source_locators",
)


def connections(conversions: list[dict]) -> str:
    """The substrates, products, connection type and enzyme of each connection, as JSON."""
    return json.dumps(
        [
            {
                "inputs": row["inputs"],
                "outputs": row["outputs"],
                "connection_type": row["connection_type"],
                "enzyme": row.get("enzyme"),
            }
            for row in conversions
        ],
        sort_keys=True,
    )


def requirements(conversions: list[dict]) -> str:
    """The unresolved requirements of each connection, as JSON."""
    return json.dumps(
        [row.get("unresolved", []) for row in conversions],
        sort_keys=True,
    )


def record_rows(config: dict, inventories: dict) -> list[dict]:
    """One row per literature record, with the connections that enter each broader network."""
    output = []
    for record in config["records"]:
        conversions = record["conversions"]
        supported = [
            row for row in conversions if row["connection_type"] == "supported"
        ]
        conditional = [
            row for row in conversions if row["connection_type"] == "conditional"
        ]
        disposition = (
            "supported_and_conditional"
            if supported and conditional
            else (
                "supported"
                if supported
                else "conditional"
                if conditional
                else "observation_only"
            )
        )
        edge_ids = {}
        # Inventory keys as stored in the network change records.
        for version, inventory in inventories.items():
            entry = inventory[record["record_id"]]
            edge_ids[version] = sorted(
                set(entry["kept_edge_ids"] + entry["retained_native_evidence_edge_ids"])
            )
        output.append(
            {
                "record_id": record["record_id"],
                "disposition": disposition,
                "connections_primary_graph": len(edge_ids["primary"]),
                "connections_with_conditional": len(edge_ids["with_conditional"]),
                "conditional_connections": len(conditional),
                "substrate_product_connections_json": connections(conversions),
                "finding": record["finding"],
                "remaining_scientific_uncertainty": "; ".join(
                    record.get("remaining", [])
                ),
                "unresolved_requirements_json": requirements(conversions),
                "source_dois": "; ".join(
                    sorted(
                        {
                            source["doi"]
                            for source in record["sources"]
                            if source.get("doi")
                        }
                    )
                ),
                "source_locators": "; ".join(
                    source["locator"] for source in record["sources"]
                ),
            }
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--networks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads((args.networks / "corrections.json").read_text())
    inventories = {}
    for version in ("primary", "with_conditional"):
        changes = json.loads(
            (args.networks / FOLDERS[BROADER, version] / "changes.json").read_text()
        )
        inventories[version] = {row["record_id"]: row for row in changes["inventory"]}
    rows = record_rows(config, inventories)
    write_csv(args.out / "TableS8_literature_reactions.csv", rows, COLUMNS)
    print(f"Wrote {len(rows)} literature records", flush=True)


if __name__ == "__main__":
    main()
