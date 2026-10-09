"""Optional: read compound structures from a MetaCyc database through Pathway Tools.

For every instance of the MetaCyc Compounds class, reads the SMILES, InChI,
common-name and molecular-weight slots through the PythonCyc interface of a
running Pathway Tools server, and writes a tab-separated file with the columns
Compound, SMILES, inchi, common-name and mw. Compounds with neither a SMILES nor
an InChI are left out. The structures are MetaCyc data and are covered by its
license; they are not distributed with this repository. compute_descriptors.py
reads this file.

Needs Pathway Tools with the MetaCyc 29.5 database, started with its Python API
enabled, and the pythoncyc package.

  python analyses/02_projected_paths/extract_structures.py --out structures.tsv
"""

from __future__ import annotations

import argparse
import signal
from pathlib import Path

SLOTS = ("SMILES", "inchi", "common-name", "molecular-weight")
SLOT_TIMEOUT_SECONDS = 5
COLUMNS = ("Compound", "SMILES", "inchi", "common-name", "mw")


class SlotTimeout(Exception):
    pass


def slot_value(database, compound: str, slot: str):
    """One slot of a compound, or None when it is absent or takes longer than the limit."""

    def timeout(signum, frame):
        raise SlotTimeout()

    previous = signal.signal(signal.SIGALRM, timeout)
    signal.alarm(SLOT_TIMEOUT_SECONDS)
    try:
        return database.get_slot_value(compound, slot)
    except Exception:
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def compound_record(database, compound_id: str) -> dict | None:
    """The structure record of one compound, or None without SMILES and InChI."""
    compound = compound_id.replace("|", "")
    values = {slot: slot_value(database, compound, slot) for slot in SLOTS}
    if not values["SMILES"] and not values["inchi"]:
        return None
    return {
        "Compound": compound,
        "SMILES": values["SMILES"],
        "inchi": values["inchi"],
        "common-name": values["common-name"] or "",
        "mw": values["molecular-weight"],
    }


def main(out: Path, organism: str) -> None:
    import pythoncyc

    database = pythoncyc.select_organism(organism)
    compounds = [
        str(getattr(instance, "frameid", instance)).replace("|", "")
        for instance in database["Compounds"].instances
    ]
    written = 0
    with open(out, "w", encoding="utf-8") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for compound in compounds:
            record = compound_record(database, compound)
            if record is None:
                continue
            handle.write("\t".join(str(record[column]) for column in COLUMNS) + "\n")
            written += 1
    print(f"{written} of {len(compounds)} compounds have a structure")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--organism", default="meta", help="Pathway Tools database (default: meta)"
    )
    args = parser.parse_args()
    main(args.out, args.organism)
