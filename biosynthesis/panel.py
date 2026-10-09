"""The biomarker panel and the compound identities considered for each biomarker."""

from __future__ import annotations

import json
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
IDENTITIES = DATA / "biomarker_identities.jsonl"


def load_identities(path: Path = IDENTITIES) -> list[dict]:
    """Read the identity records, one per biomarker identity, in panel order.

    Each record gives the biomarker number, name and class, the panel compound
    identifier, whether the identity is accepted, a candidate that was not
    accepted, or absent, and the MetaCyc compounds that represent it.
    """
    with open(path, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def load_biomarkers(path: Path = IDENTITIES) -> list[dict]:
    """One record per biomarker: number, name, class and panel compound identifier."""
    seen = {}
    for row in load_identities(path):
        seen.setdefault(
            row["biomarker_number"],
            {
                "biomarker_number": row["biomarker_number"],
                "name": row["name"],
                "class": row["class"],
                "compound_id": row["compound_id"],
            },
        )
    return [seen[number] for number in sorted(seen)]
