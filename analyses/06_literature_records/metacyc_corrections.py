"""Write the Table S8 files on MetaCyc reactions: corrections and the single-reaction test.

TableS8_metacyc_corrections.csv lists every MetaCyc reaction that the curated
decisions replace, exclude, or keep with unresolved requirements recorded
(reactions kept unchanged are not listed): the equation as recorded in MetaCyc
(compound common names, in the recorded direction), the action, the reason and
the DOIs of the sources.

TableS8_native_control.csv describes the three MetaCyc hopanoid reactions
RXN-13529, RXN-13530 and RXN-13531, how each is treated in the networks, and
the DOIs and locators of the sources.

The equations come from the user's MetaCyc flat files, so the first file is
subject to the MetaCyc license.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))

from biosynthesis.tables import write_csv  # noqa: E402

DATA = REPOSITORY / "data"
ACTIONS = {
    "replace": "Replaced",
    "exclude": "Excluded",
    "annotate": "Retained; unresolved requirements recorded",
}
CORRECTION_COLUMNS = [
    "metacyc_reaction",
    "recorded_equation",
    "action",
    "reason",
    "source_dois",
]
TEST_COLUMNS = [
    "reaction_id",
    "chemical_evidence",
    "network_treatment",
    "source_dois",
    "source_locators",
]
TEST_REACTIONS = ("RXN-13529", "RXN-13530", "RXN-13531")


def metacyc_records(metacyc: Path, name: str) -> dict[str, dict[str, list[str]]]:
    """Read the fields of each record of a flat file, keyed by UNIQUE-ID."""
    records = {}
    for block in (metacyc / name).read_text(encoding="latin-1").split("\n//\n"):
        fields: dict[str, list[str]] = {}
        for line in block.splitlines():
            if " - " in line and not line.startswith(("#", "/")):
                key, value = line.split(" - ", 1)
                fields.setdefault(key, []).append(value)
        if fields.get("UNIQUE-ID"):
            records[fields["UNIQUE-ID"][0]] = fields
    return records


def equation(reaction: dict, names: dict) -> str:
    """Write a reaction as 'left → right' in its recorded direction."""

    def side(ids: list[str]) -> str:
        return " + ".join(names.get(i, i) for i in ids)

    left, right = reaction.get("LEFT", []), reaction.get("RIGHT", [])
    if reaction.get("REACTION-DIRECTION", [""])[0].endswith("RIGHT-TO-LEFT"):
        left, right = right, left
    return f"{side(left)} → {side(right)}"


def correction_rows(metacyc: Path, decisions: list[dict]) -> list[dict]:
    """One row per replaced, excluded or annotated MetaCyc reaction."""
    reactions = metacyc_records(metacyc, "reactions.dat")
    names = {
        key: html.unescape(re.sub("<[^>]+>", "", record.get("COMMON-NAME", [key])[0]))
        for key, record in metacyc_records(metacyc, "compounds.dat").items()
    }
    listed = sorted(
        (row for row in decisions if row["action"] in ACTIONS),
        key=lambda row: (list(ACTIONS).index(row["action"]), row["native_id"]),
    )
    return [
        {
            "metacyc_reaction": row["native_id"],
            "recorded_equation": equation(reactions[row["native_id"]], names),
            "action": ACTIONS[row["action"]],
            "reason": row["table_reason"],
            "source_dois": "; ".join(
                sorted({s["doi"] for s in row.get("sources", []) if s.get("doi")})
            ),
        }
        for row in listed
    ]


def test_rows(decisions: list[dict]) -> list[dict]:
    """One row per reaction of the single-reaction test, with its evidence and treatment."""
    by_id = {row["native_id"]: row for row in decisions}
    rows = []
    for reaction in TEST_REACTIONS:
        decision = by_id[reaction]
        sources = decision["sources"]
        rows.append(
            {
                "reaction_id": reaction,
                "chemical_evidence": decision["single_reaction_test"][
                    "chemical_evidence"
                ],
                "network_treatment": decision["single_reaction_test"][
                    "network_treatment"
                ],
                "source_dois": "; ".join(
                    sorted({str(row["doi"]) for row in sources if row["doi"]})
                ),
                "source_locators": "; ".join(row["locator"] for row in sources),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--metacyc",
        type=Path,
        required=True,
        help="folder with the MetaCyc 29.5 flat files",
    )
    parser.add_argument(
        "--out", type=Path, required=True, help="output folder, e.g. results/literature"
    )
    parser.add_argument(
        "--data", type=Path, default=DATA, help="folder with metacyc_corrections.json"
    )
    args = parser.parse_args()
    decisions = json.loads(
        (args.data / "metacyc_corrections.json").read_text(encoding="utf-8")
    )["reactions"]
    count = write_csv(
        args.out / "TableS8_metacyc_corrections.csv",
        correction_rows(args.metacyc, decisions),
        CORRECTION_COLUMNS,
    )
    write_csv(
        args.out / "TableS8_native_control.csv", test_rows(decisions), TEST_COLUMNS
    )
    print(
        f"{count} MetaCyc corrections and {len(TEST_REACTIONS)} single-reaction test rows written to {args.out}"
    )


if __name__ == "__main__":
    main()
