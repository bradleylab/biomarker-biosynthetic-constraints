"""Print the numbers on matched controls (Table S9) quoted in the paper.

Counts are by distinct traced state: queries whose identifiers resolve to the
same compound state are counted once.

    python analyses/07_matched_controls/text_numbers.py --results results
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from biosynthesis.tables import read_csv  # noqa: E402

CONTROLS = (
    ("stop_at_starting_compounds", "stopping at starting compounds"),
    ("pathway_reactions_only", "following only pathway-assigned reactions"),
    ("pathway_reactions_and_stop", "pathway-assigned reactions and stopping"),
    ("strict_reactions_broader_compounds", "strict reactions on broader compounds"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    folder = args.results / "controls"
    queries = read_csv(folder / "TableS9_per_query_controls.csv")
    contexts = read_csv(folder / "TableS9_per_context_controls.csv")
    annotations = read_csv(folder / "TableS9_annotation_invariance.csv")
    for version in ("primary", "with_conditional"):
        print(f"Network version: {version}")
        rows = [r for r in annotations if r["network_version"] == version]
        distinct = {(r["reaction_criteria"], r["node"]) for r in rows}
        changed = {
            (r["reaction_criteria"], r["node"])
            for r in rows
            if r["annotation_changes_result"] == "true"
        }
        print(
            f"  Supplying unused annotations changes {len(changed)} of {len(distinct)} distinct "
            f"traced queries."
        )
        for criteria in ("strict", "broader"):
            chosen = [r for r in rows if r["reaction_criteria"] == criteria]
            nodes = {r["node"] for r in chosen}
            print(
                f"  {criteria}: {len(chosen)} query identifiers, {len(nodes)} distinct traced "
                f"states ({len(chosen) - len(nodes)} identifiers resolve to an already listed state)."
            )
        print(f"  Table S9 lists {len(rows)} queries.")
        for criteria in ("strict", "broader"):
            for control, description in CONTROLS:
                chosen = [
                    r
                    for r in queries
                    if r["network_version"] == version
                    and r["reaction_criteria"] == criteria
                    and r["control"] == control
                ]
                if not chosen:
                    continue
                nodes = {r["node"] for r in chosen}
                moved = {
                    r["node"]
                    for r in chosen
                    if r["membership_changed_by_control"] == "true"
                }
                print(
                    f"  {criteria}, {description}: changes {len(moved)} of {len(nodes)} distinct "
                    f"traced queries."
                )
        chosen = [r for r in contexts if r["network_version"] == version]
        distinct = {
            (r["reaction_criteria"], r["node"], r["pathway"]): r for r in chosen
        }
        removed = sum(
            int(r["reactions_removed_by_pathway_restriction"]) > 0
            for r in distinct.values()
        )
        stopped = sum(
            int(r["additional_reactions_removed_by_stopping"]) > 0
            for r in distinct.values()
        )
        print(
            f"  Single pathway: {len(chosen)} rows, {len(distinct)} distinct combinations of traced "
            f"state, pathway and criteria; tracing within the pathway removes reactions for "
            f"{removed}, and also stopping at starting compounds removes further reactions for "
            f"{stopped}."
        )


if __name__ == "__main__":
    main()
