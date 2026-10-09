"""Build the reaction networks from the MetaCyc 29.5 flat files and the curated literature records.

Writes, under --out:

- reaction_records.jsonl: one record per MetaCyc reaction;
- base_graph/: the strict and broader base networks (nodes, edges, and the
  screened reactions with the reasons any orientation was not kept) and
  identity_accounting.jsonl, the network nodes of each biomarker identity;
- corrections.json: the curated literature records and MetaCyc decisions as applied;
- strict_primary/, strict_with_conditional/, broader_primary/,
  broader_with_conditional/: the final networks, with changes.json listing what
  the literature records and MetaCyc decisions changed.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))

from biosynthesis import network  # noqa: E402
from biosynthesis.metacyc import read_file  # noqa: E402
from biosynthesis.panel import load_identities  # noqa: E402
from biosynthesis.reactions import build_reaction_records  # noqa: E402

DATA = REPOSITORY / "data"
FINAL_FOLDERS = {
    (network.STRICT, "primary"): "strict_primary",
    (network.STRICT, "with_conditional"): "strict_with_conditional",
    (network.BROADER, "primary"): "broader_primary",
    (network.BROADER, "with_conditional"): "broader_with_conditional",
}


def build(
    metacyc: Path, out: Path, data: Path = DATA, identities: list[dict] | None = None
) -> None:
    """Run the four steps: reaction records, base networks, corrections, final networks."""
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    lines = build_reaction_records(metacyc)
    (out / "reaction_records.jsonl").write_bytes(b"".join(lines))
    del lines
    print(f"reaction records: {time.perf_counter() - started:.0f} s", flush=True)

    started = time.perf_counter()
    compounds = read_file(metacyc, "compounds.dat")
    base_literature = (data / "base_literature_reactions.json").read_bytes()
    base = network.build_base_networks(
        out / "reaction_records.jsonl", compounds, base_literature
    )
    network.write_base_networks(base, out / "base_graph")
    rows = network.identity_accounting(
        identities
        if identities is not None
        else load_identities(data / "biomarker_identities.jsonl"),
        base,
    )
    network.write_rows(out / "base_graph" / "identity_accounting.jsonl", rows)
    del base, rows
    print(f"base networks: {time.perf_counter() - started:.0f} s", flush=True)

    started = time.perf_counter()
    config = network.load_corrections(data, base_literature, compounds)
    network.save_json(out / "corrections.json", config)
    for tier in network.CRITERIA:
        base_network = network.load_base_network(out / "base_graph", tier)
        for version in network.VERSIONS:
            final, changes = network.build_final_network(base_network, config, version)
            network.write_final_network(
                final, changes, out / FINAL_FOLDERS[tier, version]
            )
            print(
                f"{FINAL_FOLDERS[tier, version]}: {len(final.nodes)} nodes, {len(final.edges)} edges",
                flush=True,
            )
    print(f"final networks: {time.perf_counter() - started:.0f} s", flush=True)


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
        "--out", type=Path, required=True, help="output folder, e.g. results/networks"
    )
    parser.add_argument(
        "--data", type=Path, default=DATA, help="folder with the curated data files"
    )
    args = parser.parse_args()
    build(args.metacyc, args.out, args.data)


if __name__ == "__main__":
    main()
