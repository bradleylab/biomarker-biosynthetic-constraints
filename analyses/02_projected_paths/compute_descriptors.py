"""Optional: recalculate the descriptor columns of data/descriptors.csv from structures.

Needs the MetaCyc structures written by extract_structures.py (licensed; not
distributed) and the nine literature structures. Three steps, each in its own
environment:

  rdkit        the 12 RDKit descriptors (RDKit 2023.03.3)
  eligibility  structural eligibility (RDKit 2025.09.4; RDKit 2023.03.3 parses one
               more structure, CPD-19511, and would give 26 ineligible compounds)
  bottcher     Böttcher scores with mesomer handling disabled and a symmetry-search
               limit of 10,000; compounds whose search hit that limit are scored again
               at 3,000,000. Warnings are read from the scorer's messages. Needs Open
               Babel (openbabel-wheel 3.1.1.23) and the scorer file bottchscore3.py of
               forlilab/bottchscore (revision 6f1368cd), given with --scorer.

  python analyses/02_projected_paths/compute_descriptors.py rdkit \\
      --structures structures.tsv --literature-structures data/literature_structures.csv \\
      --out results/descriptor_calculation
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import io
import json
import math
import select
import subprocess
import sys
from contextlib import redirect_stdout
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

from biosynthesis import descriptors as desc  # noqa: E402
from biosynthesis.tables import read_csv, write_csv  # noqa: E402

DATA = REPO / "data"
RDKIT_COLUMNS = [
    "sa_score",
    "bertz_ct",
    "stereocenters",
    "fsp3",
    "ring_count",
    "rotatable_bonds",
    "hb_acceptors",
    "hb_donors",
    "heavy_atoms",
    "tpsa",
    "logp",
    "mw",
]
SCORE_TIMEOUT_SECONDS = 30
STARTUP_TIMEOUT_SECONDS = 120


def structures(
    main_file: Path, literature_file: Path
) -> tuple[dict[str, str], set[str]]:
    """SMILES by compound: MetaCyc structures, then literature structures not among them.

    Also returns the compounds whose structure comes from the literature file; they
    have no Böttcher score.
    """
    with open(main_file, newline="", encoding="utf-8") as handle:
        result = {
            row["Compound"]: row["SMILES"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    added = set()
    for row in read_csv(literature_file):
        if row["compound_id"] not in result:
            result[row["compound_id"]] = row["smiles"]
            added.add(row["compound_id"])
    return result, added


def descriptor_compounds() -> list[str]:
    """The compounds with descriptor records, in the order of data/descriptors.csv."""
    return [row["compound_id"] for row in read_csv(DATA / "descriptors.csv")]


def rdkit_step(smiles: dict[str, str], out: Path) -> None:
    sascorer = desc.sa_scorer()
    rows = []
    for compound in descriptor_compounds():
        values = None
        if compound in smiles:
            values = desc.compute_structural_descriptors(
                str(smiles[compound]), sascorer
            )
        values = values or {column: np.nan for column in RDKIT_COLUMNS}
        rows.append(
            {"compound_id": compound, **{k: float(v) for k, v in values.items()}}
        )
    write_csv(out / "rdkit_descriptors.csv", rows, ["compound_id", *RDKIT_COLUMNS])


def eligibility_step(smiles: dict[str, str], out: Path) -> None:
    rows = [
        {"compound_id": compound, **desc.structural_eligibility(smiles.get(compound))}
        for compound in descriptor_compounds()
    ]
    write_csv(
        out / "structural_eligibility.csv",
        rows,
        ["compound_id", "structurally_eligible", "eligibility_reason"],
    )


def score_worker(scorer_path: str) -> None:
    """Score SMILES read from stdin, recording the symmetry search and failed atomic terms."""
    from openbabel import openbabel as ob

    spec = importlib.util.spec_from_file_location("bottchscore3", scorer_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_find = ob.FindAutomorphisms
    original_atom = module.BottchScore._calculate_complexity
    original_score = module.BottchScore.score
    observations: dict = {}

    def find(*args):
        status = original_find(*args)
        observations["automorphism_search"] = dict(
            returned=bool(status), mappings=len(args[1])
        )
        return status

    def atom(self, idx):
        values = dict(self._indices[idx])
        try:
            math.log(values["Vi"] * values["bi"], 2)
        except (ValueError, OverflowError) as exc:
            observations["failed_terms"].append(
                dict(index=idx, terms=values, error=str(exc))
            )
        return original_atom(self, idx)

    def score(self, mol, disable_mesomeric=False):
        value = original_score(self, mol, disable_mesomeric)
        observations["atom_count"] = mol.NumAtoms()
        return value

    ob.FindAutomorphisms = find
    module.BottchScore._calculate_complexity = atom
    module.BottchScore.score = score
    print(json.dumps(dict(status="ready")), flush=True)
    for line in sys.stdin:
        request = json.loads(line)
        observations.clear()
        observations["failed_terms"] = []
        output = io.StringIO()
        try:
            with redirect_stdout(output):
                value = module.calculate_bottchscore_from_smiles(
                    request["smiles"],
                    disable_mesomer=desc.BOTTCHER_SETTINGS["disable_mesomer"],
                    automorp_memory_maxsize=request["memory_limit"],
                    verbose_response=False,
                    debug_arg=False,
                )
            if not math.isfinite(value):
                raise ValueError("Scoring returned a nonfinite value")
            result = dict(status="scored", score=value, **observations)
        except Exception as exc:
            result = dict(
                status="error", error=f"{type(exc).__name__}: {exc}", **observations
            )
        result["python_messages"] = output.getvalue()
        print(json.dumps(result, allow_nan=False), flush=True)


class Scorer:
    """Run the scorer in a separate process; keep its native messages per compound."""

    def __init__(self, scorer_path: Path, log: Path):
        self.scorer_path = scorer_path
        self.log_path = log
        self.log = open(log, "ab", buffering=0)
        self.process = None

    def close(self) -> None:
        if self.process is not None:
            self.process.kill()
            self.process.wait()
            self.process = None

    def call(self, smiles: str, memory_limit: int) -> dict:
        if self.process is None:
            self.process = subprocess.Popen(
                [
                    sys.executable,
                    "-B",
                    str(Path(__file__)),
                    "score-worker",
                    "--scorer",
                    str(self.scorer_path),
                ],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=self.log,
                text=True,
                bufsize=1,
            )
            ready, _, _ = select.select(
                [self.process.stdout], [], [], STARTUP_TIMEOUT_SECONDS
            )
            if not ready or json.loads(self.process.stdout.readline()) != {
                "status": "ready"
            }:
                self.close()
                raise RuntimeError("The scorer process did not start")
        start = self.log.tell()
        self.process.stdin.write(
            json.dumps({"smiles": smiles, "memory_limit": memory_limit}) + "\n"
        )
        self.process.stdin.flush()
        ready, _, _ = select.select(
            [self.process.stdout], [], [], SCORE_TIMEOUT_SECONDS
        )
        if ready:
            result = json.loads(self.process.stdout.readline())
        else:
            self.close()
            result = dict(status="timeout")
        end = self.log.tell()
        with open(self.log_path, "rb") as handle:
            handle.seek(start)
            result["native_messages"] = handle.read(end - start).decode(
                errors="replace"
            )
        return result


def bottcher_step(
    smiles: dict[str, str], literature: set[str], scorer_path: Path, out: Path
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    scorer = Scorer(scorer_path, out / "scorer_messages.log")
    rows, records = [], []
    try:
        for compound in descriptor_compounds():
            row = {"compound_id": compound}
            if compound not in smiles or compound in literature:
                row.update(bottchscore_unscreened="", bottcher_warnings="")
                rows.append(row)
                continue
            first = scorer.call(
                smiles[compound], desc.BOTTCHER_SETTINGS["automorp_memory_maxsize"]
            )
            records.append({"compound_id": compound, **first})
            unscreened = (
                repr(float(first["score"])) if first["status"] == "scored" else ""
            )
            messages = first["native_messages"] + first.get("python_messages", "")
            warnings = desc.classify_warnings(messages) if messages else []
            replacement = None
            if "automorphism_memory_limit" in warnings:
                larger = scorer.call(smiles[compound], desc.LARGER_SYMMETRY_SEARCH)
                records.append({"compound_id": compound, **larger})
                if (
                    larger["status"] == "scored"
                    and larger.get("automorphism_search", {}).get("returned")
                    and not larger["failed_terms"]
                    and not larger["native_messages"]
                    and not larger["python_messages"]
                ):
                    replacement = repr(float(larger["score"]))
            row.update(
                bottchscore_unscreened=unscreened,
                bottcher_warnings=";".join(warnings),
                bottchscore_larger_search=replacement,
            )
            rows.append(row)
    finally:
        scorer.close()
        scorer.log.close()
    for row in rows:
        status, score = desc.bottcher_status(
            row["bottchscore_unscreened"],
            row["bottcher_warnings"],
            row.get("bottchscore_larger_search"),
        )
        row.update(bottcher_status=status, bottchscore=score)
    write_csv(
        out / "bottcher_scores.csv",
        rows,
        [
            "compound_id",
            "bottchscore",
            "bottchscore_unscreened",
            "bottcher_status",
            "bottcher_warnings",
            "bottchscore_larger_search",
        ],
    )
    with open(out / "bottcher_records.jsonl", "w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, allow_nan=False) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "step", choices=["rdkit", "eligibility", "bottcher", "score-worker"]
    )
    parser.add_argument(
        "--structures", type=Path, help="tab-separated Compound/SMILES file"
    )
    parser.add_argument(
        "--literature-structures", type=Path, default=DATA / "literature_structures.csv"
    )
    parser.add_argument("--scorer", type=Path, help="path to bottchscore3.py")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.step == "score-worker":
        score_worker(str(args.scorer))
        return
    smiles, literature = structures(args.structures, args.literature_structures)
    if args.step in ("rdkit", "eligibility"):
        desc.quiet_rdkit()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.step == "rdkit":
        rdkit_step(smiles, args.out)
    elif args.step == "eligibility":
        eligibility_step(smiles, args.out)
    else:
        bottcher_step(smiles, literature, args.scorer, args.out)


if __name__ == "__main__":
    main()
