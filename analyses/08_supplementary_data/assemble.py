"""Collect the Supplementary Data tables written by the analyses into one folder.

Usage:
    python analyses/08_supplementary_data/assemble.py --results results --out results/supplementary_data
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

TABLES = [
    "Table1_projection_class_means.csv",
    "Table1.md",
    "TableS1_projection.csv",
    "TableS1_tracing.csv",
    "TableS1_tracing_native_control.csv",
    "TableS1_producing_reactions.csv",
    "TableS1_producing_reactions_native_control.csv",
    "TableS2_projection.csv",
    "TableS2_population_summaries.csv",
    "TableS2_zero_path_causes.csv",
    "TableS2_exemption_sensitivity.csv",
    "TableS3_pair_comparisons.csv",
    "TableS4_associations.csv",
    "TableS4_descriptor_cohort.csv",
    "TableS4_descriptor_selection.csv",
    "TableS4_score_treatment_summary.csv",
    "TableS5_source_observations.csv",
    "TableS6_analyte_precursor_map.csv",
    "TableS6_citations.csv",
    "TableS6_compatibility_links.csv",
    "TableS6_representative_targets.csv",
    "TableS6_source_compatibility.csv",
    "TableS6_source_profiles.csv",
    "TableS7_additions.csv",
    "TableS7_analysis_states.csv",
    "TableS7_compatible_combinations.csv",
    "TableS7_delegation_by_sample.csv",
    "TableS7_evidence_removal.csv",
    "TableS7_minimal_combination_assignments.csv",
    "TableS7_omissions.csv",
    "TableS7_scenarios.csv",
    "TableS7_thresholds_by_sample.csv",
    "TableS8_literature_reactions.csv",
    "TableS8_metacyc_corrections.csv",
    "TableS8_native_control.csv",
    "TableS9_annotation_invariance.csv",
    "TableS9_per_context_controls.csv",
    "TableS9_per_query_controls.csv",
]


def find_tables(results: Path, out: Path) -> dict[str, Path]:
    """Locate each table under the results folder, outside the output folder."""
    found: dict[str, list[Path]] = {name: [] for name in TABLES}
    for path in sorted(results.rglob("*")):
        if path.name in found and path.is_file() and out not in path.parents:
            found[path.name].append(path)
    missing = [name for name, paths in found.items() if not paths]
    if missing:
        raise SystemExit("Tables not found: " + ", ".join(missing))
    repeated = [name for name, paths in found.items() if len(paths) > 1]
    if repeated:
        raise SystemExit("Tables found more than once: " + ", ".join(repeated))
    return {name: paths[0] for name, paths in found.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    tables = find_tables(args.results.resolve(), out)
    out.mkdir(parents=True, exist_ok=True)
    for name, source in tables.items():
        shutil.copyfile(source, out / name)
    print(f"{len(tables)} tables written to {args.out}")


if __name__ == "__main__":
    main()
