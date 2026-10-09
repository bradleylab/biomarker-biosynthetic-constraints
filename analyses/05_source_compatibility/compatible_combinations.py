"""Compatible combinations of source categories for every sample scenario (Table S7, first four files).

Reads Table S5 and the assignment tables in data/geology/, sweeps the detection
threshold, the low-signal treatment and the strict and broader assignments for
every sample, and analyzes each distinct set of detected classes once,
including the sensitivity tests. Writes:

- sample_scenarios.jsonl: one record per sample scenario;
- detected_class_sets.jsonl: one record per set of detected classes and
  assignments, with all combinations and sensitivity tests;
- samples.jsonl: the workbook rows, labels and ages of each sample;
- TableS7_scenarios.csv, TableS7_analysis_states.csv,
  TableS7_compatible_combinations.csv and
  TableS7_minimal_combination_assignments.csv.

Usage:
  python analyses/05_source_compatibility/compatible_combinations.py --observations results/sources --out results/sources
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis.sensitivity import (  # noqa: E402
    CriteriaAssignment,
    Observation,
    RouteState,
    assign_identifiers,
    scenarios_and_analyses,
)
from biosynthesis.sources import (  # noqa: E402
    STRICT,
    BROADER,
    AssignmentOption,
    assignment_options,
    load_links,
    load_source_categories,
)
from biosynthesis.tables import read_csv, write_csv, write_jsonl  # noqa: E402

TABLE_S5 = "TableS5_source_observations.csv"
NUMERIC_TYPES = frozenset(
    {"reported_zero", "reported_positive", "reported_positive_low_signal"}
)
SAMPLE_FIELDS = (
    "workbook_id",
    "sheet_name",
    "sample_key",
    "source_sample_id_text",
    "source_sample_id_coordinate",
    "geological_era_text",
    "formation_or_basin_text",
    "age_source_coordinate",
    "age_raw_text",
    "age_kind",
    "age_parse_status",
)


def _bool(value: str) -> bool:
    if value not in {"true", "false"}:
        raise ValueError(f"expected true or false, found {value!r}")
    return value == "true"


def _optional_float(value: str) -> float | None:
    return float(value) if value.strip() else None


def read_observations(path: Path):
    """Samples, fossil classes used, numeric values and sample records from Table S5."""
    rows = read_csv(path)
    sample_rows = [row for row in rows if _bool(row["included_sample"])]
    sample_ids = tuple(sorted({row["canonical_sample_key"] for row in sample_rows}))
    source_records: dict[tuple[str, str], dict[str, object]] = {}
    for row in sample_rows:
        record = {field: row[field] for field in SAMPLE_FIELDS}
        record.update(
            {
                "source_row": int(row["source_row"]),
                "source_rows": row["source_rows"].split(";"),
                "age_lower_ma": _optional_float(row["age_lower_ma"]),
                "age_upper_ma": _optional_float(row["age_upper_ma"]),
                "age_central_ma": _optional_float(row["age_central_ma"]),
                "age_operator": row["age_operator"] or None,
            }
        )
        previous = source_records.setdefault(
            (row["canonical_sample_key"], row["sample_key"]), record
        )
        if previous != record:
            raise ValueError(
                f"sample record differs between fossil classes: {row['sample_key']}"
            )
    samples: dict[str, list[dict[str, object]]] = {
        sample_id: [] for sample_id in sample_ids
    }
    for (sample_id, _), record in sorted(source_records.items()):
        samples[sample_id].append(record)

    used_rows = [row for row in sample_rows if _bool(row["value_used"])]
    analyte_ids = tuple(sorted({row["analyte_id"] for row in used_rows}))
    observations: list[Observation] = []
    seen: set[tuple[str, str]] = set()
    for row in used_rows:
        state = row["value_type"]
        if state not in NUMERIC_TYPES:
            continue
        low_signal = _bool(row["low_signal"])
        if low_signal != (state == "reported_positive_low_signal"):
            raise ValueError(
                f"low-signal flag differs from value type: {row['source_coordinate']}"
            )
        key = (row["canonical_sample_key"], row["analyte_id"])
        if key in seen:
            raise ValueError(f"two values for one sample and fossil class: {key}")
        seen.add(key)
        observations.append(
            Observation(
                sample_id=key[0],
                analyte_id=key[1],
                measurement_channel=row["measurement_channel"],
                denominator_group=row["denominator_group"],
                reported_abundance=float(row["reported_abundance"]),
                low_signal=low_signal,
            )
        )
    observations.sort(key=lambda row: (row.sample_id, row.group_id, row.analyte_id))
    return (
        sample_ids,
        analyte_ids,
        tuple(observations),
        [
            {"sample_id": sample_id, "source_records": records}
            for sample_id, records in samples.items()
        ],
    )


def read_route_states(path: Path) -> tuple[RouteState, ...]:
    """The representative reaction target of each precursor class (data/representative_targets.csv)."""
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    states = tuple(
        RouteState(
            route_record_id=row["route_record_id"],
            precursor_id=row["precursor_id"],
            state=row["route_state"],
            reason_code=row["route_reason"],
        )
        for row in sorted(rows, key=lambda row: row["route_record_id"])
    )
    if not states or any(row.state not in {"supported", "frontier"} for row in states):
        raise ValueError(
            "each representative target needs a route state of supported or frontier"
        )
    return states


def criteria_assignments(data: Path, route_states) -> tuple[CriteriaAssignment, ...]:
    """Every broader assignment option, with its representative target and its criteria."""
    mappings = load_links(data / "geology" / "fossil_to_precursor.csv", "mapping")
    compatibilities = load_links(
        data / "geology" / "precursor_to_source.csv", "compatibility"
    )
    strict = set(assignment_options(mappings, compatibilities, criteria=STRICT))
    broader = set(assignment_options(mappings, compatibilities, criteria=BROADER))
    route_by_precursor = {row.precursor_id: row.route_record_id for row in route_states}
    assignments = []
    for option in sorted(broader):
        route_record_id = route_by_precursor.get(option.precursor_id)
        if route_record_id is None:
            raise ValueError(f"no representative target for {option.precursor_id}")
        assignments.append(
            CriteriaAssignment(
                option=AssignmentOption(
                    **{**asdict(option), "route_record_id": route_record_id}
                ),
                criteria=STRICT if option in strict else BROADER,
            )
        )
    return tuple(assignments)


def joined(values) -> str:
    return ";".join(values)


def scenario_table(scenarios, analyses_by_id):
    for scenario in scenarios:
        analysis = analyses_by_id.get(scenario["analysis_id"])
        result = analysis["result"] if analysis else {}
        yield {
            **scenario,
            "detected_analyte_ids": joined(scenario["detected_analyte_ids"]),
            "status": result.get("status"),
            "compatible_combinations": result.get("compatible_combinations"),
            "min_categories": result.get("min_categories"),
            "single_category_combinations": joined(
                result.get("single_category_combinations", ())
            ),
        }


def analysis_tables(analyses):
    states, combinations, supporting = [], [], []
    for analysis in analyses:
        result = analysis["result"]
        common = {
            "analysis_id": analysis["analysis_id"],
            "assignment_criteria": analysis["assignment_criteria"],
            "detected_analyte_ids": joined(analysis["detected_analyte_ids"]),
        }
        states.append(
            {
                **common,
                "status": result["status"],
                "compatible_combinations": result["compatible_combinations"],
                "min_categories": result["min_categories"],
                "reason_code": result["reason_code"],
            }
        )
        for combination in result["combinations"]:
            combinations.append(
                {
                    **common,
                    "source_class_ids": joined(combination),
                    "minimal": combination in result["minimal_combinations"],
                    "smallest": combination in result["smallest_combinations"],
                }
            )
        for item in result["minimal_combination_assignments"]:
            supporting.append(
                {
                    **common,
                    "presence_profile_id": analysis["presence_profile_id"],
                    "source_class_ids": joined(item["source_class_ids"]),
                    "assignments": [dict(row) for row in item["assignments"]],
                }
            )
    return states, combinations, supporting


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--observations", type=Path, required=True, help=f"Folder holding {TABLE_S5}."
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder.")
    parser.add_argument(
        "--data", type=Path, default=REPO / "data", help="The repository data folder."
    )
    args = parser.parse_args()

    settings = yaml.safe_load(
        (args.data / "geology" / "settings.yaml").read_text(encoding="utf-8")
    )
    sample_ids, analyte_ids, observations, samples = read_observations(
        args.observations / TABLE_S5
    )
    route_states = read_route_states(args.data / "representative_targets.csv")
    scenarios, analyses = scenarios_and_analyses(
        sample_ids=sample_ids,
        observations=observations,
        assignments=criteria_assignments(args.data, route_states),
        source_categories=load_source_categories(
            args.data / "geology" / "assignment_criteria.yaml"
        ),
        all_analyte_ids=analyte_ids,
        route_states=route_states,
        low_signal_treatments=settings["low_signal_treatments"],
        criteria_list=settings["assignment_criteria"],
    )
    scenarios, analyses = assign_identifiers(scenarios, analyses)
    analyses_by_id = {row["analysis_id"]: row for row in analyses}

    out = args.out
    write_jsonl(out / "sample_scenarios.jsonl", scenarios)
    write_jsonl(out / "detected_class_sets.jsonl", analyses)
    write_jsonl(out / "samples.jsonl", samples)
    write_csv(
        out / "TableS7_scenarios.csv",
        scenario_table(scenarios, analyses_by_id),
        (
            "scenario_id",
            "sample_id",
            "low_signal_mode",
            "assignment_criteria",
            "presence_profile_id",
            "threshold_regions",
            "detected_analyte_ids",
            "analysis_id",
            "status",
            "compatible_combinations",
            "min_categories",
            "single_category_combinations",
        ),
    )
    states, combinations, supporting = analysis_tables(analyses)
    common = ("analysis_id", "assignment_criteria", "detected_analyte_ids")
    write_csv(
        out / "TableS7_analysis_states.csv",
        states,
        (*common, "status", "compatible_combinations", "min_categories", "reason_code"),
    )
    write_csv(
        out / "TableS7_compatible_combinations.csv",
        combinations,
        (*common, "source_class_ids", "minimal", "smallest"),
    )
    write_csv(
        out / "TableS7_minimal_combination_assignments.csv",
        supporting,
        (
            "analysis_id",
            "presence_profile_id",
            "assignment_criteria",
            "detected_analyte_ids",
            "source_class_ids",
            "assignments",
        ),
    )
    status = Counter(
        analyses_by_id[row["analysis_id"]]["result"]["status"] for row in scenarios
    )
    print(
        f"{len(sample_ids)} samples; {len(scenarios)} scenarios "
        f"({', '.join(f'{count} {name}' for name, count in sorted(status.items()))}); "
        f"{len(analyses)} analyses of {len({row['presence_profile_id'] for row in analyses})} sets of detected classes."
    )


if __name__ == "__main__":
    main()
