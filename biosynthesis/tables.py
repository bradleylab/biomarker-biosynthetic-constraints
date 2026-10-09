"""Reading and writing the CSV and JSON files of the analyses."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

csv.field_size_limit(2**31 - 1)


def cell(value: object) -> str:
    """Format one value: blank for missing, lower-case booleans, integers without '.0'.

    NumPy scalars are converted to the equivalent Python values first.
    """
    if value is None:
        return ""
    if hasattr(value, "item") and not isinstance(value, (list, dict, tuple)):
        value = value.item()
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return repr(value)
    if isinstance(value, tuple):
        value = list(value)
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def write_csv(
    path: Path, rows: Iterable[Mapping[str, object]], columns: Sequence[str]
) -> int:
    """Write rows with the given columns and LF line endings; return the row count."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([cell(row.get(column)) for column in columns])
            count += 1
    return count


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV file into a list of dictionaries of strings."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, value: object) -> None:
    """Write JSON with sorted keys and two-space indentation."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def read_jsonl(path: Path) -> list[dict]:
    """Read a JSON-lines file."""
    with open(path, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def write_jsonl(path: Path, rows: Iterable[object]) -> int:
    """Write one compact JSON value per line, with sorted keys."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                json.dumps(
                    row, sort_keys=True, separators=(",", ":"), ensure_ascii=False
                )
            )
            handle.write("\n")
            count += 1
    return count
