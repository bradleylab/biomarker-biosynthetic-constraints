"""Reading the MetaCyc attribute-value flat files.

The reader works on bytes and keeps every field and annotation with its byte
range in the file, so that each reaction record can point back to the exact
source text. Record separators (``//``), value continuations (lines beginning
with ``/``), annotations (names beginning with ``^``) and comments (lines
beginning with ``#``) are handled as in the Pathway Tools export format.
"""

from __future__ import annotations

import base64
import hashlib
import re
from dataclasses import dataclass, field
from io import BytesIO
from pathlib import Path
from typing import BinaryIO, Iterator, Literal, Sequence

from biosynthesis.records import line_content_id

FIELD_VALUE_SEPARATOR = b" - "
RECORD_SEPARATOR = b"//"

# Stored in the network records; part of the record identifiers.
COMPOUND_RECORD_SCHEMA_VERSION = "metacyc-control-covariate-record/v1"
TIER_IX = "IX_conflicting_or_invalid_source_fields"
TIER_I1 = "I1_inchi_key_formula_complete"
TIER_I2 = "I2_structure_field_incomplete"
TIER_I3 = "I3_formula_only"
TIER_I4 = "I4_source_frame_only"

FORMULA_TOKEN_GRAMMAR = r"^\(([A-Z]{1,2}) ([1-9][0-9]*)\)$"
FORMULA_ELEMENTS = frozenset(
    {
        "AG",
        "AL",
        "AS",
        "AU",
        "B",
        "BA",
        "BE",
        "BR",
        "C",
        "CA",
        "CD",
        "CL",
        "CO",
        "CR",
        "CS",
        "CU",
        "F",
        "FE",
        "H",
        "HG",
        "I",
        "K",
        "LA",
        "LI",
        "MG",
        "MN",
        "MO",
        "N",
        "NA",
        "NB",
        "NI",
        "O",
        "P",
        "PB",
        "PD",
        "PT",
        "RB",
        "S",
        "SB",
        "SE",
        "SI",
        "SN",
        "SR",
        "TE",
        "TI",
        "TL",
        "U",
        "V",
        "W",
        "ZN",
    }
)
_FORMULA_TOKEN_RE = re.compile(FORMULA_TOKEN_GRAMMAR.encode("ascii"))
_STRUCTURE_FIELDS = (b"INCHI", b"INCHI-KEY", b"SMILES", b"CHEMICAL-FORMULA")


@dataclass(frozen=True, slots=True)
class SourceSpan:
    """A half-open byte range and an inclusive line range in one file."""

    source_member_path: str
    start_byte: int
    end_byte: int
    start_line: int
    end_line: int

    @property
    def locator(self) -> str:
        return f"{self.source_member_path}#bytes={self.start_byte}:{self.end_byte}"


@dataclass(frozen=True, slots=True)
class FieldOccurrence:
    """One field of a record, with its raw name and value bytes."""

    field_occurrence_id: str
    ordinal: int
    field_name_raw: bytes
    value_raw: bytes
    span: SourceSpan
    raw_bytes: bytes
    annotation_occurrence_ids: tuple[str, ...] = ()

    @property
    def start_byte(self) -> int:
        return self.span.start_byte

    @property
    def end_byte(self) -> int:
        return self.span.end_byte

    @property
    def start_line(self) -> int:
        return self.span.start_line

    @property
    def locator(self) -> str:
        return self.span.locator


@dataclass(frozen=True, slots=True)
class AnnotationOccurrence:
    """One annotation (``^NAME - value``) attached to the field before it."""

    annotation_occurrence_id: str
    ordinal: int
    annotation_name_raw: bytes
    value_raw: bytes
    span: SourceSpan
    raw_bytes: bytes
    attached_to_field_occurrence_id: str

    @property
    def start_byte(self) -> int:
        return self.span.start_byte

    @property
    def end_byte(self) -> int:
        return self.span.end_byte

    @property
    def start_line(self) -> int:
        return self.span.start_line


@dataclass(frozen=True, slots=True)
class FlatFileRecord:
    """One record of a flat file with its fields and annotations in file order."""

    ordinal: int
    span: SourceSpan
    raw_bytes: bytes
    field_occurrences: tuple[FieldOccurrence, ...]
    annotation_occurrences: tuple[AnnotationOccurrence, ...]

    @property
    def source_member_path(self) -> str:
        return self.span.source_member_path

    @property
    def record_locator(self) -> str:
        return self.span.locator

    @property
    def start_byte(self) -> int:
        return self.span.start_byte

    @property
    def end_byte(self) -> int:
        return self.span.end_byte

    @property
    def start_line(self) -> int:
        return self.span.start_line

    @property
    def end_line(self) -> int:
        return self.span.end_line


@dataclass(slots=True)
class _OccurrenceBuilder:
    name_raw: bytes
    value_parts: list[bytes]
    raw_parts: list[bytes]
    start_byte: int
    end_byte: int
    start_line: int
    end_line: int
    last_line_ending: bytes

    def append_continuation(
        self,
        *,
        content: bytes,
        line_ending: bytes,
        raw_line: bytes,
        end_byte: int,
        line_number: int,
    ) -> None:
        self.value_parts.extend((self.last_line_ending, content))
        self.raw_parts.append(raw_line)
        self.end_byte = end_byte
        self.end_line = line_number
        self.last_line_ending = line_ending

    @property
    def value_raw(self) -> bytes:
        return b"".join(self.value_parts)

    @property
    def raw_bytes(self) -> bytes:
        return b"".join(self.raw_parts)

    def span(self, source_member_path: str) -> SourceSpan:
        return SourceSpan(
            source_member_path,
            self.start_byte,
            self.end_byte,
            self.start_line,
            self.end_line,
        )


@dataclass(slots=True)
class _FieldBuilder:
    occurrence: _OccurrenceBuilder
    annotations: list[_OccurrenceBuilder] = field(default_factory=list)


def _split_line_ending(raw_line: bytes) -> tuple[bytes, bytes]:
    if raw_line.endswith(b"\r\n"):
        return raw_line[:-2], b"\r\n"
    if raw_line.endswith(b"\n"):
        return raw_line[:-1], b"\n"
    return raw_line, b""


def _new_occurrence(
    *,
    name_raw: bytes,
    value_raw: bytes,
    raw_line: bytes,
    line_ending: bytes,
    start_byte: int,
    end_byte: int,
    line_number: int,
) -> _OccurrenceBuilder:
    return _OccurrenceBuilder(
        name_raw=name_raw,
        value_parts=[value_raw],
        raw_parts=[raw_line],
        start_byte=start_byte,
        end_byte=end_byte,
        start_line=line_number,
        end_line=line_number,
        last_line_ending=line_ending,
    )


def _materialize_record(
    *,
    ordinal: int,
    source_member_path: str,
    start_byte: int,
    end_byte: int,
    start_line: int,
    end_line: int,
    raw_bytes: bytes,
    field_builders: list[_FieldBuilder],
) -> FlatFileRecord:
    fields: list[FieldOccurrence] = []
    annotations: list[AnnotationOccurrence] = []
    annotation_ordinal = 0
    for field_ordinal, field_builder in enumerate(field_builders, start=1):
        occurrence = field_builder.occurrence
        field_span = occurrence.span(source_member_path)
        field_id = f"field:{field_span.locator}"
        field_annotation_ids: list[str] = []
        for annotation_builder in field_builder.annotations:
            annotation_ordinal += 1
            annotation_span = annotation_builder.span(source_member_path)
            annotation_id = f"annotation:{annotation_span.locator}"
            field_annotation_ids.append(annotation_id)
            annotations.append(
                AnnotationOccurrence(
                    annotation_occurrence_id=annotation_id,
                    ordinal=annotation_ordinal,
                    annotation_name_raw=annotation_builder.name_raw,
                    value_raw=annotation_builder.value_raw,
                    span=annotation_span,
                    raw_bytes=annotation_builder.raw_bytes,
                    attached_to_field_occurrence_id=field_id,
                )
            )
        fields.append(
            FieldOccurrence(
                field_occurrence_id=field_id,
                ordinal=field_ordinal,
                field_name_raw=occurrence.name_raw,
                value_raw=occurrence.value_raw,
                span=field_span,
                raw_bytes=occurrence.raw_bytes,
                annotation_occurrence_ids=tuple(field_annotation_ids),
            )
        )
    return FlatFileRecord(
        ordinal=ordinal,
        span=SourceSpan(source_member_path, start_byte, end_byte, start_line, end_line),
        raw_bytes=raw_bytes,
        field_occurrences=tuple(fields),
        annotation_occurrences=tuple(annotations),
    )


def iter_records(
    stream: BinaryIO, *, source_member_path: str
) -> Iterator[FlatFileRecord]:
    """Yield the records of one flat file in file order.

    ``source_member_path`` is the file name used in every byte locator
    (for example ``reactions.dat``).
    """
    byte_offset = 0
    line_number = 0
    record_ordinal = 0
    record_start_byte: int | None = None
    record_start_line: int | None = None
    record_parts: list[bytes] = []
    field_builders: list[_FieldBuilder] = []
    last_field: _FieldBuilder | None = None
    pending: _OccurrenceBuilder | None = None
    pending_kind: Literal["field", "annotation"] | None = None

    def fail(problem: str) -> None:
        raise ValueError(f"{problem} at {source_member_path}:{line_number}")

    def finish_pending() -> None:
        nonlocal last_field, pending, pending_kind
        if pending is None:
            return
        if pending_kind == "field":
            last_field = _FieldBuilder(pending)
            field_builders.append(last_field)
        else:
            if last_field is None:
                fail("annotation without a field")
            last_field.annotations.append(pending)
        pending = None
        pending_kind = None

    while True:
        raw_line = stream.readline()
        if raw_line == b"":
            break
        line_number += 1
        line_start_byte = byte_offset
        byte_offset += len(raw_line)
        content, line_ending = _split_line_ending(raw_line)

        if content == RECORD_SEPARATOR:
            finish_pending()
            if record_start_byte is None or record_start_line is None:
                fail("empty record")
            record_parts.append(raw_line)
            record_ordinal += 1
            record = _materialize_record(
                ordinal=record_ordinal,
                source_member_path=source_member_path,
                start_byte=record_start_byte,
                end_byte=byte_offset,
                start_line=record_start_line,
                end_line=line_number,
                raw_bytes=b"".join(record_parts),
                field_builders=field_builders,
            )
            record_start_byte = None
            record_start_line = None
            record_parts = []
            field_builders = []
            last_field = None
            pending = None
            pending_kind = None
            yield record
            continue

        if content.startswith(b"#"):
            finish_pending()
            if record_start_byte is not None:
                record_parts.append(raw_line)
            continue

        if content.startswith(b"/"):
            if pending is None:
                fail("continuation line without a value")
            pending.append_continuation(
                content=content,
                line_ending=line_ending,
                raw_line=raw_line,
                end_byte=byte_offset,
                line_number=line_number,
            )
            record_parts.append(raw_line)
            continue

        if content.startswith(b"^"):
            annotation_name, separator, value = content[1:].partition(
                FIELD_VALUE_SEPARATOR
            )
            if not separator or not annotation_name:
                fail("malformed annotation")
            finish_pending()
            if last_field is None:
                fail("annotation without a field")
            pending = _new_occurrence(
                name_raw=annotation_name,
                value_raw=value,
                raw_line=raw_line,
                line_ending=line_ending,
                start_byte=line_start_byte,
                end_byte=byte_offset,
                line_number=line_number,
            )
            pending_kind = "annotation"
            record_parts.append(raw_line)
            continue

        field_name, separator, value = content.partition(FIELD_VALUE_SEPARATOR)
        if not separator or not field_name:
            fail("malformed field")
        finish_pending()
        if record_start_byte is None:
            record_start_byte = line_start_byte
            record_start_line = line_number
        pending = _new_occurrence(
            name_raw=field_name,
            value_raw=value,
            raw_line=raw_line,
            line_ending=line_ending,
            start_byte=line_start_byte,
            end_byte=byte_offset,
            line_number=line_number,
        )
        pending_kind = "field"
        record_parts.append(raw_line)

    if record_start_byte is not None:
        fail("unterminated record")


def read_file(folder: Path, name: str) -> bytes:
    """Return the bytes of one flat file from the MetaCyc folder."""
    return (Path(folder) / name).read_bytes()


def iter_file_records(source: bytes, name: str) -> Iterator[FlatFileRecord]:
    """Yield the records of flat-file bytes, using ``name`` in the locators."""
    return iter_records(BytesIO(source), source_member_path=name)


def base64_text(raw: bytes) -> str:
    return base64.b64encode(raw).decode("ascii")


def iter_compound_records(source: bytes) -> Iterator[dict[str, object]]:
    """Yield one record per compound of ``compounds.dat``.

    Each record keeps the raw InChI, InChIKey, SMILES and formula fields, the
    parsed formula, the completeness of the declared structure fields (tiers
    I1 to I4, or IX for conflicting or invalid fields) and the heavy-atom count.
    """
    for record in iter_file_records(source, "compounds.dat"):
        fields = record.field_occurrences
        if not fields or fields[0].field_name_raw != b"UNIQUE-ID":
            raise ValueError(
                f"record at {record.record_locator} does not begin with UNIQUE-ID"
            )
        compound_id = fields[0].value_raw.decode("utf-8")
        occurrences = {
            name.decode("ascii"): [
                _raw_occurrence(item) for item in fields if item.field_name_raw == name
            ]
            for name in _STRUCTURE_FIELDS
        }
        raw_values = {
            name: tuple(
                item.value_raw
                for item in fields
                if item.field_name_raw == name.encode("ascii")
            )
            for name in occurrences
        }
        formula = _formula_state(
            [item for item in fields if item.field_name_raw == b"CHEMICAL-FORMULA"]
        )
        nonblank_inchi = _distinct_nonblank(raw_values["INCHI"])
        nonblank_inchi_key = _distinct_nonblank(raw_values["INCHI-KEY"])
        has_structure = any(
            value
            for field_name in ("INCHI", "INCHI-KEY", "SMILES")
            for value in raw_values[field_name]
        )
        tier_issue_codes: list[str] = []
        if len(nonblank_inchi) > 1:
            tier_issue_codes.append("conflicting_inchi")
        if len(nonblank_inchi_key) > 1:
            tier_issue_codes.append("conflicting_inchi_key")
        if formula["status"] == "invalid":
            tier_issue_codes.append("invalid_chemical_formula")

        if tier_issue_codes:
            tier = TIER_IX
        elif (
            len(nonblank_inchi) == 1
            and len(nonblank_inchi_key) == 1
            and formula["status"] == "valid"
        ):
            tier = TIER_I1
        elif has_structure:
            tier = TIER_I2
        elif formula["status"] == "valid":
            tier = TIER_I3
        else:
            tier = TIER_I4

        if formula["status"] == "valid":
            heavy_atom_count = {
                "status": "evaluable",
                "value": sum(
                    int(item["count"])
                    for item in formula["element_counts"]
                    if item["element"] != "H"
                ),
                "reason": None,
            }
        else:
            heavy_atom_count = {
                "status": "not_evaluable",
                "value": None,
                "reason": "formula_absent"
                if formula["status"] == "absent"
                else "formula_invalid",
            }
        payload: dict[str, object] = {
            "schema_version": COMPOUND_RECORD_SCHEMA_VERSION,
            "compound_id": compound_id,
            "source": {
                "source_member_path": record.source_member_path,
                "record_ordinal": record.ordinal,
                "record_locator": record.record_locator,
                "record_bytes_sha256": hashlib.sha256(record.raw_bytes).hexdigest(),
                "field_occurrence_count": len(record.field_occurrences),
            },
            "raw_field_occurrences": occurrences,
            "chemical_formula": formula,
            "chemical_identity_completeness_tier": tier,
            "tier_issue_codes": sorted(tier_issue_codes),
            "heavy_atom_count": heavy_atom_count,
        }
        payload["control_covariate_record_id"] = line_content_id(payload)
        yield payload


def _formula_state(fields: Sequence[FieldOccurrence]) -> dict[str, object]:
    """Parse ``CHEMICAL-FORMULA`` fields of the form ``(C 27)``, one element each."""
    if not fields:
        return {"status": "absent", "element_counts": [], "issues": []}
    element_counts: dict[str, int] = {}
    first_ordinals: dict[str, int] = {}
    issues: list[dict[str, object]] = []
    for item in fields:
        match = _FORMULA_TOKEN_RE.fullmatch(item.value_raw)
        if match is None:
            issues.append(
                {
                    "code": "invalid_formula_token",
                    "field_ordinal": item.ordinal,
                    "value_raw_base64": base64_text(item.value_raw),
                }
            )
            continue
        element = match.group(1).decode("ascii")
        if element not in FORMULA_ELEMENTS:
            issues.append(
                {
                    "code": "formula_element_not_allowed",
                    "field_ordinal": item.ordinal,
                    "element": element,
                    "value_raw_base64": base64_text(item.value_raw),
                }
            )
            continue
        if element in element_counts:
            issues.append(
                {
                    "code": "duplicate_formula_element",
                    "field_ordinal": item.ordinal,
                    "first_field_ordinal": first_ordinals[element],
                    "element": element,
                    "value_raw_base64": base64_text(item.value_raw),
                }
            )
            continue
        element_counts[element] = int(match.group(2))
        first_ordinals[element] = item.ordinal
    if issues:
        return {"status": "invalid", "element_counts": [], "issues": issues}
    return {
        "status": "valid",
        "element_counts": [
            {"element": element, "count": element_counts[element]}
            for element in sorted(element_counts)
        ],
        "issues": [],
    }


def _raw_occurrence(item: FieldOccurrence) -> dict[str, object]:
    return {
        "field_ordinal": item.ordinal,
        "field_locator": item.locator,
        "value_raw_base64": base64_text(item.value_raw),
    }


def _distinct_nonblank(values: Sequence[bytes]) -> tuple[bytes, ...]:
    seen: set[bytes] = set()
    result: list[bytes] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            result.append(value)
    return tuple(result)
