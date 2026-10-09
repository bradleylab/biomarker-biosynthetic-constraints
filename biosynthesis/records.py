"""Canonical JSON and the SHA-256 identifiers of network records."""

from __future__ import annotations

import hashlib
import json
from typing import Any


def canonical_bytes(value: Any) -> bytes:
    """Serialize JSON with sorted keys, no spaces and UTF-8 text."""
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode()


def content_id(value: Any) -> str:
    """Identify a JSON value by the SHA-256 of its canonical bytes."""
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def canonical_line(value: Any) -> bytes:
    """Serialize JSON as ``canonical_bytes`` followed by one newline."""
    return canonical_bytes(value) + b"\n"


def line_content_id(value: Any) -> str:
    """Identify a JSON value by the SHA-256 of its canonical line (with the newline)."""
    return "sha256:" + hashlib.sha256(canonical_line(value)).hexdigest()
