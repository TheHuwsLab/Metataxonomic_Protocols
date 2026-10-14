"""Shared helpers: warning/notes registry, file hashing, naming."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import pandas as pd


class Notes:
    """Collects every warning/caveat so none is lost in stderr; written to warnings.tsv."""

    def __init__(self, verbose: bool = True):
        self.rows: list[dict] = []
        self.verbose = verbose

    def add(self, severity: str, scope: str, message: str) -> None:
        self.rows.append({"severity": severity, "scope": scope, "message": message})
        if self.verbose and severity in {"WARNING", "ERROR"}:
            print(f"{severity} [{scope}]: {message}", file=sys.stderr)

    def warn(self, scope: str, message: str) -> None:
        self.add("WARNING", scope, message)

    def info(self, scope: str, message: str) -> None:
        self.add("INFO", scope, message)

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows, columns=["severity", "scope", "message"]).drop_duplicates()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_name(value) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "x"


def fmt_p(p) -> str:
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "NA"
    if p != p:
        return "NA"
    return "<0.0001" if p < 1e-4 else f"{p:.4f}"
