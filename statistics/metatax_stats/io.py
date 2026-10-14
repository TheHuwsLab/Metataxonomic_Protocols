"""Input handling: feature tables, taxonomy files, sample metadata."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .taxonomy import looks_like_lineage

ANNOTATION_COLUMNS = {"taxonomy", "taxon", "consensus lineage", "lineage", "classification", "confidence",
                      "feature id", "featureid", "otu id", "#otu id", "sequence", "seq"}
TAXONOMY_COLUMN_NAMES = ["taxonomy", "taxon", "consensus lineage", "lineage", "classification"]


class InputError(SystemExit):
    pass


def _sep_for(path: Path):
    suf = path.suffix.lower()
    if suf == ".csv":
        return ","
    if suf in {".tsv", ".txt", ".tab"}:
        return "\t"
    return None


def read_delimited(path: Path) -> pd.DataFrame:
    """Read CSV/TSV as strings. Skips the QIIME2/BIOM '# Constructed from biom file' banner."""
    path = Path(path)
    with open(path, "r", encoding="utf-8-sig") as fh:
        first = fh.readline()
    skip = 1 if first.startswith("# Constructed from biom") else 0
    sep = _sep_for(path)
    kw = dict(dtype=str, skiprows=skip, na_values=["", "NA", "N/A", "nan", "NaN"], keep_default_na=False)
    df = pd.read_csv(path, sep=sep, engine="python" if sep is None else "c", **kw)
    df = df.dropna(how="all")
    if df.empty:
        raise InputError(f"File is empty: {path}")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def read_metadata(path: Path, sample_col: str | None = None) -> pd.DataFrame:
    df = read_delimited(path)
    sid = sample_col or df.columns[0]
    if sid not in df.columns:
        raise InputError(f"Sample-ID column '{sid}' not found in metadata. Columns: {list(df.columns)}")
    df = df.rename(columns={sid: "sample_id"})
    df = df[~df["sample_id"].astype(str).str.startswith("#q2:")]  # QIIME2 types row
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    dup = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
    if dup:
        raise InputError(f"Duplicate sample IDs in metadata: {dup[:10]}")
    for c in df.columns:
        if c != "sample_id":
            df[c] = df[c].astype("string").str.strip()
    return df.set_index("sample_id")


def read_taxonomy_file(path: Path) -> pd.Series:
    df = read_delimited(path)
    fcol = df.columns[0]
    tcol = next((c for c in df.columns if c.lower() in TAXONOMY_COLUMN_NAMES), None)
    if tcol is None:
        if df.shape[1] < 2:
            raise InputError("Taxonomy file needs a feature-ID column and a taxonomy column.")
        tcol = df.columns[1]
    s = df.set_index(df[fcol].astype(str).str.strip())[tcol]
    if s.index.duplicated().any():
        raise InputError("Duplicate feature IDs in taxonomy file.")
    return s.astype(str)


def read_feature_table(path: Path, metadata_ids: pd.Index, taxonomy_path: Path | None, notes):
    """Return (counts[features x samples], lineage Series or None).

    Handles samples-as-rows or samples-as-columns, QIIME2 banner lines, extra
    taxonomy columns, and a separate taxonomy file for ASV/OTU tables.
    """
    raw = read_delimited(path)
    raw = raw.set_index(raw.columns[0])
    raw.index = raw.index.astype(str).str.strip()
    ids = pd.Index(metadata_ids.astype(str))
    col_hits = len(raw.columns.intersection(ids))
    row_hits = len(raw.index.intersection(ids))
    if col_hits == 0 and row_hits == 0:
        raise InputError("No sample IDs shared between the abundance table and metadata (checked both "
                         "orientations). Check sample IDs.")
    if row_hits > col_hits:
        raw = raw.T
        notes.info("input", "Abundance table had samples as rows; transposed to features x samples.")
    annotations = {}
    count_cols = []
    for c in raw.columns:
        if str(c).strip().lower() in ANNOTATION_COLUMNS:
            annotations[c] = raw[c]
        else:
            count_cols.append(c)
    counts_raw = raw[count_cols]
    numeric = counts_raw.apply(pd.to_numeric, errors="coerce")
    bad = numeric.isna() & counts_raw.notna()
    # Columns that are entirely non-numeric are annotations (e.g. a lineage column); else it's an error.
    for c in list(counts_raw.columns):
        if numeric[c].isna().all() and counts_raw[c].notna().any() and c not in ids:
            annotations[c] = raw[c]
            counts_raw = counts_raw.drop(columns=c)
            numeric = numeric.drop(columns=c)
            bad = bad.drop(columns=c)
    if bad.to_numpy().any():
        ex = [(str(counts_raw.index[i]), str(counts_raw.columns[j]), str(counts_raw.iat[i, j]))
              for i, j in zip(*np.where(bad.to_numpy()))][:5]
        raise InputError(f"Non-numeric abundance values found; refusing to coerce them to zero. Examples: {ex}")
    counts = numeric.fillna(0.0).astype(float)
    if counts.columns.duplicated().any():
        raise InputError(f"Duplicate sample columns in abundance table: {counts.columns[counts.columns.duplicated()].tolist()[:10]}")
    if counts.index.duplicated().any():
        notes.warn("input", "Duplicate feature rows were summed.")
        counts = counts.groupby(level=0).sum()

    lineages = None
    if taxonomy_path is not None:
        tax = read_taxonomy_file(Path(taxonomy_path))
        lineages = tax.reindex(counts.index)
        miss = int(lineages.isna().sum())
        if miss:
            notes.warn("taxonomy", f"{miss} of {len(counts)} features are absent from the taxonomy file and are treated as unassigned.")
        lineages = lineages.fillna("")
    else:
        for name in TAXONOMY_COLUMN_NAMES:
            col = next((c for c in annotations if str(c).lower() == name), None)
            if col is not None:
                lineages = annotations[col].reindex(counts.index).fillna("").astype(str)
                notes.info("taxonomy", f"Using the table's '{col}' column as lineage.")
                break
        if lineages is None and looks_like_lineage(counts.index):
            lineages = pd.Series(counts.index, index=counts.index).astype(str)
    if lineages is None:
        notes.warn("taxonomy", "No taxonomy found; only the 'feature' level can be analysed. Supply --taxonomy "
                               "for ASV/OTU tables.")
    return counts, lineages


def validate_counts(counts: pd.DataFrame) -> None:
    arr = counts.to_numpy(float)
    if not np.isfinite(arr).all():
        raise InputError("Abundance table contains non-finite values.")
    if (arr < 0).any():
        raise InputError(f"Abundance table contains {(arr < 0).sum()} negative values.")
    frac = np.abs(arr - np.rint(arr)) > 1e-8
    if frac.any():
        raise InputError("Abundance table contains non-integer values; this toolkit needs raw integer counts "
                         f"(rarefaction/CLR/Chao1). Example values: {arr[frac][:5].tolist()}")
