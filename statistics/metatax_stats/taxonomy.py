"""Taxonomy parsing and collapsing to ranks.

Supports QIIME2/Greengenes/GTDB/UNITE style prefixes (``d__``/``k__``/``p__`` ...),
SILVA ``D_0__`` style, and legacy positional strings. Unresolved labels
(unknown/uncultured/blank ...) never become named taxa: they are collapsed into
lineage-specific ``<ancestor>_unclassified`` bins that are tracked explicitly.
Collapsing is by *lineage*, not by name, so homonyms in different lineages are
never merged (they receive a ``name [parent]`` label instead).
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy import sparse

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
LEVELS = ["feature"] + RANKS
_PREFIX_TO_RANK = {"k": "kingdom", "d": "kingdom", "p": "phylum", "c": "class",
                   "o": "order", "f": "family", "g": "genus", "s": "species"}
_PREFIX = re.compile(r"^([dkpcofgs])__(.*)$", re.I)
_SILVA = re.compile(r"^D_(\d+)__(.*)$")
_BLANK = re.compile(r"^[_\s]*$")
_UNRESOLVED = (
    re.compile(r"^(?:unknown|unclassified|uncultured|unidentified|unassigned|incertae[_ ]sedis|nan|none|na)$", re.I),
    re.compile(r"^(?:unknown|unclassified|uncultured|unidentified|unassigned|metagenome|environmental|no[_ -]?hit|ambiguous)[_ -].*$", re.I),
    re.compile(r"^sp\.?(?:$|[_ -].*)", re.I),
)
UNCLASSIFIED_SUFFIX = "_unclassified"
UNASSIGNED = "Unassigned"


def is_unresolved_value(value) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return True
    text = str(value).strip()
    if _BLANK.match(text):
        return True
    return any(p.match(text) for p in _UNRESOLVED)


def parse_lineage(lineage) -> dict:
    """Parse one lineage string into {rank: value or None}."""
    parsed = {r: None for r in RANKS}
    if lineage is None or (isinstance(lineage, float) and np.isnan(lineage)):
        return parsed
    for i, raw in enumerate(str(lineage).split(";")):
        part = raw.strip()
        rank, value = None, part
        m = _PREFIX.match(part)
        s = _SILVA.match(part)
        if m:
            rank, value = _PREFIX_TO_RANK[m.group(1).lower()], m.group(2).strip()
        elif s:
            idx = int(s.group(1))
            rank, value = (RANKS[idx] if idx < len(RANKS) else None), s.group(2).strip()
        elif i < len(RANKS):
            rank = RANKS[i]
        if rank is not None and not is_unresolved_value(value):
            parsed[rank] = value
    return parsed


def looks_like_lineage(labels) -> bool:
    labels = pd.Index(labels).astype(str)
    if len(labels) == 0:
        return False
    return bool(labels.str.contains(r";|(?:^|;)\s*[dkpcofgs]__|D_\d+__", regex=True).mean() >= 0.5)


def lineage_frame(lineages: pd.Series) -> pd.DataFrame:
    rows = {f: parse_lineage(l) for f, l in lineages.items()}
    return pd.DataFrame.from_dict(rows, orient="index", columns=RANKS)


def available_levels(lineages: pd.Series | None) -> list[str]:
    if lineages is None:
        return ["feature"]
    frame = lineage_frame(lineages)
    return ["feature"] + [r for r in RANKS if frame[r].notna().any()]


class LevelMap:
    """Mapping from input features to the labels of one taxonomic level."""

    def __init__(self, level: str, table: pd.DataFrame):
        self.level = level
        self.table = table  # index=feature; columns label,status,deepest_ancestor,reason
        self.labels = pd.Index(pd.unique(table["label"]))
        col = self.labels.get_indexer(table["label"])
        self.matrix = sparse.csr_matrix(
            (np.ones(len(table)), (col, np.arange(len(table)))), shape=(len(self.labels), len(table)))
        st = table.groupby("label")["status"].agg(lambda x: "unresolved" if (x == "unresolved").any() else "resolved")
        self.status = st.reindex(self.labels)

    def aggregate(self, counts: np.ndarray) -> np.ndarray:
        """Sum features (rows) into labels. counts: features x samples."""
        return np.asarray(self.matrix @ counts)

    @property
    def unresolved_mask(self) -> np.ndarray:
        return self.status.values == "unresolved"


def build_level_map(features: pd.Index, lineages: pd.Series | None, level: str) -> LevelMap:
    level = level.lower()
    features = pd.Index(features)
    if level == "feature":
        t = pd.DataFrame({"label": features.astype(str), "status": "resolved",
                          "deepest_ancestor": "", "reason": ""}, index=features)
        return LevelMap(level, t)
    if level not in RANKS:
        raise ValueError(f"Unknown taxonomic level: {level}")
    if lineages is None:
        raise ValueError("No taxonomy available; only the 'feature' level can be analysed.")
    pos = RANKS.index(level)
    frame = lineage_frame(lineages.reindex(features))
    keys, names, statuses, deepest, reasons = [], [], [], [], []
    for f in features:
        vals = [None if pd.isna(frame.at[f, r]) else frame.at[f, r] for r in RANKS[: pos + 1]]
        anc_idx = next((i for i in range(pos - 1, -1, -1) if vals[i] is not None), None)
        anc = vals[anc_idx] if anc_idx is not None else None
        if vals[-1] is not None:
            keys.append(tuple(vals)); names.append(vals[-1]); statuses.append("resolved")
            deepest.append(vals[-1]); reasons.append("")
        else:
            keys.append(("unresolved",) + tuple(vals[: anc_idx + 1] if anc_idx is not None else ()))
            names.append(f"{anc}{UNCLASSIFIED_SUFFIX}" if anc else UNASSIGNED)
            statuses.append("unresolved"); deepest.append(anc or "")
            reasons.append(f"{level} unresolved or missing")

    # Distinct lineages that share a name stay distinct: label them "name [parent]".
    distinct: dict[str, list] = {}
    for k, n in zip(keys, names):
        distinct.setdefault(n, [])
        if k not in distinct[n]:
            distinct[n].append(k)
    label_of = {}
    for n, ks in distinct.items():
        for i, k in enumerate(ks):
            if len(ks) == 1:
                label_of[k] = n
                continue
            anc = [v for v in (k[1:] if k[0] == "unresolved" else k[:-1]) if v is not None]
            label_of[k] = f"{n} [{anc[-1]}]" if anc else f"{n} [#{i + 1}]"
    labels = [label_of[k] for k in keys]
    owner: dict[str, tuple] = {}
    for i, (lab, k) in enumerate(zip(labels, keys)):
        n_try = 2
        while lab in owner and owner[lab] != k:
            lab = f"{labels[i]} #{n_try}"; n_try += 1
        owner[lab] = k
        labels[i] = lab
    t = pd.DataFrame({"label": labels, "status": statuses, "deepest_ancestor": deepest, "reason": reasons},
                     index=features)
    return LevelMap(level, t)


def is_unclassified_label(label: str) -> bool:
    text = str(label).strip()
    return text == UNASSIGNED or text.endswith(UNCLASSIFIED_SUFFIX) or is_unresolved_value(text)
