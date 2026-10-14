"""Rarefaction (repeated draws) and depth selection."""

from __future__ import annotations

import numpy as np
import pandas as pd


def rarefy_matrix(counts: np.ndarray, depth: int, rng) -> np.ndarray:
    """Sample without replacement to `depth` reads per sample (rows). Samples below depth must be removed first."""
    out = np.zeros_like(counts, dtype=np.int64)
    for i in range(counts.shape[0]):
        row = counts[i].astype(np.int64)
        nz = np.flatnonzero(row)
        if row.sum() == depth:
            out[i] = row
        else:
            out[i, nz] = rng.multivariate_hypergeometric(row[nz], depth)
    return out


def choose_depth(depths: pd.Series, requested, min_depth: int | None, notes):
    """Return (depth or None, kept sample ids, dropped Series of reasons)."""
    reasons = {}
    d = depths.copy()
    zero = d[d <= 0].index
    for s in zero:
        reasons[s] = "zero reads"
    d = d.drop(zero)
    if min_depth:
        low = d[d < min_depth].index
        for s in low:
            reasons[s] = f"depth {int(d[s])} < --min-depth {min_depth}"
        d = d.drop(low)
    if str(requested).lower() == "none":
        return None, list(d.index), pd.Series(reasons, dtype=object)
    if str(requested).lower() == "auto":
        depth = int(d.min()) if len(d) else None
        if depth is not None and depth < 0.1 * float(d.median()):
            notes.warn("rarefaction", f"Automatic depth ({depth}) is <10% of the median depth ({int(d.median())}); "
                       "one very shallow sample is dictating the depth. Consider --min-depth or --rarefy-depth.")
    else:
        depth = int(requested)
        low = d[d < depth].index
        for s in low:
            reasons[s] = f"depth {int(d[s])} < rarefaction depth {depth}"
        d = d.drop(low)
    return depth, list(d.index), pd.Series(reasons, dtype=object)


def rarefaction_curves(counts: np.ndarray, ids, rng, steps: int = 12, max_depth: int | None = None) -> pd.DataFrame:
    """Observed features at increasing depths (one draw per step), for QC plots."""
    tot = counts.sum(1)
    top = int(max_depth or tot.min())
    depths = sorted(set(np.linspace(max(top // steps, 1), top, steps).astype(int)))
    rows = []
    for d in depths:
        for i, sid in enumerate(ids):
            if tot[i] < d:
                continue
            row = counts[i].astype(np.int64)
            nz = np.flatnonzero(row)
            samp = rng.multivariate_hypergeometric(row[nz], d)
            rows.append(dict(sample_id=sid, depth=d, observed_features=int((samp > 0).sum())))
    return pd.DataFrame(rows)
