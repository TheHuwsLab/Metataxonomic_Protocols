"""Model-based contrast families shared by alpha-diversity and taxon models."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from .models import Fit, contrast


def adjust_p(p, method: str) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method=method)[1]
    return out


def contrast_family(fit: Fit, frame, mode: str = "pairwise", conf: float = 0.95) -> pd.DataFrame:
    """Contrasts of group means (within each time, or overall) and time means (within each group).

    mode='pairwise': all group pairs; mode='reference': each group vs the reference (first level).
    Direction is always first minus second. Each row carries its family label for multiplicity control.
    """
    fx = fit.fixed
    gl, tl = frame.group_levels, frame.time_levels if frame.has_time else []
    rows = []

    def add(family, a, b, ca, cb, **lab):
        r = contrast(fit, ca - cb, conf)
        rows.append(dict(family=family, **lab, level_a=a, level_b=b, **r))

    pairs = list(itertools.combinations(gl, 2)) if mode == "pairwise" else [(g, gl[0]) for g in gl[1:]]
    if tl:
        for t in tl:
            for a, b in pairs:
                add("group within time", a, b, fx.row(a, t), fx.row(b, t), time=t, group="")
        tpairs = list(itertools.combinations(tl, 2)) if mode == "pairwise" else [(t, tl[0]) for t in tl[1:]]
        for g in gl:
            for a, b in tpairs:
                add("time within group", a, b, fx.row(g, a), fx.row(g, b), time="", group=g)
    else:
        for a, b in pairs:
            add("group", a, b, fx.row(a, None), fx.row(b, None), time="", group="")
    df = pd.DataFrame(rows)
    return df
