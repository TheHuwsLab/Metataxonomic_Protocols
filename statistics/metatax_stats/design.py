"""Study-design specification, fixed-effect design matrices, estimability and degrees of freedom.

The design is generic:

* ``group``   primary factor (treatment, diet, site, disease state, ...)
* ``time``    optional second categorical factor (time point, phase, batch, ...)
* ``subject`` optional repeated-measures unit (vessel, animal, patient, ...)
* ``cluster`` optional unit of allocation containing subjects (pen, cage, hospital, ...)
* ``covariates`` optional extra numeric/categorical adjustment terms

Fixed effects use sum (effect) coding so that main effects are averages over the
other factor rather than simple effects at a reference level.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import lgamma, log

import numpy as np
import pandas as pd
from scipy import linalg

DESIGN_COLS = ["group", "time", "subject", "cluster", "upper", "lower"]


@dataclass
class DesignSpec:
    group_col: str
    time_col: str | None = None
    subject_col: str | None = None
    cluster_col: str | None = None
    covariates: tuple = ()
    group_order: list | None = None
    time_order: list | None = None
    reference_group: str | None = None

    def required(self) -> list[str]:
        cols = [self.group_col, self.time_col, self.subject_col, self.cluster_col, *self.covariates]
        return [c for c in cols if c]

    def describe(self) -> str:
        parts = [f"outcome ~ {self.group_col}"]
        if self.time_col:
            parts.append(f"* {self.time_col}")
        parts.extend(f"+ {c}" for c in self.covariates)
        if self.cluster_col and self.subject_col:
            parts.append(f"+ (1 | {self.cluster_col}) + (1 | {self.cluster_col}:{self.subject_col})")
        elif self.cluster_col or self.subject_col:
            parts.append(f"+ (1 | {self.cluster_col or self.subject_col})")
        return " ".join(parts)


@dataclass
class Frame:
    """Analysis-ready metadata with standardised column names."""
    df: pd.DataFrame
    spec: DesignSpec
    group_levels: list
    time_levels: list
    cov_kinds: dict = field(default_factory=dict)  # name -> 'numeric' | 'categorical'

    @property
    def has_time(self) -> bool:
        return self.spec.time_col is not None and len(self.time_levels) > 1

    @property
    def has_units(self) -> bool:
        return self.df["upper"].notna().any()

    def subset(self, ids) -> "Frame":
        df = self.df.loc[list(ids)].copy()
        gl = [g for g in self.group_levels if g in set(df["group"])]
        tl = [t for t in self.time_levels if t in set(df["time"])] if "time" in df else []
        return Frame(df, self.spec, gl, tl, self.cov_kinds)


def _blank_to_nan(s: pd.Series) -> pd.Series:
    s = s.astype("string").str.strip()
    return s.mask(s.isin(["", "nan", "NaN", "None", "NA"]))


def build_frame(meta: pd.DataFrame, spec: DesignSpec, notes):
    """Return (Frame, excluded[(sample, reason)])."""
    missing_cols = [c for c in spec.required() if c not in meta.columns]
    if missing_cols:
        raise SystemExit(f"Metadata columns not found: {missing_cols}. Available: {list(meta.columns)}")
    df = pd.DataFrame(index=meta.index)
    df["group"] = _blank_to_nan(meta[spec.group_col])
    df["time"] = _blank_to_nan(meta[spec.time_col]) if spec.time_col else pd.Series("all", index=meta.index, dtype="string")
    subj = _blank_to_nan(meta[spec.subject_col]) if spec.subject_col else pd.Series(pd.NA, index=meta.index, dtype="string")
    clus = _blank_to_nan(meta[spec.cluster_col]) if spec.cluster_col else pd.Series(pd.NA, index=meta.index, dtype="string")
    required = ["group"] + (["time"] if spec.time_col else [])
    covs = {}
    kinds = {}
    for c in spec.covariates:
        raw = _blank_to_nan(meta[c])
        num = pd.to_numeric(raw, errors="coerce")
        numeric = bool(((~num.isna()) | raw.isna()).all()) and num.notna().sum() > 0
        covs[c] = num.astype(float) if numeric else raw
        kinds[c] = "numeric" if numeric else "categorical"
    complete = pd.Series(True, index=meta.index)
    reasons = {}
    check = {"group": df["group"], "time": df["time"]}
    if spec.subject_col:
        check["subject"] = subj
    if spec.cluster_col:
        check["cluster"] = clus
    check.update({f"covariate {c}": v for c, v in covs.items()})
    for name, series in check.items():
        if name == "time" and not spec.time_col:
            continue
        bad = series.isna()
        for s in meta.index[bad & complete]:
            reasons[s] = f"missing {name}"
        complete &= ~bad
    excluded = [(s, reasons[s]) for s in meta.index[~complete]]
    df = df.loc[complete].copy()
    df["subject"] = subj.loc[complete].astype(object)
    df["cluster"] = clus.loc[complete].astype(object)
    for c in ("group", "time"):
        df[c] = df[c].astype(str)
    if spec.cluster_col and spec.subject_col:
        df["subject"] = df["cluster"].astype(str) + "::" + df["subject"].astype(str)
        df["upper"], df["lower"] = df["cluster"].astype(str), df["subject"]
    elif spec.cluster_col:
        df["upper"], df["lower"] = df["cluster"].astype(str), None
    elif spec.subject_col:
        df["upper"], df["lower"] = df["subject"].astype(str), None
    else:
        df["upper"], df["lower"] = None, None
    for c, v in covs.items():
        df[c] = v.loc[complete]
    if df.empty:
        raise SystemExit("No samples remain after removing rows with missing design metadata.")

    def levels(col, order, name):
        present = list(pd.unique(df[col]))
        if order:
            order = [str(o) for o in order]
            unknown = [o for o in order if o not in present]
            if unknown:
                notes.warn("design", f"{name} order lists levels not present: {unknown}")
            present = [o for o in order if o in present] + [p for p in present if p not in order]
        return present

    gl = levels("group", spec.group_order, "group")
    if spec.reference_group is not None:
        ref = str(spec.reference_group)
        if ref not in gl:
            raise SystemExit(f"--reference-group '{ref}' not found among {gl}")
        gl = [ref] + [g for g in gl if g != ref]
    tl = levels("time", spec.time_order, "time") if spec.time_col else ["all"]
    if len(gl) < 2:
        raise SystemExit(f"Group column '{spec.group_col}' needs at least 2 levels after filtering; found {gl}.")
    if spec.time_col and len(tl) < 2:
        notes.warn("design", f"Time column '{spec.time_col}' has a single level; analysed without a time factor.")
    return Frame(df, spec, gl, tl, kinds), excluded


# -----------------------------------------------------------------------------
# Fixed-effect design matrix (sum coding)
# -----------------------------------------------------------------------------

class FixedDesign:
    """Sum-coded design for group (* time) + covariates."""

    def __init__(self, frame: Frame):
        self.frame = frame
        df = frame.df
        self.gl = frame.group_levels
        self.tl = frame.time_levels if frame.has_time else []
        self.cov_means = {c: float(df[c].mean()) for c, k in frame.cov_kinds.items() if k == "numeric"}
        self.cov_levels = {c: sorted(df[c].astype(str).unique()) for c, k in frame.cov_kinds.items() if k == "categorical"}
        cols, terms = [("int",)], {}
        terms["group"] = [("g", lv) for lv in self.gl[:-1]]
        cols += terms["group"]
        if self.tl:
            terms["time"] = [("t", lv) for lv in self.tl[:-1]]
            cols += terms["time"]
            terms["group:time"] = [("gt", g[1], t[1]) for g in terms["group"] for t in terms["time"]]
            cols += terms["group:time"]
        for c, kind in frame.cov_kinds.items():
            if kind == "numeric":
                terms[c] = [("cov", c)]
            else:
                terms[c] = [("cat", c, lv) for lv in self.cov_levels[c][:-1]]
            cols += terms[c]
        self.specs = cols
        self.names = [self._name(s) for s in cols]
        self.terms = {t: [self._name(s) for s in ss] for t, ss in terms.items()}
        self.X = pd.DataFrame(self.matrix(df), index=df.index, columns=self.names)

    @staticmethod
    def _name(s) -> str:
        return {"int": "Intercept"}.get(s[0]) or ("G[%s]" % s[1] if s[0] == "g" else "T[%s]" % s[1] if s[0] == "t"
                else "G[%s]:T[%s]" % (s[1], s[2]) if s[0] == "gt" else "C[%s]" % s[1] if s[0] == "cov"
                else "K[%s=%s]" % (s[1], s[2]))

    def _sum(self, values, lv, levels):
        v = np.asarray(values, dtype=object)
        return (v == lv).astype(float) - (v == levels[-1]).astype(float)

    def matrix(self, df: pd.DataFrame) -> np.ndarray:
        out = np.empty((len(df), len(self.specs)))
        g, t = df["group"].to_numpy(object), df["time"].to_numpy(object)
        for j, s in enumerate(self.specs):
            if s[0] == "int":
                out[:, j] = 1.0
            elif s[0] == "g":
                out[:, j] = self._sum(g, s[1], self.gl)
            elif s[0] == "t":
                out[:, j] = self._sum(t, s[1], self.tl)
            elif s[0] == "gt":
                out[:, j] = self._sum(g, s[1], self.gl) * self._sum(t, s[2], self.tl)
            elif s[0] == "cov":
                out[:, j] = df[s[1]].to_numpy(float) - self.cov_means[s[1]]
            else:
                out[:, j] = self._sum(df[s[1]].astype(str).to_numpy(object), s[2], self.cov_levels[s[1]])
        return out

    def row(self, group=None, time=None) -> np.ndarray:
        """Contrast row for a cell mean; None averages over that factor (covariates at their mean)."""
        v = np.zeros(len(self.specs))
        for j, s in enumerate(self.specs):
            gval = 0.0 if group is None else float(group == s[1]) - float(group == self.gl[-1]) if s[0] in ("g", "gt") else 0.0
            if s[0] == "int":
                v[j] = 1.0
            elif s[0] == "g":
                v[j] = gval
            elif s[0] == "t":
                v[j] = 0.0 if time is None else float(time == s[1]) - float(time == self.tl[-1])
            elif s[0] == "gt":
                tval = 0.0 if time is None else float(time == s[2]) - float(time == self.tl[-1])
                v[j] = gval * tval
        return v


# -----------------------------------------------------------------------------
# Estimability and containment degrees of freedom
# -----------------------------------------------------------------------------

@dataclass
class RankInfo:
    kept: list
    rank: int
    rowspace: np.ndarray  # orthonormal basis (p x r)
    dropped_terms: list

    def estimable(self, c: np.ndarray, tol: float = 1e-8) -> bool:
        c = np.asarray(c, float)
        resid = c - self.rowspace @ (self.rowspace.T @ c)
        return bool(np.linalg.norm(resid) <= tol * (1.0 + np.linalg.norm(c)))


def analyse_rank(fixed: FixedDesign) -> RankInfo:
    X = fixed.X.to_numpy(float)
    _, R, piv = linalg.qr(X, pivoting=True, mode="economic")
    d = np.abs(np.diag(R))
    tol = max(X.shape) * np.finfo(float).eps * (d[0] if len(d) else 1.0) * 10
    rank = int((d > tol).sum())
    kept = sorted(piv[:rank].tolist())
    _, s, Vt = np.linalg.svd(X, full_matrices=False)
    r = int((s > max(X.shape) * np.finfo(float).eps * s[0] * 10).sum())
    names = fixed.names
    dropped = sorted({t for t, cols in fixed.terms.items() for c in cols if names.index(c) not in kept})
    return RankInfo(kept, rank, Vt[:r].T, dropped)


def _const_within(x: np.ndarray, codes: np.ndarray) -> bool:
    means = np.bincount(codes, x) / np.bincount(codes)
    return bool(np.abs(x - means[codes]).max() < 1e-9)


def _center(M: np.ndarray, codes: np.ndarray) -> np.ndarray:
    cnt = np.bincount(codes)
    out = M.copy()
    for j in range(M.shape[1]):
        out[:, j] = M[:, j] - (np.bincount(codes, M[:, j]) / cnt)[codes]
    return out


def _rank(M: np.ndarray) -> int:
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > max(M.shape) * np.finfo(float).eps * (s[0] if len(s) else 1) * 10).sum())


def containment_ddf(fixed: FixedDesign, kept: list) -> dict:
    """Containment denominator df for every fixed-effect column.

    Strata: columns constant within the upper unit (between-unit), constant within the
    lower unit (between-subject within cluster), or varying within the lowest unit.
    Between-unit: n_U - rank(between cols); within-unit: n_lower - n_U - rank(centred cols);
    residual: n - n_lowest - rank(centred cols). Floored at 1.
    """
    df = fixed.frame.df
    X = fixed.X.to_numpy(float)
    n = len(X)
    if df["upper"].isna().all():
        r = _rank(X[:, kept])
        return {j: max(n - r, 1) for j in kept}
    u = pd.factorize(df["upper"])[0]
    nU = u.max() + 1
    has_lower = df["lower"].notna().any()
    s = pd.factorize(df["lower"])[0] if has_lower else None
    nS = s.max() + 1 if has_lower else None
    strata = {}
    for j in kept:
        x = X[:, j]
        if _const_within(x, u):
            strata[j] = 0
        elif has_lower and _const_within(x, s):
            strata[j] = 1
        else:
            strata[j] = 2 if has_lower else 1
    out = {}
    s0 = [j for j in kept if strata[j] == 0]
    d0 = max(nU - _rank(X[:, s0]), 1)
    s1 = [j for j in kept if strata[j] == 1]
    if has_lower:
        d1 = max(nS - nU - _rank(_center(X[:, s1], u)), 1) if s1 else 1
        s2 = [j for j in kept if strata[j] == 2]
        d2 = max(n - nS - _rank(_center(X[:, s2], s)), 1) if s2 else 1
        dd = {0: d0, 1: d1, 2: d2}
    else:
        dd = {0: d0, 1: max(n - nU - _rank(_center(X[:, s1], u)), 1) if s1 else 1}
    return {j: float(dd[strata[j]]) for j in kept}


def min_attainable_p(counts_per_group: list[int]) -> tuple[float, float]:
    """(number of distinct label arrangements, smallest attainable permutation p)."""
    n = sum(counts_per_group)
    lg = lgamma(n + 1) - sum(lgamma(c + 1) for c in counts_per_group)
    arrangements = float(np.exp(lg)) if lg < 700 else float("inf")
    return arrangements, (1.0 / arrangements if np.isfinite(arrangements) else 0.0)


# -----------------------------------------------------------------------------
# Design validation report
# -----------------------------------------------------------------------------

def validate_design(frame: Frame, fixed: FixedDesign, rank: RankInfo, excluded, notes):
    """Return (checks DataFrame, cell-count DataFrame, repeated-measures audit DataFrame)."""
    df = frame.df
    rows = []

    def add(check, status, detail):
        rows.append({"check": check, "status": status, "detail": detail})
        if status in ("WARN", "FAIL"):
            notes.warn("design", f"{check}: {detail}")

    add("design_formula", "INFO", frame.spec.describe())
    add("samples_retained", "INFO", f"{len(df)}")
    if excluded:
        add("samples_excluded", "WARN", f"{len(excluded)} excluded: " + "; ".join(f"{s} ({w})" for s, w in excluded[:20]))
    gcount = df.groupby("group").size().reindex(frame.group_levels)
    add("group_levels", "INFO", ", ".join(f"{g} (n={n})" for g, n in gcount.items()))
    if frame.has_time:
        tcount = df.groupby("time").size().reindex(frame.time_levels)
        add("time_levels", "INFO", ", ".join(f"{t} (n={n})" for t, n in tcount.items()))

    # cell counts
    key = ["group", "time"] if frame.has_time else ["group"]
    cells = df.groupby(key).agg(n_samples=("group", "size"),
                                n_units=("upper", lambda x: x.nunique() if x.notna().any() else np.nan)).reset_index()
    if frame.has_time:
        full = pd.MultiIndex.from_product([frame.group_levels, frame.time_levels], names=key).to_frame(index=False)
        cells = full.merge(cells, on=key, how="left").fillna({"n_samples": 0})
        empty = cells[cells["n_samples"] == 0]
        if len(empty):
            add("empty_cells", "WARN", f"{len(empty)} group x time cell(s) have no samples; related effects are not estimable.")
    low = gcount[gcount < 3]
    if len(low):
        add("small_groups", "WARN", f"Groups with <3 samples: {low.to_dict()}")

    audit = pd.DataFrame()
    if frame.has_units:
        units_per_group = df.groupby("group")["upper"].nunique().reindex(frame.group_levels)
        add("units_per_group", "INFO", ", ".join(f"{g}: {n}" for g, n in units_per_group.items()))
        n_assignment = units_per_group.astype(int)
        unit_name = frame.spec.cluster_col or frame.spec.subject_col or "sample"
        add("independent_treatment_units",
            "WARN" if int(n_assignment.min()) < 5 else "INFO",
            f"Independent treatment-assignment units ({unit_name}) per group: "
            + ", ".join(f"{g}: {n}" for g, n in n_assignment.items())
            + ("; small numbers of independent units can make treatment inference unstable and limit permutation "
               "resolution. The warning threshold is a screening heuristic, not a validity or adequacy cutoff."
               if int(n_assignment.min()) < 5 else "."))
        const = (df.groupby("upper")["group"].nunique() == 1).all()
        add("group_constant_within_unit", "PASS" if const else "WARN",
            "Group is constant within each unit (between-unit factor)." if const else
            "Group varies within at least one unit: restricted-permutation beta tests are skipped (alpha/taxon models still apply).")
        if frame.spec.cluster_col and frame.spec.subject_col:
            lower_counts = df.groupby("lower").size()
            if not (lower_counts > 1).any():
                add("lower_random_effect", "INFO",
                    "Each nested subject has one observation; a separate subject-level random-intercept variance "
                    "cannot be identified from residual variance and is omitted from the REML model.")
        if const:
            arrangements, minp = min_attainable_p(units_per_group.astype(int).tolist())
            status = "PASS" if minp <= 0.05 else "WARN"
            add("min_attainable_p_group", status,
                f"{arrangements:.4g} distinct unit-label arrangements; smallest attainable permutation p = {minp:.3g}"
                + ("" if status == "PASS" else " (the between-unit permutation test can never reach p<0.05)"))
        # QC completeness is evaluated at the upper independent unit used for
        # treatment-level restricted permutations, including nested subjects.
        rep = df["upper"]
        if frame.has_time:
            rows_audit = []
            tl = set(frame.time_levels)
            for u, sub in df.assign(_rep=rep).groupby("_rep"):
                obs = list(pd.unique(sub["time"]))
                complete = len(obs) == len(tl) and set(obs) == tl
                rows_audit.append({"repeated_unit": u, "group": "; ".join(sub["group"].unique()),
                                   "n_samples": len(sub), "time_levels_present": "; ".join(obs),
                                   "complete": complete})
            audit = pd.DataFrame(rows_audit)
            n_inc = int((~audit["complete"]).sum())
            add("repeated_measures_complete", "PASS" if n_inc == 0 else "WARN",
            "All upper independent units have every time level." if n_inc == 0 else
            f"{n_inc} of {len(audit)} upper independent units are incomplete (restricted-permutation beta tests use "
            "complete upper units only; mixed models use all data).")
    else:
        add("units", "INFO", "No subject/cluster column: all samples treated as independent.")
        if int(gcount.min()) < 5:
            add("independent_treatment_units", "WARN",
                "Independent treatment-assignment units (samples) per group: "
                + ", ".join(f"{g}: {int(n)}" for g, n in gcount.items())
                + "; small numbers of independent units can make treatment inference unstable and limit permutation "
                "resolution. The warning threshold is a screening heuristic, not a validity or adequacy cutoff.")
    if rank.dropped_terms:
        add("estimability", "WARN", f"Non-estimable/aliased terms (columns dropped): {rank.dropped_terms}")
    else:
        add("estimability", "PASS", f"Design matrix has full column rank ({rank.rank}).")
    return pd.DataFrame(rows), cells, audit
