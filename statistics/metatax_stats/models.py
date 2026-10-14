"""Gaussian linear (mixed) models fitted by REML with explicit estimability and df handling.

Supported random-effect structures (all random intercepts):

* none                       -> OLS
* upper unit                 -> (1 | subject)  or  (1 | cluster)
* upper + repeatedly observed nested lower unit
                             -> (1 | cluster) + (1 | cluster:subject)
* upper + singleton lower unit
                             -> (1 | cluster); lower variance is confounded with residual variance

REML is profiled over the variance ratios, so boundary solutions (a variance
component equal to zero) are handled exactly rather than failing with a singular
matrix. Denominator df are containment df (see design.containment_ddf); they are
approximate for unbalanced designs and are reported as such.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize
from scipy.stats import f as f_dist, levene, shapiro, t as t_dist

from .design import FixedDesign, RankInfo, analyse_rank, containment_ddf


@dataclass
class Fit:
    beta: np.ndarray
    cov: np.ndarray
    fixed: FixedDesign
    rank: RankInfo
    ddf: dict
    sigma2: float
    variance: dict
    model: str
    note: str
    resid_marginal: np.ndarray
    resid_conditional: np.ndarray
    fitted: np.ndarray
    n: int
    boundary: bool = False
    ok: bool = True
    ddf_note: str = ""

    @property
    def icc(self) -> float:
        re = sum(v for k, v in self.variance.items() if k != "residual")
        tot = re + self.variance.get("residual", np.nan)
        return float(re / tot) if tot and np.isfinite(tot) and tot > 0 else np.nan


def _indicator(codes: np.ndarray, n_levels: int) -> np.ndarray:
    Z = np.zeros((len(codes), n_levels))
    Z[np.arange(len(codes)), codes] = 1.0
    return Z


class _Reml:
    def __init__(self, y, X, Zs: list[np.ndarray]):
        self.y, self.X, self.n, self.p = y, X, len(y), X.shape[1]
        self.Z = np.hstack(Zs) if Zs else np.zeros((len(y), 0))
        self.sizes = [z.shape[1] for z in Zs]
        self.XtX, self.Xty, self.yty = X.T @ X, X.T @ y, float(y @ y)
        self.ZtX, self.Zty, self.ZtZ = self.Z.T @ X, self.Z.T @ y, self.Z.T @ self.Z
        self.diag_only = len(Zs) <= 1

    def _scale(self, lam):
        return np.concatenate([np.full(s, np.sqrt(l)) for s, l in zip(self.sizes, lam)]) if self.sizes else np.zeros(0)

    def pieces(self, lam):
        d = self._scale(lam)
        if len(d) == 0:
            return self.XtX, self.Xty, self.yty, 0.0, d, None
        ZtX, Zty = d[:, None] * self.ZtX, d * self.Zty
        if self.diag_only:
            kd = 1.0 + d * d * np.diag(self.ZtZ)
            logdet = float(np.log(kd).sum())
            solve = lambda M: M / (kd[:, None] if M.ndim == 2 else kd)
        else:
            K = np.eye(len(d)) + (d[:, None] * self.ZtZ) * d[None, :]
            sign, logdet = np.linalg.slogdet(K)
            if sign <= 0:
                return None
            solve = lambda M: np.linalg.solve(K, M)
        XVX = self.XtX - ZtX.T @ solve(ZtX)
        XVy = self.Xty - ZtX.T @ solve(Zty)
        yVy = self.yty - float(Zty @ solve(Zty))
        return XVX, XVy, yVy, logdet, d, solve

    def objective(self, lam):
        pc = self.pieces(lam)
        if pc is None:
            return np.inf
        XVX, XVy, yVy, logdet, _, _ = pc
        s, ld_a = np.linalg.slogdet(XVX)
        if s <= 0:
            return np.inf
        q = yVy - float(XVy @ np.linalg.solve(XVX, XVy))
        dfr = self.n - self.p
        if q <= 0 or dfr <= 0:
            return np.inf
        return dfr * np.log(q / dfr) + logdet + ld_a

    def optimise(self):
        k = len(self.sizes)
        if k == 0:
            return np.zeros(0)
        if k == 1:
            grid = np.concatenate([[0.0], 10.0 ** np.arange(-6, 4.01, 0.5)])
            vals = np.array([self.objective([g]) for g in grid])
            i = int(np.argmin(vals))
            best, bestv = grid[i], vals[i]
            lo, hi = grid[max(i - 1, 1)], grid[min(i + 1, len(grid) - 1)]
            if 0 < i < len(grid) - 1 and lo > 0:
                r = optimize.minimize_scalar(lambda lg: self.objective([10 ** lg]),
                                             bounds=(np.log10(lo), np.log10(hi)), method="bounded",
                                             options={"xatol": 1e-6})
                if r.fun < bestv:
                    best = 10 ** r.x
            return np.array([best])
        best, bestv = np.zeros(k), self.objective(np.zeros(k))
        starts = [(a, b) for a in (0.0, 0.1, 1.0, 10.0) for b in (0.0, 0.1, 1.0, 10.0)]
        for st in starts:
            r = optimize.minimize(lambda v: self.objective(np.maximum(v, 0)), np.array(st, float),
                                  method="L-BFGS-B", bounds=[(0, 1e4)] * k)
            if r.fun < bestv - 1e-9:
                best, bestv = np.maximum(r.x, 0), r.fun
        return best


def fit_model(y: np.ndarray, fixed: FixedDesign, rank: RankInfo | None = None) -> Fit:
    """Fit y ~ fixed design with random intercepts implied by the frame's units."""
    y = np.asarray(y, float)
    rank = rank or analyse_rank(fixed)
    frame = fixed.frame
    Xf = fixed.X.to_numpy(float)
    Xk = Xf[:, rank.kept]
    n, p_full = len(y), Xf.shape[1]
    note = ""
    blocks, names = [], []
    df = frame.df
    if df["upper"].notna().any():
        u = pd.factorize(df["upper"])[0]
        blocks.append(_indicator(u, u.max() + 1)); names.append("upper")
        if df["lower"].notna().any():
            lower_counts = df.groupby("lower").size()
            if (lower_counts > 1).any():
                s = pd.factorize(df["lower"])[0]
                blocks.append(_indicator(s, s.max() + 1)); names.append("lower")
            else:
                note = "lower-level unit has one observation per level; its variance is not separately identifiable from residual variance, so it was omitted"
    if not np.isfinite(y).all():
        return _failed(fixed, rank, n, "non-finite response")
    if np.ptp(y) <= 1e-12 * max(1.0, np.abs(y).max()):
        return _failed(fixed, rank, n, "response is constant; effects not estimable")
    if n - len(rank.kept) <= 0:
        return _failed(fixed, rank, n, "no residual degrees of freedom (saturated design)")
    reml = _Reml(y, Xk, blocks)
    lam = reml.optimise()
    pc = reml.pieces(lam)
    if pc is None or not np.isfinite(reml.objective(lam)):
        return _failed(fixed, rank, n, "REML optimisation failed")
    XVX, XVy, yVy, logdet, d, solve = pc
    A_inv = np.linalg.inv(XVX)
    bk = A_inv @ XVy
    dfr = n - len(rank.kept)
    s2 = float((yVy - XVy @ bk) / dfr)
    beta = np.zeros(p_full); beta[rank.kept] = bk
    cov = np.zeros((p_full, p_full)); cov[np.ix_(rank.kept, rank.kept)] = s2 * A_inv
    r = y - Xk @ bk
    if len(d):
        Zr = d * (reml.Z.T @ r)
        rc = r - reml.Z @ (d * solve(Zr))  # V^-1 r = conditional residual
    else:
        rc = r
    sizes_names = dict(zip(names, lam))
    variance = {k: float(v * s2) for k, v in sizes_names.items()}
    variance["residual"] = s2
    boundary = bool(len(lam) and (lam < 1e-8).any())
    if boundary:
        zero = [k for k, v in sizes_names.items() if v < 1e-8]
        boundary_note = f"variance component at boundary (0): {', '.join(zero)} (REML solution equals the model without it)"
        note = "; ".join(part for part in (note, boundary_note) if part)
    if "lower" in frame.df and df["lower"].notna().any() and "lower" not in names and not note:
        note = "lower-level unit adds no separately identifiable structure and was omitted"
    ddf = containment_ddf(fixed, rank.kept)
    model = ("OLS" if not names else "REML random intercept" if names == ["upper"]
             else "REML nested random intercepts")
    return Fit(beta, cov, fixed, rank, ddf, s2, variance, model, note.strip("; "), r, rc, Xf @ beta, n, boundary, True)


def _failed(fixed, rank, n, msg) -> Fit:
    p = len(fixed.names)
    return Fit(np.full(p, np.nan), np.full((p, p), np.nan), fixed, rank, {}, np.nan, {}, "not fitted", msg,
               np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), n, False, False)


def wald_term(fit: Fit, term: str) -> dict:
    out = dict(term=term, test="Wald F", statistic=np.nan, df=np.nan, ddf=np.nan, p_value=np.nan, note="")
    if not fit.ok:
        out["test"], out["note"] = "not fitted", fit.note
        return out
    cols = fit.fixed.terms.get(term, [])
    idx = [fit.fixed.names.index(c) for c in cols if fit.fixed.names.index(c) in fit.rank.kept]
    if not idx:
        out["test"], out["note"] = "not estimable", "all columns of this term are aliased with others"
        return out
    b, V = fit.beta[idx], fit.cov[np.ix_(idx, idx)]
    try:
        F = float(b @ np.linalg.solve(V, b)) / len(idx)
    except np.linalg.LinAlgError:
        out["note"] = "singular covariance"
        return out
    ddf = min(fit.ddf[j] for j in idx)
    out.update(statistic=F, df=len(idx), ddf=ddf, p_value=float(f_dist.sf(F, len(idx), ddf)))
    if len(idx) < len(cols):
        out["note"] = f"reduced df ({len(idx)} of {len(cols)} columns estimable)"
    return out


def contrast(fit: Fit, c: np.ndarray, conf: float = 0.95) -> dict:
    out = dict(estimate=np.nan, se=np.nan, t=np.nan, df=np.nan, p_value=np.nan, ci_low=np.nan, ci_high=np.nan,
               estimable=False)
    if not fit.ok or not fit.rank.estimable(c):
        return out
    est = float(c @ fit.beta)
    var = float(c @ fit.cov @ c)
    out.update(estimate=est, estimable=True)
    if not np.isfinite(var) or var <= 0:
        return out
    cols = [j for j in fit.rank.kept if abs(c[j]) > 1e-12]
    ddf = min((fit.ddf[j] for j in cols), default=np.nan)
    se = float(np.sqrt(var))
    tt = est / se
    q = t_dist.ppf(0.5 + conf / 2, ddf)
    out.update(se=se, t=tt, df=ddf, p_value=float(2 * t_dist.sf(abs(tt), ddf)), ci_low=est - q * se, ci_high=est + q * se)
    return out


def residual_diagnostics(fit: Fit, frame_df: pd.DataFrame) -> dict:
    """Shapiro-Wilk and Brown-Forsythe on conditional residuals + flags."""
    rc = fit.resid_conditional
    out = dict(shapiro_p=np.nan, brown_forsythe_p=np.nan, max_abs_std_resid=np.nan, flag="")
    if not fit.ok or len(rc) < 4 or not np.isfinite(rc).all():
        return out
    sd = np.sqrt(fit.sigma2) if fit.sigma2 and fit.sigma2 > 0 else np.nan
    if len(rc) <= 5000 and np.ptp(rc) > 0:
        out["shapiro_p"] = float(shapiro(rc).pvalue)
    key = frame_df["group"].astype(str) + ("|" + frame_df["time"].astype(str) if "time" in frame_df else "")
    groups = [rc[(key == k).to_numpy()] for k in key.unique()]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) >= 2 and all(np.ptp(g) >= 0 for g in groups):
        try:
            out["brown_forsythe_p"] = float(levene(*groups, center="median").pvalue)
        except Exception:
            pass
    out["max_abs_std_resid"] = float(np.abs(rc).max() / sd) if sd and np.isfinite(sd) and sd > 0 else np.nan
    flags = []
    if np.isfinite(out["shapiro_p"]) and out["shapiro_p"] < 0.05:
        flags.append("non-normal residuals")
    if np.isfinite(out["brown_forsythe_p"]) and out["brown_forsythe_p"] < 0.05:
        flags.append("unequal variances")
    if np.isfinite(out["max_abs_std_resid"]) and out["max_abs_std_resid"] > 3.5:
        flags.append("possible outlier")
    out["flag"] = "; ".join(flags)
    return out
