"""Alpha diversity metrics and repeated-measures-aware models."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .contrasts import adjust_p, contrast_family
from .models import fit_model, residual_diagnostics, wald_term

RICHNESS = {"observed_features", "chao1"}
ALL_METRICS = ["observed_features", "shannon", "simpson", "pielou_evenness", "chao1", "goods_coverage"]
FEATURE_ONLY = {"chao1", "goods_coverage"}


def compute_alpha(counts: np.ndarray, ids, metrics: list[str], feature_level: bool) -> pd.DataFrame:
    """counts: samples x taxa (integers). Natural-log Shannon; Simpson is 1-D."""
    c = counts.astype(float)
    n = c.sum(1)
    obs = (c > 0).sum(1)
    out = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        p = c / n[:, None]
        h = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
        out["observed_features"] = obs.astype(float)
        out["shannon"] = h
        out["simpson"] = 1 - (p ** 2).sum(1)
        out["pielou_evenness"] = np.where(obs > 1, h / np.log(np.maximum(obs, 2)), np.nan)
        f1, f2 = (c == 1).sum(1), (c == 2).sum(1)
        out["chao1"] = obs + f1 * (f1 - 1) / (2.0 * (f2 + 1))
        out["goods_coverage"] = 1 - f1 / n
    df = pd.DataFrame(out, index=pd.Index(ids))
    keep = [m for m in metrics if m in df and (feature_level or m not in FEATURE_ONLY)]
    return df[keep]


def transform_name(metric: str, mode: str) -> str:
    if mode == "log" or (mode == "auto" and metric in RICHNESS):
        return "log"
    return "none"


def _tf(y: pd.Series, name: str) -> np.ndarray:
    v = y.to_numpy(float)
    return np.log(v) if name == "log" else v


def run_alpha(draws: list[pd.DataFrame], frame, fixed, rank, mode: str, adjust: str, notes, scope: str,
              external: pd.DataFrame | None = None) -> dict:
    mean = sum(draws) / len(draws) if draws else pd.DataFrame(index=frame.df.index)
    if external is not None and not external.empty:
        mean = mean.join(external, how="left")
    idx = frame.df.index
    res = {"omnibus": [], "contrasts": [], "variance": [], "diagnostics": [], "descriptives": [], "stability": []}
    resid = {}
    for m in mean.columns:
        y = mean[m].reindex(idx)
        tf = transform_name(m, mode)
        ok = y.notna() & np.isfinite(y)
        if tf == "log" and (y[ok] <= 0).any():
            tf = "none"
            notes.warn(scope, f"{m}: non-positive values; log transform disabled.")
        if ok.sum() < len(idx):
            notes.warn(scope, f"{m}: {int((~ok).sum())} sample(s) without a value dropped from this model.")
        if ok.sum() < 4:
            continue
        fr, fx = (frame, fixed) if ok.all() else (frame.subset(idx[ok]), None)
        if fx is None:
            from .design import FixedDesign, analyse_rank
            fx = FixedDesign(fr)
            rk = analyse_rank(fx)
        else:
            rk = rank
        fit = fit_model(_tf(y[ok], tf), fx, rk)
        base = dict(metric=m, transform=tf, model=fit.model)
        for term in fx.terms:
            w = wald_term(fit, term)
            w.update(base, n=int(ok.sum()))
            # sensitivity to rarefaction draw
            ps = []
            for d in (draws if len(draws) > 1 and m not in {"chao1", "goods_coverage"} else []):
                if m in d:
                    yy = d[m].reindex(idx)[ok]
                    if np.isfinite(yy).all() and not (tf == "log" and (yy <= 0).any()):
                        f2 = fit_model(_tf(yy, tf), fx, rk)
                        ps.append(wald_term(f2, term)["p_value"])
            ps = [p for p in ps if np.isfinite(p)]
            w["p_min_over_draws"] = min(ps) if ps else np.nan
            w["p_max_over_draws"] = max(ps) if ps else np.nan
            res["omnibus"].append(w)
        fam = contrast_family(fit, frame)
        if not fam.empty:
            fam.insert(0, "metric", m)
            fam["transform"] = tf
            if tf == "log":
                fam["ratio"] = np.exp(fam["estimate"])
                fam["ratio_ci_low"] = np.exp(fam["ci_low"])
                fam["ratio_ci_high"] = np.exp(fam["ci_high"])
            res["contrasts"].append(fam)
        var = dict(metric=m, **{f"var_{k}": v for k, v in fit.variance.items()}, icc=fit.icc, boundary=fit.boundary,
                   note=fit.note)
        res["variance"].append(var)
        dg = residual_diagnostics(fit, fr.df)
        res["diagnostics"].append(dict(metric=m, transform=tf, **dg))
        resid[m] = (fit.fitted, fit.resid_conditional)
        d = fr.df.assign(v=y[ok].to_numpy())
        keys = ["group"] + (["time"] if frame.has_time else [])
        desc = d.groupby(keys, observed=True)["v"].agg(n="count", mean="mean", sd="std", median="median").reset_index()
        desc.insert(0, "metric", m)
        res["descriptives"].append(desc)
    out = {k: pd.concat(v, ignore_index=True) if v and isinstance(v[0], pd.DataFrame) else pd.DataFrame(v)
           for k, v in res.items() if k != "stability"}
    om = out["omnibus"]
    if not om.empty:
        om["p_adj_across_metrics"] = np.nan
        for _, ix in om.groupby("term").groups.items():
            om.loc[ix, "p_adj_across_metrics"] = adjust_p(om.loc[ix, "p_value"].to_numpy(float), "fdr_bh")
    ct = out["contrasts"]
    if not ct.empty:
        ct["p_adj"] = np.nan
        for _, ix in ct.groupby(["metric", "family"]).groups.items():
            ct.loc[ix, "p_adj"] = adjust_p(ct.loc[ix, "p_value"].to_numpy(float), adjust)
        ct["p_adjust_method"] = f"{adjust} within metric x family"
        ct["df_method"] = "containment (min over involved coefficients)"
    out["residuals"] = resid
    out["mean_values"] = mean
    return out
