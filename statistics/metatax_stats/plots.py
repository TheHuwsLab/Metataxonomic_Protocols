"""Matplotlib figures (headless)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

DPI = 200


def set_dpi(v: int) -> None:
    global DPI
    DPI = v


def _colors(levels):
    cmap = plt.get_cmap("tab10" if len(levels) <= 10 else "tab20")
    return {lv: cmap(i % cmap.N) for i, lv in enumerate(levels)}


def _save(fig, path, rect=None):
    from pathlib import Path
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=rect)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def alpha_plot(values: pd.Series, frame, metric: str, path, connect: bool = True):
    df = frame.df.assign(v=values.reindex(frame.df.index))
    col = _colors(frame.group_levels)
    times = frame.time_levels if frame.has_time else [None]
    fig, ax = plt.subplots(figsize=(max(4, 1.2 * len(frame.group_levels) * len(times) + 1), 4.2))
    pos, rng = {}, np.random.default_rng(0)
    w = 0.8 / len(frame.group_levels)
    for ti, t in enumerate(times):
        for gi, g in enumerate(frame.group_levels):
            d = df[(df["group"] == g) & ((df["time"] == t) if t is not None else True)]
            x = ti + (gi - (len(frame.group_levels) - 1) / 2) * w
            if len(d):
                ax.boxplot(d["v"].dropna(), positions=[x], widths=w * 0.8, showfliers=False,
                           boxprops=dict(color=col[g]), medianprops=dict(color="black"),
                           whiskerprops=dict(color=col[g]), capprops=dict(color=col[g]))
                xs = x + rng.uniform(-w * 0.2, w * 0.2, len(d))
                ax.scatter(xs, d["v"], color=col[g], s=18, zorder=3, edgecolor="white", linewidth=0.4)
                for sid, xx in zip(d.index, xs):
                    pos[sid] = (xx, d.loc[sid, "v"])
    if connect and frame.has_time and frame.has_units:
        rep = df["lower"].where(df["lower"].notna(), df["upper"])
        for _, ids in df.groupby(rep).groups.items():
            ids = [i for i in sorted(ids, key=lambda i: frame.time_levels.index(df.loc[i, "time"])) if i in pos]
            if len(ids) > 1:
                ax.plot([pos[i][0] for i in ids], [pos[i][1] for i in ids], color="grey", lw=0.5, alpha=0.5, zorder=1)
    ax.set_xticks(range(len(times)))
    ax.set_xticklabels([str(t) if t is not None else "" for t in times])
    ax.set_ylabel(metric)
    ax.set_title(metric)
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=col[g], label=g) for g in frame.group_levels],
              fontsize=7, frameon=False)
    _save(fig, path)


def residual_plot(fitted, resid, metric: str, path):
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    stats.probplot(resid, plot=axes[0])
    axes[0].set_title(f"{metric}: normal Q-Q")
    axes[1].scatter(fitted, resid, s=14)
    axes[1].axhline(0, color="grey", lw=0.8)
    axes[1].set_xlabel("fitted"); axes[1].set_ylabel("conditional residual"); axes[1].set_title("residuals vs fitted")
    _save(fig, path)


def pcoa_plot(coords: pd.DataFrame, var: pd.Series, frame, path, title: str, connect: bool = True, ellipses: bool = True):
    df = frame.df.loc[coords.index]
    col = _colors(frame.group_levels)
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    tmark = {t: markers[i % len(markers)] for i, t in enumerate(frame.time_levels)} if frame.has_time else {}
    fig, ax = plt.subplots(figsize=(6, 5))
    xa, ya = coords.columns[0], coords.columns[1] if coords.shape[1] > 1 else coords.columns[0]
    if connect and frame.has_time and frame.has_units:
        rep = df["lower"].where(df["lower"].notna(), df["upper"])
        for _, ids in df.groupby(rep).groups.items():
            ids = sorted(ids, key=lambda i: frame.time_levels.index(df.loc[i, "time"]))
            if len(ids) > 1:
                ax.plot(coords.loc[ids, xa], coords.loc[ids, ya], color="grey", lw=0.5, alpha=0.5, zorder=1)
    for g in frame.group_levels:
        for t in (frame.time_levels if frame.has_time else [None]):
            m = (df["group"] == g) & ((df["time"] == t) if t is not None else True)
            if m.any():
                ax.scatter(coords.loc[m.values, xa], coords.loc[m.values, ya], color=col[g], marker=tmark.get(t, "o"),
                           s=36, edgecolor="white", linewidth=0.4, zorder=3, label=f"{g} {t}" if t is not None else g)
        if ellipses:
            _ellipse(ax, coords.loc[(df["group"] == g).values, [xa, ya]].to_numpy(), col[g])
    ax.set_xlabel(f"{xa} ({var.get(xa, np.nan) * 100:.1f}%)")
    ax.set_ylabel(f"{ya} ({var.get(ya, np.nan) * 100:.1f}%)")
    ax.set_title(title)
    ax.legend(fontsize=6, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    _save(fig, path)


def _ellipse(ax, pts, color):
    if len(pts) < 4:
        return
    from matplotlib.patches import Ellipse
    cov = np.cov(pts.T)
    w, v = np.linalg.eigh(cov)
    if (w <= 0).any():
        return
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(pts.mean(0), 2 * 1.5 * np.sqrt(w[1]), 2 * 1.5 * np.sqrt(w[0]), angle=ang, fill=False,
                         edgecolor=color, lw=0.8, ls="--", zorder=2))


def volcano_plot(df: pd.DataFrame, title: str, path, alpha: float = 0.05):
    d = df.dropna(subset=["log2_estimate", "p_value"])
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sig = d["q_value"] < alpha
    ax.scatter(d.loc[~sig, "log2_estimate"], -np.log10(d.loc[~sig, "p_value"]), s=12, color="lightgrey")
    ax.scatter(d.loc[sig, "log2_estimate"], -np.log10(d.loc[sig, "p_value"]), s=16, color="firebrick")
    for _, r in d[sig].nsmallest(8, "p_value").iterrows():
        ax.annotate(str(r["taxon"])[:30], (r["log2_estimate"], -np.log10(r["p_value"])), fontsize=6)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("log2 CLR difference"); ax.set_ylabel("-log10 p"); ax.set_title(title, fontsize=9)
    _save(fig, path)


def stacked_bar(rel: pd.DataFrame, frame, path, top_n: int = 15, by: str = "group"):
    """Mean relative abundance per group[/time] (by='group') or per sample (by='sample')."""
    top = rel.mean(0).nlargest(top_n).index
    plot = rel[top].copy()
    plot["Other"] = 1 - plot.sum(1)
    df = frame.df.loc[plot.index]
    if by == "group":
        key = df["group"].astype(str) + ((" " + df["time"].astype(str)) if frame.has_time else "")
        plot = plot.groupby(key.values, sort=False).mean()
        order = [f"{g}" + (f" {t}" if frame.has_time else "") for t in (frame.time_levels or [None]) for g in frame.group_levels]
        plot = plot.reindex([o for o in order if o in plot.index])
    else:
        order = df.sort_values(["group"] + (["time"] if frame.has_time else [])).index
        plot = plot.loc[order]
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(plot) + 4), 5))
    bottom = np.zeros(len(plot))
    for i, c in enumerate(plot.columns):
        ax.bar(range(len(plot)), plot[c], bottom=bottom, color="lightgrey" if c == "Other" else cmap(i % 20),
               label=str(c)[:40], width=0.85)
        bottom += plot[c].to_numpy()
    ax.set_xticks(range(len(plot))); ax.set_xticklabels(plot.index, rotation=70, ha="right", fontsize=7)
    ax.set_ylabel("relative abundance")
    ax.legend(fontsize=6, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    _save(fig, path)


def heatmap(clr: pd.DataFrame, frame, taxa: list, path, title: str = ""):
    if not taxa:
        return
    df = frame.df.loc[clr.index]
    order = df.sort_values(["group"] + (["time"] if frame.has_time else [])).index
    m = clr.loc[order, taxa].T
    z = (m.sub(m.mean(1), axis=0)).div(m.std(1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(max(6, 0.22 * len(order) + 4), max(3, 0.22 * len(taxa) + 2)))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    ax.set_yticks(range(len(taxa))); ax.set_yticklabels([str(t)[:40] for t in taxa], fontsize=6)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{df.loc[i, 'group']}" + (f" {df.loc[i, 'time']}" if frame.has_time else "") for i in order],
                       rotation=90, fontsize=5)
    fig.colorbar(im, ax=ax, label="CLR z-score")
    ax.set_title(title, fontsize=9)
    _save(fig, path)


def rarefaction_plot(curves: pd.DataFrame, depth, path):
    if curves.empty:
        return
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for _, d in curves.groupby("sample_id"):
        ax.plot(d["depth"], d["observed_features"], color="steelblue", lw=0.6, alpha=0.6)
    if depth:
        ax.axvline(depth, color="firebrick", ls="--", lw=1, label=f"rarefaction depth {depth}")
        ax.legend(fontsize=7, frameon=False)
    ax.set_xlabel("sequencing depth"); ax.set_ylabel("observed features"); ax.set_title("Rarefaction curves")
    _save(fig, path)


def rarefaction_alpha_comparison(raw: pd.DataFrame, draws: list[pd.DataFrame], frame, metric: str, path):
    """Compare unrarefied group means with draw-to-draw rarefaction variation."""
    meta = frame.df
    strata = [(g, t) for t in (frame.time_levels if frame.has_time else [None])
              for g in frame.group_levels]
    raw_means, draw_means = [], []
    for group, time in strata:
        mask = meta["group"].eq(group)
        if time is not None:
            mask &= meta["time"].eq(time)
        sample_ids = meta.index[mask]
        raw_means.append(float(raw.reindex(sample_ids)[metric].mean()) if len(sample_ids) else np.nan)
        draw_means.append(
            [float(draw.reindex(sample_ids)[metric].mean()) for draw in draws if len(sample_ids)]
        )
    centers = np.array([np.nanmean(x) if x else np.nan for x in draw_means])
    lows = np.array([np.nanquantile(x, 0.05) if x else np.nan for x in draw_means])
    highs = np.array([np.nanquantile(x, 0.95) if x else np.nan for x in draw_means])
    x = np.arange(len(strata))
    fig, ax = plt.subplots(figsize=(max(6, 0.75 * len(strata) + 2), 4.2))
    ax.vlines(x, lows, highs, color="steelblue", linewidth=3, alpha=0.7,
              label="Rarefied draw means (5th–95th percentile)")
    ax.scatter(x, centers, color="steelblue", marker="o", s=35, zorder=3,
               label="Mean across rarefied draws")
    ax.scatter(x, raw_means, color="darkorange", marker="D", s=32, zorder=4,
               label="Unrarefied mean")
    ax.set_xticks(x)
    ax.set_xticklabels([str(g) if t is None else f"{g} / {t}" for g, t in strata], rotation=45, ha="right")
    ax.set_ylabel(metric)
    ax.set_title(f"{metric}: rarefaction sensitivity (descriptive means)")
    ax.legend(fontsize=7, frameon=False)
    _save(fig, path)


def rarefaction_taxa_prevalence(summary: pd.DataFrame, path, top_n: int = 15):
    """Plot independent-unit taxon detection prevalence across rarefaction draws."""
    if summary.empty:
        return
    widths = summary.assign(interval_width=summary["upper_95pct"] - summary["lower_5pct"])
    top = (widths.groupby("taxon")["interval_width"].mean()
           .nlargest(top_n).index.tolist())
    data = summary[summary["taxon"].isin(top)].copy()
    strata = list(dict.fromkeys(data["stratum"]))
    taxa_order = list(reversed(top))
    ncols = min(3, max(1, len(strata)))
    nrows = int(np.ceil(len(strata) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, max(3.0, 0.36 * len(top) + 1.5)),
                             squeeze=False, sharex=True)
    for ax, stratum in zip(axes.flat, strata):
        d = data[data["stratum"] == stratum].set_index("taxon").reindex(taxa_order)
        x = d["mean_prevalence"].to_numpy(float)
        lo = d["lower_5pct"].to_numpy(float)
        hi = d["upper_95pct"].to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi)
        y = np.arange(len(d))
        ax.hlines(y[valid], lo[valid], hi[valid], color="lightseagreen", linewidth=1.5)
        ax.vlines(np.concatenate([lo[valid], hi[valid]]),
                  np.concatenate([y[valid] - 0.12, y[valid] - 0.12]),
                  np.concatenate([y[valid] + 0.12, y[valid] + 0.12]),
                  color="lightseagreen", linewidth=1)
        if "unrarefied_prevalence" in d:
            baseline = d["unrarefied_prevalence"].to_numpy(float)
            ax.scatter(baseline[valid], y[valid], marker="D", s=18, color="darkorange", zorder=3)
        ax.scatter(x[valid], y[valid], marker="o", s=18, color="teal",
                   zorder=4)
        ax.set_yticks(y)
        ax.set_yticklabels([str(t)[:42] for t in taxa_order], fontsize=7)
        ax.set_xlim(0, 1)
        ax.set_title(stratum, fontsize=9)
        ax.set_xlabel("Independent-unit prevalence")
    for ax in axes.flat[len(strata):]:
        ax.set_visible(False)
    fig.suptitle("Taxa with most variable detection prevalence after rarefaction", fontsize=10, y=0.995)
    fig.legend(
        handles=[
            plt.Line2D([], [], marker="o", ls="", color="teal", label="Rarefied mean (range: 5th–95th percentile)"),
            plt.Line2D([], [], marker="D", ls="", color="darkorange", label="Unrarefied"),
        ],
        loc="upper center", bbox_to_anchor=(0.5, 0.96), ncol=2, fontsize=7, frameon=False,
    )
    _save(fig, path, rect=(0, 0, 1, 0.90))


def depth_plot(depths: pd.Series, depth, path):
    s = depths.sort_values()
    fig, ax = plt.subplots(figsize=(max(5, 0.2 * len(s) + 2), 3.6))
    ax.bar(range(len(s)), s.values, color="steelblue")
    if depth:
        ax.axhline(depth, color="firebrick", ls="--", lw=1)
    ax.set_xticks(range(len(s))); ax.set_xticklabels(s.index, rotation=90, fontsize=5)
    ax.set_ylabel("reads"); ax.set_title("Library sizes")
    _save(fig, path)


def first_pass_ordination(coords: pd.DataFrame, frame, path, title: str, variance=None):
    if coords.shape[1] < 2:
        return
    meta = frame.df.loc[coords.index]
    colors = _colors(frame.group_levels)
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    time_markers = {t: markers[i % len(markers)] for i, t in enumerate(frame.time_levels)}
    fig, ax = plt.subplots(figsize=(6, 5))
    for group in frame.group_levels:
        for time in (frame.time_levels if frame.has_time else [None]):
            mask = meta["group"] == group
            if time is not None:
                mask &= meta["time"] == time
            if not mask.any():
                continue
            ax.scatter(coords.loc[mask, coords.columns[0]], coords.loc[mask, coords.columns[1]],
                       color=colors[group], marker=time_markers.get(time, "o"), s=38, edgecolor="white",
                       linewidth=0.4, label=f"{group} / {time}" if time is not None else group)
    xlab, ylab = coords.columns[:2]
    if variance is not None and len(variance) >= 2:
        xlab += f" ({variance[0] * 100:.1f}%)"
        ylab += f" ({variance[1] * 100:.1f}%)"
    ax.set_xlabel(xlab); ax.set_ylabel(ylab); ax.set_title(title)
    ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    _save(fig, path)


def mean_cpm_barplot(means: pd.DataFrame, group_cols: list[str], path, top_n: int = 15):
    if means.empty or not group_cols:
        return
    taxa_cols = [c for c in means.columns if c not in {*group_cols, "n_samples", "n_upper"}]
    if not taxa_cols:
        return
    top = means[taxa_cols].mean(0).nlargest(top_n).index
    data = means.set_index(group_cols)[list(top)]
    labels = ["/".join(map(str, ix if isinstance(ix, tuple) else (ix,))) for ix in data.index]
    fig, ax = plt.subplots(figsize=(max(7, 0.5 * len(data) + 3), 5))
    bottom = np.zeros(len(data))
    cmap = plt.get_cmap("tab20")
    for i, taxon in enumerate(top):
        values = data[taxon].to_numpy(float)
        ax.bar(np.arange(len(data)), values, bottom=bottom, width=0.82, color=cmap(i % 20),
               label=str(taxon)[:38])
        bottom += values
    ax.set_xticks(np.arange(len(data))); ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Mean CPM (counts per million)"); ax.set_title(f"Top {len(top)} taxa by mean CPM")
    ax.legend(fontsize=6, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    _save(fig, path)
