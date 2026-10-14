"""Beta diversity: distances, PCoA, PERMANOVA (design-aware) and PERMDISP.

Design-aware permutation tests
------------------------------
* No repeated units (independent samples): Freedman-Lane permutation of the
  residuals of the reduced model; type-II marginal tests of group, time,
  group:time and covariates.
* Repeated/clustered units (split-plot): the multivariate analogue of a
  split-plot ANOVA on principal coordinates.
    - group (between units): unit means, group labels permuted among whole units;
    - time (within units): within-unit centred data, time labels permuted within units
      (exact enumeration of all swaps when 2 levels and <= 15 units);
    - group:time: within-unit centred data with the pooled time profile removed,
      then group labels permuted among whole treatment-assignment units. This is a
      custom directional interaction test; it is not a drop-in implementation of
      a named package's PERMANOVA and requires validation for each design family.
  Replicates within a unit x time cell are averaged first; only units with every time
  level are used for time/interaction tests. With nested subjects, observations are
  averaged at assignment-unit x time, so inference is at the assignment-unit level.

Non-Euclidean distances are Lingoes-corrected before inference. This means tests
are performed on the corrected Euclidean representation, not the unmodified
distance matrix; the correction constant and original negative-eigenvalue fraction
are recorded.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata

from .design import FixedDesign, min_attainable_p

BUILTIN_METRICS = {"braycurtis", "jaccard", "aitchison"}
TOL = 1e-10


# -----------------------------------------------------------------------------
# Distances and ordination
# -----------------------------------------------------------------------------

def distance_from_table(table: pd.DataFrame, metric: str) -> pd.DataFrame:
    """table: taxa x samples (relative abundance for braycurtis, 0/1 for jaccard, CLR for aitchison)."""
    x = table.T.to_numpy(float)
    labels = table.columns.astype(str)
    metric = metric.lower()
    if metric == "braycurtis":
        d = pdist(x, "braycurtis")
    elif metric == "jaccard":
        d = pdist(x > 0, "jaccard")
    elif metric == "aitchison":
        d = pdist(x, "euclidean")
    else:
        raise ValueError(f"Unknown beta metric '{metric}'")
    d = np.nan_to_num(d, nan=0.0)  # two empty samples: identical
    return pd.DataFrame(squareform(d), index=labels, columns=labels)


def read_distance_matrix(path) -> pd.DataFrame:
    from .io import read_delimited
    df = read_delimited(path)
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    d = df.apply(pd.to_numeric, errors="coerce")
    if d.isna().any().any() or d.shape[0] != d.shape[1]:
        raise SystemExit(f"External distance matrix {path} must be a square numeric matrix with sample IDs.")
    return d.loc[d.index, d.index]


def gower_coords(D: np.ndarray):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    B = (B + B.T) / 2
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    tol = TOL * max(abs(w).max(), 1e-300)
    pos = w > tol
    Y = V[:, pos] * np.sqrt(w[pos])
    neg_frac = float(np.abs(w[w < -tol]).sum() / max(np.abs(w).sum(), 1e-300))
    return Y, {"eigenvalues": w, "n_axes": int(pos.sum()), "negative_fraction": neg_frac}


def lingoes_correct(D: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Return a Euclidean Lingoes-corrected distance matrix, correction c, and original negative mass."""
    _, info = gower_coords(D)
    w = info["eigenvalues"]
    minimum = float(w.min()) if len(w) else 0.0
    c = max(0.0, -minimum)
    if c == 0:
        return D.copy(), 0.0, info["negative_fraction"]
    corrected_sq = D ** 2 + 2.0 * c * (1.0 - np.eye(D.shape[0]))
    np.fill_diagonal(corrected_sq, 0.0)
    return np.sqrt(np.maximum(corrected_sq, 0.0)), c, info["negative_fraction"]


def pcoa(dist: pd.DataFrame, n_axes: int = 3):
    if dist.shape[0] < 3:
        raise ValueError("PCoA needs at least 3 samples")
    corrected, _, negative_fraction = lingoes_correct(dist.to_numpy(float))
    Y, info = gower_coords(corrected)
    w = info["eigenvalues"]
    posw = w[w > TOL * max(abs(w).max(), 1e-300)]
    if len(posw) == 0:
        raise ValueError("PCoA produced no positive eigenvalues")
    k = min(n_axes, Y.shape[1])
    cols = [f"PC{i + 1}" for i in range(k)]
    coords = pd.DataFrame(Y[:, :k], index=dist.index, columns=cols)
    var = pd.Series(posw[:k] / posw.sum(), index=cols, name="variance_explained")
    return coords, var, negative_fraction


# -----------------------------------------------------------------------------
# Vectorised permutation helpers
# -----------------------------------------------------------------------------

def _chunks(total: int, size: int):
    done = 0
    while done < total:
        b = min(size, total - done)
        yield b
        done += b


def _p_value(obs: float, null: np.ndarray, exact: bool) -> float:
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(obs):
        return np.nan
    ge = float(np.sum(null >= obs - 1e-9 * max(abs(obs), 1.0)))
    return ge / len(null) if exact else (1.0 + ge) / (1.0 + len(null))


def _ss_between(G: np.ndarray, L: np.ndarray, ng: int) -> np.ndarray:
    """Between-group sums of squares for each row of label matrix L (B x n) given Gram matrix G."""
    out = np.zeros(L.shape[0])
    for g in range(ng):
        Z = (L == g).astype(float)
        cnt = Z.sum(1)
        q = ((Z @ G) * Z).sum(1)
        out += np.divide(q, cnt, out=np.zeros_like(q), where=cnt > 0)
    return out


def _permuted_labels(L0: np.ndarray, B: int, rng) -> np.ndarray:
    return rng.permuted(np.tile(L0, (B, 1)), axis=1)


def _gram(Y: np.ndarray) -> np.ndarray:
    Yc = Y - Y.mean(0, keepdims=True)
    return Yc @ Yc.T


def oneway_perm(Yu: np.ndarray, labels: np.ndarray, nperm: int, rng, size: int = 1000) -> dict | None:
    """One-way distance-based test on unit-level rows; labels permuted freely among units."""
    codes, uniq = pd.factorize(pd.Series(labels))
    n, g = len(codes), len(uniq)
    if g < 2 or n <= g:
        return None
    G = _gram(Yu)
    sst = float(np.trace(G))
    ssb = float(_ss_between(G, codes[None, :], g)[0])
    ssw = sst - ssb
    F = (ssb / (g - 1)) / (ssw / (n - g)) if ssw > TOL * max(sst, 1) else np.inf
    null = []
    for b in _chunks(nperm, size):
        sb = _ss_between(G, _permuted_labels(codes, b, rng), g)
        sw = sst - sb
        null.append(np.where(sw > TOL * max(sst, 1), (sb / (g - 1)) / (np.where(sw > 0, sw, 1) / (n - g)), np.inf))
    null = np.concatenate(null)
    arrangements, minp = min_attainable_p(list(np.bincount(codes)))
    return dict(n=n, groups=g, pseudo_F=F, R2=ssb / sst if sst > 0 else np.nan, df_num=g - 1, df_den=n - g,
                permutations=nperm, p_value=_p_value(F, null, False),
                min_attainable_p=max(minp, 1.0 / (nperm + 1)), arrangements=arrangements)


def permdisp(Yu: np.ndarray, labels: np.ndarray, nperm: int, rng, size: int = 1000) -> dict | None:
    codes, uniq = pd.factorize(pd.Series(labels))
    n, g = len(codes), len(uniq)
    if g < 2 or n <= g:
        return None
    def distances_to_centroid(Y, L):
        z = np.empty_like(L, dtype=float)
        for k in range(g):
            m = L == k
            if m.sum() < 2:
                z[m] = np.nan
                continue
            centroid = Y[m].mean(0)
            z[m] = np.linalg.norm(Y[m] - centroid, axis=1)
        return z

    def fstat(zv, L):
        # L: (B x n) labels
        B = L.shape[0]
        tot = np.nanmean(zv, axis=1)
        ssb = np.zeros(L.shape[0]); ssw = np.zeros(L.shape[0])
        for k in range(g):
            Zk = (L == k).astype(float)
            cnt = Zk.sum(1)
            valid = np.isfinite(zv)
            Zk *= valid
            cnt = Zk.sum(1)
            mean = np.divide(np.nansum(Zk * zv, axis=1), cnt, out=np.zeros(L.shape[0]), where=cnt > 0)
            ssb += cnt * (mean - tot) ** 2
            ssw += np.nansum(Zk * (zv - mean[:, None]) ** 2, axis=1)
        return (ssb / (g - 1)) / np.where(ssw > 0, ssw / (n - g), np.nan)

    z0 = distances_to_centroid(Yu, codes)
    F = float(fstat(z0[None, :], codes[None, :])[0])
    null_batches = []
    for b in _chunks(nperm, size):
        L = _permuted_labels(codes, b, rng)
        z = np.vstack([distances_to_centroid(Yu, labels_b) for labels_b in L])
        null_batches.append(fstat(z, L))
    null = np.concatenate(null_batches)
    means = {str(u): float(z0[codes == k].mean()) for k, u in enumerate(uniq)}
    return dict(dispersion_F=F, dispersion_p=_p_value(F, null, False),
                mean_distance_to_centroid="; ".join(f"{k}={v:.4g}" for k, v in means.items()))


def _orth(X: np.ndarray) -> np.ndarray:
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    r = int((s > max(X.shape) * np.finfo(float).eps * (s[0] if len(s) else 1) * 10).sum())
    return U[:, :r]


def fl_term_test(Y: np.ndarray, Xa: np.ndarray, Xb: np.ndarray, Xall: np.ndarray, nperm: int, rng,
                 size: int = 200) -> dict | None:
    """Freedman-Lane test of the extra terms in Xa over reduced model Xb; error from Xall."""
    n = Y.shape[0]
    Qa, Qb, Qf = _orth(Xa), _orth(Xb), _orth(Xall)
    df_t, df_r = Qa.shape[1] - Qb.shape[1], n - Qf.shape[1]
    if df_t <= 0 or df_r <= 0:
        return None

    def stat(Yb):  # (B, n, k)
        sa = (np.einsum("np,bnk->bpk", Qa, Yb) ** 2).sum((1, 2))
        sb = (np.einsum("np,bnk->bpk", Qb, Yb) ** 2).sum((1, 2))
        sf = (np.einsum("np,bnk->bpk", Qf, Yb) ** 2).sum((1, 2))
        tot = (Yb ** 2).sum((1, 2))
        ssr = tot - sf
        return (sa - sb) / df_t / np.where(ssr > TOL * np.maximum(tot, 1e-300), ssr / df_r, np.nan), sa - sb, tot

    F0, ss_t, tot = stat(Y[None])
    fitted = Qb @ (Qb.T @ Y)
    R = Y - fitted
    nulls = []
    for b in _chunks(nperm, size):
        idx = np.argsort(rng.random((b, n)), axis=1)
        nulls.append(stat(fitted[None] + R[idx])[0])
    null = np.concatenate(nulls)
    return dict(n=n, df_num=df_t, df_den=df_r, pseudo_F=float(F0[0]),
                R2=float(ss_t[0] / tot[0]) if tot[0] > 0 else np.nan, permutations=nperm,
                p_value=_p_value(float(F0[0]), null, False), min_attainable_p=1.0 / (nperm + 1))


def independent_term_tests(Y: np.ndarray, fixed: FixedDesign, nperm: int, rng) -> pd.DataFrame:
    X = fixed.X.to_numpy(float)
    names = fixed.names
    terms = fixed.terms

    def cols(*ts):
        idx = [0] + [names.index(c) for t in ts for c in terms[t]]
        return X[:, sorted(set(idx))]

    all_terms = list(terms)
    covs = [t for t in all_terms if t not in ("group", "time", "group:time")]
    rows = []
    for t in all_terms:
        if t == "group":
            red = [x for x in ["time", *covs] if x in terms]
            a, b = cols("group", *red), cols(*red)
        elif t == "time":
            red = ["group", *covs]
            a, b = cols("time", *red), cols(*red)
        elif t == "group:time":
            red = ["group", "time", *covs]
            a, b = cols("group:time", *red), cols(*red)
        else:
            red = [x for x in all_terms if x != t]
            a, b = cols(*all_terms), cols(*red)
        res = fl_term_test(Y, a, b, cols(*all_terms), nperm, rng)
        if res is None:
            rows.append(dict(term=t, note="not estimable"))
            continue
        res.update(term=t, design="independent samples; Freedman-Lane residual permutation (type II)", note="")
        rows.append(res)
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Split-plot (repeated / clustered) restricted permutation tests
# -----------------------------------------------------------------------------

class UnitData:
    """Cell-level data (one row per repeated unit x time) for restricted permutation tests."""

    def __init__(self, Y: np.ndarray, df: pd.DataFrame, has_time: bool, time_levels: list):
        self.ok, self.reason = True, ""
        d = df.assign(_i=np.arange(len(df)))
        # Permutation units must match treatment assignment. When subjects are
        # nested inside a treatment-assignment cluster, average at cluster x
        # time rather than treating nested subjects as independent clusters.
        d["_rep"] = d["upper"]
        if (d.groupby("upper")["group"].nunique() > 1).any():
            self.ok, self.reason = False, "group varies within units"
            return
        d["_t"] = d["time"] if has_time else "all"
        grp = d.groupby(["_rep", "_t"], sort=False)
        cells = grp.agg(upper=("upper", "first"), group=("group", "first"), idx=("_i", list)).reset_index()
        self.n_replicate_cells = int((cells["idx"].map(len) > 1).sum())
        self.Y = np.vstack([Y[i].mean(0) for i in cells["idx"]])
        self.cells = cells
        self.has_time, self.tl = has_time, time_levels
        self.M = len(time_levels) if has_time else 1
        counts = cells.groupby("_rep")["_t"].nunique()
        self.complete_reps = counts.index[counts == self.M] if has_time else counts.index
        self.n_incomplete = int(len(counts) - len(self.complete_reps))

    def complete(self):
        m = self.cells["_rep"].isin(self.complete_reps).to_numpy()
        return self.cells.loc[m].reset_index(drop=True), self.Y[m]

    def at_time(self, t):
        m = (self.cells["_t"] == t).to_numpy()
        return self.cells.loc[m].reset_index(drop=True), self.Y[m]


def _unit_means(cells: pd.DataFrame, Y: np.ndarray):
    codes, uniq = pd.factorize(cells["upper"])
    Yu = np.vstack([Y[codes == k].mean(0) for k in range(len(uniq))])
    g = cells.groupby("upper", sort=False)["group"].first().reindex(uniq).to_numpy()
    return Yu, g, codes


def group_test(cells, Y, nperm, rng):
    Yu, g, _ = _unit_means(cells, Y)
    return oneway_perm(Yu, g, nperm, rng)


def _time_perm_labels(rep_rows: np.ndarray, tcodes: np.ndarray, nperm: int, rng, M: int):
    """Label matrix with time labels permuted within each repeated unit. Returns (labels, exact)."""
    nR = rep_rows.shape[0]
    n = tcodes.shape[0]
    if M == 2 and nR <= 15:
        bits = ((np.arange(2 ** nR)[:, None] >> np.arange(nR)) & 1).astype(bool)
        L = np.tile(tcodes, (len(bits), 1))
        a, b = rep_rows[:, 0], rep_rows[:, 1]
        L[:, a] = np.where(bits, tcodes[b], tcodes[a])
        L[:, b] = np.where(bits, tcodes[a], tcodes[b])
        return L, True
    L = np.tile(tcodes, (nperm, 1))
    for r in range(nR):
        rows = rep_rows[r]
        L[:, rows] = rng.permuted(np.tile(tcodes[rows], (nperm, 1)), axis=1)
    return L, False


def _rep_matrix(cells: pd.DataFrame) -> np.ndarray:
    return np.vstack([np.asarray(ix) for _, ix in cells.groupby("_rep", sort=False).indices.items()])


def _within_center(Y: np.ndarray, cells: pd.DataFrame) -> np.ndarray:
    codes = pd.factorize(cells["_rep"])[0]
    cnt = np.bincount(codes)
    means = np.vstack([np.bincount(codes, Y[:, j]) / cnt for j in range(Y.shape[1])]).T
    return Y - means[codes]


def time_test(cells, Y, tl, nperm, rng):
    """Within-unit time test; all units must be complete (each time level exactly once)."""
    M = len(tl)
    tcodes = pd.Categorical(cells["_t"], categories=tl).codes.astype(int)
    nR = cells["_rep"].nunique()
    if nR < 3:
        return None
    Yw = _within_center(Y, cells)
    G = Yw @ Yw.T
    sst = float(np.trace(G))
    rep_rows = _rep_matrix(cells)
    ssb0 = float(_ss_between(G, tcodes[None, :], M)[0])
    dfn, dfd = M - 1, (nR - 1) * (M - 1)

    def F(ssb):
        ssw = sst - ssb
        return np.where(ssw > TOL * max(sst, 1), (ssb / dfn) / (np.where(ssw > 0, ssw, 1) / dfd), np.inf)

    L, exact = _time_perm_labels(rep_rows, tcodes, nperm, rng, M)
    null = np.concatenate([F(_ss_between(G, L[i:i + 1000], M)) for i in range(0, len(L), 1000)])
    F0 = float(F(np.array([ssb0]))[0])
    return dict(n=len(cells), n_units=nR, pseudo_F=F0, R2=ssb0 / sst if sst > 0 else np.nan, df_num=dfn, df_den=dfd,
                permutations=len(L), exact=exact,
                min_attainable_p=(2.0 / len(L)) if exact else 1.0 / (nperm + 1),
                p_value=_p_value(F0, null, exact))


def interaction_test(cells, Y, tl, nperm, rng):
    M = len(tl)
    tcodes = pd.Categorical(cells["_t"], categories=tl).codes.astype(int)
    Yw = _within_center(Y, cells)
    for m in range(M):  # remove the (pooled) time effect
        Yw[tcodes == m] -= Yw[tcodes == m].mean(0)
    upper_levels = list(pd.unique(cells["upper"]))
    profile_rows, upper_groups = [], []
    for upper in upper_levels:
        rows = cells["upper"].to_numpy() == upper
        profile = []
        for m in range(M):
            ix = np.flatnonzero(rows & (tcodes == m))
            if len(ix) == 0:
                return None
            profile.append(Yw[ix].mean(0))
        profile_rows.append(np.concatenate(profile))
        upper_groups.append(cells.loc[rows, "group"].iloc[0])
    # Remove the pooled time profile so this tests group differences in trajectories,
    # not a common time shift. Permutations occur at the treatment-assignment unit.
    profiles = np.asarray(profile_rows).reshape(len(upper_levels), M, -1)
    profiles -= profiles.mean(axis=0, keepdims=True)
    profiles = profiles.reshape(len(upper_levels), -1)
    gcodes, guniq = pd.factorize(pd.Series(upper_groups))
    g, n_upper = len(guniq), len(upper_levels)
    dfn, dfd = g - 1, n_upper - g
    if dfn <= 0 or dfd <= 0:
        return None
    G = _gram(profiles)
    sst = float(np.trace(G))
    ssb0 = float(_ss_between(G, gcodes[None, :], g)[0])

    def F(ssb):
        ssw = sst - ssb
        return np.where(ssw > TOL * max(sst, 1), (ssb / dfn) / (np.where(ssw > 0, ssw, 1) / dfd), np.inf)

    nulls = []
    for b in _chunks(nperm, 500):
        gp = _permuted_labels(gcodes, b, rng)
        nulls.append(F(_ss_between(G, gp, g)))
    null = np.concatenate(nulls)
    F0 = float(F(np.array([ssb0]))[0])
    return dict(n=len(cells), n_units=n_upper, pseudo_F=F0, R2=ssb0 / sst if sst > 0 else np.nan, df_num=dfn, df_den=dfd,
                permutations=nperm, exact=False, min_attainable_p=1.0 / (nperm + 1),
                p_value=_p_value(F0, null, False))


def _kruskal_batch(values, L, ng):
    N = len(values)
    r = rankdata(values)
    H = np.zeros(L.shape[0])
    for k in range(ng):
        Z = (L == k).astype(float)
        cnt = Z.sum(1)
        s = Z @ r
        H += np.divide(s ** 2, cnt, out=np.zeros_like(s), where=cnt > 0)
    return 12.0 / (N * (N + 1)) * H - 3 * (N + 1)


def change_magnitude(cells, Y, tl, nperm, rng):
    """Magnitude of the within-unit change between two time levels, compared across groups (Kruskal H)."""
    if len(tl) != 2:
        return None
    rep_rows = _rep_matrix(cells)
    dist = np.linalg.norm(Y[rep_rows[:, 0]] - Y[rep_rows[:, 1]], axis=1)
    reps = list(cells.groupby("_rep", sort=False).groups.keys())
    rep_info = cells.groupby("_rep", sort=False)[["upper", "group"]].first().reindex(reps)
    ucodes, uuniq = pd.factorize(rep_info["upper"])
    g_u = rep_info.groupby("upper", sort=False)["group"].first().reindex(uuniq).to_numpy()
    gcodes_u, guniq = pd.factorize(g_u)

    def run(sel_groups):
        sel = np.isin(g_u[ucodes], sel_groups)
        vals = dist[sel]
        if sel.sum() < 3:
            return None
        gu_mask = np.isin(g_u, sel_groups)
        lab_u = pd.factorize(g_u[gu_mask])[0]
        map_u = -np.ones(len(g_u), int); map_u[gu_mask] = np.arange(gu_mask.sum())
        u_of_row = map_u[ucodes[sel]]
        ng = len(np.unique(lab_u))
        if ng < 2:
            return None
        H0 = float(_kruskal_batch(vals, lab_u[u_of_row][None, :], ng)[0])
        nulls = []
        for b in _chunks(nperm, 2000):
            lp = rng.permuted(np.tile(lab_u, (b, 1)), axis=1)[:, u_of_row]
            nulls.append(_kruskal_batch(vals, lp, ng))
        return dict(n_units=int(sel.sum()), H=H0, p_value=_p_value(H0, np.concatenate(nulls), False))

    allg = list(guniq)
    glob = run(allg)
    pair_rows = []
    for a, b in itertools.combinations(allg, 2):
        r = run([a, b])
        if r:
            meds = {k: float(np.median(dist[g_u[ucodes] == k])) for k in (a, b)}
            pair_rows.append(dict(group_a=a, group_b=b, median_change_a=meds[a], median_change_b=meds[b], **r))
    per = pd.DataFrame({"repeated_unit": reps, "upper": rep_info["upper"].values,
                        "group": rep_info["group"].values, "change_distance": dist})
    return glob, pd.DataFrame(pair_rows), per


# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------

def _holm(p: np.ndarray, method: str) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method=method)[1]
    return out


def run_beta(D: pd.DataFrame, frame, fixed: FixedDesign, nperm: int, rng, adjust: str, notes, scope: str,
             include_change_magnitude: bool = False) -> dict:
    ids = [i for i in frame.df.index if i in D.index]
    df = frame.df.loc[ids]
    Dm = D.loc[ids, ids].to_numpy(float)
    Dm, correction, original_negative_fraction = lingoes_correct(Dm)
    Y, info = gower_coords(Dm)
    if Y.shape[1] == 0:
        notes.warn(scope, "All samples are identical under this distance (no variation); tests skipped.")
        return {"terms": pd.DataFrame([{"term": "all", "note": "not run: zero distance variation"}])}
    if original_negative_fraction > 0.10:
        notes.info(scope, f"Lingoes correction applied to the non-Euclidean distance (c={correction:.6g}; "
                          f"original negative eigenvalue mass {original_negative_fraction:.0%}). Tests use all corrected axes.")
    out: dict[str, pd.DataFrame] = {}
    has_time = frame.has_time
    fx = FixedDesign(frame.subset(ids)) if len(ids) != len(frame.df) else fixed
    repeated = frame.has_units and df.groupby("upper").size().max() > 1
    pair_adj = adjust

    if not repeated:
        out["terms"] = independent_term_tests(Y, fx, nperm, rng)
        out["design_note"] = pd.DataFrame([{"note": "independent samples" if not frame.has_units else
                                            "units have a single sample each; treated as independent"}])
        strata = [("all", np.ones(len(df), bool))] if not has_time else [(t, (df["time"] == t).to_numpy()) for t in frame.time_levels]
        gw, pw, pdisp = [], [], []
        for label, mask in strata:
            sub = df[mask]
            Ys = Y[mask]
            _group_stratum(Ys, sub["group"].to_numpy(), label, nperm, rng, gw, pw, pdisp, True)
        _finish(out, gw, pw, pdisp, pair_adj)
        _record_distance_correction(out, correction, original_negative_fraction)
        return out

    ud = UnitData(Y, df, has_time, frame.time_levels)
    if not ud.ok:
        notes.warn(scope, f"Restricted-permutation beta tests skipped: {ud.reason}. Use time-stratified tests on unit "
                          "means is not valid here; rely on alpha/taxon models.")
        out["terms"] = pd.DataFrame([{"term": "all", "note": f"not run: {ud.reason}"}])
        return out
    if getattr(ud, "n_replicate_cells", 0):
        notes.info(scope, f"{ud.n_replicate_cells} unit x time cell(s) had replicate samples; replicates were averaged "
                          "at upper-unit x time cells for restricted-permutation tests.")
    cells_c, Yc = ud.complete()
    design_txt = "repeated/clustered units; restricted permutation (split-plot on principal coordinates)"
    rows = []
    if has_time and ud.n_incomplete:
        notes.warn(scope, f"{ud.n_incomplete} repeated unit(s) lack a time level; time and group:time beta tests use "
                          f"the {cells_c['_rep'].nunique()} complete units only.")
    # group (between units): complete units if time exists, else all
    gc, gY = (cells_c, Yc) if has_time else (ud.cells, ud.Y)
    r = group_test(gc, gY, nperm, rng)
    rows.append({**(r or {}), "term": "group", "design": design_txt + "; unit means, labels permuted among units",
                 "note": "" if r else "not estimable"})
    if has_time and cells_c["_rep"].nunique() >= 3:
        r = time_test(cells_c, Yc, frame.time_levels, nperm, rng)
        rows.append({**(r or {}), "term": "time", "design": design_txt + "; time labels permuted within units",
                     "note": "" if r else "not estimable"})
        r = interaction_test(cells_c, Yc, frame.time_levels, nperm, rng)
        rows.append({**(r or {}), "term": "group:time",
                     "design": design_txt + "; time effect removed, group labels permuted among units",
                     "note": "" if r else "not estimable"})
    out["terms"] = pd.DataFrame(rows)
    if frame.cov_kinds:
        notes.warn(scope, "Covariates are not used in restricted-permutation (repeated-unit) beta tests.")

    gw, pw, pdisp = [], [], []
    if has_time:
        for t in frame.time_levels:
            cs, Ys = ud.at_time(t)
            Yu, g, _ = _unit_means(cs, Ys)
            _group_stratum(Yu, g, t, nperm, rng, gw, pw, pdisp, False)
        # time within each group + pooled
        tw = []
        for gname in frame.group_levels + ["ALL (pooled)"]:
            cs = cells_c if gname.startswith("ALL") else cells_c[cells_c["group"] == gname]
            Ys = Yc if gname.startswith("ALL") else Yc[(cells_c["group"] == gname).to_numpy()]
            cs = cs.reset_index(drop=True)
            if cs["_rep"].nunique() < 3:
                continue
            r = time_test(cs, Ys, frame.time_levels, nperm, rng)
            if r:
                tw.append(dict(group=gname, **r))
        twdf = pd.DataFrame(tw)
        if not twdf.empty:
            mask = ~twdf["group"].str.startswith("ALL")
            twdf["p_adj"] = np.nan
            twdf.loc[mask, "p_adj"] = _holm(twdf.loc[mask, "p_value"].to_numpy(float), pair_adj)
            twdf["note"] = "labels swapped within units; p_adj covers per-group rows only; min_attainable_p shows test power limit"
        out["time_within_group"] = twdf
        cm = change_magnitude(cells_c, Yc, frame.time_levels, nperm, rng) if include_change_magnitude else None
        if cm:
            glob, pair, per = cm
            out["change_magnitude"] = pd.DataFrame([dict(
                test="group effect on magnitude of within-unit change (Kruskal H, labels permuted among units)",
                **glob, note="magnitude only; direction is tested by the group:time term")]) if glob else pd.DataFrame()
            if not pair.empty:
                pair["p_adj"] = _holm(pair["p_value"].to_numpy(float), pair_adj)
            out["change_magnitude_pairwise"] = pair
            out["change_magnitude_units"] = per
    else:
        Yu, g, _ = _unit_means(ud.cells, ud.Y)
        _group_stratum(Yu, g, "all", nperm, rng, gw, pw, pdisp, False)
    _finish(out, gw, pw, pdisp, pair_adj)
    _record_distance_correction(out, correction, original_negative_fraction)
    return out


def _record_distance_correction(out, constant, negative_fraction):
    for key, df in list(out.items()):
        if isinstance(df, pd.DataFrame) and not df.empty:
            df["lingoes_constant"] = constant
            df["original_negative_eigenvalue_fraction"] = negative_fraction


def _group_stratum(Yu, labels, stratum, nperm, rng, gw, pw, pdisp, adjust_cov_note):
    r = oneway_perm(Yu, np.asarray(labels), nperm, rng)
    if r:
        gw.append(dict(stratum=stratum, **r))
    d = permdisp(Yu, np.asarray(labels), nperm, rng)
    if d:
        pdisp.append(dict(stratum=stratum, **d))
    labs = np.asarray(labels)
    for a, b in itertools.combinations(list(pd.unique(labs)), 2):
        m = np.isin(labs, [a, b])
        rr = oneway_perm(Yu[m], labs[m], nperm, rng)
        if rr:
            pw.append(dict(stratum=stratum, group_a=a, group_b=b, **rr))


def _finish(out, gw, pw, pdisp, adjust):
    gwdf = pd.DataFrame(gw)
    pwdf = pd.DataFrame(pw)
    pddf = pd.DataFrame(pdisp)
    if not gwdf.empty and not pddf.empty:
        gwdf = gwdf.merge(pddf, on="stratum", how="left")
        interpretations = []
        for p, q in zip(gwdf["p_value"], gwdf["dispersion_p"]):
            if not np.isfinite(p):
                text = "PERMANOVA unavailable"
            elif p < 0.05 and np.isfinite(q) and q < 0.05:
                text = "significant PERMANOVA and PERMDISP; centroid and dispersion effects cannot be separated"
            elif p < 0.05 and np.isfinite(q):
                text = "significant PERMANOVA; no dispersion difference detected (not proof of equal dispersion)"
            elif p < 0.05:
                text = "significant PERMANOVA; dispersion test unavailable"
            else:
                text = "no significant PERMANOVA group difference"
            interpretations.append(text)
        gwdf["interpretation"] = interpretations
    if not gwdf.empty:
        gwdf["p_adj_across_strata"] = _holm(gwdf["p_value"].to_numpy(float), adjust)
    if not pwdf.empty:
        pwdf["p_adj"] = np.nan
        for _, idx in pwdf.groupby("stratum").groups.items():
            pwdf.loc[idx, "p_adj"] = _holm(pwdf.loc[idx, "p_value"].to_numpy(float), adjust)
    out["group_within_stratum"] = gwdf
    out["pairwise_group"] = pwdf
