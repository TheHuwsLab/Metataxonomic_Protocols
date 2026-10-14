import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

GENERA = [("Bacteria", "Firmicutes", "Clostridia", "Lachnospirales", "Lachnospiraceae", f"Gen{i}") for i in range(12)] + \
         [("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Prevotellaceae", f"Prev{i}") for i in range(8)]


def lineage(i, style="gtdb"):
    k, p, c, o, f, g = GENERA[i % len(GENERA)]
    if style == "gtdb":
        return f"d__{k};p__{p};c__{c};o__{o};f__{f};g__{g};s__{g} sp{i}"
    return f"{k};{p};{c};{o};{f};{g};{g} sp{i}"


def make_counts(meta, n_tax=40, effect=None, seed=0, depth=3000, unit_sd=0.4):
    rng = np.random.default_rng(seed)
    base = rng.gamma(0.5, 2, n_tax) + 0.05
    unit = {u: rng.normal(0, unit_sd, n_tax) for u in meta["unit"].unique()}
    rows = {}
    for sid, r in meta.iterrows():
        lam = np.log(base) + unit[r["unit"]]
        if effect is not None:
            lam = lam + effect(r)
        p = np.exp(lam); p /= p.sum()
        rows[sid] = rng.multinomial(int(depth * rng.uniform(0.8, 1.2)), p)
    return pd.DataFrame(rows, index=[f"asv{i}" for i in range(n_tax)])


def rusitec_like(seed=0, effect=True):
    rows = []
    for v in range(16):
        for t in ("0h", "4h"):
            rows.append(dict(sample_id=f"V{v}_{t}", Vessel=f"V{v}", Treatment=["CTL", "BCP", "TAN", "SAP"][v % 4],
                             Timepoint=t, unit=f"V{v}"))
    meta = pd.DataFrame(rows).set_index("sample_id")
    eff = None
    if effect:
        def eff(r):
            e = np.zeros(40)
            if r["Timepoint"] == "4h":
                e[:5] += 1.2
            return e
    return meta, make_counts(meta, effect=eff, seed=seed)


def piglet_like(seed=0):
    rows = []
    for pen in range(12):
        for a in range(3):
            rows.append(dict(sample_id=f"P{pen}A{a}", Pen=f"Pen{pen}", Animal=f"A{a}", Treatment=["C", "T1", "T2"][pen % 3],
                             unit=f"Pen{pen}"))
    meta = pd.DataFrame(rows).set_index("sample_id")
    return meta, make_counts(meta, seed=seed, unit_sd=0.5)


@pytest.fixture
def rus(tmp_path):
    meta, counts = rusitec_like()
    mp, tp = tmp_path / "meta.csv", tmp_path / "table.csv"
    meta.drop(columns="unit").reset_index().to_csv(mp, index=False)
    t = counts.T
    t.columns = [f"asv{i}" for i in range(t.shape[1])]
    # lineage columns as in a GTDB/QIIME2 collapsed table
    t.columns = [lineage(i) + f"_{i}" if False else lineage(i) for i in range(t.shape[1])]
    t = t.T.groupby(level=0).sum().T
    t.index.name = "index"
    t.to_csv(tp)
    return tp, mp, tmp_path
