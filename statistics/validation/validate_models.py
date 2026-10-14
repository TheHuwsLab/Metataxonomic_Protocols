"""Benchmark custom REML fits against statsmodels MixedLM.

Run from statistics/:
    python validation/validate_models.py --replicates 100

This compares estimates and model-based standard errors for the same fixed and
random-intercept structures. Containment denominator degrees of freedom are a
separate approximation and are not expected to match statsmodels' asymptotic
tests.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from metatax_stats.common import Notes
from metatax_stats.design import DesignSpec, FixedDesign, analyse_rank, build_frame
from metatax_stats.models import fit_model


def _case(seed: int, nested: bool, unbalanced: bool, boundary: bool):
    rng = np.random.default_rng(seed)
    groups = (("Control", 8), ("Treatment", 5 if unbalanced else 8))
    rows = []
    for group, n_units in groups:
        for unit_i in range(n_units):
            unit = f"{group}_U{unit_i}"
            animals = range(3) if nested else [None]
            for animal_i in animals:
                animal = f"A{animal_i}" if animal_i is not None else None
                for time in ("T0", "T1"):
                    rows.append({
                        "sample_id": f"{unit}_{animal}_{time}",
                        "group": group,
                        "time": time,
                        "subject": animal if nested else unit,
                        "cluster": unit if nested else None,
                        "unit": unit,
                        "animal_key": f"{unit}_{animal}" if animal is not None else None,
                    })
    meta = pd.DataFrame(rows).set_index("sample_id")
    spec = (DesignSpec("group", "time", "subject", "cluster") if nested
            else DesignSpec("group", "time", "subject"))
    frame, _ = build_frame(meta, spec, Notes(verbose=False))
    fixed = FixedDesign(frame)
    rank = analyse_rank(fixed)
    upper_sd = 0.0 if boundary else 0.8
    lower_sd = 0.0 if boundary else 0.45
    upper_effect = {u: rng.normal(0, upper_sd) for u in frame.df["upper"].unique()}
    lower_effect = ({a: rng.normal(0, lower_sd) for a in meta["animal_key"].dropna().unique()}
                    if nested else {})
    beta = np.linspace(-0.2, 0.5, fixed.X.shape[1])
    mean = fixed.X.to_numpy() @ beta
    y = []
    for sid, r in frame.df.iterrows():
        value = mean[frame.df.index.get_loc(sid)] + upper_effect[r["upper"]]
        if nested:
            value += lower_effect[meta.loc[sid, "animal_key"]]
        y.append(value + rng.normal(0, 0.7))
    return frame, fixed, rank, np.asarray(y), nested, beta


def _fit_reference(frame, fixed, y, nested):
    import statsmodels.api as sm

    d = frame.df.copy()
    for j in range(fixed.X.shape[1]):
        d[f"x{j}"] = fixed.X.iloc[:, j].to_numpy()
    d["y"] = y
    formula = "y ~ 0 + " + " + ".join(f"x{j}" for j in range(fixed.X.shape[1]))
    kwargs = {"groups": d["upper"], "re_formula": "1"}
    if nested:
        kwargs["vc_formula"] = {"lower": "0 + C(lower)"}
    model = sm.MixedLM.from_formula(formula, d, **kwargs)
    failures = []
    converged_fits = []
    last_fit = None
    for method in ("lbfgs", "bfgs", "cg", "powell"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = model.fit(reml=True, method=method, disp=False, maxiter=1000)
            last_fit = (fit, method)
            if fit.converged:
                converged_fits.append(last_fit)
        except (np.linalg.LinAlgError, ValueError, RuntimeError) as exc:
            failures.append(f"{method}: {exc}")
    if converged_fits:
        return max(converged_fits, key=lambda item: item[0].llf)
    if last_fit is not None:
        return last_fit
    raise RuntimeError("; ".join(failures))


def validate(replicates: int, seed: int):
    scenarios = [
        ("repeated_balanced", False, False, False),
        ("repeated_unbalanced", False, True, False),
        ("nested_balanced", True, False, False),
        ("nested_unbalanced", True, True, False),
        ("boundary_variances", False, False, True),
    ]
    rows = []
    for scenario_index, (name, nested, unbalanced, boundary) in enumerate(scenarios):
        for replicate in range(replicates):
            frame, fixed, rank, y, is_nested, _ = _case(
                seed + scenario_index * 100_000 + replicate, nested, unbalanced, boundary
            )
            custom = fit_model(y, fixed, rank)
            row = {"scenario": name, "replicate": replicate, "custom_ok": custom.ok}
            if not custom.ok:
                row["reference_ok"] = False
                row["failure"] = custom.note
                rows.append(row)
                continue
            try:
                ref, ref_optimizer = _fit_reference(frame, fixed, y, is_nested)
            except (np.linalg.LinAlgError, ValueError, RuntimeError) as exc:
                row["reference_ok"] = False
                row["failure"] = str(exc)
                rows.append(row)
                continue
            row.update({
                "reference_ok": bool(ref.converged),
                "reference_optimizer": ref_optimizer,
                "max_abs_fixed_effect_difference": float(np.max(np.abs(custom.beta - ref.fe_params.to_numpy()))),
                "max_abs_fixed_se_difference": float(np.max(np.abs(
                    np.sqrt(np.diag(custom.cov)) - ref.bse_fe.to_numpy()
                ))),
                "custom_upper_variance": custom.variance.get("upper", 0.0),
                "reference_upper_variance": float(ref.cov_re.iloc[0, 0]),
                "custom_lower_variance": custom.variance.get("lower", 0.0),
                "reference_lower_variance": float(ref.vcomp[0]) if is_nested else np.nan,
                "custom_residual_variance": custom.sigma2,
                "reference_residual_variance": float(ref.scale),
                "custom_ddf_method": "containment approximation",
            })
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=8117)
    parser.add_argument("--out", type=Path, default=Path("validation/reml_benchmark.tsv"))
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("replicates must be positive")
    result = validate(args.replicates, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    valid = result[result["reference_ok"].fillna(False)]
    summary_path = args.out.with_suffix(".md")
    if len(valid):
        metrics = ["max_abs_fixed_effect_difference", "max_abs_fixed_se_difference"]
        stats = valid.groupby("scenario")[metrics].agg(["median", lambda x: x.quantile(0.95)])
        stats.columns = ["_".join((col, "p95" if stat.startswith("<lambda") else stat))
                         for col, stat in stats.columns]
        summary_text = "```\n" + stats.to_string() + "\n```"
    else:
        summary_text = "_No reference fits converged._"
    summary_path.write_text(
        "# REML benchmark summary\n\n"
        f"Seed: {args.seed}. Replicates per scenario: {args.replicates}. "
        f"Converged reference fits: {len(valid)} of {len(result)}.\n\n"
        "The TSV contains fit-level fixed-effect, standard-error and variance-component comparisons. "
        "Boundary or nonconverged fits are retained and should be reviewed rather than silently dropped. "
        "Containment denominator degrees of freedom remain an approximation and are not benchmarked against "
        "statsmodels' asymptotic tests.\n\n"
        "Median and 95th-percentile absolute differences for converged fits:\n\n"
        f"{summary_text}\n"
    )
    print(f"Wrote {len(result)} benchmark fits ({len(valid)} converged references) to {args.out}")
    if len(valid):
        print(valid.groupby("scenario")[[
            "max_abs_fixed_effect_difference", "max_abs_fixed_se_difference"
        ]].median().to_string())


if __name__ == "__main__":
    main()
