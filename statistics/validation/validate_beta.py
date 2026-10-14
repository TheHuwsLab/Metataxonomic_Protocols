"""Simulation study for the custom repeated/clustered beta permutation tests.

Run from statistics/:
    python validation/validate_beta.py --replicates 200 --permutations 999

The script reports estimated rejection rates and Wilson intervals; it does not
apply a pass/fail cutoff or imply validation outside the simulated designs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from metatax_stats.beta import (
    UnitData, group_test, independent_term_tests, interaction_test, time_test,
)
from metatax_stats.common import Notes
from metatax_stats.design import DesignSpec, FixedDesign, build_frame


def _simulate(n_per_group: int, nested_subjects: int, effect: str, seed: int):
    rng = np.random.default_rng(seed)
    groups, times = ("Control", "Treatment"), ("T0", "T1")
    rows, values = [], []
    dim = 6
    for gi, group in enumerate(groups):
        for unit_i in range(n_per_group):
            unit = f"{group}_U{unit_i}"
            unit_effect = rng.normal(0, 0.7, dim)
            subjects = range(nested_subjects) if nested_subjects > 0 else [None]
            for subject_i in subjects:
                subject = f"S{subject_i}" if subject_i is not None else None
                subject_effect = rng.normal(0, 0.35, dim) if subject is not None else 0
                for ti, time in enumerate(times):
                    sid = f"{unit}_{subject}_{time}" if subject is not None else f"{unit}_{time}"
                    rows.append({
                        "sample_id": sid,
                        "group": group,
                        "time": time,
                        "unit": unit if nested_subjects == 0 else None,
                        "cluster": unit if nested_subjects > 0 else None,
                        "subject": subject,
                    })
                    x = unit_effect + subject_effect + rng.normal(0, 1, dim)
                    if effect == "group" and gi == 1:
                        x[0] += 0.9
                    if effect == "time" and ti == 1:
                        x[1] += 0.9
                    if effect == "interaction" and gi == 1 and ti == 1:
                        x[2] += 1.0
                    values.append(x)
    meta = pd.DataFrame(rows).set_index("sample_id")
    spec = (DesignSpec("group", "time", "subject", "cluster") if nested_subjects > 0
            else DesignSpec("group", "time", "unit"))
    frame, _ = build_frame(meta, spec, Notes(verbose=False))
    data = UnitData(np.asarray(values), frame.df, True, list(times))
    cells, cell_values = data.complete()
    return cells, cell_values


def _simulate_independent(n_per_group: int, effect: str, seed: int):
    rng = np.random.default_rng(seed)
    rows, values = [], []
    for group_i, group in enumerate(("Control", "Treatment")):
        for time_i, time in enumerate(("T0", "T1")):
            n_cell = n_per_group // 2 if time_i == 0 else n_per_group - n_per_group // 2
            for replicate in range(n_cell):
                rows.append({
                    "sample_id": f"{group}_{time}_{replicate}",
                    "group": group,
                    "time": time,
                })
                x = rng.normal(0, 1, 6)
                if effect == "group" and group_i == 1:
                    x[0] += 0.9
                if effect == "time" and time_i == 1:
                    x[1] += 0.9
                if effect == "interaction" and group_i == 1 and time_i == 1:
                    x[2] += 1.0
                values.append(x)
    meta = pd.DataFrame(rows).set_index("sample_id")
    frame, _ = build_frame(meta, DesignSpec("group", "time"), Notes(verbose=False))
    return FixedDesign(frame), np.asarray(values)


def _wilson(successes: int, trials: int, confidence: float = 0.95):
    if not trials:
        return np.nan, np.nan
    z = float(norm.ppf(0.5 + confidence / 2))
    p = successes / trials
    den = 1 + z * z / trials
    centre = (p + z * z / (2 * trials)) / den
    half = z * np.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / den
    return centre - half, centre + half


def validate(replicates: int, permutations: int, unit_counts: list[int], nested_subjects: int, seed: int):
    effects = ("no_effect", "group", "time", "interaction")
    terms = ("group", "time", "group:time")
    rows = []
    rng = np.random.default_rng(seed)
    index = 0
    for n_units in unit_counts:
        for design, nested in (("independent", 0), ("repeated", 0), ("nested", nested_subjects)):
            for effect in effects:
                pvals = {term: [] for term in terms}
                for _ in range(replicates):
                    if design == "independent":
                        fixed, Y = _simulate_independent(n_units, effect, seed + index)
                    else:
                        cells, Y = _simulate(n_units, nested if design == "nested" else 0,
                                             effect, seed + index)
                    index += 1
                    if design == "independent":
                        terms_frame = independent_term_tests(Y, fixed, permutations, rng)
                        tests = terms_frame.set_index("term").to_dict("index")
                    else:
                        tests = {
                            "group": group_test(cells, Y, permutations, rng),
                            "time": time_test(cells, Y, ["T0", "T1"], permutations, rng),
                            "group:time": interaction_test(cells, Y, ["T0", "T1"], permutations, rng),
                        }
                    for term, result in tests.items():
                        if result is not None and pd.notna(result.get("p_value", np.nan)):
                            pvals[term].append(float(result["p_value"]))
                for term, values in pvals.items():
                    sig = int(np.sum(np.asarray(values) <= 0.05))
                    lo, hi = _wilson(sig, len(values))
                    rows.append({
                        "units_per_group": n_units,
                        "design": design,
                        "nested_subjects_per_unit": nested,
                        "effect_scenario": effect,
                        "tested_term": term,
                        "replicates": len(values),
                        "permutations": permutations,
                        "rejection_rate_at_0.05": sig / len(values) if values else np.nan,
                        "wilson_95_low": lo,
                        "wilson_95_high": hi,
                    })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--permutations", type=int, default=199)
    parser.add_argument("--unit-counts", default="4,8,12,20")
    parser.add_argument("--nested-subjects", type=int, default=3)
    parser.add_argument("--seed", type=int, default=7319)
    parser.add_argument("--out", type=Path, default=Path("validation/beta_simulation_results.tsv"))
    args = parser.parse_args()
    if args.replicates < 1 or args.permutations < 1 or args.nested_subjects < 1:
        parser.error("replicates, permutations and nested-subjects must be positive")
    units = [int(v) for v in args.unit_counts.split(",")]
    if not units or any(v < 2 for v in units):
        parser.error("unit-counts must be comma-separated integers >= 2")
    result = validate(args.replicates, args.permutations, units, args.nested_subjects, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    report_rows = []
    for (n_units, design, nested, effect), group in result.groupby(
            ["units_per_group", "design", "nested_subjects_per_unit", "effect_scenario"], sort=False):
        row = {"units/group": n_units, "design": design, "nested subjects/unit": nested, "effect": effect}
        for _, tested in group.iterrows():
            row[tested["tested_term"]] = (
                f"{tested['rejection_rate_at_0.05']:.3f} "
                f"[{tested['wilson_95_low']:.3f}, {tested['wilson_95_high']:.3f}]"
            )
        report_rows.append(row)
    report_table = pd.DataFrame(report_rows).to_string(index=False)
    summary_path = args.out.with_suffix(".md")
    summary_path.write_text(
        "# Repeated/clustered beta permutation simulation\n\n"
        f"Seed: {args.seed}. Simulations per scenario: {args.replicates}. "
        f"Monte Carlo permutations: {args.permutations}. Unit counts/group: {', '.join(map(str, units))}.\n\n"
        "The table gives the fraction of simulations with p ≤ 0.05 and its 95% Wilson interval. "
        "`no_effect` scenarios estimate type-I error; the other rows show power only for the named effect "
        "and also show rejection of non-target terms (which should be interpreted as empirical error). "
        "Independent designs use Freedman-Lane term tests; repeated and nested designs use the toolkit's "
        "restricted group/time/interaction procedures. "
        "These estimates are uncertain and apply only to these simulated structures; there is no automatic "
        "pass/fail cutoff or guarantee for other designs.\n\n"
        "```\n" + report_table + "\n```\n\n"
        f"Fit-level results are in `{args.out.name}`.\n"
    )
    print(f"Wrote simulation rejection-rate estimates and Wilson intervals to {args.out}")


if __name__ == "__main__":
    main()
