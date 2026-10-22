#!/usr/bin/env python3
"""
Read-only review of the reference time cuts in CALIBRATION/set_reftimes/reftime_qa.

Nothing in PARAM/ or DBASE/ is written. The configuration definition is the SAME as in
parameter_generator.py (special_change_runs and manual_groups below must be kept in sync
with that file). The report answers three questions.

 1. How much do the cuts scatter inside one current Configuration (the noise floor)?
 2. How much do they shift between Configurations of the same kinematics setting?
    (If small, those Configurations could be merged into one.)
 3. How much do they shift from one kinematics setting to the next?

Sign convention (see _REFTIME_HEADER in parameter_generator.py): the first reference time
greater than abs(refcut) is used. A lower abs(refcut) never rejects the true reference time,
so the conservative cut of a group is the one with the SMALLEST abs value, whatever sign the
files store. All spreads are computed on abs values, in module channels.

Usage
  python parameter_report.py
  python parameter_report.py --tol-tdc 50 --tol-adc 50 --out-dir param_report
"""
import argparse
import os
import re
from glob import glob

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--rcdb", default="rcdb_draft.csv")
ap.add_argument("--reftime-dir", default="CALIBRATION/set_reftimes/reftime_qa")
ap.add_argument("--tol-tdc", type=float, default=100.0, help="tolerance for *_tdcrefcut spread (channels)")
ap.add_argument("--tol-adc", type=float, default=100.0, help="tolerance for *_adcrefcut spread (channels)")
ap.add_argument("--out-dir", default="param_report")
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration definition (copied from parameter_generator.py, keep in sync)
# ---------------------------------------------------------------------------
run_db = pd.read_csv(args.rcdb)
config_cols = ["Experiment", "Kinematics Setting", "Target"]

special_change_runs = [26085, 26088, 26132, 26148, 26155, 26165, 26161, 26190, 26271,
                       26284, 26302, 26361, 26448, 26489, 26593, 26612, 26613, 26641,
                       26652, 26769, 26724, 26784]

run_db = run_db.sort_values(config_cols + ["Run Number"]).reset_index(drop=True)
run_db["is_change_point"] = run_db["Run Number"].isin(special_change_runs)
run_db["setting_group"] = run_db.groupby(config_cols)["is_change_point"].cumsum()
run_db["Configuration"] = run_db.groupby(config_cols + ["setting_group"]).ngroup()

manual_groups = [
    [26148, 26149, 26150],
    [26151, 26152, 26153, 26154],
    [26155, 26156, 26157],
]
run_db.loc[run_db.Configuration == 1, "Configuration"] = 0
run_db["kin_sub"] = 0  # the manual groups are different kinematics under one label
for target_value, group in enumerate(manual_groups):
    mask = run_db["Run Number"].isin(group)
    run_db.loc[mask, "Configuration"] = target_value + 1
    run_db.loc[mask, "kin_sub"] = target_value + 1

run_db = run_db.loc[(run_db["Run Number"] != 26084) & (run_db["Run Status"] != "Junk") &
                    (run_db["Run Status"] != "Short") & (run_db["Target"] != "Carbon Hole") &
                    (run_db["Kinematics Setting"] != "Beam Checkout") & (run_db["Target"] != "Home"), :]
run_db = run_db.sort_values(by="Run Number")
run_db["Configuration"] = pd.factorize(run_db["Configuration"])[0]

run_db["setting"] = (run_db["Experiment"] + " " + run_db["Kinematics Setting"] +
                     np.where(run_db["kin_sub"] > 0, " [" + run_db["kin_sub"].astype(str) + "]", ""))
cfg_tbl = (run_db.groupby("Configuration")
           .agg(setting=("setting", "first"), target=("Target", lambda s: ",".join(sorted(set(s)))),
                first_run=("Run Number", "min"), last_run=("Run Number", "max"))
           .sort_index())

# ---------------------------------------------------------------------------
# Read the cuts
# ---------------------------------------------------------------------------
CUT_PARAMS = [
    "hdc_tdcrefcut", "hhodo_tdcrefcut", "hhodo_adcrefcut", "hcer_adcrefcut", "hcal_adcrefcut",
    "pdc_tdcrefcut", "phodo_tdcrefcut", "phodo_adcrefcut", "pngcer_adcrefcut", "phgcer_adcrefcut",
    "paero_adcrefcut", "pcal_adcrefcut", "t_coin_trig_tdcrefcut", "t_coin_trig_adcrefcut",
]
def tol_for(p):
    return args.tol_tdc if p.endswith("_tdcrefcut") else args.tol_adc

def parse_scalars(path):
    values = {}
    for line in open(path).read().splitlines():
        line = line.split(";", 1)[0].strip()
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([-+]?\d+\.?\d*)\s*$", line)
        if m:
            values[m.group(1)] = float(m.group(2))
    return values

d = args.reftime_dir
rows = []
for f in sorted(glob(os.path.join(d, "tcoin_*.param"))):
    run = int(re.search(r"tcoin_(\d+)\.param$", f).group(1))
    scalars = parse_scalars(f)
    for other in (f"h_reftime_cut_coindaq_{run}.param", f"p_reftime_cut_{run}.param"):
        path = os.path.join(d, other)
        if os.path.exists(path):
            scalars.update(parse_scalars(path))
        else:
            print(f"WARNING: missing {path}")
    rows.append({"run": run, **{p: scalars.get(p, np.nan) for p in CUT_PARAMS}})
if not rows:
    raise SystemExit(f"No tcoin_*.param found in {d}")
per_run = pd.DataFrame(rows).set_index("run").sort_index()

known = per_run.index.isin(run_db["Run Number"])
if (~known).any():
    print("WARNING: surveyed runs not in the run database (ignored):", list(per_run.index[~known]))
per_run = per_run[known]
rmap = run_db.set_index("Run Number")
per_run.insert(0, "configuration", rmap.loc[per_run.index, "Configuration"].to_numpy())
per_run.insert(1, "setting", rmap.loc[per_run.index, "setting"].to_numpy())
per_run.to_csv(os.path.join(args.out_dir, "report_per_run.csv"))  # raw values

# A reference time cut of exactly 0 is not a valid cut (typically a broken or empty QA file).
# It is excluded from every spread below and listed here. parameter_generator.py has no such
# check, so a surviving 0 would win the "smallest abs value" merge and be written to the output.
zero_cuts = [{"run": run, "configuration": row["configuration"], "setting": row["setting"], "param": p}
             for run, row in per_run.iterrows() for p in CUT_PARAMS if row[p] == 0]
zero_cuts = pd.DataFrame(zero_cuts, columns=["run", "configuration", "setting", "param"])
zero_cuts.to_csv(os.path.join(args.out_dir, "report_zero_cuts.csv"), index=False)
per_run[CUT_PARAMS] = per_run[CUT_PARAMS].mask(per_run[CUT_PARAMS] == 0)

def abs_spread(s):
    a = s.dropna().abs()
    return float(a.max() - a.min()) if len(a) else np.nan

def conservative(s):
    s = s.dropna()
    return float(s.iloc[np.argmin(np.abs(s.to_numpy()))]) if len(s) else np.nan

# ---------------------------------------------------------------------------
# 1. Inside each current Configuration
# ---------------------------------------------------------------------------
w_cfg = []
for cfg, g in per_run.groupby("configuration"):
    for p in CUT_PARAMS:
        s = g[p].dropna()
        if s.empty:
            continue
        w_cfg.append({"configuration": cfg, "setting": cfg_tbl.loc[cfg, "setting"], "param": p,
                      "n_runs": len(s), "abs_min": s.abs().min(), "abs_max": s.abs().max(),
                      "spread": abs_spread(s)})
w_cfg = pd.DataFrame(w_cfg)
w_cfg.to_csv(os.path.join(args.out_dir, "report_within_configuration.csv"), index=False)

# ---------------------------------------------------------------------------
# 2. Between Configurations of the same setting
# ---------------------------------------------------------------------------
w_set = []
for setting, g in per_run.groupby("setting", sort=False):
    cfgs = sorted(g["configuration"].unique())
    for p in CUT_PARAMS:
        s = g[p].dropna()
        if s.empty:
            continue
        per_cfg = g.groupby("configuration")[p].apply(conservative).dropna()
        w_set.append({"setting": setting, "param": p, "n_configurations": len(cfgs),
                      "n_runs": len(s), "spread_all_runs": abs_spread(s),
                      "spread_between_configs": abs_spread(per_cfg) if len(per_cfg) > 1 else 0.0,
                      "tolerance": tol_for(p)})
w_set = pd.DataFrame(w_set)
w_set["exceeds_tol"] = w_set["spread_all_runs"] > w_set["tolerance"]
w_set.to_csv(os.path.join(args.out_dir, "report_within_setting.csv"), index=False)

# ---------------------------------------------------------------------------
# 3. Between settings (conservative cut per setting, shift to the previous setting)
# ---------------------------------------------------------------------------
first_run = run_db.groupby("setting")["Run Number"].min()
vals = pd.DataFrame({st: {p: conservative(g[p]) for p in CUT_PARAMS}
                     for st, g in per_run.groupby("setting", sort=False)}).T
vals.insert(0, "first_run", first_run.reindex(vals.index).to_numpy())
vals = vals.sort_values("first_run")
shift = vals[CUT_PARAMS].abs().diff().add_prefix("shift_prev_")
between = pd.concat([vals.add_prefix("value_").rename(columns={"value_first_run": "first_run"}), shift], axis=1)
between.index.name = "setting"
between.to_csv(os.path.join(args.out_dir, "report_between_settings.csv"))

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
if len(zero_cuts):
    print("Cuts equal to 0 (excluded from all spreads below, fix or drop these runs):")
    print(zero_cuts.groupby(["run", "configuration", "setting"]).param.apply(", ".join).to_string())
    print()

neg = {p: round(float((per_run[p].dropna() < 0).mean()), 2) for p in CUT_PARAMS}
print("Fraction of negative stored cuts per parameter (1.0 = all negative):")
print("  ", neg)

print(f"\nSurveyed runs {len(per_run)}, configurations with a survey {per_run.configuration.nunique()} of {len(cfg_tbl)}")
missing = cfg_tbl.loc[~cfg_tbl.index.isin(per_run.configuration.unique())]
if len(missing):
    print("Configurations without any survey (they need a donor when the database is generated):")
    print(missing[["setting", "target", "first_run", "last_run"]].to_string())

multi = w_cfg[w_cfg.n_runs > 1]
if len(multi):
    q = multi["spread"].quantile([0.5, 0.9, 1.0])
    print("\nNoise floor, spread inside one current Configuration (channels), median / 90% / max:",
          f"{q.iloc[0]:.0f} / {q.iloc[1]:.0f} / {q.iloc[2]:.0f}")
    print("A tolerance well above this noise floor and well below the real shifts is a sensible choice.")

bad = w_set[w_set.exceeds_tol]
n_set = w_set.setting.nunique()
print(f"\nSettings with all cuts inside tolerance (tdc {args.tol_tdc:g}, adc {args.tol_adc:g}): {n_set - bad.setting.nunique()} of {n_set}")
if len(bad):
    worst = (bad.sort_values("spread_all_runs", ascending=False)
             .groupby("setting", sort=False).first()[["param", "spread_all_runs", "n_configurations"]])
    print("Settings above tolerance (worst parameter shown):")
    print(worst.to_string())

multi_cfg = w_set[w_set.n_configurations > 1].groupby("setting").agg(
    n_configurations=("n_configurations", "first"),
    max_spread_all_runs=("spread_all_runs", "max"),
    max_between_configs=("spread_between_configs", "max"))
if len(multi_cfg):
    print("\nSettings that currently have more than one Configuration")
    print(multi_cfg.assign(mergeable=lambda x: [
        not (w_set[(w_set.setting == s)].exceeds_tol.any()) for s in x.index]).to_string())

print(f"\nLargest jumps between consecutive settings (channels):")
jump = shift.max(axis=1).sort_values(ascending=False).head(8)
print(pd.DataFrame({"max_shift_prev": jump, "worst_param": shift.loc[jump.index].idxmax(axis=1).str.replace("shift_prev_", "")}).to_string())
print(f"\nCSV files written to {args.out_dir}/")