import pandas as pd
import numpy as np
import re
from glob import glob
import re
from datetime import datetime
import os
import math
from copy import copy
def build_config_row(group):
    kin_setting = group["Kinematics Setting"].iloc[0]
    is_production = bool(production_pattern.match(str(kin_setting).strip()))

    n_target = 1 if is_production else 4
    n = min(n_target, len(group))
    sampled = rng.choice(group["Run Number"].values, size=n, replace=False)
    sampled = np.sort(sampled).astype(int).tolist()

    return pd.Series({
        "Experiment": group["Experiment"].iloc[0],
        "Kinematics Setting": kin_setting,
        "Target": group["Target"].iloc[0],
        "sample_run_numbers": sampled,
    })


run_db = pd.read_csv("rcdb_draft.csv")
config_cols = ["Experiment", "Kinematics Setting", "Target"]

special_change_runs = [26085, 26088, 26132, 26148, 26155, 26165, 26161, 26190, 26271, 
                        26284, 26302, 26361, 26448, 26489, 26593, 26612, 26613, 26641, 
                        26652, 26769, 26724, 26784]  # See notes_for_some_runs.md

# Sort so cumsum respects run-number order within each config
run_db = run_db.sort_values(config_cols + ["Run Number"]).reset_index(drop=True)

# Mark the exact runs where a setting change was noted
run_db["is_change_point"] = run_db["Run Number"].isin(special_change_runs)

# Cumulative count of change points within each kinematic config
# -> increments every time a change point is hit, giving each "regime" its own id
run_db["setting_group"] = run_db.groupby(config_cols)["is_change_point"].cumsum()

# Runs with setting_group > 0 belong to a non-default (post-change) setting
run_db["flagged_changed_setting"] = run_db["setting_group"] > 0

# Unique integer label per (config_cols + setting_group) combination
run_db["Configuration"] = run_db.groupby(config_cols + ["setting_group"]).ngroup()

# Manual override: force these runs to share one Configuration,
# regardless of what the automatic setting_group logic produced
manual_groups = [
    [26148, 26149, 26150],
    [26151, 26152, 26153, 26154],
    [26155, 26156, 26157],
]
# manually sort out
run_db.loc[run_db.Configuration == 1, "Configuration"] = 0
for target_value, group in enumerate(manual_groups):
    mask = run_db["Run Number"].isin(group)
    run_db.loc[mask, "Configuration"] = target_value + 1

    
run_db = run_db.loc[ (run_db["Run Number"] != 26084) & (run_db["Run Status"]!="Junk") & (run_db["Run Status"]!="Short") & (run_db["Target"]!="Carbon Hole") & (run_db["Kinematics Setting"]!="Beam Checkout") & (run_db["Target"]!="Home"), :]
run_db.sort_values(by = 'Run Number', inplace = True)
run_db["Configuration"] = pd.factorize(run_db["Configuration"])[0]

# rng = np.random.default_rng(42)

# # Matches kinematics settings like "12a", "3b", "104e", etc.
# production_pattern = re.compile(r"^\d+[a-e]$")


# normalized_db = (
#     run_db.groupby("Configuration")
#     .apply(build_config_row, include_groups=False)
# )
# normalized_db.index.name = "Configuration"

# # Manual override for Configuration 52
# manual_sample_run_numbers = {
#     52: [26483, 26484, 26485, 26486, 26487, 26488],
# }

# for config_num, runs in manual_sample_run_numbers.items():
#     normalized_db.at[config_num, "sample_run_numbers"] = runs

# normalized_db.to_csv("calib_run_list_by_run_period.csv")
# all_runs = np.concatenate(normalized_db["sample_run_numbers"].to_numpy())
# np.savetxt("calib_run_list.csv", all_runs, fmt="%d", delimiter=",")

DEFAULT_TDC_NAMES = (    "h1X h1Y h2X h2Y h1T h2T hT1 hASUM hBSUM hCSUM hDSUM hPRLO hPRHI hSHWR hEDTM hCER hT2 "    "hDCREF1 hDCREF2 hDCREF3 hDCREF4 "    "hTRIG1_ROC1 hTRIG2_ROC1 hTRIG3_ROC1 hTRIG4_ROC1 hTRIG5_ROC1 hTRIG6_ROC1 "    "pTRIG1_ROC1 pTRIG2_ROC1 pTRIG3_ROC1 pTRIG4_ROC1 pTRIG5_ROC1 pTRIG6_ROC1 "    "pT1 pT2 p1X p1Y p2X p2Y p1T p2T pT3 pAER pHGCER pNGCER "    "pDCREF1 pDCREF2 pDCREF3 pDCREF4 pDCREF5 pDCREF6 pDCREF7 pDCREF8 pDCREF9 pDCREF10 "    "pEDTM pPRLO pPRHI "    "pTRIG1_ROC2 pTRIG2_ROC2 pTRIG3_ROC2 pTRIG4_ROC2 pTRIG5_ROC2 pTRIG6_ROC2 "    "hTRIG1_ROC2 hTRIG2_ROC2 hTRIG3_ROC2 hTRIG4_ROC2 hTRIG5_ROC2 hTRIG6_ROC2 "    "pSTOF_ROC2 pEL_LO_LO_ROC2 pEL_LO_ROC2 pEL_HI_ROC2 pEL_REAL_ROC2 pEL_CLEAN_ROC2 "    "hSTOF_ROC2 hEL_LO_LO_ROC2 hEL_LO_ROC2 hEL_HI_ROC2 hEL_REAL_ROC2 hEL_CLEAN_ROC2 "    "pSTOF_ROC1 pEL_LO_LO_ROC1 pEL_LO_ROC1 pEL_HI_ROC1 pEL_REAL_ROC1 pEL_CLEAN_ROC1 "    "hSTOF_ROC1 hEL_LO_LO_ROC1 hEL_LO_ROC1 hEL_HI_ROC1 hEL_REAL_ROC1 hEL_CLEAN_ROC1 "    "pPRE40_ROC1 pPRE100_ROC1 pPRE150_ROC1 pPRE200_ROC1 "    "hPRE40_ROC1 hPRE100_ROC1 hPRE150_ROC1 hPRE200_ROC1 "    "pPRE40_ROC2 pPRE100_ROC2 pPRE150_ROC2 pPRE200_ROC2 "    "hPRE40_ROC2 hPRE100_ROC2 hPRE150_ROC2 hPRE200_ROC2 "    "hDCREF5 hRF pRF hHODO_RF pHODO_RF").split()
ADC_CHANNELS = ["pFADC_TREF_ROC2", "hFADC_TREF_ROC1"]
PARAM_MAP = {
    **{f"pDCREF{i}": ["pdc_tdcrefcut"] for i in range(1, 11)},
    **{f"hDCREF{i}": ["hdc_tdcrefcut"] for i in range(1, 6)},
    "pT1": ["phodo_tdcrefcut"],
    "pT2": ["t_coin_trig_tdcrefcut"],
    "hT2": ["hhodo_tdcrefcut"],
    "pFADC_TREF_ROC2": ["phodo_adcrefcut", "pngcer_adcrefcut", "phgcer_adcrefcut",
                         "paero_adcrefcut", "pcal_adcrefcut", "t_coin_trig_adcrefcut"],
    "hFADC_TREF_ROC1": ["hhodo_adcrefcut", "hcer_adcrefcut", "hcal_adcrefcut"],
}
NEEDED_CHANNELS = list(PARAM_MAP)  # ordered list of the ~20 channels that feed a .param value
NEEDED_SET = set(NEEDED_CHANNELS)

TRIG_NAMES = "pTRIG1_ROC1 pTRIG4_ROC1 pTRIG1_ROC2 pTRIG4_ROC2"
ADC_NAMES_FULL = ("hASUM hBSUM hCSUM hDSUM hPSHWR hSHWR hAER hCER hFADC_TREF_ROC1 pAER "
                   "pHGCER pNGCER pPSHWR pFADC_TREF_ROC2 pHGCER_MOD pNGCER_MOD pHEL_NEG "
                   "pHEL_POS pHEL_MPS")
NUM_ADC = 19
TDC_OFFSET, ADC_TDC_OFFSET, TDC_CHANPERNS, EHADCOINTIME_OFFSET = 300.0, 200.0, 0.09766, 0.0

_REFTIME_HEADER = """; Cut to select the Reference time when multiple hits in reference time
; The units in channels for the module (CAEN tdc or FADC)
; negative value refcut means that the first reference time greater than the abs(refcut)
;     is used as reftime. If no ref time is found  greater than the abs(refcut) then first
;     reference time is used.
; positive value refcut means that the the first reference time greater than the abs(refcut)
;     is used as reftime. If no ref time is found  greater than the abs(refcut) then no
;     reference time is used and warning message is produced.
; Cut is on reference time per detector.
"""

def parse_run_num(tcoin_param):
  return int(tcoin_param.split('_')[-1].split('.')[0])

def run_to_hms_param(run):
  return 'CALIBRATION/set_reftimes/reftime_qa/h_reftime_cut_coindaq_{}.param'.format(run)

def run_to_shms_param(run):
  return 'CALIBRATION/set_reftimes/reftime_qa/p_reftime_cut_{}.param'.format(run)

def run_to_tcoin_param(run):
  return 'CALIBRATION/set_reftimes/reftime_qa/tcoin_{}.param'.format(run)

def run_to_configuration(run):
  return int(run_db.loc[run_db["Run Number"] == run, "Configuration"].to_numpy()[0])

tcoin_params = glob("CALIBRATION/set_reftimes/reftime_qa/tcoin_*.param")
run_numbers = np.sort(list(map(parse_run_num, tcoin_params)))
hms_params  = list(map(run_to_hms_param, run_numbers))
shms_params = list(map(run_to_shms_param, run_numbers))
tcoin_params = list(map(run_to_tcoin_param, run_numbers))
configuration_params = list(map(run_to_configuration, run_numbers))

def _extract_block(text, key):
    m = re.search(rf"{re.escape(key)}\s*=\s*(.*?)(?=\n\s*\n|\Z)", text, re.S)
    return m.group(1) if m else None

def parse_scalar_params(text):
    """key = value assignments only. Skips ';' comment lines, blank lines,
    and multi-line array assignments (e.g. t_coin_TdcTimeWindowMin)."""
    values = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(';'):
            continue
        line = line.split(';', 1)[0].strip()  # drop trailing inline comments
        m = re.match(r'^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([-+]?\d+\.?\d*)\s*$', line)
        if m:
            values[m.group(1)] = float(m.group(2))
    return values

def parse_tcoin_window(text):
    """Parse t_coin_tdcNames + t_coin_TdcTimeWindowMin/Max into
    {channel_name: (min, max)}. Only tcoin.param has these arrays."""
    names_block = _extract_block(text, "t_coin_tdcNames")
    names_match = re.search(r'"([^"]*)"', names_block) if names_block else None
    names = names_match.group(1).split() if names_match else []

    min_block = _extract_block(text, "t_coin_TdcTimeWindowMin")
    max_block = _extract_block(text, "t_coin_TdcTimeWindowMax")
    mins = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*", min_block)] if min_block else []
    maxs = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*", max_block)] if max_block else []

    if not (names and len(names) == len(mins) == len(maxs)):
        return {}
    return dict(zip(names, zip(mins, maxs)))

def parse_param_file(path):
    """Parse one Hall C .param file into:
        {"scalars": {param_name: value}, "window": {channel_name: (min, max)}}
    'window' is only non-empty for tcoin.param. Kept as a *separate*
    namespace from 'scalars' deliberately -- a channel name (e.g. 'hT2')
    and the param name it feeds (e.g. 'hhodo_tdcrefcut', via PARAM_MAP)
    are different quantities and would otherwise collide if merged into
    one flat dict.
    """
    text = open(path).read()
    return {"scalars": parse_scalar_params(text), "window": parse_tcoin_window(text)}

def survey_row(run, config, tcoin_path, hms_path, shms_path):
    """One row: {name}_min/{name}_max for every TDC channel (230 cols)
    + {name}_min for every ADC channel (2 cols) = 232 data columns, plus run.

    For a NEEDED channel, 'min' comes from its scalar cut (reversing
    PARAM_MAP's first/most-exclusive param -- same convention the app
    itself uses); otherwise 'min' falls back to the raw window array.
    'max' always comes straight from the window array -- there's no other
    source for it. No boilerplate/offset scalars (t_coin_numAdc,
    t_coin_tdcoffset, eHadCoinTime_Offset, etc.) are stored at all.
    """
    tcoin = parse_param_file(tcoin_path)
    hms = parse_param_file(hms_path)
    shms = parse_param_file(shms_path)
    scalars = {**tcoin["scalars"], **hms["scalars"], **shms["scalars"]}
    window = tcoin["window"]

    row = {"run": run, "configuration": config}
    for name in DEFAULT_TDC_NAMES:
        win = window.get(name, (0.0, 100000.0))
        lo = None
        if name in NEEDED_SET:
            p0 = PARAM_MAP[name][0]
            if p0 in scalars:
                lo = -scalars[p0]
        if lo is None:
            lo = win[0]
        row[f"{name}_min"] = lo
        row[f"{name}_max"] = win[1] if (win[1]>0) & (win[1]<=100000.0) else 100000.0
        assert row[f"{name}_max"] > row[f"{name}_min"]

    for name in ADC_CHANNELS:
        p0 = PARAM_MAP[name][0]
        row[f"{name}_min"] = -scalars[p0] if p0 in scalars else None

    return row




def format_array_floor(values, per_line=10, indent="\t\t\t\t  "):
    def fmt(v):
        return f"{v:.0f}" if float(v).is_integer() else f"{math.floor(v):d}"
    rows = []
    for i in range(0, len(values), per_line):
        rows.append(", ".join(fmt(v) for v in values[i:i + per_line]))
    return (",\n" + indent).join(rows)

def format_array_ceiling(values, per_line=10, indent="\t\t\t\t  "):
  def fmt(v):
      return f"{v:.0f}" if float(v).is_integer() else f"{math.ceil(v):d}"
  rows = []
  for i in range(0, len(values), per_line):
      rows.append(", ".join(fmt(v) for v in values[i:i + per_line]))
  return (",\n" + indent).join(rows)

def generate_tcoin_param(run, tdc_min, tdc_max, trig_tdcrefcut, trig_adcrefcut):
    tdcnames_str = " ".join(DEFAULT_TDC_NAMES)
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""; Copied from reftime_qa {run} on {now}

t_coin_numAdc = {NUM_ADC}
t_coin_numTdc = {len(DEFAULT_TDC_NAMES)}

t_coin_tdcoffset = {TDC_OFFSET}
t_coin_adc_tdc_offset = {ADC_TDC_OFFSET}

t_coin_tdcchanperns = {TDC_CHANPERNS}
eHadCoinTime_Offset = {EHADCOINTIME_OFFSET}

t_coin_trigNames="{TRIG_NAMES}"

; tdc cut is on pTRef2
; adc cut on pFADC_ROC2
t_coin_trig_tdcrefcut = {math.floor(trig_tdcrefcut):d}
t_coin_trig_adcrefcut = {math.floor(trig_adcrefcut):d}

t_coin_adcNames = "{ADC_NAMES_FULL}"

t_coin_tdcNames = "{tdcnames_str}"

t_coin_TdcTimeWindowMin = {format_array_floor(tdc_min)}

t_coin_TdcTimeWindowMax = {format_array_ceiling(tdc_max)}
"""

def generate_hms_param(run, hdc, hhodo_tdc, hhodo_adc, hcer_adc, hcal_adc):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return _REFTIME_HEADER + f"""
; Copied from reftime_qa {run} on {now}
hdc_tdcrefcut={math.floor(hdc):d}
hhodo_tdcrefcut={math.floor(hhodo_tdc):d}
hhodo_adcrefcut={math.floor(hhodo_adc):d}
hcer_adcrefcut={math.floor(hcer_adc):d}
hcal_adcrefcut={math.floor(hcal_adc):d}
"""

def generate_shms_param(run, pdc, phodo_tdc, phodo_adc, pngcer_adc, phgcer_adc, paero_adc, pcal_adc):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return _REFTIME_HEADER + f"""
; Copied from reftime_qa {run} on {now}
pdc_tdcrefcut={math.floor(pdc):d}
phodo_tdcrefcut={math.floor(phodo_tdc):d}
phodo_adcrefcut={math.floor(phodo_adc):d}
pngcer_adcrefcut={math.floor(pngcer_adc):d}
phgcer_adcrefcut={math.floor(phgcer_adc):d}
paero_adcrefcut={math.floor(paero_adc):d}
pcal_adcrefcut={math.floor(pcal_adc):d}
"""

def write_param_files_for_configuration(df, configuration, out_run=None, out_dir="."):
    """Generate tcoin_<label>.param, h_reftime_cut_coindaq_<label>.param, and
    p_reftime_cut_<label>.param from EVERY row of the survey DataFrame `df`
    with df.configuration == configuration, taking the most conservative
    window across all matching runs: tdc_min = min of each channel's _min
    column, tdc_max = max of each channel's _max column -- wide enough to
    cover every reference run for that configuration. `out_run` sets the
    label used in the output filenames (defaults to the smallest run number
    among the matching rows, since run is df's index)."""
    rows = df.loc[df.configuration == configuration, :]
    if rows.empty:
        raise ValueError(f"No rows found for configuration={configuration!r}")

    label = configuration#out_run if out_run is not None else int(rows.index.min())

    tdc_min = [rows[f"{name}_min"].min() for name in DEFAULT_TDC_NAMES]
    tdc_max = [rows[f"{name}_max"].max() for name in DEFAULT_TDC_NAMES]

    # most conservative (smallest) lo across all channels feeding each param,
    # AND across every matching run
    param_value = {}
    for channel, params in PARAM_MAP.items():
        col = f"{channel}_min"
        if col not in rows.columns:
            continue
        lo = rows[col].min()  # most conservative across all matching runs
        if pd.isna(lo):
            continue
        for p in params:
            if p not in param_value or lo < param_value[p]:
                param_value[p] = lo

    def v(name):
        return -param_value.get(name, 0.0)

    tcoin_text = generate_tcoin_param(",".join(list(rows.index.astype(str))), tdc_min, tdc_max,
                                       v("t_coin_trig_tdcrefcut"), v("t_coin_trig_adcrefcut"))
    hms_text = generate_hms_param(",".join(list(rows.index.astype(str))), v("hdc_tdcrefcut"), v("hhodo_tdcrefcut"),
                                   v("hhodo_adcrefcut"), v("hcer_adcrefcut"), v("hcal_adcrefcut"))
    shms_text = generate_shms_param(",".join(list(rows.index.astype(str))), v("pdc_tdcrefcut"), v("phodo_tdcrefcut"),
                                     v("phodo_adcrefcut"), v("pngcer_adcrefcut"),
                                     v("phgcer_adcrefcut"), v("paero_adcrefcut"), v("pcal_adcrefcut"))

    paths = {
        os.path.join(out_dir, f"tcoin_{label}.param"): tcoin_text,
        os.path.join(out_dir, f"h_reftime_cut_coindaq_{label}.param"): hms_text,
        os.path.join(out_dir, f"p_reftime_cut_{label}.param"): shms_text,
    }
    for path, text in paths.items():
        with open(path, "w") as f:
            f.write(text)
    print(f"Wrote (from {len(rows)} run(s) matching configuration={configuration!r}: "
          f"{list(rows.index)}) as label={label}:")
    for p in paths:
        print(" ", p)
    return list(paths.keys())
rows = [survey_row(run, config, tcoin_p, hms_p, shms_p)
        for run, config, tcoin_p, hms_p, shms_p in zip(run_numbers, configuration_params, tcoin_params, hms_params, shms_params)]
df = pd.DataFrame(rows).set_index("run").sort_index()
# df.reset_index(inplace = True)

write_param_files_for_configuration(df, 51)

run_configurations          = run_db["Configuration"].to_numpy()
run_configurations_next     = np.array([0] + list(run_configurations[:-1]))
changing_point_min          = run_configurations != run_configurations_next
run_configurations_previous = list(run_configurations[1:]) + [99]
changing_point_max          = run_configurations != run_configurations_previous

run_period_mins     = [26088] +  list(run_db["Run Number"].to_numpy()[changing_point_min])
run_period_maxs     = list(run_db["Run Number"].astype(int).to_numpy()[changing_point_max])  + [int(np.max(run_db["Run Number"]))]
for i in range(len(run_period_mins)):
  run_period_min = run_period_mins[i]
  run_period_max = run_period_maxs[i]
  this_configuration = run_db.loc[(run_db["Run Number"] >= run_period_min) & (run_db["Run Number"] <= run_period_max), "Configuration"].unique()
  assert len(this_configuration) == 1
  this_configuration = this_configuration[0]
  this_experiment    = run_db.loc[(run_db["Run Number"] >= run_period_min) & (run_db["Run Number"] <= run_period_max), "Experiment"].unique()[0]
  this_kinematics    = run_db.loc[(run_db["Run Number"] >= run_period_min) & (run_db["Run Number"] <= run_period_max), "Kinematics Setting"].unique()[0]
  
  if ~np.isin(this_configuration, df.configuration):
    continue
  if i == 0:
    standard_database_txt = "#{}. {} {}\n".format(this_configuration, this_experiment, this_kinematics)
  else:
    standard_database_txt = standard_database_txt+ "#{}. {} {}\n".format(this_configuration, this_experiment, this_kinematics)

  standard_database_txt  = standard_database_txt+ "{}--{}\n".format(run_period_min, run_period_max)
  standard_database_txt  = standard_database_txt + 'g_ctp_parm_filename       = "DBASE/COIN/general.param"\n'
  standard_database_txt  = standard_database_txt + 'g_ctp_kinematics_filename = "DBASE/COIN/standard.kinematics"\n'
  standard_database_txt  = standard_database_txt + 'g_ctp_map_filename        = "MAPS/COIN/DETEC/coin.map"\n'
  standard_database_txt  = standard_database_txt + 'g_ctp_trigdet_filename    = "PARAM/TRIG/tcoin_{}.param"\n\n'.format(this_configuration)
  

  print(standard_database_txt)

  with open("DBASE/COIN/general.param", "r") as general_param_template:
    general_param_lines = general_param_template.readlines()
    
  general_param_this_run = copy(general_param_lines)
  general_param_this_run[19-1] = general_param_this_run[19-1].replace('.param', '_{}.param'.format(this_configuration))
  general_param_this_run[51-1] = general_param_this_run[51-1].replace('.param', '_{}.param'.format(this_configuration))
  with open("DBASE/COIN/general_{}.param".format(this_configuration), "w") as general_param_file_this_run:
    general_param_file_this_run.writelines(general_param_this_run)