#!/usr/bin/env python3
"""
reftime_cut_app.py

Interactive matplotlib app for setting Hall C reference-time cuts and TRIG
TDC time-window cuts by eye, and generating the three per-run param files
from scratch:

    tcoin.param
    h_reftime_cut_coindaq.param
    p_reftime_cut.param

--------------------------------------------------------------------------
Everything is set by hand -- there is no algorithmic suggestion. Instead,
every channel is seeded from a reference (see below) so that "don't touch
it" means "keep the reference value," not "silently disable the cut."

1. The 14 single-value reference-time cuts (pdc_tdcrefcut, phodo_tdcrefcut,
   phodo_adcrefcut, pngcer_adcrefcut, phgcer_adcrefcut, paero_adcrefcut,
   pcal_adcrefcut, t_coin_trig_tdcrefcut, t_coin_trig_adcrefcut,
   hdc_tdcrefcut, hhodo_tdcrefcut, hhodo_adcrefcut, hcer_adcrefcut,
   hcal_adcrefcut). These come from the ~20 channels tagged "NEEDED" in
   the app. The right-hand panel flags any that are still unreviewed
   (tagged [ref]) so you don't miss one before saving.

2. The t_coin_TdcTimeWindowMin / t_coin_TdcTimeWindowMax arrays in
   tcoin.param, one min/max pair per entry of t_coin_tdcNames (115
   channels). Looking at a real tcoin.param, most of these are legitimately
   left wide open (0, 100000) -- only a handful of channels actually carry
   a real window, and that's exactly what a reference (or the vanilla
   defaults) already encodes.

Reference: pass --reference-run to seed every channel from a previously
generated tcoin_<run>.param / h_reftime_cut_coindaq_<run>.param /
p_reftime_cut_<run>.param. If you don't, the app falls back to the
standard hallc_replay layout under --param-dir (default ./PARAM):
    <param-dir>/TRIG/tcoin.param
    <param-dir>/HMS/GEN/h_reftime_cut_coindaq.param
    <param-dir>/SHMS/GEN/p_reftime_cut.param
and seeds from those instead, so there's always a sane baseline. Use the
"Use Reference" button to explicitly (re)accept a channel's reference value.

Click twice on a TDC channel's histogram to set (lo, hi); click once on an
ADC channel's histogram to set lo only. A third click starts over.

Navigation: << Prev / Accept & Next >> step through all 117 channels in
order (also Left/Right arrow keys). "Skip to next NEEDED channel >>" jumps
straight between the ~20 channels that actually feed one of the 14 param
values above, skipping over the ~94 other TRIG TDC channels (scalers,
trigger words, prescales, RF, ...) that usually don't need a window at all
-- use this if you don't care about reviewing every single TRIG channel.

Zoom: scroll the mouse wheel over the plot to zoom in/out centered on the
cursor, or use the Zoom In / Zoom Out / Reset View buttons.
--------------------------------------------------------------------------

Branch naming convention:

    TDC channel <name>:  T.coin.<name>_tdcTimeRaw / _tdcMultiplicity
    ADC channel <name>:  T.coin.<name>_adcPulseTimeRaw / _adcMultiplicity

Usage:
    python reftime_cut_app.py 26092
    python reftime_cut_app.py 26092 --out-dir ./reftime_qa

By default the ROOT file is:
    /volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_<run>_2000000_0.root
Pass --root-file to override, and --tree to use a Tree name other than "T".
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import uproot
except ImportError:
    sys.exit("This script requires uproot: pip install uproot --break-system-packages")

import matplotlib

# The toolbar (NavigationToolbar2Tk in particular) can fail over X11 forwarding even when
# the base backend works fine (known Tk/X11 font issue) -- disabling it sidesteps that, and
# we don't need it for this click-to-set-cuts workflow anyway.
matplotlib.rcParams["toolbar"] = "None"

plt = None  # bound by _init_backend() once a working backend is found


def _init_backend(preferred=None):
    global plt
    candidates = [preferred] if preferred else ["MacOSX", "QtAgg", "Qt5Agg", "TkAgg"]
    tried = []
    for candidate in candidates:
        try:
            matplotlib.use(candidate, force=True)
            import matplotlib.pyplot as _plt
            fig = _plt.figure()
            _plt.close(fig)
            plt = _plt
            return candidate
        except Exception as e:
            tried.append((candidate, str(e)))
    print("Could not initialize any interactive matplotlib backend. Tried:")
    for name, err in tried:
        print(f"  {name}: {err}")
    print("Try installing a GUI toolkit (pip install PyQt5 --break-system-packages, or "
          "sudo apt-get install python3-tk), or pass --backend to force one explicitly.")
    sys.exit(1)


from matplotlib.widgets import Button, TextBox


# ==========================================================================
# Channel lists
# ==========================================================================

DEFAULT_TDC_NAMES = (
    "h1X h1Y h2X h2Y h1T h2T hT1 hASUM hBSUM hCSUM hDSUM hPRLO hPRHI hSHWR hEDTM hCER hT2 "
    "hDCREF1 hDCREF2 hDCREF3 hDCREF4 "
    "hTRIG1_ROC1 hTRIG2_ROC1 hTRIG3_ROC1 hTRIG4_ROC1 hTRIG5_ROC1 hTRIG6_ROC1 "
    "pTRIG1_ROC1 pTRIG2_ROC1 pTRIG3_ROC1 pTRIG4_ROC1 pTRIG5_ROC1 pTRIG6_ROC1 "
    "pT1 pT2 p1X p1Y p2X p2Y p1T p2T pT3 pAER pHGCER pNGCER "
    "pDCREF1 pDCREF2 pDCREF3 pDCREF4 pDCREF5 pDCREF6 pDCREF7 pDCREF8 pDCREF9 pDCREF10 "
    "pEDTM pPRLO pPRHI "
    "pTRIG1_ROC2 pTRIG2_ROC2 pTRIG3_ROC2 pTRIG4_ROC2 pTRIG5_ROC2 pTRIG6_ROC2 "
    "hTRIG1_ROC2 hTRIG2_ROC2 hTRIG3_ROC2 hTRIG4_ROC2 hTRIG5_ROC2 hTRIG6_ROC2 "
    "pSTOF_ROC2 pEL_LO_LO_ROC2 pEL_LO_ROC2 pEL_HI_ROC2 pEL_REAL_ROC2 pEL_CLEAN_ROC2 "
    "hSTOF_ROC2 hEL_LO_LO_ROC2 hEL_LO_ROC2 hEL_HI_ROC2 hEL_REAL_ROC2 hEL_CLEAN_ROC2 "
    "pSTOF_ROC1 pEL_LO_LO_ROC1 pEL_LO_ROC1 pEL_HI_ROC1 pEL_REAL_ROC1 pEL_CLEAN_ROC1 "
    "hSTOF_ROC1 hEL_LO_LO_ROC1 hEL_LO_ROC1 hEL_HI_ROC1 hEL_REAL_ROC1 hEL_CLEAN_ROC1 "
    "pPRE40_ROC1 pPRE100_ROC1 pPRE150_ROC1 pPRE200_ROC1 "
    "hPRE40_ROC1 hPRE100_ROC1 hPRE150_ROC1 hPRE200_ROC1 "
    "pPRE40_ROC2 pPRE100_ROC2 pPRE150_ROC2 pPRE200_ROC2 "
    "hPRE40_ROC2 hPRE100_ROC2 hPRE150_ROC2 hPRE200_ROC2 "
    "hDCREF5 hRF pRF hHODO_RF pHODO_RF"
).split()

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

# ---- fixed boilerplate for tcoin.param (structural, doesn't change run to run) ----
TRIG_NAMES = "pTRIG1_ROC1 pTRIG4_ROC1 pTRIG1_ROC2 pTRIG4_ROC2"
ADC_NAMES_FULL = ("hASUM hBSUM hCSUM hDSUM hPSHWR hSHWR hAER hCER hFADC_TREF_ROC1 pAER "
                   "pHGCER pNGCER pPSHWR pFADC_TREF_ROC2 pHGCER_MOD pNGCER_MOD pHEL_NEG "
                   "pHEL_POS pHEL_MPS")
NUM_ADC = 19
NUM_TDC = 115
TDC_OFFSET = 300.0
ADC_TDC_OFFSET = 200.0
TDC_CHANPERNS = 0.09766
EHADCOINTIME_OFFSET = 0.0


def tdc_branches(name: str) -> tuple[str, str]:
    return f"T.coin.{name}_tdcTimeRaw", f"T.coin.{name}_tdcMultiplicity"


def adc_branches(name: str) -> tuple[str, str]:
    return f"T.coin.{name}_adcPulseTimeRaw", f"T.coin.{name}_adcMultiplicity"


# ==========================================================================
# Per-channel data (no cut-finding, just the histogram inputs)
# ==========================================================================

def dominant_multiplicity(mult: np.ndarray, max_mult: int = 3):
    mult = np.asarray(mult)
    mult = mult[(mult >= 1) & (mult <= max_mult)]
    if mult.size == 0:
        return None
    vals, counts = np.unique(mult, return_counts=True)
    return int(vals[np.argmax(counts)])


def load_channel_data(raw, mult):
    """Return dict(raw, sel, sel_mult2, sel_mult3, dominant_mult, data_min, data_max)
    or None if the channel has no usable data. No cut is computed here -- the
    person sets lo/hi entirely by clicking on the plot."""
    raw = np.asarray(raw, dtype=float)
    mult = np.asarray(mult)
    finite = np.isfinite(raw)
    raw = raw[finite]
    mult = mult[finite]
    if raw.size == 0:
        return None

    dom_mult = dominant_multiplicity(mult)
    sel = raw if dom_mult is None else raw[mult == dom_mult]
    if sel.size < 20:
        sel = raw
        dom_mult = None

    sel_mult2 = raw[mult == 2]
    sel_mult3 = raw[mult == 3]

    return dict(raw=raw, sel=sel, dominant_mult=dom_mult,
                sel_mult2=sel_mult2 if sel_mult2.size > 0 else None,
                sel_mult3=sel_mult3 if sel_mult3.size > 0 else None,
                data_min=float(raw.min()), data_max=float(raw.max()))


# ==========================================================================
# Param file generation (written from scratch, not by editing a template)
# ==========================================================================

def format_array(values, per_line=10, indent="\t\t\t\t  "):
    def fmt(v):
        return f"{v:.0f}" if float(v).is_integer() else f"{v:.2f}"
    rows = []
    for i in range(0, len(values), per_line):
        rows.append(", ".join(fmt(v) for v in values[i:i + per_line]))
    return (",\n" + indent).join(rows)


def generate_tcoin_param(run, tdc_min, tdc_max, trig_tdcrefcut, trig_adcrefcut):
    tdcnames_str = " ".join(DEFAULT_TDC_NAMES)
    min_str = format_array(tdc_min)
    max_str = format_array(tdc_max)
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""; Auto-generated by reftime_cut_app.py for run {run} on {now}

t_coin_numAdc = {NUM_ADC}
t_coin_numTdc = {NUM_TDC}

t_coin_tdcoffset = {TDC_OFFSET}
t_coin_adc_tdc_offset = {ADC_TDC_OFFSET}

t_coin_tdcchanperns = {TDC_CHANPERNS}
eHadCoinTime_Offset = {EHADCOINTIME_OFFSET}

t_coin_trigNames="{TRIG_NAMES}"

; tdc cut is on pTRef2
; adc cut on pFADC_ROC2
t_coin_trig_tdcrefcut = {trig_tdcrefcut:.1f}
t_coin_trig_adcrefcut = {trig_adcrefcut:.1f}

;NOTE: pTRIG1_ROC1, pTRIG4_ROC1, pTRIG1_ROC2, pTRIG4_ROC2 (tdcTimeRaw), are hard-coded in THcTrigDet.cxx class
;      These variables and their indices are used for the coincidence time calculation, and should NOT change here.
;      If you do change them, please make sure to modify the THcTrigDet.cxx class.

t_coin_adcNames = "{ADC_NAMES_FULL}"

t_coin_tdcNames = "{tdcnames_str}"

t_coin_TdcTimeWindowMin = {min_str}

t_coin_TdcTimeWindowMax = {max_str}
"""


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


def generate_hms_param(run, hdc, hhodo_tdc, hhodo_adc, hcer_adc, hcal_adc):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return (_REFTIME_HEADER + f"""
; Auto-generated by reftime_cut_app.py for run {run} on {now}
hdc_tdcrefcut={hdc:.1f}
hhodo_tdcrefcut={hhodo_tdc:.1f}
hhodo_adcrefcut={hhodo_adc:.1f}
hcer_adcrefcut={hcer_adc:.1f}
hcal_adcrefcut={hcal_adc:.1f}
""")


def generate_shms_param(run, pdc, phodo_tdc, phodo_adc, pngcer_adc, phgcer_adc, paero_adc, pcal_adc):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return (_REFTIME_HEADER + f"""
; Auto-generated by reftime_cut_app.py for run {run} on {now}
pdc_tdcrefcut={pdc:.1f}
phodo_tdcrefcut={phodo_tdc:.1f}
phodo_adcrefcut={phodo_adc:.1f}
pngcer_adcrefcut={pngcer_adc:.1f}
phgcer_adcrefcut={phgcer_adc:.1f}
paero_adcrefcut={paero_adc:.1f}
pcal_adcrefcut={pcal_adc:.1f}
""")


# ==========================================================================
# Loading a previous run's generated param files as a visual reference
# ==========================================================================

def _parse_scalar_params(text):
    """key = number (or key=number) assignments, one per line -- deliberately
    excludes array lines like 't_coin_TdcTimeWindowMin = 0, 0, ...' since
    those don't end the line with just a bare number."""
    result = {}
    for line in text.splitlines():
        m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([-+]?\d+\.?\d*)\s*$", line)
        if m:
            result[m.group(1)] = float(m.group(2))
    return result


def _extract_block(text, key):
    """Everything after 'key = ' up to the next blank line (matches the
    blank-line-separated layout generate_tcoin_param() writes)."""
    m = re.search(rf"{re.escape(key)}\s*=\s*(.*?)(?=\n\s*\n|\Z)", text, re.S)
    return m.group(1) if m else None


def load_reference(tcoin_path, hms_path, shms_path):
    """Parse a previously-generated set of the three param files and return
    dict(window={name: (lo, hi)}, channel_lo={channel: lo}) for use as
    reference lines. Returns None if the tcoin.param's arrays can't be
    parsed (e.g. an unrelated/hand-edited file)."""
    tcoin_text = Path(tcoin_path).read_text()
    hms_text = Path(hms_path).read_text()
    shms_text = Path(shms_path).read_text()

    scalars = {}
    scalars.update(_parse_scalar_params(tcoin_text))
    scalars.update(_parse_scalar_params(hms_text))
    scalars.update(_parse_scalar_params(shms_text))

    names_block = _extract_block(tcoin_text, "t_coin_tdcNames")
    names_match = re.search(r'"([^"]*)"', names_block) if names_block else None
    names = names_match.group(1).split() if names_match else []

    min_block = _extract_block(tcoin_text, "t_coin_TdcTimeWindowMin")
    max_block = _extract_block(tcoin_text, "t_coin_TdcTimeWindowMax")
    mins = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*", min_block)] if min_block else []
    maxs = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*", max_block)] if max_block else []

    window = {}
    if names and len(names) == len(mins) == len(maxs):
        for nm, lo, hi in zip(names, mins, maxs):
            window[nm] = (lo, hi)

    # reverse PARAM_MAP -> per-channel lo, using each channel's first (most
    # exclusive) param as the representative value
    channel_lo = {}
    for ch, params in PARAM_MAP.items():
        p0 = params[0]
        if p0 in scalars:
            channel_lo[ch] = -scalars[p0]  # undo the param file's sign flip

    if not window and not channel_lo:
        return None
    return {"window": window, "channel_lo": channel_lo}


# ==========================================================================
# Interactive app
# ==========================================================================

@dataclass
class ChannelState:
    name: str
    kind: str  # "tdc" or "adc"
    lo: float = None
    hi: float = None  # stays None for "adc" kind
    source: str = "manual"  # "manual" (you clicked it) or "reference" (auto-seeded, unreviewed)


class ReftimeCutApp:
    def __init__(self, run, channels: list[tuple[str, str]], arrays, want_branches,
                 progress_path: Path = None, reference: dict = None):
        self.run = run
        self.channels = channels  # [(name, kind), ...] in order
        self.arrays = arrays
        self.want_branches = want_branches
        self.progress_path = progress_path
        self.reference = reference  # dict(run=int, window={...}, channel_lo={...}) or None
        self.clicks: list[float] = []
        self.results: dict[str, ChannelState] = {}
        self.data: dict[str, dict] = {}
        self._full_xlim = (0.0, 1.0)
        self._load_all_channel_data()

        self.idx = 0
        loaded = self.load_progress() if self.progress_path else None
        if loaded is not None:
            self.idx = min(loaded["idx"], len(channels) - 1)
            self.results = loaded["results"]
            n_set = len(self.results)
            print(f"Resuming saved progress from {self.progress_path}: "
                  f"{n_set} channel(s) already set, continuing at "
                  f"channel {self.idx + 1}/{len(channels)}.")

        if self.reference is not None:
            n_seeded = 0
            for name, kind in self.channels:
                if name in self.results:
                    continue  # saved progress always takes priority over the reference
                needed = name in NEEDED_SET
                ref_lo, ref_hi = self.get_reference(name, needed)
                if ref_lo is None and ref_hi is None:
                    continue
                self.results[name] = ChannelState(name, kind, ref_lo, ref_hi, source="reference")
                n_seeded += 1
            print(f"Seeded {n_seeded} channel(s) from the {self.reference['label']} reference "
                  f"-- anything you don't click keeps that value. Use 'Use Reference' to "
                  f"explicitly (re)accept one, or just click to override it.")

        self.fig = plt.figure(figsize=(15, 7.5))
        self.fig.suptitle(f"Run {run} -- Hall C Reference-Time / TDC Window Cut Tool",
                          fontsize=13, fontweight="bold")

        self.ax = self.fig.add_axes([0.06, 0.30, 0.62, 0.58])
        self.info_ax = self.fig.add_axes([0.72, 0.30, 0.26, 0.58])
        self.info_ax.axis("off")

        ax_goto_label = self.fig.add_axes([0.72, 0.90, 0.02, 0.045])
        ax_goto_label.axis("off")
        ax_goto_label.text(0, 0.5, "Go to:", transform=ax_goto_label.transAxes,
                           va="center", ha="left", fontsize=9)
        ax_goto_box = self.fig.add_axes([0.775, 0.90, 0.13, 0.045])
        self.goto_box = TextBox(ax_goto_box, "", initial="")
        self.goto_box.on_submit(self.on_goto_submit)
        ax_goto_btn = self.fig.add_axes([0.91, 0.90, 0.07, 0.045])
        self.goto_button = Button(ax_goto_btn, "Go")
        self.goto_button.on_clicked(lambda event: self.on_goto_submit(self.goto_box.text))

        row1 = [
            ("<< Prev", self.on_prev),
            ("Accept & Next >>", self.on_next),
            ("Reset clicks", self.on_reset),
            ("Use Reference", self.on_use_reference),
            ("Skip to next NEEDED >>", self.on_next_required),
        ]
        row2 = [
            ("Zoom In (+)", self.on_zoom_in),
            ("Zoom Out (-)", self.on_zoom_out),
            ("Reset View", self.on_reset_view),
            ("Pause (save && exit)", self.on_pause),
            ("Save && Finish", self.on_finish),
        ]
        self.buttons = self._make_button_row(row1, y=0.16)
        self.buttons += self._make_button_row(row2, y=0.06)

        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("scroll_event", self.on_scroll)

        self.finished = False
        self._current_source = "manual"
        self._fresh_visit = True
        self.draw_channel()

    def _make_button_row(self, specs, y, height=0.07):
        n = len(specs)
        margin, gap = 0.06, 0.015
        total_w = 0.92 - margin
        w = (total_w - gap * (n - 1)) / n
        buttons = []
        for i, (label, cb) in enumerate(specs):
            x = margin + i * (w + gap)
            ax_b = self.fig.add_axes([x, y, w, height])
            b = Button(ax_b, label)
            b.on_clicked(cb)
            buttons.append(b)
        return buttons

    # ---- data ----

    def _load_all_channel_data(self):
        for name, kind in self.channels:
            branches = self.want_branches.get(name)
            if branches is None:
                self.data[name] = None
                continue
            raw_b, mult_b = branches
            self.data[name] = load_channel_data(self.arrays[raw_b], self.arrays[mult_b])

    def current(self):
        return self.channels[self.idx]

    def get_reference(self, name, needed):
        """(ref_lo, ref_hi), either possibly None, from self.reference for this
        channel. lo prefers the single-value param (most authoritative for a
        NEEDED channel); hi comes from the window array only if it's a real
        (non-default) window, since most TDC window entries are just (0, 1e5)."""
        if self.reference is None:
            return None, None
        ref_lo = self.reference["channel_lo"].get(name) if needed else None
        win = self.reference["window"].get(name)
        ref_hi = None
        if win is not None and win != (0.0, 100000.0):
            if ref_lo is None:
                ref_lo = win[0]
            ref_hi = win[1]
        return ref_lo, ref_hi

    # ---- drawing ----

    def draw_channel(self):
        self.ax.clear()
        name, kind = self.current()
        d = self.data.get(name)
        needed = name in NEEDED_SET

        if d is None:
            self.ax.set_title(f"{name} ({kind.upper()}) -- missing branch(es), no data", fontsize=11)
            self.ax.axis("off")
            self.clicks = []
            self._full_xlim = (0.0, 1.0)
            self.update_info_panel()
            self.fig.canvas.draw_idle()
            return

        rng = (d["data_min"], d["data_max"])
        self.ax.hist(d["raw"], bins=300, range=rng, histtype="stepfilled",
                     color="tab:blue", alpha=0.3, lw=0.8, label="raw")
        mlabel = f"mult={d['dominant_mult']} (dominant)" if d["dominant_mult"] is not None else "raw (no mult sel.)"
        self.ax.hist(d["sel"], bins=300, range=rng, histtype="step",
                     color="tab:red", lw=1.0, label=mlabel)
        if d.get("sel_mult2") is not None and d["dominant_mult"] != 2:
            self.ax.hist(d["sel_mult2"], bins=300, range=rng, histtype="step",
                         color="tab:green", lw=0.9, label="mult=2")
        if d.get("sel_mult3") is not None and d["dominant_mult"] != 3:
            self.ax.hist(d["sel_mult3"], bins=300, range=rng, histtype="step",
                         color="tab:purple", lw=0.9, label="mult=3")

        self.ax.set_yscale("log")
        margin = 0.02 * max(rng[1] - rng[0], 1.0)
        self._full_xlim = (rng[0] - margin, rng[1] + margin)
        self.ax.set_xlim(*self._full_xlim)

        ref_lo, ref_hi = self.get_reference(name, needed)
        if ref_lo is not None or ref_hi is not None:
            ylim = self.ax.get_ylim()
            y_marker = ylim[1] * 0.6
            ref_label = f"ref ({self.reference['label']})"
            for x in (ref_lo, ref_hi):
                if x is None:
                    continue
                self.ax.axvline(x, color="dimgray", ls="--", lw=1.0, zorder=1)
                # a marker at a fixed height stays visible even when the cut line below
                # lands exactly on top of the reference (the norm now that channels are
                # seeded from it) -- a plain axvline-on-axvline would be invisible then.
                self.ax.plot([x], [y_marker], marker="v", color="dimgray", markersize=9,
                            markeredgecolor="black", markeredgewidth=0.4, zorder=6,
                            linestyle="None", label=ref_label)
                ref_label = None
            self.ax.set_ylim(*ylim)

        # restore any previously-recorded clicks for this channel
        stored = self.results.get(name)
        if stored is not None and stored.lo is not None:
            self.clicks = [stored.lo] if stored.hi is None else [stored.lo, stored.hi]
            self._current_source = stored.source
        else:
            self.clicks = []
            self._current_source = "manual"
        self._fresh_visit = True  # first click on this visit starts a new pair,
                                   # instead of appending to a restored/seeded value

        tag = "[NEEDED -- feeds a .param value]" if needed else "[optional -- TDC window only]"
        n_clicks_needed = 1 if kind == "adc" else 2
        self.ax.set_title(
            f"[{self.idx + 1}/{len(self.channels)}] {name} ({kind.upper()}) {tag}\n"
            f"click {n_clicks_needed} point(s): {'lo' if n_clicks_needed == 1 else 'lo, hi'} "
            f"| scroll to zoom",
            fontsize=10)
        self.ax.legend(fontsize=7, loc="upper right")
        self.redraw_clicks()

    def redraw_clicks(self):
        for ln in list(self.ax.lines):
            if getattr(ln, "_manual_cut_line", False):
                ln.remove()
        for x in self.clicks:
            ln = self.ax.axvline(x, color="orange", ls="-", lw=1.5)
            ln._manual_cut_line = True
        self.update_info_panel()
        self.fig.canvas.draw_idle()

    def update_info_panel(self):
        self.info_ax.clear()
        self.info_ax.axis("off")
        name, kind = self.current()
        needed = name in NEEDED_SET

        lines = [f"Run {self.run}", ""]
        lines.append(f"Channel: {name} ({kind.upper()})")
        lines.append("Status: NEEDED" if needed else "Status: optional (window only)")
        lines.append("")
        lo_str = f"{self.clicks[0]:.1f}" if len(self.clicks) >= 1 else "-- not set --"
        hi_str = (f"{self.clicks[1]:.1f}" if len(self.clicks) >= 2
                  else ("(n/a)" if kind == "adc" else "-- not set --"))
        src_str = f" [{self._current_source}]" if self.clicks else ""
        lines.append(f"lo = {lo_str}{src_str}")
        lines.append(f"hi = {hi_str}")
        if self.reference is not None:
            ref_lo, ref_hi = self.get_reference(name, needed)
            if ref_lo is not None or ref_hi is not None:
                ref_lo_str = f"{ref_lo:.1f}" if ref_lo is not None else "--"
                ref_hi_str = f"{ref_hi:.1f}" if ref_hi is not None else "--"
                lines.append(f"(ref {self.reference['label']}: lo={ref_lo_str}, hi={ref_hi_str})")
        lines.append("")
        lines.append("-" * 30)
        lines.append(f"NEEDED channels ({len(NEEDED_CHANNELS)} total):")
        for n in NEEDED_CHANNELS:
            r = self.results.get(n)
            mark = ">" if n == name else " "
            if r is not None and r.lo is not None:
                tag = " [ref]" if r.source == "reference" else ""
                lines.append(f"{mark} {n}: {r.lo:.1f}{tag}")
            else:
                lines.append(f"{mark} {n}: UNSET")

        self.info_ax.text(0.0, 1.0, "\n".join(lines), transform=self.info_ax.transAxes,
                          fontsize=7.5, va="top", ha="left", family="monospace")

    # ---- interaction: clicks / navigation ----

    def on_click(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        name, kind = self.current()
        if self.data.get(name) is None:
            return
        needs_two = kind == "tdc"
        if self._fresh_visit:
            self.clicks = [event.xdata]
            self._fresh_visit = False
        elif needs_two:
            if len(self.clicks) >= 2:
                self.clicks = [event.xdata]
            else:
                self.clicks.append(event.xdata)
        else:
            self.clicks = [event.xdata]
        self._current_source = "manual"
        self.redraw_clicks()

    def on_key(self, event):
        if getattr(self.goto_box, "capturekeystrokes", False):
            return  # typing in the "Go to" box -- don't also trigger shortcuts
        if event.key == "right":
            self.on_next(None)
        elif event.key == "left":
            self.on_prev(None)
        elif event.key in ("+", "="):
            self.on_zoom_in(None)
        elif event.key == "-":
            self.on_zoom_out(None)
        elif event.key == "0":
            self.on_reset_view(None)

    def on_goto_submit(self, text):
        """Jump directly to a channel by 1-based index (matching the [i/117]
        shown in the title) or by name (exact match, or unique prefix)."""
        text = text.strip()
        if not text:
            return
        target_idx = None
        if text.lstrip("-").isdigit():
            i = int(text)
            if 1 <= i <= len(self.channels):
                target_idx = i - 1
            else:
                print(f"Index {i} out of range (1-{len(self.channels)}).")
                return
        else:
            lower = text.lower()
            exact = [i for i, (n, _) in enumerate(self.channels) if n.lower() == lower]
            if exact:
                target_idx = exact[0]
            else:
                prefix = [i for i, (n, _) in enumerate(self.channels) if n.lower().startswith(lower)]
                if len(prefix) == 1:
                    target_idx = prefix[0]
                elif len(prefix) > 1:
                    names = [self.channels[i][0] for i in prefix[:10]]
                    print(f"'{text}' matches multiple channels: {names}"
                          f"{' ...' if len(prefix) > 10 else ''}. Be more specific.")
                    return
                else:
                    print(f"No channel matches '{text}'.")
                    return
        self._store_current()
        self.idx = target_idx
        self.draw_channel()
        self.goto_box.set_val("")

    def _store_current(self):
        name, kind = self.current()
        if self.data.get(name) is None:
            return
        needs_two = kind == "tdc"
        if needs_two:
            if len(self.clicks) == 2:
                lo, hi = sorted(self.clicks)
                self.results[name] = ChannelState(name, kind, lo, hi, source=self._current_source)
            elif len(self.clicks) == 0:
                self.results.pop(name, None)  # explicit "unset" -> default at generation
            # a single pending click with no second click: leave whatever was
            # previously stored untouched rather than discarding it
        else:
            if len(self.clicks) == 1:
                self.results[name] = ChannelState(name, kind, self.clicks[0], None,
                                                   source=self._current_source)
            elif len(self.clicks) == 0:
                self.results.pop(name, None)
        self.save_progress()

    # ---- pause / resume ----

    def save_progress(self):
        if self.progress_path is None:
            return
        data = {
            "run": self.run,
            "idx": self.idx,
            "results": {
                name: {"kind": cs.kind, "lo": cs.lo, "hi": cs.hi, "source": cs.source}
                for name, cs in self.results.items()
            },
        }
        self.progress_path.parent.mkdir(parents=True, exist_ok=True)
        self.progress_path.write_text(json.dumps(data, indent=2))

    def load_progress(self):
        if not self.progress_path.exists():
            return None
        try:
            data = json.loads(self.progress_path.read_text())
        except (json.JSONDecodeError, OSError) as e:
            print(f"Could not read saved progress at {self.progress_path} ({e}); "
                  f"starting fresh.")
            return None
        results = {
            name: ChannelState(name, v["kind"], v.get("lo"), v.get("hi"),
                              source=v.get("source", "manual"))
            for name, v in data.get("results", {}).items()
        }
        return {"idx": data.get("idx", 0), "results": results}

    def on_pause(self, event):
        self._store_current()
        print(f"\nProgress saved to {self.progress_path}. "
              f"Re-run the same command to resume where you left off.")
        self.finished = False
        plt.close(self.fig)

    def on_next(self, event):
        self._store_current()
        if self.idx < len(self.channels) - 1:
            self.idx += 1
            self.draw_channel()
        else:
            print("Reached the last channel. Click 'Save & Finish' to write the param files.")

    def on_prev(self, event):
        self._store_current()
        if self.idx > 0:
            self.idx -= 1
            self.draw_channel()

    def on_next_required(self, event):
        """Jump to the next channel that feeds one of the 14 .param values
        (i.e. skip the ~94 TRIG TDC channels that are usually left at the
        (0, 100000) 'no cut' default)."""
        self._store_current()
        j = self.idx + 1
        while j < len(self.channels) and self.channels[j][0] not in NEEDED_SET:
            j += 1
        if j < len(self.channels):
            self.idx = j
            self.draw_channel()
        else:
            print("No more NEEDED channels after this one.")

    def on_reset(self, event):
        self.clicks = []
        self._fresh_visit = False  # explicit reset already clears state; a following click
                                    # should accumulate normally, not re-trigger fresh-start logic
        self.redraw_clicks()

    def on_use_reference(self, event):
        """Explicitly (re)accept this channel's reference value(s), marking it
        as reviewed rather than just carried over silently."""
        name, kind = self.current()
        needed = name in NEEDED_SET
        ref_lo, ref_hi = self.get_reference(name, needed)
        if ref_lo is None and ref_hi is None:
            print(f"No reference value available for {name}.")
            return
        if kind == "tdc":
            if ref_lo is not None and ref_hi is not None:
                self.clicks = [ref_lo, ref_hi]
            elif ref_lo is not None:
                self.clicks = [ref_lo]  # only lo available -- second click still needed for hi
            else:
                print(f"Reference for {name} only has an hi value, which isn't enough to set "
                      f"a TDC window on its own -- click lo manually.")
                return
        else:
            self.clicks = [ref_lo] if ref_lo is not None else []
        self._current_source = "reference"
        self._fresh_visit = False
        self.redraw_clicks()

    def on_finish(self, event):
        self._store_current()
        self.finished = True
        plt.close(self.fig)

    def run_interactive(self):
        try:
            plt.show()  # blocks until the window is closed
        except Exception as e:
            print(f"\nWindow hit an error while running: {e}")
            print("If this looks like an X11/font error ('BadName', 'X_OpenFont'), it's a "
                  "known Tk-over-X11-forwarding issue. Try: python reftime_cut_app.py "
                  "<run> --backend QtAgg (after pip install PyQt5 --break-system-packages)")
            raise
        self._store_current()

    # ---- interaction: zoom ----

    def on_zoom_in(self, event):
        self._zoom(factor=0.7)

    def on_zoom_out(self, event):
        self._zoom(factor=1 / 0.7)

    def _zoom(self, factor):
        lo, hi = self.ax.get_xlim()
        center = 0.5 * (lo + hi)
        half = 0.5 * (hi - lo) * factor
        self.ax.set_xlim(center - half, center + half)
        self.fig.canvas.draw_idle()

    def on_reset_view(self, event):
        self.ax.set_xlim(*self._full_xlim)
        self.fig.canvas.draw_idle()

    def on_scroll(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        factor = 0.8 if event.button == "up" else 1.25
        lo, hi = self.ax.get_xlim()
        center = event.xdata
        new_lo = center - (center - lo) * factor
        new_hi = center + (hi - center) * factor
        self.ax.set_xlim(new_lo, new_hi)
        self.fig.canvas.draw_idle()

    # ---- final value resolution ----

    def resolve_window_arrays(self):
        """(min_list, max_list) for t_coin_TdcTimeWindowMin/Max, in DEFAULT_TDC_NAMES order.
        Default per channel, if never clicked, is (0, 100000) -- 'no cut'."""
        mins, maxs = [], []
        for name in DEFAULT_TDC_NAMES:
            r = self.results.get(name)
            if r is not None and r.lo is not None and r.hi is not None:
                mins.append(r.lo)
                maxs.append(r.hi)
            else:
                mins.append(0.0)
                maxs.append(100000.0)
        return mins, maxs

    def resolve_channel_lo(self, name):
        """lo value for a NEEDED channel: clicked value if set, else 0.0
        (no algorithm -- an unset NEEDED channel just defaults to 0.0, which
        is effectively 'no cut' under the sign convention used here)."""
        r = self.results.get(name)
        if r is not None and r.lo is not None:
            return r.lo
        return 0.0

    def resolve_param_values(self):
        by_param: dict[str, list[tuple[str, float]]] = {}
        for name, params in PARAM_MAP.items():
            lo = self.resolve_channel_lo(name)
            for p in params:
                by_param.setdefault(p, []).append((name, lo))
        summary = {}
        for p, items in by_param.items():
            best_name, best_lo = min(items, key=lambda kv: kv[1])
            summary[p] = {"param_value": -best_lo, "raw_lo": best_lo,
                          "from_channel": best_name, "n_channels_feeding": len(items)}
        return summary

    # ---- summary PDF ----

    def save_summary_pdf(self, out_path: Path, per_page=10, ncols=2):
        """Multi-page PDF, `per_page` channel histograms per page, showing raw
        + mult overlays, the reference line (if any), and the final resolved
        cut (orange) for every channel -- a static record of the whole session."""
        from matplotlib.backends.backend_pdf import PdfPages

        tdc_min, tdc_max = self.resolve_window_arrays()
        final_values = dict(zip(DEFAULT_TDC_NAMES, zip(tdc_min, tdc_max)))
        for name in ADC_CHANNELS:
            final_values[name] = (self.resolve_channel_lo(name), None)

        nrows = -(-per_page // ncols)  # ceil
        n_pages = -(-len(self.channels) // per_page)

        with PdfPages(out_path) as pdf:
            for page in range(n_pages):
                chunk = self.channels[page * per_page:(page + 1) * per_page]
                fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 3.2 * nrows))
                axes = np.atleast_1d(axes).ravel()
                for ax, (name, kind) in zip(axes, chunk):
                    self._plot_channel_static(ax, name, kind, final_values)
                for ax in axes[len(chunk):]:
                    ax.axis("off")
                fig.suptitle(f"Run {self.run} -- Reference-time / TDC window cuts "
                             f"(page {page + 1}/{n_pages})", fontsize=12, fontweight="bold")
                fig.tight_layout(rect=[0, 0, 1, 0.96])
                pdf.savefig(fig)
                plt.close(fig)
        print(f"Wrote summary PDF: {out_path} ({len(self.channels)} channels, {n_pages} pages)")

    def _plot_channel_static(self, ax, name, kind, final_values):
        d = self.data.get(name)
        needed = name in NEEDED_SET
        if d is None:
            ax.set_title(f"{name} ({kind.upper()}) -- missing", fontsize=8, color="gray")
            ax.axis("off")
            return

        rng = (d["data_min"], d["data_max"])
        ax.hist(d["raw"], bins=200, range=rng, histtype="stepfilled",
                color="tab:blue", alpha=0.3, lw=0.6)
        ax.hist(d["sel"], bins=200, range=rng, histtype="step", color="tab:red", lw=0.8)
        if d.get("sel_mult2") is not None and d["dominant_mult"] != 2:
            ax.hist(d["sel_mult2"], bins=200, range=rng, histtype="step", color="tab:green", lw=0.6)
        if d.get("sel_mult3") is not None and d["dominant_mult"] != 3:
            ax.hist(d["sel_mult3"], bins=200, range=rng, histtype="step", color="tab:purple", lw=0.6)
        ax.set_yscale("log")

        ref_lo, ref_hi = self.get_reference(name, needed)
        if ref_lo is not None or ref_hi is not None:
            ylim = ax.get_ylim()
            y_marker = ylim[1] * 0.6
            for x in (ref_lo, ref_hi):
                if x is None:
                    continue
                ax.axvline(x, color="dimgray", ls="--", lw=0.8, zorder=1)
                ax.plot([x], [y_marker], marker="v", color="dimgray", markersize=6,
                       markeredgecolor="black", markeredgewidth=0.3, zorder=6, linestyle="None")
            ax.set_ylim(*ylim)

        lo, hi = final_values.get(name, (None, None))
        if lo is not None:
            ax.axvline(lo, color="orange", lw=1.3)
        if hi is not None:
            ax.axvline(hi, color="orange", lw=1.3)

        r = self.results.get(name)
        src_tag = {"reference": " [ref]", "manual": ""}.get(r.source if r else None, " [default]")
        label_lines = []
        if ref_lo is not None or ref_hi is not None:
            ref_lo_s = f"{ref_lo:.1f}" if ref_lo is not None else "-"
            ref_hi_s = f"{ref_hi:.1f}" if ref_hi is not None else "-"
            label_lines.append(f"ref: {ref_lo_s} / {ref_hi_s}")
        lo_s = f"{lo:.1f}" if lo is not None else "-"
        hi_s = f"{hi:.1f}" if hi is not None else "-"
        label_lines.append(f"cut: {lo_s} / {hi_s}{src_tag}")
        ax.text(0.98, 0.03, "\n".join(label_lines), transform=ax.transAxes,
               fontsize=6, va="bottom", ha="right",
               bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="lightgray"))

        tag = "NEEDED" if needed else "optional"
        ax.set_title(f"{name} ({kind.upper()}, {tag})", fontsize=8)
        ax.tick_params(labelsize=6)


# ==========================================================================
# Main
# ==========================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=int, help="Run number, e.g. 26092")
    ap.add_argument("--root-file", default=None,
                     help="Override the default ROOT file path (default: "
                          "/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/"
                          "coin_replay_production_<run>_2000000_0.root)")
    ap.add_argument("--tree", default="T", help="Tree name (default: T)")
    ap.add_argument("--out-dir", default="./reftime_qa",
                     help="Directory for the generated per-run param files "
                          "(default ./reftime_qa)")
    ap.add_argument("--backend", default=None,
                     help="Force a specific matplotlib backend (e.g. QtAgg, TkAgg, MacOSX) "
                          "instead of auto-detecting. Use this if auto-detection picks a "
                          "backend that later hits a Tk/X11 font error over remote X forwarding.")
    ap.add_argument("--fresh", action="store_true",
                     help="Ignore any saved progress for this run and start over from scratch.")
    ap.add_argument("--reference-run", type=int, default=None,
                     help="Load a previously-generated set of "
                          "tcoin_<run>.param / h_reftime_cut_coindaq_<run>.param / "
                          "p_reftime_cut_<run>.param (from --reference-dir, default --out-dir) "
                          "and draw them as gray reference lines -- handy when a new run isn't "
                          "expected to differ much from one you already reviewed.")
    ap.add_argument("--reference-dir", default=None,
                     help="Directory to look for --reference-run's param files in "
                          "(default: --out-dir).")
    ap.add_argument("--param-dir", default="PARAM",
                     help="hallc_replay-style PARAM directory used for the vanilla default "
                          "reference when --reference-run isn't given -- expects "
                          "<param-dir>/TRIG/tcoin.param, <param-dir>/HMS/GEN/"
                          "h_reftime_cut_coindaq.param, and <param-dir>/SHMS/GEN/"
                          "p_reftime_cut.param. Default: ./PARAM")
    args = ap.parse_args()

    default_path = Path(
        "/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/"
        f"coin_replay_production_{args.run}_2000000_0.root"
    )
    root_path = Path(args.root_file) if args.root_file else default_path
    if not root_path.exists():
        sys.exit(f"File not found: {root_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    channels = [(n, "tdc") for n in DEFAULT_TDC_NAMES] + [(n, "adc") for n in ADC_CHANNELS]

    tree = uproot.open(root_path)[args.tree]
    available = set(tree.keys())
    want_branches = {}
    for name, kind in channels:
        raw_b, mult_b = tdc_branches(name) if kind == "tdc" else adc_branches(name)
        if raw_b in available and mult_b in available:
            want_branches[name] = (raw_b, mult_b)

    all_branch_names = sorted({b for pair in want_branches.values() for b in pair})
    print(f"Loading {len(all_branch_names)} branches from {root_path.name} "
          f"({len(channels)} channels requested, {len(want_branches)} found) ...")
    arrays = tree.arrays(all_branch_names, library="np")

    backend = _init_backend(args.backend)
    print(f"Using matplotlib backend: {backend}")

    progress_path = out_dir / f"{args.run}_progress.json"
    if args.fresh and progress_path.exists():
        progress_path.unlink()
        print(f"--fresh given: removed existing progress file {progress_path}")

    reference = None
    if args.reference_run is not None:
        ref_dir = Path(args.reference_dir) if args.reference_dir else out_dir
        tcoin_ref = ref_dir / f"tcoin_{args.reference_run}.param"
        hms_ref = ref_dir / f"h_reftime_cut_coindaq_{args.reference_run}.param"
        shms_ref = ref_dir / f"p_reftime_cut_{args.reference_run}.param"
        missing = [p for p in (tcoin_ref, hms_ref, shms_ref) if not p.exists()]
        if missing:
            print(f"WARNING: --reference-run {args.reference_run} given but missing file(s): "
                  f"{[str(p) for p in missing]}. Continuing without a reference.")
        else:
            reference = load_reference(tcoin_ref, hms_ref, shms_ref)
            if reference is None:
                print(f"WARNING: could not parse reference files for run {args.reference_run}. "
                      f"Continuing without a reference.")
            else:
                reference["label"] = f"run {args.reference_run}"
                print(f"Loaded reference from run {args.reference_run}: "
                      f"{len(reference['window'])} window entries, "
                      f"{len(reference['channel_lo'])} cut values.")

    if reference is None:
        # fall back to the standard hallc_replay param layout -- a fixed baseline that's
        # always there, so there's always something sane to seed from even if you haven't
        # reviewed a specific prior run yet.
        param_dir = Path(args.param_dir)
        tcoin_v = param_dir / "TRIG" / "tcoin.param"
        hms_v = param_dir / "HMS" / "GEN" / "h_reftime_cut_coindaq.param"
        shms_v = param_dir / "SHMS" / "GEN" / "p_reftime_cut.param"
        missing_v = [p for p in (tcoin_v, hms_v, shms_v) if not p.exists()]
        if missing_v:
            print(f"No --reference-run given, and vanilla default(s) not found: "
                  f"{[str(p) for p in missing_v]}. Continuing without a reference.")
        else:
            reference = load_reference(tcoin_v, hms_v, shms_v)
            if reference is not None:
                reference["label"] = "defaults (PARAM/)"
                print(f"No --reference-run given; using vanilla defaults from {param_dir}: "
                      f"{len(reference['window'])} window entries, "
                      f"{len(reference['channel_lo'])} cut values.")

    app = ReftimeCutApp(args.run, channels, arrays, want_branches,
                        progress_path=progress_path, reference=reference)
    print("\nInteractive window opened.")
    print("  Click twice per TDC channel (lo, hi), once per ADC channel (lo).")
    print("  'Skip to next NEEDED >>' jumps between the ~20 channels that feed a .param value.")
    print("  Scroll to zoom, or use the Zoom In/Out/Reset View buttons.")
    if reference is not None:
        print(f"  Gray dashed lines show the reference values from {reference['label']}.")
    print("  Progress auto-saves after every step. 'Pause' (or just closing the window) exits")
    print("  without writing the final .param files -- re-run this same command to resume.")
    print("  Click 'Save & Finish' when you're done reviewing to write the final files.\n")
    app.run_interactive()

    if not app.finished:
        print(f"\nSession ended without clicking 'Save & Finish' -- progress is saved at "
              f"{progress_path}, but no .param files were written this time.")
        print(f"Re-run the same command (python reftime_cut_app.py {args.run} ...) to resume.")
        return

    # ---- resolve final values and write files ----
    tdc_min, tdc_max = app.resolve_window_arrays()
    param_summary = app.resolve_param_values()

    print("\n" + "=" * 78)
    print("Suggested .param values")
    print("=" * 78)
    for p, info in sorted(param_summary.items()):
        n = info["n_channels_feeding"]
        extra = f"  (most conservative of {n} channels, from {info['from_channel']})" if n > 1 else ""
        print(f"  {p} = {info['param_value']:.1f}{extra}")

    def v(name):
        return param_summary[name]["param_value"]

    tcoin_path = out_dir / f"tcoin_{args.run}.param"
    hms_path = out_dir / f"h_reftime_cut_coindaq_{args.run}.param"
    shms_path = out_dir / f"p_reftime_cut_{args.run}.param"

    tcoin_path.write_text(generate_tcoin_param(
        args.run, tdc_min, tdc_max, v("t_coin_trig_tdcrefcut"), v("t_coin_trig_adcrefcut")))
    hms_path.write_text(generate_hms_param(
        args.run, v("hdc_tdcrefcut"), v("hhodo_tdcrefcut"), v("hhodo_adcrefcut"),
        v("hcer_adcrefcut"), v("hcal_adcrefcut")))
    shms_path.write_text(generate_shms_param(
        args.run, v("pdc_tdcrefcut"), v("phodo_tdcrefcut"), v("phodo_adcrefcut"),
        v("pngcer_adcrefcut"), v("phgcer_adcrefcut"), v("paero_adcrefcut"), v("pcal_adcrefcut")))

    pdf_path = out_dir / f"{args.run}_reftime_summary.pdf"
    app.save_summary_pdf(pdf_path, per_page=10, ncols=2)

    print(f"\nWrote:\n  {tcoin_path}\n  {hms_path}\n  {shms_path}\n  {pdf_path}")
    print(f"(Progress file {progress_path} left in place -- delete it, or pass --fresh next "
          f"time, if you want to start this run over from scratch.)")


if __name__ == "__main__":
    main()
