#!/usr/bin/env python3
"""
reftime_cut_app.py

Interactive matplotlib app for setting Hall C reference-time cuts and TRIG
TDC time-window cuts by eye, and generating the three per-run param files
from scratch:

    tcoin.param
    h_reftime_cut_coindaq.param
    p_reftime_cut.param
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

# Disable toolbar to avoid Tk/X11 font forwarding crashes
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
    print("Try installing a GUI toolkit or pass --backend to force one explicitly.")
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

# which generated file each single-value param lives in, for the "used by" labeling
PARAM_FILE_OF = {
    "hdc_tdcrefcut": "hms", "hhodo_tdcrefcut": "hms", "hhodo_adcrefcut": "hms",
    "hcer_adcrefcut": "hms", "hcal_adcrefcut": "hms",
    "pdc_tdcrefcut": "shms", "phodo_tdcrefcut": "shms", "phodo_adcrefcut": "shms",
    "pngcer_adcrefcut": "shms", "phgcer_adcrefcut": "shms", "paero_adcrefcut": "shms",
    "pcal_adcrefcut": "shms",
    "t_coin_trig_tdcrefcut": "tcoin", "t_coin_trig_adcrefcut": "tcoin",
}


def describe_channel_usage(name, kind):
    groups: dict[str, list[str]] = {}
    if kind == "tdc":
        groups.setdefault("tcoin", []).extend(
            ["t_coin_TdcTimeWindowMin", "t_coin_TdcTimeWindowMax"])
    for p in PARAM_MAP.get(name, []):
        groups.setdefault(PARAM_FILE_OF.get(p, "tcoin"), []).append(p)
    parts = [f"{f}\u2014{', '.join(groups[f])}" for f in ("hms", "shms", "tcoin") if f in groups]
    return ", ".join(parts)


# ---- fixed boilerplate for tcoin.param ----
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
# Per-channel data
# ==========================================================================

def dominant_multiplicity(mult: np.ndarray, max_mult: int = 3):
    mult = np.asarray(mult)
    mult = mult[(mult >= 1) & (mult <= max_mult)]
    if mult.size == 0:
        return None
    vals, counts = np.unique(mult, return_counts=True)
    return int(vals[np.argmax(counts)])


def load_channel_data(raw, mult):
    raw = np.asarray(raw, dtype=float)
    mult = np.asarray(mult)
    finite = np.isfinite(raw)
    raw = raw[finite]
    mult = mult[finite]
    if raw.size == 0:
        return None

    dom_mult = dominant_multiplicity(mult)
    sel_mult1 = raw[mult == 1]
    sel_mult2 = raw[mult == 2]
    sel_mult3 = raw[mult == 3]

    return dict(raw=raw, dominant_mult=dom_mult,
                sel_mult1=sel_mult1 if sel_mult1.size > 0 else None,
                sel_mult2=sel_mult2 if sel_mult2.size > 0 else None,
                sel_mult3=sel_mult3 if sel_mult3.size > 0 else None,
                data_min=float(raw.min()), data_max=float(raw.max()))


MULT_COLORS = {1: "tab:green", 2: "tab:purple", 3: "tab:brown"}


def draw_mult_histograms(ax, d, rng, bins, raw_lw=0.8, mult_lw=0.9, fontsize=None):
    ax.hist(d["raw"], bins=bins, range=rng, histtype="stepfilled",
            color="tab:blue", alpha=0.3, lw=raw_lw, label="raw")
    for m in (1, 2, 3):
        sel = d.get(f"sel_mult{m}")
        if sel is None:
            continue
        label = f"mult={m}" + (" (dominant)" if d["dominant_mult"] == m else "")
        lw = mult_lw * (1.3 if d["dominant_mult"] == m else 1.0)
        ax.hist(sel, bins=bins, range=rng, histtype="step", color=MULT_COLORS[m], lw=lw, label=label)
    ax.set_yscale("log")


# ==========================================================================
# Param file generation
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

t_coin_trig_tdcrefcut = {trig_tdcrefcut:.1f}
t_coin_trig_adcrefcut = {trig_adcrefcut:.1f}

t_coin_adcNames = "{ADC_NAMES_FULL}"

t_coin_tdcNames = "{tdcnames_str}"

t_coin_TdcTimeWindowMin = {min_str}

t_coin_TdcTimeWindowMax = {max_str}
"""


_REFTIME_HEADER = """; Cut to select the Reference time when multiple hits in reference time
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
# Loading previous run's generated param files
# ==========================================================================

def _parse_scalar_params(text):
    result = {}
    for line in text.splitlines():
        m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([-+]?\d+\.?\d*)\s*$", line)
        if m:
            result[m.group(1)] = float(m.group(2))
    return result


def _extract_block(text, key):
    m = re.search(rf"{re.escape(key)}\s*=\s*(.*?)(?=\n\s*\n|\Z)", text, re.S)
    return m.group(1) if m else None


def load_reference(tcoin_path, hms_path, shms_path):
    try:
        tcoin_text = Path(tcoin_path).read_text()
        hms_text = Path(hms_path).read_text()
        shms_text = Path(shms_path).read_text()
    except OSError:
        return None

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

    channel_lo = {}
    for ch, params in PARAM_MAP.items():
        p0 = params[0]
        if p0 in scalars:
            channel_lo[ch] = -scalars[p0]

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
    source: str = "manual"


class ReftimeCutApp:
    def __init__(self, run, channels: list[tuple[str, str]], arrays, want_branches,
                 progress_path: Path = None, baseline_ref: dict = None, this_run_ref: dict = None, interactive: bool = True):
        self.run = run
        self.channels = channels
        self.arrays = arrays
        self.want_branches = want_branches
        self.progress_path = progress_path
        self.baseline_ref = baseline_ref  # Target baseline parameter cuts
        self.this_run_ref = this_run_ref  # Saved progress cuts for *this* run context
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
            print(f"Resuming saved progress from {self.progress_path}: {n_set} channels set.")

        # Seed initial session states using local progress/historical this-run data if available
        if self.this_run_ref is not None:
            for name, kind in self.channels:
                if name in self.results:
                    continue
                needed = name in NEEDED_SET
                ref_lo, ref_hi = self.get_reference_values(self.this_run_ref, name, needed)
                if ref_lo is None and ref_hi is None:
                    continue
                if ref_lo is not None:
                    ref_lo = max(0.0, ref_lo)
                if ref_hi is not None:
                    ref_hi = max(0.0, ref_hi)
                self.results[name] = ChannelState(name, kind, ref_lo, ref_hi, source="reference")

        self.interactive = interactive
        if not interactive:
            self.fig = None
            self.ax = None
            self.info_ax = None
            self.finished = True
            self._current_source = "manual"
            self._fresh_visit = True
            return

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

    def get_reference_values(self, ref_dict, name, needed):
        if ref_dict is None:
            return None, None
        ref_lo = ref_dict["channel_lo"].get(name) if needed else None
        win = ref_dict["window"].get(name)
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
        draw_mult_histograms(self.ax, d, rng, bins=300)
        margin = 0.02 * max(rng[1] - rng[0], 1.0)
        self._full_xlim = (rng[0] - margin, rng[1] + margin)
        self.ax.set_xlim(*self._full_xlim)

        # Plot Comparative Baseline Reference-Run Lines
        ref_lo, ref_hi = self.get_reference_values(self.baseline_ref, name, needed)
        if ref_lo is not None or ref_hi is not None:
            ylim = self.ax.get_ylim()
            y_marker = ylim[1] * 0.6
            ref_label_base = f"cuts-ref run {self.baseline_ref.get('run_num', 'default')}"
            ref_lo_str = f"{ref_lo:.1f}" if ref_lo is not None else "-"
            ref_hi_str = f"{ref_hi:.1f}" if ref_hi is not None else "-"
            full_ref_label = f"{ref_label_base} (lo={ref_lo_str}, hi={ref_hi_str})"

            first = True
            for x in (ref_lo, ref_hi):
                if x is None:
                    continue
                self.ax.axvline(x, color="dimgray", ls="--", lw=1.0, zorder=1)
                self.ax.plot([x], [y_marker], marker="v", color="dimgray", markersize=9,
                            markeredgecolor="black", markeredgewidth=0.4, zorder=6,
                            linestyle="None", label=full_ref_label if first else None)
                first = False
            self.ax.set_ylim(*ylim)

        # Restore / set current session's working cuts
        stored = self.results.get(name)
        if stored is not None and stored.lo is not None:
            self.clicks = [stored.lo] if stored.hi is None else [stored.lo, stored.hi]
            self._current_source = stored.source
        else:
            self.clicks = []
            self._current_source = "manual"
        self._fresh_visit = True

        usage = describe_channel_usage(name, kind)
        n_clicks_needed = 1 if kind == "adc" else 2
        self.ax.set_title(
            f"[{self.idx + 1}/{len(self.channels)}] {name} ({kind.upper()}) -- used by {usage}\n"
            f"click {n_clicks_needed} point(s): {'lo' if n_clicks_needed == 1 else 'lo, hi'} "
            f"| scroll to zoom",
            fontsize=9)
        self.redraw_clicks()

    def redraw_clicks(self):
        for ln in list(self.ax.lines):
            if getattr(ln, "_manual_cut_line", False):
                ln.remove()

        name, kind = self.current()
        if len(self.clicks) > 0:
            lo = self.clicks[0]
            hi = self.clicks[1] if len(self.clicks) >= 2 else None
            lo_str = f"{lo:.1f}"
            hi_str = f"{hi:.1f}" if hi is not None else ("(n/a)" if kind == "adc" else "-")
            legend_label = f"cuts-this run {self.run} (lo={lo_str}, hi={hi_str})"

            first = True
            for x in self.clicks:
                ln = self.ax.axvline(x, color="orange", ls="-", lw=1.5, label=legend_label if first else None)
                ln._manual_cut_line = True
                first = False

        self.ax.legend(fontsize=8, loc="upper right")
        self.update_info_panel()
        self.fig.canvas.draw_idle()

    def update_info_panel(self):
        self.info_ax.clear()
        self.info_ax.axis("off")
        name, _ = self.current()

        lines = [
            f"Run {self.run}",
            f"Channel: {name}",
            "",
            "-" * 30,
            f"NEEDED channels ({len(NEEDED_CHANNELS)} total):"
        ]
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

    # ---- interaction ----

    def on_click(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        name, kind = self.current()
        if self.data.get(name) is None:
            return
        
        # Clamp to zero so no lower/upper bound parameters can drop negative
        click_val = max(0.0, event.xdata)

        needs_two = kind == "tdc"
        if self._fresh_visit:
            self.clicks = [click_val]
            self._fresh_visit = False
        elif needs_two:
            if len(self.clicks) >= 2:
                self.clicks = [click_val]
            else:
                self.clicks.append(click_val)
        else:
            self.clicks = [click_val]
        self._current_source = "manual"
        self.redraw_clicks()

    def on_key(self, event):
        if getattr(self.goto_box, "capturekeystrokes", False):
            return
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
        text = text.strip()
        if not text:
            return
        target_idx = None
        if text.lstrip("-").isdigit():
            i = int(text)
            if 1 <= i <= len(self.channels):
                target_idx = i - 1
            else:
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
                else:
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
                self.results.pop(name, None)
        else:
            if len(self.clicks) == 1:
                self.results[name] = ChannelState(name, kind, self.clicks[0], None,
                                                   source=self._current_source)
            elif len(self.clicks) == 0:
                self.results.pop(name, None)
        self.save_progress()

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
        except (json.JSONDecodeError, OSError):
            return None
        results = {
            name: ChannelState(name, v["kind"], v.get("lo"), v.get("hi"),
                              source=v.get("source", "manual"))
            for name, v in data.get("results", {}).items()
        }
        return {"idx": data.get("idx", 0), "results": results}

    def on_pause(self, event):
        self._store_current()
        self.finished = False
        plt.close(self.fig)

    def on_next(self, event):
        self._store_current()
        if self.idx < len(self.channels) - 1:
            self.idx += 1
            self.draw_channel()

    def on_prev(self, event):
        self._store_current()
        if self.idx > 0:
            self.idx -= 1
            self.draw_channel()

    def on_next_required(self, event):
        self._store_current()
        j = self.idx + 1
        while j < len(self.channels) and self.channels[j][0] not in NEEDED_SET:
            j += 1
        if j < len(self.channels):
            self.idx = j
            self.draw_channel()

    def on_reset(self, event):
        self.clicks = []
        self._fresh_visit = False
        self.redraw_clicks()

    def on_use_reference(self, event):
        name, kind = self.current()
        needed = name in NEEDED_SET
        ref_lo, ref_hi = self.get_reference_values(self.baseline_ref, name, needed)
        if ref_lo is None and ref_hi is None:
            return
        
        if ref_lo is not None:
            ref_lo = max(0.0, ref_lo)
        if ref_hi is not None:
            ref_hi = max(0.0, ref_hi)

        if kind == "tdc":
            if ref_lo is not None and ref_hi is not None:
                self.clicks = [ref_lo, ref_hi]
            elif ref_lo is not None:
                self.clicks = [ref_lo]
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
        if not self.interactive:
            return
        plt.show()
        self._store_current()

    # ---- zoom ----

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

    def resolve_window_arrays(self):
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


# ==========================================================================
# Main
# ==========================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=int, help="Run number, e.g. 26092")
    ap.add_argument("--root-file", default=None)
    ap.add_argument("--tree", default="T")
    ap.add_argument("--out-dir", default="./reftime_qa")
    ap.add_argument("--backend", default=None)
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--reference-run", type=int, default=None)
    ap.add_argument("--reference-dir", default=None)
    ap.add_argument("--param-dir", default="../../PARAM")
    ap.add_argument("--non-interactive", action="store_true")
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
    arrays = tree.arrays(all_branch_names, library="np")

    if args.non_interactive:
        _init_backend("Agg")
    else:
        _init_backend(args.backend)

    progress_path = out_dir / f"{args.run}_progress.json"
    if args.fresh and progress_path.exists():
        progress_path.unlink()

    # ----------------------------------------------------------------------
    # 1. Determine "This Run" context configurations
    # ----------------------------------------------------------------------
    this_run_ref = None
    self_tcoin = out_dir / f"tcoin_{args.run}.param"
    self_hms = out_dir / f"h_reftime_cut_coindaq_{args.run}.param"
    self_shms = out_dir / f"p_reftime_cut_{args.run}.param"
    
    if all(p.exists() for p in (self_tcoin, self_hms, self_shms)):
        this_run_ref = load_reference(self_tcoin, self_hms, self_shms)
        if this_run_ref is not None:
            this_run_ref["label"] = f"this run {args.run}"
            this_run_ref["run_num"] = args.run

    if this_run_ref is None and progress_path.exists():
        try:
            p_data = json.loads(progress_path.read_text())
            prog_results = p_data.get("results", {})
            prog_window = {}
            prog_channel_lo = {}
            for ch, cs in prog_results.items():
                if cs.get("lo") is not None:
                    prog_channel_lo[ch] = cs.get("lo")
                if cs.get("lo") is not None or cs.get("hi") is not None:
                    prog_window[ch] = (cs.get("lo", 0.0), cs.get("hi", 100000.0))
            
            this_run_ref = {
                "window": prog_window, 
                "channel_lo": prog_channel_lo, 
                "label": f"progress {args.run}", 
                "run_num": args.run
            }
        except Exception:
            pass

    # ----------------------------------------------------------------------
    # 2. Determine "Baseline Reference Target" cuts (independent tracking)
    # ----------------------------------------------------------------------
    baseline_ref = None

    # Explicit baseline reference run requested via CLI
    if args.reference_run is not None:
        ref_dir = Path(args.reference_dir) if args.reference_dir else out_dir
        tcoin_ref = ref_dir / f"tcoin_{args.reference_run}.param"
        hms_ref = ref_dir / f"h_reftime_cut_coindaq_{args.reference_run}.param"
        shms_ref = ref_dir / f"p_reftime_cut_{args.reference_run}.param"
        missing = [p for p in (tcoin_ref, hms_ref, shms_ref) if not p.exists()]
        
        if missing:
            sys.exit(f"Error: reference-run {args.reference_run} requested but missing files: {[str(p) for p in missing]}")
        
        baseline_ref = load_reference(tcoin_ref, hms_ref, shms_ref)
        if baseline_ref is None:
            sys.exit(f"Error: could not parse reference files for run {args.reference_run}")
        
        baseline_ref["label"] = f"run {args.reference_run}"
        baseline_ref["run_num"] = args.reference_run

    # Cascade down to standard/default global environment values
    else:
        param_dir = Path(args.param_dir)
        tcoin_v = param_dir / "TRIG" / "tcoin.param"
        hms_v = param_dir / "HMS" / "GEN" / "h_reftime_cut_coindaq.param"
        shms_v = param_dir / "SHMS" / "GEN" / "p_reftime_cut.param"
        missing_v = [p for p in (tcoin_v, hms_v, shms_v) if not p.exists()]
        if not missing_v:
            baseline_ref = load_reference(tcoin_v, hms_v, shms_v)
            if baseline_ref is not None:
                baseline_ref["label"] = "default parameters"
                baseline_ref["run_num"] = "default"

    # Initialize visual workspace app
    app = ReftimeCutApp(args.run, channels, arrays, want_branches,
                        progress_path=None if args.non_interactive else progress_path,
                        baseline_ref=baseline_ref, this_run_ref=this_run_ref, 
                        interactive=not args.non_interactive)

    if args.non_interactive:
# Generate the diagnostic PDF report
        pdf_path = out_dir / f"reftime_cuts_{args.run}.pdf"
        print(f"Generating diagnostic PDF: {pdf_path}")
        
        from matplotlib.backends.backend_pdf import PdfPages
        
        # We temporarily enable interactive plotting or use a local figure 
        # so we don't run into GUI backend issues on a headless node.
        with PdfPages(pdf_path) as pdf:
            fig, ax = plt.subplots(figsize=(11, 8.5))
            
            for idx, (name, kind) in enumerate(channels):
                ax.clear()
                d = app.data.get(name)
                needed = name in NEEDED_SET
                
                if d is None:
                    ax.text(0.5, 0.5, f"{name} ({kind.upper()})\nNo Data Available", 
                            ha='center', va='center', fontsize=14)
                    ax.set_title(f"{name} ({kind.upper()}) - Missing Branch")
                    pdf.savefig(fig)
                    continue
                
                # Draw the multi-hit histograms
                rng = (d["data_min"], d["data_max"])
                draw_mult_histograms(ax, d, rng, bins=300)
                
                # Set range with a slight margin
                margin = 0.02 * max(rng[1] - rng[0], 1.0)
                ax.set_xlim(rng[0] - margin, rng[1] + margin)
                
                # Retrieve and plot baseline/reference cuts (grey dashed)
                ref_lo, ref_hi = app.get_reference_values(baseline_ref, name, needed)
                ref_label = None
                if ref_lo is not None or ref_hi is not None:
                    ref_lo_str = f"{ref_lo:.1f}" if ref_lo is not None else "-"
                    ref_hi_str = f"{ref_hi:.1f}" if ref_hi is not None else "-"
                    ref_label = f"cuts-ref run {baseline_ref.get('run_num', 'default')} (lo={ref_lo_str}, hi={ref_hi_str})"
                    
                    for x in (ref_lo, ref_hi):
                        if x is not None:
                            ax.axvline(x, color="dimgray", ls="--", lw=1.0)
                
                # Retrieve and plot current run cuts (orange solid)
                # Headless mode seeds results from 'this_run_ref' inside ReftimeCutApp initialization
                stored = app.results.get(name)
                this_label = None
                if stored is not None and stored.lo is not None:
                    lo, hi = stored.lo, stored.hi
                    lo_str = f"{lo:.1f}"
                    hi_str = f"{hi:.1f}" if hi is not None else "-"
                    this_label = f"cuts-this run {args.run} (lo={lo_str}, hi={hi_str})"
                    
                    for x in (lo, hi):
                        if x is not None:
                            ax.axvline(x, color="orange", ls="-", lw=1.5)
                
                # Add dummy points to force legend entries if they exist
                if ref_label:
                    ax.plot([], [], color="dimgray", ls="--", label=ref_label)
                if this_label:
                    ax.plot([], [], color="orange", ls="-", label=this_label)
                
                usage = describe_channel_usage(name, kind)
                ax.set_title(f"[{idx + 1}/{len(channels)}] {name} ({kind.upper()}) - {usage}", fontsize=10)
                ax.legend(fontsize=8, loc="upper right")
                
                pdf.savefig(fig)
            
            plt.close(fig)
        return

    app.run_interactive()

    if not app.finished:
        return

    # Process and write resulting modifications
    tdc_min, tdc_max = app.resolve_window_arrays()
    param_summary = app.resolve_param_values()

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


if __name__ == "__main__":
    main()
