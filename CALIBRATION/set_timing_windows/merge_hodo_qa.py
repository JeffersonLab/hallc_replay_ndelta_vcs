#!/usr/bin/env python3
"""
merge_hodo_qa.py

Merges every student's hodo_diff_cuts_<run>.json under a QA directory into
one final SHMS+HMS hodoscope cut file, applying known SHMS dead-PMT rules
along the way so cuts accidentally set on channels that were physically
turned off don't end up in the final result. HMS is never filtered -- only
SHMS is subject to any of this.

--------------------------------------------------------------------------
RULES
--------------------------------------------------------------------------
1) UNIVERSALLY dead, every run including calibration -- SHMS 2Y quartz bars:
     2y Pos: PMTs 1, 2, 5, 19, 20, 21
     2y Neg: PMTs 1, 2, 19, 20, 21
   Any cut reported for these is ignored outright. Falls back to the
   vanilla PARAM/SHMS/HODO/phodo_cuts.param default.

2) Production-run-dependent -- SHMS 1X and 2X on/off PMTs. The calibration
   run (26088) is EXEMPT from this filter (its cuts are always trusted for
   every 1X/2X PMT). For any other run:
     26165 <= run < 26284:  1X on = {5,6,7}          2X on = {5,6,7,8}
     26284 <= run < 26302:  1X on = {5,6,7,8}         2X on = {5,6,7,8,9}
     run >= 26302:          1X on = {4,5,6,7,8,9}     2X on = {4,5,6,7,8,9,10}
   A cut reported by a NON-calibration run's file for a 1X/2X PMT NOT in
   that run's on-list is ignored (that channel was physically off, so the
   "cut" is presumably noise/background, not a real hit distribution).

3) Because calibration is exempt from (2) and just gets folded into the
   same union pool as valid production entries, a 1X/2X PMT that's off in
   EVERY production bracket naturally ends up using only the calibration
   run's value once filtering removes everything else -- no separate case
   needed for "other 1x/2x PMTs must come from 26088".

4) Final fallback, anything still completely unresolved (no valid
   contribution survives filtering, and no calibration entry exists
   either): the vanilla PARAM/{SHMS,HMS}/HODO/{p,h}hodo_cuts.param default,
   parsed the same PMT-major/plane-minor way the C++ tools do, with (0,0)
   treated as "disabled paddle, no reference" -- same convention.

Every channel that DOES have at least one valid (post-filter) contribution
is merged via widest-union (most inclusive): lo=min(lo), hi=max(hi) across
every file/run that reported a valid cut for it. 1Y and non-dead 2Y bars,
and all of HMS, are never filtered at all -- plain inclusive union.

In addition to the JSON report, this also writes real hcana-format
phodo_cuts.param / hhodo_cuts.param (PMT-major/plane-minor flat arrays,
same layout the vanilla files use): built by starting from the vanilla
arrays and overlaying every resolved channel on top, so untouched channels
correctly keep their existing vanilla value instead of being zeroed out.
These are STAGED output -- review before copying into the real
PARAM/{SHMS,HMS}/HODO/ directories. Use --no-write-param to skip this and
get only the JSON report.

Usage:
    python3 merge_hodo_qa.py [qa_dir] [-v VERSION] [-o OUTPUT] [--param-dir PARAM_DIR] [-r]

qa_dir defaults to ./hodo_diff_qa and --param-dir defaults to ../../PARAM,
so a bare invocation with no arguments at all covers the common case. With
no -o given, output defaults to <qa_dir>/hodo_diff_cuts_<version>.json
(--version defaults to 0), and phodo_cuts.param/hhodo_cuts.param default to
that same directory. A same-content compatibility copy is also written as
<qa_dir>/hodo_diff_cuts_<version>.json -- the exact filename
hodo_timediff_cut_app.C's --reference-run looks for -- so the merge is
directly usable for visual verification with no manual renaming:

    python3 merge_hodo_qa.py
    root -l -b -q 'hodo_timediff_cut_app.C(26088, 0, true)'   # 0 = default --version
"""
import argparse
import glob
import json
import os
import re
import sys
from datetime import datetime

PLANES = ["1x", "1y", "2x", "2y"]
PLANE_IDX = {p: i for i, p in enumerate(PLANES)}
CALIBRATION_RUN = 26088

# Rule 1: universally-dead SHMS 2Y quartz bars
DEAD_2Y = {
    "Pos": {1, 2, 5, 19, 20, 21},
    "Neg": {1, 2, 19, 20, 21},
}

# Rule 2: SHMS 1X/2X on-lists by production-run bracket
# (lo_inclusive, hi_exclusive_or_None, {"1x": set, "2x": set})
BRACKETS = [
    (26165, 26284, {"1x": {5, 6, 7},           "2x": {5, 6, 7, 8}}),
    (26284, 26302, {"1x": {5, 6, 7, 8},        "2x": {5, 6, 7, 8, 9}}),
    (26302, None,  {"1x": {4, 5, 6, 7, 8, 9},  "2x": {4, 5, 6, 7, 8, 9, 10}}),
]


def bracket_for_run(run):
    for lo, hi, onlists in BRACKETS:
        if run >= lo and (hi is None or run < hi):
            return onlists
    return None  # run < 26165, or otherwise outside every known bracket


def is_dead_2y(spec, plane, side, ipmt):
    if spec != "p" or plane != "2y":
        return False
    return ipmt in DEAD_2Y.get(side, set())


def is_1x2x_off_for_run(spec, plane, run, ipmt):
    """True if this (spec,plane,ipmt) contribution should be dropped
    because the PMT was off for THIS run's production bracket. Always
    False for HMS, non-1x/2x planes, or the calibration run (exempt)."""
    if spec != "p" or plane not in ("1x", "2x"):
        return False
    if run == CALIBRATION_RUN:
        return False
    onlists = bracket_for_run(run)
    if onlists is None:
        # run isn't in any known bracket -- can't validate it, so don't
        # silently trust it either; treat as off/ignore and flag it
        return True
    return ipmt not in onlists[plane]


# ---- vanilla param file parsing -- mirrors the C++ tools' ExtractBlock/
# FlatIndex exactly (same line-based block-end detection, same PMT-major/
# plane-minor indexing, same (0,0)-means-disabled convention) ----

def read_file(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return ""


def extract_block(text, key):
    start = text.find(key)
    if start == -1:
        return None
    eq = text.find("=", start)
    if eq == -1:
        return None
    pos = eq + 1
    n = len(text)
    end = n
    while pos < n:
        nl = text.find("\n", pos)
        line_end = n if nl == -1 else nl
        line = text[pos:line_end]
        trimmed = line.strip(" \t\r")
        if trimmed == "" or trimmed[0] == ';' or trimmed[0].isalpha() or trimmed[0] == '_':
            end = pos
            break
        if nl == -1:
            end = n
            break
        pos = nl + 1
    return text[eq + 1:end]


def split_numbers(block):
    if block is None:
        return []
    toks = re.split(r'[,\s]+', block.strip())
    out = []
    for t in toks:
        if t == '':
            continue
        try:
            out.append(float(t))
        except ValueError:
            pass
    return out


def load_vanilla_hodo(param_dir, spec):
    """dict: (plane,side,ipmt) -> (lo,hi) from the vanilla hcana param
    file. (0,0) is treated as a disabled paddle (no reference), matching
    the C++ tools' convention."""
    path = os.path.join(param_dir, "SHMS" if spec == "p" else "HMS", "HODO",
                         "phodo_cuts.param" if spec == "p" else "hhodo_cuts.param")
    text = read_file(path)
    if not text:
        print(f"WARNING: vanilla defaults not found at {path}", file=sys.stderr)
        return {}
    n_pmt = 21 if spec == "p" else 16
    expected = 4 * n_pmt
    pos_min = split_numbers(extract_block(text, "hodo_PosAdcTimeWindowMin"))
    pos_max = split_numbers(extract_block(text, "hodo_PosAdcTimeWindowMax"))
    neg_min = split_numbers(extract_block(text, "hodo_NegAdcTimeWindowMin"))
    neg_max = split_numbers(extract_block(text, "hodo_NegAdcTimeWindowMax"))
    if any(len(arr) != expected for arr in (pos_min, pos_max, neg_min, neg_max)):
        print(f"WARNING: {path} array length(s) != {expected} -- vanilla baseline unavailable for spec '{spec}'.",
              file=sys.stderr)
        return {}
    out = {}
    for plane in PLANES:
        pi = PLANE_IDX[plane]
        for ipmt in range(1, n_pmt + 1):
            idx = (ipmt - 1) * 4 + pi
            for side, lo_arr, hi_arr in (("Pos", pos_min, pos_max), ("Neg", neg_min, neg_max)):
                lo, hi = lo_arr[idx], hi_arr[idx]
                if lo == 0.0 and hi == 0.0:
                    continue
                out[(plane, side, ipmt)] = (lo, hi)
    return out


# ---- load + merge ----

def load_files(paths):
    for path in paths:
        try:
            with open(path) as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
            continue
        # This script's own outputs (hodo_diff_cuts_<v>.json and its
        # hodo_diff_cuts_<v>.json compatibility copy) now live inside
        # qa_dir by default, so a re-run's *.json glob will see them too.
        # Identify them by content (a "calibration_run" key, which only
        # this script's own output ever has) rather than by filename
        # pattern, so this is robust across any version number/renaming
        # -- otherwise a re-run would silently try to re-merge its own
        # prior merged output as if it were a fresh student submission.
        if "calibration_run" in data:
            print(f"NOTE: {path} looks like merge_hodo_qa.py's own output -- skipping to avoid double-counting.",
                  file=sys.stderr)
            continue
        channels = data.get("channels")
        if not channels:
            print(f"NOTE: {path} has no channels -- skipping.", file=sys.stderr)
            continue
        run = data.get("run")
        if run is None:
            print(f"WARNING: {path} has no \"run\" field -- skipping (can't classify calibration vs production).",
                  file=sys.stderr)
            continue
        yield path, run, channels


# ---- writing real hcana-format param files back out ----
# Same PMT-major/plane-minor flat-array layout the vanilla files use.
# Starts from the vanilla arrays (so untouched channels correctly keep
# their existing vanilla value rather than getting zeroed out) and
# overlays every resolved merged channel on top at the correct index.

def n_pmt_for(spec):
    return 21 if spec == "p" else 16


def build_flat_arrays(spec, merged, vanilla):
    n_pmt = n_pmt_for(spec)
    size = 4 * n_pmt
    arrays = {
        ("Pos", "Min"): [0.0] * size, ("Pos", "Max"): [0.0] * size,
        ("Neg", "Min"): [0.0] * size, ("Neg", "Max"): [0.0] * size,
    }
    for (plane, side, ipmt), (lo, hi) in vanilla.items():
        pi = PLANE_IDX.get(plane)
        if pi is None:
            continue
        idx = (ipmt - 1) * 4 + pi
        if idx >= size:
            continue
        arrays[(side, "Min")][idx] = lo
        arrays[(side, "Max")][idx] = hi
    for (s, plane, side, ipmt), entry in merged.items():
        if s != spec:
            continue
        pi = PLANE_IDX.get(plane)
        if pi is None:
            continue
        idx = (ipmt - 1) * 4 + pi
        if idx >= size:
            continue
        arrays[(side, "Min")][idx] = entry["lo"]
        arrays[(side, "Max")][idx] = entry["hi"]
    return arrays


def fmt_val(v):
    return f"{v:.0f}" if float(v).is_integer() else f"{v:.2f}"


def fmt_block(values):
    rows = []
    for i in range(0, len(values), 4):
        rows.append(", ".join(fmt_val(v) for v in values[i:i + 4]))
    return ",\n\t\t\t\t\t  ".join(rows)


def write_param_file(path, prefix, spec_label, arrays, n_pmt, note_lines):
    with open(path, "w") as f:
        f.write(f"; {os.path.basename(path)} -- STAGED by merge_hodo_qa.py, review before copying into PARAM/\n")
        for note in note_lines:
            f.write(f"; {note}\n")
        f.write(f"; layout: PMT-major, plane-minor (row = PMT 1..{n_pmt}, columns = 1x,1y,2x,2y)\n\n")
        f.write(f"{prefix}hodo_PosAdcTimeWindowMin = {fmt_block(arrays[('Pos', 'Min')])}\n\n")
        f.write(f"{prefix}hodo_PosAdcTimeWindowMax = {fmt_block(arrays[('Pos', 'Max')])}\n\n")
        f.write(f"{prefix}hodo_NegAdcTimeWindowMin = {fmt_block(arrays[('Neg', 'Min')])}\n\n")
        f.write(f"{prefix}hodo_NegAdcTimeWindowMax = {fmt_block(arrays[('Neg', 'Max')])}\n")


def channel_to_json_line(ch):
    """Serializes one channel dict as a SINGLE physical line of JSON.
    Required for compatibility with hodo_timediff_cut_app.C's
    LoadJsonAsMap(), which is a naive line-scanner (not a real JSON
    parser) that expects spec/plane/side/ipmt/lo/hi all on the same line
    -- json.dump(..., indent=2)'s one-field-per-line style would make
    every channel silently fail to parse there."""
    parts = []
    for k, v in ch.items():
        parts.append(f"{json.dumps(k)}: {json.dumps(v)}")
    return "{" + ", ".join(parts) + "}"


def write_hybrid_json(path, output):
    """Metadata fields nicely indented; each entry in "channels" forced
    onto one line each (see channel_to_json_line). Still fully valid,
    re-parseable JSON -- verified against json.load in testing -- just
    formatted for compatibility with the C++ tool's line-based parser."""
    meta = {k: v for k, v in output.items() if k != "channels"}
    lines = ["{"]
    for k, v in meta.items():
        lines.append(f"  {json.dumps(k)}: {json.dumps(v)},")
    lines.append('  "channels": [')
    ch_lines = [f"    {channel_to_json_line(ch)}" for ch in output["channels"]]
    lines.append(",\n".join(ch_lines))
    lines.append("  ]")
    lines.append("}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("qa_dir", nargs="?", default="./hodo_diff_qa",
                     help="Directory containing hodo_diff_cuts_<run>.json files (default: ./hodo_diff_qa)")
    ap.add_argument("-v", "--version", type=int, default=0,
                     help="Version/label number for this merge (default: 0). Used to name the default "
                          "output files, and doubles as the pseudo run-number for a hodo_diff_cuts_<version>.json "
                          "compatibility copy -- see below.")
    ap.add_argument("-o", "--output", default=None,
                     help="Output JSON path (default: <qa_dir>/hodo_diff_cuts_<version>.json)")
    ap.add_argument("-r", "--recursive", action="store_true", help="Search qa_dir recursively for *.json")
    ap.add_argument("--pattern", default="*.json", help="Glob pattern for input files (default: *.json)")
    ap.add_argument("--param-dir", default="../../PARAM", help="Path to PARAM/ for the vanilla fallback (default: ../../PARAM)")
    ap.add_argument("--no-write-param", action="store_true",
                     help="Skip writing phodo_cuts.param/hhodo_cuts.param (JSON output only)")
    ap.add_argument("--phodo-out", default=None, help="Output path for the staged SHMS param file (default: alongside -o, phodo_cuts.param)")
    ap.add_argument("--hhodo-out", default=None, help="Output path for the staged HMS param file (default: alongside -o, hhodo_cuts.param)")
    ap.add_argument("--no-compat-copy", action="store_true",
                     help="Skip writing the hodo_diff_cuts_<version>.json compatibility copy "
                          "(the name hodo_timediff_cut_app.C's --reference-run actually looks for)")
    args = ap.parse_args()

    if args.output is None:
        args.output = os.path.join(args.qa_dir, f"hodo_diff_cuts_{args.version}.json")

    pattern = os.path.join(args.qa_dir, "**", args.pattern) if args.recursive else os.path.join(args.qa_dir, args.pattern)
    files = sorted(glob.glob(pattern, recursive=args.recursive))
    if not files:
        print(f"No files matched {pattern}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} file(s) matching {pattern}:")
    for fpath in files:
        print(f"  {fpath}")

    files_data = list(load_files(files))
    if not files_data:
        print("No usable *.json files found.", file=sys.stderr)
        sys.exit(1)

    runs_seen = sorted(set(run for _, run, _ in files_data))
    print(f"\nRuns present: {runs_seen}")
    print(f"Calibration run: {CALIBRATION_RUN}"
          f"{' (present)' if CALIBRATION_RUN in runs_seen else ' (NOT present in this directory!)'}")

    # ---- pass 1: collect valid (post-filter) contributions and a record
    # of everything that got ignored, per channel ----
    contributions = {}  # (spec,plane,side,ipmt) -> [(lo,hi,file,run), ...]
    ignored = []         # (key, file, run, lo, hi, reason)

    for path, run, channels in files_data:
        fname = os.path.basename(path)
        for ch in channels:
            spec, plane, side, ipmt = ch.get("spec"), ch.get("plane"), ch.get("side"), ch.get("ipmt")
            lo, hi = ch.get("lo"), ch.get("hi")
            if None in (spec, plane, side, ipmt, lo, hi):
                continue
            key = (spec, plane, side, ipmt)

            if is_dead_2y(spec, plane, side, ipmt):
                ignored.append((key, fname, run, lo, hi, "universally-dead SHMS 2Y quartz bar"))
                continue
            if is_1x2x_off_for_run(spec, plane, run, ipmt):
                bracket = bracket_for_run(run)
                reason = (f"PMT {ipmt} off for run {run}'s production bracket"
                          if bracket is not None else f"run {run} outside all known brackets -- unvalidated, ignored")
                ignored.append((key, fname, run, lo, hi, reason))
                continue

            contributions.setdefault(key, []).append((lo, hi, fname, run))

    # ---- pass 2: union-merge whatever survived, per channel ----
    merged = {}
    for key, contribs in contributions.items():
        lo = min(c[0] for c in contribs)
        hi = max(c[1] for c in contribs)
        merged[key] = {"lo": lo, "hi": hi, "source": "merged", "contributors": contribs}

    # ---- pass 3: vanilla fallback for anything that ended up with zero
    # valid contributions but was reported (and filtered out) somewhere ----
    vanilla_cache = {}

    def vanilla_lookup(spec, plane, side, ipmt):
        if spec not in vanilla_cache:
            vanilla_cache[spec] = load_vanilla_hodo(args.param_dir, spec)
        return vanilla_cache[spec].get((plane, side, ipmt))

    fell_back_to_vanilla = []
    unresolved = []
    ignored_keys = set(k for k, *_ in ignored)
    for key in sorted(ignored_keys - set(merged.keys())):
        spec, plane, side, ipmt = key
        vanilla = vanilla_lookup(spec, plane, side, ipmt)
        if vanilla is not None:
            lo, hi = vanilla
            merged[key] = {"lo": lo, "hi": hi, "source": "vanilla-fallback", "contributors": []}
            fell_back_to_vanilla.append(key)
        else:
            unresolved.append(key)

    # ---- write output ----
    out_channels = []
    for key, entry in sorted(merged.items()):
        spec, plane, side, ipmt = key
        out_channels.append({
            "spec": spec, "det": "hod", "plane": plane, "side": side, "ipmt": ipmt,
            "lo": entry["lo"], "hi": entry["hi"],
            "source": entry["source"],
            "n_contributors": len(entry["contributors"]),
            "contributor_files": sorted(set(c[2] for c in entry["contributors"])),
        })

    output = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "calibration_run": CALIBRATION_RUN,
        "runs_present": runs_seen,
        "source_files": [os.path.basename(f) for f in files],
        "n_source_files": len(files),
        "n_channels": len(out_channels),
        "n_ignored_contributions": len(ignored),
        "n_vanilla_fallback": len(fell_back_to_vanilla),
        "n_unresolved": len(unresolved),
        "channels": out_channels,
    }
    write_hybrid_json(args.output, output)

    print(f"\nWrote {args.output}: {len(out_channels)} channel(s) "
          f"({len(fell_back_to_vanilla)} via vanilla fallback, {len(unresolved)} unresolved)")

    # hodo_timediff_cut_app.C's --reference-run looks specifically for
    # "hodo_diff_cuts_<run>.json" -- write a same-content copy under that
    # exact name (using --version as the pseudo run-number) so this merge
    # is directly usable as `hodo_timediff_cut_app.C(<some_run>, <version>, true)`
    # without any manual renaming.
    if not args.no_compat_copy:
        compat_path = os.path.join(args.qa_dir, f"hodo_diff_cuts_{args.version}.json")
        write_hybrid_json(compat_path, output)
        print(f"Wrote {compat_path} (compatibility copy -- use referenceRun={args.version} "
              f"with hodo_timediff_cut_app.C)")

    # ---- write real hcana-format phodo_cuts.param / hhodo_cuts.param ----
    if not args.no_write_param:
        out_dir = os.path.dirname(os.path.abspath(args.output))
        phodo_path = args.phodo_out or os.path.join(out_dir, "phodo_cuts.param")
        hhodo_path = args.hhodo_out or os.path.join(out_dir, "hhodo_cuts.param")

        # make sure both specs' vanilla arrays are loaded, not just whichever
        # ones happened to be needed for a vanilla-fallback lookup above
        for spec in ("p", "h"):
            if spec not in vanilla_cache:
                vanilla_cache[spec] = load_vanilla_hodo(args.param_dir, spec)

        notes = [
            f"generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} from {len(files)} file(s) in {args.qa_dir}",
            f"runs merged: {runs_seen} (calibration run: {CALIBRATION_RUN})",
            "SHMS 2Y dead quartz bars and off-for-run SHMS 1X/2X PMTs were excluded per known dead-channel rules --",
            f"see {os.path.basename(args.output)} for the full per-channel audit trail (source/contributors/ignored list)",
        ]

        p_arrays = build_flat_arrays("p", merged, vanilla_cache["p"])
        write_param_file(phodo_path, "p", "SHMS", p_arrays, n_pmt_for("p"), notes)
        print(f"Wrote {phodo_path} (STAGED -- review before copying into PARAM/SHMS/HODO/)")

        h_arrays = build_flat_arrays("h", merged, vanilla_cache["h"])
        write_param_file(hhodo_path, "h", "HMS", h_arrays, n_pmt_for("h"), notes)
        print(f"Wrote {hhodo_path} (STAGED -- review before copying into PARAM/HMS/HODO/)")

    if ignored:
        print(f"\n{len(ignored)} contribution(s) were IGNORED (dead/off SHMS PMT rules):")
        for key, fname, run, lo, hi, reason in ignored[:80]:
            spec, plane, side, ipmt = key
            print(f"  {spec}.hod.{plane}.{side}[{ipmt}] in {fname} (run {run}): "
                  f"({lo:.1f},{hi:.1f}) -- {reason}")
        if len(ignored) > 80:
            print(f"  ... and {len(ignored) - 80} more")

    if unresolved:
        print(f"\n{len(unresolved)} channel(s) have NO value at all -- ignored by the dead-PMT rules, "
              f"AND no vanilla default was found either:")
        for spec, plane, side, ipmt in unresolved:
            print(f"  {spec}.hod.{plane}.{side}[{ipmt}]")


if __name__ == "__main__":
    main()
