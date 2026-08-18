// hodo_timediff_cut_app.C
//
// Click-to-cut tool for the Hodoscope AdcTdcDiffTime window cuts, i.e.
//   P.hod.[pl].GoodPosAdcTdcDiffTime / GoodNegAdcTdcDiffTime   (SHMS)
//   H.hod.[pl].GoodPosAdcTdcDiffTime / GoodNegAdcTdcDiffTime   (HMS)
// for pl = 1x, 1y, 2x, 2y. Ported from reftime_cut_app.C -- same argument
// signature, same two-tier reference logic, same Prev/Next/click/Use
// Reference/Pause/Save workflow.
//
// Usage (must run interactively -- NOT with -q, it waits for clicks):
//   root -l hodo_timediff_cut_app.C                                  // spec="p", run=26107
//   root -l 'hodo_timediff_cut_app.C("p", 26107)'                    // one run, no baseline reference
//   root -l 'hodo_timediff_cut_app.C("p", 26107, 26092)'             // + comparison to run 26092
//   root -l 'hodo_timediff_cut_app.C("both", 26107, 26092)'          // SHMS + HMS together
//
// Non-interactive (batch), 4th arg = true -- just renders the summary PDF
// from whatever's already on record and exits, no display needed:
//   root -l -b -q 'hodo_timediff_cut_app.C("p", 26107, 26092, true)'
//
// Controls (same as reftime_cut_app.C):
//   Click twice on the histogram: 1st click = lo, 2nd click = hi (orange).
//   << Prev / Next >>      -- step through (plane, side, PMT) channels
//   Reset clicks           -- clear the current channel's clicks
//   Use Reference           -- copy the baseline reference's window for
//                             this channel into the current cut
//   Pause (save)            -- save progress to disk, keep going
//   Save && Finish            -- write hodo_diff_cuts_<spec>_<run>.json
//                             under outDir, close
//   Keys: n/b/r/p/s do the same as the buttons above.
//
// --------------------------------------------------------------------------
// HISTOGRAM SOURCE (built live from the raw replay tree, not a pre-built
// "golden file" 2D histogram):
//
//   For each (spec, plane, side) we build the 2D hist ourselves with a
//   single TTree::Draw(), mirroring exactly the hcana histogram-definition
//   rule you gave:
//     TH2F phodo_1x_good_adctdc_diff_time_vs_pmt_pos '...' [I+1]
//          P.hod.1x.GoodPosAdcTdcDiffTime  13 0.5 13.5  400 -200 200
//   translates to
//     tree->Draw("P.hod.1x.GoodPosAdcTdcDiffTime:Iteration$+1"
//                ">>h2d_p_1x_Pos(13,0.5,13.5,400,-200,200)", "", "goff")
//   -- ROOT's Iteration$ is the native equivalent of hcana's def-file
//   [I+1] loop index, so X = PMT slot index (1-based), Y = that slot's
//   diff-time value. Per-PMT diff-time distributions then come from
//   TH2::ProjectionY(ipmt, ipmt), same as before.
//
//   PMT COUNT PER PLANE IS AUTO-DETECTED, NOT HARDCODED: your example
//   shows SHMS 1x has 13 PMTs (not the 21 that timing_window_setup.C's
//   fHodoScin=nPlanes*nPMT padding constant would suggest -- that 21 is
//   only the max, used to pad hcana's own flat-array param file to a
//   single common size across all 4 planes). Since real per-plane counts
//   differ, DetectNPmt() reads the branch's own array length straight off
//   the tree (TLeaf::GetLenStatic(), falling back to GetLen() after
//   GetEntry(0) for leafcount-declared arrays) before building each 2D
//   hist, and the channel list itself is only built AFTER the tree is
//   open, from these detected per-plane counts.
//
// --------------------------------------------------------------------------
// OUTPUT: Save & Finish writes a plain JSON file per spec,
//   <outDir>/hodo_diff_cuts_<spec>_<run>.json
// with one self-describing object per cut channel (plane/side/ipmt/lo/hi,
// plus ref_lo/ref_hi from the baseline reference when available) -- not
// hcana's own flat-array param format. See a prior revision's comments (or
// just read ToJson() below) if you want the exact shape.
//
// Ported from reftime_cut_app.C: same (channel/spec, run, referenceRun,
// nonInteractive, rootFile, outDir, paramDir, ...) argument order; same
// baseline_ref (gray dashed line, comparison only -- --reference-run's
// staged JSON, or hcana's own vanilla PARAM/{SHMS,HMS}/HODO/
// {p,h}hodo_cuts.param read in ITS native flat, PMT-major/plane-minor
// array layout since that file isn't ours to change) vs. this_run_ref
// (silently seeds gResults from this run's own previously-staged JSON).
//
// Simplifications relative to reftime_cut_app.C (say the word if you want
// any of these ported over too): no cross-run "most conservative wins"
// merge (session's clicks just win outright over a prior save on
// conflict), no retained-fraction QA pass.
//
// MULTIPLICITY-CUT OVERLAY: like reftime_cut_app.C's "NO Cut" (blue) vs.
// "Multiplicity Cut == N" (red) overlay, each channel now draws two
// histograms -- the base (no cut, transparent blue outline) and a red
// overlay built from the same branch with a cut on the parallel
// Good[side]AdcMult array (Good[side]AdcMult == multCutValue, default 1),
// via the same Iteration$+1 trick since both arrays share the same
// per-PMT-slot indexing. If a plane/side has no companion Mult branch, the
// overlay is silently skipped for it (base histogram still shown alone).
// --------------------------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TControlBar.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TDatime.h>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cctype>
#include <cstdio>

// ==========================================================================
// Channel definition -- built AFTER the tree is open, from detected nPMT
// ==========================================================================

struct HodoChannel {
  TString spec;     // "p" or "h"
  TString plane;    // "1x","1y","2x","2y"
  int     planeIdx; // 0..3, matches offset in hcana's own padded param file
  TString side;     // "Pos" or "Neg"
  int     ipmt;     // 1-based PMT number
};

std::vector<TString> Planes() { return {"1x", "1y", "2x", "2y"}; }

// max PMT count used ONLY for reading hcana's own vanilla param file, which
// pads every plane's array out to a single common size (21 for SHMS, 16 for
// HMS -- see timing_window_setup.C's fHodoScin=nPlanes*nPMT comment). Real
// per-plane counts (used for everything else) come from DetectNPmt().
int PaddedNPmtFor(const TString& spec) { return spec == "p" ? 21 : 16; }

TString ChanLabel(const HodoChannel& c) {
  return Form("%s.hod.%s.Good%sAdcTdcDiffTime[%d]",
              c.spec == "p" ? "P" : "H", c.plane.Data(), c.side.Data(), c.ipmt);
}

TString BranchName(const TString& spec, const TString& plane, const TString& side) {
  return Form("%s.hod.%s.Good%sAdcTdcDiffTime", spec == "p" ? "P" : "H", plane.Data(), side.Data());
}

// flat-array index into hcana's own e.g. phodo_PosAdcTimeWindowMin -- PMT-
// major, plane-minor, using the PADDED max nPMT, matching
// calc_timing_windows()'s minArr[ipmt*4+offset]. Only used to read the
// vanilla baseline file.
int PaddedFlatIndex(const HodoChannel& c) { return (c.ipmt - 1) * 4 + c.planeIdx; }

// ==========================================================================
// Build the 2D vs-PMT histogram live from the tree, per your rule:
//   TH2F <name> '<title>' [I+1] <branch> nbinsX 0.5 nPmt+0.5 nbinsY ymin ymax
// i.e. tree->Draw("<branch>:Iteration$+1>>h2d(...)", "", "goff")
// ==========================================================================

static TTree* gTree = nullptr;
static int gDiffNbinsY = 400;
static double gDiffYmin = -200, gDiffYmax = 200;

// per-plane/side PMT count, auto-detected once and cached
static std::map<TString, int> gNPmtCache;   // key = "<spec>_<plane>_<side>"

// GoodPosAdcTdcDiffTime is a Double_t[N] branch (one slot per physical
// PMT, N constant across events), declared via a leaflist/leafcount rather
// than a std::vector -- so SetBranchAddress needs a bare double buffer,
// not vector<double>*. TLeaf::GetLen() (after GetEntry so any leafcount is
// resolved for that entry) gives N directly.
int DetectNPmt(const TString& spec, const TString& plane, const TString& side) {
  TString key = Form("%s_%s_%s", spec.Data(), plane.Data(), side.Data());
  auto it = gNPmtCache.find(key);
  if (it != gNPmtCache.end()) return it->second;

  int n = 0;
  TString branch = BranchName(spec, plane, side);
  TLeaf* leaf = gTree ? gTree->GetLeaf(branch) : nullptr;
  if (leaf && gTree->GetEntries() > 0) {
    gTree->GetEntry(0);
    n = leaf->GetLen();
  }
  if (n <= 0) printf("WARNING: could not detect PMT count for branch '%s' -- skipping this plane/side.\n",
                      branch.Data());
  gNPmtCache[key] = n;
  return n;
}

std::vector<HodoChannel> BuildChannelsFromTree(const TString& spec) {
  std::vector<HodoChannel> out;
  int planeIdx = 0;
  for (auto& pl : Planes()) {
    for (auto& side : { TString("Pos"), TString("Neg") }) {
      int nPmt = DetectNPmt(spec, pl, side);
      for (int ipmt = 1; ipmt <= nPmt; ++ipmt) out.push_back({ spec, pl, planeIdx, side, ipmt });
    }
    ++planeIdx;
  }
  return out;
}

TString Hist2DTag(const TString& spec, const TString& plane, const TString& side) {
  return Form("h2d_%s_%s_%s", spec.Data(), plane.Data(), side.Data());
}

TString MultBranchName(const TString& spec, const TString& plane, const TString& side) {
  return Form("%s.hod.%s.Good%sAdcMult", spec == "p" ? "P" : "H", plane.Data(), side.Data());
}

static std::map<TString, TH2F*> gHist2DCache;         // no-cut 2D hist, one build per (spec,plane,side)
static std::map<TString, TH2F*> gHist2DMultCache;     // multiplicity-cut 2D hist, same keying
static std::map<TString, bool>  gHasMultBranch;        // whether the companion Mult branch exists at all
static int gMultCutValue = 1;                          // "Good...AdcMult == this" -- overlay selection

TH2F* GetHist2D(const HodoChannel& c) {
  TString tag = Hist2DTag(c.spec, c.plane, c.side);
  auto it = gHist2DCache.find(tag);
  if (it != gHist2DCache.end()) return it->second;

  int nPmt = DetectNPmt(c.spec, c.plane, c.side);
  TH2F* h = nullptr;
  if (gTree && nPmt > 0) {
    TString branch = BranchName(c.spec, c.plane, c.side);
    TString expr = Form("%s:Iteration$+1>>%s(%d,0.5,%.1f,%d,%.1f,%.1f)",
                         branch.Data(), tag.Data(), nPmt, nPmt + 0.5, gDiffNbinsY, gDiffYmin, gDiffYmax);
    gTree->Draw(expr, "", "goff");
    h = (TH2F*)gDirectory->Get(tag);
    if (h) { h->SetDirectory(nullptr); h->SetStats(kFALSE); }
  }
  gHist2DCache[tag] = h;  // cache the miss too
  if (!h) printf("WARNING: could not build 2D hist for %s.hod.%s.Good%sAdcTdcDiffTime\n",
                  c.spec.Data(), c.plane.Data(), c.side.Data());
  return h;
}

// Same as GetHist2D, but with a cut on the parallel Good[side]AdcMult array
// (Iteration$ synchronizes across both arrays since they're the same
// length -- one multiplicity value per PMT slot, same indexing as the
// diff-time array). Returns nullptr if the companion Mult branch doesn't
// exist for this plane/side, in which case the overlay is just skipped.
TH2F* GetHist2DMultCut(const HodoChannel& c) {
  TString tag = Hist2DTag(c.spec, c.plane, c.side) + "_multcut";
  auto it = gHist2DMultCache.find(tag);
  if (it != gHist2DMultCache.end()) return it->second;

  TString multKey = Form("%s_%s_%s", c.spec.Data(), c.plane.Data(), c.side.Data());
  TString multBranch = MultBranchName(c.spec, c.plane, c.side);
  bool hasMult = gTree && gTree->GetBranch(multBranch) != nullptr;
  gHasMultBranch[multKey] = hasMult;

  int nPmt = DetectNPmt(c.spec, c.plane, c.side);
  TH2F* h = nullptr;
  if (hasMult && nPmt > 0) {
    TString diffBranch = BranchName(c.spec, c.plane, c.side);
    TString expr = Form("%s:Iteration$+1>>%s(%d,0.5,%.1f,%d,%.1f,%.1f)",
                         diffBranch.Data(), tag.Data(), nPmt, nPmt + 0.5, gDiffNbinsY, gDiffYmin, gDiffYmax);
    TString cutExpr = Form("%s==%d", multBranch.Data(), gMultCutValue);
    gTree->Draw(expr, cutExpr, "goff");
    h = (TH2F*)gDirectory->Get(tag);
    if (h) { h->SetDirectory(nullptr); h->SetStats(kFALSE); }
  } else if (!hasMult) {
    printf("NOTE: no '%s' branch -- multiplicity-cut overlay skipped for %s.hod.%s.Good%s*\n",
           multBranch.Data(), c.spec.Data(), c.plane.Data(), c.side.Data());
  }
  gHist2DMultCache[tag] = h;
  return h;
}

// fresh 1D projection for one PMT -- caller owns the returned histogram.
// Transparent outline (no fill) so a multiplicity-cut overlay drawn on top
// stays visible instead of being hidden underneath a solid fill.
TH1D* GetOrBuildProjection(const HodoChannel& c) {
  TH2F* h2 = GetHist2D(c);
  if (!h2) return nullptr;
  if (c.ipmt < 1 || c.ipmt > h2->GetNbinsX()) return nullptr;
  TString tag = Form("h_%s_%s_%s_pmt%d_proj", c.spec.Data(), c.plane.Data(), c.side.Data(), c.ipmt);
  TH1D* h1 = (TH1D*)h2->ProjectionY(tag, c.ipmt, c.ipmt);
  h1->SetDirectory(nullptr);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2);
  h1->SetFillStyle(0);  // transparent -- was SetFillColorAlpha(kBlue, 0.3)
  h1->SetStats(kFALSE);
  return h1;
}

// same idea, but from the multiplicity-cut 2D hist. Returns nullptr if no
// companion Mult branch was available for this plane/side.
TH1D* GetOrBuildProjectionMultCut(const HodoChannel& c) {
  TH2F* h2 = GetHist2DMultCut(c);
  if (!h2) return nullptr;
  if (c.ipmt < 1 || c.ipmt > h2->GetNbinsX()) return nullptr;
  TString tag = Form("h_%s_%s_%s_pmt%d_multcut_proj", c.spec.Data(), c.plane.Data(), c.side.Data(), c.ipmt);
  TH1D* h1 = (TH1D*)h2->ProjectionY(tag, c.ipmt, c.ipmt);
  h1->SetDirectory(nullptr);
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->SetFillStyle(0);
  h1->SetStats(kFALSE);
  return h1;
}

// ==========================================================================
// Minimal, regex-free text helpers (same approach as reftime_cut_app.C)
// ==========================================================================

TString ReadFile(const TString& path) {
  std::ifstream in(path.Data());
  if (!in.is_open()) return "";
  std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  return TString(content.c_str());
}

TString NowString() { TDatime now; return now.AsSQLString(); }

// ---- vanilla hcana .param parsing (only for the un-touched PARAM/ defaults) ----

std::vector<double> SplitNumbers(const TString& s) {
  std::vector<double> out;
  TObjArray* arr = s.Tokenize(", \t\r\n");
  for (int i = 0; i < arr->GetEntries(); ++i) {
    TString tok = ((TObjString*)arr->At(i))->GetString();
    if (tok.IsFloat()) out.push_back(tok.Atof());
  }
  delete arr;
  return out;
}

TString ExtractBlock(const TString& text, const TString& key) {
  Ssiz_t start = text.Index(key);
  if (start == kNPOS) return "";
  Ssiz_t eqPos = text.Index("=", start);
  if (eqPos == kNPOS) return "";
  Ssiz_t blockStart = eqPos + 1;
  Ssiz_t blankPos = text.Index("\n\n", blockStart);
  Ssiz_t end = (blankPos == kNPOS) ? text.Length() : blankPos;
  return text(blockStart, end - blockStart);
}

// ---- shared in-memory result type: channel key -> cut window ----

struct ChannelState { double lo = 0, hi = 0; bool hasLo = false, hasHi = false; };
using HodoResultsMap = std::map<TString, ChannelState>;  // key = ChanLabel()

// Loads hcana's native flat-array hhodo_cuts.param / phodo_cuts.param
// (PMT-major/plane-minor, padded to PaddedNPmtFor(spec)) and converts it
// into a HodoResultsMap restricted to the REAL channels this run actually
// has (gChannels), dropping any padded-but-nonexistent higher-PMT slots.
HodoResultsMap LoadVanillaParamAsMap(const TString& spec, const TString& path,
                                      const std::vector<HodoChannel>& realChannels) {
  HodoResultsMap out;
  TString text = ReadFile(path);
  if (text.Length() == 0) {
    printf("No vanilla defaults found at %s.\n", path.Data());
    return out;
  }
  std::vector<double> posMin = SplitNumbers(ExtractBlock(text, "hodo_PosAdcTimeWindowMin"));
  std::vector<double> posMax = SplitNumbers(ExtractBlock(text, "hodo_PosAdcTimeWindowMax"));
  std::vector<double> negMin = SplitNumbers(ExtractBlock(text, "hodo_NegAdcTimeWindowMin"));
  std::vector<double> negMax = SplitNumbers(ExtractBlock(text, "hodo_NegAdcTimeWindowMax"));
  int expected = 4 * PaddedNPmtFor(spec);
  if ((int)posMin.size() != expected || (int)posMax.size() != expected ||
      (int)negMin.size() != expected || (int)negMax.size() != expected) {
    printf("NOTE: %s array length(s) != %d (nPlanes*paddedNPMT) -- treating vanilla baseline as unavailable.\n",
           path.Data(), expected);
    return out;
  }
  int nSkippedOff = 0;
  for (auto& c : realChannels) {
    if (c.spec != spec) continue;
    int idx = PaddedFlatIndex(c);
    if (idx < 0 || idx >= expected) continue;
    ChannelState cs;
    if (c.side == "Pos") { cs.lo = posMin[idx]; cs.hi = posMax[idx]; }
    else                 { cs.lo = negMin[idx]; cs.hi = negMax[idx]; }
    // (0.00, 0.00) is hcana's convention for a disabled/turned-off paddle,
    // not a real (degenerate) cut window -- skip it so it doesn't draw a
    // misleading reference line at 0, and so "Use Reference" can't copy a
    // fake cut onto a channel that was never actually enabled.
    if (cs.lo == 0.0 && cs.hi == 0.0) { ++nSkippedOff; continue; }
    cs.hasLo = cs.hasHi = true;
    out[ChanLabel(c)] = cs;
  }
  if (nSkippedOff) printf("Skipped %d channel(s) with a (0.00, 0.00) window in %s (treated as disabled paddles).\n",
                           nSkippedOff, path.Data());
  return out;
}

// ---- our own JSON output (write + read-back) ----

TString ToJson(const TString& spec, int run, const std::vector<HodoChannel>& realChannels,
                const HodoResultsMap& merged, const HodoResultsMap& baseline) {
  TString out;
  out += "{\n";
  out += Form("  \"run\": %d,\n", run);
  out += Form("  \"spec\": \"%s\",\n", spec.Data());
  out += Form("  \"generated\": \"%s\",\n", NowString().Data());
  out += "  \"channels\": [\n";
  bool first = true;
  for (auto& c : realChannels) {
    if (c.spec != spec) continue;
    TString key = ChanLabel(c);
    auto it = merged.find(key);
    if (it == merged.end() || !it->second.hasLo || !it->second.hasHi) continue;
    if (!first) out += ",\n";
    first = false;
    out += Form("    {\"plane\": \"%s\", \"side\": \"%s\", \"ipmt\": %d, \"lo\": %.3f, \"hi\": %.3f",
                 c.plane.Data(), c.side.Data(), c.ipmt, it->second.lo, it->second.hi);
    auto rit = baseline.find(key);
    if (rit != baseline.end() && rit->second.hasLo && rit->second.hasHi) {
      out += Form(", \"ref_lo\": %.3f, \"ref_hi\": %.3f", rit->second.lo, rit->second.hi);
    }
    out += "}";
  }
  out += "\n  ]\n}\n";
  return out;
}

TString ExtractJsonString(const TString& line, const TString& key) {
  TString pat = TString("\"") + key + "\"";
  Ssiz_t k = line.Index(pat);
  if (k == kNPOS) return "";
  Ssiz_t colon = line.Index(":", k);
  if (colon == kNPOS) return "";
  Ssiz_t q1 = line.Index("\"", colon + 1);
  if (q1 == kNPOS) return "";
  Ssiz_t q2 = line.Index("\"", q1 + 1);
  if (q2 == kNPOS) return "";
  return line(q1 + 1, q2 - q1 - 1);
}

bool ExtractJsonNumber(const TString& line, const TString& key, double& val) {
  TString pat = TString("\"") + key + "\"";
  Ssiz_t k = line.Index(pat);
  if (k == kNPOS) return false;
  Ssiz_t colon = line.Index(":", k);
  if (colon == kNPOS) return false;
  Ssiz_t start = colon + 1;
  while (start < line.Length() && line[start] == ' ') ++start;
  Ssiz_t end = start;
  while (end < line.Length() && (isdigit((unsigned char)line[end]) || line[end] == '-' || line[end] == '.'))
    ++end;
  if (end <= start) return false;
  TString numStr = line(start, end - start);
  if (!numStr.IsFloat()) return false;
  val = numStr.Atof();
  return true;
}

// Parses back exactly the shape ToJson() writes: one channel object per
// line inside "channels". Not a general JSON parser -- tailored to our own
// output, same hand-rolled approach reftime_cut_app.C uses for .param files.
HodoResultsMap LoadJsonAsMap(const TString& spec, const TString& path) {
  HodoResultsMap out;
  TString text = ReadFile(path);
  if (text.Length() == 0) return out;
  std::stringstream ss(text.Data());
  std::string lineStd;
  while (std::getline(ss, lineStd)) {
    TString line(lineStd.c_str());
    if (line.Index("\"plane\"") == kNPOS) continue;
    TString plane = ExtractJsonString(line, "plane");
    TString side = ExtractJsonString(line, "side");
    double ipmtD, lo, hi;
    if (plane.Length() == 0 || side.Length() == 0) continue;
    if (!ExtractJsonNumber(line, "ipmt", ipmtD)) continue;
    if (!ExtractJsonNumber(line, "lo", lo)) continue;
    if (!ExtractJsonNumber(line, "hi", hi)) continue;
    HodoChannel c{ spec, plane, 0, side, (int)ipmtD };
    ChannelState cs; cs.lo = lo; cs.hi = hi; cs.hasLo = cs.hasHi = true;
    out[ChanLabel(c)] = cs;
  }
  return out;
}

// ==========================================================================
// Interactive state
// ==========================================================================

static std::vector<HodoChannel> gChannels;
static int gIdx = 0;
static int gRun = 0;
static TString gOutDir, gParamDir, gRootPath;

static HodoResultsMap gResults;   // this session's clicks
static std::vector<double> gClicks;
static bool gFreshVisit = true;
static std::vector<TObject*> gLines, gRefLines;
static TCanvas* gCanvas = nullptr;
static TControlBar* gControlBar = nullptr;
static TH1D* gHcur = nullptr;
static TH1D* gHcurMult = nullptr;  // multiplicity-cut overlay for the current channel, may be null
static TString gProgressPath;

static std::map<TString, HodoResultsMap> gBaselineRef;  // gray dashed comparison line, per spec
static std::map<TString, HodoResultsMap> gThisRunRef;   // seeds gResults, per spec

TString CurKey() { return ChanLabel(gChannels[gIdx]); }

void SaveProgress() {
  if (gProgressPath.Length() == 0) return;
  std::ofstream out(gProgressPath.Data());
  for (const auto& kv : gResults)
    out << kv.first << "|" << (kv.second.hasLo ? Form("%.4f", kv.second.lo) : "") << "|"
        << (kv.second.hasHi ? Form("%.4f", kv.second.hi) : "") << "\n";
}

int LoadProgress() {
  TString text = ReadFile(gProgressPath);
  if (text.Length() == 0) return 0;
  int n = 0;
  std::stringstream ss(text.Data());
  std::string lineStd;
  while (std::getline(ss, lineStd)) {
    TObjArray* parts = TString(lineStd.c_str()).Tokenize("|");
    if (parts->GetEntries() < 1) { delete parts; continue; }
    TString key = ((TObjString*)parts->At(0))->GetString();
    ChannelState cs;
    if (parts->GetEntries() > 1) {
      TString loStr = ((TObjString*)parts->At(1))->GetString();
      if (loStr.IsFloat()) { cs.lo = loStr.Atof(); cs.hasLo = true; }
    }
    if (parts->GetEntries() > 2) {
      TString hiStr = ((TObjString*)parts->At(2))->GetString();
      if (hiStr.IsFloat()) { cs.hi = hiStr.Atof(); cs.hasHi = true; }
    }
    if (cs.hasLo || cs.hasHi) { gResults[key] = cs; ++n; }
    delete parts;
  }
  return n;
}

bool GetWindowFrom(const std::map<TString, HodoResultsMap>& src, const HodoChannel& c,
                    double& lo, double& hi) {
  auto specIt = src.find(c.spec.Data());
  if (specIt == src.end()) return false;
  auto it = specIt->second.find(ChanLabel(c));
  if (it == specIt->second.end() || !it->second.hasLo || !it->second.hasHi) return false;
  lo = it->second.lo; hi = it->second.hi;
  return true;
}

bool GetBaselineWindow(const HodoChannel& c, double& lo, double& hi) { return GetWindowFrom(gBaselineRef, c, lo, hi); }

// ---- drawing ----

void RedrawCutLines() {
  for (auto* o : gLines) delete o;
  gLines.clear();
  if (!gHcur) return;
  double ymax = gHcur->GetMaximum() * 2;
  for (double x : gClicks) {
    TLine* ln = new TLine(x, 0, x, ymax);
    ln->SetLineColor(kOrange + 1);
    ln->SetLineWidth(2);
    ln->Draw();
    gLines.push_back(ln);
  }
  gPad->Modified();
  gPad->Update();
}

void DrawReferenceLines() {
  for (auto* o : gRefLines) delete o;
  gRefLines.clear();
  if (!gHcur) return;
  double lo, hi;
  if (!GetBaselineWindow(gChannels[gIdx], lo, hi)) return;
  double ymax = gHcur->GetMaximum() * 2;
  for (double x : { lo, hi }) {
    TLine* ln = new TLine(x, 0, x, ymax);
    ln->SetLineColor(kGray + 2);
    ln->SetLineStyle(2);
    ln->Draw();
    gRefLines.push_back(ln);
  }
}

void UpdateLegend() {
  TLegend* leg = new TLegend(0.60, 0.62, 0.99, 0.92);
  leg->SetTextSize(0.024);
  leg->AddEntry(gHcur, "no cut", "l");
  if (gHcurMult) leg->AddEntry(gHcurMult, Form("mult == %d", gMultCutValue), "l");
  double refLo, refHi;
  if (GetBaselineWindow(gChannels[gIdx], refLo, refHi)) {
    TLine* proxy = new TLine(0, 0, 0, 0);
    proxy->SetLineColor(kGray + 2); proxy->SetLineStyle(2);
    leg->AddEntry(proxy, Form("cuts-ref (%.1f, %.1f)", refLo, refHi), "l");
  }
  if (!gClicks.empty()) {
    TString loS = Form("%.1f", gClicks[0]);
    TString hiS = gClicks.size() >= 2 ? Form("%.1f", gClicks[1]) : "-";
    TLine* proxy = new TLine(0, 0, 0, 0);
    proxy->SetLineColor(kOrange + 1);
    leg->AddEntry(proxy, Form("cuts-this run %d (%s, %s)", gRun, loS.Data(), hiS.Data()), "l");
  }
  leg->Draw();
  gPad->Modified();
  gPad->Update();
}

void StoreCurrent() {
  TString key = CurKey();
  if (gClicks.size() == 2) {
    ChannelState cs;
    cs.lo = std::min(gClicks[0], gClicks[1]);
    cs.hi = std::max(gClicks[0], gClicks[1]);
    cs.hasLo = cs.hasHi = true;
    gResults[key] = cs;
  } else if (gClicks.empty()) {
    gResults.erase(key);
  }
  SaveProgress();
}

void DrawChannel() {
  const HodoChannel& c = gChannels[gIdx];
  gPad->Clear();
  gPad->SetRightMargin(0.02);

  if (gHcur) { delete gHcur; gHcur = nullptr; }
  if (gHcurMult) { delete gHcurMult; gHcurMult = nullptr; }
  gHcur = GetOrBuildProjection(c);
  if (!gHcur) {
    TLatex* msg = new TLatex(0.5, 0.5, Form("%s -- could not build histogram", ChanLabel(c).Data()));
    msg->SetNDC(); msg->SetTextAlign(22); msg->SetTextSize(0.03);
    msg->Draw();
    gClicks.clear();
    gPad->Modified(); gPad->Update();
    return;
  }
  gHcur->SetTitle(Form("[%d/%d] %s;AdcTdcDiffTime (channels);Counts",
                        gIdx + 1, (int)gChannels.size(), ChanLabel(c).Data()));
  gHcur->Draw("HIST");
  gHcurMult = GetOrBuildProjectionMultCut(c);
  if (gHcurMult) gHcurMult->Draw("HIST SAME");
  gPad->SetLogy();

  auto it = gResults.find(CurKey());
  if (it != gResults.end() && it->second.hasLo) {
    gClicks.clear();
    gClicks.push_back(it->second.lo);
    if (it->second.hasHi) gClicks.push_back(it->second.hi);
  } else {
    gClicks.clear();
  }
  gFreshVisit = true;

  DrawReferenceLines();
  RedrawCutLines();
  UpdateLegend();
}

void DrawChannelStatic(TVirtualPad* pad, const HodoChannel& c) {
  pad->cd();
  pad->SetLogy();
  pad->SetLeftMargin(0.12); pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.16); pad->SetBottomMargin(0.12);

  TH1D* h = GetOrBuildProjection(c);
  if (!h) {
    TLatex* msg = new TLatex(0.5, 0.5, "missing");
    msg->SetNDC(); msg->SetTextAlign(22); msg->SetTextSize(0.09);
    msg->SetTextColor(kGray + 2);
    msg->Draw();
    return;
  }
  h->SetTitle("");
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->Draw("HIST");
  TH1D* hMult = GetOrBuildProjectionMultCut(c);
  if (hMult) hMult->Draw("HIST SAME");

  double refLo, refHi;
  bool hasRef = GetBaselineWindow(c, refLo, refHi);
  double lo = 0, hi = 0; bool hasLo = false, hasHi = false;
  auto it = gResults.find(ChanLabel(c));
  if (it != gResults.end()) { hasLo = it->second.hasLo; lo = it->second.lo;
                               hasHi = it->second.hasHi; hi = it->second.hi; }
  double ymax = h->GetMaximum();
  if (hasRef) {
    for (double x : { refLo, refHi }) {
      TLine* l = new TLine(x, 0, x, ymax * 1.8);
      l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->Draw();
    }
  }
  if (hasLo) { TLine* l = new TLine(lo, 0, lo, ymax * 1.8); l->SetLineColor(kOrange + 1); l->SetLineWidth(2); l->Draw(); }
  if (hasHi) { TLine* l = new TLine(hi, 0, hi, ymax * 1.8); l->SetLineColor(kOrange + 1); l->SetLineWidth(2); l->Draw(); }

  TLatex* title = new TLatex(0.02, 0.90, Form("%s.%s.%s[%d]", c.spec.Data(), c.plane.Data(), c.side.Data(), c.ipmt));
  title->SetNDC(); title->SetTextSize(0.09); title->Draw();
  delete h;
  if (hMult) delete hMult;
}

// ---- interaction ----

void UseReference() {
  double lo, hi;
  if (!GetBaselineWindow(gChannels[gIdx], lo, hi)) {
    printf("No reference available for %s.\n", CurKey().Data());
    return;
  }
  gClicks = { lo, hi };
  gFreshVisit = false;
  RedrawCutLines();
  UpdateLegend();
  StoreCurrent();
  printf("[Use Reference] %s -> lo=%.2f hi=%.2f\n", CurKey().Data(), lo, hi);
}

void GoNext() { StoreCurrent(); if (gIdx < (int)gChannels.size() - 1) { ++gIdx; DrawChannel(); } }
void GoPrev() { StoreCurrent(); if (gIdx > 0) { --gIdx; DrawChannel(); } }
void PauseSave() { StoreCurrent(); printf("[progress saved]\n"); }

void WriteJsonResults() {
  std::set<TString> specsPresent;
  for (auto& c : gChannels) specsPresent.insert(c.spec);

  for (auto& s : specsPresent) {
    HodoResultsMap merged = gThisRunRef.count(s.Data()) ? gThisRunRef[s.Data()] : HodoResultsMap();
    for (auto& c : gChannels) {
      if (c.spec != s) continue;
      auto it = gResults.find(ChanLabel(c));
      if (it == gResults.end() || !it->second.hasLo || !it->second.hasHi) continue;
      merged[ChanLabel(c)] = it->second;
    }
    if (merged.empty()) { printf("Nothing set for spec '%s' -- skipping.\n", s.Data()); continue; }

    TString path = Form("%s/hodo_diff_cuts_%s_%d.json", gOutDir.Data(), s.Data(), gRun);
    HodoResultsMap baseline = gBaselineRef.count(s.Data()) ? gBaselineRef[s.Data()] : HodoResultsMap();
    std::ofstream(path.Data()) << ToJson(s, gRun, gChannels, merged, baseline);
    printf("Wrote %s (%d channel(s))\n", path.Data(), (int)merged.size());
  }
}

void SaveAndFinish() {
  StoreCurrent();
  SaveProgress();
  WriteJsonResults();
  if (gCanvas) gCanvas->Close();
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "Hodo DiffTime Cut Controls", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Previous PMT/plane/side channel");
  gControlBar->AddButton("Next >>", "GoNext();", "Next PMT/plane/side channel");
  gControlBar->AddButton("Reset clicks", "gClicks.clear(); RedrawCutLines(); UpdateLegend();",
                          "Clear this channel's clicks");
  gControlBar->AddButton("Use Reference", "UseReference();",
                          "Copy the baseline reference's window into this cut");
  gControlBar->AddButton("Pause (save)", "PauseSave();", "Save progress without finishing");
  gControlBar->AddButton("Save && Finish", "SaveAndFinish();",
                          "Write hodo_diff_cuts_<spec>_<run>.json under outDir and close");
  gControlBar->Show();
}

void OnClick() {
  int event = gPad->GetEvent();
  if (event == kKeyPress) {
    int key = gPad->GetEventX();
    if (key == 'n' || key == 'N') GoNext();
    else if (key == 'b' || key == 'B') GoPrev();
    else if (key == 'r' || key == 'R') UseReference();
    else if (key == 'p' || key == 'P') PauseSave();
    else if (key == 's' || key == 'S') SaveAndFinish();
    return;
  }
  if (event != kButton1Down || !gHcur) return;

  double x = gPad->AbsPixeltoX(gPad->GetEventX());
  if (gFreshVisit) { gClicks.clear(); gClicks.push_back(x); gFreshVisit = false; }
  else if (gClicks.size() >= 2) { gClicks.clear(); gClicks.push_back(x); }
  else { gClicks.push_back(x); }

  RedrawCutLines();
  UpdateLegend();
  if (gClicks.size() == 1) printf("lo = %.2f\n", gClicks[0]);
  else printf("lo = %.2f   hi = %.2f\n", std::min(gClicks[0], gClicks[1]), std::max(gClicks[0], gClicks[1]));
  if (gClicks.size() == 2) StoreCurrent();
}

// ==========================================================================
// Main entry point -- argument order matches reftime_cut_app.C
// ==========================================================================

void hodo_timediff_cut_app(
    const char* spec = "p",            // "p", "h", or "both" -- like reftime_cut_app.C's channel="all"
    int run = 26107,
    int referenceRun = -1,
    bool nonInteractive = false,
    const char* rootFile = nullptr,    // the run's own replay ROOT file (has the raw branches, not a golden 2D hist)
    const char* outDir = "./hodo_diff_qa",
    const char* paramDir = "../../PARAM",
    int gridCols = 4,
    int gridRows = 4,
    const char* treeName = "T",
    int diffNbinsY = 400,               // -- these three match your booking rule's "400 -200 200"
    double diffYmin = -200,
    double diffYmax = 200,
    int multCutValue = 1)               // red overlay = Good...AdcMult == this value
{
  if (nonInteractive) gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  gHist2DCache.clear();
  gHist2DMultCache.clear();
  gHasMultBranch.clear();
  gNPmtCache.clear();
  gRun = run;
  gOutDir = outDir;
  gParamDir = paramDir;
  gDiffNbinsY = diffNbinsY; gDiffYmin = diffYmin; gDiffYmax = diffYmax;
  gMultCutValue = multCutValue;
  gSystem->mkdir(gOutDir, true);

  gRootPath = rootFile ? TString(rootFile)
      : Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/"
             "coin_replay_production_%d_2000000_0.root", run);
  TFile* f = TFile::Open(gRootPath);
  if (!f || f->IsZombie()) { printf("Could not open %s\n", gRootPath.Data()); return; }
  gTree = (TTree*)f->Get(treeName);
  if (!gTree) { printf("Tree '%s' not found in %s\n", treeName, gRootPath.Data()); return; }

  gChannels.clear();
  TString specStr(spec);
  if (specStr == "both") {
    for (auto& c : BuildChannelsFromTree("p")) gChannels.push_back(c);
    for (auto& c : BuildChannelsFromTree("h")) gChannels.push_back(c);
  } else {
    gChannels = BuildChannelsFromTree(specStr);
  }
  if (gChannels.empty()) { printf("No channels found -- check branch names/tree contents.\n"); return; }

  std::set<TString> specsPresent;
  for (auto& c : gChannels) specsPresent.insert(c.spec);

  // ---- baseline_ref (gray dashed line): --reference-run's JSON, or vanilla PARAM/ ----
  for (auto& s : specsPresent) {
    if (referenceRun >= 0) {
      TString rPath = Form("%s/hodo_diff_cuts_%s_%d.json", outDir, s.Data(), referenceRun);
      HodoResultsMap m = LoadJsonAsMap(s, rPath);
      if (m.empty()) printf("WARNING: --reference-run %d requested but %s missing/empty for spec '%s'.\n",
                             referenceRun, rPath.Data(), s.Data());
      gBaselineRef[s.Data()] = m;
    } else {
      TString vPath = s == "p" ? Form("%s/SHMS/HODO/phodo_cuts.param", paramDir)
                                : Form("%s/HMS/HODO/hhodo_cuts.param", paramDir);
      gBaselineRef[s.Data()] = LoadVanillaParamAsMap(s, vPath, gChannels);
    }
  }

  // ---- this_run_ref (seeds gResults silently): this run's own prior JSON save ----
  for (auto& s : specsPresent) {
    TString path = Form("%s/hodo_diff_cuts_%s_%d.json", outDir, s.Data(), run);
    gThisRunRef[s.Data()] = LoadJsonAsMap(s, path);
  }

  gProgressPath = Form("%s/%s_%d_progress.txt", outDir, specStr.Data(), run);
  int nLoaded = LoadProgress();
  if (nLoaded > 0) printf("Resuming saved progress from %s: %d channel(s) set.\n", gProgressPath.Data(), nLoaded);

  int nSeeded = 0;
  for (auto& c : gChannels) {
    TString key = ChanLabel(c);
    if (gResults.count(key)) continue;
    double lo, hi;
    if (!GetWindowFrom(gThisRunRef, c, lo, hi)) continue;
    ChannelState cs; cs.lo = lo; cs.hi = hi; cs.hasLo = cs.hasHi = true;
    gResults[key] = cs;
    ++nSeeded;
  }
  if (nSeeded) printf("Seeded %d channel(s) from this run's own previously-saved JSON.\n", nSeeded);

  gIdx = 0;
  gCanvas = new TCanvas("c1", Form("hodo diff-time cuts, run %d (%s)", run, spec), 1000, 650);
  DrawChannel();

  if (nonInteractive) {
    TString pdfPath = Form("%s/%d_%s_all_summary.pdf", outDir, run, specStr.Data());
    int perPage = gridCols * gridRows;
    int nPages = ((int)gChannels.size() + perPage - 1) / perPage;
    gCanvas->SetCanvasSize(gridCols * 380, gridRows * 300);
    gCanvas->Print(pdfPath + "[");
    for (int page = 0; page < nPages; ++page) {
      gCanvas->Clear();
      gCanvas->Divide(gridCols, gridRows, 0.001, 0.001);
      for (int cell = 0; cell < perPage; ++cell) {
        int idx = page * perPage + cell;
        if (idx >= (int)gChannels.size()) break;
        DrawChannelStatic(gCanvas->cd(cell + 1), gChannels[idx]);
      }
      gCanvas->Print(pdfPath);
    }
    gCanvas->Print(pdfPath + "]");
    printf("\nWrote %s (%d channels, %dx%d grid, %d pages -- inspection only, no .json written)\n",
           pdfPath.Data(), (int)gChannels.size(), gridCols, gridRows, nPages);
    delete gCanvas;
    return;
  }

  gCanvas->AddExec("dynamic", "OnClick()");
  gCanvas->Update();
  MakeControlBar();

  printf("\n=== run %d -- %d channel(s) loaded (spec=%s) ===\n", run, (int)gChannels.size(), spec);
  printf("Click on the histogram: 1st = lo, 2nd = hi (orange).\n");
  printf("Use the 'Hodo DiffTime Cut Controls' panel, or keys n/b/r/p/s.\n");
}
