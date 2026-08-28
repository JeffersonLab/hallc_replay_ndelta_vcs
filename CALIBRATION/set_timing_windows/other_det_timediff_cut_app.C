// other_det_timediff_cut_app.C
//
// Generalization of hodo_timediff_cut_app.C to the OTHER per-PMT/per-block
// AdcTdcDiffTime timing cuts documented in your analysis notes (Tables 4/5)
// and timing_window_setup.C's calc_timing_windows(): Calorimeter and
// Cherenkov. Hodoscope is intentionally EXCLUDED here -- that's
// hodo_timediff_cut_app.C's job; keeping it in both tools would risk two
// independent JSON/progress files disagreeing on the same channels. Same
// architecture throughout -- live TTree::Draw()+Iteration$+1 histogram
// building, mult=1/2/3 overlay, two-tier reference, single JSON output
// with source provenance, batch-PDF with per-channel zoom and leak-safe
// deferred deletion, Zoom In/Out, "Go to #..." input dialog, and the
// never-reviewed/no-value reports.
//
// NOT COVERED: Drift Chamber (per-plane raw TDC hit-lists, not per-PMT-
// slot arrays -- structurally different) and TRIG (a different indexed-
// lookup param scheme entirely). Both would need a separate tool.
//
// --------------------------------------------------------------------------
// DETECTOR GROUPS (pass one of these as the `group` argument, or "all"):
//
//   cal_h     HMS  Calorimeter H.cal.{1pr,2ta,3ta,4ta}.good{Pos,Neg}AdcTdcDiffTime
//   cal_pr_p  SHMS Preshower   P.cal.pr.good{Pos,Neg}AdcTdcDiffTime
//   cal_fly_p SHMS Shower      P.cal.fly.goodAdcTdcDiffTime      (single-ended)
//   cer_h     HMS  Cherenkov   H.cer.goodAdcTdcDiffTime          (single-ended)
//   hgcer_p   SHMS HG Cherenkov P.hgcer.goodAdcTdcDiffTime       (single-ended)
//
// NGCER is intentionally NOT a group -- not used for this experiment
// (replaced by a vacuum chamber).
//
// Branch capitalization ("Good" vs "good") is taken literally from Tables
// 4/5. This is exactly the kind of detail that's easy to get subtly wrong
// without a live tree to check against, so if a group comes up with
// "could not build histogram" for everything, that capitalization is the
// first thing to check.
//
// Usage (interactive, must NOT use -q):
//   root -l 'other_det_timediff_cut_app.C("cal_h", 26107)'
//   root -l 'other_det_timediff_cut_app.C("cer_h", 26107, 26092)'   // + reference run
//   root -l 'other_det_timediff_cut_app.C("all", 26107)'            // every non-hodo group in one session
//
// Batch PDF (3rd arg = true):
//   root -l -b -q 'other_det_timediff_cut_app.C("cal_pr_p", 26107, -1, true)'
//
// Controls: identical to hodo_timediff_cut_app.C -- click twice (lo/hi),
// << Prev / Next >>, Reset clicks, Use Reference, Pause (save),
// Save && Finish, Go to #... (input dialog), Zoom In/Out (starts zoomed
// in, -30 to 70 by default). Keys: n/b/g/z/r/p/s.
//
// --------------------------------------------------------------------------
// OUTPUT: <outDir>/det_diff_cuts_<run>.json, one file covering every group
// touched this session, each channel tagged with its group/spec/detBase/
// plane/side/ipmt plus lo/hi/source (manual, reference-confirmed, or
// reference-auto) and ref_lo/ref_hi when a baseline reference exists.
// Save & Finish also fills any untouched-but-referenced channel from the
// baseline automatically (tagged reference-auto), and reports both real
// gaps (no value at all) and never-reviewed channels (reference-auto) --
// same behavior as hodo_timediff_cut_app.C.
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
#include <TGClient.h>
#include <TGInputDialog.h>
#include <TRootCanvas.h>
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
// Detector group configuration
// ==========================================================================

enum class RefIndexScheme {
  kFlat,          // index = ipmt-1 -- single flat array, no plane concept
  kHodoPmtMajor,  // index = (ipmt-1)*4 + planeIdx, padded to a common max nPMT
  kCalPlaneMajor  // index = planeIdx*blockPerPlane + (ipmt-1) -- HMS Cal's row layout
};

struct DetGroup {
  TString key;                    // selector, e.g. "cal_h"
  TString label;                  // human-readable, for titles/printouts
  TString spec;                   // "p" or "h"
  TString detBase;                // "cal", "cer", "hgcer"
  std::vector<TString> planes;    // e.g. {"1x","1y","2x","2y"}, {"pr"}, {"fly"}, or {} (no plane segment)
  std::vector<TString> sides;     // {"Pos","Neg"} or {} (single-ended, no side segment)
  TString goodWord;               // "Good" or "good" -- literal per Tables 4/5
  // vanilla PARAM/ baseline reference -- key names/array sizes/index
  // schemes verified against real uploaded param files for every group
  // except NGCER (removed, not used for this experiment)
  TString vanillaPath;            // relative to paramDir
  TString vanillaKeyBase;         // e.g. "cal" -> keys "{vanillaKeyBase}_pos_AdcTimeWindowMin" etc.,
                                   // or "" if this group has no side-specific keys (see vanillaKeyFlat)
  TString vanillaKeyFlat;         // for single-ended/no-side groups: bare key prefix, e.g. "cer_adcTimeWindowMin"
  int vanillaLen;                 // expected total flat-array length
  RefIndexScheme scheme;
  int blockPerPlane;              // only used by kCalPlaneMajor (13 for HMS Cal)
};

std::vector<DetGroup> BuildGroupTable() {
  std::vector<DetGroup> g;
  // Hodoscope intentionally excluded -- that's hodo_timediff_cut_app.C's
  // job. Keeping it here too would mean two tools could independently
  // write conflicting JSON/progress files for the same channels.
  g.push_back({ "cal_h", "HMS Calorimeter", "h", "cal", {"1pr","2ta","3ta","4ta"}, {"Pos","Neg"}, "good",
                "HMS/CAL/hcal_cuts.param", "cal", "", 4 * 13, RefIndexScheme::kCalPlaneMajor, 13 });
  g.push_back({ "cal_pr_p", "SHMS Preshower", "p", "cal", {"pr"}, {"Pos","Neg"}, "good",
                "SHMS/CAL/pcal_cuts.param", "cal", "", 14, RefIndexScheme::kFlat, 0 });
  g.push_back({ "cal_fly_p", "SHMS Shower", "p", "cal", {"fly"}, {}, "good",
                "SHMS/CAL/pcal_cuts.param", "", "cal_arr_AdcTimeWindow", 224, RefIndexScheme::kFlat, 0 });
  g.push_back({ "cer_h", "HMS Cherenkov", "h", "cer", {}, {}, "good",
                "HMS/CER/hcer_cuts.param", "", "cer_adcTimeWindow", 2, RefIndexScheme::kFlat, 0 });
  g.push_back({ "hgcer_p", "SHMS HG Cherenkov", "p", "hgcer", {}, {}, "good",
                "SHMS/HGCER/phgcer_cuts.param", "", "hgcer_adcTimeWindow", 4, RefIndexScheme::kFlat, 0 });
  // NGCER intentionally excluded -- not used for this experiment (replaced
  // by a vacuum chamber). Don't add it back without checking first.
  return g;
}

static std::vector<DetGroup> gGroupTable = BuildGroupTable();

const DetGroup* FindGroup(const TString& key) {
  for (auto& g : gGroupTable) if (g.key == key) return &g;
  return nullptr;
}

// ==========================================================================
// Channel definition
// ==========================================================================

struct DetChannel {
  TString groupKey;
  TString spec;
  TString detBase;
  int     planeIdx;   // index into the group's planes list, or -1 if no plane segment
  TString plane;       // "" if no plane segment
  TString side;        // "" if no side segment
  TString goodWord;
  int     ipmt;         // 1-based
};

// {spec}.{detBase}{.plane}.{good}{side}AdcTdcDiffTime -- verified against
// Tables 4/5 for every group above (see BuildGroupTable() call sites).
// {SPEC}.{detBase}{.plane}.{good}{side}AdcTdcDiffTime -- verified against
// Tables 4/5 for every group above (see BuildGroupTable() call sites).
// NOTE: tree branches use uppercase P./H., but c.spec is lowercase "p"/"h"
// (the internal convention also used for param-file prefixes like
// hcal_/pcal_, which genuinely ARE lowercase) -- so the spec letter must be
// upper-cased here specifically for branch construction, same as
// hodo_timediff_cut_app.C's `c.spec == "p" ? "P" : "H"`. Using c.spec
// directly would silently build "p.cal..." instead of "P.cal...", making
// every single histogram fail to build.
TString BranchBase(const DetChannel& c) {
  TString specUp = (c.spec == "p") ? "P" : "H";
  TString s = specUp + "." + c.detBase;
  if (c.plane.Length()) s += "." + c.plane;
  return s;
}
TString BranchName(const DetChannel& c) { return BranchBase(c) + "." + c.goodWord + c.side + "AdcTdcDiffTime"; }
TString MultBranchName(const DetChannel& c) { return BranchBase(c) + "." + c.goodWord + c.side + "AdcMult"; }
TString ChanLabel(const DetChannel& c) { return BranchName(c) + Form("[%d]", c.ipmt); }  // DISPLAY ONLY

// Internal map/storage key. Deliberately does NOT depend on goodWord ("Good"
// vs "good") -- groupKey already fixes spec+detBase+goodWord uniquely, so
// this is fully reconstructable from the JSON fields alone (group/plane/
// side/ipmt), unlike ChanLabel(). Using ChanLabel() as the map key would
// silently break cross-session lookups: LoadJsonAsMap() has no "good" field
// to read back (it was never serialized), so a channel reloaded from JSON
// would reconstruct with goodWord="" and never match the key a live session
// computes via BuildChannelsForGroup() (which does know goodWord).
TString CanonicalKey(const DetChannel& c) { return Form("%s|%s|%s|%d", c.groupKey.Data(), c.plane.Data(), c.side.Data(), c.ipmt); }

// ==========================================================================
// Live histogram building, identical mechanism to hodo_timediff_cut_app.C:
// tree->Draw("<branch>:Iteration$+1>>h2d(...)", "<cut>", "goff")
// ==========================================================================

static TTree* gTree = nullptr;
static int gDiffNbinsY = 400;
static double gDiffYmin = -200, gDiffYmax = 200;
static double gZoomInXMin = -30, gZoomInXMax = 70;
static bool gZoomedIn = true;

double CurViewXMin() { return gZoomedIn ? gZoomInXMin : gDiffYmin; }
double CurViewXMax() { return gZoomedIn ? gZoomInXMax : gDiffYmax; }

static std::map<TString, int> gNPmtCache;  // key = BranchName-ish tag

// Same as hodo's DetectNPmt(): these are Double_t[N] branches (one slot
// per physical PMT/block), N constant across events -- read via
// TLeaf::GetLen() after GetEntry(0).
int DetectNPmt(const DetChannel& proto) {  // proto.ipmt is ignored here
  TString key = BranchName(proto);
  auto it = gNPmtCache.find(key);
  if (it != gNPmtCache.end()) return it->second;
  int n = 0;
  TLeaf* leaf = gTree ? gTree->GetLeaf(key) : nullptr;
  if (leaf && gTree->GetEntries() > 0) { gTree->GetEntry(0); n = leaf->GetLen(); }
  if (n <= 0) printf("WARNING: could not detect element count for branch '%s' -- skipping.\n", key.Data());
  gNPmtCache[key] = n;
  return n;
}

std::vector<DetChannel> BuildChannelsForGroup(const DetGroup& g) {
  std::vector<DetChannel> out;
  std::vector<TString> planes = g.planes.empty() ? std::vector<TString>{ "" } : g.planes;
  std::vector<TString> sides = g.sides.empty() ? std::vector<TString>{ "" } : g.sides;
  for (size_t pi = 0; pi < planes.size(); ++pi) {
    for (auto& side : sides) {
      DetChannel proto{ g.key, g.spec, g.detBase, (int)pi, planes[pi], side, g.goodWord, 1 };
      int n = DetectNPmt(proto);
      for (int ipmt = 1; ipmt <= n; ++ipmt) {
        DetChannel c = proto; c.ipmt = ipmt;
        out.push_back(c);
      }
    }
  }
  return out;
}

TString Hist2DTag(const DetChannel& c) { return "h2d_" + BranchBase(c) + "_" + c.side; }

static std::map<TString, TH2F*> gHist2DCache;
static std::map<TString, TH2F*> gHist2DMultCache;
static const int kMultColors[3] = { kGreen + 2, kMagenta + 1, kOrange + 3 };

TH2F* GetHist2D(const DetChannel& c) {
  TString tag = Hist2DTag(c);
  auto it = gHist2DCache.find(tag);
  if (it != gHist2DCache.end()) return it->second;
  int nPmt = DetectNPmt(c);
  TH2F* h = nullptr;
  if (gTree && nPmt > 0) {
    TString branch = BranchName(c);
    TString expr = Form("%s:Iteration$+1>>%s(%d,0.5,%.1f,%d,%.1f,%.1f)",
                         branch.Data(), tag.Data(), nPmt, nPmt + 0.5, gDiffNbinsY, gDiffYmin, gDiffYmax);
    gTree->Draw(expr, "", "goff");
    h = (TH2F*)gDirectory->Get(tag);
    if (h) { h->SetDirectory(nullptr); h->SetStats(kFALSE); }
  }
  gHist2DCache[tag] = h;
  if (!h) printf("WARNING: could not build 2D hist for %s\n", BranchName(c).Data());
  return h;
}

TH2F* GetHist2DMultCut(const DetChannel& c, int multVal) {
  TString tag = Hist2DTag(c) + Form("_m%d", multVal);
  auto it = gHist2DMultCache.find(tag);
  if (it != gHist2DMultCache.end()) return it->second;
  TString multBranch = MultBranchName(c);
  bool hasMult = gTree && gTree->GetBranch(multBranch) != nullptr;
  int nPmt = DetectNPmt(c);
  TH2F* h = nullptr;
  if (hasMult && nPmt > 0) {
    TString diffBranch = BranchName(c);
    TString expr = Form("%s:Iteration$+1>>%s(%d,0.5,%.1f,%d,%.1f,%.1f)",
                         diffBranch.Data(), tag.Data(), nPmt, nPmt + 0.5, gDiffNbinsY, gDiffYmin, gDiffYmax);
    TString cutExpr = Form("%s==%d", multBranch.Data(), multVal);
    gTree->Draw(expr, cutExpr, "goff");
    h = (TH2F*)gDirectory->Get(tag);
    if (h) { h->SetDirectory(nullptr); h->SetStats(kFALSE); }
  } else if (!hasMult && multVal == 1) {
    printf("NOTE: no '%s' branch -- multiplicity-cut overlays skipped for %s\n",
           multBranch.Data(), BranchName(c).Data());
  }
  gHist2DMultCache[tag] = h;
  return h;
}

TH1D* GetOrBuildProjection(const DetChannel& c) {
  TH2F* h2 = GetHist2D(c);
  if (!h2) return nullptr;
  if (c.ipmt < 1 || c.ipmt > h2->GetNbinsX()) return nullptr;
  TString tag = Form("h_%s_pmt%d_proj", Hist2DTag(c).Data(), c.ipmt);
  TH1D* h1 = (TH1D*)h2->ProjectionY(tag, c.ipmt, c.ipmt);
  h1->SetDirectory(nullptr);
  h1->SetLineColor(kBlue);
  h1->SetFillColorAlpha(kBlue, 0.3);
  h1->SetStats(kFALSE);
  h1->GetXaxis()->SetRangeUser(CurViewXMin(), CurViewXMax());
  return h1;
}

TH1D* GetOrBuildProjectionMultCut(const DetChannel& c, int multVal) {
  TH2F* h2 = GetHist2DMultCut(c, multVal);
  if (!h2) return nullptr;
  if (c.ipmt < 1 || c.ipmt > h2->GetNbinsX()) return nullptr;
  TString tag = Form("h_%s_pmt%d_m%d_proj", Hist2DTag(c).Data(), c.ipmt, multVal);
  TH1D* h1 = (TH1D*)h2->ProjectionY(tag, c.ipmt, c.ipmt);
  h1->SetDirectory(nullptr);
  h1->SetLineColor(kMultColors[multVal - 1]);
  h1->SetLineWidth(2);
  h1->SetStats(kFALSE);
  h1->GetXaxis()->SetRangeUser(CurViewXMin(), CurViewXMax());
  return h1;
}

// ==========================================================================
// Minimal, regex-free text helpers (same approach as hodo_timediff_cut_app.C)
// ==========================================================================

TString ReadFile(const TString& path) {
  std::ifstream in(path.Data());
  if (!in.is_open()) return "";
  std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  return TString(content.c_str());
}
TString NowString() { TDatime now; return now.AsSQLString(); }

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
// Stops at the first line that ISN'T a pure numeric continuation of the
// value list: blank, starts with ';' (comment), or starts with a letter/
// underscore (the next key's assignment). NOT simply "next blank line" --
// several real param files stack multiple single-line keys back-to-back
// with no blank line between them (e.g. pcal_cuts.param's four Preshower
// window keys, phgcer_cuts.param's Min/Max pair at EOF with no trailing
// blank line), and the old blank-line-only version would silently swallow
// the next several keys' values into one oversized block.
TString ExtractBlock(const TString& text, const TString& key) {
  Ssiz_t start = text.Index(key);
  if (start == kNPOS) return "";
  Ssiz_t eqPos = text.Index("=", start);
  if (eqPos == kNPOS) return "";
  Ssiz_t pos = eqPos + 1;
  Ssiz_t n = text.Length();
  Ssiz_t end = n;
  while (pos < n) {
    Ssiz_t nl = text.Index("\n", pos);
    Ssiz_t lineEnd = (nl == kNPOS) ? n : nl;
    TString line = text(pos, lineEnd - pos);
    Ssiz_t i = 0;
    while (i < line.Length() && (line[i] == ' ' || line[i] == '\t' || line[i] == '\r')) ++i;
    if (i >= line.Length() || line[i] == ';' || isalpha((unsigned char)line[i]) || line[i] == '_') {
      end = pos;
      break;
    }
    if (nl == kNPOS) { end = n; break; }
    pos = nl + 1;
  }
  return text(eqPos + 1, end - (eqPos + 1));
}

struct ChannelState { double lo = 0, hi = 0; bool hasLo = false, hasHi = false; TString source = "manual"; };
using DetResultsMap = std::map<TString, ChannelState>;  // key = CanonicalKey()

int RefIndex(const DetGroup& g, const DetChannel& c) {
  switch (g.scheme) {
    case RefIndexScheme::kFlat: return c.ipmt - 1;
    case RefIndexScheme::kHodoPmtMajor: return (c.ipmt - 1) * 4 + c.planeIdx;
    case RefIndexScheme::kCalPlaneMajor: return c.planeIdx * g.blockPerPlane + (c.ipmt - 1);
  }
  return -1;
}

// Loads the vanilla PARAM/ file for one group and converts it into a
// DetResultsMap, restricted to the group's REAL channels. (0,0) windows
// are treated as "disabled/no reference", same convention as hodo.
DetResultsMap LoadVanillaParamAsMap(const DetGroup& g, const TString& paramDir,
                                     const std::vector<DetChannel>& realChannels) {
  DetResultsMap out;
  TString path = Form("%s/%s", paramDir.Data(), g.vanillaPath.Data());
  TString text = ReadFile(path);
  if (text.Length() == 0) { printf("No vanilla defaults found at %s.\n", path.Data()); return out; }

  TString specPrefix = g.spec;  // e.g. "h" or "p" -- param keys are prefixed same as hodo's

  std::vector<double> posMin, posMax, negMin, negMax, flatMin, flatMax;
  bool haveSide = g.vanillaKeyBase.Length() > 0;
  if (haveSide) {
    posMin = SplitNumbers(ExtractBlock(text, specPrefix + g.vanillaKeyBase + "_pos_AdcTimeWindowMin"));
    posMax = SplitNumbers(ExtractBlock(text, specPrefix + g.vanillaKeyBase + "_pos_AdcTimeWindowMax"));
    negMin = SplitNumbers(ExtractBlock(text, specPrefix + g.vanillaKeyBase + "_neg_AdcTimeWindowMin"));
    negMax = SplitNumbers(ExtractBlock(text, specPrefix + g.vanillaKeyBase + "_neg_AdcTimeWindowMax"));
  } else {
    flatMin = SplitNumbers(ExtractBlock(text, specPrefix + g.vanillaKeyFlat + "Min"));
    flatMax = SplitNumbers(ExtractBlock(text, specPrefix + g.vanillaKeyFlat + "Max"));
  }

  int nSkippedOff = 0;
  for (auto& c : realChannels) {
    if (c.groupKey != g.key) continue;
    int idx = RefIndex(g, c);
    ChannelState cs;
    if (haveSide) {
      const std::vector<double>& lo = (c.side == "Pos") ? posMin : negMin;
      const std::vector<double>& hi = (c.side == "Pos") ? posMax : negMax;
      if ((int)lo.size() != g.vanillaLen || (int)hi.size() != g.vanillaLen) continue;
      if (idx < 0 || idx >= g.vanillaLen) continue;
      cs.lo = lo[idx]; cs.hi = hi[idx];
    } else {
      if ((int)flatMin.size() != g.vanillaLen || (int)flatMax.size() != g.vanillaLen) continue;
      if (idx < 0 || idx >= g.vanillaLen) continue;
      cs.lo = flatMin[idx]; cs.hi = flatMax[idx];
    }
    if (cs.lo == 0.0 && cs.hi == 0.0) { ++nSkippedOff; continue; }
    cs.hasLo = cs.hasHi = true;
    out[CanonicalKey(c)] = cs;
  }
  if (nSkippedOff) printf("Skipped %d channel(s) with a (0.00, 0.00) window in %s (treated as disabled).\n",
                           nSkippedOff, path.Data());
  return out;
}

// ---- our own JSON output (write + read-back), same shape as hodo's but
// with group/detBase/plane added per channel ----

TString ToJson(int run, const std::vector<DetChannel>& realChannels,
                const DetResultsMap& merged, const DetResultsMap& baseline) {
  TString out;
  out += "{\n";
  out += Form("  \"run\": %d,\n", run);
  out += Form("  \"generated\": \"%s\",\n", NowString().Data());
  out += "  \"channels\": [\n";
  bool first = true;
  for (auto& c : realChannels) {
    TString key = CanonicalKey(c);
    auto it = merged.find(key);
    if (it == merged.end() || !it->second.hasLo || !it->second.hasHi) continue;
    if (!first) out += ",\n";
    first = false;
    out += Form("    {\"group\": \"%s\", \"spec\": \"%s\", \"det\": \"%s\", \"plane\": \"%s\", "
                "\"side\": \"%s\", \"ipmt\": %d, \"lo\": %.3f, \"hi\": %.3f, \"source\": \"%s\"",
                 c.groupKey.Data(), c.spec.Data(), c.detBase.Data(), c.plane.Data(), c.side.Data(),
                 c.ipmt, it->second.lo, it->second.hi, it->second.source.Data());
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
  while (end < line.Length() && (isdigit((unsigned char)line[end]) || line[end] == '-' || line[end] == '.')) ++end;
  if (end <= start) return false;
  TString numStr = line(start, end - start);
  if (!numStr.IsFloat()) return false;
  val = numStr.Atof();
  return true;
}

DetResultsMap LoadJsonAsMap(const TString& path) {
  DetResultsMap out;
  TString text = ReadFile(path);
  if (text.Length() == 0) return out;
  std::stringstream ss(text.Data());
  std::string lineStd;
  while (std::getline(ss, lineStd)) {
    TString line(lineStd.c_str());
    if (line.Index("\"group\"") == kNPOS) continue;
    TString groupKey = ExtractJsonString(line, "group");
    TString spec = ExtractJsonString(line, "spec");
    TString detBase = ExtractJsonString(line, "det");
    TString plane = ExtractJsonString(line, "plane");
    TString side = ExtractJsonString(line, "side");
    double ipmtD, lo, hi;
    if (groupKey.Length() == 0 || spec.Length() == 0 || detBase.Length() == 0) continue;
    if (!ExtractJsonNumber(line, "ipmt", ipmtD)) continue;
    if (!ExtractJsonNumber(line, "lo", lo)) continue;
    if (!ExtractJsonNumber(line, "hi", hi)) continue;
    DetChannel c{ groupKey, spec, detBase, 0, plane, side, "", (int)ipmtD };
    ChannelState cs; cs.lo = lo; cs.hi = hi; cs.hasLo = cs.hasHi = true;
    TString source = ExtractJsonString(line, "source");
    cs.source = source.Length() > 0 ? source : "manual";
    out[CanonicalKey(c)] = cs;
  }
  return out;
}

// ==========================================================================
// Interactive state
// ==========================================================================

static std::vector<DetChannel> gChannels;
static int gIdx = 0;
static int gRun = 0;
static TString gOutDir, gParamDir, gRootPath;

static DetResultsMap gResults;
static std::vector<double> gClicks;
static bool gFreshVisit = true;
static std::vector<TObject*> gLines, gRefLines;
static TCanvas* gCanvas = nullptr;
static TControlBar* gControlBar = nullptr;
static TH1D* gHcur = nullptr;
static TH1D* gHcurMult[3] = { nullptr, nullptr, nullptr };
static TString gProgressPath;
static TLegend* gCurLegend = nullptr;
static std::vector<TObject*> gLegendProxies;
static std::vector<TObject*> gStaticPagePrimitives;

static DetResultsMap gBaselineRef;
static DetResultsMap gThisRunRef;

TString CurKey() { return CanonicalKey(gChannels[gIdx]); }

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

bool GetWindowFrom(const DetResultsMap& src, const DetChannel& c, double& lo, double& hi) {
  auto it = src.find(CanonicalKey(c));
  if (it == src.end() || !it->second.hasLo || !it->second.hasHi) return false;
  lo = it->second.lo; hi = it->second.hi;
  return true;
}
bool GetBaselineWindow(const DetChannel& c, double& lo, double& hi) { return GetWindowFrom(gBaselineRef, c, lo, hi); }

// ---- drawing ----

void RedrawCutLines() {
  for (auto* o : gLines) delete o;
  gLines.clear();
  if (!gHcur) return;
  double ymax = gHcur->GetMaximum() * 2;
  for (double x : gClicks) {
    TLine* ln = new TLine(x, 0, x, ymax);
    ln->SetLineColor(kOrange + 1); ln->SetLineWidth(2); ln->Draw();
    gLines.push_back(ln);
  }
  gPad->Modified(); gPad->Update();
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
    ln->SetLineColor(kGray + 2); ln->SetLineStyle(2); ln->Draw();
    gRefLines.push_back(ln);
  }
}
void UpdateLegend() {
  if (gCurLegend) { delete gCurLegend; gCurLegend = nullptr; }
  for (auto* o : gLegendProxies) delete o;
  gLegendProxies.clear();

  TLegend* leg = new TLegend(0.785, 0.55, 0.995, 0.92);
  leg->SetTextSize(0.024);
  leg->AddEntry(gHcur, "no cut", "l");
  for (int m = 0; m < 3; ++m) if (gHcurMult[m]) leg->AddEntry(gHcurMult[m], Form("mult=%d", m + 1), "l");
  double refLo, refHi;
  if (GetBaselineWindow(gChannels[gIdx], refLo, refHi)) {
    TLine* proxy = new TLine(0, 0, 0, 0);
    proxy->SetLineColor(kGray + 2); proxy->SetLineStyle(2);
    leg->AddEntry(proxy, Form("cuts-ref (%.1f, %.1f)", refLo, refHi), "l");
    gLegendProxies.push_back(proxy);
  }
  if (!gClicks.empty()) {
    TString loS = Form("%.1f", gClicks[0]);
    TString hiS = gClicks.size() >= 2 ? Form("%.1f", gClicks[1]) : "-";
    TLine* proxy = new TLine(0, 0, 0, 0);
    proxy->SetLineColor(kOrange + 1);
    leg->AddEntry(proxy, Form("cuts-this run %d (%s, %s)", gRun, loS.Data(), hiS.Data()), "l");
    gLegendProxies.push_back(proxy);
  }
  leg->Draw();
  gCurLegend = leg;
  gPad->Modified(); gPad->Update();
}

void StoreCurrent(bool isFreshAction = false, const TString& source = "manual") {
  TString key = CurKey();
  if (gClicks.size() == 2) {
    ChannelState cs;
    cs.lo = std::min(gClicks[0], gClicks[1]); cs.hi = std::max(gClicks[0], gClicks[1]);
    cs.hasLo = cs.hasHi = true;
    if (isFreshAction) cs.source = source;
    else { auto existing = gResults.find(key); cs.source = (existing != gResults.end()) ? existing->second.source : "manual"; }
    gResults[key] = cs;
  } else if (gClicks.empty()) {
    gResults.erase(key);
  }
  SaveProgress();
}

void DrawChannel() {
  const DetChannel& c = gChannels[gIdx];
  gPad->Clear();
  gPad->SetRightMargin(0.24);

  if (gHcur) { delete gHcur; gHcur = nullptr; }
  for (int m = 0; m < 3; ++m) { if (gHcurMult[m]) { delete gHcurMult[m]; gHcurMult[m] = nullptr; } }
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
  for (int m = 0; m < 3; ++m) {
    gHcurMult[m] = GetOrBuildProjectionMultCut(c, m + 1);
    if (gHcurMult[m]) gHcurMult[m]->Draw("HIST SAME");
  }
  gPad->SetLogy();

  auto it = gResults.find(CurKey());
  if (it != gResults.end() && it->second.hasLo) {
    gClicks.clear(); gClicks.push_back(it->second.lo);
    if (it->second.hasHi) gClicks.push_back(it->second.hi);
  } else {
    gClicks.clear();
  }
  gFreshVisit = true;

  DrawReferenceLines();
  RedrawCutLines();
  UpdateLegend();
}

void DrawChannelStatic(TVirtualPad* pad, const DetChannel& c) {
  pad->cd();
  pad->SetLogy();
  pad->SetLeftMargin(0.12); pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.16); pad->SetBottomMargin(0.12);

  TH1D* h = GetOrBuildProjection(c);
  if (!h) {
    TLatex* msg = new TLatex(0.5, 0.5, "missing");
    msg->SetNDC(); msg->SetTextAlign(22); msg->SetTextSize(0.09); msg->SetTextColor(kGray + 2);
    msg->Draw();
    gStaticPagePrimitives.push_back(msg);
    return;
  }
  h->SetTitle("");
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);

  double refLo, refHi;
  bool hasRef = GetBaselineWindow(c, refLo, refHi);
  double lo = 0, hi = 0; bool hasLo = false, hasHi = false;
  auto it = gResults.find(CanonicalKey(c));
  if (it != gResults.end()) { hasLo = it->second.hasLo; lo = it->second.lo;
                               hasHi = it->second.hasHi; hi = it->second.hi; }
  if (hasLo && hasHi) h->GetXaxis()->SetRangeUser(lo - 10, hi + 10);
  else if (hasRef) h->GetXaxis()->SetRangeUser(refLo - 10, refHi + 10);

  h->Draw("HIST");
  TH1D* hMult[3];
  for (int m = 0; m < 3; ++m) {
    hMult[m] = GetOrBuildProjectionMultCut(c, m + 1);
    if (hMult[m]) hMult[m]->Draw("HIST SAME");
  }

  double ymax = h->GetMaximum();
  if (hasRef) {
    for (double x : { refLo, refHi }) {
      TLine* l = new TLine(x, 0, x, ymax * 1.8); l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->Draw();
      gStaticPagePrimitives.push_back(l);
    }
  }
  if (hasLo) { TLine* l = new TLine(lo, 0, lo, ymax * 1.8); l->SetLineColor(kOrange + 1); l->SetLineWidth(2); l->Draw(); gStaticPagePrimitives.push_back(l); }
  if (hasHi) { TLine* l = new TLine(hi, 0, hi, ymax * 1.8); l->SetLineColor(kOrange + 1); l->SetLineWidth(2); l->Draw(); gStaticPagePrimitives.push_back(l); }

  TLatex* title = new TLatex(0.02, 0.90, ChanLabel(c));
  title->SetNDC(); title->SetTextSize(0.075); title->Draw();
  gStaticPagePrimitives.push_back(title);

  gStaticPagePrimitives.push_back(h);
  for (int m = 0; m < 3; ++m) { if (hMult[m]) gStaticPagePrimitives.push_back(hMult[m]); }
}

// ---- interaction ----

void UseReference() {
  double lo, hi;
  if (!GetBaselineWindow(gChannels[gIdx], lo, hi)) { printf("No reference available for %s.\n", CurKey().Data()); return; }
  gClicks = { lo, hi };
  gFreshVisit = false;
  RedrawCutLines();
  UpdateLegend();
  StoreCurrent(true, "reference-confirmed");
  printf("[Use Reference] %s -> lo=%.2f hi=%.2f\n", CurKey().Data(), lo, hi);
}

void GoNext() { StoreCurrent(); if (gIdx < (int)gChannels.size() - 1) { ++gIdx; DrawChannel(); } }
void GoPrev() { StoreCurrent(); if (gIdx > 0) { --gIdx; DrawChannel(); } }
void PauseSave() { StoreCurrent(); printf("[progress saved]\n"); }

void GoToChannelPrompt() {
  StoreCurrent();
  char buf[256] = "";
  const TGWindow* parent = gClient->GetRoot();
  if (gCanvas && gCanvas->GetCanvasImp()) {
    TRootCanvas* rc = dynamic_cast<TRootCanvas*>(gCanvas->GetCanvasImp());
    if (rc) parent = rc;
  }
  new TGInputDialog(gClient->GetRoot(), parent,
                     Form("Enter channel # (1-%d), or a label like P.hod.1x.GoodPosAdcTdcDiffTime[3]:",
                          (int)gChannels.size()), "", buf);
  TString input(buf);
  input = input.Strip(TString::kBoth);
  if (input.Length() == 0) { printf("(cancelled)\n"); return; }
  if (input.IsDigit()) {
    int idx = input.Atoi();
    if (idx >= 1 && idx <= (int)gChannels.size()) { gIdx = idx - 1; DrawChannel(); return; }
    printf("Index %d out of range (1-%d)\n", idx, (int)gChannels.size());
    return;
  }
  for (size_t i = 0; i < gChannels.size(); ++i) {
    if (ChanLabel(gChannels[i]) == input) { gIdx = (int)i; DrawChannel(); return; }
  }
  printf("No channel matches '%s'\n", input.Data());
}

void ToggleZoom() {
  gZoomedIn = !gZoomedIn;
  printf("[view: %s (%.0f, %.0f)]\n", gZoomedIn ? "zoomed in" : "zoomed out", CurViewXMin(), CurViewXMax());
  DrawChannel();
}

DetResultsMap WriteJsonResults() {
  DetResultsMap merged = gThisRunRef;
  for (auto& c : gChannels) {
    auto it = gResults.find(CanonicalKey(c));
    if (it == gResults.end() || !it->second.hasLo || !it->second.hasHi) continue;
    merged[CanonicalKey(c)] = it->second;
  }
  int nFromRef = 0;
  for (auto& c : gChannels) {
    TString key = CanonicalKey(c);
    if (merged.count(key)) continue;
    double lo, hi;
    if (!GetWindowFrom(gBaselineRef, c, lo, hi)) continue;
    ChannelState cs; cs.lo = lo; cs.hi = hi; cs.hasLo = cs.hasHi = true; cs.source = "reference-auto";
    merged[key] = cs;
    ++nFromRef;
  }
  if (nFromRef) printf("Filled %d channel(s) from the baseline reference (no manual cut was ever set).\n", nFromRef);
  if (merged.empty()) { printf("Nothing set -- skipping write.\n"); return merged; }

  TString path = Form("%s/det_diff_cuts_%d.json", gOutDir.Data(), gRun);
  std::ofstream(path.Data()) << ToJson(gRun, gChannels, merged, gBaselineRef);
  printf("Wrote %s (%d channel(s))\n", path.Data(), (int)merged.size());
  return merged;
}

void ReportUnsavedChannels(const DetResultsMap& merged) {
  std::vector<std::pair<int, TString>> noValue, autoFilled;
  for (size_t i = 0; i < gChannels.size(); ++i) {
    const DetChannel& c = gChannels[i];
    auto it = merged.find(CanonicalKey(c));
    if (it == merged.end() || !it->second.hasLo || !it->second.hasHi) noValue.push_back({ (int)i + 1, ChanLabel(c) });
    else if (it->second.source == "reference-auto") autoFilled.push_back({ (int)i + 1, ChanLabel(c) });
  }
  if (noValue.empty() && autoFilled.empty()) { printf("All %d channel(s) have a manually-reviewed cut.\n", (int)gChannels.size()); return; }
  if (!noValue.empty()) {
    printf("\n%d of %d channel(s) have NO value at all (no manual cut, no reference):\n", (int)noValue.size(), (int)gChannels.size());
    for (auto& e : noValue) printf("%d) %s\n", e.first, e.second.Data());
  }
  if (!autoFilled.empty()) {
    printf("\n%d of %d channel(s) were NEVER REVIEWED -- filled from the baseline reference automatically:\n",
           (int)autoFilled.size(), (int)gChannels.size());
    for (auto& e : autoFilled) printf("%d) %s\n", e.first, e.second.Data());
  }
}

void SaveAndFinish() {
  StoreCurrent();
  SaveProgress();
  DetResultsMap merged = WriteJsonResults();
  ReportUnsavedChannels(merged);
  if (gCanvas) gCanvas->Close();
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "Detector DiffTime Cut Controls", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Previous channel");
  gControlBar->AddButton("Next >>", "GoNext();", "Next channel");
  gControlBar->AddButton("Go to #...", "GoToChannelPrompt();", "Jump to a channel by index or label (input dialog)");
  gControlBar->AddButton("Zoom In/Out", "ToggleZoom();", "Toggle between full range and the zoomed-in window");
  gControlBar->AddButton("Reset clicks", "gClicks.clear(); RedrawCutLines(); UpdateLegend();", "Clear this channel's clicks");
  gControlBar->AddButton("Use Reference", "UseReference();", "Copy the baseline reference's window into this cut");
  gControlBar->AddButton("Pause (save)", "PauseSave();", "Save progress without finishing");
  gControlBar->AddButton("Save && Finish", "SaveAndFinish();", "Write det_diff_cuts_<run>.json under outDir and close");
  gControlBar->Show();
}

void OnClick() {
  int event = gPad->GetEvent();
  if (event == kKeyPress) {
    int key = gPad->GetEventX();
    if (key == 'n' || key == 'N') GoNext();
    else if (key == 'b' || key == 'B') GoPrev();
    else if (key == 'g' || key == 'G') GoToChannelPrompt();
    else if (key == 'z' || key == 'Z') ToggleZoom();
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
  if (gClicks.size() == 2) StoreCurrent(true, "manual");
}

// ==========================================================================
// Main entry point
// ==========================================================================

void other_det_timediff_cut_app(
    const char* group = "cal_h",        // one group key above, or "all"
    int run = 26107,
    int referenceRun = -1,
    bool nonInteractive = false,
    const char* rootFile = nullptr,
    const char* outDir = "./other_det_timediff_qa",
    const char* paramDir = "../../PARAM",
    int gridCols = 4,
    int gridRows = 4,
    const char* treeName = "T",
    int diffNbinsY = 400,
    double diffYmin = -200,
    double diffYmax = 200,
    double zoomInXMin = -30,
    double zoomInXMax = 70,
    bool startZoomedIn = true)
{
  if (nonInteractive) gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  gHist2DCache.clear(); gHist2DMultCache.clear(); gNPmtCache.clear();
  gRun = run; gOutDir = outDir; gParamDir = paramDir;
  gDiffNbinsY = diffNbinsY; gDiffYmin = diffYmin; gDiffYmax = diffYmax;
  gZoomInXMin = zoomInXMin; gZoomInXMax = zoomInXMax; gZoomedIn = startZoomedIn;
  gSystem->mkdir(gOutDir, true);

  gRootPath = rootFile ? TString(rootFile)
      : Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
  TFile* f = TFile::Open(gRootPath);
  if (!f || f->IsZombie()) { printf("Could not open %s\n", gRootPath.Data()); return; }
  gTree = (TTree*)f->Get(treeName);
  if (!gTree) { printf("Tree '%s' not found in %s\n", treeName, gRootPath.Data()); return; }

  TString groupStr(group);
  std::vector<TString> groupKeys;
  if (groupStr == "all") { for (auto& g : gGroupTable) groupKeys.push_back(g.key); }
  else if (!FindGroup(groupStr)) { printf("Unknown group '%s'. Valid: cal_h, cal_pr_p, cal_fly_p, cer_h, hgcer_p, all\n", group); return; }
  else groupKeys.push_back(groupStr);

  gChannels.clear();
  for (auto& gk : groupKeys) {
    const DetGroup* g = FindGroup(gk);
    for (auto& c : BuildChannelsForGroup(*g)) gChannels.push_back(c);
  }
  if (gChannels.empty()) { printf("No channels found -- check branch names/tree contents.\n"); return; }

  // ---- baseline_ref: --reference-run's JSON, or vanilla PARAM/ per group ----
  if (referenceRun >= 0) {
    TString rPath = Form("%s/det_diff_cuts_%d.json", outDir, referenceRun);
    gBaselineRef = LoadJsonAsMap(rPath);
    if (gBaselineRef.empty()) printf("WARNING: --reference-run %d requested but %s missing/empty.\n", referenceRun, rPath.Data());
  } else {
    gBaselineRef.clear();
    for (auto& gk : groupKeys) {
      const DetGroup* g = FindGroup(gk);
      DetResultsMap m = LoadVanillaParamAsMap(*g, paramDir, gChannels);
      for (auto& kv : m) gBaselineRef[kv.first] = kv.second;
    }
  }

  TString thisRunPath = Form("%s/det_diff_cuts_%d.json", outDir, run);
  gThisRunRef = LoadJsonAsMap(thisRunPath);

  gProgressPath = Form("%s/%s_%d_progress.txt", outDir, groupStr.Data(), run);
  int nLoaded = LoadProgress();
  if (nLoaded > 0) printf("Resuming saved progress from %s: %d channel(s) set.\n", gProgressPath.Data(), nLoaded);

  int nSeeded = 0;
  for (auto& c : gChannels) {
    TString key = CanonicalKey(c);
    if (gResults.count(key)) continue;
    auto it = gThisRunRef.find(key);
    if (it == gThisRunRef.end() || !it->second.hasLo || !it->second.hasHi) continue;
    gResults[key] = it->second;
    ++nSeeded;
  }
  if (nSeeded) printf("Seeded %d channel(s) from this run's own previously-saved JSON.\n", nSeeded);

  gIdx = 0;
  gCanvas = new TCanvas("c1", Form("detector diff-time cuts, run %d (%s)", run, group), 1000, 650);
  DrawChannel();

  if (nonInteractive) {
    TString pdfPath = Form("%s/%d_%s_summary.pdf", outDir, run, groupStr.Data());
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
      for (auto* p : gStaticPagePrimitives) delete p;
      gStaticPagePrimitives.clear();
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

  printf("\n=== run %d -- %d channel(s) loaded (group=%s) ===\n", run, (int)gChannels.size(), group);
  printf("Click on the histogram: 1st = lo, 2nd = hi (orange).\n");
  printf("Use the 'Detector DiffTime Cut Controls' panel, or keys n/b/g/z/r/p/s.\n");
}