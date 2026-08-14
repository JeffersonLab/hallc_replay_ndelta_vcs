// reftime_cut_app.C
//
// ROOT port of reftime_cut_app.py's workflow: reference loading, click-to-
// cut, progress persistence, usage summary, and (new) real param-file
// generation. Started as a single-channel speed test; pass channel="all"
// to run the full 117-channel tool with Prev/Next navigation, exactly
// mirroring the Python app's channel list and output.
//
// Usage (must run interactively -- NOT with -q, since it waits for clicks):
//   root -l reftime_cut_app.C                                   // pDCREF1 only
//   root -l 'reftime_cut_app.C("hT1", 26107)'                   // one channel
//   root -l 'reftime_cut_app.C("all", 26107, 26092)'            // all 117 channels
//
// Non-interactive (batch), 4th arg = true -- just renders the summary PDF
// from whatever's already on record and exits, no display needed:
//   root -l -b -q 'reftime_cut_app.C("pDCREF3", 26107, 26092, true)'
//   root -l -b -q 'reftime_cut_app.C("all", 26107, 26092, true)'          // 3x4 grid/page (default)
//   root -l -b -q 'reftime_cut_app.C("all", 26107, 26092, true, nullptr,
//       "./reftime_qa", "../../PARAM", "T", 300, 2, 6)'                     // 2x6 grid/page instead
//
// Controls:
//   Click twice on the histogram: 1st click = lo, 2nd click = hi (orange).
//   A 3rd click starts over.
//   A separate "Reftime Cut Controls" button panel (a TControlBar window)
//   opens alongside the canvas, mirroring the Python app's buttons:
//     << Prev / Next >>          -- step through channels
//     Skip to next NEEDED >>     -- jump to the next channel that feeds a
//                                    .param value, skipping scaler/trigger-
//                                    word channels
//     Reset clicks               -- clear the current channel's clicks
//     Use Reference               -- copy the baseline reference's value
//                                    into this channel's cut
//     Pause (save)                -- save progress to disk, keep going
//     Save && Finish               -- write tcoin.param / h_reftime_cut_
//                                    coindaq.param / p_reftime_cut.param
//                                    (channel="all") or a single-channel
//                                    summary (single-channel mode), close
//   The same keyboard shortcuts (n/b/j/r/p/s) still work too, if you prefer
//   them with the canvas focused, but the buttons are the primary interface.
//
// --------------------------------------------------------------------------
// Ported from reftime_cut_app.py:
//   - DEFAULT_TDC_NAMES / ADC_CHANNELS / PARAM_MAP / describe_channel_usage
//   - baseline_ref (axvline 2): --reference-run's saved files, hard error if
//     missing; else the vanilla PARAM/TRIG,HMS/GEN,SHMS/GEN defaults. Parsed
//     ONCE for all channels (not re-parsed per channel) -- same idea as
//     Python's load_reference(). Purely a comparison line, never seeds a cut.
//   - this_run_ref (seeds the cut): this run's own already-saved param
//     files, if present. A whole-run progress file (channel|kind|lo|hi|
//     source per line, not JSON, but the same idea) always wins per channel.
//   - Same click semantics: first click after navigating to a channel always
//     starts a fresh pair; non-negative clamp; legend labels formatted like
//     "cuts-this run <run> (lo=.., hi=..)" / "cuts-ref <label> (lo=.., hi=..)".
//   - resolve_window_arrays / resolve_param_values ("most conservative"
//     value across all channels feeding the same param) and the real
//     generate_tcoin_param / generate_hms_param / generate_shms_param
//     output -- channel="all" writes the actual three files now.
//
// New (not in the Python version): cut quality quantification. For each
// channel, RetainedFraction() computes what fraction of the dominant-
// multiplicity "good hit" peak actually survives the resolved cut window
// (manual > reference > none, same fallback as the write-time resolution).
// A channel below --qaThreshold (default 0.80) is flagged. Runs
// automatically on Save & Finish (channel="all") and at the end of
// non-interactive mode, or on demand via the "Check Cut Quality" button.
// Writes <run>_cut_qa.csv (every channel's fraction, for comparing across
// runs) and <run>_flagged.pdf (just the flagged channels, so you don't have
// to open the full 117-channel PDF to find the handful that need attention).
//
// Single-channel mode (channel != "all") patches the real tcoin.param /
// h_reftime_cut_coindaq.param / p_reftime_cut.param directly instead of
// writing a standalone summary: this channel's window entry and whichever
// single-value param(s) it feeds get updated in place, taking the more
// conservative of the existing file value vs. this channel's new one (so
// reviewing channels one at a time across separate sessions never regresses
// a value another channel already set). If the files don't exist yet at
// all, they're first initialized from the reference for every other
// channel/param, so a lone channel's save doesn't leave everything else
// zeroed out.
// --------------------------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TControlBar.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TDatime.h>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdio>

// ==========================================================================
// Channel lists -- ported directly from reftime_cut_app.py
// ==========================================================================

std::vector<TString> BuildDefaultTdcNames() {
  TString all =
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
    "hDCREF5 hRF pRF hHODO_RF pHODO_RF";
  std::vector<TString> out;
  TObjArray* arr = all.Tokenize(" ");
  for (int i = 0; i < arr->GetEntries(); ++i) out.push_back(((TObjString*)arr->At(i))->GetString());
  delete arr;
  return out;
}

std::vector<TString> BuildAdcChannels() {
  return { "pFADC_TREF_ROC2", "hFADC_TREF_ROC1" };
}

struct ParamRef { TString file; TString name; };  // file: "hms" | "shms" | "tcoin"

std::map<TString, std::vector<ParamRef>> BuildParamMap() {
  std::map<TString, std::vector<ParamRef>> m;
  for (int i = 1; i <= 10; ++i) m[Form("pDCREF%d", i)] = { {"shms", "pdc_tdcrefcut"} };
  for (int i = 1; i <= 5; ++i)  m[Form("hDCREF%d", i)] = { {"hms", "hdc_tdcrefcut"} };
  m["pT1"] = { {"shms", "phodo_tdcrefcut"} };
  m["pT2"] = { {"tcoin", "t_coin_trig_tdcrefcut"} };
  m["hT2"] = { {"hms", "hhodo_tdcrefcut"} };
  m["pFADC_TREF_ROC2"] = { {"shms", "phodo_adcrefcut"}, {"shms", "pngcer_adcrefcut"},
                            {"shms", "phgcer_adcrefcut"}, {"shms", "paero_adcrefcut"},
                            {"shms", "pcal_adcrefcut"}, {"tcoin", "t_coin_trig_adcrefcut"} };
  m["hFADC_TREF_ROC1"] = { {"hms", "hhodo_adcrefcut"}, {"hms", "hcer_adcrefcut"},
                            {"hms", "hcal_adcrefcut"} };
  return m;
}

bool IsTdcChannel(const TString& channel) { return !channel.Contains("FADC"); }

TString DescribeChannelUsage(const TString& channel, bool isTdc,
                              const std::map<TString, std::vector<ParamRef>>& paramMap) {
  std::map<TString, std::vector<TString>> groups;
  if (isTdc) groups["tcoin"] = { "t_coin_TdcTimeWindowMin", "t_coin_TdcTimeWindowMax" };
  auto it = paramMap.find(channel);
  if (it != paramMap.end())
    for (const auto& pr : it->second) groups[pr.file].push_back(pr.name);

  TString out;
  for (const char* f : {"hms", "shms", "tcoin"}) {
    auto g = groups.find(f);
    if (g == groups.end()) continue;
    if (out.Length()) out += ", ";
    out += TString(f) + "-";
    for (size_t i = 0; i < g->second.size(); ++i) {
      if (i) out += ", ";
      out += g->second[i];
    }
  }
  return out;
}

// ==========================================================================
// Minimal, regex-free .param file parsing
// ==========================================================================

TString ReadFile(const TString& path) {
  std::ifstream in(path.Data());
  if (!in.is_open()) return "";
  std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  return TString(content.c_str());
}

std::vector<TString> SplitTokens(const TString& s) {
  std::vector<TString> out;
  TObjArray* arr = s.Tokenize(" \t\r\n");
  for (int i = 0; i < arr->GetEntries(); ++i) out.push_back(((TObjString*)arr->At(i))->GetString());
  delete arr;
  return out;
}

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

TString ExtractQuoted(const TString& text, const TString& key) {
  Ssiz_t start = text.Index(key);
  if (start == kNPOS) return "";
  Ssiz_t q1 = text.Index("\"", start);
  if (q1 == kNPOS) return "";
  Ssiz_t q2 = text.Index("\"", q1 + 1);
  if (q2 == kNPOS) return "";
  return text(q1 + 1, q2 - q1 - 1);
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

std::map<TString, double> ParseScalarParams(const TString& text) {
  std::map<TString, double> result;
  std::stringstream ss(text.Data());
  std::string lineStd;
  while (std::getline(ss, lineStd)) {
    TString line(lineStd.c_str());
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0 || line.BeginsWith(";")) continue;
    Ssiz_t semi = line.Index(";");
    if (semi != kNPOS) line = TString(line(0, semi)).Strip(TString::kBoth);
    Ssiz_t eq = line.Index("=");
    if (eq == kNPOS) continue;
    TString key = TString(line(0, eq)).Strip(TString::kBoth);
    TString val = TString(line(eq + 1, line.Length())).Strip(TString::kBoth);
    if (val.Length() == 0 || !val.IsFloat()) continue;
    result[key] = val.Atof();
  }
  return result;
}

// ==========================================================================
// Reference data -- parsed ONCE for all channels (not re-parsed per
// channel), same idea as Python's load_reference().
// ==========================================================================

struct ReferenceData {
  bool valid = false;
  TString label;
  std::map<TString, std::pair<double, double>> window;  // channel -> (min, max)
  std::map<TString, double> channelLo;                   // NEEDED channel -> lo (sign-flipped)
};

ReferenceData LoadReferenceData(const TString& tcoinPath, const TString& hmsPath,
                                 const TString& shmsPath, const TString& label,
                                 const std::map<TString, std::vector<ParamRef>>& paramMap) {
  ReferenceData rd;
  TString tcoinText = ReadFile(tcoinPath);
  TString hmsText = ReadFile(hmsPath);
  TString shmsText = ReadFile(shmsPath);
  if (tcoinText.Length() == 0 && hmsText.Length() == 0 && shmsText.Length() == 0) return rd;

  std::map<TString, double> scalars = ParseScalarParams(tcoinText);
  auto hmsScalars = ParseScalarParams(hmsText);
  auto shmsScalars = ParseScalarParams(shmsText);
  scalars.insert(hmsScalars.begin(), hmsScalars.end());
  scalars.insert(shmsScalars.begin(), shmsScalars.end());

  TString namesStr = ExtractQuoted(tcoinText, "t_coin_tdcNames");
  std::vector<TString> names = SplitTokens(namesStr);
  std::vector<double> mins = SplitNumbers(ExtractBlock(tcoinText, "t_coin_TdcTimeWindowMin"));
  std::vector<double> maxs = SplitNumbers(ExtractBlock(tcoinText, "t_coin_TdcTimeWindowMax"));
  if (!names.empty() && names.size() == mins.size() && names.size() == maxs.size()) {
    for (size_t i = 0; i < names.size(); ++i) rd.window[names[i]] = { mins[i], maxs[i] };
  }

  for (const auto& kv : paramMap) {
    if (kv.second.empty()) continue;
    const TString& p0 = kv.second[0].name;
    auto it = scalars.find(p0);
    if (it != scalars.end()) rd.channelLo[kv.first] = -it->second;
  }

  rd.valid = !rd.window.empty() || !rd.channelLo.empty();
  rd.label = label;
  return rd;
}

// Same derivation as Python's get_reference_values(): prefer the scalar cut
// for a NEEDED channel, fall back to / supplement with the window array.
// The (0, 100000) window entry is the "no cut" sentinel -- not a real window.
bool GetReferenceValues(const ReferenceData& rd, const TString& channel, bool needed,
                         double& lo, bool& hasLo, double& hi, bool& hasHi) {
  hasLo = hasHi = false;
  if (!rd.valid) return false;
  if (needed) {
    auto it = rd.channelLo.find(channel);
    if (it != rd.channelLo.end()) { lo = it->second; hasLo = true; }
  }
  auto wIt = rd.window.find(channel);
  if (wIt != rd.window.end() && !(wIt->second.first == 0.0 && wIt->second.second == 100000.0)) {
    if (!hasLo) { lo = wIt->second.first; hasLo = true; }
    hi = wIt->second.second;
    hasHi = true;
  }
  return hasLo || hasHi;
}

// ==========================================================================
// Param file generation -- ported from generate_tcoin_param/hms/shms
// ==========================================================================

const char* kTrigNames = "pTRIG1_ROC1 pTRIG4_ROC1 pTRIG1_ROC2 pTRIG4_ROC2";
const char* kAdcNamesFull =
    "hASUM hBSUM hCSUM hDSUM hPSHWR hSHWR hAER hCER hFADC_TREF_ROC1 pAER "
    "pHGCER pNGCER pPSHWR pFADC_TREF_ROC2 pHGCER_MOD pNGCER_MOD pHEL_NEG pHEL_POS pHEL_MPS";
const int kNumAdc = 19;
const double kTdcOffset = 300.0, kAdcTdcOffset = 200.0, kTdcChanPerNs = 0.09766, kEHadOffset = 0.0;

TString FormatArray(const std::vector<double>& values, int perLine = 10) {
  TString out;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i > 0) out += ((i % perLine == 0) ? ",\n\t\t\t\t  " : ", ");
    double v = values[i];
    out += (v == (long)v) ? Form("%.0f", v) : Form("%.2f", v);
  }
  return out;
}

TString NowString() {
  TDatime now;
  return now.AsSQLString();
}

TString GenerateTcoinParam(int run, const std::vector<TString>& tdcNames,
                            const std::vector<double>& tdcMin, const std::vector<double>& tdcMax,
                            double trigTdcRefCut, double trigAdcRefCut) {
  TString namesStr;
  for (size_t i = 0; i < tdcNames.size(); ++i) { if (i) namesStr += " "; namesStr += tdcNames[i]; }

  return Form(
    "; Auto-generated by reftime_cut_app.C for run %d on %s\n\n"
    "t_coin_numAdc = %d\n"
    "t_coin_numTdc = %d\n\n"
    "t_coin_tdcoffset = %.1f\n"
    "t_coin_adc_tdc_offset = %.1f\n\n"
    "t_coin_tdcchanperns = %.5f\n"
    "eHadCoinTime_Offset = %.1f\n\n"
    "t_coin_trigNames=\"%s\"\n\n"
    "; tdc cut is on pTRef2\n"
    "; adc cut on pFADC_ROC2\n"
    "t_coin_trig_tdcrefcut = %.1f\n"
    "t_coin_trig_adcrefcut = %.1f\n\n"
    "t_coin_adcNames = \"%s\"\n\n"
    "t_coin_tdcNames = \"%s\"\n\n"
    "t_coin_TdcTimeWindowMin = %s\n\n"
    "t_coin_TdcTimeWindowMax = %s\n",
    run, NowString().Data(), kNumAdc, (int)tdcNames.size(), kTdcOffset, kAdcTdcOffset,
    kTdcChanPerNs, kEHadOffset, kTrigNames, trigTdcRefCut, trigAdcRefCut, kAdcNamesFull,
    namesStr.Data(), FormatArray(tdcMin).Data(), FormatArray(tdcMax).Data());
}

const char* kReftimeHeader = "; Cut to select the Reference time when multiple hits in reference time\n";

TString GenerateHmsParam(int run, double hdc, double hhodoTdc, double hhodoAdc,
                          double hcerAdc, double hcalAdc) {
  return TString(kReftimeHeader) + Form(
    "\n; Auto-generated by reftime_cut_app.C for run %d on %s\n"
    "; cut variable = hDCREF2\n"
    "hdc_tdcrefcut=%.1f\n"
    "; cut variable = hT2\n"
    "hhodo_tdcrefcut=%.1f\n"
    "; cut variable = hFADC_TREF_ROC1\n"
    "hhodo_adcrefcut=%.1f\n"
    "hcer_adcrefcut=%.1f\n"
    "hcal_adcrefcut=%.1f\n",
    run, NowString().Data(), hdc, hhodoTdc, hhodoAdc, hcerAdc, hcalAdc);
}

TString GenerateShmsParam(int run, double pdc, double phodoTdc, double phodoAdc,
                           double pngcerAdc, double phgcerAdc, double paeroAdc, double pcalAdc) {
  return TString(kReftimeHeader) + Form(
    "\n; Auto-generated by reftime_cut_app.C for run %d on %s\n"
    "; cut variable = pDCREF(min)\n"
    "pdc_tdcrefcut=%.1f\n"
    "; cut variable = pT1\n"
    "phodo_tdcrefcut=%.1f\n"
    "; cut variable = pFADC_TREF_ROC2\n"
    "phodo_adcrefcut=%.1f\n"
    "pngcer_adcrefcut=%.1f\n"
    "phgcer_adcrefcut=%.1f\n"
    "paero_adcrefcut=%.1f\n"
    "pcal_adcrefcut=%.1f\n",
    run, NowString().Data(), pdc, phodoTdc, phodoAdc, pngcerAdc, phgcerAdc, paeroAdc, pcalAdc);
}

// ==========================================================================
// Module-level interactive state
// ==========================================================================

struct ChannelState {
  double lo = 0, hi = 0;
  bool hasLo = false, hasHi = false;
  TString source = "manual";  // "manual" | "reference"
};

struct HistSet { TH1F *raw = nullptr, *m1 = nullptr, *m2 = nullptr, *m3 = nullptr; };

static std::vector<std::pair<TString, TString>> gChannels;  // (name, kind)
static std::map<TString, std::vector<ParamRef>> gParamMap;
static std::set<TString> gNeededSet;
static int gIdx = 0;
static int gRun = 0;
static TString gOutDir, gRootPath, gTreeName;
static TTree* gTree = nullptr;
static int gNbins = 300;

static std::map<TString, HistSet> gHistCache;   // built lazily, kept for the whole session
static std::map<TString, ChannelState> gResults;
static std::vector<double> gClicks;             // live edit buffer for the current channel
static bool gFreshVisit = true;
static std::vector<TObject*> gLines, gRefLines;
static TCanvas* gCanvas = nullptr;
static TControlBar* gControlBar = nullptr;
static TH1F* gHraw = nullptr;

static ReferenceData gBaselineRef;   // axvline(2): --reference-run or vanilla PARAM/
static ReferenceData gThisRunRef;    // axvline(1) seed: this run's own saved files
static TString gProgressPath;
static bool gAllMode = false;
static double gQaThreshold = 0.80;  // flag a channel if less than this fraction of its
                                     // dominant-multiplicity peak survives the cut window

TString CurName() { return gChannels[gIdx].first; }
TString CurKind() { return gChannels[gIdx].second; }
bool CurNeeded() { return gNeededSet.count(CurName()) > 0; }

// ---- whole-run progress persistence: "channel|kind|lo|hi|source" per line ----

void SaveProgress() {
  if (gProgressPath.Length() == 0) return;
  std::ofstream out(gProgressPath.Data());
  for (const auto& kv : gResults) {
    const ChannelState& cs = kv.second;
    out << kv.first << "|" << gChannels[0].second << "|"
        << (cs.hasLo ? Form("%.4f", cs.lo) : "") << "|"
        << (cs.hasHi ? Form("%.4f", cs.hi) : "") << "|" << cs.source << "\n";
  }
  out.close();
}

int LoadProgress() {
  TString text = ReadFile(gProgressPath);
  if (text.Length() == 0) return 0;
  int n = 0;
  std::stringstream ss(text.Data());
  std::string lineStd;
  while (std::getline(ss, lineStd)) {
    TObjArray* parts = TString(lineStd.c_str()).Tokenize("|");
    if (parts->GetEntries() < 5) { delete parts; continue; }
    TString name = ((TObjString*)parts->At(0))->GetString();
    ChannelState cs;
    TString loStr = ((TObjString*)parts->At(2))->GetString();
    TString hiStr = ((TObjString*)parts->At(3))->GetString();
    cs.source = ((TObjString*)parts->At(4))->GetString();
    if (loStr.Length() && loStr.IsFloat()) { cs.lo = loStr.Atof(); cs.hasLo = true; }
    if (hiStr.Length() && hiStr.IsFloat()) { cs.hi = hiStr.Atof(); cs.hasHi = true; }
    if (cs.hasLo || cs.hasHi) { gResults[name] = cs; ++n; }
    delete parts;
  }
  return n;
}

// ---- histogram loading (lazy, cached) ----

HistSet& GetOrBuildHist(const TString& name, const TString& kind) {
  auto it = gHistCache.find(name);
  if (it != gHistCache.end()) return it->second;

  HistSet hs;
  TString rawBranch = kind == "tdc" ? Form("T.coin.%s_tdcTimeRaw", name.Data())
                                     : Form("T.coin.%s_adcPulseTimeRaw", name.Data());
  TString multBranch = kind == "tdc" ? Form("T.coin.%s_tdcMultiplicity", name.Data())
                                      : Form("T.coin.%s_adcMultiplicity", name.Data());
  if (!gTree->GetBranch(rawBranch) || !gTree->GetBranch(multBranch)) {
    gHistCache[name] = hs;  // empty -- "missing branch"
    return gHistCache[name];
  }

  TString tag = name; tag.ReplaceAll(".", "_");
  const int colors[3] = {kGreen + 2, kMagenta + 1, kOrange + 3};

  // (re)builds raw + mult=1/2/3 into hs over [lo, hi], discarding any
  // previous partial build under the same name first
  auto build = [&](double lo, double hi) {
    delete hs.raw; delete hs.m1; delete hs.m2; delete hs.m3;
    hs.raw = hs.m1 = hs.m2 = hs.m3 = nullptr;

    hs.raw = new TH1F("h_" + tag + "_raw", "", gNbins, lo, hi);
    hs.raw->SetLineColor(kBlue);
    hs.raw->SetFillColorAlpha(kBlue, 0.3);
    hs.raw->SetStats(kFALSE);
    gTree->Draw(Form("%s>>h_%s_raw", rawBranch.Data(), tag.Data()), "", "goff");

    TH1F** slots[3] = { &hs.m1, &hs.m2, &hs.m3 };
    for (int m = 1; m <= 3; ++m) {
      TH1F* h = new TH1F(Form("h_%s_m%d", tag.Data(), m), "", gNbins, lo, hi);
      h->SetLineColor(colors[m - 1]);
      h->SetLineWidth(2);
      h->SetStats(kFALSE);
      gTree->Draw(Form("%s>>h_%s_m%d", rawBranch.Data(), tag.Data(), m),
                  Form("%s==%d", multBranch.Data(), m), "goff");
      *slots[m - 1] = h;
    }
  };

  const double kDefaultLo = 0.0, kDefaultHi = 100000.0;  // same "no cut" sentinel used elsewhere

  double xmin = gTree->GetMinimum(rawBranch);
  double xmax = gTree->GetMaximum(rawBranch);
  bool usedDefaultRange = (xmax <= xmin);
  if (usedDefaultRange) { xmin = kDefaultLo; xmax = kDefaultHi; }  // degenerate -> default span
  
// Small margin beyond the observed data extent (same idea as the earlier
  // Python version's 0.02x-range padding): without this, if a channel's real
  // peak sits right at the data's observed max, the histogram's plotted
  // range ends exactly there -- leaving no room to click a "max" cut
  // boundary past the peak at all. This is purely display/binning padding;
  // it doesn't touch the separate (0, 100000) "no cut" sentinel used
  // elsewhere for stored cut values.
  double margin = 0.02 * std::max(xmax - xmin, 1.0);
  xmax += margin;

  build(xmin, xmax);

  // Genuinely empty (zero entries) even though GetMinimum/GetMaximum looked
  // fine -- rebuild over the full (0, 100000) default span instead of
  // whatever narrow/degenerate range came out of an empty selection, so an
  // empty channel is visually unambiguous (flat & empty across the whole
  // plausible range) rather than looking like a rendering glitch.
  if (hs.raw->GetEntries() == 0 && !usedDefaultRange) {
    build(kDefaultLo, kDefaultHi);
  }

  gHistCache[name] = hs;
  return gHistCache[name];
}

// Builds raw + mult=1/2/3 fresh, over [rangeLo, rangeHi] with `nbins` bins --
// NOT cached in gHistCache (that cache is for each channel's one full-range
// display copy; this is a throwaway, zoom-specific set). Used by both
// DrawChannelStatic() and ComputeQaRows() so what's plotted and what's
// measured for QA come from the exact same binning, not just "the same
// parameters" -- e.g. run==referenceRun's QA fraction should visually match
// what the plot shows, not just agree in principle.
HistSet BuildZoomedHistSet(const TString& name, const TString& kind,
                            double rangeLo, double rangeHi, int nbins) {
  HistSet hs;
  TString rawBranch = kind == "tdc" ? Form("T.coin.%s_tdcTimeRaw", name.Data())
                                     : Form("T.coin.%s_adcPulseTimeRaw", name.Data());
  TString multBranch = kind == "tdc" ? Form("T.coin.%s_tdcMultiplicity", name.Data())
                                      : Form("T.coin.%s_adcMultiplicity", name.Data());
  TString tag = name; tag.ReplaceAll(".", "_"); tag += "_zoom";

  hs.raw = new TH1F("h_" + tag + "_raw", "", nbins, rangeLo, rangeHi);
  hs.raw->SetLineColor(kBlue);
  hs.raw->SetFillColorAlpha(kBlue, 0.3);
  hs.raw->SetStats(kFALSE);
  gTree->Draw(Form("%s>>h_%s_raw", rawBranch.Data(), tag.Data()), "", "goff");

  const int colors[3] = {kGreen + 2, kMagenta + 1, kOrange + 3};
  TH1F** slots[3] = { &hs.m1, &hs.m2, &hs.m3 };
  for (int m = 1; m <= 3; ++m) {
    TH1F* h = new TH1F(Form("h_%s_m%d", tag.Data(), m), "", nbins, rangeLo, rangeHi);
    h->SetLineColor(colors[m - 1]);
    h->SetLineWidth(2);
    h->SetStats(kFALSE);
    gTree->Draw(Form("%s>>h_%s_m%d", rawBranch.Data(), tag.Data(), m),
                Form("%s==%d", multBranch.Data(), m), "goff");
    *slots[m - 1] = h;
  }
  return hs;
}

// DrawChannelStatic() registers its zoomed histograms here (they must
// survive until the whole PDF page is printed, since Divide() puts several
// channels on one page together) -- cleared once per page, between pages,
// to avoid unbounded accumulation across a full 117-channel/many-page PDF.
static std::vector<TH1F*> gZoomedHistPool;

void ClearZoomedHistPool() {
  for (auto* h : gZoomedHistPool) delete h;
  gZoomedHistPool.clear();
}

// ---- drawing ----

void RedrawCutLines() {
  for (auto* obj : gLines) delete obj;
  gLines.clear();
  if (!gHraw) return;
  double ymin = 0.5, ymax = gHraw->GetMaximum() * 2;
  for (double x : gClicks) {
    TLine* ln = new TLine(x, ymin, x, ymax);
    ln->SetLineColor(kOrange + 1);
    ln->SetLineWidth(2);
    ln->Draw();
    gLines.push_back(ln);
  }
  gPad->Modified();
  gPad->Update();
}

void UpdateLegend() {
  // sits in the right margin (see SetRightMargin in DrawChannel), outside the
  // histogram frame itself rather than overlapping the plotted data
  TLegend* leg = new TLegend(0.745, 0.15, 0.995, 0.90);
  leg->SetTextSize(0.024);
  leg->AddEntry(gHraw, "raw", "f");
  HistSet& hs = gHistCache[CurName()];
  if (hs.m1) leg->AddEntry(hs.m1, "mult=1", "l");
  if (hs.m2) leg->AddEntry(hs.m2, "mult=2", "l");
  if (hs.m3) leg->AddEntry(hs.m3, "mult=3", "l");

  if (gBaselineRef.valid) {
    double lo, hi; bool hasLo, hasHi;
    GetReferenceValues(gBaselineRef, CurName(), CurNeeded(), lo, hasLo, hi, hasHi);
    TString loS = hasLo ? Form("%.1f", lo) : "-";
    TString hiS = hasHi ? Form("%.1f", hi) : "-";
    TLine* proxy = new TLine(0, 0, 0, 0);
    proxy->SetLineColor(kGray + 2);
    proxy->SetLineStyle(2);
    // #splitline{}{} is ROOT's multi-line text syntax -- plain "\n" does not
    // render as a line break here. Nested to get three lines (label/lo/hi).
    // TLegend gives every entry equal row height regardless of line count,
    // so this entry's text is shrunk to keep all three lines inside its row.
    auto* refEntry = leg->AddEntry(proxy,
        Form("#splitline{cuts-ref %s}{#splitline{lo=%s}{hi=%s}}",
             gBaselineRef.label.Data(), loS.Data(), hiS.Data()), "l");
    refEntry->SetTextSize(0.017);
  }
  if (!gClicks.empty()) {
    TString loS = Form("%.1f", gClicks[0]);
    TString hiS = gClicks.size() >= 2 ? Form("%.1f", gClicks[1])
                 : (CurKind() == "adc" ? "(n/a)" : "-");
    TLine* proxy = new TLine(0, 0, 0, 0);
    proxy->SetLineColor(kOrange + 1);
    auto* cutEntry = leg->AddEntry(proxy,
        Form("#splitline{cuts-this run %d}{#splitline{lo=%s}{hi=%s}}", gRun, loS.Data(), hiS.Data()), "l");
    cutEntry->SetTextSize(0.017);
  }
  leg->Draw();
  gPad->Modified();
  gPad->Update();
}

void DrawReferenceOnce() {
  for (auto* obj : gRefLines) delete obj;
  gRefLines.clear();
  if (!gBaselineRef.valid || !gHraw) return;
  double lo, hi; bool hasLo, hasHi;
  if (!GetReferenceValues(gBaselineRef, CurName(), CurNeeded(), lo, hasLo, hi, hasHi)) return;

  double ymax = gHraw->GetMaximum();
  double yMarker = ymax * 0.6;
  for (auto p : { std::make_pair(hasLo, lo), std::make_pair(hasHi, hi) }) {
    if (!p.first) continue;
    TLine* ln = new TLine(p.second, 0.5, p.second, ymax * 2);
    ln->SetLineColor(kGray + 2);
    ln->SetLineStyle(2);
    ln->Draw();
    gRefLines.push_back(ln);
    TMarker* mk = new TMarker(p.second, yMarker, 23);
    mk->SetMarkerColor(kGray + 2);
    mk->SetMarkerSize(1.3);
    mk->Draw();
    gRefLines.push_back(mk);
  }
}

void StoreCurrent() {
  TString name = CurName();
  bool needsTwo = CurKind() == "tdc";
  if (needsTwo) {
    if (gClicks.size() == 2) {
      ChannelState cs;
      cs.lo = std::min(gClicks[0], gClicks[1]);
      cs.hi = std::max(gClicks[0], gClicks[1]);
      cs.hasLo = cs.hasHi = true;
      cs.source = "manual";
      gResults[name] = cs;
    } else if (gClicks.empty()) {
      gResults.erase(name);
    }
  } else {
    if (gClicks.size() == 1) {
      ChannelState cs;
      cs.lo = gClicks[0];
      cs.hasLo = true;
      gResults[name] = cs;
    } else if (gClicks.empty()) {
      gResults.erase(name);
    }
  }
  SaveProgress();
}

// Compact renderer for one grid cell in the multi-channel PDF: same data as
// DrawChannel() (raw+mult histograms, reference line, recorded cut), but a
// small TLatex label instead of the full sidebar TLegend, since a 3-line
// #splitline entry won't fit in a cell that's 1/12th of the page.
void DrawChannelStatic(TVirtualPad* pad, const TString& name, const TString& kind) {
  pad->cd();
  pad->SetLogy();
  pad->SetRightMargin(0.03);
  pad->SetLeftMargin(0.12);
  pad->SetTopMargin(0.14);
  pad->SetBottomMargin(0.12);

  bool needed = gNeededSet.count(name) > 0;
  HistSet& wideHs = GetOrBuildHist(name, kind);
  if (!wideHs.raw) {
    TLatex* msg = new TLatex(0.5, 0.5, Form("%s -- missing", name.Data()));
    msg->SetNDC();
    msg->SetTextAlign(22);
    msg->SetTextSize(0.09);
    msg->SetTextColor(kGray + 2);
    msg->Draw();
    return;
  }

  double refLo, refHi; bool hasRefLo, hasRefHi;
  bool hasRef = GetReferenceValues(gBaselineRef, name, needed, refLo, hasRefLo, refHi, hasRefHi);

  double lo = 0, hi = 0; bool hasLo = false, hasHi = false;
  auto it = gResults.find(name);
  if (it != gResults.end()) {
    hasLo = it->second.hasLo; lo = it->second.lo;
    hasHi = it->second.hasHi; hi = it->second.hi;
  }

  // Zoom to 0.9x/1.1x around whichever cut values (reference and/or this
  // run's) actually exist for this channel -- and, when there's something to
  // zoom to, REBUILD raw+mult=1/2/3 fresh over just that window at fine
  // (100-bin) resolution via the same BuildZoomedHistSet() ComputeQaRows()
  // uses, rather than just narrowing the view on the wide 300-bin cached
  // histogram. This is what's plotted here AND what QA measures against --
  // same object, same binning, so the panel and the QA number agree.
  std::vector<double> anchors;
  if (hasRefLo) anchors.push_back(refLo);
  if (hasRefHi) anchors.push_back(refHi);
  if (hasLo) anchors.push_back(lo);
  if (hasHi) anchors.push_back(hi);

  HistSet zoomedHs;
  bool usedZoom = false;
  if (!anchors.empty()) {
    double vmin = *std::min_element(anchors.begin(), anchors.end());
    double vmax = *std::max_element(anchors.begin(), anchors.end());
    double rangeLo = 0.9 * vmin, rangeHi = 1.1 * vmax;
    if (rangeHi > rangeLo) {
      zoomedHs = BuildZoomedHistSet(name, kind, rangeLo, rangeHi, 100);
      if (zoomedHs.raw) {
        for (TH1F* h : { zoomedHs.raw, zoomedHs.m1, zoomedHs.m2, zoomedHs.m3 })
          if (h) gZoomedHistPool.push_back(h);  // must outlive this call, until the page prints
        usedZoom = true;
      }
    }
  }
  HistSet& hs = usedZoom ? zoomedHs : wideHs;

  hs.raw->SetTitle("");
  hs.raw->GetXaxis()->SetLabelSize(0.055);
  hs.raw->GetYaxis()->SetLabelSize(0.055);
  hs.raw->Draw("HIST");
  if (hs.m1) hs.m1->Draw("HIST SAME");
  if (hs.m2) hs.m2->Draw("HIST SAME");
  if (hs.m3) hs.m3->Draw("HIST SAME");

  // ymax from whichever histogram is actually being shown, so a taller peak
  // sitting outside the (possibly zoomed) window doesn't inflate the
  // reference/cut line height beyond what's visible here.
  double ymax = hs.raw->GetMaximum();
  hs.raw->SetMaximum(ymax * 3);  // matching headroom for the panel's own y-axis too,
                                  // so there's no leftover empty space from the old range
  gPad->Modified();

  if (hasRef) {
    if (hasRefLo) {
      TLine* l = new TLine(refLo, 0.5, refLo, ymax * 2);
      l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->Draw();
    }
    if (hasRefHi) {
      TLine* l = new TLine(refHi, 0.5, refHi, ymax * 2);
      l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->Draw();
    }
  }
  if (hasLo) { TLine* l = new TLine(lo, 0.5, lo, ymax * 2); l->SetLineColor(kOrange + 1); l->SetLineWidth(2); l->Draw(); }
  if (hasHi) { TLine* l = new TLine(hi, 0.5, hi, ymax * 2); l->SetLineColor(kOrange + 1); l->SetLineWidth(2); l->Draw(); }

  // Anchored well below the pad's top edge (1.0), not right at it: TLatex's
  // default alignment draws text upward from its y-anchor, so a title placed
  // too close to y=1.0 gets its top clipped by the pad boundary once text
  // height is accounted for.
  TLatex* title = new TLatex(0.02, 0.88, name);
  title->SetNDC();
  title->SetTextSize(0.075);
  title->Draw();

  TString label;
  if (hasRef) label += Form("ref: %s/%s  ", hasRefLo ? Form("%.0f", refLo) : "-",
                             hasRefHi ? Form("%.0f", refHi) : "-");
  label += Form("cut: %s/%s", hasLo ? Form("%.0f", lo) : "-", hasHi ? Form("%.0f", hi) : "-");
  TLatex* lab = new TLatex(0.02, 0.02, label);
  lab->SetNDC();
  lab->SetTextSize(0.065);
  lab->Draw();
}

void DrawChannel() {
  TString name = CurName(), kind = CurKind();

  gPad->Clear();
  gPad->SetRightMargin(0.28);  // leave room to the right of the histogram frame for the legend
  HistSet& hs = GetOrBuildHist(name, kind);
  if (!hs.raw) {
    gPad->cd();
    TLatex* msg = new TLatex(0.5, 0.5, Form("%s (%s) -- missing branch(es), no data", name.Data(), kind.Data()));
    msg->SetTextAlign(22);
    msg->Draw();
    gClicks.clear();
    gPad->Modified();
    gPad->Update();
    return;
  }

  gHraw = hs.raw;
  gHraw->Draw("HIST");
  if (hs.m1) hs.m1->Draw("HIST SAME");
  if (hs.m2) hs.m2->Draw("HIST SAME");
  if (hs.m3) hs.m3->Draw("HIST SAME");
  gPad->SetLogy();

  TString usage = DescribeChannelUsage(name, kind == "tdc", gParamMap);
  gHraw->SetTitle(Form("[%d/%d] %s (%s) -- used by %s;Raw TDC/ADC (channel);Counts",
                        gIdx + 1, (int)gChannels.size(), name.Data(), kind.Data(), usage.Data()));

  // restore this channel's stored cut, if any
  auto it = gResults.find(name);
  if (it != gResults.end() && it->second.hasLo) {
    gClicks.clear();
    gClicks.push_back(it->second.lo);
    if (it->second.hasHi) gClicks.push_back(it->second.hi);
  } else {
    gClicks.clear();
  }
  gFreshVisit = true;

  DrawReferenceOnce();
  RedrawCutLines();
  UpdateLegend();
}

// ---- interaction ----

void UseReference() {
  if (!gBaselineRef.valid) { printf("No reference available for %s.\n", CurName().Data()); return; }
  double lo, hi; bool hasLo, hasHi;
  if (!GetReferenceValues(gBaselineRef, CurName(), CurNeeded(), lo, hasLo, hi, hasHi)) return;

  gClicks.clear();
  if (hasLo) gClicks.push_back(lo);
  if (CurKind() == "tdc" && hasHi) gClicks.push_back(hi);
  gFreshVisit = false;
  RedrawCutLines();
  UpdateLegend();

  // tag this result as reference-sourced once committed
  StoreCurrent();
  auto it = gResults.find(CurName());
  if (it != gResults.end()) it->second.source = "reference";
  SaveProgress();

  printf("[Use Reference] %s -> lo=%s hi=%s\n", CurName().Data(),
         hasLo ? Form("%.2f", lo) : "-", hasHi ? Form("%.2f", hi) : "-");
}

void GoNext() { StoreCurrent(); if (gIdx < (int)gChannels.size() - 1) { ++gIdx; DrawChannel(); } }
void GoPrev() { StoreCurrent(); if (gIdx > 0) { --gIdx; DrawChannel(); } }

void GoNextNeeded() {
  StoreCurrent();
  int j = gIdx + 1;
  while (j < (int)gChannels.size() && !gNeededSet.count(gChannels[j].first)) ++j;
  if (j < (int)gChannels.size()) { gIdx = j; DrawChannel(); }
  else printf("No more NEEDED channels after this one.\n");
}

// resolve a channel's lo, same fallback as Python's resolve_channel_lo(): 0.0 if never set
// resolve a channel's lo for writing the final file(s): the manual cut if
// set, else the reference value (so saving without visiting every channel
// still writes real reference-derived values instead of silently zeroing
// them out), else 0.0 as the last resort if no reference exists either
double ResolveChannelLo(const TString& name) {
  auto it = gResults.find(name);
  if (it != gResults.end() && it->second.hasLo) return it->second.lo;

  bool needed = gNeededSet.count(name) > 0;
  double lo, hi; bool hasLo, hasHi;
  if (GetReferenceValues(gBaselineRef, name, needed, lo, hasLo, hi, hasHi) && hasLo) return lo;
  return 0.0;
}

// ==========================================================================
// Single-channel mode: patch the real .param file(s) for just this channel's
// contribution, rather than writing a standalone summary. "Most conservative
// wins" against whatever's already there, so reviewing channels one at a
// time across separate sessions never regresses a value another channel
// already set -- except when the existing value is the "unset" sentinel
// (0 for a scalar param, (0, 100000) for a window entry), which always
// just gets replaced outright rather than compared against.
// ==========================================================================

// Patches one "key = value" / "key=value" line's numeric value with
// min(existing_lo, newLo) (converting from the file's negated-value
// convention), preserving every other line untouched. Appends a new line if
// the key isn't present at all. existing == 0.0 is treated as "not really
// set" and simply replaced with newLo, not compared against.
TString PatchScalarParam(const TString& text, const TString& key, double newLo) {
  std::stringstream ss(text.Data());
  std::string lineStd;
  std::vector<TString> outLines;
  bool found = false;
  while (std::getline(ss, lineStd)) {
    TString line(lineStd.c_str());
    TString stripped = TString(line).Strip(TString::kBoth);
    bool matched = false;
    if (!stripped.BeginsWith(";") && stripped.Length() > 0) {
      Ssiz_t eq = stripped.Index("=");
      if (eq != kNPOS) {
        TString k = TString(stripped(0, eq)).Strip(TString::kBoth);
        if (k == key) {
          TString valPart = TString(stripped(eq + 1, stripped.Length())).Strip(TString::kBoth);
          double existingLo = valPart.IsFloat() ? -valPart.Atof() : newLo;
          double finalLo = (existingLo == 0.0) ? newLo : std::min(existingLo, newLo);
          outLines.push_back(Form("%s=%.1f", key.Data(), -finalLo));
          found = true;
          matched = true;
        }
      }
    }
    if (!matched) outLines.push_back(line);
  }
  if (!found) outLines.push_back(Form("%s=%.1f", key.Data(), -newLo));

  TString result;
  for (auto& l : outLines) { result += l; result += "\n"; }
  return result;
}

// Patches this channel's single entry in the t_coin_TdcTimeWindowMin/Max
// arrays with min(existing, newLo)/max(existing, newHi) -- unless the
// existing entry is still the (0, 100000) "no cut" sentinel, in which case
// it's just replaced outright. Everything else in the file is untouched;
// only the two array blocks get reformatted (values unchanged except this
// one channel's).
TString PatchTcoinWindow(const TString& text, const TString& channel, double newLo, double newHi) {
  TString namesStr = ExtractQuoted(text, "t_coin_tdcNames");
  std::vector<TString> names = SplitTokens(namesStr);
  int idx = -1;
  for (size_t i = 0; i < names.size(); ++i) if (names[i] == channel) { idx = (int)i; break; }
  if (idx < 0) return text;  // channel isn't a TDC channel in this file -- nothing to patch

  TString minBlockOld = ExtractBlock(text, "t_coin_TdcTimeWindowMin");
  TString maxBlockOld = ExtractBlock(text, "t_coin_TdcTimeWindowMax");
  std::vector<double> mins = SplitNumbers(minBlockOld);
  std::vector<double> maxs = SplitNumbers(maxBlockOld);
  if ((int)mins.size() != (int)names.size() || (int)maxs.size() != (int)names.size()) return text;

  bool existingIsSentinel = (mins[idx] == 0.0 && maxs[idx] == 100000.0);
  if (existingIsSentinel) {
    mins[idx] = newLo;
    maxs[idx] = newHi;
  } else {
    mins[idx] = std::min(mins[idx], newLo);
    maxs[idx] = std::max(maxs[idx], newHi);
  }

  TString result = text;
  result.ReplaceAll(minBlockOld, FormatArray(mins));
  result.ReplaceAll(maxBlockOld, FormatArray(maxs));
  return result;
}

// If the run's param files don't exist yet at all, initialize them from the
// reference (or vanilla defaults) for every channel/param -- so a lone
// channel's save doesn't produce a file that's mostly zeros for everything
// it isn't reviewing right now.
void InitializeParamFilesFromReference(const TString& tcoinPath, const TString& hmsPath,
                                        const TString& shmsPath) {
  std::vector<TString> tdcNames = BuildDefaultTdcNames();

  std::vector<double> tdcMin, tdcMax;
  for (const auto& name : tdcNames) {
    bool needed = gNeededSet.count(name) > 0;
    double lo, hi; bool hasLo, hasHi;
    if (GetReferenceValues(gBaselineRef, name, needed, lo, hasLo, hi, hasHi) && hasLo && hasHi) {
      tdcMin.push_back(lo); tdcMax.push_back(hi);
    } else {
      tdcMin.push_back(0.0); tdcMax.push_back(100000.0);
    }
  }

  std::map<TString, double> paramValue;
  for (const auto& kv : gParamMap) {
    double lo, hi; bool hasLo, hasHi;
    double chLo = 0.0;
    if (GetReferenceValues(gBaselineRef, kv.first, true, lo, hasLo, hi, hasHi) && hasLo) chLo = lo;
    for (const auto& pr : kv.second) {
      auto it = paramValue.find(pr.name);
      if (it == paramValue.end() || chLo < it->second) paramValue[pr.name] = chLo;
    }
  }
  auto v = [&](const char* name) -> double {
    auto it = paramValue.find(name);
    return it != paramValue.end() ? -it->second : 0.0;
  };

  std::ofstream(tcoinPath.Data()) << GenerateTcoinParam(
      gRun, tdcNames, tdcMin, tdcMax, v("t_coin_trig_tdcrefcut"), v("t_coin_trig_adcrefcut"));
  std::ofstream(hmsPath.Data()) << GenerateHmsParam(
      gRun, v("hdc_tdcrefcut"), v("hhodo_tdcrefcut"), v("hhodo_adcrefcut"),
      v("hcer_adcrefcut"), v("hcal_adcrefcut"));
  std::ofstream(shmsPath.Data()) << GenerateShmsParam(
      gRun, v("pdc_tdcrefcut"), v("phodo_tdcrefcut"), v("phodo_adcrefcut"),
      v("pngcer_adcrefcut"), v("phgcer_adcrefcut"), v("paero_adcrefcut"), v("pcal_adcrefcut"));

  printf("Initialized fresh param files from %s, since none existed yet for run %d.\n",
         gBaselineRef.valid ? gBaselineRef.label.Data() : "hardcoded defaults (no reference available)",
         gRun);
}

void SaveSingleChannelToParamFiles() {
  TString name = CurName(), kind = CurKind();
  auto it = gResults.find(name);
  if (it == gResults.end() || !it->second.hasLo) {
    printf("\nNo cut set for %s -- nothing to save.\n", name.Data());
    return;
  }
  double newLo = it->second.lo;
  double newHi = it->second.hasHi ? it->second.hi : 100000.0;  // window array only; ADC has no hi

  TString tcoinPath = Form("%s/tcoin_%d.param", gOutDir.Data(), gRun);
  TString hmsPath = Form("%s/h_reftime_cut_coindaq_%d.param", gOutDir.Data(), gRun);
  TString shmsPath = Form("%s/p_reftime_cut_%d.param", gOutDir.Data(), gRun);

  if (gSystem->AccessPathName(tcoinPath) || gSystem->AccessPathName(hmsPath) ||
      gSystem->AccessPathName(shmsPath)) {
    InitializeParamFilesFromReference(tcoinPath, hmsPath, shmsPath);
  }

  TString tcoinText = ReadFile(tcoinPath);
  TString hmsText = ReadFile(hmsPath);
  TString shmsText = ReadFile(shmsPath);

  if (kind == "tdc") tcoinText = PatchTcoinWindow(tcoinText, name, newLo, newHi);

  auto pmIt = gParamMap.find(name);
  if (pmIt != gParamMap.end()) {
    for (const auto& pr : pmIt->second) {
      if (pr.file == "tcoin") tcoinText = PatchScalarParam(tcoinText, pr.name, newLo);
      else if (pr.file == "hms") hmsText = PatchScalarParam(hmsText, pr.name, newLo);
      else if (pr.file == "shms") shmsText = PatchScalarParam(shmsText, pr.name, newLo);
    }
  }

  std::ofstream(tcoinPath.Data()) << tcoinText;
  std::ofstream(hmsPath.Data()) << hmsText;
  std::ofstream(shmsPath.Data()) << shmsText;

  TString usage = DescribeChannelUsage(name, kind == "tdc", gParamMap);
  printf("\nUpdated:\n  %s\n  %s\n  %s\n", tcoinPath.Data(), hmsPath.Data(), shmsPath.Data());
  printf("(%s feeds: %s)\n", name.Data(), usage.Data());
}

void WriteAllParamFiles() {
  std::vector<TString> tdcNames;
  for (auto& ch : gChannels) if (ch.second == "tdc") tdcNames.push_back(ch.first);

  std::vector<double> tdcMin, tdcMax;
  for (const auto& name : tdcNames) {
    auto it = gResults.find(name);
    if (it != gResults.end() && it->second.hasLo && it->second.hasHi) {
      tdcMin.push_back(it->second.lo);
      tdcMax.push_back(it->second.hi);
    } else {
      // fall back to the reference window (if a full lo+hi pair exists there)
      // rather than the hardcoded (0, 100000) "no cut" default -- same reasoning
      // as ResolveChannelLo()
      bool needed = gNeededSet.count(name) > 0;
      double lo, hi; bool hasLo, hasHi;
      if (GetReferenceValues(gBaselineRef, name, needed, lo, hasLo, hi, hasHi) && hasLo && hasHi) {
        tdcMin.push_back(lo);
        tdcMax.push_back(hi);
      } else {
        tdcMin.push_back(0.0);
        tdcMax.push_back(100000.0);
      }
    }
  }

  // resolve_param_values: most conservative (smallest lo) across all channels feeding each param
  std::map<TString, double> paramValue;  // param name -> smallest raw lo seen
  for (const auto& kv : gParamMap) {
    double lo = ResolveChannelLo(kv.first);
    for (const auto& pr : kv.second) {
      auto it = paramValue.find(pr.name);
      if (it == paramValue.end() || lo < it->second) paramValue[pr.name] = lo;
    }
  }
  auto v = [&](const char* name) -> double {
    auto it = paramValue.find(name);
    return it != paramValue.end() ? -it->second : 0.0;
  };

  TString tcoinText = GenerateTcoinParam(gRun, tdcNames, tdcMin, tdcMax,
                                          v("t_coin_trig_tdcrefcut"), v("t_coin_trig_adcrefcut"));
  TString hmsText = GenerateHmsParam(gRun, v("hdc_tdcrefcut"), v("hhodo_tdcrefcut"),
                                      v("hhodo_adcrefcut"), v("hcer_adcrefcut"), v("hcal_adcrefcut"));
  TString shmsText = GenerateShmsParam(gRun, v("pdc_tdcrefcut"), v("phodo_tdcrefcut"),
                                        v("phodo_adcrefcut"), v("pngcer_adcrefcut"),
                                        v("phgcer_adcrefcut"), v("paero_adcrefcut"), v("pcal_adcrefcut"));

  TString tcoinPath = Form("%s/tcoin_%d.param", gOutDir.Data(), gRun);
  TString hmsPath = Form("%s/h_reftime_cut_coindaq_%d.param", gOutDir.Data(), gRun);
  TString shmsPath = Form("%s/p_reftime_cut_%d.param", gOutDir.Data(), gRun);

  std::ofstream(tcoinPath.Data()) << tcoinText;
  std::ofstream(hmsPath.Data()) << hmsText;
  std::ofstream(shmsPath.Data()) << shmsText;

  printf("\nWrote:\n  %s\n  %s\n  %s\n", tcoinPath.Data(), hmsPath.Data(), shmsPath.Data());
}

// ==========================================================================
// Cut quality: for each channel, what fraction of its dominant-multiplicity
// "good hit" peak actually survives the cut window? A cut clipping most of
// the peak (low retained fraction) gets flagged, so instead of eyeballing
// all 117 channels every run, you only need to look at the flagged ones.
// ==========================================================================

TH1F* DominantMultHist(const HistSet& hs) {
  TH1F* best = nullptr;
  double bestIntegral = -1;
  for (TH1F* h : { hs.m1, hs.m2, hs.m3 }) {
    if (!h) continue;
    double integ = h->Integral();
    if (integ > bestIntegral) { bestIntegral = integ; best = h; }
  }
  return best;
}

// Fraction of `h`'s integral that falls within [lo, hi]. An unset bound is
// treated as open-ended (matches the sign-convention note: a TDC cut with no
// hi, or an ADC cut which never has one, means "keep anything past lo").
// Returns -1 if there's nothing to judge against (empty histogram).
double RetainedFraction(TH1F* h, double lo, bool hasLo, double hi, bool hasHi) {
  if (!h) return -1.0;
  double total = h->Integral();
  if (total <= 0) return -1.0;
  int b1 = hasLo ? h->GetXaxis()->FindBin(lo) : 1;
  int b2 = hasHi ? h->GetXaxis()->FindBin(hi) : h->GetNbinsX();
  // strictly interior bins: b1 and b2 themselves may only be partially inside
  // the true [lo, hi] window (the boundary can fall mid-bin), so they're
  // excluded rather than counted as fully in or fully out. Falls back to the
  // inclusive range if the window is too narrow to leave anything after that.
  int loBin = b1 + 1, hiBin = b2 - 1;
  if (loBin > hiBin) { loBin = b1; hiBin = b2; }
  return h->Integral(loBin, hiBin) / total;
}

// Same manual > reference > none fallback as ResolveChannelLo()/
// WriteAllParamFiles(), but returns the full (lo, hi) window with hasLo/hasHi
// flags intact, since RetainedFraction() needs to know which bound is real.
void ResolveChannelWindow(const TString& name, bool needed,
                           double& lo, bool& hasLo, double& hi, bool& hasHi, TString& source) {
  auto it = gResults.find(name);
  if (it != gResults.end() && (it->second.hasLo || it->second.hasHi)) {
    lo = it->second.lo; hasLo = it->second.hasLo;
    hi = it->second.hi; hasHi = it->second.hasHi;
    source = "manual";
    return;
  }
  if (GetReferenceValues(gBaselineRef, name, needed, lo, hasLo, hi, hasHi)) {
    source = "reference";
    return;
  }
  hasLo = hasHi = false;
  source = "none";
}

struct QaRow {
  TString name, kind;
  double lo = 0, hi = 0;
  bool hasLo = false, hasHi = false;
  TString source = "none";   // "manual" | "reference" | "none" -- where (lo,hi) came from
  double frac = -1;          // retained fraction for the resolved (effective/written) cut
  double refFrac = -1;       // retained fraction for the reference cut alone, regardless of
                              // whether this run overrode it -- lets you tell "I broke this"
                              // apart from "the reference itself is already bad here"
  bool flagged = false;
};

std::vector<QaRow> ComputeQaRows() {
  const int kQaNbins = 100;
  std::vector<QaRow> rows;
  for (auto& ch : gChannels) {
    QaRow row;
    row.name = ch.first;
    row.kind = ch.second;
    bool needed = gNeededSet.count(row.name) > 0;
    ResolveChannelWindow(row.name, needed, row.lo, row.hasLo, row.hi, row.hasHi, row.source);

    double refLo, refHi; bool hasRefLo, hasRefHi;
    GetReferenceValues(gBaselineRef, row.name, needed, refLo, hasRefLo, refHi, hasRefHi);

    // Same union-of-anchors zoom range, and the same BuildZoomedHistSet()
    // call, DrawChannelStatic() uses for the plotted view -- frac and
    // refFrac are measured against exactly what's on the page, not a
    // separately-constructed (if similarly-parameterized) histogram.
    std::vector<double> anchors;
    if (row.hasLo) anchors.push_back(row.lo);
    if (row.hasHi) anchors.push_back(row.hi);
    if (hasRefLo) anchors.push_back(refLo);
    if (hasRefHi) anchors.push_back(refHi);
    if (anchors.empty()) { rows.push_back(row); continue; }  // nothing to judge against

    double vmin = *std::min_element(anchors.begin(), anchors.end());
    double vmax = *std::max_element(anchors.begin(), anchors.end());
    double rangeLo = 0.9 * vmin, rangeHi = 1.1 * vmax;
    if (rangeHi <= rangeLo) { rows.push_back(row); continue; }

    HistSet qhs = BuildZoomedHistSet(row.name, row.kind, rangeLo, rangeHi, kQaNbins);
    if (!qhs.raw) { rows.push_back(row); continue; }

    // dominant multiplicity WITHIN the zoomed window, not the wide range --
    // more physically meaningful here since the wide range's dominance can
    // be driven by background far from the region that actually matters
    TH1F* dom = DominantMultHist(qhs);
    if (dom) {
      row.frac = RetainedFraction(dom, row.lo, row.hasLo, row.hi, row.hasHi);
      row.refFrac = RetainedFraction(dom, refLo, hasRefLo, refHi, hasRefHi);
    }
    for (TH1F* h : { qhs.raw, qhs.m1, qhs.m2, qhs.m3 }) delete h;  // transient -- not pooled

    row.flagged = (row.frac >= 0 && row.frac < gQaThreshold);
    rows.push_back(row);
  }
  return rows;
}

void WriteQaSummary(bool writeFlaggedPdf) {
  std::vector<QaRow> rows = ComputeQaRows();

  TString csvPath = Form("%s/%d_cut_qa.csv", gOutDir.Data(), gRun);
  std::ofstream csv(csvPath.Data());
  csv << "channel,kind,lo,hi,source,retained_fraction,ref_retained_fraction,flagged\n";
  for (const auto& r : rows) {
    csv << r.name << "," << r.kind << ","
        << (r.hasLo ? Form("%.1f", r.lo) : "") << ","
        << (r.hasHi ? Form("%.1f", r.hi) : "") << ","
        << r.source << ","
        << (r.frac >= 0 ? Form("%.3f", r.frac) : "") << ","
        << (r.refFrac >= 0 ? Form("%.3f", r.refFrac) : "") << ","
        << (r.flagged ? "1" : "0") << "\n";
  }
  csv.close();

  std::vector<QaRow> flagged;
  for (const auto& r : rows) if (r.flagged) flagged.push_back(r);
  std::sort(flagged.begin(), flagged.end(),
            [](const QaRow& a, const QaRow& b) { return a.frac < b.frac; });

  printf("\n=== Cut quality (run %d): %d/%d channel(s) below %.0f%% retained ===\n",
         gRun, (int)flagged.size(), (int)rows.size(), gQaThreshold * 100);
  for (const auto& r : flagged) {
    // when source == "manual", refFrac (if available) shows what the
    // reference alone would have given, so you can tell whether your click
    // made things better or worse than just leaving the reference in place
    TString note = r.source == "manual" && r.refFrac >= 0
        ? Form(" [manual; reference alone gives %.1f%%]", r.refFrac * 100)
        : Form(" [%s]", r.source.Data());
    printf("  %-20s %5.1f%% retained  (lo=%s, hi=%s)%s\n", r.name.Data(), r.frac * 100,
           r.hasLo ? Form("%.1f", r.lo) : "-", r.hasHi ? Form("%.1f", r.hi) : "-", note.Data());
  }
  printf("Full table: %s\n", csvPath.Data());

  if (writeFlaggedPdf && !flagged.empty()) {
    TString pdfPath = Form("%s/%d_flagged.pdf", gOutDir.Data(), gRun);
    int cols = std::min((int)flagged.size(), 3);
    int rowsPerPage = 4;
    int perPage = cols * rowsPerPage;
    int nPages = ((int)flagged.size() + perPage - 1) / perPage;

    TCanvas* qc = new TCanvas("qc", "flagged", cols * 500, rowsPerPage * 380);
    qc->Print(pdfPath + "[");
    for (int page = 0; page < nPages; ++page) {
      qc->Clear();
      qc->Divide(cols, rowsPerPage, 0.001, 0.001);
      for (int cell = 0; cell < perPage; ++cell) {
        int idx = page * perPage + cell;
        if (idx >= (int)flagged.size()) break;
        TVirtualPad* sub = qc->cd(cell + 1);
        DrawChannelStatic(sub, flagged[idx].name, flagged[idx].kind);
      }
      qc->Print(pdfPath);
      ClearZoomedHistPool();  // safe now that this page is fully printed
    }
    qc->Print(pdfPath + "]");
    delete qc;
    printf("Flagged-only PDF (%d channel(s)): %s\n", (int)flagged.size(), pdfPath.Data());
  }
}

void PauseSave() {
  StoreCurrent();
  printf("[progress saved]\n");
}

void SaveAndFinish() {
  StoreCurrent();
  SaveProgress();

  if (gAllMode) {
    WriteAllParamFiles();
    WriteQaSummary(true);
  } else {
    SaveSingleChannelToParamFiles();
  }

  gCanvas->Close();
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "Reftime Cut Controls", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Go to the previous channel");
  gControlBar->AddButton("Next >>", "GoNext();", "Go to the next channel");
  gControlBar->AddButton("Skip to next NEEDED >>", "GoNextNeeded();",
                          "Jump to the next channel that feeds a .param value");
  gControlBar->AddButton("Reset clicks", "gClicks.clear(); RedrawCutLines(); UpdateLegend();",
                          "Clear the current channel's clicks");
  gControlBar->AddButton("Use Reference", "UseReference();",
                          "Copy the reference value into this channel's cut");
  gControlBar->AddButton("Check Cut Quality", "WriteQaSummary(true);",
                          "Flag channels whose cut clips most of the peak, without saving/finishing");
  gControlBar->AddButton("Pause (save)", "PauseSave();",
                          "Save progress to disk without finishing");
  gControlBar->AddButton("Save && Finish", "SaveAndFinish();",
                          "Write the final file(s) and close");
  gControlBar->Show();
}

void OnClick() {
  int event = gPad->GetEvent();

  if (event == kKeyPress) {
    int key = gPad->GetEventX();
    if (key == 'n' || key == 'N') GoNext();
    else if (key == 'b' || key == 'B') GoPrev();
    else if (key == 'j' || key == 'J') GoNextNeeded();
    else if (key == 'r' || key == 'R') UseReference();
    else if (key == 'p' || key == 'P') PauseSave();
    else if (key == 's' || key == 'S') SaveAndFinish();
    return;
  }

  if (event != kButton1Down) return;
  if (!gHraw) return;

  double x = gPad->AbsPixeltoX(gPad->GetEventX());
  x = std::max(0.0, x);

  bool needsTwo = CurKind() == "tdc";
  if (gFreshVisit) {
    gClicks.clear();
    gClicks.push_back(x);
    gFreshVisit = false;
  } else if (needsTwo) {
    if (gClicks.size() >= 2) gClicks.clear();
    gClicks.push_back(x);
  } else {
    gClicks.clear();
    gClicks.push_back(x);
  }

  RedrawCutLines();
  UpdateLegend();

  if (gClicks.size() == 1) printf("lo = %.2f\n", gClicks[0]);
  else printf("lo = %.2f   hi = %.2f\n", std::min(gClicks[0], gClicks[1]), std::max(gClicks[0], gClicks[1]));

  StoreCurrent();
}

// ==========================================================================
// Main entry point
// ==========================================================================

void reftime_cut_app(
    const char* channel = "pDCREF1",   // pass "all" to run every channel, like the Python app
    int run = 26107,
    int referenceRun = -1,
    bool nonInteractive = false,
    const char* rootFile = nullptr,
    const char* outDir = "./reftime_qa",
    const char* paramDir = "../../PARAM",
    const char* treeName = "T",
    int nbins = 300,
    int gridCols = 3,    // channel="all" + nonInteractive: channels per PDF page = gridCols*gridRows
    int gridRows = 4,    // default 3x4 = 12/page; try 2,6 for a 2x6 layout instead
    double qaThreshold = 0.80)  // flag a channel if less than this fraction of its dominant-
                                 // multiplicity peak survives the cut window (0-1)
{
  if (nonInteractive) gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);  // no per-histogram stats box -- the legend covers what's needed
  gQaThreshold = qaThreshold;

  // gHistCache/gZoomedHistPool are plain static globals that otherwise persist
  // across repeated calls within the same ROOT session -- without clearing
  // them here, editing GetOrBuildHist()/BuildZoomedHistSet() and re-running
  // this macro without a full ROOT restart would silently keep using the
  // OLD cached histograms, since GetOrBuildHist() returns early on a cache
  // hit with no way to know the code that built them has changed.
  gHistCache.clear();
  ClearZoomedHistPool();

  gRun = run;
  gOutDir = outDir;
  gTreeName = treeName;
  gNbins = nbins;
  gParamMap = BuildParamMap();
  for (const auto& kv : gParamMap) gNeededSet.insert(kv.first);
  gSystem->mkdir(gOutDir, true);

  gAllMode = (TString(channel) == "all");
  gChannels.clear();
  if (gAllMode) {
    for (auto& n : BuildDefaultTdcNames()) gChannels.push_back({n, "tdc"});
    for (auto& n : BuildAdcChannels()) gChannels.push_back({n, "adc"});
  } else {
    gChannels.push_back({ TString(channel), IsTdcChannel(channel) ? "tdc" : "adc" });
  }

  gRootPath = rootFile ? TString(rootFile)
      : Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/"
             "coin_replay_production_%d_2000000_0.root", run);
  TFile* f = TFile::Open(gRootPath);
  if (!f || f->IsZombie()) { printf("Could not open %s\n", gRootPath.Data()); return; }
  gTree = (TTree*)f->Get(gTreeName);
  if (!gTree) { printf("Tree '%s' not found in %s\n", gTreeName.Data(), gRootPath.Data()); return; }

  // ---- baseline_ref (axvline 2) ----
  if (referenceRun >= 0) {
    TString rTcoin = Form("%s/tcoin_%d.param", outDir, referenceRun);
    TString rHms = Form("%s/h_reftime_cut_coindaq_%d.param", outDir, referenceRun);
    TString rShms = Form("%s/p_reftime_cut_%d.param", outDir, referenceRun);
    if (gSystem->AccessPathName(rTcoin) || gSystem->AccessPathName(rHms) || gSystem->AccessPathName(rShms)) {
      printf("ERROR: --reference-run %d requested but missing param files in %s\n", referenceRun, outDir);
      return;
    }
    gBaselineRef = LoadReferenceData(rTcoin, rHms, rShms, Form("run %d", referenceRun), gParamMap);
  } else {
    TString vTcoin = Form("%s/TRIG/tcoin.param", paramDir);
    TString vHms = Form("%s/HMS/GEN/h_reftime_cut_coindaq.param", paramDir);
    TString vShms = Form("%s/SHMS/GEN/p_reftime_cut.param", paramDir);
    if (!gSystem->AccessPathName(vTcoin) && !gSystem->AccessPathName(vHms) && !gSystem->AccessPathName(vShms)) {
      gBaselineRef = LoadReferenceData(vTcoin, vHms, vShms, "vanilla defaults", gParamMap);
    } else {
      printf("No --reference-run given, and vanilla defaults not found under %s.\n", paramDir);
    }
  }

  // ---- this_run_ref: this run's own already-saved param files ----
  TString sTcoin = Form("%s/tcoin_%d.param", outDir, run);
  TString sHms = Form("%s/h_reftime_cut_coindaq_%d.param", outDir, run);
  TString sShms = Form("%s/p_reftime_cut_%d.param", outDir, run);
  if (!gSystem->AccessPathName(sTcoin) && !gSystem->AccessPathName(sHms) && !gSystem->AccessPathName(sShms)) {
    gThisRunRef = LoadReferenceData(sTcoin, sHms, sShms, Form("this run %d", run), gParamMap);
  }

  gProgressPath = Form("%s/%d_progress.txt", outDir, run);
  int nLoaded = LoadProgress();
  if (nLoaded > 0) {
    printf("Resuming saved progress from %s: %d channel(s) set.\n", gProgressPath.Data(), nLoaded);
  }
  if (gThisRunRef.valid) {
    int nSeeded = 0;
    for (auto& ch : gChannels) {
      if (gResults.count(ch.first)) continue;
      bool needed = gNeededSet.count(ch.first) > 0;
      double lo, hi; bool hasLo, hasHi;
      if (!GetReferenceValues(gThisRunRef, ch.first, needed, lo, hasLo, hi, hasHi)) continue;
      ChannelState cs; cs.lo = lo; cs.hasLo = hasLo; cs.hi = hi; cs.hasHi = hasHi; cs.source = "reference";
      gResults[ch.first] = cs;
      ++nSeeded;
    }
    if (nSeeded) printf("Seeded %d channel(s) from this run's own previously-saved param files.\n", nSeeded);
  }

  gIdx = 0;
  gCanvas = new TCanvas("c1", Form("run %d", run), 1000, 650);
  DrawChannel();

  if (nonInteractive) {
    gSystem->mkdir(outDir, true);
    TString pdfPath = gAllMode ? Form("%s/%d_all_summary.pdf", outDir, run)
                                : Form("%s/%d_%s_summary.pdf", outDir, run, channel);
    if (gAllMode) {
      int perPage = gridCols * gridRows;
      int nPages = ((int)gChannels.size() + perPage - 1) / perPage;
      gCanvas->SetCanvasSize(gridCols * 500, gridRows * 380);  // room per cell to stay legible

      gCanvas->Print(pdfPath + "[");
      for (int page = 0; page < nPages; ++page) {
        gCanvas->Clear();
        gCanvas->Divide(gridCols, gridRows, 0.001, 0.001);
        for (int cell = 0; cell < perPage; ++cell) {
          int idx = page * perPage + cell;
          if (idx >= (int)gChannels.size()) break;
          TVirtualPad* sub = gCanvas->cd(cell + 1);
          DrawChannelStatic(sub, gChannels[idx].first, gChannels[idx].second);
        }
        // No page-level suptitle: Divide() fills the full canvas height with
        // no reserved margin above row 1, so a suptitle here collided with
        // (and got overlapped by) the top row's own per-cell titles. The run
        // number/page number are cheap to get from the PDF filename/viewer
        // page indicator instead of fighting for space here.
        gCanvas->Print(pdfPath);
        ClearZoomedHistPool();  // safe now that this page is fully printed
      }
      gCanvas->Print(pdfPath + "]");
      printf("\nWrote %s (%d channels, %dx%d grid, %d pages -- inspection only, no .param files)\n",
             pdfPath.Data(), (int)gChannels.size(), gridCols, gridRows, nPages);
      WriteQaSummary(true);
    } else {
      gCanvas->Print(pdfPath);
      printf("\nWrote %s (inspection only -- no .param files written)\n", pdfPath.Data());
    }
    delete gCanvas;
    return;
  }

  gCanvas->AddExec("dynamic", "OnClick()");
  gCanvas->Update();
  MakeControlBar();

  printf("\n=== run %d -- %d channel(s) loaded ===\n", run, (int)gChannels.size());
  printf("Click on the histogram: 1st = lo, 2nd = hi (orange). 3rd click resets.\n");
  printf("Use the 'Reftime Cut Controls' button panel for Prev/Next/Use Reference/Pause/Save.\n");
  printf("(Keys still work too, if you prefer: n=next b=back j=next NEEDED r=Use Reference p=Pause s=Save&Finish)\n");
}
