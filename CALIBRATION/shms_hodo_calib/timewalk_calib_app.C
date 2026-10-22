// timewalk_calib_app.C
//
// Interactive click-to-fit tool for the SHMS Hodoscope time-walk
// correction (Eq. 2 of hodo_calib.pdf):
//
//   f_TW(a) = c1 + 1 / (a / TDC_Thrs)^c2
//
// fitted to TDC-ADC Time vs. ADC Pulse Amplitude, per PMT (plane, side,
// paddle). Replaces the old two-script batch flow (timeWalkHistos.C then
// timeWalkCalib.C) with one interactive session: the raw histograms are
// built once from the coincidence replay tree (same event-loop logic as
// shms_hodo_calib/timeWalkHistos.C), then you step through every PMT,
// click twice on the amplitude axis to set the fit range (1st click = lo,
// 2nd click = hi), see the fit update live, and move on.
//
// The reference/trigger apparatus branches (FADC_TREF, T1/T2/T3) live
// under "T.coin.*" in this experiment's coincidence replay tree, NOT
// "T.shms.*" -- this app always reads them with the "coin" tree prefix.
//
// A handful of PMTs are permanently off for this experiment (see
// gPermanentlyOff -- currently SHMS 2y Pos {1,2,5,19,20,21} and 2y Neg
// {1,2,19,20,21}) and are excluded from the channel list entirely, so
// no histogram is booked and no fit is attempted for them; the staged
// .param file writes (0,0) for those channels, same convention hcana
// uses for an uncalibrated/disabled paddle. Note 2y's paddle count was
// bumped 18->21 to match: paddles 19-21 physically exist but are always
// off.
//
// Usage (must run interactively -- NOT with -q, it waits for clicks):
//   root -l 'timewalk_calib_app.C(26107)'                        // no reference overlay
//   root -l 'timewalk_calib_app.C(26107, 26092)'                 // + reference run 26092
//   root -l 'timewalk_calib_app.C(26483, 0, "26483-26488")'      // + reference = a chained-run tag
//
// Non-interactive (batch), 4th arg = true -- fits every channel with its
// default/seeded range, writes the summary PDF and param/json files, exits:
//   root -l -b -q 'timewalk_calib_app.C(26107, 26092, "", true)'
//
// Compare-only batch mode: fits every channel fresh (red curve) against a
// loaded reference (black dotted curve) in the summary PDF, but NEVER
// writes JSON/.param output -- use this to sanity-check a new fit against
// the currently-stored calibration without touching it. Needs a
// referenceRun/referenceTag to have anything to compare against:
//   root -l -b -q 'timewalk_calib_app.C(26107, 26092, "", true, "", 1.0e9, "", "./timewalk_qa", "../../PARAM", 1200.0, 20.0, 300.0, true)'
//
// Combining several runs into one TChain (e.g. adding statistics across a
// set of clean calibration runs), with an optional P.gtr.dp cut to keep
// away from elastic-peak bias -- pass a comma-separated run list as the
// "runs" argument (5th positional arg, right after nonInteractive) and it
// overrides the single "run" arg for building histograms; "run" is still
// used as the output tag base when "runs" is empty. Everything after
// dpCut (rootFile, outDir, paramDir, tdcThresh, fitRangeLow/High) stays
// at its default unless you need to override it:
//   root -l 'timewalk_calib_app.C(0, 0, "", false, "26483,26484,26485,26486,26487,26488", 10.0)'
//
// Non-interactive batch fit of the same combined runs, straight to PDF:
//   root -l -b -q 'timewalk_calib_app.C(0, 0, "", true, "26483,26484,26485,26486,26487,26488", 10.0)'
//
// Verifying a chained batch fit interactively, run by run -- referenceTag
// points "Use Reference" / the gray overlay curve at the chain's saved
// result instead of a single reference run:
//   root -l 'timewalk_calib_app.C(26483, 0, "26483-26488")'
//   root -l 'timewalk_calib_app.C(26484, 0, "26483-26488")'
//   ... etc for each run in the chain
//
// Controls:
//   Click twice on the histogram: 1st click = fit range lo, 2nd = hi.
//     The fit re-runs automatically after the 2nd click.
//   << Prev / Next >>   -- step through (plane, side, paddle) channels
//   Refit                -- re-run the fit with the current range
//   Reset Range           -- back to the default [twFitRangeLow, twFitRangeHigh]
//   Use Reference          -- copy the reference run's range+result for
//                            this channel into the current fit, then refit
//   Go to channel...        -- console prompt, index or label "1x.pos[7]"
//   Pause (save)            -- write progress to outDir, keep going
//   Save && Finish            -- write JSON + phodo_TWcalib_<run>.param
//                            under paramDir, close
//
// Outputs (under outDir, default "./timewalk_qa") -- filenames use a "tag":
// a single run number normally, or "<first>-<last>" when "runs" chains
// several runs together (the JSON's "runs" array lists every run actually
// used):
//   timewalk_shms_<tag>.json          -- full per-channel results + provenance
//   timewalk_shms_<tag>_summary.pdf   -- one page per channel, fit overlaid
// Outputs (under paramDir, default "../../PARAM"):
//   SHMS/HODO/phodo_TWcalib_<tag>.param -- staged with the tag in the
//   filename, same convention as reftime_cut_app.C's tcoin_<run>.param;
//   rename to phodo_TWcalib.param for hcana to read it.
//
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TControlBar.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

// ===========================================================================
// Section 1: SHMS constants (from shms_hodo_calib/timeWalkHistos.C)
// ===========================================================================

static const TString gPlaneNames[4] = {"1x", "1y", "2x", "2y"};
static const TString gSideNames[2]  = {"pos", "neg"};
static const Int_t   gNbars[4]      = {13, 13, 14, 21}; // 2y extended 18->21: paddles 19-21 exist but are permanently off (see gPermanentlyOff)
static const Int_t   nBarsMax       = 21; // fixed column padding, matches legacy .param layout

// PMTs permanently off for this experiment -- excluded from the channel
// list entirely (no histogram booked, no fit attempted; WriteLegacyParam
// then naturally writes (0,0) for these, same convention as a paddle
// that was never calibrated).
struct OffPmt { Int_t plane; Int_t side; Int_t paddle; };
static const std::vector<OffPmt> gPermanentlyOff = {
  {3, 0, 1}, {3, 0, 2}, {3, 0, 5}, {3, 0, 19}, {3, 0, 20}, {3, 0, 21}, // SHMS 2y Pos
  {3, 1, 1}, {3, 1, 2},            {3, 1, 19}, {3, 1, 20}, {3, 1, 21}, // SHMS 2y Neg
};

Bool_t IsPermanentlyOff(Int_t plane, Int_t side, Int_t paddle) {
  for (auto &o : gPermanentlyOff)
    if (o.plane == plane && o.side == side && o.paddle == paddle) return kTRUE;
  return kFALSE;
}

static const TString kRefAdcName = "pFADC_TREF_ROC2";
static const std::vector<TString> kRefTdcNames = {"pT1", "pT2", "pT3"};
static const TString kPidBranch = "P.cal.etracknorm"; static const Double_t kPidCutLow = 0.7;
static const TString kCerBranch = "P.hgcer.npeSum";   static const Double_t kCerCutLow = 0.5;
static const Double_t kHodoAmpLo = 15.0, kHodoAmpHi = 1000.0;
static const Double_t kRefAmpLo = 40.0, kRefAmpHi = 70.0;
static const Double_t kRefTimeLo = 300.0, kRefTimeHi = 370.0;

static const Double_t tdcChanToTime = 0.09766; // ns
static const Double_t adcChanToTime = 0.0625;  // ns

// ===========================================================================
// Section 2: channel bookkeeping
// ===========================================================================

struct Channel { Int_t plane; Int_t side; Int_t paddle; };

struct FitRec {
  Double_t lo = 20.0, hi = 300.0;
  Double_t c1 = 1.0, c2 = 1.0;
  Double_t c1err = 0.0, c2err = 0.0;
  Double_t chi2ndf = 0.0;
  Bool_t   fitted = kFALSE;
  TString  source = "default"; // default | manual | reference-confirmed | reference-file | resumed
};

TString ChanLabel(const Channel &c) {
  return Form("%s.%s[%d]", gPlaneNames[c.plane].Data(), gSideNames[c.side].Data(), c.paddle);
}

// ===========================================================================
// Section 3: globals
// ===========================================================================

Int_t    gRun = 0, gReferenceRun = 0;
TString  gReferenceTag; // if non-empty, overrides gReferenceRun for loading a reference (e.g. a chained-run tag)
TString  gRootFile, gOutDir = "./timewalk_qa", gParamDir = "../../PARAM";
Double_t gTdcThresh = 1200.0; // TW fit threshold, units of FADC channels (see timeWalkCalib.C)
Double_t gTwFitRangeLow = 20.0, gTwFitRangeHigh = 300.0;

std::vector<Int_t> gRunList;  // non-empty -> BuildHistos chains all these runs together
TString  gRunTag;             // output-file tag: single run number, or "first-last" when chained
Double_t gDpCut = 1.0e9;      // |P.gtr.dp| < gDpCut cut; default effectively off (elastic-peak avoidance)
Bool_t   gCompareOnly = kFALSE; // batch mode: write the summary PDF (new fit vs. reference) but never touch JSON/.param output

std::vector<Channel> gChannels;
std::map<TString, TH2F*> gHist;        // key = ChanLabel -> raw h2(amp, tdc-adc diff)
std::map<TString, FitRec> gResults;    // key = ChanLabel -> current fit
std::map<TString, FitRec> gRefResults; // key = ChanLabel -> reference run's fit (for overlay/seed)

Int_t     gIdx = 0;
TCanvas  *gCanvas = nullptr;
TControlBar *gControlBar = nullptr;
TF1      *gCurrentFit = nullptr;
TF1      *gRefFit = nullptr;
Int_t     gClickStage = 0; // 0 = waiting for lo click, 1 = waiting for hi click
Double_t  gPendingLo = 0.0;

// Time-walk fit function: c1 + 1/(a/thresh)^c2
Double_t TwFitFunc(Double_t *x, Double_t *p) {
  return p[0] + 1.0 / TMath::Power(x[0] / gTdcThresh, p[1]);
}

// Default ROOT-file path for a single run number, same convention used
// when rootFile isn't explicitly supplied.
TString ResolvedRunPath(Int_t run) {
  return Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
}

// ===========================================================================
// Section 4: build the channel list
// ===========================================================================

void BuildChannelList() {
  gChannels.clear();
  Int_t nSkipped = 0;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t isd = 0; isd < 2; isd++)
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        if (IsPermanentlyOff(ipl, isd, ipad)) { ++nSkipped; continue; }
        gChannels.push_back(Channel{ipl, isd, ipad});
      }
  if (nSkipped) printf("Excluded %d permanently-off PMT(s) from the channel list.\n", nSkipped);
}

// ===========================================================================
// Section 5: build raw histograms (single pass over the coin replay tree)
//   -- same event-selection logic as shms_hodo_calib/timeWalkHistos.C
// ===========================================================================

void BuildHistos() {
  TFile *f = nullptr;
  TChain *chain = nullptr;
  TTree *T = nullptr;

  if (!gRunList.empty()) {
    chain = new TChain("T");
    for (auto r : gRunList) {
      TString path = ResolvedRunPath(r);
      Int_t nAdded = chain->Add(path);
      if (nAdded == 0) printf("WARNING: could not add %s to chain\n", path.Data());
      else printf("  chained: %s\n", path.Data());
    }
    if (chain->GetNtrees() == 0) { printf("ERROR: no files added to chain for runs %s\n", gRunTag.Data()); return; }
    T = chain;
  } else {
    f = new TFile(gRootFile, "READ");
    if (!f || f->IsZombie()) { printf("ERROR: could not open %s\n", gRootFile.Data()); return; }
    T = (TTree*)f->Get("T");
    if (!T) { printf("ERROR: no tree 'T' in %s\n", gRootFile.Data()); return; }
  }

  // Book histograms
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    gHist[key] = new TH2F("h2_" + key,
        Form("TDC-ADC Time vs. Pulse Amp %s %s Paddle %d;Pulse Amplitude (mV);TDC-ADC Time (ns)",
             gPlaneNames[ch.plane].Data(), gSideNames[ch.side].Data(), ch.paddle),
        500, 0, 500, 1500, -50, 50);
    gHist[key]->SetDirectory(nullptr);
  }

  // Reference apparatus branches -- ALWAYS under T.coin.*
  Double_t refAdcTime = 0, refAdcAmp = 0, refAdcMult = 0;
  std::vector<Double_t> refTdc(kRefTdcNames.size(), 0.0);
  T->SetBranchAddress(Form("T.coin.%s_adcPulseTimeRaw", kRefAdcName.Data()), &refAdcTime);
  T->SetBranchAddress(Form("T.coin.%s_adcPulseAmp", kRefAdcName.Data()),     &refAdcAmp);
  T->SetBranchAddress(Form("T.coin.%s_adcMultiplicity", kRefAdcName.Data()), &refAdcMult);
  for (size_t i = 0; i < kRefTdcNames.size(); i++)
    T->SetBranchAddress(Form("T.coin.%s_tdcTimeRaw", kRefTdcNames[i].Data()), &refTdc[i]);

  // PID branches
  Double_t pidVal = 0, cerVal = 0;
  T->SetBranchAddress(kPidBranch, &pidVal);
  T->SetBranchAddress(kCerBranch, &cerVal);

  // Optional momentum acceptance cut (e.g. |P.gtr.dp| < 10 to stay off the
  // elastic peak when combining runs). Off by default (gDpCut ~ 1e9).
  Double_t dpVal = 0;
  T->SetBranchAddress("P.gtr.dp", &dpVal);

  // Hodoscope raw arrays: fixed C-style arrays with parallel Ndata
  // counters, exactly as in timeWalkHistos.C -- one set per (plane, side).
  static const Int_t kMaxHits = 128;
  struct HodoArr {
    Int_t nAdc = 0, nTdc = 0;
    Double_t adcPad[kMaxHits], adcErr[kMaxHits], adcTimeRaw[kMaxHits], adcAmp[kMaxHits];
    Double_t tdcPad[kMaxHits], tdcTimeRaw[kMaxHits];
  };
  std::map<TString, HodoArr*> hodo; // key = "plane.side"

  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t isd = 0; isd < 2; isd++) {
      TString key = Form("%s.%s", gPlaneNames[ipl].Data(), gSideNames[isd].Data());
      HodoArr *h = new HodoArr();
      hodo[key] = h;
      TString base = "P.hod." + gPlaneNames[ipl] + "." + gSideNames[isd];
      T->SetBranchAddress("Ndata." + base + "AdcCounter", &h->nAdc);
      T->SetBranchAddress(base + "AdcCounter",            h->adcPad);
      T->SetBranchAddress(base + "AdcErrorFlag",           h->adcErr);
      T->SetBranchAddress(base + "AdcPulseTimeRaw",        h->adcTimeRaw);
      T->SetBranchAddress(base + "AdcPulseAmp",            h->adcAmp);
      T->SetBranchAddress("Ndata." + base + "TdcCounter", &h->nTdc);
      T->SetBranchAddress(base + "TdcCounter",             h->tdcPad);
      T->SetBranchAddress(base + "TdcTimeRaw",              h->tdcTimeRaw);
    }

  Long64_t nentries = T->GetEntries();
  printf("\n=== Building SHMS time-walk histograms: %lld events ===\n", nentries);
  for (Long64_t ievt = 0; ievt < nentries; ievt++) {
    T->GetEntry(ievt);
    if (ievt % 200000 == 0 && ievt != 0) printf("  %lld / %lld events...\n", ievt, nentries);

    // PID / reference apparatus cuts, same logic as timeWalkHistos.C
    if (pidVal < kPidCutLow) continue;
    if (cerVal < kCerCutLow) continue;
    if (TMath::Abs(dpVal) >= gDpCut) continue; // e.g. |dp| < 10 to avoid elastic-peak bias
    Bool_t refMultCut = (refAdcMult < 1.0);
    Bool_t refAmpCut  = (refAdcAmp < kRefAmpLo || refAdcAmp > kRefAmpHi);
    Double_t refTimeNs = refAdcTime * adcChanToTime;
    Bool_t refTimeCut = (refTimeNs < kRefTimeLo || refTimeNs > kRefTimeHi);
    if (refMultCut || refAmpCut || refTimeCut) continue;

    for (Int_t ipl = 0; ipl < 4; ipl++) {
      for (Int_t isd = 0; isd < 2; isd++) {
        TString hkey = Form("%s.%s", gPlaneNames[ipl].Data(), gSideNames[isd].Data());
        HodoArr *h = hodo[hkey];
        // Plane 2y, pos side historically uses T2 instead of T1 as the
        // TDC reference (see shms_hodo_calib/timeWalkHistos.C).
        Double_t refTdcNs = (ipl == 3 && isd == 0) ? refTdc[1] * tdcChanToTime : refTdc[0] * tdcChanToTime;

        for (Int_t ia = 0; ia < h->nAdc && ia < kMaxHits; ia++) {
          if (h->adcErr[ia] != 0) continue;
          Double_t amp = h->adcAmp[ia];
          if (amp < kHodoAmpLo || amp > kHodoAmpHi) continue;
          Double_t adcPadNum = h->adcPad[ia];
          Double_t adcTimeRawNs = h->adcTimeRaw[ia] * adcChanToTime;
          Double_t adcPulseTime = adcTimeRawNs - refTimeNs;

          for (Int_t it = 0; it < h->nTdc && it < kMaxHits; it++) {
            if ((Int_t)h->tdcPad[it] != (Int_t)adcPadNum) continue;
            Int_t padNum = (Int_t)h->tdcPad[it];
            if (padNum < 1 || padNum > gNbars[ipl]) continue;
            Double_t tdcTimeNs = h->tdcTimeRaw[it] * tdcChanToTime - refTdcNs;
            Double_t diff = tdcTimeNs - adcPulseTime;
            TString ckey = ChanLabel(Channel{ipl, isd, padNum});
            if (gHist.count(ckey)) gHist[ckey]->Fill(amp, diff);
          }
        }
      }
    }
  }
  printf("=== Done: %lld events processed ===\n\n", nentries);
  if (f) f->Close();
  if (chain) delete chain; // histograms already filled and detached (SetDirectory(nullptr))
}

// ===========================================================================
// Section 6: fitting
// ===========================================================================

void DoFit(const Channel &ch, Double_t lo, Double_t hi, const TString &source) {
  TString key = ChanLabel(ch);
  TH2F *h2 = gHist[key];
  FitRec &r = gResults[key];
  r.lo = lo; r.hi = hi; r.source = source;
  if (!h2 || h2->GetEntries() == 0) {
    printf("WARNING: no entries for %s -- cannot fit\n", key.Data());
    r.fitted = kFALSE;
    return;
  }
  TF1 fit("twFit", TwFitFunc, lo, hi, 2);
  fit.SetParameters(1.0, 1.0);
  TFitResultPtr res = h2->Fit(&fit, "SREQ0"); // 0 = don't draw here, caller draws
  Int_t status = res;
  r.fitted = (status == 0);
  r.c1 = fit.GetParameter(0); r.c1err = fit.GetParError(0);
  r.c2 = fit.GetParameter(1); r.c2err = fit.GetParError(1);
  r.chi2ndf = (fit.GetNDF() > 0) ? fit.GetChisquare() / fit.GetNDF() : 0.0;
  if (status != 0) printf("WARNING: fit did not converge cleanly for %s (status=%d)\n", key.Data(), status);
}

// ===========================================================================
// Section 7: legacy hcana .param I/O (for reference-run loading + final save)
// ===========================================================================

// Parses a block like:
//   pc1_Pos = v1x, v1y, v2x, v2y,
//              v1x, v1y, v2x, v2y,
//              ...
// into out[plane][paddleIdx] for paddleIdx = 0..nBarsMax-1.
void ParseParamBlock(std::ifstream &in, const TString &key, Double_t out[4][nBarsMax]) {
  std::string line;
  std::streampos startPos = in.tellg();
  bool found = false;
  while (std::getline(in, line)) {
    if (line.find(key.Data()) != std::string::npos && line.find('=') != std::string::npos) { found = true; break; }
  }
  if (!found) { in.clear(); in.seekg(startPos); return; }
  std::string rest = line.substr(line.find('=') + 1);
  std::vector<Double_t> vals;
  auto pushNums = [&](const std::string &s) {
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      try { vals.push_back(std::stod(tok)); } catch (...) {}
    }
  };
  pushNums(rest);
  while ((Int_t)vals.size() < 4 * nBarsMax) {
    std::streampos p = in.tellg();
    if (!std::getline(in, line)) break;
    TString t(line.c_str()); t = t.Strip(TString::kBoth);
    if (t.Length() == 0 || t.BeginsWith(";") || t.Contains("=")) { in.seekg(p); break; }
    pushNums(line);
  }
  for (Int_t i = 0; i < (Int_t)vals.size() && i < 4 * nBarsMax; i++)
    out[i % 4][i / 4] = vals[i];
}

// Loads a legacy phodo_TWcalib.param (vanilla or staged-with-runnumber)
// into gRefResults. Same convention as hodo_timediff_cut_app.C's cuts
// loader: a (0.00, 0.00) pair for c1/c2 means the paddle was never
// calibrated (turned off / disabled), not a real degenerate fit -- skip
// those so they don't show up as a fake reference curve or get copied
// onto a channel via "Use Reference".
void LoadLegacyParam(const TString &path) {
  std::ifstream in(path.Data());
  if (!in.is_open()) { printf("  (no legacy param file at %s)\n", path.Data()); return; }
  Double_t c1Pos[4][nBarsMax] = {{0}}, c1Neg[4][nBarsMax] = {{0}};
  Double_t c2Pos[4][nBarsMax] = {{0}}, c2Neg[4][nBarsMax] = {{0}};
  ParseParamBlock(in, "pc1_Pos", c1Pos);
  in.clear(); in.seekg(0); ParseParamBlock(in, "pc1_Neg", c1Neg);
  in.clear(); in.seekg(0); ParseParamBlock(in, "pc2_Pos", c2Pos);
  in.clear(); in.seekg(0); ParseParamBlock(in, "pc2_Neg", c2Neg);
  Int_t nSkippedOff = 0, nSkippedPermanent = 0;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      for (Int_t isd = 0; isd < 2; isd++) {
        if (IsPermanentlyOff(ipl, isd, ipad)) { ++nSkippedPermanent; continue; }
        Channel ch{ipl, isd, ipad};
        Double_t c1 = (isd == 0) ? c1Pos[ipl][ipad - 1] : c1Neg[ipl][ipad - 1];
        Double_t c2 = (isd == 0) ? c2Pos[ipl][ipad - 1] : c2Neg[ipl][ipad - 1];
        if (c1 == 0.0 && c2 == 0.0) { ++nSkippedOff; continue; }
        FitRec r;
        r.c1 = c1; r.c2 = c2;
        r.lo = gTwFitRangeLow; r.hi = gTwFitRangeHigh;
        r.fitted = kTRUE; r.source = "reference-file";
        gRefResults[ChanLabel(ch)] = r;
      }
    }
  printf("  loaded reference param: %s", path.Data());
  if (nSkippedOff || nSkippedPermanent)
    printf(" (skipped %d disabled-in-file + %d permanently-off channel(s))", nSkippedOff, nSkippedPermanent);
  printf("\n");
}

TString LegacyParamPath(const TString &tag, Bool_t staged) {
  if (staged) return Form("%s/SHMS/HODO/phodo_TWcalib_%s.param", gParamDir.Data(), tag.Data());
  return Form("%s/SHMS/HODO/phodo_TWcalib.param", gParamDir.Data());
}

void LoadReference() {
  // A referenceTag (e.g. "26483-26488", a chained-run tag) takes priority
  // over a plain referenceRun number when both are given. The literal
  // string "vanilla" explicitly loads PARAM/SHMS/HODO/phodo_TWcalib.param
  // (the untagged file hcana actually reads) rather than relying on a
  // staged-tag lookup incidentally not being found.
  if (gReferenceTag == "vanilla") { LoadLegacyParam(LegacyParamPath("", kFALSE)); return; }
  if (gReferenceTag.Length() == 0 && gReferenceRun == 0) return;
  TString refTag = (gReferenceTag.Length() > 0) ? gReferenceTag : Form("%d", gReferenceRun);
  TString staged = LegacyParamPath(refTag, kTRUE);
  std::ifstream test(staged.Data());
  if (test.is_open()) { test.close(); LoadLegacyParam(staged); }
  else LoadLegacyParam(LegacyParamPath("", kFALSE)); // fall back to vanilla PARAM/ default
}

// Seed gResults for any channel not yet touched this session from this
// run tag's own previously-staged param file (a resumed session).
void SeedFromOwnPriorSave() {
  TString path = LegacyParamPath(gRunTag, kTRUE);
  std::ifstream test(path.Data());
  if (!test.is_open()) return;
  test.close();
  std::map<TString, FitRec> saved = gRefResults; // don't clobber reference map
  gRefResults.clear();
  LoadLegacyParam(path);
  for (auto &kv : gRefResults)
    if (!gResults.count(kv.first)) { gResults[kv.first] = kv.second; gResults[kv.first].source = "resumed"; }
  gRefResults = saved;
}

void WriteLegacyParam(const TString &tag) {
  TString outPath = LegacyParamPath(tag, kTRUE);
  gSystem->mkdir(gSystem->DirName(outPath), kTRUE);
  std::ofstream out(outPath.Data());
  TString runsNote = gRunList.empty() ? tag : TString("runs ") + tag + Form(" (%d chained files)", (int)gRunList.size());
  out << Form(";SHMS Hodoscopes Time Walk Output Parameter File: %s", runsNote.Data()) << std::endl << std::endl;
  out << "pTDC_threshold = " << gTdcThresh << " ;units of FADC channels" << std::endl << std::endl;

  auto writeBlock = [&](const TString &name, Bool_t isC1, Int_t side) {
    out << ";Param " << name << std::endl;
    out << ";           1x             1y             2x             2y " << std::endl;
    out << "p" << name << " = ";
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) {
      for (Int_t ipl = 0; ipl < 4; ipl++) {
        Double_t v = 0.0;
        if (ipad < gNbars[ipl]) {
          Channel ch{ipl, side, ipad + 1};
          TString key = ChanLabel(ch);
          // Only write a real value for channels that actually got a
          // converged fit; anything untouched/unfitted (disabled paddle,
          // empty histogram, or simply never visited) writes (0, 0) --
          // the same "turned off" convention LoadLegacyParam skips on
          // read-back, rather than leaking FitRec's 1.0/1.0 default.
          if (gResults.count(key) && gResults[key].fitted) v = isC1 ? gResults[key].c1 : gResults[key].c2;
        }
        out << std::fixed << v;
        if (ipl != 3) out << ", ";
      }
      out << "," << std::endl;
    }
    out << std::endl;
  };
  writeBlock("c1_Pos", kTRUE, 0);
  writeBlock("c1_Neg", kTRUE, 1);
  writeBlock("c2_Pos", kFALSE, 0);
  writeBlock("c2_Neg", kFALSE, 1);
  out.close();
  printf("Wrote %s\n", outPath.Data());
}

// ===========================================================================
// Section 8: JSON I/O (full state, with provenance)
// ===========================================================================

void SaveJSON(const TString &tag) {
  gSystem->mkdir(gOutDir, kTRUE);
  TString path = Form("%s/timewalk_shms_%s.json", gOutDir.Data(), tag.Data());
  std::ofstream out(path.Data());
  out << "{\n  \"run_tag\": \"" << tag << "\",\n  \"runs\": [";
  for (size_t i = 0; i < gRunList.size(); i++) out << (i ? ", " : "") << gRunList[i];
  if (gRunList.empty()) out << gRun;
  out << "],\n  \"dp_cut\": " << gDpCut << ",\n  \"reference_run\": " << gReferenceRun
      << ",\n  \"tdc_threshold\": " << gTdcThresh << ",\n  \"channels\": [\n";
  Bool_t first = kTRUE;
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    if (!gResults.count(key)) continue;
    FitRec &r = gResults[key];
    if (!first) out << ",\n";
    first = kFALSE;
    out << "    {\"channel\": \"" << key << "\", \"plane\": \"" << gPlaneNames[ch.plane]
        << "\", \"side\": \"" << gSideNames[ch.side] << "\", \"paddle\": " << ch.paddle
        << ", \"lo\": " << r.lo << ", \"hi\": " << r.hi
        << ", \"c1\": " << r.c1 << ", \"c1err\": " << r.c1err
        << ", \"c2\": " << r.c2 << ", \"c2err\": " << r.c2err
        << ", \"chi2ndf\": " << r.chi2ndf << ", \"fitted\": " << (r.fitted ? "true" : "false")
        << ", \"source\": \"" << r.source << "\"}";
  }
  out << "\n  ]\n}\n";
  out.close();
  printf("Wrote %s\n", path.Data());
}

// ===========================================================================
// Section 9: drawing
// ===========================================================================

void DrawChannel() {
  if (gChannels.empty()) return;
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  TH2F *h2 = gHist[key];
  gCanvas->cd();
  gCanvas->Clear();
  if (!h2) { printf("No histogram for %s\n", key.Data()); return; }
  h2->SetStats(kFALSE);
  h2->Draw("COLZ");

  FitRec &r = gResults.count(key) ? gResults[key] : (gResults[key] = FitRec{gTwFitRangeLow, gTwFitRangeHigh});

  // Curves are drawn across the full x-axis (extrapolated beyond the fit
  // range) so you can see where the function is heading outside where it
  // was actually constrained -- the fit itself still only used [lo,hi].
  // Start just above 0 since the function diverges as amplitude -> 0.
  // Drawing "SAME" keeps the pad's existing y-range (set by the histogram),
  // so the divergent part near x=0 is simply clipped by the frame, not
  // rescaled into it.
  Double_t dispLo = TMath::Max(1.0, h2->GetXaxis()->GetXmin());
  Double_t dispHi = h2->GetXaxis()->GetXmax();

  if (gRefResults.count(key)) {
    FitRec &rr = gRefResults[key];
    if (gRefFit) delete gRefFit;
    gRefFit = new TF1("refFit", TwFitFunc, dispLo, dispHi, 2);
    gRefFit->SetNpx(2000);
    gRefFit->SetParameters(rr.c1, rr.c2);
    gRefFit->SetLineColor(kBlack);
    gRefFit->SetLineStyle(3); // dotted
    gRefFit->Draw("SAME");
  }

  if (r.fitted) {
    if (gCurrentFit) delete gCurrentFit;
    gCurrentFit = new TF1("curFit", TwFitFunc, dispLo, dispHi, 2);
    gCurrentFit->SetNpx(2000);
    gCurrentFit->SetParameters(r.c1, r.c2);
    gCurrentFit->SetLineColor(kRed);
    gCurrentFit->SetLineWidth(2);
    gCurrentFit->Draw("SAME");
  }

  TLine loLine(r.lo, h2->GetYaxis()->GetXmin(), r.lo, h2->GetYaxis()->GetXmax());
  TLine hiLine(r.hi, h2->GetYaxis()->GetXmin(), r.hi, h2->GetYaxis()->GetXmax());
  loLine.SetLineColor(kOrange + 2); loLine.SetLineStyle(2); loLine.DrawClone();
  hiLine.SetLineColor(kOrange + 2); hiLine.SetLineStyle(2); hiLine.DrawClone();

  TPaveText pt(0.5, 0.12, 0.89, 0.36, "NDC");
  pt.SetFillColor(kWhite); pt.SetTextAlign(12); pt.SetFillStyle(1001);
  pt.AddText(Form("Channel %d / %d: %s", gIdx + 1, (int)gChannels.size(), key.Data()));
  pt.AddText(Form("Entries = %.0f", h2->GetEntries()));
  pt.AddText(Form("Range = [%.1f, %.1f] mV", r.lo, r.hi));
  if (r.fitted) {
    pt.AddText(Form("c1 = %.3f #pm %.3f", r.c1, r.c1err));
    pt.AddText(Form("c2 = %.3f #pm %.3f", r.c2, r.c2err));
    pt.AddText(Form("#chi^{2}/NDF = %.2f", r.chi2ndf));
  } else pt.AddText("(not fitted)");
  pt.AddText(Form("source: %s", r.source.Data()));
  if (gRefResults.count(key)) pt.AddText("red = new fit, black dotted = stored/reference");
  pt.DrawClone();

  gPad->Modified(); gPad->Update();
}

// ===========================================================================
// Section 10: interaction
// ===========================================================================

void GoNext() { if (gIdx < (int)gChannels.size() - 1) { ++gIdx; gClickStage = 0; DrawChannel(); } }
void GoPrev() { if (gIdx > 0) { --gIdx; gClickStage = 0; DrawChannel(); } }

void Refit() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  FitRec &r = gResults[key];
  DoFit(ch, r.lo, r.hi, "manual");
  DrawChannel();
}

void ResetRange() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  gResults[key].lo = gTwFitRangeLow;
  gResults[key].hi = gTwFitRangeHigh;
  gClickStage = 0;
  Refit();
}

void UseReference() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  if (!gRefResults.count(key)) { printf("No reference result for %s\n", key.Data()); return; }
  FitRec &rr = gRefResults[key];
  DoFit(ch, rr.lo, rr.hi, "reference-confirmed");
  DrawChannel();
}

void GoToChannelPrompt() {
  printf("\nEnter channel # (1-%d), or a label like 1x.pos[7]: ", (int)gChannels.size());
  fflush(stdout);
  std::string line;
  std::getline(std::cin, line);
  TString input(line.c_str()); input = input.Strip(TString::kBoth);
  if (input.Length() == 0) { printf("(cancelled)\n"); return; }
  if (input.IsDigit()) {
    Int_t idx = input.Atoi();
    if (idx >= 1 && idx <= (int)gChannels.size()) { gIdx = idx - 1; DrawChannel(); return; }
    printf("Index out of range\n"); return;
  }
  for (size_t i = 0; i < gChannels.size(); i++)
    if (ChanLabel(gChannels[i]) == input) { gIdx = (int)i; DrawChannel(); return; }
  printf("No channel matches '%s'\n", input.Data());
}

void PauseSave() { SaveJSON(gRunTag); WriteLegacyParam(gRunTag); printf("[progress saved]\n"); }

void SaveFinish() {
  PauseSave();
  printf("[saved and finished -- rename %s/SHMS/HODO/phodo_TWcalib_%s.param to phodo_TWcalib.param for hcana]\n",
         gParamDir.Data(), gRunTag.Data());
}

void OnClick() {
  Int_t event = gPad->GetEvent();
  if (event != 11) return; // 11 = button-1 down
  Int_t px = gPad->GetEventX();
  Double_t x = gPad->AbsPixeltoX(px);
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  FitRec &r = gResults.count(key) ? gResults[key] : (gResults[key] = FitRec{gTwFitRangeLow, gTwFitRangeHigh});
  if (gClickStage == 0) {
    gPendingLo = x;
    gClickStage = 1;
    printf("  lo = %.1f mV (click again for hi)\n", x);
  } else {
    Double_t lo = TMath::Min(gPendingLo, x);
    Double_t hi = TMath::Max(gPendingLo, x);
    gClickStage = 0;
    DoFit(ch, lo, hi, "manual");
    DrawChannel();
  }
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "SHMS Time-Walk Fit Controls", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Previous channel");
  gControlBar->AddButton("Next >>", "GoNext();", "Next channel");
  gControlBar->AddButton("Refit", "Refit();", "Re-run the fit with the current range");
  gControlBar->AddButton("Reset Range", "ResetRange();", "Back to default fit range");
  gControlBar->AddButton("Use Reference", "UseReference();", "Copy reference run's range+result here");
  gControlBar->AddButton("Go to channel...", "GoToChannelPrompt();", "Jump to a channel by index or label");
  gControlBar->AddButton("Pause (save)", "PauseSave();", "Write progress, keep going");
  gControlBar->AddButton("Save && Finish", "SaveFinish();", "Write final outputs");
  gControlBar->Show();
}

// ===========================================================================
// Section 11: non-interactive batch mode
// ===========================================================================

void RunNonInteractive() {
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    Double_t lo = gTwFitRangeLow, hi = gTwFitRangeHigh;
    TString src = "default";
    if (gRefResults.count(key)) { lo = gRefResults[key].lo; hi = gRefResults[key].hi; src = "reference-auto"; }
    DoFit(ch, lo, hi, src);
  }
  gSystem->mkdir(gOutDir, kTRUE);
  TString pdfPath = Form("%s/timewalk_shms_%s_summary.pdf", gOutDir.Data(), gRunTag.Data());
  TCanvas c("c", "c", 900, 700);
  gCanvas = &c;
  Bool_t firstPage = kTRUE;
  for (size_t i = 0; i < gChannels.size(); i++) {
    gIdx = (int)i;
    DrawChannel();
    TString pagePath = pdfPath;
    if (firstPage) { pagePath += "("; firstPage = kFALSE; }
    c.Print(pagePath);
  }
  TString closePath = pdfPath + ")";
  c.Print(closePath);
  printf("Wrote %s\n", pdfPath.Data());
  if (!gCompareOnly) {
    SaveJSON(gRunTag);
    WriteLegacyParam(gRunTag);
  } else {
    printf("[compareOnly: PDF written for comparison, nothing saved -- current fit results and .param file left untouched]\n");
  }
}

// ===========================================================================
// Section 12: main entry point
// ===========================================================================

void timewalk_calib_app(Int_t run = 0, Int_t referenceRun = 0, TString referenceTag = "",
                              Bool_t nonInteractive = kFALSE,
                              TString runs = "", Double_t dpCut = 1.0e9,
                              TString rootFile = "", TString outDir = "./timewalk_qa",
                              TString paramDir = "../../PARAM", Double_t tdcThresh = 1200.0,
                              Double_t fitRangeLow = 20.0, Double_t fitRangeHigh = 300.0,
                              Bool_t compareOnly = kTRUE) {
  if (run == 0 && runs.Length() == 0) {
    printf("ERROR: must supply a run number, e.g. timewalk_calib_app(26107)\n"); return;
  }
  gReferenceRun = referenceRun;
  gReferenceTag = referenceTag;
  gOutDir = outDir; gParamDir = paramDir;
  gTdcThresh = tdcThresh; gTwFitRangeLow = fitRangeLow; gTwFitRangeHigh = fitRangeHigh;
  gDpCut = dpCut;
  gCompareOnly = compareOnly;
  if (compareOnly && referenceTag.Length() == 0 && referenceRun == 0)
    printf("WARNING: compareOnly=true but no referenceRun/referenceTag given -- there will be nothing to compare against.\n");

  gRunList.clear();
  if (runs.Length() > 0) {
    std::stringstream ss(runs.Data());
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      TString t(tok.c_str()); t = t.Strip(TString::kBoth);
      if (t.Length() > 0) gRunList.push_back(t.Atoi());
    }
    if (gRunList.empty()) { printf("ERROR: could not parse any run numbers from runs=\"%s\"\n", runs.Data()); return; }
    gRun = gRunList.front();
    gRunTag = (gRunList.size() > 1) ? Form("%d-%d", gRunList.front(), gRunList.back()) : Form("%d", gRunList.front());
  } else {
    gRun = run;
    gRunTag = Form("%d", run);
  }

  if (gRunList.empty()) {
    if (rootFile.Length() == 0) gRootFile = ResolvedRunPath(run);
    else gRootFile = rootFile;
  }

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(nonInteractive);

  BuildChannelList();
  BuildHistos();
  LoadReference();
  SeedFromOwnPriorSave();

  if (nonInteractive) { RunNonInteractive(); return; }

  gIdx = 0;
  gCanvas = new TCanvas("timewalkCanvas", "SHMS Time-Walk Calibration", 900, 700);
  DrawChannel();
  gCanvas->AddExec("dynamic", "OnClick()");
  gCanvas->Update();
  MakeControlBar();

  printf("\n=== run(s) %s -- %d channel(s) loaded (SHMS) ===\n", gRunTag.Data(), (int)gChannels.size());
  if (gDpCut < 1.0e8) printf("|P.gtr.dp| < %.2f cut applied.\n", gDpCut);
  printf("Click twice on the histogram: 1st = fit lo, 2nd = fit hi (auto-refits).\n");
  printf("Use the 'SHMS Time-Walk Fit Controls' panel to navigate and save.\n");
}
