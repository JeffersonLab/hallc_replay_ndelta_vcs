// timewalk_calib_app.C (HMS)
//
// Interactive click-to-fit tool for the HMS Hodoscope time-walk
// correction (Eq. 2 of hodo_calib.pdf):
//
//   f_TW(a) = c1 + 1 / (a / TDC_Thrs)^c2
//
// fitted to TDC-ADC Time vs. ADC Pulse Amplitude, per PMT (plane, side,
// paddle). Replaces the old two-script batch flow (timeWalkHistos.C then
// timeWalkCalib.C) with one interactive session: the raw histograms are
// built once from the coincidence replay tree (same event-loop logic as
// hms_hodo_calib/timeWalkHistos.C), then you step through every PMT,
// click twice on the amplitude axis to set the fit range (1st click = lo,
// 2nd click = hi), see the fit update live, and move on.
//
// The reference/trigger apparatus branches (FADC_TREF, T1/T2) live under
// "T.coin.*" in this experiment's coincidence replay tree, NOT "T.hms.*"
// -- this app always reads them with the "coin" tree prefix.
//
// NOTE on cuts (mirrors hms_hodo_calib/timeWalkHistos.C exactly): that
// script computes PID (H.cal.etot, H.cer.npeSum) and reference-apparatus
// (multiplicity/amplitude/time) cuts but leaves them COMMENTED OUT --
// only the error-flag cut, the hodo pulse-amplitude range, the ADC/TDC
// paddle-match, and a broad [-100,100] ns diff-time window are actually
// applied. This app defaults to that same (weaker) cut set for fidelity;
// set applyPidRefCuts=true to turn the PID+ref-apparatus cuts on if you
// decide you want them (see the main entry point's arguments).
//
// Default fit range is [50, 300] mV, and the fit itself (not the display)
// is additionally restricted to y (TDC-ADC diff time) < 40 ns via
// fitYLow/fitYHigh -- a sparse high-|y| noise streak, mostly at low
// amplitude, otherwise pulls a bad fit (wrong-sign c2, huge chi2/NDF).
// The full histogram is still drawn/extrapolated in the display.
//
// Usage (must run interactively -- NOT with -q, it waits for clicks):
//   root -l 'timewalk_calib_app.C(26138)'                        // -50V setting, no reference
//   root -l 'timewalk_calib_app.C(26138, 26927)'                 // + reference run 26927
//   root -l 'timewalk_calib_app.C(26138, 0, "0")'                // + reference = tag "0"
//
// Non-interactive (batch), 4th arg = true -- fits every channel with its
// default/seeded range, writes the summary PDF and param/json files, and
// exits. outputTag (5th arg) controls the .param filename directly --
// use this to get hhodo_TWcalib_0.param for the -50V run and
// hhodo_TWcalib_1.param for the nominal-HV run, decoupled from the run
// number itself:
//   root -l -b -q 'timewalk_calib_app.C(26138, 0, "", true, "0")'   // -50V setting  -> hhodo_TWcalib_0.param
//   root -l -b -q 'timewalk_calib_app.C(26927, 0, "", true, "1")'   // nominal HV    -> hhodo_TWcalib_1.param
//
// Compare-only batch mode: fits every channel fresh (red curve) against a
// loaded reference (black dotted curve) in the summary PDF, but NEVER
// writes JSON/.param output -- use this to sanity-check a new fit against
// the currently-stored calibration without touching it. compareOnly is
// the last positional arg; needs a referenceRun/referenceTag to compare
// against:
//   root -l -b -q 'timewalk_calib_app.C(26138, 0, "0", true, "", "", 1.0e9, "", "./timewalk_qa", "../../PARAM", 120.0, 50.0, 300.0, false, -1.0e8, 40.0, true)'
//
// Combining several runs into one TChain, with an optional H.gtr.dp cut
// to keep away from elastic-peak bias -- pass a comma-separated run list
// as the "runs" argument (6th positional arg) and it overrides the
// single "run" arg for building histograms; outputTag (if given) still
// controls the .param filename, otherwise it falls back to "<first>-<last>":
//   root -l -b -q 'timewalk_calib_app.C(0, 0, "", true, "0", "26138,26139", 10.0)'
//
// Verifying a saved batch fit interactively -- referenceTag points "Use
// Reference" / the gray overlay curve at a saved tag's result instead of
// a single reference run (useful once you've saved under "0"/"1" and
// want to sanity-check individual runs against it):
//   root -l 'timewalk_calib_app.C(26138, 0, "0")'
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
//   Save && Finish            -- write JSON + hhodo_TWcalib_<tag>.param
//                            under paramDir, close
//
// Outputs (under outDir, default "./timewalk_qa") -- filenames use a
// "tag": outputTag if given, else a single run number, or "<first>-<last>"
// when "runs" chains several runs together (the JSON's "runs" array
// lists every run actually used):
//   timewalk_hms_<tag>.json          -- full per-channel results + provenance
//   timewalk_hms_<tag>_summary.pdf   -- one page per channel, fit overlaid
// Outputs (under paramDir, default "../../PARAM"):
//   HMS/HODO/hhodo_TWcalib_<tag>.param -- staged with the tag in the
//   filename, same convention as reftime_cut_app.C's tcoin_<run>.param;
//   rename to hhodo_TWcalib.param for hcana to read it.
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
// Section 1: HMS constants (from hms_hodo_calib/timeWalkHistos.C and
//   hms_hodo_calib/timeWalkCalib.C)
// ===========================================================================

static const TString gPlaneNames[4] = {"1x", "1y", "2x", "2y"};
static const TString gSideNames[2]  = {"pos", "neg"};
static const Int_t   gNbars[4]      = {16, 10, 16, 10};
static const Int_t   nBarsMax       = 16; // fixed column padding, matches legacy .param layout (HMS: 16, not SHMS's 21)

// PMTs permanently off for this experiment -- excluded from the channel
// list entirely (no histogram booked, no fit attempted; WriteLegacyParam
// then naturally writes (0,0) for these, same convention as a paddle
// that was never calibrated). None known yet for HMS -- add entries here
// the same way the SHMS app's gPermanentlyOff was populated, e.g.
// {plane, side, paddle} with plane 0=1x,1=1y,2=2x,3=2y and side 0=pos,1=neg.
struct OffPmt { Int_t plane; Int_t side; Int_t paddle; };
static const std::vector<OffPmt> gPermanentlyOff = {
  // (none yet for HMS)
};

Bool_t IsPermanentlyOff(Int_t plane, Int_t side, Int_t paddle) {
  for (auto &o : gPermanentlyOff)
    if (o.plane == plane && o.side == side && o.paddle == paddle) return kTRUE;
  return kFALSE;
}

static const TString kRefAdcName = "hFADC_TREF_ROC1";
// hT1 is read (for completeness/monitoring) but the source script's fill
// loop always references hT2 as the actual TDC reference, for every
// plane and side (no per-plane exception like SHMS's 2y-pos T2 special
// case) -- index 1 (hT2) is what gets used below.
static const std::vector<TString> kRefTdcNames = {"hT1", "hT2"};
static const Int_t kRefTdcIdx = 1; // hT2
static const TString kPidBranch = "H.cal.etot";   static const Double_t kPidCutLow = 0.100;
static const TString kCerBranch = "H.cer.npeSum"; static const Double_t kCerCutLow = 1.5;
static const Double_t kHodoAmpLo = 25.0, kHodoAmpHi = 1000.0;
static const Double_t kRefAmpLo = 50.0, kRefAmpHi = 60.0;
static const Double_t kRefTimeLo = 210.0, kRefTimeHi = 225.0;
// Broad ADC-TDC diff-time sanity window -- IS actually applied in the
// source script (unlike the PID/ref-apparatus cuts above, which are
// computed there but left commented out).
static const Double_t kDiffTimeLo = -100.0, kDiffTimeHi = 100.0;

static const Double_t tdcChanToTime = 0.09766; // ns
static const Double_t adcChanToTime = 0.0625;  // ns

// ===========================================================================
// Section 2: channel bookkeeping
// ===========================================================================

struct Channel { Int_t plane; Int_t side; Int_t paddle; };

struct FitRec {
  Double_t lo = 50.0, hi = 300.0;
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
Double_t gTdcThresh = 120.0; // TW fit threshold, units of mV per hcana's convention for this param (see timeWalkCalib.C)
Double_t gTwFitRangeLow = 50.0, gTwFitRangeHigh = 300.0;
// y (TDC-ADC diff time) window applied ONLY to the fit itself, via a
// temporary TH2 Y-axis SetRangeUser around the Fit() call -- excludes a
// sparse high-|y| noise streak (mostly at low amplitude) from pulling the
// fit off the real, tightly-clustered population. Sentinels (+-1e8) mean
// "use the histogram's own booked axis limit", i.e. no extra cut.
// Display (DrawChannel) always restores the full axis afterward.
Double_t gFitYLo = -1.0e8, gFitYHi = 40.0;

std::vector<Int_t> gRunList;  // non-empty -> BuildHistos chains all these runs together
TString  gRunTag;             // output-file tag: single run number, or "first-last" when chained
TString  gOutputTag;          // if non-empty, overrides gRunTag for output filenames (e.g. "0"/"1" for HV setting)
Double_t gDpCut = 1.0e9;      // |H.gtr.dp| < gDpCut cut; default effectively off (elastic-peak avoidance / acceptance sanity cut)
Bool_t   gCompareOnly = kFALSE; // batch mode: write the summary PDF (new fit vs. reference) but never touch JSON/.param output
Bool_t   gApplyPidRefCuts = kFALSE; // off by default, matching the uploaded HMS script's actual behavior

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
//   -- same event-selection logic as hms_hodo_calib/timeWalkHistos.C
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

  // Book histograms. x: pulse amplitude (mV), wide enough for the HMS
  // amp cut ceiling (1000 mV) and default fit range up to 600 mV.
  // y: TDC-ADC diff time, sized to the [-100,100] ns window that's
  // actually applied below.
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    gHist[key] = new TH2F("h2_" + key,
        Form("TDC-ADC Time vs. Pulse Amp %s %s Paddle %d;Pulse Amplitude (mV);TDC-ADC Time (ns)",
             gPlaneNames[ch.plane].Data(), gSideNames[ch.side].Data(), ch.paddle),
        1000, 0, 1000, 2000, kDiffTimeLo, kDiffTimeHi);
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

  // PID branches (bound regardless; only enforced if gApplyPidRefCuts)
  Double_t pidVal = 0, cerVal = 0;
  T->SetBranchAddress(kPidBranch, &pidVal);
  T->SetBranchAddress(kCerBranch, &cerVal);

  // Optional momentum acceptance cut (e.g. |H.gtr.dp| < 10 to stay off the
  // elastic peak when combining runs). Off by default (gDpCut ~ 1e9).
  Double_t dpVal = 0;
  T->SetBranchAddress("H.gtr.dp", &dpVal);

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
      TString base = "H.hod." + gPlaneNames[ipl] + "." + gSideNames[isd];
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
  printf("\n=== Building HMS time-walk histograms: %lld events ===\n", nentries);
  if (gApplyPidRefCuts) printf("  (PID + reference-apparatus cuts ENABLED -- not the uploaded script's default)\n");
  for (Long64_t ievt = 0; ievt < nentries; ievt++) {
    T->GetEntry(ievt);
    if (ievt % 200000 == 0 && ievt != 0) printf("  %lld / %lld events...\n", ievt, nentries);

    if (TMath::Abs(dpVal) >= gDpCut) continue; // e.g. |dp| < 10 as an acceptance/garbage-data sanity cut

    // PID / reference apparatus cuts -- computed but NOT applied by the
    // uploaded hms_hodo_calib/timeWalkHistos.C; only enforced here if
    // gApplyPidRefCuts is explicitly turned on.
    if (gApplyPidRefCuts) {
      if (pidVal < kPidCutLow) continue;
      if (cerVal < kCerCutLow) continue;
      Bool_t refMultCut = (refAdcMult != 1.0);
      Bool_t refAmpCut  = (refAdcAmp < kRefAmpLo || refAdcAmp > kRefAmpHi);
      Double_t refTimeNs0 = refAdcTime * adcChanToTime;
      Bool_t refTimeCut = (refTimeNs0 < kRefTimeLo || refTimeNs0 > kRefTimeHi);
      if (refMultCut || refAmpCut || refTimeCut) continue;
    }
    Double_t refTimeNs = refAdcTime * adcChanToTime;
    // HMS always uses hT2 as the TDC reference, for every plane/side --
    // no per-plane exception like SHMS's 2y-pos special case.
    Double_t refTdcNs = refTdc[kRefTdcIdx] * tdcChanToTime;

    for (Int_t ipl = 0; ipl < 4; ipl++) {
      for (Int_t isd = 0; isd < 2; isd++) {
        TString hkey = Form("%s.%s", gPlaneNames[ipl].Data(), gSideNames[isd].Data());
        HodoArr *h = hodo[hkey];

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
            // Broad sanity window, actually applied in the source script
            // (unlike the PID/ref-apparatus cuts above).
            if (diff < kDiffTimeLo || diff > kDiffTimeHi) continue;
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
  // Temporarily restrict the histogram's Y (diff-time) axis range so the
  // fit ignores a sparse high-|y| noise streak -- ROOT's TH1/TH2::Fit only
  // uses bins within the current axis range. Restored right after, so
  // display (DrawChannel) always sees the full axis.
  Double_t yLoFull = h2->GetYaxis()->GetXmin(), yHiFull = h2->GetYaxis()->GetXmax();
  Double_t yLoFit = (gFitYLo > -1.0e7) ? gFitYLo : yLoFull;
  Double_t yHiFit = (gFitYHi <  1.0e7) ? gFitYHi : yHiFull;
  h2->GetYaxis()->SetRangeUser(yLoFit, yHiFit);
  TF1 fit("twFit", TwFitFunc, lo, hi, 2);
  fit.SetParameters(1.0, 1.0);
  TFitResultPtr res = h2->Fit(&fit, "SREQ0"); // 0 = don't draw here, caller draws
  h2->GetYaxis()->SetRangeUser(yLoFull, yHiFull); // restore full axis for display
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
//   hc1_Pos = v1x, v1y, v2x, v2y,
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

// Loads a legacy hhodo_TWcalib.param (vanilla or staged-with-tag)
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
  ParseParamBlock(in, "hc1_Pos", c1Pos);
  in.clear(); in.seekg(0); ParseParamBlock(in, "hc1_Neg", c1Neg);
  in.clear(); in.seekg(0); ParseParamBlock(in, "hc2_Pos", c2Pos);
  in.clear(); in.seekg(0); ParseParamBlock(in, "hc2_Neg", c2Neg);
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
  if (staged) return Form("%s/HMS/HODO/hhodo_TWcalib_%s.param", gParamDir.Data(), tag.Data());
  return Form("%s/HMS/HODO/hhodo_TWcalib.param", gParamDir.Data());
}

// The tag actually used for THIS session's output filenames: outputTag
// if given (e.g. "0"/"1" for HV setting), else the run-based gRunTag.
TString EffectiveTag() { return gOutputTag.Length() > 0 ? gOutputTag : gRunTag; }

void LoadReference() {
  // A referenceTag (e.g. "0"/"1", a saved-output tag) takes priority
  // over a plain referenceRun number when both are given. The literal
  // string "vanilla" explicitly loads PARAM/HMS/HODO/hhodo_TWcalib.param
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
// session's own previously-staged param file (a resumed session).
void SeedFromOwnPriorSave() {
  TString path = LegacyParamPath(EffectiveTag(), kTRUE);
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
  TString runsNote = gRunList.empty() ? Form("%d", gRun) : Form("%d chained files", (int)gRunList.size());
  out << Form(";HMS Hodoscopes Output Parameter File: %s (tag %s)", runsNote.Data(), tag.Data()) << std::endl << std::endl;
  out << "hTDC_threshold=" << gTdcThresh << ". ;units of mV " << std::endl << std::endl;

  auto writeBlock = [&](const TString &name, Bool_t isC1, Int_t side) {
    out << ";Param " << name << std::endl;
    out << ";           1x             1y             2x             2y " << std::endl;
    out << "h" << name << " = ";
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
  TString path = Form("%s/timewalk_hms_%s.json", gOutDir.Data(), tag.Data());
  std::ofstream out(path.Data());
  out << "{\n  \"tag\": \"" << tag << "\",\n  \"runs\": [";
  for (size_t i = 0; i < gRunList.size(); i++) out << (i ? ", " : "") << gRunList[i];
  if (gRunList.empty()) out << gRun;
  out << "],\n  \"dp_cut\": " << gDpCut << ",\n  \"apply_pid_ref_cuts\": " << (gApplyPidRefCuts ? "true" : "false")
      << ",\n  \"reference_run\": " << gReferenceRun << ",\n  \"reference_tag\": \"" << gReferenceTag
      << "\",\n  \"tdc_threshold\": " << gTdcThresh << ",\n  \"channels\": [\n";
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

  // The y-window actually used to restrict the FIT (not the display) --
  // see DoFit()'s temporary GetYaxis()->SetRangeUser(). Shown here as
  // horizontal lines so it's visible even though the fit itself doesn't
  // touch the drawn axis range.
  Double_t yLoFull = h2->GetYaxis()->GetXmin(), yHiFull = h2->GetYaxis()->GetXmax();
  Double_t yLoFit = (gFitYLo > -1.0e7) ? gFitYLo : yLoFull;
  Double_t yHiFit = (gFitYHi <  1.0e7) ? gFitYHi : yHiFull;
  TLine yLoLine(h2->GetXaxis()->GetXmin(), yLoFit, h2->GetXaxis()->GetXmax(), yLoFit);
  TLine yHiLine(h2->GetXaxis()->GetXmin(), yHiFit, h2->GetXaxis()->GetXmax(), yHiFit);
  yLoLine.SetLineColor(kGreen + 2); yLoLine.SetLineStyle(2); yLoLine.DrawClone();
  yHiLine.SetLineColor(kGreen + 2); yHiLine.SetLineStyle(2); yHiLine.DrawClone();

  TPaveText pt(0.5, 0.12, 0.89, 0.40, "NDC");
  pt.SetFillColor(kWhite); pt.SetTextAlign(12); pt.SetFillStyle(1001);
  pt.AddText(Form("Channel %d / %d: %s", gIdx + 1, (int)gChannels.size(), key.Data()));
  pt.AddText(Form("Entries = %.0f", h2->GetEntries()));
  pt.AddText(Form("x-range (orange) = [%.1f, %.1f] mV", r.lo, r.hi));
  pt.AddText(Form("y-window (green, fit only) = [%.1f, %.1f] ns", yLoFit, yHiFit));
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

void PauseSave() { SaveJSON(EffectiveTag()); WriteLegacyParam(EffectiveTag()); printf("[progress saved]\n"); }

void SaveFinish() {
  PauseSave();
  printf("[saved and finished -- rename %s/HMS/HODO/hhodo_TWcalib_%s.param to hhodo_TWcalib.param for hcana]\n",
         gParamDir.Data(), EffectiveTag().Data());
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
  gControlBar = new TControlBar("vertical", "HMS Time-Walk Fit Controls", 20, 20);
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
  TString pdfPath = Form("%s/timewalk_hms_%s_summary.pdf", gOutDir.Data(), EffectiveTag().Data());
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
    SaveJSON(EffectiveTag());
    WriteLegacyParam(EffectiveTag());
  } else {
    printf("[compareOnly: PDF written for comparison, nothing saved -- current fit results and .param file left untouched]\n");
  }
}

// ===========================================================================
// Section 12: main entry point
// ===========================================================================

void timewalk_calib_app(Int_t run = 0, Int_t referenceRun = 0, TString referenceTag = "",
                              Bool_t nonInteractive = kFALSE, TString outputTag = "",
                              TString runs = "", Double_t dpCut = 1.0e9,
                              TString rootFile = "", TString outDir = "./timewalk_qa",
                              TString paramDir = "../../PARAM", Double_t tdcThresh = 120.0,
                              Double_t fitRangeLow = 50.0, Double_t fitRangeHigh = 300.0,
                              Bool_t applyPidRefCuts = kFALSE,
                              Double_t fitYLow = 10.0, Double_t fitYHigh = 30.0,
                              Bool_t compareOnly = kTRUE) {//change it to kFalse to do the new calibration
  if (run == 0 && runs.Length() == 0) {
    printf("ERROR: must supply a run number, e.g. timewalk_calib_app(26138)\n"); return;
  }
  gReferenceRun = referenceRun;
  gReferenceTag = referenceTag;
  gOutputTag = outputTag;
  gOutDir = outDir; gParamDir = paramDir;
  gTdcThresh = tdcThresh; gTwFitRangeLow = fitRangeLow; gTwFitRangeHigh = fitRangeHigh;
  gFitYLo = fitYLow; gFitYHi = fitYHigh;
  gDpCut = dpCut;
  gApplyPidRefCuts = applyPidRefCuts;
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
  gCanvas = new TCanvas("timewalkCanvas", "HMS Time-Walk Calibration", 900, 700);
  DrawChannel();
  gCanvas->AddExec("dynamic", "OnClick()");
  gCanvas->Update();
  MakeControlBar();

  printf("\n=== run(s) %s (output tag: %s) -- %d channel(s) loaded (HMS) ===\n",
         gRunTag.Data(), EffectiveTag().Data(), (int)gChannels.size());
  if (gDpCut < 1.0e8) printf("|H.gtr.dp| < %.2f cut applied.\n", gDpCut);
  printf("Click twice on the histogram: 1st = fit lo, 2nd = fit hi (auto-refits).\n");
  printf("Use the 'HMS Time-Walk Fit Controls' panel to navigate and save.\n");
}
