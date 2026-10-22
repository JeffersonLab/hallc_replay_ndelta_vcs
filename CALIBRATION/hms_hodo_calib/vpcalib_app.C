// vpcalib_app.C (HMS)
//
// Interactive click-to-fit tool for the HMS Hodoscope cable-length /
// effective-velocity correction, replicating hms_hodo_calib/fitHodoCalib.C's
// velocity/cable section:
//
//   TW-Corr. Time Diff.(y_trk) = (1/v_p) * y_trk + b0
//
// IMPORTANT -- the HMS script is structurally SIMPLER than the SHMS one
// this app's sibling (../shms_hodo_calib/vpcalib_app.C) was built from:
//   * NO cross-plane single-hit-in-all-4-planes requirement for this
//     stage -- each paddle is selected independently: PID + its own
//     <100ns hard time cut + its own +-1sigma mean-time window. (The
//     "require all 4 planes" logic only shows up later, in the lambda
//     stage, via hod_nhits==1 -- see lambda_calib_app.C.)
//   * NO velocity sanity-check / constrained-refit fallback in the
//     original -- a failed or degenerate fit (slope exactly 1, i.e. the
//     TF1's untouched default parameter) just falls back to a FLAT
//     v_p=15.0 cm/ns for every plane, cable=0, no bounds checking at
//     all. The sample hhodo_Vpcalib.param you get without any sanity
//     check can come out with wildly unphysical velocities (negative,
//     hundreds of cm/ns) -- this app still defaults to matching that
//     free-fit-only behavior, but ALSO offers the same interactive
//     "Toggle Fit Mode" / constrained-fit safety net the SHMS app has,
//     as a genuinely useful addition given how easily the free fit goes
//     wrong here. It's opt-in, not automatic, so it doesn't change what
//     batch mode produces relative to the original script unless you
//     turn on autoConstrainOutOfRange.
//   * Fit range is PER PLANE TYPE, not uniform: X-planes (1x, 2x) use
//     [-35, 35] cm, Y-planes (1y, 2y) use [-50, 50] cm.
//
// Depends on GoodPosTdcTimeWalkCorr/GoodNegTdcTimeWalkCorr, i.e. a
// replay done with an up-to-date hhodo_TWcalib.param.
//
// PID cuts (H.cal.etracknorm > 0.7, H.cer.npeSum > 0.7, H.dc.ntrack > 0)
// match fitHodoCalib.C, and can be turned off independently: applyCalCut,
// applyCerCut, applyTrackCut. If the calorimeter isn't reconstructed for
// a given run (H.cal.etrack itself reads 0, not just the momentum-
// normalized etracknorm), set applyCalCut=false and leave the other two
// on -- don't drop all PID just because one branch is dead.
//
// SIGMA (hhodo_PosSigma/NegSigma, a per-paddle timing resolution) is
// computed here the same way as in the SHMS app. The original HMS
// script's own sigma code references DiffDistTWCorr before it is ever
// assigned, so there is nothing reliable to replicate; the SHMS
// definition is used instead. PosSigma and NegSigma are the same array
// written twice. For each paddle, with v_p/cable as written to the
// .param file:
//   DiffDistTWCorr = v_p * (0.5*(TW_neg - TW_pos) - cable)     [cm]
//   residual       = DiffDistTWCorr - TrackPos                  [cm]
//   sigma          = StdDev(residual) / (2 * v_p)               [ns]
// using the same events that fill the fit histograms. A paddle whose
// v_p AND cable are identical (to the 6 decimals of the .param file) to
// the loaded reference (referenceTag, "vanilla" for the untagged file)
// is not recomputed; the reference's own sigma is copied instead. If no
// reference is loaded, nothing can be copied and every paddle is
// recomputed.
//
// SIGMAONLY: like compareOnly, loads this tag's already-saved
// hhodo_Vpcalib_<tag>.param and trusts its v_p/cable (no fitting), but
// DOES save: sigma is recomputed (or copied from the reference where
// v_p/cable match it) and written together with the unchanged
// v_p/cable. Use this to update sigma alone after fitting.
// (compareOnly and sigmaOnly together is nonsensical -- compareOnly wins.)
//
// COMPAREONLY: a pure review mode, not a fitting mode. compareOnly=true
// never calls DoFit() at all -- it asserts a calibration for this exact
// run/tag was ALREADY committed (hhodo_Vpcalib_<tag>.param must already
// exist under paramDir; the app aborts with an error if it doesn't) and
// draws that saved calibration in red alongside a loaded reference
// (referenceTag, "vanilla" for the untagged file) in black dotted, using
// this run's own histograms for context. Nothing is fit, nothing is
// saved -- Pause/Save&&Finish are disabled for the whole session.
//
// Usage (must run interactively -- NOT with -q, it waits for clicks):
//   root -l 'vpcalib_app.C(23859)'
//   root -l 'vpcalib_app.C(23859, 0, "0")'     // + reference = a saved tag
//
// Recompute sigma only, leaving velFit/cableFit as already saved (last arg):
//   root -l -b -q 'vpcalib_app.C(26138, 0, "vanilla", true, "0", "", 10.0, "", "./vpcalib_qa", "../../PARAM", true, 0.0, false, true, true, false, true)'
//
// Non-interactive (batch), 4th arg = true:
//   root -l -b -q 'vpcalib_app.C(23859, 0, "", true)'
//
// Review an already-saved calibration against vanilla, no fitting:
//   root -l -b -q 'vpcalib_app.C(26138, 0, "vanilla", true, "0", "", 10.0, "", "./vpcalib_qa", "../../PARAM", true, 0.0, false, true, true, true)'
//
// Combining runs / |H.gtr.dp| cut, "runs" is the 6th positional arg:
//   root -l 'vpcalib_app.C(0, 0, "", false, "", "23859,23860", 10.0)'
//
// Controls:
//   Click twice on the histogram: 1st click = fit range lo, 2nd = hi
//     (track position, cm). Fits FREE (slope, intercept) by default.
//   << Prev / Next >>     -- step through (plane, paddle) channels
//   Refit                  -- re-run the current fit mode with the current range
//   Toggle Fit Mode         -- FREE <-> CONSTRAINED (slope fixed to 1/v_p,
//                            only the offset floats) -- not in the
//                            original script, offered because the free
//                            fit has no safety net here (see note above)
//   Reset Range              -- back to the plane-type default range, FREE mode
//   Use Reference             -- RE-FIT this channel's own histogram,
//                            seeded with the reference's range+mode
//   Copy Reference Value       -- ADOPT the reference's vp/cable directly,
//                            no local fit (for low-statistics paddles;
//                            referenceTag="vanilla" loads the existing
//                            PARAM/HMS/HODO/hhodo_Vpcalib.param)
//   Go to channel...           -- console prompt, index or label "1x[7]"
//   Pause (save) / Save && Finish
//
// Outputs (under outDir, default "./vpcalib_qa"):
//   vpcalib_hms_<tag>.json          -- per-paddle results + provenance
//   vpcalib_hms_<tag>_summary.pdf   -- one page per paddle, fit overlaid
// Outputs (under paramDir, default "../../PARAM"):
//   HMS/HODO/hhodo_Vpcalib_<tag>.param -- staged with the tag in the
//   filename; rename to hhodo_Vpcalib.param for hcana to read it.
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
// Section 1: HMS constants (from hms_hodo_calib/fitHodoCalib.C)
// ===========================================================================

static const TString gPlaneNames[4] = {"1x", "1y", "2x", "2y"};
static const Int_t   gNbars[4]      = {16, 10, 16, 10}; // matches maxPMT[] in fitHodoCalib.C
static const Int_t   nBarsMax       = 16;

// Per-plane-type fit range, matching fit1x/fit1y/fit2x/fit2y's own
// hardcoded TF1 ranges in the original script.
Double_t DefaultFitLo(Int_t plane) { return (plane == 0 || plane == 2) ? -35.0 : -50.0; }
Double_t DefaultFitHi(Int_t plane) { return (plane == 0 || plane == 2) ?  35.0 :  50.0; }

// No permanently-off list supplied for HMS yet -- same mechanism as the
// other HMS apps, ready whenever one is given.
struct OffPmt { Int_t plane; Int_t side; Int_t paddle; };
static const std::vector<OffPmt> gPermanentlyOff = {
  // (none yet for HMS)
};
Bool_t IsPermanentlyOff(Int_t plane, Int_t side, Int_t paddle) {
  for (auto &o : gPermanentlyOff)
    if (o.plane == plane && o.side == side && o.paddle == paddle) return kTRUE;
  return kFALSE;
}

static const TString kPidBranch1 = "H.cal.etracknorm"; static const Double_t kPidCut1Low = 0.7;
static const TString kPidBranch2 = "H.cer.npeSum";     static const Double_t kPidCut2Low = 0.7;
static const TString kPidBranch3 = "H.dc.ntrack";      static const Double_t kPidCut3Low = 0.0;
static const Double_t kHardTimeCutHi = 100.0; // TdcTimeWalkCorr < 100 ns, both sides (uniform for HMS, unlike SHMS's 200/125 split)
static const Double_t kNSig = 1.0;

// Flat fallback velocity, matching the original's literal "hhodo_velArr = 15.0"
// (no per-plane velSet, no bounds window, in the original).
static const Double_t kVelFallback = 15.0;
// Not part of the original script -- used only for the OPTIONAL
// sanity-check flag / constrained-fit target this app adds (see note above).
static const Double_t kVelVar = 3.0;

TString ResolvedRunPath(Int_t run) {
  // See note in shms_hodo_calib/vpcalib_app.C -- switch back once the
  // normal replay dir is current:
  return Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
  // return Form("/volatile/hallc/alphaE/ndelta_vcs2/selected_runs_before_hodo_calibration/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
}

// ===========================================================================
// Section 2: channel bookkeeping
// ===========================================================================

struct Channel { Int_t plane; Int_t side; Int_t paddle; }; // side unused (per-paddle quantity), kept for label consistency

struct FitRec {
  Double_t lo = -35.0, hi = 35.0; // overwritten per-plane-type when seeded
  Bool_t   constrained = kFALSE;
  Double_t vp = 0.0, cable = 0.0;
  Double_t vpErr = 0.0, cableErr = 0.0;
  Double_t chi2ndf = 0.0;
  Bool_t   fitted = kFALSE;
  TString  source = "default";
};

TString ChanLabel(const Channel &c) { return Form("%s[%d]", gPlaneNames[c.plane].Data(), c.paddle); }

// ===========================================================================
// Section 3: globals
// ===========================================================================

Int_t    gRun = 0, gReferenceRun = 0;
TString  gReferenceTag;
TString  gRootFile, gOutDir = "./vpcalib_qa", gParamDir = "../../PARAM";

std::vector<Int_t> gRunList;
TString  gRunTag, gOutputTag;
Double_t gDpCut = 1.0e9; // |H.gtr.dp| < gDpCut, off by default (not part of the original script's cuts)
Bool_t   gCompareOnly = kFALSE; // batch mode: write the summary PDF (new fit vs. reference) but never touch JSON/.param output
Bool_t   gApplyCalCut = kTRUE;   // H.cal.etracknorm cut -- turn off if the calorimeter isn't reconstructed for this run
Bool_t   gApplyCerCut = kTRUE;   // H.cer.npeSum cut
Bool_t   gApplyTrackCut = kTRUE; // H.dc.ntrack cut
Double_t gMinEntriesForFit = 0.0;
Bool_t   gAutoConstrainOutOfRange = kFALSE; // opt-in; original script has no such fallback at all
Bool_t   gSigmaOnly = kFALSE;   // load this tag's ALREADY-SAVED velFit/cableFit (no fitting, like compareOnly),
                                // but DOES save -- recomputes+writes sigma, leaving velFit/cableFit as loaded

// Sigma (timing resolution) support. BuildHistos() caches, for every
// event that fills a paddle's fit histogram, the half time difference
// and track position; ComputeSigma() turns them into sigma right before
// any save, using whatever v_p/cable the session holds at that point.
struct SigCand { Float_t halfDiff; Float_t pos; };
std::vector<SigCand> gSigCands[4][nBarsMax];
Double_t gSigma[4][nBarsMax] = {{0}};
Bool_t   gSigmaComputed = kFALSE;
// The reference file's own PosSigma, loaded alongside velFit/cableFit.
Double_t gReferenceSigma[4][nBarsMax] = {{0}};
Bool_t   gHaveReferenceSigma = kFALSE;
// Sigma blocks carried over from an existing hhodo_Vpcalib.param -- only
// used as a fallback if ComputeSigma() has nothing to write.
TString  gSigmaPosBlockRaw, gSigmaNegBlockRaw;
Bool_t   gHaveCarriedSigma = kFALSE;

std::vector<Channel> gChannels;
std::map<TString, TH2F*> gHist;
std::map<TString, FitRec> gResults;
std::map<TString, FitRec> gRefResults;

Int_t     gIdx = 0;
TCanvas  *gCanvas = nullptr;
TControlBar *gControlBar = nullptr;
TF1      *gCurrentFit = nullptr;
TF1      *gRefFit = nullptr;
Int_t     gClickStage = 0;
Double_t  gPendingLo = 0.0;

TString EffectiveTag() { return gOutputTag.Length() > 0 ? gOutputTag : gRunTag; }

// ===========================================================================
// Section 4: channel list
// ===========================================================================

void BuildChannelList() {
  gChannels.clear();
  Int_t nSkipped = 0;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad) || IsPermanentlyOff(ipl, 1, ipad)) { ++nSkipped; continue; }
      gChannels.push_back(Channel{ipl, 0, ipad});
    }
  if (nSkipped) printf("Excluded %d permanently-off paddle(s) from the channel list.\n", nSkipped);
}

// ===========================================================================
// Section 5: build raw histograms -- single pass, per-paddle selection
//   (NO cross-plane single-hit requirement, unlike the SHMS app)
// ===========================================================================

void BuildHistos() {
  TFile *f = nullptr;
  TChain *chain = nullptr;
  TTree *T = nullptr;

  if (!gRunList.empty()) {
    chain = new TChain("T");
    for (auto r : gRunList) {
      Int_t nAdded = chain->Add(ResolvedRunPath(r));
      if (nAdded == 0) printf("WARNING: could not add run %d to chain\n", r);
    }
    if (chain->GetNtrees() == 0) { printf("ERROR: no files added to chain\n"); return; }
    T = chain;
  } else {
    f = new TFile(gRootFile, "READ");
    if (!f || f->IsZombie()) { printf("ERROR: could not open %s\n", gRootFile.Data()); return; }
    T = (TTree*)f->Get("T");
    if (!T) { printf("ERROR: no tree 'T' in %s\n", gRootFile.Data()); return; }
  }

  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    gHist[key] = new TH2F("h2_" + key,
        Form("TW-Corr. Time Diff vs. Track Pos %s Paddle %d;Hodoscope Track Position from Center (cm);Time-Walk Corrected TDC Time Difference (ns)",
             gPlaneNames[ch.plane].Data(), ch.paddle),
        400, -60, 60, 400, -15, 15);
    gHist[key]->SetDirectory(nullptr);
  }

  Double_t pid1 = 0, pid2 = 0, pid3 = 0, dpVal = 0;
  T->SetBranchAddress(kPidBranch1, &pid1);
  T->SetBranchAddress(kPidBranch2, &pid2);
  T->SetBranchAddress(kPidBranch3, &pid3);
  T->SetBranchAddress("H.gtr.dp", &dpVal);

  Double_t twPos[4][nBarsMax], twNeg[4][nBarsMax];
  Double_t trackX[4], trackY[4];
  for (Int_t ipl = 0; ipl < 4; ipl++) {
    TString base = "H.hod." + gPlaneNames[ipl];
    T->SetBranchAddress(base + ".GoodPosTdcTimeWalkCorr", twPos[ipl]);
    T->SetBranchAddress(base + ".GoodNegTdcTimeWalkCorr", twNeg[ipl]);
    T->SetBranchAddress(base + ".TrackXPos", &trackX[ipl]);
    T->SetBranchAddress(base + ".TrackYPos", &trackY[ipl]);
  }

  Long64_t nentries = T->GetEntries();
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) gSigCands[ipl][ipad].clear();

  // Pass 1: per-(plane,paddle) mean/stddev of the TW-average time
  // (fitHodoCalib.C's "FIRST PASS"). Note the original ALSO has a
  // redundant inline "hcal_etrkNorm>0.7" check here identical to the
  // outer PID cut -- harmless, not replicated separately.
  std::map<TString, TH1F*> twAvg;
  for (auto &ch : gChannels)
    twAvg[ChanLabel(ch)] = new TH1F("twavg_" + ChanLabel(ch), "", 4000, -200, 200);

  printf("\n=== HMS vpcalib pass 1/2 (mean/stddev): %lld events ===\n", nentries);
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    Bool_t pidOk = (!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low);
    if (!pidOk) continue;
    for (Int_t ipl = 0; ipl < 4; ipl++)
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        TString key = ChanLabel(Channel{ipl, 0, ipad});
        if (!twAvg.count(key)) continue;
        Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
        if (p < kHardTimeCutHi && n < kHardTimeCutHi) twAvg[key]->Fill(0.5 * (p + n));
      }
  }
  std::map<TString, Double_t> meanOf, stdOf;
  for (auto &kv : twAvg) { meanOf[kv.first] = kv.second->GetMean(); stdOf[kv.first] = kv.second->GetStdDev(); delete kv.second; }

  // Pass 2: fill directly -- EACH PADDLE INDEPENDENTLY, no requirement
  // that the other 3 planes also have a clean single hit (that
  // requirement only appears later, in the lambda stage).
  printf("=== HMS vpcalib pass 2/2 (fill): %lld events ===\n", nentries);
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    Bool_t pidOk = (!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low);
    if (!pidOk) continue;

    for (Int_t ipl = 0; ipl < 4; ipl++) {
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        TString key = ChanLabel(Channel{ipl, 0, ipad});
        if (!gHist.count(key)) continue;
        Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
        if (!(p < kHardTimeCutHi && n < kHardTimeCutHi)) continue;
        Double_t avg = 0.5 * (p + n);
        if (!meanOf.count(key)) continue;
        if (!(avg > meanOf[key] - kNSig * stdOf[key] && avg < meanOf[key] + kNSig * stdOf[key])) continue;
        Double_t halfDiff = 0.5 * (n - p);
        Double_t pos = (ipl == 0 || ipl == 2) ? trackY[ipl] : trackX[ipl]; // 1x/2x -> TrackYPos, 1y/2y -> TrackXPos
        gHist[key]->Fill(pos, halfDiff);
        gSigCands[ipl][ipad - 1].push_back(SigCand{(Float_t)halfDiff, (Float_t)pos});
      }
    }
  }
  printf("=== Done: %lld events processed ===\n\n", nentries);

  if (f) f->Close();
  if (chain) delete chain;
}

// ===========================================================================
// Section 6: fitting
// ===========================================================================

void DoFit(const Channel &ch, Double_t lo, Double_t hi, Bool_t constrained, const TString &source) {
  TString key = ChanLabel(ch);
  TH2F *h2 = gHist[key];
  FitRec &r = gResults[key];
  r.lo = lo; r.hi = hi; r.constrained = constrained; r.source = source;
  if (!h2 || h2->GetEntries() == 0) {
    printf("WARNING: no entries for %s -- cannot fit\n", key.Data());
    r.fitted = kFALSE;
    return;
  }
  Int_t status;
  if (!constrained) {
    TF1 fit("vpFit", "[0]*x + [1]", lo, hi);
    fit.SetParameters(1.0, 0.0); // matches the original's TF1 default start (slope=1) -- see the "==1" degenerate-fit check below
    TFitResultPtr res = h2->Fit(&fit, "SREQ0");
    status = res;
    Double_t slope = fit.GetParameter(0);
    // Same degenerate-fit detection as fitHodoCalib.C: status==-1 OR the
    // slope came back exactly at the untouched default (1.0, i.e. v_p=1)
    // both mean "this didn't really fit".
    if (status == -1 || slope == 1.0) {
      r.fitted = kFALSE; // caller (DoFit's invoker) decides the flat fallback at write-time, same as WriteLegacyParam below
      r.vp = kVelFallback; r.vpErr = 0.0; r.cable = 0.0; r.cableErr = 0.0; r.chi2ndf = 0.0;
      printf("Could not fit %s (status=%d, slope=%.3f) -- falling back to v_p=%.1f, cable=0\n", key.Data(), status, slope, kVelFallback);
      return;
    }
    r.fitted = kTRUE;
    r.vp = (slope != 0.0) ? 1.0 / slope : 0.0;
    r.vpErr = (slope != 0.0) ? TMath::Abs(fit.GetParError(0) / (slope * slope)) : 0.0;
    r.cable = fit.GetParameter(1); r.cableErr = fit.GetParError(1);
    r.chi2ndf = (fit.GetNDF() > 0) ? fit.GetChisquare() / fit.GetNDF() : 0.0;
  } else {
    // Not in the original script -- opt-in safety net (see header note).
    TF1 fit("vpFitConstrained", Form("(1./%f)*x + [0]", kVelFallback), lo, hi);
    fit.SetParameter(0, 0.0);
    TFitResultPtr res = h2->Fit(&fit, "SREQ0");
    status = res;
    r.fitted = (status == 0);
    r.vp = kVelFallback; r.vpErr = 0.0;
    r.cable = fit.GetParameter(0); r.cableErr = fit.GetParError(0);
    r.chi2ndf = (fit.GetNDF() > 0) ? fit.GetChisquare() / fit.GetNDF() : 0.0;
    if (status != 0) printf("WARNING: constrained fit did not converge cleanly for %s (status=%d)\n", key.Data(), status);
  }
}

Bool_t VelocityInRange(const Channel &ch, Double_t vp) {
  return vp > (kVelFallback - kVelVar) && vp < (kVelFallback + kVelVar);
}

// ===========================================================================
// Section 7: legacy hcana .param I/O
// ===========================================================================

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
    while (std::getline(ss, tok, ',')) { try { vals.push_back(std::stod(tok)); } catch (...) {} }
  };
  pushNums(rest);
  while ((Int_t)vals.size() < 4 * nBarsMax) {
    std::streampos p = in.tellg();
    if (!std::getline(in, line)) break;
    TString t(line.c_str()); t = t.Strip(TString::kBoth);
    if (t.Length() == 0 || t.BeginsWith(";") || t.Contains("=")) { in.seekg(p); break; }
    pushNums(line);
  }
  for (Int_t i = 0; i < (Int_t)vals.size() && i < 4 * nBarsMax; i++) out[i % 4][i / 4] = vals[i];
}

TString ExtractRawBlock(const TString &path, const TString &key) {
  std::ifstream in(path.Data());
  if (!in.is_open()) return "";
  std::string line;
  bool found = false;
  while (std::getline(in, line)) {
    if (line.find(key.Data()) != std::string::npos && line.find('=') != std::string::npos) { found = true; break; }
  }
  if (!found) return "";
  TString block = line.c_str(); block += "\n";
  while (true) {
    std::streampos p = in.tellg();
    if (!std::getline(in, line)) break;
    TString t(line.c_str()); t = t.Strip(TString::kBoth);
    if (t.Length() == 0 || t.BeginsWith(";") || t.Contains("=")) break;
    block += line.c_str(); block += "\n";
  }
  return block;
}

TString LegacyParamPath(const TString &tag, Bool_t staged) {
  if (staged) return Form("%s/HMS/HODO/hhodo_Vpcalib_%s.param", gParamDir.Data(), tag.Data());
  return Form("%s/HMS/HODO/hhodo_Vpcalib.param", gParamDir.Data());
}

void LoadLegacyParam(const TString &path) {
  std::ifstream in(path.Data());
  if (!in.is_open()) { printf("  (no legacy param file at %s)\n", path.Data()); return; }
  Double_t vel[4][nBarsMax] = {{0}}, cable[4][nBarsMax] = {{0}};
  Double_t sig[4][nBarsMax] = {{0}};
  ParseParamBlock(in, "hhodo_velFit", vel);
  in.clear(); in.seekg(0); ParseParamBlock(in, "hhodo_cableFit", cable);
  in.clear(); in.seekg(0); ParseParamBlock(in, "hhodo_PosSigma", sig);
  gHaveReferenceSigma = (ExtractRawBlock(path, "hhodo_PosSigma").Length() > 0);
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) gReferenceSigma[ipl][ipad] = sig[ipl][ipad];
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad) || IsPermanentlyOff(ipl, 1, ipad)) continue;
      Channel ch{ipl, 0, ipad};
      FitRec r;
      r.vp = vel[ipl][ipad - 1]; r.cable = cable[ipl][ipad - 1];
      r.lo = DefaultFitLo(ipl); r.hi = DefaultFitHi(ipl);
      r.constrained = !VelocityInRange(ch, r.vp);
      r.fitted = kTRUE; r.source = "reference-file";
      gRefResults[ChanLabel(ch)] = r;
    }
  printf("  loaded reference param: %s\n", path.Data());
}

// compareOnly mode: loads THIS run/tag's own already-staged calibration
// (not the reference/vanilla) directly into gResults, with no fitting at
// all. Returns false if no staged file exists for EffectiveTag() -- the
// caller aborts in that case.
Bool_t LoadOwnSavedResults() {
  TString path = LegacyParamPath(EffectiveTag(), kTRUE);
  std::ifstream test(path.Data());
  if (!test.is_open()) return kFALSE;
  test.close();
  Double_t vel[4][nBarsMax] = {{0}}, cable[4][nBarsMax] = {{0}};
  std::ifstream in(path.Data());
  ParseParamBlock(in, "hhodo_velFit", vel);
  in.clear(); in.seekg(0); ParseParamBlock(in, "hhodo_cableFit", cable);
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad) || IsPermanentlyOff(ipl, 1, ipad)) continue;
      Channel ch{ipl, 0, ipad};
      FitRec r;
      r.vp = vel[ipl][ipad - 1]; r.cable = cable[ipl][ipad - 1];
      r.lo = DefaultFitLo(ipl); r.hi = DefaultFitHi(ipl);
      r.constrained = !VelocityInRange(ch, r.vp);
      r.fitted = kTRUE; r.source = "saved";
      gResults[ChanLabel(ch)] = r;
    }
  printf("  loaded this run's own saved calibration: %s\n", path.Data());
  return kTRUE;
}

void LoadReference() {
  if (gReferenceTag.Length() == 0 && gReferenceRun == 0) return;
  if (gReferenceTag == "vanilla") { LoadLegacyParam(LegacyParamPath("", kFALSE)); return; }
  TString refTag = (gReferenceTag.Length() > 0) ? gReferenceTag : Form("%d", gReferenceRun);
  TString staged = LegacyParamPath(refTag, kTRUE);
  std::ifstream test(staged.Data());
  if (test.is_open()) { test.close(); LoadLegacyParam(staged); }
  else LoadLegacyParam(LegacyParamPath("", kFALSE));
}

void SeedFromOwnPriorSave() {
  TString path = LegacyParamPath(EffectiveTag(), kTRUE);
  std::ifstream test(path.Data());
  if (!test.is_open()) return;
  test.close();
  std::map<TString, FitRec> saved = gRefResults;
  // LoadLegacyParam() also loads PosSigma into the reference-sigma slots;
  // keep the real reference's sigma, not this tag's own prior save.
  Double_t savedSigma[4][nBarsMax];
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) savedSigma[ipl][ipad] = gReferenceSigma[ipl][ipad];
  Bool_t savedHaveSigma = gHaveReferenceSigma;
  gRefResults.clear();
  LoadLegacyParam(path);
  for (auto &kv : gRefResults)
    if (!gResults.count(kv.first)) { gResults[kv.first] = kv.second; gResults[kv.first].source = "resumed"; }
  gRefResults = saved;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) gReferenceSigma[ipl][ipad] = savedSigma[ipl][ipad];
  gHaveReferenceSigma = savedHaveSigma;
}

// Carry over PosSigma/NegSigma from whatever param file already exists
// (this tag's own staged file if present, else vanilla). Must run BEFORE
// WriteLegacyParam() opens (and truncates) the output file.
void LoadCarriedSigma() {
  TString ownPath = LegacyParamPath(EffectiveTag(), kTRUE);
  std::ifstream ownTest(ownPath.Data());
  TString sourcePath = ownTest.is_open() ? ownPath : LegacyParamPath("", kFALSE);
  ownTest.close();
  gSigmaPosBlockRaw = ExtractRawBlock(sourcePath, "hhodo_PosSigma");
  gSigmaNegBlockRaw = ExtractRawBlock(sourcePath, "hhodo_NegSigma");
  gHaveCarriedSigma = (gSigmaPosBlockRaw.Length() > 0 && gSigmaNegBlockRaw.Length() > 0);
  if (gHaveCarriedSigma) printf("Carried over PosSigma/NegSigma from %s\n", sourcePath.Data());
  else printf("No existing PosSigma/NegSigma found -- will write a 1.0 placeholder unless sigma is computed.\n");
}

// The v_p/cable that WriteLegacyParam() will write for this channel
// (flat fallback if the channel was never touched).
void EffectiveVpCable(const TString &key, Double_t &vp, Double_t &cable) {
  vp = kVelFallback; cable = 0.0;
  if (gResults.count(key)) { vp = gResults[key].vp; cable = gResults[key].cable; }
}

// True if this channel's v_p AND cable match the loaded reference to
// within the precision of the .param files (6 decimals). Such a paddle
// has the reference's timing calibration, so the reference's sigma is
// carried over instead of being recomputed.
const Double_t kRefMatchTol = 1.0e-6;
Bool_t SameAsReference(const TString &key) {
  if (!gHaveReferenceSigma || !gRefResults.count(key)) return kFALSE;
  Double_t vp, cable;
  EffectiveVpCable(key, vp, cable);
  const FitRec &rr = gRefResults[key];
  return (TMath::Abs(vp - rr.vp) < kRefMatchTol && TMath::Abs(cable - rr.cable) < kRefMatchTol);
}

// Computes hhodo_PosSigma/NegSigma from the events cached in
// BuildHistos() (no second tree read), using the v_p/cable this session
// holds for each paddle right now (fresh fit, resumed, saved, or
// reference-copied). Paddles identical to the reference take the
// reference's sigma. See the header for the formula.
void ComputeSigma() {
  Int_t nComputed = 0, nCopiedFromRef = 0, nUndetermined = 0;
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    Int_t ipl = ch.plane, ipad = ch.paddle - 1;
    gSigma[ipl][ipad] = 0.0;
    Bool_t copyRef = (gResults.count(key) && gResults[key].source.BeginsWith("reference-copied")) || SameAsReference(key);
    if (copyRef && gHaveReferenceSigma) {
      gSigma[ipl][ipad] = gReferenceSigma[ipl][ipad];
      ++nCopiedFromRef;
      continue;
    }
    Double_t vp, cable;
    EffectiveVpCable(key, vp, cable);
    if (vp == 0.0 || gSigCands[ipl][ipad].empty()) { ++nUndetermined; continue; }
    TH1F h("resid", "", 2000, -120, 80);
    h.SetDirectory(nullptr);
    for (auto &c : gSigCands[ipl][ipad]) h.Fill(vp * (c.halfDiff - cable) - c.pos);
    gSigma[ipl][ipad] = h.GetStdDev() / (2.0 * vp);
    ++nComputed;
  }
  gSigmaComputed = (nComputed > 0 || nCopiedFromRef > 0);
  printf("Sigma: computed for %d paddle(s) from this run's own data, copied from the reference for %d paddle(s) (v_p/cable identical to it), %d paddle(s) left undetermined (0.0), out of %d total.\n",
         nComputed, nCopiedFromRef, nUndetermined, (int)gChannels.size());
}

void WriteLegacyParam(const TString &tag) {
  TString outPath = LegacyParamPath(tag, kTRUE);
  gSystem->mkdir(gSystem->DirName(outPath), kTRUE);
  std::ofstream out(outPath.Data());
  out << "; HMS Hodoscope Parameter File Containing propagation velocities per paddle " << std::endl;
  out << "; and signal cable time diff. offsets per paddle " << std::endl;
  out << Form("; tag %s ", tag.Data()) << std::endl << std::endl;

  auto writeBlock = [&](const TString &name, Bool_t isVp) {
    out << ";" << name << std::endl;
    out << name << " = ";
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) {
      for (Int_t ipl = 0; ipl < 4; ipl++) {
        Double_t v = isVp ? kVelFallback : 0.0; // matches fitHodoCalib.C's flat fallback
        if (ipad < gNbars[ipl]) {
          Channel ch{ipl, 0, ipad + 1};
          TString key = ChanLabel(ch);
          // Unlike the SHMS app, a DEGENERATE fit here already has
          // r.vp/r.cable set to the correct fallback by DoFit() itself
          // (r.fitted==false in that case) -- so write whatever's in
          // gResults whenever the channel was touched at all, not only
          // when fitted==true.
          if (gResults.count(key)) v = isVp ? gResults[key].vp : gResults[key].cable;
        }
        out << std::fixed << v;
        if (ipl != 3) out << ", ";
      }
      out << "," << std::endl;
    }
    out << std::endl;
  };
  writeBlock("hhodo_velFit", kTRUE);
  writeBlock("hhodo_cableFit", kFALSE);

  if (gSigmaComputed) {
    out << ";PMTs Time Diff. Sigma Parameters (computed this session)" << std::endl;
    auto writeSigmaBlock = [&](const TString &name) {
      out << name << " = ";
      for (Int_t ipad = 0; ipad < nBarsMax; ipad++) {
        for (Int_t ipl = 0; ipl < 4; ipl++) {
          Double_t v = (ipad < gNbars[ipl] && !IsPermanentlyOff(ipl, 0, ipad + 1)) ? gSigma[ipl][ipad] : 0.0;
          out << std::fixed << v;
          if (ipl != 3) out << ", ";
        }
        out << "," << std::endl;
      }
      out << std::endl;
    };
    // PosSigma and NegSigma are the same array written twice.
    writeSigmaBlock("hhodo_PosSigma");
    writeSigmaBlock("hhodo_NegSigma");
  } else if (gHaveCarriedSigma) {
    out << ";PMTs Time Diff. Sigma Parameters (carried over, not recomputed this session)" << std::endl;
    out << gSigmaPosBlockRaw.Data() << std::endl << gSigmaNegBlockRaw.Data() << std::endl;
  } else {
    out << ";PMTs Time Diff. Sigma Parameters (placeholder)" << std::endl;
    out << "hhodo_PosSigma = ";
    for (Int_t i = 0; i < nBarsMax; i++) out << "1.000000, 1.000000, 1.000000, 1.000000," << std::endl;
    out << std::endl << "hhodo_NegSigma = ";
    for (Int_t i = 0; i < nBarsMax; i++) out << "1.000000, 1.000000, 1.000000, 1.000000," << std::endl;
  }
  out.close();
  printf("Wrote %s\n", outPath.Data());
}

// ===========================================================================
// Section 8: JSON I/O
// ===========================================================================

void SaveJSON(const TString &tag) {
  gSystem->mkdir(gOutDir, kTRUE);
  TString path = Form("%s/vpcalib_hms_%s.json", gOutDir.Data(), tag.Data());
  std::ofstream out(path.Data());
  out << "{\n  \"tag\": \"" << tag << "\",\n  \"runs\": [";
  for (size_t i = 0; i < gRunList.size(); i++) out << (i ? ", " : "") << gRunList[i];
  if (gRunList.empty()) out << gRun;
  out << "],\n  \"dp_cut\": " << gDpCut << ",\n  \"reference_run\": " << gReferenceRun
      << ",\n  \"reference_tag\": \"" << gReferenceTag << "\",\n  \"channels\": [\n";
  Bool_t first = kTRUE;
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    if (!gResults.count(key)) continue;
    FitRec &r = gResults[key];
    if (!first) out << ",\n";
    first = kFALSE;
    out << "    {\"channel\": \"" << key << "\", \"plane\": \"" << gPlaneNames[ch.plane]
        << "\", \"paddle\": " << ch.paddle << ", \"lo\": " << r.lo << ", \"hi\": " << r.hi
        << ", \"constrained\": " << (r.constrained ? "true" : "false")
        << ", \"vp\": " << r.vp << ", \"vp_err\": " << r.vpErr
        << ", \"cable\": " << r.cable << ", \"cable_err\": " << r.cableErr
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

  FitRec &r = gResults.count(key) ? gResults[key] : (gResults[key] = FitRec{DefaultFitLo(ch.plane), DefaultFitHi(ch.plane)});

  Double_t xLoFull = h2->GetXaxis()->GetXmin(), xHiFull = h2->GetXaxis()->GetXmax();

  if (gRefResults.count(key)) {
    FitRec &rr = gRefResults[key];
    if (gRefFit) delete gRefFit;
    gRefFit = new TF1("refFit", "[0]*x + [1]", xLoFull, xHiFull);
    gRefFit->SetNpx(1000);
    Double_t refSlope = (rr.vp != 0.0) ? 1.0 / rr.vp : 0.0;
    gRefFit->SetParameters(refSlope, rr.cable);
    gRefFit->SetLineColor(kGray + 2);
    gRefFit->SetLineStyle(2);
    gRefFit->Draw("SAME");
  }

  // Skip the "new" curve when it's a verbatim copy of the reference
  // (compareOnly's no-fitting-at-all mode, or the low-stat fallback) --
  // it would sit exactly on top of the black dotted reference curve.
  if ((r.vp != 0.0 || r.cable != 0.0) && !r.source.BeginsWith("reference-copied")) {
    if (gCurrentFit) delete gCurrentFit;
    gCurrentFit = new TF1("curFit", "[0]*x + [1]", xLoFull, xHiFull);
    gCurrentFit->SetNpx(1000);
    Double_t slope = (r.vp != 0.0) ? 1.0 / r.vp : 0.0;
    gCurrentFit->SetParameters(slope, r.cable);
    gCurrentFit->SetLineColor(r.constrained ? kMagenta + 1 : (r.fitted ? kRed : kGray + 1));
    gCurrentFit->SetLineWidth(2);
    gCurrentFit->Draw("SAME");
  }

  TLine loLine(r.lo, h2->GetYaxis()->GetXmin(), r.lo, h2->GetYaxis()->GetXmax());
  TLine hiLine(r.hi, h2->GetYaxis()->GetXmin(), r.hi, h2->GetYaxis()->GetXmax());
  loLine.SetLineColor(kOrange + 2); loLine.SetLineStyle(2); loLine.DrawClone();
  hiLine.SetLineColor(kOrange + 2); hiLine.SetLineStyle(2); hiLine.DrawClone();

  TPaveText pt(0.5, 0.12, 0.89, 0.42, "NDC");
  pt.SetFillColor(kWhite); pt.SetTextAlign(12); pt.SetFillStyle(1001);
  pt.AddText(Form("Channel %d / %d: %s", gIdx + 1, (int)gChannels.size(), key.Data()));
  pt.AddText(Form("Entries = %.0f", h2->GetEntries()));
  pt.AddText(Form("Range = [%.1f, %.1f] cm", r.lo, r.hi));
  pt.AddText(r.constrained ? "Mode: CONSTRAINED (v_p fixed)" : "Mode: FREE");
  if (r.source == "saved") {
    pt.AddText(Form("v_p = %.3f cm/ns (this run's saved calibration)", r.vp));
    pt.AddText(Form("cable offset = %.3f ns (saved)", r.cable));
  } else if (r.source.BeginsWith("reference-copied")) {
    pt.AddText(Form("v_p = %.3f cm/ns (copied, no fit performed)", r.vp));
    pt.AddText(Form("cable offset = %.3f ns (copied)", r.cable));
  } else if (r.fitted) {
    if (r.constrained) pt.AddText(Form("v_p = %.3f cm/ns (fixed)", r.vp));
    else pt.AddText(Form("v_p = %.3f #pm %.3f cm/ns", r.vp, r.vpErr));
    pt.AddText(Form("cable offset = %.3f #pm %.3f ns", r.cable, r.cableErr));
    pt.AddText(Form("#chi^{2}/NDF = %.2f", r.chi2ndf));
    if (!r.constrained && !VelocityInRange(ch, r.vp)) {
      TText *warn = pt.AddText(Form("v_p far from %.1f #pm %.1f cm/ns -- try Constrained", kVelFallback, kVelVar));
      warn->SetTextColor(kRed);
    }
  } else {
    pt.AddText(Form("(fit failed -- flat fallback v_p=%.1f, cable=0)", kVelFallback));
  }
  pt.AddText(Form("source: %s", r.source.Data()));
  if (gCompareOnly && gRefResults.count(key)) pt.AddText("red = this run's saved calibration, black dotted = reference");
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
  DoFit(ch, r.lo, r.hi, r.constrained, "manual");
  DrawChannel();
}

void ToggleFitMode() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  FitRec &r = gResults[key];
  DoFit(ch, r.lo, r.hi, !r.constrained, "manual");
  DrawChannel();
}

void ResetRange() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  gResults[key].lo = DefaultFitLo(ch.plane); gResults[key].hi = DefaultFitHi(ch.plane);
  gClickStage = 0;
  DoFit(ch, DefaultFitLo(ch.plane), DefaultFitHi(ch.plane), kFALSE, "manual");
  DrawChannel();
}

void UseReference() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  if (!gRefResults.count(key)) { printf("No reference result for %s\n", key.Data()); return; }
  FitRec &rr = gRefResults[key];
  DoFit(ch, rr.lo, rr.hi, rr.constrained, "reference-confirmed");
  DrawChannel();
}

void CopyReferenceValue() {
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  if (!gRefResults.count(key)) { printf("No reference result for %s\n", key.Data()); return; }
  FitRec &rr = gRefResults[key];
  FitRec &r = gResults[key];
  r.lo = rr.lo; r.hi = rr.hi; r.constrained = rr.constrained;
  r.vp = rr.vp; r.vpErr = 0.0; r.cable = rr.cable; r.cableErr = 0.0;
  r.chi2ndf = 0.0; r.fitted = kTRUE; r.source = "reference-copied";
  DrawChannel();
}

void GoToChannelPrompt() {
  printf("\nEnter channel # (1-%d), or a label like 1x[7]: ", (int)gChannels.size());
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

void PauseSave() {
  if (gCompareOnly) { printf("compareOnly is on -- saving is disabled for this session.\n"); return; }
  ComputeSigma();
  SaveJSON(EffectiveTag()); WriteLegacyParam(EffectiveTag()); printf("[progress saved]\n");
}

void SaveFinish() {
  if (gCompareOnly) { printf("compareOnly is on -- saving is disabled for this session.\n"); return; }
  PauseSave();
  printf("[saved and finished -- rename %s/HMS/HODO/hhodo_Vpcalib_%s.param to hhodo_Vpcalib.param for hcana]\n",
         gParamDir.Data(), EffectiveTag().Data());
}

void OnClick() {
  Int_t event = gPad->GetEvent();
  if (event != 11) return;
  Int_t px = gPad->GetEventX();
  Double_t x = gPad->AbsPixeltoX(px);
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  FitRec &r = gResults.count(key) ? gResults[key] : (gResults[key] = FitRec{DefaultFitLo(ch.plane), DefaultFitHi(ch.plane)});
  if (gClickStage == 0) {
    gPendingLo = x;
    gClickStage = 1;
    printf("  lo = %.1f cm (click again for hi)\n", x);
  } else {
    Double_t lo = TMath::Min(gPendingLo, x);
    Double_t hi = TMath::Max(gPendingLo, x);
    gClickStage = 0;
    DoFit(ch, lo, hi, r.constrained, "manual");
    DrawChannel();
  }
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "HMS Vp/Cable Fit Controls", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Previous channel");
  gControlBar->AddButton("Next >>", "GoNext();", "Next channel");
  gControlBar->AddButton("Refit", "Refit();", "Re-run the current fit mode with the current range");
  gControlBar->AddButton("Toggle Fit Mode", "ToggleFitMode();", "Switch FREE <-> CONSTRAINED (fixed v_p)");
  gControlBar->AddButton("Reset Range", "ResetRange();", "Back to the plane-type default range, FREE mode");
  gControlBar->AddButton("Use Reference", "UseReference();", "Re-fit here, seeded with reference's range+mode");
  gControlBar->AddButton("Copy Reference Value", "CopyReferenceValue();", "Adopt reference's vp/cable directly, no local fit");
  gControlBar->AddButton("Go to channel...", "GoToChannelPrompt();", "Jump to a channel by index or label");
  gControlBar->AddButton("Pause (save)", "PauseSave();", "Write progress, keep going");
  gControlBar->AddButton("Save && Finish", "SaveFinish();", "Write final outputs");
  gControlBar->Show();
}

// ===========================================================================
// Section 11: non-interactive batch mode
// ===========================================================================

void RunNonInteractive() {
  if (!gCompareOnly && !gSigmaOnly) {
    for (auto &ch : gChannels) {
      TString key = ChanLabel(ch);
      TH2F *h2 = gHist[key];
      if (gMinEntriesForFit > 0.0 && gRefResults.count(key) && h2 && h2->GetEntries() < gMinEntriesForFit) {
        FitRec &rr = gRefResults[key];
        FitRec &r = gResults[key];
        r.lo = rr.lo; r.hi = rr.hi; r.constrained = rr.constrained;
        r.vp = rr.vp; r.vpErr = 0.0; r.cable = rr.cable; r.cableErr = 0.0;
        r.chi2ndf = 0.0; r.fitted = kTRUE; r.source = "reference-copied-lowstat";
        continue;
      }
      Double_t lo = DefaultFitLo(ch.plane), hi = DefaultFitHi(ch.plane);
      Bool_t constrained = kFALSE;
      TString src = "default";
      if (gRefResults.count(key)) { lo = gRefResults[key].lo; hi = gRefResults[key].hi; constrained = gRefResults[key].constrained; src = "reference-auto"; }
      DoFit(ch, lo, hi, constrained, src);
      // Opt-in only -- the original script has NO automatic fallback like
      // this; it just accepts whatever the free fit gives (or the flat
      // 15.0 default on an outright degenerate fit).
      if (gAutoConstrainOutOfRange && !constrained && gResults[key].fitted && !VelocityInRange(ch, gResults[key].vp)) {
        DoFit(ch, lo, hi, kTRUE, src + "+auto-constrained");
      }
    }
  }
  // else: gCompareOnly or gSigmaOnly -- gResults was already fully
  // populated by LoadOwnSavedResults() before this function was ever
  // called; no fitting happens here at all.
  gSystem->mkdir(gOutDir, kTRUE);
  TString pdfPath = Form("%s/vpcalib_hms_%s_summary.pdf", gOutDir.Data(), EffectiveTag().Data());
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
    ComputeSigma();
    SaveJSON(EffectiveTag());
    WriteLegacyParam(EffectiveTag());
  } else {
    printf("[compareOnly: PDF written from this run's saved calibration vs. reference -- no fitting performed, .param file left untouched]\n");
  }
}

// ===========================================================================
// Section 12: main entry point
// ===========================================================================

void vpcalib_app(Int_t run = 0, Int_t referenceRun = 0, TString referenceTag = "",
                  Bool_t nonInteractive = kFALSE, TString outputTag = "",
                  TString runs = "", Double_t dpCut = 1.0e9,
                  TString rootFile = "", TString outDir = "./vpcalib_qa",
                  TString paramDir = "../../PARAM",
                  Bool_t applyCalCut = kTRUE, Double_t minEntriesForFit = 0.0,
                  Bool_t autoConstrainOutOfRange = kFALSE,
                  Bool_t applyCerCut = kTRUE, Bool_t applyTrackCut = kTRUE,
                  Bool_t compareOnly = kFALSE, Bool_t sigmaOnly = kFALSE) {
  if (run == 0 && runs.Length() == 0) {
    printf("ERROR: must supply a run number, e.g. vpcalib_app(23859)\n"); return;
  }
  gReferenceRun = referenceRun;
  gReferenceTag = referenceTag;
  gOutputTag = outputTag;
  gOutDir = outDir; gParamDir = paramDir;
  gDpCut = dpCut;
  gApplyCalCut = applyCalCut;
  gApplyCerCut = applyCerCut;
  gApplyTrackCut = applyTrackCut;
  gMinEntriesForFit = minEntriesForFit;
  gAutoConstrainOutOfRange = autoConstrainOutOfRange;
  gCompareOnly = compareOnly;
  gSigmaOnly = sigmaOnly;
  if (compareOnly && sigmaOnly)
    printf("WARNING: both compareOnly and sigmaOnly are set -- compareOnly wins, nothing will be saved.\n");
  if (sigmaOnly && !compareOnly && referenceTag.Length() == 0 && referenceRun == 0)
    printf("WARNING: sigmaOnly without a reference (e.g. referenceTag=\"vanilla\") -- no sigma can be copied, every paddle will be recomputed.\n");
  if (compareOnly && referenceTag.Length() == 0 && referenceRun == 0) {
    gReferenceTag = "vanilla";
    printf("compareOnly=true with no referenceRun/referenceTag given -- defaulting to referenceTag=\"vanilla\" "
           "(the currently-saved ../../PARAM/HMS/HODO/hhodo_Vpcalib.param) so there's something to compare against.\n");
  }

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

  if (compareOnly || sigmaOnly) {
    // Assert the calibration was already committed -- compareOnly and
    // sigmaOnly never fit from scratch; both load this tag's existing
    // velFit/cableFit. compareOnly never saves; sigmaOnly saves the
    // recomputed sigma with velFit/cableFit unchanged.
    TString ownPath = LegacyParamPath(EffectiveTag(), kTRUE);
    std::ifstream ownTest(ownPath.Data());
    if (!ownTest.is_open()) {
      printf("ERROR: compareOnly/sigmaOnly requires an existing staged calibration for tag \"%s\", but none found at %s.\n"
             "Run vpcalib_app interactively, or in batch mode with compareOnly=false, to produce it first.\n",
             EffectiveTag().Data(), ownPath.Data());
      return;
    }
    ownTest.close();
  }

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(nonInteractive);

  BuildChannelList();
  BuildHistos();
  LoadReference();
  if (compareOnly || sigmaOnly) LoadOwnSavedResults(); // already confirmed to exist above; no fitting, ever
  else SeedFromOwnPriorSave();
  LoadCarriedSigma();

  if (nonInteractive) { RunNonInteractive(); return; }

  gIdx = 0;
  gCanvas = new TCanvas("vpcalibCanvas", "HMS Vp/Cable Calibration", 900, 700);
  DrawChannel();
  gCanvas->AddExec("dynamic", "OnClick()");
  gCanvas->Update();
  MakeControlBar();

  printf("\n=== run(s) %s (output tag: %s) -- %d channel(s) loaded (HMS) ===\n",
         gRunTag.Data(), EffectiveTag().Data(), (int)gChannels.size());
  printf("Click twice on the histogram: 1st = fit lo, 2nd = fit hi (auto-refits, keeps current mode).\n");
  printf("Use the 'HMS Vp/Cable Fit Controls' panel to navigate, toggle fit mode, and save.\n");
}
