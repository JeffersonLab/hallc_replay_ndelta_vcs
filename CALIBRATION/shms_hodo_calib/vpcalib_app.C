// vpcalib_app.C (SHMS)
//
// Interactive click-to-fit tool for the SHMS Hodoscope cable-length /
// effective-velocity correction (the "Cable Time Corrections" step of
// hodo_calib.pdf, Eq. 10-11 and Fig. 5-6):
//
//   Delta t_Corr.TW(y_trk) = (1/v_p) * y_trk + b0
//
// fitted per paddle to (track position at the plane, half the TW-corrected
// pos-minus-neg TDC time difference), using the SAME two-pass single-hit
// selection logic as shms_hodo_calib/fitHodoCalib.C's velocity/cable
// section: v_p = 1/slope, b0 = intercept.
//
// This depends on GoodPosTdcTimeWalkCorr / GoodNegTdcTimeWalkCorr, i.e. a
// replay done with an up-to-date phodo_TWcalib.param -- run this AFTER
// the time-walk-corrected replay, not on raw data.
//
// SIGMA (phodo_PosSigma/NegSigma, a per-paddle timing resolution) is
// computed here too -- NOT as a separate app/tree-read. It was
// originally built as its own sigma_calib_app.C, but that just
// re-derived the identical single-hit-per-plane candidates this app's
// own BuildHistos() already computes, then re-opened the ROOT files a
// second time for a quantity that fitHodoCalib.C itself gets almost for
// free as a byproduct of the same event loop -- so it's folded in
// directly: BuildHistos()'s existing pass caches those same candidates
// (EventCand), and ComputeSigma() turns them into PosSigma/NegSigma
// right before any save, using whatever v_p/cable this session currently
// holds per paddle (fresh fit, resumed, or reference-copied). Formula,
// confirmed from fitHodoCalib.C's actual Part-3 block (not recalled from
// memory -- an earlier read of this had gotten tangled up with stale
// commented-out code from an older revision):
//   DiffDistTWCorr = v_p * 0.5*(TW_neg - 2*cable - TW_pos)   [cm]
//   residual       = DiffDistTWCorr - TrackPos               [cm]
//   sigma          = StdDev(residual) / (2 * v_p)            [ns]
// phodo_PosSigma and phodo_NegSigma are the SAME array written twice in
// the original script -- there is no separate pos/neg sigma quantity,
// despite the naming. Skipped entirely under compareOnly (see below).
// If an existing phodo_Vpcalib.param has PosSigma/NegSigma and this
// session computes nothing new for some reason, those are carried over
// unchanged rather than overwritten with a 1.0 placeholder.
//
// PID cuts (P.cal.etracknorm, P.hgcer.npeSum, P.dc.ntrack) are applied by
// default, matching fitHodoCalib.C, and can be turned off independently:
// applyCalCut, applyCerCut, applyTrackCut. Calorimeter and Cherenkov are
// electron-ID cuts -- turn both off (applyCalCut=false, applyCerCut=false)
// for a proton-arm run, where those signals don't mean what the fixed
// thresholds assume; applyTrackCut is species-independent and usually
// fine to leave on. BuildHistos() also prints a per-cut event-survival count and
// flags any empty channel with whether it had standalone candidates that
// got cut only by the cross-plane single-hit requirement (see below) --
// use that before assuming PID is the culprit for an empty paddle.
//
// NOTE ON EMPTY PADDLES: this app requires exactly one good paddle in
// ALL FOUR planes simultaneously (fitHodoCalib.C's "track matching lite"
// technique), which the time-walk calibration never required -- a paddle
// with plenty of TW-calibration statistics can still come up empty here
// if events that hit it rarely also produce a clean single hit in the
// other three planes (common at the geometric edges of a plane's
// acceptance). That's a structural difference from PID cuts, not a bug;
// the pass-2 diagnostic above tells you which case you're looking at.
//
// COMPAREONLY: a pure review mode, not a fitting mode. compareOnly=true
// never calls DoFit() at all -- it asserts a calibration for this exact
// run/tag was ALREADY committed (phodo_Vpcalib_<tag>.param must already
// exist under paramDir; the app aborts with an error if it doesn't) and
// draws that saved calibration in red alongside a loaded reference
// (referenceTag, "vanilla" for the untagged file) in black dotted, using
// this run's own histograms for context. Nothing is fit, nothing is
// saved -- Pause/Save&&Finish are disabled for the whole session.
//
// SIGMAONLY: like compareOnly, loads this tag's already-saved
// velFit/cableFit instead of fitting -- but, UNLIKE compareOnly, DOES
// save: ComputeSigma() runs against the cached candidates using those
// loaded (not freshly fit) v_p/cable values, and the merged output is
// written with velFit/cableFit exactly as loaded (unchanged) plus a
// freshly recomputed PosSigma/NegSigma. Use this to update sigma alone
// after the fact, without re-touching the velocity/cable calibration.
// (compareOnly and sigmaOnly together is nonsensical -- compareOnly wins.)
//
// Usage (must run interactively -- NOT with -q, it waits for clicks):
//   root -l 'vpcalib_app.C(26107)'
//   root -l 'vpcalib_app.C(26107, 0, "0")'     // + reference = a saved tag
//
// Non-interactive (batch), 4th arg = true:
//   root -l -b -q 'vpcalib_app.C(26107, 0, "", true)'
//
// Review an already-saved calibration against vanilla, no fitting:
//   root -l -b -q 'vpcalib_app.C(0, 0, "vanilla", true, "", "26483,26484,26485,26486,26487,26488", 10.0, "", "./vpcalib_qa", "../../PARAM", -40.0, 40.0, true, 0.0, true, true, true)'
//
// Recompute sigma only, leaving velFit/cableFit as already saved:
//   root -l -b -q 'vpcalib_app.C(0, 0, "", true, "", "26483,26484,26485,26486,26487,26488", 10.0, "", "./vpcalib_qa", "../../PARAM", -40.0, 40.0, true, 0.0, true, true, false, true)'
//
// Combining several runs into one TChain, with an optional |P.gtr.dp| cut
// to keep away from elastic-peak bias -- pass a comma-separated run list
// as the "runs" argument (6th positional arg, right after outputTag) and
// it overrides the single "run" arg for building histograms; the output
// tag defaults to "<first>-<last>" when "runs" has more than one entry
// (referenceTag can point at that same tag to self-seed/overlay a prior
// save of this combined run, or at any other saved tag):
//   root -l 'vpcalib_app.C(0, 0, "26483-26488", false, "", "26483,26484,26485,26486,26487,26488", 10.0)'
//
// Non-interactive batch fit of the same combined runs, straight to PDF:
//   root -l -b -q 'vpcalib_app.C(0, 0, "26483-26488", true, "", "26483,26484,26485,26486,26487,26488", 10.0)'
//
// Controls:
//   Click twice on the histogram: 1st click = fit range lo, 2nd = hi
//     (track position, cm). Fits FREE (slope, intercept) by default.
//   << Prev / Next >>     -- step through (plane, side, paddle) channels
//   Refit                  -- re-run the current fit mode with the current range
//   Toggle Fit Mode         -- switch between FREE (2-param) and
//                            CONSTRAINED (slope fixed to the plane's
//                            nominal 1/v_p, only the offset floats) --
//                            use this when the free fit's v_p lands
//                            outside the plane's [vSet-3, vSet+3] cm/ns
//                            sanity window (flagged in the stat box)
//   Reset Range              -- back to the default [-40, 40] cm, FREE mode
//   Use Reference             -- RE-FIT this channel's own histogram,
//                            seeded with the reference's range+mode
//   Copy Reference Value       -- ADOPT the reference's vp/cable directly,
//                            no local fit at all -- use this for a
//                            low-statistics paddle where your own data
//                            can't be trusted (e.g. an existing
//                            PARAM/SHMS/HODO/phodo_Vpcalib.param from a
//                            prior experiment/setup: point referenceTag
//                            at the literal string "vanilla" to load it)
//   Go to channel...           -- console prompt, index or label "1x.pos[7]"
//   Pause (save) / Save && Finish
//
// For BATCH mode, minEntriesForFit does the same "copy, don't fit" thing
// automatically: any channel with fewer than that many histogram entries
// is copied from the loaded reference instead of fit (needs a reference
// loaded via referenceRun/referenceTag; default 0 = off, always fit).
//
// Outputs (under outDir, default "./vpcalib_qa"):
//   vpcalib_shms_<tag>.json          -- full per-channel results + provenance
//   vpcalib_shms_<tag>_summary.pdf   -- one page per channel, fit overlaid
// Outputs (under paramDir, default "../../PARAM"):
//   SHMS/HODO/phodo_Vpcalib_<tag>.param -- staged with the tag in the
//   filename; rename to phodo_Vpcalib.param for hcana to read it.
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
// Section 1: SHMS constants (from shms_hodo_calib/fitHodoCalib.C)
// ===========================================================================

static const TString gPlaneNames[4] = {"1x", "1y", "2x", "2y"};
static const TString gSideNames[2]  = {"pos", "neg"};
static const Int_t   gNbars[4]      = {13, 13, 14, 21}; // matches maxPMT[] in fitHodoCalib.C
static const Int_t   nBarsMax       = 21;

// Same permanently-off list as timewalk_calib_app.C.
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

static const TString kPidBranch1 = "P.cal.etracknorm"; static const Double_t kPidCut1Low = 0.7;
static const TString kPidBranch2 = "P.hgcer.npeSum";   static const Double_t kPidCut2Low = 0.5;
static const TString kPidBranch3 = "P.dc.ntrack";      static const Double_t kPidCut3Low = 0.0;
static const Double_t kHardTimeCutHi = 200.0; // TdcTimeWalkCorr < 200 ns, both sides
static const Double_t kNSig = 1.0;            // +-nSig*StdDev window on the mean TW-avg time

// Nominal per-plane propagation velocity (cm/ns) and the +-tolerance
// window used to sanity-check the free-fit's 1/slope.
static const Double_t kVelSet[4] = {15.75, 15.75, 15.75, 14.0};
static const Double_t kVelVar = 3.0;

// ===========================================================================
// Section 2: channel bookkeeping
// ===========================================================================

struct Channel { Int_t plane; Int_t side; Int_t paddle; }; // side kept for symmetry w/ other apps; fit is per-paddle (side-independent)

struct FitRec {
  Double_t lo = -40.0, hi = 40.0;
  Bool_t   constrained = kFALSE; // FALSE = free 2-param fit, TRUE = slope fixed to 1/kVelSet[plane]
  Double_t vp = 0.0, cable = 0.0;
  Double_t vpErr = 0.0, cableErr = 0.0;
  Double_t chi2ndf = 0.0;
  Bool_t   fitted = kFALSE;
  TString  source = "default"; // default | manual | reference-confirmed | reference-file | resumed
};

TString ChanLabel(const Channel &c) {
  // Per-paddle (not per-side) quantity, but keep the same "plane.side[pad]"
  // label shape as the other apps for consistency; side is always "pos"
  // here since only one channel per paddle exists.
  return Form("%s.pos[%d]", gPlaneNames[c.plane].Data(), c.paddle);
}

// ===========================================================================
// Section 3: globals
// ===========================================================================

Int_t    gRun = 0, gReferenceRun = 0;
TString  gReferenceTag;
TString  gRootFile, gOutDir = "./vpcalib_qa", gParamDir = "../../PARAM";
Double_t gFitRangeLow = -40.0, gFitRangeHigh = 40.0;
Double_t gMinEntriesForFit = 0.0; // if >0 and a reference is loaded, a channel with fewer histogram entries
                                   // than this is copied from the reference (no local fit) in batch mode --
                                   // same idea as the interactive "Copy Reference Value" button

std::vector<Int_t> gRunList;
TString  gRunTag;
TString  gOutputTag;
Double_t gDpCut = 1.0e9; // |P.gtr.dp| < gDpCut, off by default (not part of the original script's cuts)
Bool_t   gApplyCalCut = kTRUE;   // P.cal.etracknorm cut -- turn off for a proton-arm run (electron-ID cuts don't apply to protons)
Bool_t   gApplyCerCut = kTRUE;   // P.hgcer.npeSum cut -- same reasoning
Bool_t   gApplyTrackCut = kTRUE; // P.dc.ntrack cut (species-independent track quality, usually fine to leave on)
Bool_t   gCompareOnly = kFALSE; // batch mode: write the summary PDF (new fit vs. reference) but never touch JSON/.param output
Bool_t   gSigmaOnly = kFALSE;   // load this tag's ALREADY-SAVED velFit/cableFit (no fitting, like compareOnly),
                                 // but DOES save -- recomputes+writes sigma, leaving velFit/cableFit as loaded

std::vector<Channel> gChannels;

// Sigma (timing resolution) support -- merged in from what used to be a
// separate sigma_calib_app.C. Cached during BuildHistos()'s existing
// single-hit-per-plane pass (same candidates, same cuts, no second tree
// read), then turned into PosSigma/NegSigma right before any save, once
// this session's v_p/cable values are known.
struct EventCand { Int_t pad[4]; Double_t twPos[4], twNeg[4]; Double_t trackX[4], trackY[4]; };
std::vector<EventCand> gCandidates;
Double_t gSigma[4][nBarsMax] = {{0}};
Bool_t   gSigmaComputed = kFALSE;
std::map<TString, TH2F*> gHist;        // key = ChanLabel -> h2(track pos, half TW-corr time diff)
std::map<TString, FitRec> gResults;
std::map<TString, FitRec> gRefResults;
// The reference file's own PosSigma, loaded alongside velFit/cableFit --
// used for any channel whose vp/cable were themselves copied from this
// same reference (source "reference-copied"/"reference-copied-lowstat"):
// that paddle has too little of THIS run's own data to fit vp/cable, so
// it has just as little to compute a meaningful sigma residual from --
// use the reference's own sigma there instead of a near-empty StdDev().
Double_t gReferenceSigma[4][nBarsMax] = {{0}};
Bool_t   gHaveReferenceSigma = kFALSE;

Int_t     gIdx = 0;
TCanvas  *gCanvas = nullptr;
TControlBar *gControlBar = nullptr;
TF1      *gCurrentFit = nullptr;
TF1      *gRefFit = nullptr;
Int_t     gClickStage = 0;
Double_t  gPendingLo = 0.0;

// Sigma blocks carried over from an existing phodo_Vpcalib.param -- only
// used as a fallback if ComputeSigma() (see below) has nothing new.
TString gSigmaPosBlockRaw, gSigmaNegBlockRaw; // raw comma-separated rows, verbatim
Bool_t  gHaveCarriedSigma = kFALSE;

TString ResolvedRunPath(Int_t run) {
  // Original location -- switch back once the slurm-queued replay (with
  // the updated time-walk correction) has landed in the normal ROOTfiles
  // directory:
  return Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
  // return Form("/volatile/hallc/alphaE/ndelta_vcs2/selected_runs_before_hodo_calibration/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
}

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
// Section 5: build raw histograms (two-pass, single-hit-per-plane logic
//   from shms_hodo_calib/fitHodoCalib.C's velocity/cable section)
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

  // Book histograms: x = track position at the plane (cm), y = half the
  // TW-corrected pos-minus-neg TDC time diff (ns).
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    gHist[key] = new TH2F("h2_" + key,
        Form("TW-Corr. Time Diff vs. Track Pos %s Paddle %d;Hodoscope Track Position from Center (cm);Time-Walk Corrected TDC Time Difference (ns)",
             gPlaneNames[ch.plane].Data(), ch.paddle),
        400, -60, 60, 400, -15, 15);
    gHist[key]->SetDirectory(nullptr);
  }

  // PID branches
  Double_t pid1 = 0, pid2 = 0, pid3 = 0;
  T->SetBranchAddress(kPidBranch1, &pid1);
  T->SetBranchAddress(kPidBranch2, &pid2);
  T->SetBranchAddress(kPidBranch3, &pid3);

  Double_t dpVal = 0;
  T->SetBranchAddress("P.gtr.dp", &dpVal);

  // Per-plane, fixed-size (nBarsMax) arrays: hcana's "Good" TW-corrected
  // branches are already per-paddle indexed (paddle i -> slot i-1), with
  // an unhit paddle carrying a large sentinel that fails the <200ns cut
  // -- no Ndata-style hit counting needed here, unlike the raw ADC/TDC
  // branches in timewalk_calib_app.C.
  Double_t twPos[4][nBarsMax], twNeg[4][nBarsMax];
  Double_t trackX[4], trackY[4];
  for (Int_t ipl = 0; ipl < 4; ipl++) {
    TString base = "P.hod." + gPlaneNames[ipl];
    T->SetBranchAddress(base + ".GoodPosTdcTimeWalkCorr", twPos[ipl]);
    T->SetBranchAddress(base + ".GoodNegTdcTimeWalkCorr", twNeg[ipl]);
    T->SetBranchAddress(base + ".TrackXPos", &trackX[ipl]);
    T->SetBranchAddress(base + ".TrackYPos", &trackY[ipl]);
  }

  Long64_t nentries = T->GetEntries();

  // --- Pass 1: per-(plane,paddle) mean/stddev of the TW-average time ---
  std::map<TString, TH1F*> twAvg;
  for (auto &ch : gChannels)
    twAvg[ChanLabel(ch)] = new TH1F("twavg_" + ChanLabel(ch), "", 4000, -200, 200);

  printf("\n=== vpcalib pass 1/2 (mean/stddev): %lld events ===\n", nentries);
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    Bool_t pidOk = (!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low);
    if (!pidOk) continue;
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        TString key = ChanLabel(Channel{ipl, 0, ipad});
        if (!twAvg.count(key)) continue; // permanently-off, excluded
        Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
        if (p < kHardTimeCutHi && n < kHardTimeCutHi) twAvg[key]->Fill(0.5 * (p + n));
      }
    }
  }

  std::map<TString, Double_t> meanOf, stdOf;
  for (auto &kv : twAvg) { meanOf[kv.first] = kv.second->GetMean(); stdOf[kv.first] = kv.second->GetStdDev(); delete kv.second; }

  // --- Pass 2: single-hit-per-plane selection, fill the fit histograms ---
  printf("=== vpcalib pass 2/2 (single-hit selection + fill): %lld events ===\n", nentries);
  Long64_t nAfterDp = 0, nAfterPid = 0, nSingleHit4Plane = 0;
  std::map<TString, Long64_t> candidateCount; // times this paddle alone passed its own plane's window,
                                               // regardless of whether the OTHER 3 planes also had a single hit
  for (auto &ch : gChannels) candidateCount[ChanLabel(ch)] = 0;
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    ++nAfterDp;
    Bool_t pidOk = (!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low);
    if (!pidOk) continue;
    ++nAfterPid;

    Int_t goodCount[4] = {0, 0, 0, 0};
    Int_t goodPaddle[4] = {-1, -1, -1, -1};
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        TString key = ChanLabel(Channel{ipl, 0, ipad});
        if (!meanOf.count(key)) continue;
        Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
        if (!(p < kHardTimeCutHi && n < kHardTimeCutHi)) continue;
        Double_t avg = 0.5 * (p + n);
        if (avg > meanOf[key] - kNSig * stdOf[key] && avg < meanOf[key] + kNSig * stdOf[key]) {
          ++goodCount[ipl]; goodPaddle[ipl] = ipad;
        }
      }
    }
    // Record standalone candidacy for whichever paddle was its plane's
    // sole flagged hit, independent of the other 3 planes -- this is
    // what "empty despite good TW stats" needs to distinguish from.
    for (Int_t ipl = 0; ipl < 4; ipl++)
      if (goodCount[ipl] == 1) ++candidateCount[ChanLabel(Channel{ipl, 0, goodPaddle[ipl]})];

    Bool_t singleHit = (goodCount[0] == 1 && goodCount[1] == 1 && goodCount[2] == 1 && goodCount[3] == 1);
    if (!singleHit) continue;
    ++nSingleHit4Plane;

    // Same candidate used for the Vp/cable scatter fill below also goes
    // into the sigma cache -- by construction every plane here already
    // passed the <200ns hard cut (that's what made it "good" above), so
    // no extra filtering is needed later when sigma is computed from this.
    EventCand cand;
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      Int_t ipad = goodPaddle[ipl];
      cand.pad[ipl] = ipad;
      cand.twPos[ipl] = twPos[ipl][ipad - 1]; cand.twNeg[ipl] = twNeg[ipl][ipad - 1];
      cand.trackX[ipl] = trackX[ipl]; cand.trackY[ipl] = trackY[ipl];
    }
    gCandidates.push_back(cand);

    for (Int_t ipl = 0; ipl < 4; ipl++) {
      Int_t ipad = goodPaddle[ipl];
      TString key = ChanLabel(Channel{ipl, 0, ipad});
      if (!gHist.count(key)) continue;
      Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
      Double_t halfDiff = 0.5 * (n - p);
      Double_t pos = (ipl == 0 || ipl == 2) ? trackY[ipl] : trackX[ipl]; // 1x/2x -> TrackYPos, 1y/2y -> TrackXPos, per fitHodoCalib.C
      gHist[key]->Fill(pos, halfDiff);
    }
  }
  printf("  survival: %lld total -> %lld after |dp| cut -> %lld after PID -> %lld after single-hit-in-all-4-planes\n",
         nentries, nAfterDp, nAfterPid, nSingleHit4Plane);
  Int_t nEmptyWithCandidates = 0, nEmptyNoCandidates = 0;
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    if (gHist[key]->GetEntries() > 0) continue;
    if (candidateCount[key] > 0) {
      ++nEmptyWithCandidates;
      printf("  EMPTY but had %lld standalone candidate(s), cut by the OTHER 3 planes' single-hit requirement: %s\n", candidateCount[key], key.Data());
    } else {
      ++nEmptyNoCandidates;
      printf("  EMPTY with ZERO standalone candidates (own-plane time window or PID/dp cuts): %s\n", key.Data());
    }
  }
  if (nEmptyWithCandidates || nEmptyNoCandidates)
    printf("  %d empty channel(s) total: %d killed by cross-plane coincidence, %d killed upstream (own cuts)\n",
           nEmptyWithCandidates + nEmptyNoCandidates, nEmptyWithCandidates, nEmptyNoCandidates);
  printf("=== Done: %lld events processed ===\n\n", nentries);

  if (f) f->Close();
  if (chain) delete chain;
}

// ===========================================================================
// Section 6: fitting (free 2-param, or constrained slope-fixed 1-param)
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
    fit.SetParameters(1.0 / kVelSet[ch.plane], 0.0);
    TFitResultPtr res = h2->Fit(&fit, "SREQ0");
    status = res;
    r.fitted = (status == 0);
    Double_t slope = fit.GetParameter(0), slopeErr = fit.GetParError(0);
    r.vp = (slope != 0.0) ? 1.0 / slope : 0.0;
    r.vpErr = (slope != 0.0) ? TMath::Abs(slopeErr / (slope * slope)) : 0.0;
    r.cable = fit.GetParameter(1); r.cableErr = fit.GetParError(1);
    r.chi2ndf = (fit.GetNDF() > 0) ? fit.GetChisquare() / fit.GetNDF() : 0.0;
  } else {
    TF1 fit("vpFitConstrained", Form("(1./%f)*x + [0]", kVelSet[ch.plane]), lo, hi);
    fit.SetParameter(0, 0.0);
    TFitResultPtr res = h2->Fit(&fit, "SREQ0");
    status = res;
    r.fitted = (status == 0);
    r.vp = kVelSet[ch.plane]; r.vpErr = 0.0;
    r.cable = fit.GetParameter(0); r.cableErr = fit.GetParError(0);
    r.chi2ndf = (fit.GetNDF() > 0) ? fit.GetChisquare() / fit.GetNDF() : 0.0;
  }
  if (status != 0) printf("WARNING: fit did not converge cleanly for %s (status=%d)\n", key.Data(), status);
}

Bool_t VelocityInRange(const Channel &ch, Double_t vp) {
  return vp > (kVelSet[ch.plane] - kVelVar) && vp < (kVelSet[ch.plane] + kVelVar);
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

// Grabs a whole param block's raw text verbatim (from "= " through the
// last comma-row) -- used as a fallback carry-over for PosSigma/NegSigma
// if ComputeSigma() has nothing new to write for some reason.
TString ExtractRawBlock(const TString &path, const TString &key) {
  std::ifstream in(path.Data());
  if (!in.is_open()) return "";
  std::string line;
  bool found = false;
  while (std::getline(in, line)) {
    if (line.find(key.Data()) != std::string::npos && line.find('=') != std::string::npos) { found = true; break; }
  }
  if (!found) return "";
  TString block = line.c_str();
  block += "\n";
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
  if (staged) return Form("%s/SHMS/HODO/phodo_Vpcalib_%s.param", gParamDir.Data(), tag.Data());
  return Form("%s/SHMS/HODO/phodo_Vpcalib.param", gParamDir.Data());
}

TString EffectiveTag() { return gOutputTag.Length() > 0 ? gOutputTag : gRunTag; }

void LoadLegacyParam(const TString &path) {
  std::ifstream in(path.Data());
  if (!in.is_open()) { printf("  (no legacy param file at %s)\n", path.Data()); return; }
  Double_t vel[4][nBarsMax] = {{0}}, cable[4][nBarsMax] = {{0}}, sig[4][nBarsMax] = {{0}};
  ParseParamBlock(in, "phodo_velFit", vel);
  in.clear(); in.seekg(0); ParseParamBlock(in, "phodo_cableFit", cable);
  in.clear(); in.seekg(0); ParseParamBlock(in, "phodo_PosSigma", sig);
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) gReferenceSigma[ipl][ipad] = sig[ipl][ipad];
  gHaveReferenceSigma = (ExtractRawBlock(path, "phodo_PosSigma").Length() > 0);
  Int_t nSkippedPermanent = 0;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad) || IsPermanentlyOff(ipl, 1, ipad)) { ++nSkippedPermanent; continue; }
      Channel ch{ipl, 0, ipad};
      FitRec r;
      r.vp = vel[ipl][ipad - 1]; r.cable = cable[ipl][ipad - 1];
      r.lo = gFitRangeLow; r.hi = gFitRangeHigh;
      r.constrained = !VelocityInRange(ch, r.vp); // best guess at which mode produced this saved value
      r.fitted = kTRUE; r.source = "reference-file";
      gRefResults[ChanLabel(ch)] = r;
    }
  printf("  loaded reference param: %s", path.Data());
  if (nSkippedPermanent) printf(" (skipped %d permanently-off channel(s))", nSkippedPermanent);
  printf("\n");
}

// compareOnly mode: loads THIS run/tag's own already-staged calibration
// (not the reference/vanilla) directly into gResults, with no fitting at
// all. Returns false if no staged file exists for EffectiveTag() -- the
// caller aborts in that case rather than silently fitting or defaulting,
// since compareOnly asserts a calibration was already committed (via a
// prior interactive session or a non-compareOnly batch run).
Bool_t LoadOwnSavedResults() {
  TString path = LegacyParamPath(EffectiveTag(), kTRUE);
  std::ifstream test(path.Data());
  if (!test.is_open()) return kFALSE;
  test.close();
  Double_t vel[4][nBarsMax] = {{0}}, cable[4][nBarsMax] = {{0}};
  std::ifstream in(path.Data());
  ParseParamBlock(in, "phodo_velFit", vel);
  in.clear(); in.seekg(0); ParseParamBlock(in, "phodo_cableFit", cable);
  Int_t nSkippedPermanent = 0;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad) || IsPermanentlyOff(ipl, 1, ipad)) { ++nSkippedPermanent; continue; }
      Channel ch{ipl, 0, ipad};
      FitRec r;
      r.vp = vel[ipl][ipad - 1]; r.cable = cable[ipl][ipad - 1];
      r.lo = gFitRangeLow; r.hi = gFitRangeHigh;
      r.constrained = !VelocityInRange(ch, r.vp);
      r.fitted = kTRUE; r.source = "saved";
      gResults[ChanLabel(ch)] = r;
    }
  printf("  loaded this run's own saved calibration: %s", path.Data());
  if (nSkippedPermanent) printf(" (skipped %d permanently-off channel(s))", nSkippedPermanent);
  printf("\n");
  return kTRUE;
}

void LoadReference() {
  if (gReferenceTag.Length() == 0 && gReferenceRun == 0) return;
  // Explicit request for the vanilla (untagged) file, e.g. an existing
  // PARAM/SHMS/HODO/phodo_Vpcalib.param from a prior experiment/setup --
  // skip the staged-tag lookup entirely rather than relying on it
  // incidentally not being found.
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

// Try to carry over PosSigma/NegSigma from whatever Vpcalib.param already
// exists (vanilla, then this tag's own prior save takes priority if it
// has one) so the staged output file stays complete for hcana.
void LoadCarriedSigma() {
  TString ownPath = LegacyParamPath(EffectiveTag(), kTRUE);
  std::ifstream ownTest(ownPath.Data());
  TString sourcePath = ownTest.is_open() ? ownPath : LegacyParamPath("", kFALSE);
  ownTest.close();
  gSigmaPosBlockRaw = ExtractRawBlock(sourcePath, "phodo_PosSigma");
  gSigmaNegBlockRaw = ExtractRawBlock(sourcePath, "phodo_NegSigma");
  gHaveCarriedSigma = (gSigmaPosBlockRaw.Length() > 0 && gSigmaNegBlockRaw.Length() > 0);
  if (gHaveCarriedSigma) printf("Carried over PosSigma/NegSigma from %s\n", sourcePath.Data());
  else printf("No existing PosSigma/NegSigma found -- will write a 1.0 placeholder (a proper sigma pass is still needed).\n");
}

// True if this channel's current vp AND cable match the loaded reference
// (vanilla) values to within the precision the .param files are written
// with (6 decimals). Such a channel has the same timing calibration as
// the reference, so the reference's sigma is carried over instead of
// being recomputed. Covers explicit reference copies as well as fits or
// resumed values that end up identical to the reference.
const Double_t kRefMatchTol = 1.0e-6;
Bool_t SameAsReference(const TString &key) {
  if (!gHaveReferenceSigma) return kFALSE;
  if (!gResults.count(key) || !gRefResults.count(key)) return kFALSE;
  const FitRec &r = gResults[key];
  const FitRec &rr = gRefResults[key];
  if (!r.fitted) return kFALSE;
  return (TMath::Abs(r.vp - rr.vp) < kRefMatchTol && TMath::Abs(r.cable - rr.cable) < kRefMatchTol);
}

// Computes phodo_PosSigma/phodo_NegSigma from the SAME candidates cached
// during BuildHistos() (no second tree read), using whatever v_p/cable
// this session currently holds in gResults for each paddle -- fresh fit,
// resumed, or reference-copied, whatever's actually there at save time.
// phodo_PosSigma and phodo_NegSigma are the same array written twice in
// the original script; there's no separate pos/neg sigma quantity.
//   DiffDistTWCorr = v_p * 0.5*(TW_neg - 2*cable - TW_pos)   [cm]
//   residual       = DiffDistTWCorr - TrackPos               [cm]
//   sigma          = StdDev(residual) / (2 * v_p)            [ns]
// Skipped entirely under compareOnly (nothing gets saved anyway).
void ComputeSigma() {
  if (gCandidates.empty()) { printf("No cached candidates -- cannot compute sigma.\n"); return; }
  std::map<TString, TH1F*> resid;
  for (auto &ch : gChannels) resid[ChanLabel(ch)] = new TH1F("resid", "", 2000, -120, 80);

  for (auto &c : gCandidates) {
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      Int_t ipad = c.pad[ipl];
      TString key = ChanLabel(Channel{ipl, 0, ipad});
      if (!resid.count(key)) continue;
      if (!(gResults.count(key) && gResults[key].fitted)) continue; // need this paddle's v_p/cable to be known
      // A channel whose vp/cable were themselves copied from the
      // reference (not fit from this run's own data) has just as little
      // of this run's data to compute a meaningful residual from -- skip
      // filling it here; the final loop below uses the reference's own
      // sigma for these instead of a near-empty StdDev().
      if (gResults[key].source.BeginsWith("reference-copied") || SameAsReference(key)) continue;
      Double_t vp = gResults[key].vp, cable = gResults[key].cable;
      if (vp == 0.0) continue;
      Double_t diffDist = vp * 0.5 * (c.twNeg[ipl] - 2.0 * cable - c.twPos[ipl]);
      Double_t trackPos = (ipl == 0 || ipl == 2) ? c.trackY[ipl] : c.trackX[ipl];
      resid[key]->Fill(diffDist - trackPos);
    }
  }

  Int_t nComputed = 0, nCopiedFromRef = 0;
  for (auto &ch : gChannels) {
    TString key = ChanLabel(ch);
    TH1F *h = resid[key];
    if (gResults.count(key) && (gResults[key].source.BeginsWith("reference-copied") || SameAsReference(key))) {
      if (gHaveReferenceSigma) {
        gSigma[ch.plane][ch.paddle - 1] = gReferenceSigma[ch.plane][ch.paddle - 1];
        ++nCopiedFromRef;
      }
      // else: no reference sigma available either -- leaves gSigma at 0,
      // same "not yet determined" convention as an unfit channel.
    } else if (h->GetEntries() > 0 && gResults.count(key) && gResults[key].fitted && gResults[key].vp != 0.0) {
      gSigma[ch.plane][ch.paddle - 1] = h->GetStdDev() / (2.0 * gResults[key].vp);
      ++nComputed;
    }
    delete h;
  }
  gSigmaComputed = (nComputed > 0 || nCopiedFromRef > 0);
  printf("Computed sigma for %d paddle(s) from this run's own data, copied reference sigma for %d paddle(s) (vp/cable itself was reference-copied), out of %d total.\n",
         nComputed, nCopiedFromRef, (int)gChannels.size());
}

void WriteLegacyParam(const TString &tag) {
  TString outPath = LegacyParamPath(tag, kTRUE);
  gSystem->mkdir(gSystem->DirName(outPath), kTRUE);
  std::ofstream out(outPath.Data());
  out << "; SHMS Hodoscope Parameter File Containing propagation velocities per paddle " << std::endl;
  out << "; and signal cable time diff. offsets per paddle " << std::endl;
  out << Form("; tag %s ", tag.Data()) << std::endl << std::endl;

  auto writeBlock = [&](const TString &name, Bool_t isVp) {
    out << ";" << name << std::endl;
    out << name << " = ";
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++) {
      for (Int_t ipl = 0; ipl < 4; ipl++) {
        Double_t v = kVelSet[ipl]; // default fallback for vp; overwritten below for cable
        if (!isVp) v = 0.0;
        if (ipad < gNbars[ipl]) {
          Channel ch{ipl, 0, ipad + 1};
          TString key = ChanLabel(ch);
          if (gResults.count(key) && gResults[key].fitted) v = isVp ? gResults[key].vp : gResults[key].cable;
          // else: permanently-off or never-fitted paddle keeps the
          // fallback (kVelSet[plane] for vp, 0.0 for cable), matching
          // fitHodoCalib.C's own "fit failed twice" convention.
        }
        out << std::fixed << v;
        if (ipl != 3) out << ", ";
      }
      out << "," << std::endl;
    }
    out << std::endl;
  };
  writeBlock("phodo_velFit", kTRUE);
  writeBlock("phodo_cableFit", kFALSE);

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
    // PosSigma and NegSigma are the same array written twice, matching
    // the original script -- there's no separate pos/neg sigma quantity.
    writeSigmaBlock("phodo_PosSigma");
    writeSigmaBlock("phodo_NegSigma");
  } else if (gHaveCarriedSigma) {
    out << ";PMTs Time Diff. Sigma Parameters (carried over, not recomputed this session)" << std::endl;
    out << gSigmaPosBlockRaw.Data() << std::endl;
    out << gSigmaNegBlockRaw.Data() << std::endl;
  } else {
    out << ";PMTs Time Diff. Sigma Parameters (placeholder)" << std::endl;
    out << "phodo_PosSigma = ";
    for (Int_t i = 0; i < nBarsMax; i++) out << "1.000000, 1.000000, 1.000000, 1.000000," << std::endl;
    out << std::endl;
    out << "phodo_NegSigma = ";
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
  TString path = Form("%s/vpcalib_shms_%s.json", gOutDir.Data(), tag.Data());
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

  FitRec &r = gResults.count(key) ? gResults[key] : (gResults[key] = FitRec{gFitRangeLow, gFitRangeHigh});

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

  // Skip the "new fit" curve when it's really just a verbatim copy of
  // the reference (compareOnly's no-fitting-at-all mode, or the
  // low-stat fallback) -- it would sit exactly on top of the black
  // dotted reference curve, drawing it twice adds nothing.
  if (r.fitted && !r.source.BeginsWith("reference-copied")) {
    if (gCurrentFit) delete gCurrentFit;
    gCurrentFit = new TF1("curFit", "[0]*x + [1]", xLoFull, xHiFull); // extrapolated display, same as timewalk_calib_app.C
    gCurrentFit->SetNpx(1000);
    Double_t slope = (r.vp != 0.0) ? 1.0 / r.vp : 0.0;
    gCurrentFit->SetParameters(slope, r.cable);
    gCurrentFit->SetLineColor(r.constrained ? kMagenta + 1 : kRed);
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
      TText *warn = pt.AddText(Form("OUT OF RANGE [%.2f, %.2f] -- try Constrained", kVelSet[ch.plane]-kVelVar, kVelSet[ch.plane]+kVelVar));
      warn->SetTextColor(kRed);
    }
  } else pt.AddText("(not fitted)");
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
  gResults[key].lo = gFitRangeLow; gResults[key].hi = gFitRangeHigh;
  gClickStage = 0;
  DoFit(ch, gFitRangeLow, gFitRangeHigh, kFALSE, "manual");
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

// Directly ADOPTS the reference's own vp/cable numbers, with no local
// fit at all -- unlike UseReference() (which re-fits THIS paddle's own
// histogram, just seeded with the reference's range/mode). Use this for
// a low-statistics paddle where your own data can't be trusted to
// produce a reliable fit, but an existing calibration's value can.
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

void PauseSave() {
  if (gCompareOnly) { printf("compareOnly is on -- saving is disabled for this session.\n"); return; }
  ComputeSigma();
  SaveJSON(EffectiveTag()); WriteLegacyParam(EffectiveTag()); printf("[progress saved]\n");
}

void SaveFinish() {
  if (gCompareOnly) { printf("compareOnly is on -- saving is disabled for this session.\n"); return; }
  PauseSave();
  printf("[saved and finished -- rename %s/SHMS/HODO/phodo_Vpcalib_%s.param to phodo_Vpcalib.param for hcana]\n",
         gParamDir.Data(), EffectiveTag().Data());
}

void OnClick() {
  Int_t event = gPad->GetEvent();
  if (event != 11) return;
  Int_t px = gPad->GetEventX();
  Double_t x = gPad->AbsPixeltoX(px);
  Channel &ch = gChannels[gIdx];
  TString key = ChanLabel(ch);
  FitRec &r = gResults.count(key) ? gResults[key] : (gResults[key] = FitRec{gFitRangeLow, gFitRangeHigh});
  if (gClickStage == 0) {
    gPendingLo = x;
    gClickStage = 1;
    printf("  lo = %.1f cm (click again for hi)\n", x);
  } else {
    Double_t lo = TMath::Min(gPendingLo, x);
    Double_t hi = TMath::Max(gPendingLo, x);
    gClickStage = 0;
    DoFit(ch, lo, hi, r.constrained, "manual"); // keeps whatever mode was active
    DrawChannel();
  }
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "SHMS Vp/Cable Fit Controls", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Previous channel");
  gControlBar->AddButton("Next >>", "GoNext();", "Next channel");
  gControlBar->AddButton("Refit", "Refit();", "Re-run the current fit mode with the current range");
  gControlBar->AddButton("Toggle Fit Mode", "ToggleFitMode();", "Switch FREE <-> CONSTRAINED (fixed v_p)");
  gControlBar->AddButton("Reset Range", "ResetRange();", "Back to [-40,40] cm, FREE mode");
  gControlBar->AddButton("Use Reference", "UseReference();", "Re-fit here, seeded with reference's range+mode");
  gControlBar->AddButton("Copy Reference Value", "CopyReferenceValue();", "Adopt reference's vp/cable directly, no local fit (for low-stat paddles)");
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
      // Low-statistics fallback: copy the reference's vp/cable directly
      // rather than fitting noise, same as the interactive
      // "Copy Reference Value" button.
      if (gMinEntriesForFit > 0.0 && gRefResults.count(key) && h2 && h2->GetEntries() < gMinEntriesForFit) {
        FitRec &rr = gRefResults[key];
        FitRec &r = gResults[key];
        r.lo = rr.lo; r.hi = rr.hi; r.constrained = rr.constrained;
        r.vp = rr.vp; r.vpErr = 0.0; r.cable = rr.cable; r.cableErr = 0.0;
        r.chi2ndf = 0.0; r.fitted = kTRUE; r.source = "reference-copied-lowstat";
        continue;
      }
      Double_t lo = gFitRangeLow, hi = gFitRangeHigh;
      Bool_t constrained = kFALSE;
      TString src = "default";
      if (gRefResults.count(key)) { lo = gRefResults[key].lo; hi = gRefResults[key].hi; constrained = gRefResults[key].constrained; src = "reference-auto"; }
      DoFit(ch, lo, hi, constrained, src);
      // Auto fallback to constrained mode if the free fit's v_p is out of range,
      // mirroring fitHodoCalib.C's own automatic behavior in batch mode.
      if (!constrained && gResults[key].fitted && !VelocityInRange(ch, gResults[key].vp)) {
        DoFit(ch, lo, hi, kTRUE, src + "+auto-constrained");
      }
    }
  }
  // else: gCompareOnly or gSigmaOnly -- gResults was already fully
  // populated by LoadOwnSavedResults() before this function was ever
  // called; no fitting happens here at all in either mode.
  gSystem->mkdir(gOutDir, kTRUE);
  TString pdfPath = Form("%s/vpcalib_shms_%s_summary.pdf", gOutDir.Data(), EffectiveTag().Data());
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
                  Double_t fitRangeLow = -40.0, Double_t fitRangeHigh = 40.0,
                  Bool_t applyCalCut = kTRUE, Double_t minEntriesForFit = 0.0,
                  Bool_t applyCerCut = kTRUE, Bool_t applyTrackCut = kTRUE,
                  Bool_t compareOnly = kFALSE, Bool_t sigmaOnly = kFALSE) {
  if (run == 0 && runs.Length() == 0) {
    printf("ERROR: must supply a run number, e.g. vpcalib_app(26107)\n"); return;
  }
  gReferenceRun = referenceRun;
  gReferenceTag = referenceTag;
  gOutputTag = outputTag;
  gOutDir = outDir; gParamDir = paramDir;
  gFitRangeLow = fitRangeLow; gFitRangeHigh = fitRangeHigh;
  gDpCut = dpCut;
  gMinEntriesForFit = minEntriesForFit;
  gApplyCalCut = applyCalCut;
  gApplyCerCut = applyCerCut;
  gApplyTrackCut = applyTrackCut;
  gCompareOnly = compareOnly;
  if (compareOnly && referenceTag.Length() == 0 && referenceRun == 0) {
    gReferenceTag = "vanilla";
    printf("compareOnly=true with no referenceRun/referenceTag given -- defaulting to referenceTag=\"vanilla\" "
           "(the currently-saved ../../PARAM/SHMS/HODO/phodo_Vpcalib.param) so there's something to compare against.\n");
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

  gSigmaOnly = sigmaOnly;
  if (compareOnly && sigmaOnly)
    printf("WARNING: both compareOnly and sigmaOnly are set -- compareOnly wins, nothing will be saved.\n");

  if (compareOnly || sigmaOnly) {
    // Assert the calibration was already committed (interactive session,
    // or a prior batch run) -- neither mode fits from scratch; both load
    // this tag's existing velFit/cableFit instead. compareOnly never
    // saves; sigmaOnly DOES save (recomputed sigma + velFit/cableFit as
    // loaded, unchanged).
    TString ownPath = LegacyParamPath(EffectiveTag(), kTRUE);
    std::ifstream ownTest(ownPath.Data());
    if (!ownTest.is_open()) {
      printf("ERROR: compareOnly/sigmaOnly requires an existing staged calibration for tag \"%s\", but none found at %s.\n"
             "Run vpcalib_app in a normal (fitting) session first to produce it.\n",
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
  gCanvas = new TCanvas("vpcalibCanvas", "SHMS Vp/Cable Calibration", 900, 700);
  DrawChannel();
  gCanvas->AddExec("dynamic", "OnClick()");
  gCanvas->Update();
  MakeControlBar();

  printf("\n=== run(s) %s (output tag: %s) -- %d channel(s) loaded (SHMS) ===\n",
         gRunTag.Data(), EffectiveTag().Data(), (int)gChannels.size());
  printf("Click twice on the histogram: 1st = fit lo, 2nd = fit hi (auto-refits, keeps current mode).\n");
  printf("Use the 'SHMS Vp/Cable Fit Controls' panel to navigate, toggle fit mode, and save.\n");
}
