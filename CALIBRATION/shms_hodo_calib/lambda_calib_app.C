// lambda_calib_app.C (SHMS)
//
// Fits the SHMS Hodoscope's per-paddle inter-plane timing offset
// (phodo_LCoeff) via the "Hodoscope Planes Time Difference Corrections"
// step of hodo_calib.pdf (Eq. 3-8), replicating
// shms_hodo_calib/fitHodoCalib.C's final section.
//
// UNLIKE the time-walk and velocity/cable apps, this is NOT 61
// independent per-paddle fits -- it's a SINGLE global least-squares
// solve across all 61 paddles simultaneously:
//
//   For each clean single-track event, and each of the 6 hodoscope
//   plane-pairs (1x-1y, 1x-2x, 1x-2y, 1y-2x, 1y-2y, 2x-2y), the paddle
//   pair (i,j) that fired contributes one equation:
//     lambda_i - lambda_j = D_ij/c - (t_i - t_j) == b_ij
//   where D_ij is the 3D distance between the two planes' track
//   crossing points, c = 30 cm/ns (speed of light -- NOT your fitted
//   v_p, which plays no role in this stage), and t_i is the average
//   (cable-corrected) TW-corrected TDC time at plane i.
//
// All such equations across the whole dataset are accumulated into the
// normal-equations system Ay*lambda = bVec (Ay is 61x61, symmetric) and
// solved via TDecompSVD -- SVD is used because Ay can be poorly
// conditioned (one paddle, 1x-paddle-7, is fixed as the reference with
// lambda=0, an explicit degeneracy; some paddle-pairs may also simply
// have few contributing events).
//
// Because there's no natural "fit range" to click per channel, this app
// is mostly BATCH: build histograms/accumulate Ay,bVec, solve, done.
// The "interactive" part is a lightweight per-plane review of the
// resulting lambda values (bar chart + table), with the one genuinely
// tunable knob (the hard TDC time cut used to admit a paddle into the
// system, default 125 ns per the original script) adjustable and
// re-solvable without re-reading the tree.
//
// REQUIRES an existing phodo_Vpcalib_<vpcalibTag>.param (from
// vpcalib_app.C) for phodo_cableFit -- this stage does not compute cable
// offsets itself. phodo_velFit and the sigma blocks are carried over
// verbatim from that same file into the merged output (matching
// fitHodoCalib.C's convention of one combined param file with all four
// blocks: velFit, cableFit, sigma, LCoeff).
//
// Usage:
//   root -l 'lambda_calib_app.C(26107)'                        // interactive review
//   root -l -b -q 'lambda_calib_app.C(26107, "", "", true)'    // batch: solve + write, no review
//
// Comparing a fresh solve against an existing phodo_LCoeff -- referenceTag
// loads it (as with the other apps, "vanilla" means the untagged
// PARAM/SHMS/HODO/phodo_Vpcalib.param) and overlays it on DrawPlane() as
// black open circles alongside the new solve's blue filled squares.
// compareOnly solves and draws normally but never writes JSON/the merged
// param file, so you can sanity-check a new solve without touching what's
// currently saved:
//   root -l -b -q 'lambda_calib_app.C(0, "26483,26484,26485,26486,26487,26488", "", true, "", 10.0, 125.0, "", "./lambda_qa", "../../PARAM", "vanilla", true, true, true, true)'
//
// Controls (interactive):
//   << Prev / Next >>       -- step through the 4 planes' bar charts
//   Redo with new cut...     -- console prompt for a new hard TDC time
//                              cut (ns), re-solves without re-reading the
//                              tree (the event-level selection itself,
//                              i.e. which paddle is "good" per plane per
//                              event, does NOT depend on this cut -- see
//                              note in BuildAndSolve)
//   Pause (save) / Save && Finish -- write the merged 4-block param file
//
// Outputs (under outDir, default "./lambda_qa"):
//   lambda_shms_<tag>.json          -- per-paddle LCoeff + provenance
//   lambda_shms_<tag>_summary.pdf   -- one page per plane, bar chart
// Outputs (under paramDir, default "../../PARAM"):
//   SHMS/HODO/phodo_Vpcalib_<tag>.param -- MERGED: velFit/cableFit/sigma
//   carried over from vpcalibTag's file + the new phodo_LCoeff block.
//
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TControlBar.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
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
static const Int_t   gNbars[4]      = {13, 13, 14, 21};
static const Int_t   nBarsMax       = 21;
static const Int_t   kRefPad[4]     = {0, 13, 26, 40}; // cumulative global-index offset per plane
static const Int_t   kNpar          = 61;              // total paddles = sum(gNbars)
static const Int_t   kRefPlane = 0, kRefPaddle = 7;     // 1x paddle 7 fixed, lambda = 0
static const Int_t   kRefGlobalIdx = kRefPad[kRefPlane] + kRefPaddle - 1; // = 6

static const Double_t kPlaneZ0[4] = {52.1, 61.7, 271.4, 282.4}; // cm, plane central distance from focal plane
static const Double_t kDz[4] = {2.12, -2.12, -2.12, -5.4}; // cm, per-paddle-parity z offset
static const Double_t kLightSpeed = 30.0; // cm/ns, D_ij/c term -- NOT your fitted v_p

struct OffPmt { Int_t plane; Int_t side; Int_t paddle; };
static const std::vector<OffPmt> gPermanentlyOff = {
  {3, 0, 1}, {3, 0, 2}, {3, 0, 5}, {3, 0, 19}, {3, 0, 20}, {3, 0, 21},
  {3, 1, 1}, {3, 1, 2},            {3, 1, 19}, {3, 1, 20}, {3, 1, 21},
};
Bool_t IsPermanentlyOff(Int_t plane, Int_t side, Int_t paddle) {
  for (auto &o : gPermanentlyOff)
    if (o.plane == plane && o.side == side && o.paddle == paddle) return kTRUE;
  return kFALSE;
}

static const TString kPidBranch1 = "P.cal.etracknorm"; static const Double_t kPidCut1Low = 0.7;
static const TString kPidBranch2 = "P.hgcer.npeSum";   static const Double_t kPidCut2Low = 0.5;
static const TString kPidBranch3 = "P.dc.ntrack";      static const Double_t kPidCut3Low = 0.0;
static const Double_t kMeanTimeCutHi = 200.0; // same as vpcalib_app.C's pre-selection (mean/stddev + single-hit-per-plane)
static const Double_t kNSig = 1.0;
static const Double_t kTrackPosSanity = 200.0; // TrackXPos/TrackYPos < this, else treated as kBig/invalid

TString ResolvedRunPath(Int_t run) {
  // See note in vpcalib_app.C -- switch back once the normal replay dir is current:
  return Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
  // return Form("/volatile/hallc/alphaE/ndelta_vcs2/selected_runs_before_hodo_calibration/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
}

// ===========================================================================
// Section 2: globals
// ===========================================================================

Int_t    gRun = 0;
TString  gRootFile, gOutDir = "./lambda_qa", gParamDir = "../../PARAM";
std::vector<Int_t> gRunList;
TString  gRunTag, gOutputTag, gVpcalibTag; // gVpcalibTag: which phodo_Vpcalib_<tag>.param to read cableFit from (defaults to EffectiveTag())
Double_t gDpCut = 1.0e9;
Double_t gHardTimeCut = 125.0; // the one tunable knob -- redoable without re-reading the tree
Bool_t   gApplyCalCut = kTRUE;   // P.cal.etracknorm cut -- turn off for a proton-arm run
Bool_t   gApplyCerCut = kTRUE;   // P.hgcer.npeSum cut
Bool_t   gApplyTrackCut = kTRUE; // P.dc.ntrack cut

// Carried-over blocks from the input Vpcalib file (verbatim text, written
// back unchanged into the merged output).
TString gVelBlockRaw, gCableBlockRaw, gPosSigmaBlockRaw, gNegSigmaBlockRaw;
Double_t gCableFit[4][nBarsMax] = {{0}}; // parsed, needed for the actual math

// Per-event candidate cache built ONCE from the tree (single-hit-per-plane
// selection + raw TW-corrected times + track position), so "Redo with new
// cut" only re-does the cheap accumulation step, not the tree read.
struct EventCand {
  Int_t plane[4];       // always {0,1,2,3} -- kept for clarity
  Int_t pad[4];          // the single good paddle per plane (-1 if none)
  Double_t twPos[4], twNeg[4]; // both sides' TW-corrected times for that paddle
  Double_t trackX[4], trackY[4];
};
std::vector<EventCand> gCandidates; // one entry per event with all 4 planes single-hit

TVectorD gLambdaSolved(kNpar);
Bool_t   gSolved = kFALSE;
Long64_t gNGoodEvents = 0;
Bool_t   gDecomposeOk = kFALSE, gSolveOk = kFALSE;

// Optional reference: an existing phodo_LCoeff (from a prior solve, or
// the vanilla PARAM/SHMS/HODO/phodo_Vpcalib.param) loaded for visual
// comparison against a fresh solve. Nothing to do with fitting -- purely
// a second series drawn on DrawPlane().
Double_t gLambdaReference[kNpar] = {0};
Bool_t   gHaveReferenceLambda = kFALSE;
TString  gReferenceTag; // "" = no reference; "vanilla" = the untagged file; else a staged tag
Bool_t   gCompareOnly = kFALSE; // batch mode: write the summary PDF (new solve vs. reference) but never touch JSON/.param output

Int_t gPlaneIdx = 0; // which of the 4 planes the review UI is showing
TCanvas *gCanvas = nullptr;
TControlBar *gControlBar = nullptr;

TString EffectiveTag() { return gOutputTag.Length() > 0 ? gOutputTag : gRunTag; }

// ===========================================================================
// Section 3: legacy param I/O -- reading cableFit (+ carrying other blocks)
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

TString VpcalibPath(const TString &tag, Bool_t staged) {
  if (staged) return Form("%s/SHMS/HODO/phodo_Vpcalib_%s.param", gParamDir.Data(), tag.Data());
  return Form("%s/SHMS/HODO/phodo_Vpcalib.param", gParamDir.Data());
}

Bool_t LoadInputVpcalib() {
  TString tag = gVpcalibTag.Length() > 0 ? gVpcalibTag : EffectiveTag();
  TString path = VpcalibPath(tag, kTRUE);
  std::ifstream test(path.Data());
  if (!test.is_open()) {
    test.close();
    path = VpcalibPath("", kFALSE);
    std::ifstream test2(path.Data());
    if (!test2.is_open()) {
      printf("ERROR: no phodo_Vpcalib param file found (tried tag '%s' and the vanilla file). Run vpcalib_app.C first.\n", tag.Data());
      return kFALSE;
    }
    test2.close();
  } else test.close();

  Double_t cable[4][nBarsMax] = {{0}};
  std::ifstream in(path.Data());
  ParseParamBlock(in, "phodo_cableFit", cable);
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++)
      gCableFit[ipl][ipad] = cable[ipl][ipad];

  gVelBlockRaw      = ExtractRawBlock(path, "phodo_velFit");
  gCableBlockRaw    = ExtractRawBlock(path, "phodo_cableFit");
  gPosSigmaBlockRaw = ExtractRawBlock(path, "phodo_PosSigma");
  gNegSigmaBlockRaw = ExtractRawBlock(path, "phodo_NegSigma");
  printf("Loaded cableFit (and carried velFit/sigma) from %s\n", path.Data());
  return kTRUE;
}

// Loads an existing phodo_LCoeff block for visual comparison only -- has
// no effect on the fresh solve. gReferenceTag == "vanilla" loads the
// untagged file; otherwise it's treated as a staged tag (falls back to
// vanilla if that staged file doesn't exist, same convention as the
// other apps' reference loading).
void LoadReferenceLambda() {
  if (gReferenceTag.Length() == 0) return;
  TString path;
  if (gReferenceTag == "vanilla") {
    path = VpcalibPath("", kFALSE);
  } else {
    path = VpcalibPath(gReferenceTag, kTRUE);
    std::ifstream test(path.Data());
    if (!test.is_open()) path = VpcalibPath("", kFALSE);
    else test.close();
  }
  std::ifstream in(path.Data());
  if (!in.is_open()) { printf("  (no reference LCoeff file at %s)\n", path.Data()); return; }
  Double_t lcoeff[4][nBarsMax] = {{0}};
  ParseParamBlock(in, "phodo_LCoeff", lcoeff);
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < gNbars[ipl]; ipad++) {
      Int_t gidx = kRefPad[ipl] + ipad;
      gLambdaReference[gidx] = lcoeff[ipl][ipad];
    }
  gHaveReferenceLambda = kTRUE;
  printf("  loaded reference LCoeff from %s\n", path.Data());
}

void WriteMergedParam(const TString &tag) {
  TString outPath = VpcalibPath(tag, kTRUE);
  gSystem->mkdir(gSystem->DirName(outPath), kTRUE);
  std::ofstream out(outPath.Data());
  out << "; SHMS Hodoscope Parameter File Containing propagation velocities, cable offsets, " << std::endl;
  out << "; sigma resolution, and inter-plane timing (LCoeff) parameters per paddle " << std::endl;
  out << Form("; tag %s ", tag.Data()) << std::endl << std::endl;

  out << ";phodo_velFit" << std::endl;
  out << (gVelBlockRaw.Length() ? gVelBlockRaw.Data() : "phodo_velFit = ; (none carried over)\n") << std::endl;
  out << ";phodo_cableFit" << std::endl;
  out << (gCableBlockRaw.Length() ? gCableBlockRaw.Data() : "phodo_cableFit = ; (none carried over)\n") << std::endl;
  out << ";phodo_PosSigma" << std::endl;
  out << (gPosSigmaBlockRaw.Length() ? gPosSigmaBlockRaw.Data() : "phodo_PosSigma = ; (none carried over)\n") << std::endl;
  out << ";phodo_NegSigma" << std::endl;
  out << (gNegSigmaBlockRaw.Length() ? gNegSigmaBlockRaw.Data() : "phodo_NegSigma = ; (none carried over)\n") << std::endl;

  out << ";Timing Corrections Per Paddle, where 1X Paddle 7 has been set as the reference paddle" << std::endl;
  out << "phodo_LCoeff = ";
  for (Int_t ipad = 0; ipad < nBarsMax; ipad++) {
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      Double_t v = 0.0;
      if (ipad < gNbars[ipl] && gSolved && !IsPermanentlyOff(ipl, 0, ipad + 1)) {
        Int_t gidx = kRefPad[ipl] + ipad;
        v = (gidx == kRefGlobalIdx) ? 0.0 : gLambdaSolved[gidx];
      }
      out << std::fixed << v;
      if (ipl != 3) out << ", ";
    }
    out << "," << std::endl;
  }
  out.close();
  printf("Wrote %s\n", outPath.Data());
}

// ===========================================================================
// Section 4: build the candidate list (single tree pass -- reused across
//   any number of "Redo with new cut" re-solves)
// ===========================================================================

void BuildCandidates() {
  gCandidates.clear();
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

  Double_t pid1 = 0, pid2 = 0, pid3 = 0, dpVal = 0;
  T->SetBranchAddress(kPidBranch1, &pid1);
  T->SetBranchAddress(kPidBranch2, &pid2);
  T->SetBranchAddress(kPidBranch3, &pid3);
  T->SetBranchAddress("P.gtr.dp", &dpVal);

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

  // Pass 1: per-(plane,paddle) mean/stddev of the TW-average time (same
  // as vpcalib_app.C's pre-selection).
  std::map<TString, TH1F*> twAvg;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad) || IsPermanentlyOff(ipl, 1, ipad)) continue;
      twAvg[Form("%s[%d]", gPlaneNames[ipl].Data(), ipad)] = new TH1F("twavg", "", 4000, -200, 200);
    }

  printf("\n=== lambda pass 1/2 (mean/stddev): %lld events ===\n", nentries);
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    if (!((!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low))) continue;
    for (Int_t ipl = 0; ipl < 4; ipl++)
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        TString key = Form("%s[%d]", gPlaneNames[ipl].Data(), ipad);
        if (!twAvg.count(key)) continue;
        Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
        if (p < kMeanTimeCutHi && n < kMeanTimeCutHi) twAvg[key]->Fill(0.5 * (p + n));
      }
  }
  std::map<TString, Double_t> meanOf, stdOf;
  for (auto &kv : twAvg) { meanOf[kv.first] = kv.second->GetMean(); stdOf[kv.first] = kv.second->GetStdDev(); delete kv.second; }

  // Pass 2: single-hit-per-plane selection; cache one EventCand per event
  // that passes, with BOTH sides' raw TW-corrected times so "Redo with a
  // new hard cut" can re-threshold without re-reading the tree. The hard
  // cut itself is applied later, in Solve() -- it doesn't affect WHICH
  // paddle is picked as the plane's single candidate (that's governed by
  // the +-1sigma window above), only whether that candidate is admitted
  // into the matrix system.
  printf("=== lambda pass 2/2 (single-hit selection + cache): %lld events ===\n", nentries);
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    if (!((!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low))) continue;

    Int_t goodCount[4] = {0, 0, 0, 0};
    Int_t goodPaddle[4] = {-1, -1, -1, -1};
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        TString key = Form("%s[%d]", gPlaneNames[ipl].Data(), ipad);
        if (!meanOf.count(key)) continue;
        Double_t p = twPos[ipl][ipad - 1], n = twNeg[ipl][ipad - 1];
        if (!(p < kMeanTimeCutHi && n < kMeanTimeCutHi)) continue;
        Double_t avg = 0.5 * (p + n);
        if (avg > meanOf[key] - kNSig * stdOf[key] && avg < meanOf[key] + kNSig * stdOf[key]) {
          ++goodCount[ipl]; goodPaddle[ipl] = ipad;
        }
      }
    }
    if (!(goodCount[0] == 1 && goodCount[1] == 1 && goodCount[2] == 1 && goodCount[3] == 1)) continue;

    Bool_t hodTrk = kTRUE;
    for (Int_t ipl = 0; ipl < 4; ipl++)
      if (!(trackX[ipl] < kTrackPosSanity && trackY[ipl] < kTrackPosSanity)) hodTrk = kFALSE;
    if (!hodTrk) continue;

    EventCand c;
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      c.plane[ipl] = ipl;
      c.pad[ipl] = goodPaddle[ipl];
      c.twPos[ipl] = twPos[ipl][goodPaddle[ipl] - 1];
      c.twNeg[ipl] = twNeg[ipl][goodPaddle[ipl] - 1];
      c.trackX[ipl] = trackX[ipl]; c.trackY[ipl] = trackY[ipl];
    }
    gCandidates.push_back(c);
  }
  printf("=== Done: %lld events processed, %zu single-hit-all-4-planes candidates cached ===\n\n", nentries, gCandidates.size());

  if (f) f->Close();
  if (chain) delete chain;
}

// ===========================================================================
// Section 5: accumulate normal equations from the cached candidates and
//   solve via SVD -- this is what "Redo with new cut" re-runs, cheaply.
// ===========================================================================

Double_t ZOf(Int_t plane, Int_t paddle) { return kPlaneZ0[plane] + ((paddle - 1) % 2) * kDz[plane]; }

void Solve() {
  TMatrixD Ay(kNpar, kNpar);
  TVectorD bVec(kNpar);
  Ay.Zero(); bVec.Zero();

  Long64_t nGood = 0;
  for (auto &c : gCandidates) {
    // Tighter hard cut, admitted into the matrix system only if BOTH
    // sides of ALL 4 planes' selected paddle pass it (matches
    // fitHodoCalib.C's "HARD CUT NOTICE" for this section).
    Bool_t allGood = kTRUE;
    for (Int_t ipl = 0; ipl < 4; ipl++)
      if (!(c.twPos[ipl] < gHardTimeCut && c.twNeg[ipl] < gHardTimeCut)) allGood = kFALSE;
    if (!allGood) continue;
    ++nGood;

    Int_t gidx[4]; Double_t t[4], x[4], y[4], z[4];
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      gidx[ipl] = kRefPad[ipl] + c.pad[ipl] - 1;
      Double_t negCorr = c.twNeg[ipl] - 2.0 * gCableFit[ipl][c.pad[ipl] - 1];
      t[ipl] = 0.5 * (c.twPos[ipl] + negCorr);
      x[ipl] = c.trackX[ipl]; y[ipl] = c.trackY[ipl]; z[ipl] = ZOf(ipl, c.pad[ipl]);
    }

    // The 6 plane-pairs: (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
    static const Int_t kPairI[6] = {0, 0, 0, 1, 1, 2};
    static const Int_t kPairJ[6] = {1, 2, 3, 2, 3, 3};
    for (Int_t r = 0; r < 6; r++) {
      Int_t pi = kPairI[r], pj = kPairJ[r];
      Double_t D = -TMath::Sqrt(TMath::Power(x[pi]-x[pj],2) + TMath::Power(y[pi]-y[pj],2) + TMath::Power(z[pi]-z[pj],2));
      Double_t b = D / kLightSpeed - (t[pi] - t[pj]);
      Int_t a = gidx[pi], bIdx = gidx[pj];
      Double_t ca = (a == kRefGlobalIdx) ? 0.0 : 1.0;   // reference paddle: coefficient forced to 0 (lambda_ref fixed at 0)
      Double_t cb = (bIdx == kRefGlobalIdx) ? 0.0 : -1.0;
      // Rank-1 update: Ay += [ca,cb]-row outer product; bVec += [ca,cb]*b
      if (ca != 0.0) { Ay[a][a] += ca*ca; bVec[a] += ca*b; }
      if (cb != 0.0) { Ay[bIdx][bIdx] += cb*cb; bVec[bIdx] += cb*b; }
      if (ca != 0.0 && cb != 0.0) { Ay[a][bIdx] += ca*cb; Ay[bIdx][a] += cb*ca; }
    }
  }
  gNGoodEvents = nGood;
  printf("Solve(): %lld / %zu cached candidates pass the %.1f ns hard cut\n", nGood, gCandidates.size(), gHardTimeCut);

  if (nGood == 0) { gSolved = kFALSE; printf("ERROR: no events survive -- cannot solve\n"); return; }

  TDecompSVD lamD(Ay);
  gDecomposeOk = lamD.Decompose();
  TVectorD x = bVec;
  gSolveOk = lamD.Solve(x);
  gLambdaSolved = x;
  gSolved = gSolveOk;
  printf("Decompose OK? %d, Solve OK? %d\n", gDecomposeOk, gSolveOk);
}

// ===========================================================================
// Section 6: JSON output
// ===========================================================================

void SaveJSON(const TString &tag) {
  gSystem->mkdir(gOutDir, kTRUE);
  TString path = Form("%s/lambda_shms_%s.json", gOutDir.Data(), tag.Data());
  std::ofstream out(path.Data());
  out << "{\n  \"tag\": \"" << tag << "\",\n  \"vpcalib_tag\": \"" << (gVpcalibTag.Length() ? gVpcalibTag : EffectiveTag())
      << "\",\n  \"hard_time_cut\": " << gHardTimeCut << ",\n  \"n_candidates\": " << (Long64_t)gCandidates.size()
      << ",\n  \"n_good_events\": " << gNGoodEvents << ",\n  \"decompose_ok\": " << (gDecomposeOk ? "true" : "false")
      << ",\n  \"solve_ok\": " << (gSolveOk ? "true" : "false") << ",\n  \"paddles\": [\n";
  Bool_t first = kTRUE;
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
      if (IsPermanentlyOff(ipl, 0, ipad)) continue;
      Int_t gidx = kRefPad[ipl] + ipad - 1;
      Double_t v = gSolved ? ((gidx == kRefGlobalIdx) ? 0.0 : gLambdaSolved[gidx]) : 0.0;
      if (!first) out << ",\n";
      first = kFALSE;
      out << "    {\"plane\": \"" << gPlaneNames[ipl] << "\", \"paddle\": " << ipad
          << ", \"global_index\": " << gidx << ", \"is_reference\": " << (gidx == kRefGlobalIdx ? "true" : "false")
          << ", \"lambda\": " << v << "}";
    }
  out << "\n  ]\n}\n";
  out.close();
  printf("Wrote %s\n", path.Data());
}

// ===========================================================================
// Section 7: per-plane review drawing
// ===========================================================================

void DrawPlane() {
  gCanvas->cd();
  gCanvas->Clear();
  // Widen the right margin so the stat box lives entirely outside the
  // plot frame -- the lambda scatter spans nearly the whole panel (top
  // to bottom, left to right), so there's no free interior region to
  // tuck a box into the way the TDC-vs-amplitude plots had.
  gPad->SetRightMargin(0.32);
  Int_t ipl = gPlaneIdx;
  Int_t n = gNbars[ipl];
  std::vector<Double_t> xs, ys;
  for (Int_t ipad = 1; ipad <= n; ipad++) {
    xs.push_back(ipad);
    if (IsPermanentlyOff(ipl, 0, ipad)) { ys.push_back(0.0); continue; }
    Int_t gidx = kRefPad[ipl] + ipad - 1;
    ys.push_back(gSolved ? ((gidx == kRefGlobalIdx) ? 0.0 : gLambdaSolved[gidx]) : 0.0);
  }
  std::vector<Double_t> refYs;
  if (gHaveReferenceLambda)
    for (Int_t ipad = 1; ipad <= n; ipad++) {
      Int_t gidx = kRefPad[ipl] + ipad - 1;
      refYs.push_back(IsPermanentlyOff(ipl, 0, ipad) ? 0.0 : gLambdaReference[gidx]);
    }

  TH1F frame("frame", Form("Inter-Plane Timing Offset (#lambda) -- Plane %s;Paddle;#lambda (ns)", gPlaneNames[ipl].Data()),
             n, 0.5, n + 0.5);
  Double_t ymin = 0, ymax = 0;
  for (auto v : ys) { ymin = TMath::Min(ymin, v); ymax = TMath::Max(ymax, v); }
  for (auto v : refYs) { ymin = TMath::Min(ymin, v); ymax = TMath::Max(ymax, v); }
  Double_t pad = 0.2 * TMath::Max(0.1, ymax - ymin);
  frame.SetMinimum(ymin - pad); frame.SetMaximum(ymax + pad);
  frame.SetStats(kFALSE);
  frame.DrawClone();

  // Draw the new solve FIRST, reference SECOND (on top) -- same
  // draw-order reasoning as the other apps' curve overlays, so the
  // reference markers are never hidden by coincident new-solve markers.
  TGraph g(n, xs.data(), ys.data());
  g.SetMarkerStyle(21); g.SetMarkerColor(kBlue + 1); g.SetMarkerSize(1.2);
  g.DrawClone("P SAME");

  if (gHaveReferenceLambda) {
    TGraph gref(n, xs.data(), refYs.data());
    gref.SetMarkerStyle(24); // open circle
    gref.SetMarkerColor(kBlack);
    gref.SetMarkerSize(1.4);
    gref.DrawClone("P SAME");
  }

  TPaveText pt(0.70, 0.35, 0.99, 0.92, "NDC");
  pt.SetFillColor(kWhite); pt.SetTextAlign(12); pt.SetFillStyle(1001);
  pt.AddText(Form("Plane %d / 4: %s", gPlaneIdx + 1, gPlaneNames[ipl].Data()));
  pt.AddText(Form("Hard time cut: %.1f ns", gHardTimeCut));
  pt.AddText(Form("Good events in solve: %lld / %zu cached", gNGoodEvents, gCandidates.size()));
  pt.AddText(Form("Decompose/Solve OK: %d / %d", gDecomposeOk, gSolveOk));
  if (gHaveReferenceLambda) pt.AddText("blue filled square = new solve, black open circle = reference");
  if (ipl == kRefPlane) pt.AddText(Form("Paddle %d is the fixed reference (#lambda #equiv 0)", kRefPaddle));
  pt.DrawClone();

  gPad->Modified(); gPad->Update();
}

// ===========================================================================
// Section 8: interaction
// ===========================================================================

void GoNext() { if (gPlaneIdx < 3) { ++gPlaneIdx; DrawPlane(); } }
void GoPrev() { if (gPlaneIdx > 0) { --gPlaneIdx; DrawPlane(); } }

void RedoWithNewCut() {
  printf("\nCurrent hard time cut: %.1f ns. Enter a new value: ", gHardTimeCut);
  fflush(stdout);
  std::string line;
  std::getline(std::cin, line);
  TString input(line.c_str()); input = input.Strip(TString::kBoth);
  if (input.Length() == 0) { printf("(cancelled)\n"); return; }
  try { gHardTimeCut = std::stod(input.Data()); } catch (...) { printf("Could not parse '%s'\n", input.Data()); return; }
  Solve();
  DrawPlane();
}

void PauseSave() { SaveJSON(EffectiveTag()); WriteMergedParam(EffectiveTag()); printf("[progress saved]\n"); }

void SaveFinish() {
  PauseSave();
  printf("[saved and finished -- rename %s/SHMS/HODO/phodo_Vpcalib_%s.param to phodo_Vpcalib.param for hcana]\n",
         gParamDir.Data(), EffectiveTag().Data());
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "SHMS Lambda (Inter-Plane) Review", 20, 20);
  gControlBar->AddButton("<< Prev", "GoPrev();", "Previous plane");
  gControlBar->AddButton("Next >>", "GoNext();", "Next plane");
  gControlBar->AddButton("Redo with new cut...", "RedoWithNewCut();", "Re-solve with a different hard TDC time cut");
  gControlBar->AddButton("Pause (save)", "PauseSave();", "Write progress, keep going");
  gControlBar->AddButton("Save && Finish", "SaveFinish();", "Write the merged 4-block param file");
  gControlBar->Show();
}

// ===========================================================================
// Section 9: non-interactive batch mode
// ===========================================================================

void RunNonInteractive() {
  Solve();
  gSystem->mkdir(gOutDir, kTRUE);
  TString pdfPath = Form("%s/lambda_shms_%s_summary.pdf", gOutDir.Data(), EffectiveTag().Data());
  TCanvas c("c", "c", 1150, 700);
  gCanvas = &c;
  Bool_t firstPage = kTRUE;
  for (Int_t ipl = 0; ipl < 4; ipl++) {
    gPlaneIdx = ipl;
    DrawPlane();
    TString pagePath = pdfPath;
    if (firstPage) { pagePath += "("; firstPage = kFALSE; }
    c.Print(pagePath);
  }
  TString closePath = pdfPath + ")";
  c.Print(closePath);
  printf("Wrote %s\n", pdfPath.Data());
  if (!gCompareOnly) {
    SaveJSON(EffectiveTag());
    WriteMergedParam(EffectiveTag());
  } else {
    printf("[compareOnly: PDF written for comparison, nothing saved -- merged .param file left untouched]\n");
  }
}

// ===========================================================================
// Section 10: main entry point
// ===========================================================================

void lambda_calib_app(Int_t run = 0, TString runs = "", TString vpcalibTag = "",
                       Bool_t nonInteractive = kFALSE, TString outputTag = "",
                       Double_t dpCut = 1.0e9, Double_t hardTimeCut = 125.0,
                       TString rootFile = "", TString outDir = "./lambda_qa",
                       TString paramDir = "../../PARAM",
                       TString referenceTag = "",
                       Bool_t applyCalCut = kTRUE, Bool_t applyCerCut = kTRUE,
                       Bool_t applyTrackCut = kTRUE, Bool_t compareOnly = kFALSE) {
  if (run == 0 && runs.Length() == 0) {
    printf("ERROR: must supply a run number or runs list, e.g. lambda_calib_app(0, \"26483-26488\")\n"); return;
  }
  gVpcalibTag = vpcalibTag;
  gOutputTag = outputTag;
  gOutDir = outDir; gParamDir = paramDir;
  gDpCut = dpCut; gHardTimeCut = hardTimeCut;
  gApplyCalCut = applyCalCut;
  gApplyCerCut = applyCerCut;
  gApplyTrackCut = applyTrackCut;
  gReferenceTag = referenceTag;
  gCompareOnly = compareOnly;
  if (compareOnly && referenceTag.Length() == 0)
    printf("WARNING: compareOnly=true but no referenceTag given -- there will be nothing to compare against.\n");

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

  gStyle->SetOptStat(0);
  gROOT->SetBatch(nonInteractive);

  if (!LoadInputVpcalib()) return;
  LoadReferenceLambda();
  BuildCandidates();

  if (nonInteractive) { RunNonInteractive(); return; }

  Solve();
  gPlaneIdx = 0;
  gCanvas = new TCanvas("lambdaCanvas", "SHMS Lambda (Inter-Plane) Calibration", 1150, 700);
  DrawPlane();
  MakeControlBar();

  printf("\n=== run(s) %s (output tag: %s), vpcalib source tag: %s ===\n",
         gRunTag.Data(), EffectiveTag().Data(), gVpcalibTag.Length() ? gVpcalibTag.Data() : EffectiveTag().Data());
  printf("Use the 'SHMS Lambda (Inter-Plane) Review' panel to step through planes, adjust the cut, and save.\n");
}
