// lambda_calib_app.C (HMS)
//
// Fits the HMS Hodoscope's per-paddle inter-plane timing offset
// (hhodo_LCoeff) via the "Hodoscope Planes Time Difference Corrections"
// step, replicating hms_hodo_calib/fitHodoCalib.C's final section.
//
// Same single global least-squares + SVD solve as the SHMS sibling app
// (see ../shms_hodo_calib/lambda_calib_app.C for the general explanation
// of the Ay*lambda=bVec system) -- differences specific to HMS:
//   * Candidate selection is SIMPLER: the original doesn't reuse any
//     mean-time-window derivation here at all (unlike SHMS). It gates
//     purely on the tree's own "H.hod.<plane>.nhits" branch == 1 for
//     ALL 4 planes, plus a track-position sanity check -- so this app
//     needs only ONE tree pass, not two.
//   * npar = 52 (16+10+16+10), refPad = {0,16,26,42}.
//   * z0 = {77.83, 97.52, 298.82, 318.51} cm, dz = {2.12,2.12,2.12,2.12}
//     (all the same sign, unlike SHMS's alternating-sign dz).
//   * Hard TDC time cut default 100 ns (matches the uniform 100ns used
//     throughout the HMS script, vs SHMS's 125ns for this stage).
//   * Reference paddle is still 1x-paddle-7 (global index 6, same as
//     SHMS, since refPad[0]=0 for both).
//
// REQUIRES an existing hhodo_Vpcalib_<vpcalibTag>.param (from
// vpcalib_app.C) for hhodo_cableFit.
//
// Usage:
//   root -l 'lambda_calib_app.C(23859)'                        // interactive review
//   root -l -b -q 'lambda_calib_app.C(23859, "", "", true)'    // batch: solve + write
//
// NOTE on paddle-finding: given nhits==1, exactly one paddle per plane
// should carry real (non-sentinel) TW-corrected times; this app scans
// for the first paddle whose times are below a loose 1e8 sentinel (not
// the hard time cut itself) so the actual hard cut stays a re-solvable
// knob afterward, same idea as caching in the SHMS app. The original
// script instead applies the hard cut directly while searching (and,
// due to an unconditional loop overwrite, keeps the LAST matching
// paddle rather than the first) -- functionally equivalent given
// nhits==1 should guarantee at most one real candidate either way.
//
// Controls (interactive): same as the SHMS app --
//   << Prev / Next >>, Redo with new cut..., Pause (save), Save && Finish
//
// Outputs (under outDir, default "./lambda_qa"):
//   lambda_hms_<tag>.json, lambda_hms_<tag>_summary.pdf
// Outputs (under paramDir, default "../../PARAM"):
//   HMS/HODO/hhodo_Vpcalib_<tag>.param -- MERGED: velFit/cableFit/sigma
//   carried over from vpcalibTag's file + the new hhodo_LCoeff block.
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
// Section 1: HMS constants (from hms_hodo_calib/fitHodoCalib.C)
// ===========================================================================

static const TString gPlaneNames[4] = {"1x", "1y", "2x", "2y"};
static const Int_t   gNbars[4]      = {16, 10, 16, 10};
static const Int_t   nBarsMax       = 16;
static const Int_t   kRefPad[4]     = {0, 16, 26, 42};
static const Int_t   kNpar          = 52;
static const Int_t   kRefPlane = 0, kRefPaddle = 7;
static const Int_t   kRefGlobalIdx = kRefPad[kRefPlane] + kRefPaddle - 1; // = 6

static const Double_t kPlaneZ0[4] = {77.83, 97.52, 298.82, 318.51};
static const Double_t kDz[4] = {2.12, 2.12, 2.12, 2.12};
static const Double_t kLightSpeed = 30.0; // cm/ns

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
static const Double_t kTrackPosSanity = 200.0;

TString ResolvedRunPath(Int_t run) {
  // Original location -- switch back once the slurm-queued replay has
  // landed in the normal ROOTfiles directory:
  return Form("/volatile/hallc/alphaE/ndelta_vcs2/calib/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
  // return Form("/volatile/hallc/alphaE/ndelta_vcs2/selected_runs_before_hodo_calibration/ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
}

// ===========================================================================
// Section 2: globals
// ===========================================================================

Int_t    gRun = 0;
TString  gRootFile, gOutDir = "./lambda_qa", gParamDir = "../../PARAM";
std::vector<Int_t> gRunList;
TString  gRunTag, gOutputTag, gVpcalibTag;
Double_t gDpCut = 1.0e9;
Double_t gHardTimeCut = 100.0; // matches the uniform 100ns cut used throughout the HMS script
Bool_t   gApplyCalCut = kTRUE;   // H.cal.etracknorm cut -- turn off for a proton-arm run
Bool_t   gApplyCerCut = kTRUE;   // H.cer.npeSum cut
Bool_t   gApplyTrackCut = kTRUE; // H.dc.ntrack cut

TString gVelBlockRaw, gCableBlockRaw, gPosSigmaBlockRaw, gNegSigmaBlockRaw;
Double_t gCableFit[4][nBarsMax] = {{0}};

struct EventCand {
  Int_t pad[4];
  Double_t twPos[4], twNeg[4];
  Double_t trackX[4], trackY[4];
};
std::vector<EventCand> gCandidates;

TVectorD gLambdaSolved(kNpar);
Bool_t   gSolved = kFALSE;
Long64_t gNGoodEvents = 0;
Bool_t   gDecomposeOk = kFALSE, gSolveOk = kFALSE;

// Optional reference: an existing hhodo_LCoeff (from a prior solve, or
// the vanilla PARAM/HMS/HODO/hhodo_Vpcalib.param) loaded for visual
// comparison against a fresh solve.
Double_t gLambdaReference[kNpar] = {0};
Bool_t   gHaveReferenceLambda = kFALSE;
TString  gReferenceTag; // "" = no reference; "vanilla" = the untagged file; else a staged tag
Bool_t   gCompareOnly = kFALSE; // batch mode: write the summary PDF (new solve vs. reference) but never touch JSON/.param output

Int_t gPlaneIdx = 0;
TCanvas *gCanvas = nullptr;
TControlBar *gControlBar = nullptr;

TString EffectiveTag() { return gOutputTag.Length() > 0 ? gOutputTag : gRunTag; }

// ===========================================================================
// Section 3: legacy param I/O
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
  if (staged) return Form("%s/HMS/HODO/hhodo_Vpcalib_%s.param", gParamDir.Data(), tag.Data());
  return Form("%s/HMS/HODO/hhodo_Vpcalib.param", gParamDir.Data());
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
      printf("ERROR: no hhodo_Vpcalib param file found (tried tag '%s' and the vanilla file). Run vpcalib_app.C first.\n", tag.Data());
      return kFALSE;
    }
    test2.close();
  } else test.close();

  Double_t cable[4][nBarsMax] = {{0}};
  std::ifstream in(path.Data());
  ParseParamBlock(in, "hhodo_cableFit", cable);
  for (Int_t ipl = 0; ipl < 4; ipl++)
    for (Int_t ipad = 0; ipad < nBarsMax; ipad++)
      gCableFit[ipl][ipad] = cable[ipl][ipad];

  gVelBlockRaw      = ExtractRawBlock(path, "hhodo_velFit");
  gCableBlockRaw    = ExtractRawBlock(path, "hhodo_cableFit");
  gPosSigmaBlockRaw = ExtractRawBlock(path, "hhodo_PosSigma");
  gNegSigmaBlockRaw = ExtractRawBlock(path, "hhodo_NegSigma");
  printf("Loaded cableFit (and carried velFit/sigma) from %s\n", path.Data());
  return kTRUE;
}

// Loads an existing hhodo_LCoeff block for visual comparison only -- has
// no effect on the fresh solve. gReferenceTag == "vanilla" loads the
// untagged file; otherwise it's treated as a staged tag (falls back to
// vanilla if that staged file doesn't exist).
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
  ParseParamBlock(in, "hhodo_LCoeff", lcoeff);
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
  out << "; HMS Hodoscope Parameter File Containing propagation velocities, cable offsets, " << std::endl;
  out << "; sigma resolution, and inter-plane timing (LCoeff) parameters per paddle " << std::endl;
  out << Form("; tag %s ", tag.Data()) << std::endl << std::endl;

  out << ";hhodo_velFit" << std::endl;
  out << (gVelBlockRaw.Length() ? gVelBlockRaw.Data() : "hhodo_velFit = ; (none carried over)\n") << std::endl;
  out << ";hhodo_cableFit" << std::endl;
  out << (gCableBlockRaw.Length() ? gCableBlockRaw.Data() : "hhodo_cableFit = ; (none carried over)\n") << std::endl;
  out << ";hhodo_PosSigma" << std::endl;
  out << (gPosSigmaBlockRaw.Length() ? gPosSigmaBlockRaw.Data() : "hhodo_PosSigma = ; (none carried over)\n") << std::endl;
  out << ";hhodo_NegSigma" << std::endl;
  out << (gNegSigmaBlockRaw.Length() ? gNegSigmaBlockRaw.Data() : "hhodo_NegSigma = ; (none carried over)\n") << std::endl;

  out << ";Timing Corrections Per Paddle, where 1X Paddle 7 has been set as the reference paddle" << std::endl;
  out << "hhodo_LCoeff = ";
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
// Section 4: build the candidate list -- SINGLE tree pass, gated on the
//   tree's own "nhits==1" branch (no mean-time-window derivation needed
//   for HMS, unlike the SHMS app).
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
  T->SetBranchAddress("H.gtr.dp", &dpVal);

  Double_t twPos[4][nBarsMax], twNeg[4][nBarsMax];
  Double_t trackX[4], trackY[4], nhits[4];
  for (Int_t ipl = 0; ipl < 4; ipl++) {
    TString base = "H.hod." + gPlaneNames[ipl];
    T->SetBranchAddress(base + ".GoodPosTdcTimeWalkCorr", twPos[ipl]);
    T->SetBranchAddress(base + ".GoodNegTdcTimeWalkCorr", twNeg[ipl]);
    T->SetBranchAddress(base + ".TrackXPos", &trackX[ipl]);
    T->SetBranchAddress(base + ".TrackYPos", &trackY[ipl]);
    T->SetBranchAddress(base + ".nhits", &nhits[ipl]);
  }

  Long64_t nentries = T->GetEntries();
  printf("\n=== lambda (HMS) single pass: %lld events ===\n", nentries);
  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 200000 == 0 && i != 0) printf("  %lld / %lld...\n", i, nentries);
    if (TMath::Abs(dpVal) >= gDpCut) continue;
    if (!((!gApplyCalCut || pid1 > kPidCut1Low) && (!gApplyCerCut || pid2 > kPidCut2Low) && (!gApplyTrackCut || pid3 > kPidCut3Low))) continue;

    // "Require each plane to have ONLY a SINGLE HIT" -- gated directly
    // on the tree's own nhits branch, not a re-derived flag.
    Bool_t singleHitAll4 = kTRUE;
    for (Int_t ipl = 0; ipl < 4; ipl++) if (nhits[ipl] != 1) singleHitAll4 = kFALSE;
    if (!singleHitAll4) continue;

    Bool_t hodTrk = kTRUE;
    for (Int_t ipl = 0; ipl < 4; ipl++)
      if (!(trackX[ipl] < kTrackPosSanity && trackY[ipl] < kTrackPosSanity)) hodTrk = kFALSE;
    if (!hodTrk) continue;

    // Given nhits==1, exactly one paddle per plane should have real
    // (non-sentinel) TW-corrected times -- find it. The hard time cut
    // itself is applied later in Solve(), so this just locates the
    // paddle; cache raw times for both sides.
    EventCand c;
    Bool_t allFound = kTRUE;
    for (Int_t ipl = 0; ipl < 4; ipl++) {
      Int_t found = -1;
      for (Int_t ipad = 1; ipad <= gNbars[ipl]; ipad++) {
        if (twPos[ipl][ipad - 1] < 1.0e8 && twNeg[ipl][ipad - 1] < 1.0e8) { found = ipad; break; }
      }
      if (found < 0) { allFound = kFALSE; break; }
      c.pad[ipl] = found;
      c.twPos[ipl] = twPos[ipl][found - 1]; c.twNeg[ipl] = twNeg[ipl][found - 1];
      c.trackX[ipl] = trackX[ipl]; c.trackY[ipl] = trackY[ipl];
    }
    if (!allFound) continue;
    gCandidates.push_back(c);
  }
  printf("=== Done: %lld events processed, %zu single-hit-all-4-planes candidates cached ===\n\n", nentries, gCandidates.size());

  if (f) f->Close();
  if (chain) delete chain;
}

// ===========================================================================
// Section 5: accumulate normal equations and solve via SVD
// ===========================================================================

Double_t ZOf(Int_t plane, Int_t paddle) { return kPlaneZ0[plane] + ((paddle - 1) % 2) * kDz[plane]; }

void Solve() {
  TMatrixD Ay(kNpar, kNpar);
  TVectorD bVec(kNpar);
  Ay.Zero(); bVec.Zero();

  Long64_t nGood = 0;
  for (auto &c : gCandidates) {
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

    static const Int_t kPairI[6] = {0, 0, 0, 1, 1, 2};
    static const Int_t kPairJ[6] = {1, 2, 3, 2, 3, 3};
    for (Int_t r = 0; r < 6; r++) {
      Int_t pi = kPairI[r], pj = kPairJ[r];
      Double_t D = -TMath::Sqrt(TMath::Power(x[pi]-x[pj],2) + TMath::Power(y[pi]-y[pj],2) + TMath::Power(z[pi]-z[pj],2));
      Double_t b = D / kLightSpeed - (t[pi] - t[pj]);
      Int_t a = gidx[pi], bIdx = gidx[pj];
      Double_t ca = (a == kRefGlobalIdx) ? 0.0 : 1.0;
      Double_t cb = (bIdx == kRefGlobalIdx) ? 0.0 : -1.0;
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
  TString path = Form("%s/lambda_hms_%s.json", gOutDir.Data(), tag.Data());
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
  // plot frame -- see note in shms_hodo_calib/lambda_calib_app.C.
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
  printf("[saved and finished -- rename %s/HMS/HODO/hhodo_Vpcalib_%s.param to hhodo_Vpcalib.param for hcana]\n",
         gParamDir.Data(), EffectiveTag().Data());
}

void MakeControlBar() {
  if (gControlBar) { delete gControlBar; gControlBar = nullptr; }
  gControlBar = new TControlBar("vertical", "HMS Lambda (Inter-Plane) Review", 20, 20);
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
  TString pdfPath = Form("%s/lambda_hms_%s_summary.pdf", gOutDir.Data(), EffectiveTag().Data());
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
                       Double_t dpCut = 1.0e9, Double_t hardTimeCut = 100.0,
                       TString rootFile = "", TString outDir = "./lambda_qa",
                       TString paramDir = "../../PARAM",
                       TString referenceTag = "",
                       Bool_t applyCalCut = kTRUE, Bool_t applyCerCut = kTRUE,
                       Bool_t applyTrackCut = kTRUE, Bool_t compareOnly = kFALSE) {
  if (run == 0 && runs.Length() == 0) {
    printf("ERROR: must supply a run number or runs list, e.g. lambda_calib_app(23859)\n"); return;
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
  gCanvas = new TCanvas("lambdaCanvas", "HMS Lambda (Inter-Plane) Calibration", 1150, 700);
  DrawPlane();
  MakeControlBar();

  printf("\n=== run(s) %s (output tag: %s), vpcalib source tag: %s ===\n",
         gRunTag.Data(), EffectiveTag().Data(), gVpcalibTag.Length() ? gVpcalibTag.Data() : EffectiveTag().Data());
  printf("Use the 'HMS Lambda (Inter-Plane) Review' panel to step through planes, adjust the cut, and save.\n");
}
