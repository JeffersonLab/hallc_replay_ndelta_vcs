// plotHodo1x.C
//
// Usage:
//   root -l -b -q 'plotHodo1x.C(12345)'
// or from within ROOT:
//   root [0] .x plotHodo1x.C(12345)
//
// Reads ROOTfiles/coin_replay_production_<run>_2000000_0.root, draws
// 2D correlation plots of several P.hod.1x variables vs
// P.hod.1x.GoodPosAdcTdcDiffTime on a 2x5 canvas, and saves the
// result as result_<run>.pdf
//
// The x-axis binning (100, -30, 30) is fixed per your original code.
// The y-axis range can either be set manually per-variable (yRange
// array below) or determined automatically from the data itself.
//
//   autoRange = true  -> ymin/ymax for each variable are taken from
//                         T->GetMinimum()/GetMaximum() on that branch,
//                         with a small padding margin, then clamped
//                         to yHardMax (a few branches like PulseInt
//                         can have huge outlier tails). The nbins in
//                         yRange is still used for binning.
//   autoRange = false -> the ymin/ymax in yRange are used as-is
//                         (original fixed-range behavior).

void plotHodo1x(int run, bool autoRange = true, double yHardMax = 1e3)
{
    // --- Input file ---
    TString infile = Form("ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) {
        printf("ERROR: could not open %s\n", infile.Data());
        return;
    }

    TTree *T = (TTree*)f->Get("T");
    if (!T) {
        printf("ERROR: could not find tree \"T\" in %s\n", infile.Data());
        return;
    }

    // --- Variables to plot (y-axis), each vs GoodPosAdcTdcDiffTime (x-axis) ---
    std::vector<TString> yVars = {
        "P.hod.1x.GoodPosAdcHitUsed",
        "P.hod.1x.GoodPosAdcMult",
        "P.hod.1x.GoodPosAdcPed",
        "P.hod.1x.GoodPosAdcPulseAmp",
        "P.hod.1x.GoodPosAdcPulseInt",
        "P.hod.1x.GoodPosAdcPulseTime",
        "P.hod.1x.GoodPosTdcTimeCorr",
        "P.hod.1x.GoodPosTdcTimeTOFCorr",
        "P.hod.1x.GoodPosTdcTimeUnCorr",
        "P.hod.1x.GoodPosTdcTimeWalkCorr"
    };

    TString xVar = "P.hod.1x.GoodPosAdcTdcDiffTime";

    // --- Fixed x-axis binning ---
    const int    nbinsX = 100;
    const double xmin   = -30;
    const double xmax   = 30;

    // --- Per-variable y-axis binning: {nbinsY, ymin, ymax} ---
    // Defaults to (30, 0, 30) for everything, matching your original
    // code. Edit individual rows below as needed -- order matches yVars.
    struct YRange { int nbins; double ymin; double ymax; };
    std::vector<YRange> yRange = {
        {3, 0.5, 3.5},   // GoodPosAdcHitUsed
        {3, 0.5, 3.5},   // GoodPosAdcMult
        {40, 40, 80},   // GoodPosAdcPed
        {200, 0, 200},   // GoodPosAdcPulseAmp
        {200, 0, 200},   // GoodPosAdcPulseInt
        {200, 0, 200},   // GoodPosAdcPulseTime
        {200, 0, 200},   // GoodPosTdcTimeCorr
        {200, 0, 200},   // GoodPosTdcTimeTOFCorr
        {200, 0, 200},   // GoodPosTdcTimeUnCorr
        {200, 0, 200}    // GoodPosTdcTimeWalkCorr
    };

    if (yRange.size() != yVars.size()) {
        printf("ERROR: yRange and yVars size mismatch\n");
        return;
    }

    // --- Canvas ---
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 2000);
    c1->Clear();
    c1->Divide(2, 5);

    int i = 1;
    std::vector<TH2F*> hists; // keep pointers alive until printed

    for (size_t k = 0; k < yVars.size(); k++) {
        c1->cd(i);

        double ymin = yRange[k].ymin;
        double ymax = yRange[k].ymax;

        if (autoRange) {
            double dmin = T->GetMinimum(yVars[k]);
            double dmax = T->GetMaximum(yVars[k]);

            if (dmin == dmax) {
                // Constant branch -- pad by 1 so the histogram isn't degenerate
                dmin -= 1.0;
                dmax += 1.0;
            } else {
                double pad = 0.05 * (dmax - dmin);
                dmin -= pad;
                dmax += pad;
            }
            ymin = dmin;
            ymax = dmax;

            // Clamp to a hard ceiling (some ADC branches have long tails)
            if (ymax > yHardMax) ymax = yHardMax;
            if (ymin > ymax)     ymin = ymax - 1.0; // safety, keeps range valid
        }

        TString hname = Form("h%d", i);
        TString draw  = Form("%s:%s>>%s(%d,%f,%f,%d,%f,%f)",
                              yVars[k].Data(), xVar.Data(), hname.Data(),
                              nbinsX, xmin, xmax,
                              yRange[k].nbins, ymin, ymax);

        T->Draw(draw, "", "colz");

        TH2F *h = (TH2F*)gDirectory->Get(hname);
        if (h) {
            h->SetTitle(Form("Run %d: %s vs %s", run, yVars[k].Data(), xVar.Data()));
            h->GetXaxis()->SetTitle(xVar);
            h->GetYaxis()->SetTitle(yVars[k]);
            h->GetXaxis()->SetTitleSize(0.045);
            h->GetYaxis()->SetTitleSize(0.045);
            hists.push_back(h);
        }

        i++;
    }

    c1->Update();

    // --- Save output ---
    TString outfile = Form("result_%d.pdf", run);
    c1->Print(outfile);

    printf("Saved %s\n", outfile.Data());
}
