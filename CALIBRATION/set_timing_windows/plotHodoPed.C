// plotHodoPed.C
//
// Usage:
//   root -l -b -q plotHodoPed.C
// or from within ROOT:
//   root [0] .x plotHodoPed.C
//
// For each run in `runs` below, reads
// ROOTfiles/coin_replay_production_<run>_2000000_0.root and draws
// P.hod.1x.GoodPosAdcPed vs P.hod.1x.GoodPosAdcTdcDiffTime side by
// side on one canvas, then saves it as result_compare.pdf
//
// The x-axis binning (100, -30, 30) is fixed per your original code.
// The y-axis range can either be set manually (yRange below) or
// determined automatically from the data itself.
//
//   autoRange = true  -> ymin/ymax are taken from
//                         T->GetMinimum()/GetMaximum() on the branch,
//                         with a small padding margin, then clamped
//                         to yHardMax. nbins from yRange is still used.
//   autoRange = false -> the ymin/ymax in yRange are used as-is.

void plotHodoPed(bool autoRange = false, double yHardMax = 1e3)
{
    // --- Runs to compare ---
    std::vector<int> runs = {26088, 26169};

    // --- Descriptive label appended to each run's plot title ---
    std::map<int, TString> runLabel = {
        {26088, "carbon elastics"},
        {26169, "N-Delta LH2 1d"}
    };

    TString yVar = "P.hod.1x.GoodPosAdcPed";
    TString xVar = "P.hod.1x.GoodPosAdcTdcDiffTime";

    // --- Fixed x-axis binning ---
    const int    nbinsX = 100;
    const double xmin   = 10;
    const double xmax   = 30;

    // --- Y-axis binning: {nbins, ymin, ymax} (used as-is if autoRange = false,
    // or just for nbins if autoRange = true) ---
    struct YRange { int nbins; double ymin; double ymax; };
    YRange yRange = {40, 40, 80};

    // --- Canvas ---
    gStyle->SetOptStat(0); // mute the stat box on all histograms
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 700);
    c1->Clear();
    c1->Divide(2, 1);

    std::vector<TFile*> files; // keep open until printed
    std::vector<TH2F*>  hists;

    for (size_t k = 0; k < runs.size(); k++) {
        int run = runs[k];

        TString infile = Form("ROOTfiles/coin_replay_production_%d_2000000_0.root", run);
        TFile *f = TFile::Open(infile);
        if (!f || f->IsZombie()) {
            printf("ERROR: could not open %s\n", infile.Data());
            continue;
        }
        files.push_back(f);

        TTree *T = (TTree*)f->Get("T");
        if (!T) {
            printf("ERROR: could not find tree \"T\" in %s\n", infile.Data());
            continue;
        }

        c1->cd(k + 1);

        double ymin = yRange.ymin;
        double ymax = yRange.ymax;

        if (autoRange) {
            double dmin = T->GetMinimum(yVar);
            double dmax = T->GetMaximum(yVar);

            if (dmin == dmax) {
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

        TString hname = Form("h_%d", run);
        TString draw  = Form("%s:%s>>%s(%d,%f,%f,%d,%f,%f)",
                              yVar.Data(), xVar.Data(), hname.Data(),
                              nbinsX, xmin, xmax,
                              yRange.nbins, ymin, ymax);

        T->Draw(draw, "", "colz");

        TH2F *h = (TH2F*)gDirectory->Get(hname);
        if (h) {
            h->SetTitle(Form("Run %d (%s)", run, runLabel[run].Data()));
            h->SetStats(0);
            h->GetXaxis()->SetTitle(xVar);
            h->GetYaxis()->SetTitle(yVar);
            h->GetXaxis()->SetTitleSize(0.045);
            h->GetYaxis()->SetTitleSize(0.045);
            hists.push_back(h);
        }
    }

    c1->Update();

    // --- Save output ---
    TString outfile = "result_compare_carbon_to_production.pdf";
    c1->Print(outfile);

    printf("Saved %s\n", outfile.Data());
}
