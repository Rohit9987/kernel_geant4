// plotXZ_polar_contours_scaled.C
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <cmath>

static TH2D* MakeHist(const char* name, int nBins, double range_mm) {
  auto* h = new TH2D(name, "", nBins, -range_mm, range_mm, nBins, 0.0, range_mm); // upper half only (Z>=0)
  h->SetDirectory(nullptr);
  return h;
}

void contours(const char* fname = "DoseKernel_6.0MeV.root",
                                  double range_gcm2 = 10.0,   // shown radius (g cm^-2) ~ cm in water
                                  double binWidth_mm = 0.5,   // binning in mm
                                  double yTol_mm = 0.25,      // |y| slice
                                  bool logLevels = true,
                                  int  nLevels = 6,
                                  double minFrac = 1e-4,      // of max
                                  double maxFrac = 0.5,       // of max
                                  const char* units = "arb.") // label for iso values
{
  const double mm_per_gcm2 = 10.0;            // water: 1 g cm^-2 ≡ 10 mm
  const double Rmax_mm = range_gcm2 * mm_per_gcm2;

  TFile* file = TFile::Open(fname, "READ");
  if (!file || file->IsZombie()) { std::cerr << "Error opening " << fname << "\n"; return; }
  TTree* tree = static_cast<TTree*>(file->Get("DoseData"));
  if (!tree) { std::cerr << "TTree 'DoseData' not found!\n"; file->Close(); return; }

  // Fast I/O: cache + only needed branches
  tree->SetCacheSize(128*1024*1024);
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("X", 1);
  tree->SetBranchStatus("Y", 1);
  tree->SetBranchStatus("Z", 1);
  tree->SetBranchStatus("edep", 1);
  tree->AddBranchToCache("X", kTRUE);
  tree->AddBranchToCache("Y", kTRUE);
  tree->AddBranchToCache("Z", kTRUE);
  tree->AddBranchToCache("edep", kTRUE);

  double x=0, y=0, z=0, edep=0;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Z", &z);
  tree->SetBranchAddress("edep", &edep);

  const int nBins = static_cast<int>(std::floor((2.0*Rmax_mm)/binWidth_mm + 0.5));
  TH2D* h = MakeHist("hXZ_All", nBins, Rmax_mm);

  const Long64_t N = tree->GetEntries();
  const Long64_t printEvery = std::max<Long64_t>(N/50, 100000);
  for (Long64_t i=0; i<N; ++i) {
    tree->GetEntry(i);
    if ((i % printEvery) == 0) std::cout << "Processed " << i << " / " << N << '\r' << std::flush;
    if (std::abs(y) > yTol_mm) continue;
    if (z < 0.0) continue;                                  // only upper half
    if (std::hypot(x, z) > Rmax_mm) continue;               // clip to shown radius
    h->Fill(x, z, edep);
  }
  std::cout << "\nFilled.\n";

  // Build contour levels (fractions of max)
  const double hmax = h->GetMaximum();
  if (hmax <= 0) { std::cerr << "Histogram empty.\n"; file->Close(); return; }
  std::vector<double> levels; levels.reserve(nLevels);
  if (logLevels) {
    const double a = std::log(minFrac), b = std::log(maxFrac);
    for (int i=0; i<nLevels; ++i) levels.push_back(hmax * std::exp(a + (b-a)*i/(nLevels-1)));
  } else {
    for (int i=0; i<nLevels; ++i) levels.push_back(hmax * (minFrac + (maxFrac-minFrac)*i/(nLevels-1)));
  }
  h->SetContour(nLevels, levels.data());
  h->SetLineColor(kBlack);
  h->SetLineWidth(3);

  // Canvas
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cPolarScaled", "Kernel isolines (polar)", 950, 750);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.07);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetFixedAspectRatio();

  // Hide frame and axes completely (removes the “top lines”)
  gPad->SetFrameBorderMode(0);
  gPad->SetFrameLineColor(0);
  h->SetTitle("");
  h->GetXaxis()->SetLabelSize(0);
  h->GetYaxis()->SetLabelSize(0);
  h->GetXaxis()->SetTickLength(0);
  h->GetYaxis()->SetTickLength(0);
  h->Draw("CONT4");                         // contours only, no color

  // ---- Distance scale (bottom TGaxis in g cm^-2) ----
  // Place a horizontal axis just below z=0 by 6% of R
  const double axYOffset = -0.06 * Rmax_mm;
  TGaxis* ax = new TGaxis(-Rmax_mm, axYOffset, +Rmax_mm, axYOffset,
                          -range_gcm2, +range_gcm2, 510, "-"); // numeric labels in g cm^-2
  ax->SetTitle("Radius (g cm^{-2})");
  ax->CenterTitle(true);
  ax->SetTitleSize(0.04);
  ax->SetLabelSize(0.032);
  ax->SetLineColor(kBlack);
  ax->SetTextFont(42);
  ax->Draw();

  // Draw the semicircle arc (outer boundary) lightly
  TEllipse* outer = new TEllipse(0,0, Rmax_mm, Rmax_mm, 0, 180);
  outer->SetFillStyle(0);
  outer->SetLineStyle(1);
  outer->SetLineColor(kBlack);
  outer->SetLineWidth(1);
  outer->Draw("same");

  // Optional: a few radius guide arcs (dashed) — comment out if not wanted
  /*
  for (double r_gcm2 : {2.0,4.0,6.0,8.0}) {
    double R = r_gcm2 * mm_per_gcm2;
    auto* arc = new TEllipse(0,0, R, R, 0, 180);
    arc->SetFillStyle(0); arc->SetLineStyle(2); arc->SetLineColor(kGray+1); arc->SetLineWidth(1);
    arc->Draw("same");
  }
  */

  // ---- Legend for isodose lines ----
  // Show both absolute and fraction of max
  TLegend* L = new TLegend(0.73, 0.70, 0.96, 0.94); // NDC coords
  L->SetBorderSize(0);
  L->SetFillStyle(0);
  L->SetTextFont(42);
  L->SetTextSize(0.030);
  L->SetHeader("Isodose levels", "C");
  for (int i=nLevels-1; i>=0; --i) { // highest first
    const double frac = levels[i]/hmax;
    TString lab = Form("%.2g %s  (%.2g×max)", levels[i], units, frac);
    // Use a dummy line sample for display
    TLine* sample = new TLine(0,0,1,0); sample->SetLineColor(kBlack); sample->SetLineWidth(3);
    L->AddEntry(sample, lab, "l");
  }
  L->Draw();

  c->SaveAs("Kernel_isolines_polar_scaled.png");
  c->SaveAs("Kernel_isolines_polar_scaled.pdf");
  std::cout << "Saved: Kernel_isolines_polar_scaled.(png|pdf)\n";

  file->Close();
}

