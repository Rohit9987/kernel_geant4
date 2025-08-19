// plot_r2_angular.C
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

struct AngleBand { double loDeg, hiDeg; const char* label; int color; int style; };

void plot_r2_angular(const char* fname = "../physics/DoseKernel_6.0MeV.root",
                     double yTol_mm = 0.25,      // mid-plane slice thickness
                     double rmax_gcm2 = 20.0,    // radial axis extent (g cm^-2)
                     double binw_gcm2 = 0.1,     // radial bin width (g cm^-2)
                     bool   normalizeByAngle = false) // divide by angular width (deg)
{
  // ----------------- open & set up cached I/O -----------------
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) { std::cerr << "Error opening " << fname << "\n"; return; }
  TTree* T = static_cast<TTree*>(f->Get("DoseData"));
  if (!T) { std::cerr << "TTree 'DoseData' not found\n"; f->Close(); return; }

  T->SetCacheSize(128*1024*1024);
  T->SetBranchStatus("*", 0);
  T->SetBranchStatus("X", 1);
  T->SetBranchStatus("Y", 1);
  T->SetBranchStatus("Z", 1);
  T->SetBranchStatus("edep", 1);
  T->AddBranchToCache("X", kTRUE);
  T->AddBranchToCache("Y", kTRUE);
  T->AddBranchToCache("Z", kTRUE);
  T->AddBranchToCache("edep", kTRUE);

  double x=0, y=0, z=0, edep=0;
  T->SetBranchAddress("X", &x);
  T->SetBranchAddress("Y", &y);
  T->SetBranchAddress("Z", &z);
  T->SetBranchAddress("edep", &edep);

  // ----------------- configuration -----------------
  const double mm_per_gcm2 = 10.0; // water
  const double rmax_mm     = rmax_gcm2 * mm_per_gcm2;

  std::vector<AngleBand> bands = {
    { 0.00,  2.56,  "0–2.56^{#circ}",        kBlack,  1 },
    {25.84, 31.79,  "25.84–31.79^{#circ}",  kBlack,  2 }, // dashed
    {45.57, 49.46,  "45.57–49.46^{#circ}",  kBlack,  3 }  // dotted
  };

  const int nbins = std::max(1, int(std::floor(rmax_gcm2/binw_gcm2 + 0.5)));

  // Histograms, one per band
  std::vector<TH1D*> H;
  H.reserve(bands.size());
  for (size_t i=0;i<bands.size();++i) {
    auto* h = new TH1D(Form("h_r2_%zu", i), "", nbins, 0.0, rmax_gcm2);
    h->SetDirectory(nullptr);
    h->Sumw2(); // keep errors if you want
    h->SetLineColor(bands[i].color);
    h->SetLineStyle(bands[i].style);
    h->SetLineWidth(2);
    H.push_back(h);
  }

  // ----------------- loop once over TTree -----------------
  const Long64_t N = T->GetEntries();
  const Long64_t printEvery = std::max<Long64_t>(N/50, 100000);
  for (Long64_t i=0;i<N;++i) {
    T->GetEntry(i);
    if ((i % printEvery) == 0) std::cout << "Processed " << i << "/" << N << "\r" << std::flush;

    if (std::abs(y) > yTol_mm) continue;  // XZ mid-plane
    if (z <= 0.0) continue;               // forward hemisphere like the paper

    const double r_mm = std::hypot(x, z);
    if (r_mm<=0.0 || r_mm > rmax_mm) continue;

    // angle from beam axis (z-axis) in the XZ plane
    // use |x| to fold left/right symmetry, so 0° = on-axis forward.
    const double theta_deg = std::atan2(std::abs(x), z) * 180.0 / TMath::Pi();

    // radius in g cm^-2
    const double r_gcm2 = r_mm / mm_per_gcm2;

    // weight = edep * r^2  (with r in g cm^-2 as per figure)
    const double w = edep * (r_gcm2 * r_gcm2);

    // fill the matching band(s)
    for (size_t k=0;k<bands.size();++k) {
      if (theta_deg >= bands[k].loDeg && theta_deg < bands[k].hiDeg) {
        H[k]->Fill(r_gcm2, w);
      }
    }
  }
  std::cout << "\nFilling done.\n";

  // Optional: normalize by the angular width (so units per degree)
  if (normalizeByAngle) {
    for (size_t k=0;k<bands.size();++k) {
      const double width = std::max(1e-9, bands[k].hiDeg - bands[k].loDeg);
      H[k]->Scale(1.0/width);
    }
  }

  // ----------------- draw overlay -----------------
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_r2", "Dose x r^{2} vs radius", 900, 800);
  c->SetLogy();
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  // Frame histogram for axes
  TH1D* frame = (TH1D*)H[0]->Clone("frame");
  frame->Reset();
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("Radius (g cm^{-2})");
  frame->GetYaxis()->SetTitle("Dose #times r^{2} (arb.)");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.040);

  // Set a reasonable y-range
  double ymin = 1e30, ymax = -1e30;
  for (auto* h : H) { if (h->GetMaximum() > ymax) ymax = h->GetMaximum(); }
  ymin = std::max(1e-12, 0.1 * ymax * 1e-3); // heuristic floor
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax * 1.5);

  frame->Draw("axis");
  for (auto* h : H) h->Draw("hist same");

  // Legend
  TLegend* L = new TLegend(0.52, 0.70, 0.88, 0.90);
  L->SetBorderSize(0);
  L->SetFillStyle(0);
  L->SetTextFont(42);
  L->SetTextSize(0.035);
  L->SetHeader("Angular bands","C");
  for (size_t k=0;k<bands.size();++k)
    L->AddEntry(H[k], bands[k].label, "l");
  L->Draw();

  c->SaveAs("dose_r2_vs_radius_overlay.png");
  c->SaveAs("dose_r2_vs_radius_overlay.pdf");

  f->Close();
}

