// plot_edep_vs_radius_angles.C
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

struct RaySpec { double thetaDeg; const char* label; int color; int style; };

void plot_radial(const char* fname = "../physics/DoseKernel_6.0MeV.root",
                                double yTol_mm = 0.25,     // mid-plane slice thickness
                                double rmax_gcm2 = 20.0,   // radial axis extent (g cm^-2)
                                double binw_gcm2 = 0.1,    // radial bin width (g cm^-2)
                                double halfWidthDeg = 1.0, // half-width around each ray angle
                                bool   normalizeByAngle = false) // divide by 2*halfWidthDeg
{
  // ---------- open tree (fast I/O) ----------
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

  // ---------- rays to sample (center angles & styles) ----------
  std::vector<RaySpec> rays = {
    {  0.000,  "0^{#circ}",                kBlack, 1 },
    { 28.815,  "28.815^{#circ}",           kBlack, 2 }, // dashed
    { 47.515,  "47.515^{#circ}",           kBlack, 3 }  // dotted
  };

  const double mm_per_gcm2 = 10.0;     // water
  const int nbins = std::max(1, int(std::floor(rmax_gcm2/binw_gcm2 + 0.5)));

  // one histogram per ray
  std::vector<TH1D*> H; H.reserve(rays.size());
  for (size_t i=0;i<rays.size();++i) {
    auto* h = new TH1D(Form("h_edep_ray_%zu", i), "", nbins, 0.0, rmax_gcm2);
    h->SetDirectory(nullptr);
    h->Sumw2();
    h->SetLineColor(rays[i].color);
    h->SetLineStyle(rays[i].style);
    h->SetLineWidth(2);
    H.push_back(h);
  }

  // ---------- single pass over entries ----------
  const Long64_t N = T->GetEntries();
  const Long64_t printEvery = std::max<Long64_t>(N/50, 100000);
  for (Long64_t i=0;i<N;++i) {
    T->GetEntry(i);
    if ((i % printEvery) == 0) std::cout << "Processed " << i << "/" << N << "\r" << std::flush;

    if (std::abs(y) > yTol_mm) continue; // mid-plane
    if (z <= 0.0) continue;              // forward hemisphere

    const double r_mm   = std::hypot(x, z);
    const double r_gcm2 = r_mm / mm_per_gcm2;
    if (r_gcm2 <= 0.0 || r_gcm2 > rmax_gcm2) continue;

    // angle from the +z axis in XZ plane; fold L/R symmetry via |x|
    const double thetaDeg = std::atan2(std::abs(x), z) * 180.0 / TMath::Pi();

    // weight = edep (no r^2 weighting here)
    const double w = edep;

    for (size_t k=0;k<rays.size();++k) {
      if (std::abs(thetaDeg - rays[k].thetaDeg) <= halfWidthDeg) {
        H[k]->Fill(r_gcm2, w);
      }
    }
  }
  std::cout << "\nDone.\n";

  if (normalizeByAngle) {
    const double width = std::max(1e-9, 2.0*halfWidthDeg);
    for (auto* h : H) h->Scale(1.0/width);
  }

  // ---------- draw ----------
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_edep_rays", "edep vs radius along rays", 900, 800);
  c->SetLogy();
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  // frame for axes
  TH1D* frame = (TH1D*)H[0]->Clone("frame_edep_rays");
  frame->Reset();
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("Radius (g cm^{-2})");
  frame->GetYaxis()->SetTitle(normalizeByAngle ? "edep per degree (arb.)" : "edep (arb.)");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.040);

  double ymax = 0.0;
  for (auto* h : H) ymax = std::max(ymax, h->GetMaximum());
  frame->SetMinimum(std::max(1e-12, 0.1*ymax*1e-3));
  frame->SetMaximum(ymax*1.5);

  frame->Draw("axis");
  for (auto* h : H) h->Draw("hist same");

  // legend
  TLegend* L = new TLegend(0.52, 0.70, 0.88, 0.90);
  L->SetBorderSize(0);
  L->SetFillStyle(0);
  L->SetTextFont(42);
  L->SetTextSize(0.035);
  L->SetHeader(Form("|#theta - #theta_{0}| #leq %.2f^{#circ}", halfWidthDeg),"C");
  for (size_t k=0;k<rays.size();++k)
    L->AddEntry(H[k], rays[k].label, "l");
  L->Draw();

  c->SaveAs("edep_vs_radius_rays.png");
  c->SaveAs("edep_vs_radius_rays.pdf");

  f->Close();
}

