// plotXZ_fast.C
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <cmath>
#include <vector>

// Helper to build 2D hist with common settings
static TH2D* MakeHist(const char* name, const char* label,
                      int nBins, double range) {
  auto* h = new TH2D(name,
      Form("XZ Dose Projection (%s);X (mm);Z (mm)", label),
      nBins, -range, range, nBins, -range, range);
  h->SetDirectory(nullptr);   // faster; don't attach to a file
  // h->Sumw2(false);         // uncomment if you don't need bin errors
  return h;
}

void plot(const char* fname = "../physics/DoseKernel_0.5MeV.root",
                 double range = 60.0,   // mm
                 double binWidth = 0.5, // mm
                 double yTolerance = 0.25) {

  TFile* file = TFile::Open(fname, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fname << "\n";
    return;
  }

  TTree* tree = static_cast<TTree*>(file->Get("DoseData"));
  if (!tree) {
    std::cerr << "TTree 'DoseData' not found!\n";
    file->Close();
    return;
  }

  // Speed: enable cache & disable unused branches
  tree->SetCacheSize(128*1024*1024);
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("X", 1);
  tree->SetBranchStatus("Y", 1);
  tree->SetBranchStatus("Z", 1);
  tree->SetBranchStatus("edep", 1);
  tree->SetBranchStatus("Type", 1);
  tree->SetBranchStatus("Scatter", 1);
  tree->AddBranchToCache("X", kTRUE);
  tree->AddBranchToCache("Y", kTRUE);
  tree->AddBranchToCache("Z", kTRUE);
  tree->AddBranchToCache("edep", kTRUE);
  tree->AddBranchToCache("Type", kTRUE);
  tree->AddBranchToCache("Scatter", kTRUE);

  double x=0, y=0, z=0, edep=0;
  int type=9, scatter=-1;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Z", &z);
  tree->SetBranchAddress("edep", &edep);
  tree->SetBranchAddress("Type", &type);
  tree->SetBranchAddress("Scatter", &scatter);

  const int nBins = static_cast<int>(std::floor((2.0*range)/binWidth + 0.5));

  // Histograms: one pass fills all of these
  TH2D* hAll      = MakeHist("hXZ_All",      "All",      nBins, range);
  TH2D* hGamma    = MakeHist("hXZ_Gamma",    "Gamma",    nBins, range);   // Type==0
  TH2D* hComp     = MakeHist("hXZ_Compton",  "Compton",  nBins, range);   // Type==1
  TH2D* hPhot     = MakeHist("hXZ_Photo",    "Photoelectric", nBins, range); // Type==2
  TH2D* hBrem     = MakeHist("hXZ_Brem",     "BremPhoton", nBins, range); // Type==3
  TH2D* hAnni     = MakeHist("hXZ_Anni",     "AnnihilPhoton", nBins, range); // Type==4
  TH2D* hPair     = MakeHist("hXZ_Pair",     "PairProduction", nBins, range); // Type==5
  TH2D* hUnknown  = MakeHist("hXZ_Unknown",  "Unknown",  nBins, range);   // Type==9

  // Optional: Compton scatter detail
  TH2D* hCompS1   = MakeHist("hXZ_Compton_S1",   "Compton (Scatter=1)", nBins, range);
  TH2D* hCompSgt1 = MakeHist("hXZ_Compton_Sgt1", "Compton (Scatter>1)", nBins, range);

  const Long64_t N = tree->GetEntries();
  const Long64_t printEvery = std::max<Long64_t>(N/50, 100000); // ~2% or 100k
  for (Long64_t i=0; i<N; ++i) {
    tree->GetEntry(i);
    if ((i % printEvery) == 0) {
      std::cout << "Processed " << i << " / " << N << '\r' << std::flush;
    }

    // Mid-plane slice
    if (std::abs(y) > yTolerance) continue;

    // Fill "All"
    hAll->Fill(x, z, edep);

    // Fill by InteractionType (int code)
    switch (type) {
      case 0: hGamma->Fill(x, z, edep); break;
      case 1: // Compton
        hComp->Fill(x, z, edep);
        if (scatter == 1)      hCompS1->Fill(x, z, edep);
        else if (scatter > 1)  hCompSgt1->Fill(x, z, edep);
        break;
      case 2: hPhot->Fill(x, z, edep); break;
      case 3: hBrem->Fill(x, z, edep); break;
      case 4: hAnni->Fill(x, z, edep); break;
      case 5: hPair->Fill(x, z, edep); break;
      default: hUnknown->Fill(x, z, edep); break; // includes 9
    }
  }
  std::cout << "\nDone. Drawing...\n";

  gStyle->SetOptStat(0);

  auto drawSave = [](TH2D* h, const char* out){
    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), h->GetTitle(), 900, 800);
    c->SetLogz();
    h->Draw("COLZ");
    c->SaveAs(out);
    delete c;
  };

  // Save images (PNG). Add/remove as you like.
  drawSave(hAll,     "DoseXZ_y0_All.png");
  drawSave(hGamma,   "DoseXZ_y0_Gamma.png");
  drawSave(hComp,    "DoseXZ_y0_Compton.png");
  drawSave(hCompS1,  "DoseXZ_y0_Compton_S1.png");
  drawSave(hCompSgt1,"DoseXZ_y0_Compton_Sgt1.png");
  drawSave(hPhot,    "DoseXZ_y0_Photoelectric.png");
  drawSave(hBrem,    "DoseXZ_y0_BremPhoton.png");
  drawSave(hAnni,    "DoseXZ_y0_AnnihilPhoton.png");
  drawSave(hPair,    "DoseXZ_y0_PairProduction.png");
  drawSave(hUnknown, "DoseXZ_y0_Unknown.png");

  file->Close();
}

