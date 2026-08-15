#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

#include "KernelBinning.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

namespace {

G4double StandardErrorOfMean(G4double sum,
                             G4double sumSquares,
                             G4int histories)
{
  if (histories < 2) {
    return std::numeric_limits<G4double>::quiet_NaN();
  }

  const G4double n = static_cast<G4double>(histories);
  const G4double centeredSumSquares =
    std::max(0.0, sumSquares - sum * sum / n);
  return std::sqrt(centeredSumSquares / (n * (n - 1.0)));
}

} // namespace

namespace B4
{

RunAction::RunAction(G4String energyStr): G4UserRunAction(), fEnergyStr(energyStr)
{
	G4RunManager::GetRunManager()->SetPrintProgress(1);

	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetVerboseLevel(1);

	fKernelNtupleId = analysisManager->CreateNtuple(
	  "Kernel", "Radial-angular energy-deposition kernel");
	analysisManager->CreateNtupleIColumn(fKernelNtupleId, "RadialBin");
	analysisManager->CreateNtupleIColumn(fKernelNtupleId, "ThetaBin");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId, "RadialLower_mm");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId, "RadialUpper_mm");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId, "ThetaLower_rad");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId, "ThetaUpper_rad");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId, "EdepSum_MeV");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId, "EdepSum2_MeV2");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId,
	                                    "MeanEdep_MeV_per_history");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId,
	                                    "StdError_MeV_per_history");
	analysisManager->CreateNtupleDColumn(fKernelNtupleId,
	                                    "RelativeStdError");
	analysisManager->FinishNtuple(fKernelNtupleId);

	fSummaryNtupleId = analysisManager->CreateNtuple(
	  "RunSummary", "Kernel scoring summary");
	analysisManager->CreateNtupleIColumn(fSummaryNtupleId, "Histories");
	analysisManager->CreateNtupleIColumn(fSummaryNtupleId, "PrimaryPhotons");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId, "ScoredEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId, "UnbinnedEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "UnbinnedEnergyFraction");
	analysisManager->FinishNtuple(fSummaryNtupleId);
}

RunAction::~RunAction()
{

}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	if (G4Threading::IsWorkerThread()) return;

	auto analysisManager = G4AnalysisManager::Instance();
	G4String filename = "physics/DoseKernel_" + fEnergyStr + ".root";
	analysisManager->OpenFile(filename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
	if (G4Threading::IsWorkerThread()) return;

	const auto* myRun = static_cast<const MyRun*>(run);
	const G4int histories = run->GetNumberOfEvent();
	const G4double n = static_cast<G4double>(histories);
	auto analysisManager = G4AnalysisManager::Instance();

	for (std::size_t radial = 0;
	     radial < B4c::KernelBinning::NumRadialBins();
	     ++radial) {
	  for (std::size_t theta = 0;
	       theta < B4c::KernelBinning::kNumThetaBins;
	       ++theta) {
	    const auto linear = B4c::KernelBinning::LinearIndex(radial, theta);
	    const G4double sum = myRun->GetEdepSum(linear);
	    const G4double sumSquares = myRun->GetEdepSumSquares(linear);
	    const G4double mean = histories > 0 ? sum / n : 0.0;
	    const G4double standardError =
	      StandardErrorOfMean(sum, sumSquares, histories);
	    const G4double relativeError = mean > 0.0
	      ? standardError / mean
	      : std::numeric_limits<G4double>::quiet_NaN();

	    analysisManager->FillNtupleIColumn(
	      fKernelNtupleId, 0, static_cast<G4int>(radial));
	    analysisManager->FillNtupleIColumn(
	      fKernelNtupleId, 1, static_cast<G4int>(theta));
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 2, B4c::KernelBinning::RadialLowerMm(radial));
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 3, B4c::KernelBinning::RadialUpperMm(radial));
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 4, B4c::KernelBinning::ThetaLowerRad(theta));
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 5, B4c::KernelBinning::ThetaUpperRad(theta));
	    analysisManager->FillNtupleDColumn(fKernelNtupleId, 6, sum / MeV);
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 7, sumSquares / (MeV * MeV));
	    analysisManager->FillNtupleDColumn(fKernelNtupleId, 8, mean / MeV);
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 9, standardError / MeV);
	    analysisManager->FillNtupleDColumn(
	      fKernelNtupleId, 10, relativeError);
	    analysisManager->AddNtupleRow(fKernelNtupleId);
	  }
	}

	const G4double scored = myRun->GetTotalScoredEdep();
	const G4double unbinned = myRun->GetUnbinnedEdepSum();
	const G4double accounted = scored + unbinned;
	const G4double unbinnedFraction =
	  accounted > 0.0 ? unbinned / accounted : 0.0;

	analysisManager->FillNtupleIColumn(
	  fSummaryNtupleId, 0, histories);
	analysisManager->FillNtupleIColumn(
	  fSummaryNtupleId, 1, myRun->GetPhotonCount());
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, 2, scored / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, 3, unbinned / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, 4, unbinnedFraction);
	analysisManager->AddNtupleRow(fSummaryNtupleId);

	analysisManager->Write();
	analysisManager->CloseFile();

	G4cout << "Kernel output: " << B4c::KernelBinning::NumBins()
	       << " bins; unbinned energy fraction = "
	       << unbinnedFraction << G4endl;

	std::ofstream outFile(
	  "./photon_counts/physics_nphotons_" + fEnergyStr + ".txt");
	if (outFile.is_open()) {
	  outFile << myRun->GetPhotonCount() << std::endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
