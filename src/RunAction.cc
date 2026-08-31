#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

#include "KernelBinning.hh"
#include "KernelScoring.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

namespace {

enum SummaryColumn : G4int {
  kHistories = 0,
  kPrimaryPhotons,
  kPrimaryEnergy,
  kScoredEdep,
  kLocalAtOriginEdep,
  kOutsideKernelRadiusEdep,
  kMissingKernelFrameEdep,
  kInvalidDirectionOrAngleEdep,
  kUnbinnedEdep,
  kEscapedWorldEnergy,
  kAccountedEnergy,
  kUnaccountedEnergy,
  kKernelScoredFractionOfPrimary,
  kUnbinnedFractionOfPrimary,
  kEscapedFractionOfPrimary,
  kEnergyClosureFraction
};

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
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId, "PrimaryEnergy_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId, "ScoredEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "LocalAtOriginEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "OutsideKernelRadiusEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "MissingKernelFrameEdep_MeV");
	analysisManager->CreateNtupleDColumn(
	  fSummaryNtupleId, "InvalidDirectionOrAngleEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId, "UnbinnedEdep_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "EscapedWorldEnergy_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "AccountedEnergy_MeV");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "UnaccountedEnergy_MeV");
	analysisManager->CreateNtupleDColumn(
	  fSummaryNtupleId, "KernelScoredFractionOfPrimary");
	analysisManager->CreateNtupleDColumn(
	  fSummaryNtupleId, "UnbinnedFractionOfPrimary");
	analysisManager->CreateNtupleDColumn(
	  fSummaryNtupleId, "EscapedFractionOfPrimary");
	analysisManager->CreateNtupleDColumn(fSummaryNtupleId,
	                                    "EnergyClosureFraction");
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

	const G4double primaryEnergy = myRun->GetPrimaryEnergySum();
	const G4double scored = myRun->GetTotalScoredEdep();
	const G4double localAtOrigin = myRun->GetUnbinnedEdepSum(
	  B4c::UnbinnedReason::LocalAtOrigin);
	const G4double outsideKernelRadius = myRun->GetUnbinnedEdepSum(
	  B4c::UnbinnedReason::OutsideKernelRadius);
	const G4double missingKernelFrame = myRun->GetUnbinnedEdepSum(
	  B4c::UnbinnedReason::MissingKernelFrame);
	const G4double invalidDirectionOrAngle = myRun->GetUnbinnedEdepSum(
	  B4c::UnbinnedReason::InvalidDirectionOrAngle);
	const G4double unbinned = myRun->GetTotalUnbinnedEdepSum();
	const G4double escaped = myRun->GetEscapedEnergySum();
	const G4double accounted = scored + unbinned + escaped;
	const G4double unaccounted = primaryEnergy - accounted;
	const G4double scoredFraction =
	  primaryEnergy > 0.0 ? scored / primaryEnergy : 0.0;
	const G4double unbinnedFraction =
	  primaryEnergy > 0.0 ? unbinned / primaryEnergy : 0.0;
	const G4double escapedFraction =
	  primaryEnergy > 0.0 ? escaped / primaryEnergy : 0.0;
	const G4double closureFraction =
	  primaryEnergy > 0.0 ? accounted / primaryEnergy : 0.0;

	analysisManager->FillNtupleIColumn(
	  fSummaryNtupleId, kHistories, histories);
	analysisManager->FillNtupleIColumn(
	  fSummaryNtupleId, kPrimaryPhotons, myRun->GetPhotonCount());
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kPrimaryEnergy, primaryEnergy / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kScoredEdep, scored / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kLocalAtOriginEdep, localAtOrigin / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kOutsideKernelRadiusEdep,
	  outsideKernelRadius / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kMissingKernelFrameEdep,
	  missingKernelFrame / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kInvalidDirectionOrAngleEdep,
	  invalidDirectionOrAngle / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kUnbinnedEdep, unbinned / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kEscapedWorldEnergy, escaped / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kAccountedEnergy, accounted / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kUnaccountedEnergy, unaccounted / MeV);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kKernelScoredFractionOfPrimary, scoredFraction);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kUnbinnedFractionOfPrimary, unbinnedFraction);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kEscapedFractionOfPrimary, escapedFraction);
	analysisManager->FillNtupleDColumn(
	  fSummaryNtupleId, kEnergyClosureFraction, closureFraction);
	analysisManager->AddNtupleRow(fSummaryNtupleId);

	analysisManager->Write();
	analysisManager->CloseFile();

	G4cout << G4endl
	       << "Run energy accounting:" << G4endl
	       << "  Histories                         = " << histories << G4endl
	       << "  PrimaryPhotons                    = "
	       << myRun->GetPhotonCount() << G4endl
	       << "  PrimaryEnergy_MeV                 = "
	       << primaryEnergy / MeV << G4endl
	       << "  ScoredEdep_MeV                    = "
	       << scored / MeV << G4endl
	       << "  LocalAtOriginEdep_MeV             = "
	       << localAtOrigin / MeV << G4endl
	       << "  OutsideKernelRadiusEdep_MeV       = "
	       << outsideKernelRadius / MeV << G4endl
	       << "  MissingKernelFrameEdep_MeV        = "
	       << missingKernelFrame / MeV << G4endl
	       << "  InvalidDirectionOrAngleEdep_MeV   = "
	       << invalidDirectionOrAngle / MeV << G4endl
	       << "  UnbinnedEdep_MeV                  = "
	       << unbinned / MeV << G4endl
	       << "  EscapedWorldEnergy_MeV            = "
	       << escaped / MeV << G4endl
	       << "  AccountedEnergy_MeV               = "
	       << accounted / MeV << G4endl
	       << "  UnaccountedEnergy_MeV             = "
	       << unaccounted / MeV << G4endl
	       << "  KernelScoredFractionOfPrimary     = "
	       << scoredFraction << G4endl
	       << "  UnbinnedFractionOfPrimary         = "
	       << unbinnedFraction << G4endl
	       << "  EscapedFractionOfPrimary          = "
	       << escapedFraction << G4endl
	       << "  EnergyClosureFraction             = "
	       << closureFraction << G4endl
	       << "Kernel output: " << B4c::KernelBinning::NumBins()
	       << " bins" << G4endl;

	std::ofstream outFile(
	  "./photon_counts/physics_nphotons_" + fEnergyStr + ".txt");
	if (outFile.is_open()) {
	  outFile << myRun->GetPhotonCount() << std::endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
