#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
	G4RunManager::GetRunManager()->SetPrintProgress(1);

	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetNtupleMerging(true);

	analysisManager->CreateNtuple("DoseData", "Voxel Dose Data");
	analysisManager->CreateNtupleDColumn("R");
	analysisManager->CreateNtupleDColumn("Theta");
	analysisManager->CreateNtupleDColumn("edep");
	analysisManager->CreateNtupleIColumn("Scatter");
	analysisManager->FinishNtuple();
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile("voxel_dose.root");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->Write();
        analysisManager->CloseFile();  // print histogram statistics
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
