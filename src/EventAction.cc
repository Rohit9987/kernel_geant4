#include "EventAction.hh"
#include "VoxelSensitiveDetector.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

namespace B4c
{

EventAction::EventAction()
{}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}


void EventAction::EndOfEventAction(const G4Event* event)
{
	auto analysisManager = G4AnalysisManager::Instance();

	auto sdManager = G4SDManager::GetSDMpointer();
	auto voxelSD = static_cast<VoxelSensitiveDetector*>(sdManager->FindSensitiveDetector("VoxelSD"));
	if (!voxelSD) return;

	G4double halfVoxelRes = 0.5*mm;
	G4int numVoxels = 101;

	const auto& doseMap = voxelSD->GetDoseMap();
	for (const auto& entry : doseMap) {
		G4int copyNo = entry.first;
		G4double dose = entry.second / gray;

		G4int iZ = copyNo % numVoxels;
		iZ = iZ - (numVoxels-1) * halfVoxelRes; 

		G4int iY = (copyNo / numVoxels) % numVoxels;
		iY = iY - (numVoxels-1) * halfVoxelRes; 

		G4int iX = copyNo / (numVoxels * numVoxels);
		iX = iX - (numVoxels-1) * halfVoxelRes; 

		
		analysisManager->FillNtupleIColumn(0, iX);
		analysisManager->FillNtupleIColumn(1, iY);
		analysisManager->FillNtupleIColumn(2, iZ);
		analysisManager->FillNtupleDColumn(3, dose);
		analysisManager->AddNtupleRow();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
