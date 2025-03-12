#pragma once

#include "G4VSensitiveDetector.hh"
#include <map>
#include "G4SystemOfUnits.hh"

class VoxelSensitiveDetector : public G4VSensitiveDetector {
public:
    VoxelSensitiveDetector(const G4String& name) : G4VSensitiveDetector(name) {
        collectionName.insert("DoseDepositCollection");
    }

    void Initialize(G4HCofThisEvent*) override { doseMap.clear();}

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
        G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0) {
            G4int copyNo = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
            doseMap[copyNo] += edep;
        }
        return true;
    }

    void EndOfEvent(G4HCofThisEvent*) override {
    }

	std::map<G4int, G4double> GetDoseMap()
	{
		return doseMap;
	}
private:
    std::map<G4int, G4double> doseMap;
};

