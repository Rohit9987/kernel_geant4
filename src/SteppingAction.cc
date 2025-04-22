#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "MyTrackInfo.hh"
#include "Run.hh"

#define DEBUG 0
#if DEBUG
  #define debug(x) std::cout << "SA---->" << x << std::endl
#else
  #define debug(x)
#endif

namespace B4c {

SteppingAction::SteppingAction() {}
SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    auto track = step->GetTrack();
	
    bool isCompton = false;
    bool isPhoto   = false;
	bool isEBrem = false;

    // Check interaction process if it's a gamma
	
    G4StepPoint* postStep = step->GetPostStepPoint();
    const G4VProcess* process = postStep->GetProcessDefinedStep();
    G4String procName = process->GetProcessName();

	if (process) {
		MyTrackInfo* info = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());

		if (procName == "compt" || procName == "phot") {

			// Record primary interaction position once
			if (info->GetScatterOrder() == 0 && !info->IsPrimaryInteractionSet()) {
				info->SetPrimaryInteractionPosition(postStep->GetPosition());
				MyRun* run = static_cast<MyRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
				run->AddPhoton();
			}

			if (procName == "compt") {
				info->IncrementScatterOrder();
				isCompton = true;
			}
			else if (procName == "phot") {
					info->SetGamma(1);
					isPhoto = true;
			}
		}
		if (procName == "eBrem")
		{
			isEBrem = true;
		}
    }

    // Track secondaries
    const auto secondaries = step->GetSecondaryInCurrentStep();
    if (secondaries && !secondaries->empty()) {

        MyTrackInfo* parentInfo = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());

        for (auto secTrack : *secondaries) {
            const auto* secDef = secTrack->GetParticleDefinition();
            if (secDef == G4Gamma::GammaDefinition() || 
                secDef == G4Electron::ElectronDefinition() || 
                secDef == G4Positron::PositronDefinition()) {

                auto* secInfo = new MyTrackInfo(parentInfo ? parentInfo->GetScatterOrder() : -1);

                if (parentInfo && parentInfo->IsPrimaryInteractionSet()) {
                    secInfo->SetPrimaryInteractionPosition(parentInfo->GetPrimaryInteractionPosition());
					if (parentInfo->GetEBrem())
						isEBrem = 1;
                }

				if (isEBrem)   secInfo->SetEBrem(isEBrem);
				else if (isPhoto)   secInfo->SetPhotoElectron(isPhoto);
                else if (isCompton) secInfo->SetComptonElectron(isCompton);

                secTrack->SetUserInformation(secInfo);

            }
        }
    }

	
    // Handle energy deposition
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep == 0.0) {
        return;
    }

    auto* info = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());
    if (!info) return;

    G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector ppos = info->GetPrimaryInteractionPosition();

    G4double dx = pos.x() - ppos.x();
    G4double dy = pos.y() - ppos.y();
    G4double dz = pos.z() - ppos.z();

    G4double r = std::sqrt(dx*dx + dy*dy + dz*dz);
    G4double theta = (r > 0) ? std::acos(dz / r) : 0.0;
	if((dy < 0) ^ (dx < 0))
		theta *= -1;
	
    int type = 9;
    if (info->GetGammaPhoton() || track->GetParticleDefinition() == G4Gamma::GammaDefinition()) {
        type = 0;
    } else if (info->GetComptonElectron()) {
        type = 1;
    } else if (info->GetPhotoElectron()) {
        type = 2;
    } else if (info->GetEBrem())
	{
		type = 3;
	}
	else {
		type = 9;
    }

    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(0, r);
    analysisManager->FillNtupleDColumn(1, theta);
    analysisManager->FillNtupleDColumn(2, edep);
    analysisManager->FillNtupleIColumn(3, info->GetScatterOrder());
    analysisManager->FillNtupleIColumn(4, type);
    analysisManager->AddNtupleRow();
}

}  // namespace B4c

