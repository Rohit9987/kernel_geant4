#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "MyTrackInfo.hh"

namespace B4c
{

SteppingAction::SteppingAction()
{
}

SteppingAction::~SteppingAction()
{
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
	auto track = step->GetTrack();

	if(track->GetParticleDefinition() == G4Gamma::GammaDefinition())
	{
		G4StepPoint* postStepPoint = step->GetPostStepPoint();
		const G4VProcess* process = postStepPoint->GetProcessDefinedStep();

		if(process)
		{
			G4String procName = process->GetProcessName();
			if (procName == "compt" || procName == "phot") {
                MyTrackInfo* info = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());
                if (info) 
				{
					// if its the first interaction, save this position
					if(info->GetScatterOrder() == 0 && !info->IsPrimaryInteractionSet()) {
                        info->SetPrimaryInteractionPosition(postStepPoint->GetPosition());
					}
                    info->IncrementScatterOrder();
				}
            }
		}
	}

	const auto secondaries = step->GetSecondaryInCurrentStep();
    if (secondaries && !secondaries->empty()) {
        MyTrackInfo* parentInfo = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());

        for (auto secTrack : *secondaries) {
            const G4ParticleDefinition* secParticle = secTrack->GetParticleDefinition();
				
			if (secParticle == G4Gamma::GammaDefinition() ||
                secParticle == G4Electron::ElectronDefinition() ||
                secParticle == G4Positron::PositronDefinition()) {

                MyTrackInfo* secTrackInfo = new MyTrackInfo(parentInfo ? parentInfo->GetScatterOrder() : 1);

                // Copy primary interaction position from parent
                if(parentInfo && parentInfo->IsPrimaryInteractionSet()) {
                    secTrackInfo->SetPrimaryInteractionPosition(parentInfo->GetPrimaryInteractionPosition());
                }

                secTrack->SetUserInformation(secTrackInfo);
            }
        }
    }


	G4double edep = step->GetTotalEnergyDeposit();
	if(edep == 0.) return;

	MyTrackInfo* info = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());

	G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
	G4ThreeVector primarypos = info->GetPrimaryInteractionPosition();
	G4double x = position.x() - primarypos.x(),
			 y = position.y() - primarypos.y(),
			 z = position.z() - primarypos.z();


	G4double r = std::sqrt(x*x + y*y + z*z);
	G4double theta = (r > 0) ? std::acos(z/r): 0.0;

	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillNtupleDColumn(0, r);
	analysisManager->FillNtupleDColumn(1, theta);
	analysisManager->FillNtupleDColumn(2, edep);
	analysisManager->FillNtupleIColumn(3, info->GetScatterOrder());
	analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
