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

#define DEBUG 0  // Set to 0 to disable debug output

#if DEBUG
  #define debug(x) std::cout << x << std::endl
#else
  #define debug(x)
#endif

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
	debug("Particle: " << track->GetParticleDefinition()->GetParticleName());

	bool isCompton = false;
	bool isPhoto = false;

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
					debug("Photon scatter" << info->GetScatterOrder() << ", Process:" << procName);
					// if its the first interaction, save this position

					if(info->GetScatterOrder() == 0 && !info->IsPrimaryInteractionSet()) 
					{
                        info->SetPrimaryInteractionPosition(postStepPoint->GetPosition());

						// update run-level counter (thread-local)
						MyRun* run = static_cast<MyRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
						run->AddPhoton();
					}
					if(procName == "compt")					// only for compton increment the scatter
					{
						info->IncrementScatterOrder();
						isCompton = true;
					}
					if(procName == "phot")
						isPhoto = true;

					debug("Photon scatter incremented: " << info->GetScatterOrder() << ", Process:" << procName);
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

                MyTrackInfo* secTrackInfo = new MyTrackInfo(parentInfo ? parentInfo->GetScatterOrder() : -1);

                // Copy primary interaction position from parent
                if(parentInfo && parentInfo->IsPrimaryInteractionSet()) {
                    secTrackInfo->SetPrimaryInteractionPosition(parentInfo->GetPrimaryInteractionPosition());
                }

				if(secParticle == G4Electron::ElectronDefinition() || 
					secParticle == G4Positron::PositronDefinition())
				{
					if(isCompton)
						secTrackInfo->SetComptonElectron(1);
					if(isPhoto)
						secTrackInfo->SetPhotoElectron(1);
				}
				else if(secParticle == G4Gamma::GammaDefinition())
					secTrackInfo->SetGamma(1);


                secTrack->SetUserInformation(secTrackInfo);
				
				debug("    SecondaryName:" << secParticle->GetParticleName());
				debug("    ScatterOrder" << secTrackInfo->GetScatterOrder());
            }
        }
    }


	G4double edep = step->GetTotalEnergyDeposit();
	if(edep == 0.) 
	{
		debug("    No energy deposition");
		return;
	}

	MyTrackInfo* info = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());

	G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
	G4ThreeVector primarypos = info->GetPrimaryInteractionPosition();
	G4double x = position.x() - primarypos.x(),
			 y = position.y() - primarypos.y(),
			 z = position.z() - primarypos.z();


	G4double r = std::sqrt(x*x + y*y + z*z);
	G4double theta = (r > 0) ? std::acos(z/r): 0.0;

	debug("    Position: (" << position.x() << ',' << position.y() << ',' << position.z() << ")");
	debug("    Primary Position: (" << primarypos.x() << ',' << primarypos.y() << ',' << primarypos.z() << ")");
	debug("    Edep: " << edep);
	debug("    r: " << r );
	debug("    theta: " << theta);
	debug("    ScatterOrder: " << info->GetScatterOrder());
	debug("    Photo: " << info->GetPhotoElectron() );
	debug("    Compton: " << info->GetComptonElectron());
	debug("	   Gamma: " << info->GetGammaPhoton());

	if(info->GetComptonElectron() == info->GetPhotoElectron())
	{
		debug("Is Compton: " << info->GetComptonElectron() << '\n'
				  << "Is photo: " << info->GetPhotoElectron() << '\n');
	}

	int type;
	if(info->GetGammaPhoton())
		type = 0;
	else if(info->GetComptonElectron())
		type = 1;
	else if(info->GetPhotoElectron())
		type = 2;
	else 
		type = 9;
	

	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillNtupleDColumn(0, r);
	analysisManager->FillNtupleDColumn(1, theta);
	analysisManager->FillNtupleDColumn(2, edep);
	analysisManager->FillNtupleIColumn(3, info->GetScatterOrder());
	analysisManager->FillNtupleIColumn(4, type); 

	analysisManager->AddNtupleRow();
}


}
