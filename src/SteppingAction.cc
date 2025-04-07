#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Gamma.hh"
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
				if (info) std::cout << info->GetScatterOrder() << '\n';
                if (info) info->IncrementScatterOrder();
				
            }
		}
	}
	const auto secondaries = step->GetSecondaryInCurrentStep();
    if (secondaries && !secondaries->empty()) {
        MyTrackInfo* parentInfo = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());

        for (auto secTrack : *secondaries) {
            const G4ParticleDefinition* secParticle = secTrack->GetParticleDefinition();

            if (secParticle == G4Gamma::GammaDefinition())
                secTrack->SetUserInformation(new MyTrackInfo(parentInfo ? parentInfo->GetScatterOrder() : 1));
        }
    }

	// kill the primary photons when they leave the central voxel
	auto parentID = track->GetParentID();
	if(parentID == 0)
	{
		// Kill photon if it leaves the central voxel
		G4ThreeVector position = step->GetPostStepPoint()->GetPosition();

		G4double radius = 0.5 * mm; // Size of each voxel

		G4double distance_square = position.x()* position.x() +
							position.y() * position.y() +
							position.y() * position.y();
							 
		if(distance_square < (radius * radius))	
		{
		step->GetTrack()->SetTrackStatus(fStopAndKill);
		return;
		}
	}

	G4double edep = step->GetTotalEnergyDeposit();
	if(edep == 0.) return;

	G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
	G4double x = position.x(),
			 y = position.y(),
			 z = position.z();

	G4double r = std::sqrt(x*x + y*y + z*z);
	G4double theta = (r > 0) ? std::acos(z/r): 0.0;

	// Define binning parameters
    G4int r_bins = 100;
    G4double r_max = 200.0 * mm;
    G4double theta_max = CLHEP::pi;

    G4double dr = r_max / r_bins;
    G4double dtheta = theta_max / r_bins;

    // Compute volume element (in spherical coordinates)
    G4double r_inner = (r - dr / 2);
    G4double r_outer = (r + dr / 2);
    G4double volume_element = (1.0 / 3.0) * (r_outer*r_outer*r_outer - r_inner*r_inner*r_inner) * dtheta * (2 * CLHEP::pi);

    // Get material density (assuming water-like tissue for now)
    G4Material* material = step->GetTrack()->GetMaterial();
    G4double density = material->GetDensity();  // in kg/m³

    // Compute mass of the volume element
    G4double mass = density * volume_element;  // mass in kg

    // Convert edep (MeV) to dose (Gy)
    G4double edep_J = edep * 1.60218e-13;  // Convert MeV to Joules
    G4double dose = (mass > 0) ? edep_J / mass : 0.0;  // Dose in Gy

	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillNtupleDColumn(0, r);
	analysisManager->FillNtupleDColumn(1, theta);
	analysisManager->FillNtupleDColumn(2, dose);
	analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
