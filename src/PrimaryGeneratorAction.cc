
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "MyPrimaryParticleInfo.hh"

#define DEBUG 0
#if DEBUG
#define debug(x) std::cout  << "PGA--------->"<< x << std::endl
#else
#define debug(x)
#endif

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String energyStr)
	: G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
      
  G4double energy = G4UIcommand::ConvertToDimensionedDouble(energyStr);
  fParticleGun->SetParticleEnergy(energy*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // Set gun position
  fParticleGun
    ->SetParticlePosition(G4ThreeVector(0., 0., 0.));

  fParticleGun->GeneratePrimaryVertex(anEvent);
  debug("Primary generated");

  G4PrimaryVertex* vertex = anEvent->GetPrimaryVertex();
  for(G4int i =0; i < vertex->GetNumberOfParticle(); i++)
  {
	  G4PrimaryParticle* particle = vertex->GetPrimary(i);
	  particle->SetUserInformation(new MyPrimaryParticleInfo(0));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
