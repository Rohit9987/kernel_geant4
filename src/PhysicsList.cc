#include "PhysicsList.hh"

#include "PhysicsList.hh"

#include "G4EmPenelopePhysics.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

namespace B4
{
PhysicsList::PhysicsList() {
    // Register Penelope EM physics (handles gamma, e-, e+)
    RegisterPhysics(new G4EmPenelopePhysics());

    // Enable detailed atomic effects
    auto emParams = G4EmParameters::Instance();
    emParams->SetFluo(true);     // Enable fluorescence
    emParams->SetAuger(true);    // Enable Auger electron emission

    G4cout << ">> Using G4EmPenelopePhysics with fluorescence and Auger enabled" << G4endl;
}

PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
}


}

