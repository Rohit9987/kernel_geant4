/*
#include "PhysicsList.hh"


#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option3.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

namespace B4
{
PhysicsList::PhysicsList() {

  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  RegisterPhysics(new G4EmStandardPhysics_option3(1));
}

PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts()
{

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
}

}
*/
