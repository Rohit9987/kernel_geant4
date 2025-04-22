#pragma once

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

// Choose physics model here:
// 0 = Standard Option 3
// 1 = Standard Option 4
// 2 = Livermore
// 3 = Penelope
#define PHYSICS_MODEL 3

namespace B4
{

class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList() : G4VModularPhysicsList()
  {
    G4LossTableManager::Instance();
    defaultCutValue = 0.1 * mm;

#if PHYSICS_MODEL == 0
    emPhysicsList = new G4EmStandardPhysics_option3();
#elif PHYSICS_MODEL == 1
    emPhysicsList = new G4EmStandardPhysics_option4();
#elif PHYSICS_MODEL == 2
    emPhysicsList = new G4EmLivermorePhysics();
#elif PHYSICS_MODEL == 3
    emPhysicsList = new G4EmPenelopePhysics();
#else
#error "Invalid PHYSICS_MODEL selected!"
#endif

    decPhysicsList = new G4DecayPhysics();

    SetVerboseLevel(1);
  }

  ~PhysicsList() override
  {
    delete emPhysicsList;
    delete decPhysicsList;
  }

  void ConstructParticle() override
  {
    emPhysicsList->ConstructParticle();
    decPhysicsList->ConstructParticle();
  }

  void ConstructProcess() override
  {
    AddTransportation();
    emPhysicsList->ConstructProcess();
    decPhysicsList->ConstructProcess();
  }

  void SetCuts() override
  {
    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");

    if (verboseLevel > 0)
      DumpCutValuesTable();
  }

private:
  G4VPhysicsConstructor* emPhysicsList;
  G4VPhysicsConstructor* decPhysicsList;
};
}
