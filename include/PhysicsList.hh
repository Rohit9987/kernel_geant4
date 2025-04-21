#pragma once

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

namespace B4
{

class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList() : G4VModularPhysicsList()
  {
    G4LossTableManager::Instance();
    defaultCutValue = 0.1 * mm;

    // Electromagnetic Physics
    emPhysicsList = new G4EmStandardPhysics_option4();

    // Decay Physics
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
