#pragma once

#include "G4VModularPhysicsList.hh"
#include "G4EmLivermorePhysics.hh" // precise low-energy EM
#include "G4EmStandardPhysics_option4.hh" // precise high-energy EM
#include "G4DecayPhysics.hh"
#include "G4EmParameters.hh"
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
    defaultCutValue = 0.1 * mm; // Small cut for dose accuracy

    // Electromagnetic physics models
    emPhysicsList = new G4EmLivermorePhysics(1); // Priority 1

    // Decay
    decPhysicsList = new G4DecayPhysics();

    // Set verbose
    SetVerboseLevel(1);

    // Settings for atomic de-excitation
    auto emParams = G4EmParameters::Instance();
    emParams->SetFluo(true);          // Fluorescence
    emParams->SetAuger(true);          // Auger electrons
    emParams->SetAugerCascade(true);   // Full Auger cascades
    emParams->SetPixe(true);           // PIXE (proton-induced X-ray emission, if needed)
    emParams->SetDeexcitationIgnoreCut(true); // Always generate de-excitation products
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

