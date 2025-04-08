#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "MyStackingAction.hh"

using namespace B4;

namespace B4c
{


ActionInitialization::ActionInitialization(G4String energy): G4VUserActionInitialization(), fEnergyStr(energy)
{}


ActionInitialization::~ActionInitialization()
{}


void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction(fEnergyStr));
}


void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction(fEnergyStr));
  SetUserAction(new RunAction(fEnergyStr));
  SetUserAction(new EventAction);
  SetUserAction(new SteppingAction);
  SetUserAction(new MyStackingAction);
}

}
