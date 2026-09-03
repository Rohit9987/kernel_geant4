#ifndef B4cSteppingAction_h
#define B4cSteppingAction_h 1

#include "G4UserSteppingAction.hh"

namespace B4c
{

class EventAction;

class SteppingAction : public G4UserSteppingAction
{
public:
  explicit SteppingAction(EventAction* eventAction);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step* step) override;

private:
  EventAction* fEventAction;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
