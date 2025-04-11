#ifndef B4cSteppingAction_h
#define B4cSteppingAction_h 1

#include "G4UserSteppingAction.hh"

namespace B4c
{


class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step* step) override;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
