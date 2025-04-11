#ifndef B4RunAction_h
#define B4RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4String.hh"
#include "Run.hh"

class G4Run;

namespace B4
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction(G4String energyStr);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void   EndOfRunAction(const G4Run*) override;
	G4Run* GenerateRun() { return new MyRun();}

private:
	G4String fEnergyStr;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

