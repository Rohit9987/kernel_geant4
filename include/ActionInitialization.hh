
#ifndef B4cActionInitialization_h
#define B4cActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "G4String.hh"

namespace B4c
{

/// Action initialization class.

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(G4String energy);
    ~ActionInitialization() override;

    void BuildForMaster() const override;
    void Build() const override;

private:
	G4String fEnergyStr;
};

}

#endif


