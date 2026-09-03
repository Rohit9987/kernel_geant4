#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "KernelScoring.hh"

#include <array>
#include <cstddef>
#include <vector>

namespace B4c
{

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;

  void  BeginOfEventAction(const G4Event* event) override;
  void    EndOfEventAction(const G4Event* event) override;

  void AddKernelDeposit(std::size_t linearIndex, G4double edep);
  void AddUnbinnedDeposit(UnbinnedReason reason, G4double edep);
  void AddEscapedEnergy(G4double energy);
  void AddPrimaryEnergy(G4double energy);

private:
  void ResetEventBuffer();

  std::vector<G4double> fEventEdep;
  std::vector<std::size_t> fTouchedBins;
  std::array<G4double, kNumUnbinnedReasons> fUnbinnedEdep{};
  G4double fEscapedEnergy{0.0};
  G4double fPrimaryEnergy{0.0};
};

}


#endif
