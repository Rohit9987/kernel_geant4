#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <cstddef>
#include <vector>

namespace B4c
{

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  void AddKernelDeposit(std::size_t linearIndex, G4double edep);
  void AddUnbinnedDeposit(G4double edep);

private:
  void ResetEventBuffer();

  std::vector<G4double> fEventEdep;
  std::vector<std::size_t> fTouchedBins;
  G4double fUnbinnedEdep{0.0};
};

} // namespace B4c

#endif
