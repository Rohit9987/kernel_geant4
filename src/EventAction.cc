#include "EventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "KernelBinning.hh"
#include "Run.hh"

namespace B4c
{

EventAction::EventAction()
: fEventEdep(KernelBinning::NumBins(), 0.0)
{
  fTouchedBins.reserve(256);
}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  ResetEventBuffer();
}


void EventAction::EndOfEventAction(const G4Event* /*event*/)
{
  auto* run = static_cast<MyRun*>(
    G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  for (const auto linearIndex : fTouchedBins) {
    run->AddEventBin(linearIndex, fEventEdep[linearIndex]);
  }
  run->AddUnbinnedEvent(fUnbinnedEdep);
  run->AddEscapedEvent(fEscapedEnergy);
  run->AddPrimaryEnergyEvent(fPrimaryEnergy);

  ResetEventBuffer();
}

void EventAction::AddKernelDeposit(std::size_t linearIndex, G4double edep)
{
  if (linearIndex >= fEventEdep.size() || edep <= 0.0) return;

  if (fEventEdep[linearIndex] == 0.0) {
    fTouchedBins.push_back(linearIndex);
  }
  fEventEdep[linearIndex] += edep;
}

void EventAction::AddUnbinnedDeposit(UnbinnedReason reason, G4double edep)
{
  if (edep <= 0.0) return;

  const auto index = UnbinnedReasonIndex(reason);
  if (index < fUnbinnedEdep.size()) {
    fUnbinnedEdep[index] += edep;
  }
}

void EventAction::AddEscapedEnergy(G4double energy)
{
  if (energy > 0.0) fEscapedEnergy += energy;
}

void EventAction::AddPrimaryEnergy(G4double energy)
{
  if (energy > 0.0) fPrimaryEnergy += energy;
}

void EventAction::ResetEventBuffer()
{
  for (const auto linearIndex : fTouchedBins) {
    fEventEdep[linearIndex] = 0.0;
  }
  fTouchedBins.clear();
  fUnbinnedEdep.fill(0.0);
  fEscapedEnergy = 0.0;
  fPrimaryEnergy = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
