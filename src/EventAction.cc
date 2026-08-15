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

void EventAction::AddUnbinnedDeposit(G4double edep)
{
  if (edep > 0.0) fUnbinnedEdep += edep;
}

void EventAction::ResetEventBuffer()
{
  for (const auto linearIndex : fTouchedBins) {
    fEventEdep[linearIndex] = 0.0;
  }
  fTouchedBins.clear();
  fUnbinnedEdep = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
