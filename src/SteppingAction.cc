// SteppingAction.cc
#include "SteppingAction.hh"

#include "EventAction.hh"
#include "KernelBinning.hh"
#include "MyTrackInfo.hh"
#include "Run.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

#define DEBUG 0
#if DEBUG
  #define debug(x) std::cout << "SA----> " << x << std::endl
#else
  #define debug(x)
#endif

namespace {

// ----------------- Small helpers -----------------
inline bool IsPhoton(const G4ParticleDefinition* pd)   { return pd == G4Gamma::GammaDefinition(); }

// Only processes that actually CREATE delta-electrons as secondaries
inline bool IsIonizationCreator(const G4String& n) { return (n == "eIoni" || n == "ionIoni"); }

// Attach/retrieve MyTrackInfo; overload handles const secondaries
inline MyTrackInfo* EnsureInfo(G4Track* t) {
  auto* info = dynamic_cast<MyTrackInfo*>(t->GetUserInformation());
  if (!info) { info = new MyTrackInfo(); t->SetUserInformation(info); }
  return info;
}
inline MyTrackInfo* EnsureInfo(const G4Track* ct) {
  return EnsureInfo(const_cast<G4Track*>(ct));
}

// Map ONLY photon interaction processes (NOT creators like eBrem/annihil)
inline InteractionType MapPrimaryPhotonProcess(const G4String& name) {
  if (name == "compt") return InteractionType::Compton;
  if (name == "phot")  return InteractionType::Photoelectric;
  if (name == "conv")  return InteractionType::PairProduction;
  return InteractionType::Unknown;
}

} // anon

namespace B4c {

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}
SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  auto* post = step->GetPostStepPoint();
  auto* pre  = step->GetPreStepPoint();
  const G4VProcess* process = post->GetProcessDefinedStep();
  if (!process) return;

  const G4String procName = process->GetProcessName();
  const G4ParticleDefinition* pd = track->GetParticleDefinition();
  MyTrackInfo* info = EnsureInfo(track);

  // (A) Count primary photons once per history (optional)
  if (IsPhoton(pd) && track->GetCurrentStepNumber() == 1 && track->GetCreatorProcess() == nullptr) {
    auto* run = static_cast<MyRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    if (run) run->AddPhoton();
  }

	// (B) PHOTON interactions (compt/phot/conv):
	//     - set the kernel origin at the primary photon's first interaction
	//     - increment & persist scatter order ONLY for Compton
	//     - tag secondaries: preserve Brem/Annihil lineage, else map to Compton/Photo/Pair
	if (IsPhoton(pd) && (procName == "compt" || procName == "phot" || procName == "conv"))
	{
	  if (!info->IsPrimaryInteractionSet())
		info->SetPrimaryInteraction(post->GetPosition(),
		                            pre->GetMomentumDirection());

	  const bool isCompton   = (procName == "compt");
	  const int  currOrder   = info->GetScatterOrder();
	  const int  newOrder    = isCompton ? (currOrder + 1) : currOrder;

	  if (isCompton) {
		// IMPORTANT: persist on the *current photon* so future Compton events see the incremented order
		info->SetScatterOrder(newOrder);
	  }

	  const auto* secs = step->GetSecondaryInCurrentStep();
	  if (secs && !secs->empty()) {
		const InteractionType photonLineage = info->GetParentType();
		const bool preserveLineage =
		  (photonLineage == InteractionType::BremPhoton || photonLineage == InteractionType::AnnihilationPhoton);

		for (auto* s : *secs) {
		  MyTrackInfo* sinfo = EnsureInfo(s);

		  // Propagate the *photon's* current scatter order to all secondaries of this interaction
		  sinfo->SetScatterOrder(newOrder);

		  // Preserve Brem/Annihil lineage; otherwise classify by THIS photon interaction
		  sinfo->SetParentType(preserveLineage ? photonLineage : MapPrimaryPhotonProcess(procName));

		  // All secondaries inherit the fixed kernel origin and incident direction.
		  if (!sinfo->IsPrimaryInteractionSet()) {
			sinfo->SetPrimaryInteraction(
			  info->GetPrimaryInteractionPosition(),
			  info->GetPrimaryInteractionDirection());
		  }
		}
	  }
	}

  // (C) Newly created secondaries in THIS step (creator-based tagging)
  if (const auto* secs = step->GetSecondaryInCurrentStep(); secs && !secs->empty()) {
    for (auto* s : *secs) {
      MyTrackInfo* sinfo = EnsureInfo(s);

      // If scatter order not set (e.g., δ-electrons), inherit from parent
      if (sinfo->GetScatterOrder() == 0)
        sinfo->SetScatterOrder(info->GetScatterOrder());

      if (auto* cproc = s->GetCreatorProcess()) {
        const G4String cname = cproc->GetProcessName();

        // δ-electrons from ionization: inherit parent's lineage + origin
        if (IsIonizationCreator(cname)) {
          sinfo->InheritParent(info);
        }

        // Brems photons: inherit the fixed kernel origin, but retain a separate lineage tag
        if (IsPhoton(s->GetParticleDefinition()) && cname == "eBrem") {
          sinfo->InheritParent(info);
          sinfo->SetParentType(InteractionType::BremPhoton);
        }

        // Annihilation photons: inherit the same fixed origin and retain their lineage tag
        if (IsPhoton(s->GetParticleDefinition()) && cname == "annihil") {
          sinfo->InheritParent(info);
          sinfo->SetParentType(InteractionType::AnnihilationPhoton);
        }
      }
    }
  }

  // (D) Score energy deposition
  const G4double edep = step->GetTotalEnergyDeposit();	// default unit MeV	
  if (edep == 0.0) return;

  if (!info || !info->IsPrimaryInteractionSet()) {
    // Keep this energy visible in the run summary instead of silently losing it.
    fEventAction->AddUnbinnedDeposit(edep);
    return;
  }

  const G4ThreeVector displacement =
    post->GetPosition() - info->GetPrimaryInteractionPosition();
  const G4double radius = displacement.mag();

  // A direction is undefined exactly at the kernel origin. Keep such energy
  // in the unbinned counter so that we can quantify it before choosing a
  // redistribution policy.
  if (radius == 0.0) {
    fEventAction->AddUnbinnedDeposit(edep);
    return;
  }

  const int radialBin = KernelBinning::FindRadialBin(radius / mm);

  const auto& incidentDirection = info->GetPrimaryInteractionDirection();
  G4double cosTheta = displacement.dot(incidentDirection) / radius;
  cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
  const int thetaBin = KernelBinning::FindThetaBin(std::acos(cosTheta));

  if (radialBin == KernelBinning::kInvalidBin ||
      thetaBin == KernelBinning::kInvalidBin) {
    fEventAction->AddUnbinnedDeposit(edep);
    return;
  }

  const auto linearIndex = KernelBinning::LinearIndex(
    static_cast<std::size_t>(radialBin),
    static_cast<std::size_t>(thetaBin));
  fEventAction->AddKernelDeposit(linearIndex, edep);
}

} // namespace B4c
