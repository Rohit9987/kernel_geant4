// SteppingAction.cc
#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "MyTrackInfo.hh"
#include "Run.hh"

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
inline bool IsElectron(const G4ParticleDefinition* pd) { return pd == G4Electron::ElectronDefinition(); }
inline bool IsPositron(const G4ParticleDefinition* pd) { return pd == G4Positron::PositronDefinition(); }

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

SteppingAction::SteppingAction() {}
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
  //     - set origin if not set (this will set it for brem/annihil photons at their first inelastic)
  //     - increment scatter order
  //     - tag secondaries: preserve Brem/Annihil lineage, else map to Compton/Photo/Pair
  if (IsPhoton(pd) && (procName == "compt" || procName == "phot" || procName == "conv"))
  {
    if (!info->IsPrimaryInteractionSet())
      info->SetPrimaryInteractionPosition(post->GetPosition());

    const int newOrder = info->GetScatterOrder() + 1;

    const auto* secs = step->GetSecondaryInCurrentStep();
    if (secs && !secs->empty()) {
      const InteractionType photonLineage = info->GetParentType();
      const bool preserveLineage =
        (photonLineage == InteractionType::BremPhoton || photonLineage == InteractionType::AnnihilationPhoton);

      for (auto* s : *secs) {
        MyTrackInfo* sinfo = EnsureInfo(s);
        sinfo->SetScatterOrder(newOrder);

        // Preserve Brem/Annihil lineage; otherwise classify by THIS photon interaction
        sinfo->SetParentType(preserveLineage ? photonLineage : MapPrimaryPhotonProcess(procName));

        // Charged secondaries inherit the (now-set) origin
        if (!IsPhoton(s->GetParticleDefinition()) && !sinfo->IsPrimaryInteractionSet())
          sinfo->SetPrimaryInteractionPosition(info->GetPrimaryInteractionPosition());
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

        // Brems photons: mark lineage and CLEAR origin (defer origin to first inelastic)
        if (IsPhoton(s->GetParticleDefinition()) && cname == "eBrem") {
          sinfo->SetParentType(InteractionType::BremPhoton);
          sinfo->ResetPrimaryInteractionPosition();
        }

        // Annihilation photons: same handling
        if (IsPhoton(s->GetParticleDefinition()) && cname == "annihil") {
          sinfo->SetParentType(InteractionType::AnnihilationPhoton);
          sinfo->ResetPrimaryInteractionPosition();
        }
      }
    }
  }

  // (D) Score energy deposition
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep == 0.0) return;

  if (!info || !info->IsPrimaryInteractionSet()) {
    // No origin yet (e.g., brem photon before first inelastic) -> nothing to score
    return;
  }

  G4ThreeVector d = post->GetPosition() - info->GetPrimaryInteractionPosition();

  int typeCode = 0;
  int scatterForOut = -1; // Only meaningful for Compton lineages

  if (IsPhoton(pd)) {
    // photons don't deposit energy in standard EM; keep for completeness
    d = G4ThreeVector(0.,0.,0.);
    typeCode = ToCode(InteractionType::Gamma);
  } else {
    const InteractionType t = info->GetParentType();
    typeCode = ToCode(t);
    if (t == InteractionType::Compton) {
      scatterForOut = info->GetScatterOrder(); // scatter order of the photon lineage
    }
  }

  auto* ana = G4AnalysisManager::Instance();
  ana->FillNtupleDColumn(0, d.x());
  ana->FillNtupleDColumn(1, d.y());
  ana->FillNtupleDColumn(2, d.z());
  ana->FillNtupleDColumn(3, edep);
  ana->FillNtupleIColumn(4, scatterForOut); // -1 for non-Compton
  ana->FillNtupleIColumn(5, typeCode);      // 0,1,2,3,4,5,6,9...
  ana->AddNtupleRow();
}

} // namespace B4c

