// MyTrackInfo.hh
#pragma once

#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include <cstdint>

// --------------------------------------------------------------------------
// Parent interaction category used to tag the *lineage* of dose deposition.
// Keep numeric values aligned with your ntuple "type" codes.
// --------------------------------------------------------------------------
enum class InteractionType : uint8_t {
  Gamma              = 0, // photon step (logged for completeness; usually no edep)
  Compton            = 1, // lineage from Compton interaction
  Photoelectric      = 2, // lineage from photoelectric interaction
  BremPhoton         = 3, // lineage from a bremsstrahlung photon (creator: eBrem)
  AnnihilationPhoton = 4, // lineage from an annihilation photon (creator: annihil)
  PairProduction     = 5, // lineage from pair production (creator: conv)
  Unknown            = 9
};

// Optional utilities (inline so you can include this header alone)
inline int ToCode(InteractionType t) {
  return static_cast<int>(t);
}

inline const char* ToString(InteractionType t) {
  switch (t) {
    case InteractionType::Gamma:              return "Gamma";
    case InteractionType::Compton:            return "Compton";
    case InteractionType::Photoelectric:      return "Photoelectric";
    case InteractionType::BremPhoton:         return "BremPhoton";
    case InteractionType::AnnihilationPhoton: return "AnnihilationPhoton";
    case InteractionType::PairProduction:     return "PairProduction";
    default:                                  return "Unknown";
  }
}

// Map a Geant4 process name to an InteractionType (use for creator/step names)
inline InteractionType FromProcessName(const G4String& name) {
  if (name == "compt")   return InteractionType::Compton;
  if (name == "phot")    return InteractionType::Photoelectric;
  if (name == "conv")    return InteractionType::PairProduction;
  if (name == "eBrem")   return InteractionType::BremPhoton;
  if (name == "annihil") return InteractionType::AnnihilationPhoton;
  return InteractionType::Unknown;
}

// --------------------------------------------------------------------------

class MyTrackInfo : public G4VUserTrackInformation {
public:
  explicit MyTrackInfo(G4int scatterOrder = 0)
  : fScatterOrder(scatterOrder)
  , fPrimaryInteractionSet(false)
  , fPrimaryInteractionPosition()
  , fParentType(InteractionType::Unknown)
  {}

  ~MyTrackInfo() override = default;

  // -------- Scatter order (defined at photon interaction events) ----------
  void   IncrementScatterOrder()        { ++fScatterOrder; }
  void   SetScatterOrder(G4int n)       { fScatterOrder = n; }
  G4int  GetScatterOrder() const        { return fScatterOrder; }

  // -------- Parent interaction type (lineage tag) ------------------------
  void            SetParentType(InteractionType t) { fParentType = t; }
  InteractionType GetParentType() const            { return fParentType; }
  int             GetParentTypeCode() const        { return ToCode(fParentType); }

  // Copy lineage/origin from a parent track (for δ-electrons, etc.)
  void InheritParent(const MyTrackInfo* parent) {
    if (!parent) return;
    fParentType = parent->fParentType;
    // Do NOT change scatter order here unless that’s your design.
    if (parent->fPrimaryInteractionSet && !fPrimaryInteractionSet) {
      fPrimaryInteractionPosition = parent->fPrimaryInteractionPosition;
      fPrimaryInteractionSet = true;
    }
  }

  // -------- Kernel origin (first relevant interaction point) -------------
  void ResetPrimaryInteractionPosition() {
    fPrimaryInteractionSet = false;
    fPrimaryInteractionPosition = G4ThreeVector();
  }

  void SetPrimaryInteractionPosition(const G4ThreeVector& pos) {
    // If you want “first-hit wins”, keep the guard; otherwise remove it.
    if (!fPrimaryInteractionSet) {
      fPrimaryInteractionPosition = pos;
      fPrimaryInteractionSet = true;
    }
  }

  const G4ThreeVector& GetPrimaryInteractionPosition() const { return fPrimaryInteractionPosition; }
  bool                 IsPrimaryInteractionSet() const       { return fPrimaryInteractionSet; }

private:
  G4int         fScatterOrder{0};
  bool          fPrimaryInteractionSet{false};
  G4ThreeVector fPrimaryInteractionPosition{};
  InteractionType fParentType{InteractionType::Unknown};
};

