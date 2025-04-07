#pragma once

#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

class MyTrackInfo : public G4VUserTrackInformation {
public:
    MyTrackInfo(G4int scatterOrder = 0)
    : fScatterOrder(scatterOrder), fPrimaryInteractionSet(false), fPrimaryInteractionPosition(G4ThreeVector()) {}

    virtual ~MyTrackInfo() {}

    void IncrementScatterOrder() { fScatterOrder++; }
    G4int GetScatterOrder() const { return fScatterOrder; }

    void SetPrimaryInteractionPosition(const G4ThreeVector& pos) {
        if (!fPrimaryInteractionSet) {
            fPrimaryInteractionPosition = pos;
            fPrimaryInteractionSet = true;
        }
    }

    G4ThreeVector GetPrimaryInteractionPosition() const { return fPrimaryInteractionPosition; }
    bool IsPrimaryInteractionSet() const { return fPrimaryInteractionSet; }

private:
    G4int fScatterOrder;
    bool fPrimaryInteractionSet;
    G4ThreeVector fPrimaryInteractionPosition;
};

