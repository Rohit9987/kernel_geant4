#pragma once

#include "G4Run.hh"

class MyRun : public G4Run {
public:
    MyRun() : numPhotons(0) {}
    virtual ~MyRun() {}

    void AddPhoton() { numPhotons++; }
    G4int GetPhotonCount() const { return numPhotons; }

    void Merge(const G4Run* run) override {
        const MyRun* localRun = static_cast<const MyRun*>(run);
        numPhotons += localRun->numPhotons;
        G4Run::Merge(run);
    }

private:
    G4int numPhotons;
};

