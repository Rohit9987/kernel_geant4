#pragma once

#include "KernelBinning.hh"

#include "G4Run.hh"

#include <cstddef>
#include <vector>

class MyRun : public G4Run {
public:
    MyRun()
    : numPhotons(0)
    , fEdepSum(B4c::KernelBinning::NumBins(), 0.0)
    , fEdepSumSquares(B4c::KernelBinning::NumBins(), 0.0)
    {}

    ~MyRun() override = default;

    void AddPhoton() { numPhotons++; }
    G4int GetPhotonCount() const { return numPhotons; }

    void AddEventBin(std::size_t linearIndex, G4double eventEdep) {
        if (linearIndex >= fEdepSum.size() || eventEdep <= 0.0) return;

        fEdepSum[linearIndex] += eventEdep;
        fEdepSumSquares[linearIndex] += eventEdep * eventEdep;
        fTotalScoredEdep += eventEdep;
    }

    void AddUnbinnedEvent(G4double eventEdep) {
        if (eventEdep <= 0.0) return;

        fUnbinnedEdepSum += eventEdep;
        fUnbinnedEdepSumSquares += eventEdep * eventEdep;
    }

    G4double GetEdepSum(std::size_t linearIndex) const {
        return fEdepSum.at(linearIndex);
    }

    G4double GetEdepSumSquares(std::size_t linearIndex) const {
        return fEdepSumSquares.at(linearIndex);
    }

    G4double GetTotalScoredEdep() const { return fTotalScoredEdep; }
    G4double GetUnbinnedEdepSum() const { return fUnbinnedEdepSum; }
    G4double GetUnbinnedEdepSumSquares() const {
        return fUnbinnedEdepSumSquares;
    }

    void Merge(const G4Run* run) override {
        const MyRun* localRun = static_cast<const MyRun*>(run);
        numPhotons += localRun->numPhotons;

        for (std::size_t i = 0; i < fEdepSum.size(); ++i) {
            fEdepSum[i] += localRun->fEdepSum[i];
            fEdepSumSquares[i] += localRun->fEdepSumSquares[i];
        }
        fTotalScoredEdep += localRun->fTotalScoredEdep;
        fUnbinnedEdepSum += localRun->fUnbinnedEdepSum;
        fUnbinnedEdepSumSquares += localRun->fUnbinnedEdepSumSquares;

        G4Run::Merge(run);
    }

private:
    G4int numPhotons;
    std::vector<G4double> fEdepSum;
    std::vector<G4double> fEdepSumSquares;
    G4double fTotalScoredEdep{0.0};
    G4double fUnbinnedEdepSum{0.0};
    G4double fUnbinnedEdepSumSquares{0.0};
};
