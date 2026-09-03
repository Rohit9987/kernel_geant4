#pragma once

#include "KernelBinning.hh"
#include "KernelScoring.hh"

#include "G4Run.hh"

#include <array>
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

    void AddUnbinnedEvent(
        const std::array<G4double, B4c::kNumUnbinnedReasons>& eventEdep) {
        G4double eventTotal = 0.0;
        for (std::size_t i = 0; i < eventEdep.size(); ++i) {
            const G4double value = eventEdep[i];
            if (value <= 0.0) continue;

            fUnbinnedEdepSum[i] += value;
            fUnbinnedEdepSumSquares[i] += value * value;
            eventTotal += value;
        }

        fTotalUnbinnedEdepSum += eventTotal;
        fTotalUnbinnedEdepSumSquares += eventTotal * eventTotal;
    }

    void AddEscapedEvent(G4double eventEnergy) {
        if (eventEnergy <= 0.0) return;

        fEscapedEnergySum += eventEnergy;
        fEscapedEnergySumSquares += eventEnergy * eventEnergy;
    }

    void AddPrimaryEnergyEvent(G4double eventEnergy) {
        if (eventEnergy <= 0.0) return;

        fPrimaryEnergySum += eventEnergy;
        fPrimaryEnergySumSquares += eventEnergy * eventEnergy;
    }

    G4double GetEdepSum(std::size_t linearIndex) const {
        return fEdepSum.at(linearIndex);
    }

    G4double GetEdepSumSquares(std::size_t linearIndex) const {
        return fEdepSumSquares.at(linearIndex);
    }

    G4double GetTotalScoredEdep() const { return fTotalScoredEdep; }
    G4double GetUnbinnedEdepSum(B4c::UnbinnedReason reason) const {
        return fUnbinnedEdepSum.at(B4c::UnbinnedReasonIndex(reason));
    }

    G4double GetUnbinnedEdepSumSquares(B4c::UnbinnedReason reason) const {
        return fUnbinnedEdepSumSquares.at(
            B4c::UnbinnedReasonIndex(reason));
    }

    G4double GetTotalUnbinnedEdepSum() const {
        return fTotalUnbinnedEdepSum;
    }

    G4double GetTotalUnbinnedEdepSumSquares() const {
        return fTotalUnbinnedEdepSumSquares;
    }

    G4double GetEscapedEnergySum() const { return fEscapedEnergySum; }
    G4double GetEscapedEnergySumSquares() const {
        return fEscapedEnergySumSquares;
    }

    G4double GetPrimaryEnergySum() const { return fPrimaryEnergySum; }
    G4double GetPrimaryEnergySumSquares() const {
        return fPrimaryEnergySumSquares;
    }

    void Merge(const G4Run* run) override {
        const MyRun* localRun = static_cast<const MyRun*>(run);
        numPhotons += localRun->numPhotons;

        for (std::size_t i = 0; i < fEdepSum.size(); ++i) {
            fEdepSum[i] += localRun->fEdepSum[i];
            fEdepSumSquares[i] += localRun->fEdepSumSquares[i];
        }
        fTotalScoredEdep += localRun->fTotalScoredEdep;

        for (std::size_t i = 0; i < fUnbinnedEdepSum.size(); ++i) {
            fUnbinnedEdepSum[i] += localRun->fUnbinnedEdepSum[i];
            fUnbinnedEdepSumSquares[i] +=
                localRun->fUnbinnedEdepSumSquares[i];
        }
        fTotalUnbinnedEdepSum += localRun->fTotalUnbinnedEdepSum;
        fTotalUnbinnedEdepSumSquares +=
            localRun->fTotalUnbinnedEdepSumSquares;
        fEscapedEnergySum += localRun->fEscapedEnergySum;
        fEscapedEnergySumSquares += localRun->fEscapedEnergySumSquares;
        fPrimaryEnergySum += localRun->fPrimaryEnergySum;
        fPrimaryEnergySumSquares += localRun->fPrimaryEnergySumSquares;

        G4Run::Merge(run);
    }

private:
    G4int numPhotons;
    std::vector<G4double> fEdepSum;
    std::vector<G4double> fEdepSumSquares;
    G4double fTotalScoredEdep{0.0};
    std::array<G4double, B4c::kNumUnbinnedReasons> fUnbinnedEdepSum{};
    std::array<G4double, B4c::kNumUnbinnedReasons>
        fUnbinnedEdepSumSquares{};
    G4double fTotalUnbinnedEdepSum{0.0};
    G4double fTotalUnbinnedEdepSumSquares{0.0};
    G4double fEscapedEnergySum{0.0};
    G4double fEscapedEnergySumSquares{0.0};
    G4double fPrimaryEnergySum{0.0};
    G4double fPrimaryEnergySumSquares{0.0};
};
