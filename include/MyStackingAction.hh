#pragma once

#include "G4UserStackingAction.hh"
#include "G4Track.hh"
#include "G4ClassificationOfNewTrack.hh"

class MyStackingAction : public G4UserStackingAction
{
public:
    MyStackingAction();
    virtual ~MyStackingAction();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;
};
