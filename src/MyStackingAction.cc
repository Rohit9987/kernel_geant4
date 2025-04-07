#include "MyStackingAction.hh"
#include "G4Track.hh"
#include "G4PrimaryParticle.hh"
#include "MyPrimaryParticleInfo.hh"
#include "MyTrackInfo.hh"

MyStackingAction::MyStackingAction()
: G4UserStackingAction()
{}

MyStackingAction::~MyStackingAction()
{}

G4ClassificationOfNewTrack MyStackingAction::ClassifyNewTrack(const G4Track* track)
{
    // Only handle primary particles here
    if (track->GetParentID() == 0) {
        G4PrimaryParticle* primaryParticle = track->GetDynamicParticle()->GetPrimaryParticle();

        if (primaryParticle) {
            // Extract scatter order from primary particle info
            MyPrimaryParticleInfo* primaryInfo = dynamic_cast<MyPrimaryParticleInfo*>(
                primaryParticle->GetUserInformation());

            G4int scatterOrder = primaryInfo ? primaryInfo->GetScatterOrder() : 0;

            // Assign scatter order to the track using MyTrackInfo
            MyTrackInfo* trackInfo = new MyTrackInfo(scatterOrder);
            const_cast<G4Track*>(track)->SetUserInformation(trackInfo);
        } else {
            // Fallback if no primary info available
            const_cast<G4Track*>(track)->SetUserInformation(new MyTrackInfo(0));
        }
    }

    // Classify all tracks as urgent (default behavior)
    return fUrgent;
}

