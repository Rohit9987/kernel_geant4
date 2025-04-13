#pragma once

#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

class MyTrackInfo : public G4VUserTrackInformation {
public:
    MyTrackInfo(G4int scatterOrder = 0)
    :	fScatterOrder(scatterOrder), 
		fPrimaryInteractionSet(false), 
		fPrimaryInteractionPosition(G4ThreeVector()),
		isComptonElectron(0),
		isPhotoElectron(0),
		isGammaPhoton(0),
		isEBrem(0)
		{}

    virtual ~MyTrackInfo() {}

    void IncrementScatterOrder() { fScatterOrder++; }
    G4int GetScatterOrder() const { return fScatterOrder; }

	void SetComptonElectron(bool isCompton) { isComptonElectron=isCompton;}
	void SetPhotoElectron(bool isPhoto) { isPhotoElectron=isPhoto;}
	void SetGamma(bool isGamma) {isGammaPhoton = isGamma;}
	void SetEBrem(bool isEB) { isEBrem = isEB;}

 	bool GetComptonElectron() const { return isComptonElectron;}
	bool GetPhotoElectron() const { return isPhotoElectron;}
	bool GetGammaPhoton() const { return isGammaPhoton;}	
	bool GetEBrem() const { return isEBrem;}
	

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
	bool isComptonElectron;
	bool isPhotoElectron;
	bool isGammaPhoton;
	bool isEBrem;
};
