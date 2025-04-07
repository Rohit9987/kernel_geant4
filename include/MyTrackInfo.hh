#pragma once

#include "G4VUserTrackInformation.hh"

class MyTrackInfo: public G4VUserTrackInformation
{
public:
	MyTrackInfo(G4int scatterOrder = 0): fScatterOrder(scatterOrder) {}
	virtual ~MyTrackInfo() {}

	void IncrementScatterOrder() {fScatterOrder++;}
	G4int GetScatterOrder() const { return fScatterOrder;}

private:
	G4int fScatterOrder;
};
