#pragma once

#include "G4VUserPrimaryParticleInformation.hh"

class MyPrimaryParticleInfo: public G4VUserPrimaryParticleInformation
{
public:
	MyPrimaryParticleInfo(G4int scatterOrder=0): fScatterOrder(scatterOrder) {};
	virtual ~MyPrimaryParticleInfo() {};

	G4int GetScatterOrder() const { return fScatterOrder; }
	virtual void Print() const {}

private:
	G4int fScatterOrder;
};
