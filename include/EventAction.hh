#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4UserEventAction.hh"


#include "globals.hh"

namespace B4c
{

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;

  void  BeginOfEventAction(const G4Event* event) override;
  void    EndOfEventAction(const G4Event* event) override;

private:
};

}


#endif


