#ifndef NuBeamStackingActionMessenger_h
#define NuBeamStackingActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NuBeamStackingAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;

class NuBeamStackingActionMessenger: public G4UImessenger
{
public:
   NuBeamStackingActionMessenger(NuBeamStackingAction*);
  ~NuBeamStackingActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  NuBeamStackingAction* fStacker;
  
  G4UIdirectory*             fDirectory;
  G4UIcmdWithADoubleAndUnit* fLowKEThresholdCmd;
};

#endif // #ifndef NuBeamStackingActionMessenger_h
