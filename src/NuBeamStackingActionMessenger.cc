#include "NuBeamStackingActionMessenger.hh"
#include "NuBeamStackingAction.hh"

#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UIdirectory.hh"

#include <sstream>
using namespace std;

NuBeamStackingActionMessenger::NuBeamStackingActionMessenger(NuBeamStackingAction* S) :
  fStacker(S)
{
  fDirectory = new G4UIdirectory("/boone/stacking/");
  fDirectory->SetGuidance( "BooNE stacking action control commands." );

  fLowKEThresholdCmd = new G4UIcmdWithADoubleAndUnit("/boone/stacking/lowKEThreshold", this);
  fLowKEThresholdCmd->SetGuidance("Threshold of kinetic energy below which tracks are killed");
  fLowKEThresholdCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit (e.g. MeV, GeV, etc.)");
  fLowKEThresholdCmd->SetParameterName("lowKEThreshold", true, false);
  fLowKEThresholdCmd->SetDefaultUnit("MeV");
  fLowKEThresholdCmd->SetDefaultValue(50.);
}

NuBeamStackingActionMessenger::~NuBeamStackingActionMessenger()
{
  delete fLowKEThresholdCmd;
  delete fDirectory;
}

void NuBeamStackingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fLowKEThresholdCmd)
    fStacker->SetLowKEThreshold(fLowKEThresholdCmd->GetNewDoubleValue(newValue));
}
