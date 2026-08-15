#ifndef NuBeamOutputMessenger_h
#define NuBeamOutputMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NuBeamOutput;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

class NuBeamOutputMessenger: public G4UImessenger
{
  public:
    NuBeamOutputMessenger(NuBeamOutput*);
   ~NuBeamOutputMessenger();
    
    void SetNewValue(G4UIcommand*, G4String );

  private:
  NuBeamOutput*  fOutput;

  G4UIdirectory*             fDirectory;
  G4UIcmdWithAString*        fOutputFileName;
  G4UIcmdWithABool*          fSaveProductionNtuple;
  G4UIcommand*               fBoundaryNtuple;
  G4UIcmdWithADoubleAndUnit* fNuEnergyThrCmd;
  G4UIcmdWithADoubleAndUnit* fNuEnergyUpperThrCmd;
  G4UIcmdWithADoubleAndUnit* fPionMomentumThrCmd;
  G4UIcmdWithADoubleAndUnit* fMuonMomentumThrCmd;
  G4UIcmdWithADoubleAndUnit* fKaonMomentumThrCmd;
};

#endif

