#ifndef NuBeamStackingAction_H
#define NuBeamStackingAction_H 1

#include <fstream>
#include <map>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4UserStackingAction.hh"
#include "NuBeamStackingActionMessenger.hh"

class G4Track;
class G4RunManager;
class NuBeamStackingActionMessenger;

class NuBeamStackingAction : public G4UserStackingAction
{
public:
  NuBeamStackingAction();
  virtual ~NuBeamStackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

  
  void KillEMParticles(G4ClassificationOfNewTrack& classification, 
		       const G4Track* aTrack);
  void KillThresholdParticles(G4ClassificationOfNewTrack& classification, 
			       const G4Track * aTrack);
  
  G4double GetLowKEThreshold() const { return fLowKEThreshold; };
  void SetLowKEThreshold(G4double val){ fLowKEThreshold = val; };
  
private:
  G4RunManager * fRunManager; 
  
  NuBeamStackingActionMessenger* fMessenger;
  G4double fLowKEThreshold;
};

#endif

