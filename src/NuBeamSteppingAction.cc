#include "NuBeamSteppingAction.hh"
#include "NuBeamSteppingActionMessenger.hh"
#include "NuBeamGeometryConstruction.hh"
#include "globals.hh"

#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4Navigator.hh"
#include "G4UImessenger.hh"
#include "G4UImanager.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicInteraction.hh"
#include "NuBeamRunManager.hh"

#include "NuBeamOutput.hh"
#include "NuBeamTrackInformation.hh"
#include "NuBeamTrajectory.hh"

NuBeamSteppingAction::NuBeamSteppingAction()
:fMessenger(0)
{
  const NuBeamRunManager *pRunManager=
    reinterpret_cast<const NuBeamRunManager*>(G4RunManager::GetRunManager());
  fRecords = pRunManager->GetRecordPtr();
  
  fMessenger = new NuBeamSteppingActionMessenger(this);
  
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  UI->ApplyCommand("/boone/step/perfectFocusingForPositives");
  UI->ApplyCommand("/boone/step/perfectFocusingForNegatives");
  
}

NuBeamSteppingAction::~NuBeamSteppingAction()
{
  delete fMessenger;
  fMessenger=0;
}

void NuBeamSteppingAction::UserSteppingAction(const G4Step* s)
{
  G4Track* theTrack=s->GetTrack();
  G4ParticleDefinition* particleType = theTrack->GetDefinition();

  G4bool interestingTrack = 
    !(particleType==G4Electron::ElectronDefinition() ||
      particleType==G4Positron::PositronDefinition() ||
      particleType==G4Gamma::GammaDefinition() ||
      particleType==G4Geantino::GeantinoDefinition());

  if (!interestingTrack) theTrack->SetTrackStatus(fStopAndKill);
  
  if(fPerfectFocusingForPositives || fPerfectFocusingForNegatives) {
    if ( s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "HDSK"
    	 && s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "PICV")
      {
	// if perfect focusing case, make particle on-axis and with momentum parallel
	// to z-axis
	// For technical reasons, do this by killing original track and creating a new one
	// (user is not supposed to change track information in the SteppingAction, that is
	// why the pointer s is of type const G4Step* and not G4Step*)

	if(fPerfectFocusingForPositives && theTrack->GetDefinition()->GetPDGCharge() > 0.) {
	  G4Track* newTrackPos = new G4Track(*theTrack); 
	  newTrackPos->
	    SetPosition(G4ThreeVector(0.0,0.0,theTrack->GetPosition().z()));
	  newTrackPos->SetMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
	  newTrackPos->SetParentID(theTrack->GetParentID());
	  newTrackPos->SetTrackID(theTrack->GetTrackID());
	  newTrackPos->SetWeight(theTrack->GetWeight());
	  G4TrackVector* trackVec = fpSteppingManager->GetfSecondary();
	  trackVec->push_back(newTrackPos);

	  theTrack->SetTrackStatus(fStopAndKill);
	}
	
	if(fPerfectFocusingForNegatives && theTrack->GetDefinition()->GetPDGCharge() < 0.) {
	  G4Track* newTrackNeg = new G4Track(*theTrack); 
	  newTrackNeg->
	    SetPosition(G4ThreeVector(0.0,0.0,theTrack->GetPosition().z()));
	  newTrackNeg->SetMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
	  newTrackNeg->SetParentID(theTrack->GetParentID());
	  newTrackNeg->SetTrackID(theTrack->GetTrackID());
	  newTrackNeg->SetWeight(theTrack->GetWeight());
	  G4TrackVector* trackVec = fpSteppingManager->GetfSecondary();
	  trackVec->push_back(newTrackNeg);

	  theTrack->SetTrackStatus(fStopAndKill);
	}
      }
  }

  if (fRecords != NULL) fRecords->RecordStep(s);

}


