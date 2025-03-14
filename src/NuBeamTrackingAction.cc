#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4TrackingManager.hh"
#include "NuBeamTrackingAction.hh"
#include "NuBeamTrackingActionMessenger.hh"
#include "NuBeamTrackInformation.hh"

#include "globals.hh"

#include "NuBeamOutput.hh"
#include "NuBeamTrajectory.hh"
#include "NuBeamTrackInformation.hh"
#include "NuBeamRunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"

#include "G4SteppingManager.hh"
#include "G4TrackingManager.hh"

NuBeamTrackingAction::NuBeamTrackingAction() 
  :fMessenger(0)
{
  const NuBeamRunManager *pRunManager=
    reinterpret_cast<const NuBeamRunManager*>(G4RunManager::GetRunManager());
  fRecords = pRunManager->GetRecordPtr();
 
  fMessenger = new NuBeamTrackingActionMessenger(this);

  //  G4UImanager* UI = G4UImanager::GetUIpointer();  

}

NuBeamTrackingAction::~NuBeamTrackingAction() {
  delete fMessenger;
  fMessenger = 0;
}

void NuBeamTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  G4ParticleDefinition* particleType = aTrack->GetDefinition();

  NuBeamTrackInformation * aTrackInfo = 
    dynamic_cast<NuBeamTrackInformation *>(aTrack->GetUserInformation());
  
  G4bool interestingTrack = 
    !(particleType==G4Electron::ElectronDefinition() ||
      particleType==G4Positron::PositronDefinition() ||
      particleType==G4Gamma::GammaDefinition() ||
      particleType==G4Geantino::GeantinoDefinition());

  if (interestingTrack) {

    if (aTrackInfo == 0) {
      NuBeamTrackInformation* info = 
	new NuBeamTrackInformation(aTrack);
      G4Track *aNCTrack = (G4Track *) aTrack; // defeat constness..
      aNCTrack->SetUserInformation(info);
    }
    
    // trajectory object
    if ( !fpTrackingManager->GimmeTrajectory() ) {
      NuBeamTrajectory* trajectory = new NuBeamTrajectory(aTrack);
      fpTrackingManager->SetTrajectory(trajectory);
      fpTrackingManager->SetStoreTrajectory(true);

      G4String creatorProc;
      // dynamic initial trajectory information
      if (aTrack->GetTrackID()==1) 
	creatorProc="Primary";
      else
	creatorProc=aTrack->GetCreatorProcess()->GetProcessName()+":"+
	  ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorModelName();

      trajectory->SetCreatorProcessName(creatorProc);
      trajectory->SetInitialEnergy( aTrack->GetTotalEnergy() );
      trajectory->SetInitialMomentum( aTrack->GetMomentum() );
      trajectory->SetInitialPosition( aTrack->GetPosition() );
      trajectory->SetInitialPolarization(aTrack->GetPolarization());
      trajectory->SetInitialTime( aTrack->GetGlobalTime() );
      trajectory->SetInitialStepNumber( aTrack->GetCurrentStepNumber() );
      trajectory->SetInitialVolumeName(aTrack->GetVolume()->GetName());
      const G4Material *aMat = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
      trajectory->SetInitialMaterialNumber(aMat->GetIndex()); // Too early... Must be stored in 
      trajectory->SetInitialMaterialName(aMat->GetName()); // Too early... Must be stored in 
      trajectory->AddTrajectoryPoint(aTrack, creatorProc);
    }
  } else {
    fpTrackingManager->SetStoreTrajectory(false);
  }
  
  // Perform any record-keeping in NuBeamOutput
  const int pdg = particleType->GetPDGEncoding();
  // For dk2nu output..
  if ((fRecords != NULL) &&
      ((std::abs(pdg) == 12) ||  (std::abs(pdg) == 14) || (std::abs(pdg) == 16))) fRecords->RecordNeutrino(aTrack);
}

void NuBeamTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  if (fpTrackingManager->GimmeTrajectory() && 
      fpTrackingManager->GetStoreTrajectory()) {
    NuBeamTrajectory* trajectory = 
      (NuBeamTrajectory*)fpTrackingManager->GimmeTrajectory();
    
    // dynamic final track information
    trajectory->SetFinalEnergy( aTrack->GetTotalEnergy() );
    trajectory->SetFinalMomentum( aTrack->GetMomentum() );
    trajectory->SetFinalPosition( aTrack->GetPosition() );
    trajectory->SetFinalPolarization(aTrack->GetPolarization());
    trajectory->SetFinalTime( aTrack->GetGlobalTime() );
    trajectory->SetFinalStepNumber( aTrack->GetCurrentStepNumber() );
    trajectory->AddTrajectoryPoint(aTrack, "Final");
  }
  
  return;

}





