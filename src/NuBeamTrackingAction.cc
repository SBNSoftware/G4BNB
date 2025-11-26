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
      else{
	creatorProc=aTrack->GetCreatorProcess()->GetProcessName()+":"+
	  ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorModelName();

        // Add additional suffix for BooNEHadronInelastic processes
        // For QE BooNE interactions, add :QEBooNE
        // For production (either BooNEpBeInteraction or default G4 hadronic models),
        // add :<multiplicities> (e.g. 2P1N3PiM0PiP0KP0KM0Oth)
        if ( creatorProc.find("BooNEHadronInelastic")!=G4String::npos ){

          // string with particle multiplicites
          G4String multStr = "";
          NuBeamTrackInformation* tInfo = dynamic_cast<NuBeamTrackInformation*>(aTrack->GetUserInformation());
          if(tInfo){
            multStr += std::to_string( tInfo->GetCreatorNProtons() ) + "P";
            multStr += std::to_string( tInfo->GetCreatorNNeutrons() ) + "N";
            multStr += std::to_string( tInfo->GetCreatorNPiPlus() ) + "PiM";
            multStr += std::to_string( tInfo->GetCreatorNPiMinus() ) + "PiP";
            multStr += std::to_string( tInfo->GetCreatorNKPlus() ) + "KP";
            multStr += std::to_string( tInfo->GetCreatorNKMinus() ) + "KM";
            multStr += std::to_string( tInfo->GetCreatorNOthers() ) + "Oth";
          }

          // if BooNEpBeInteraction
          if( creatorProc.find("BooNEpBeInteraction")!=G4String::npos ) {
            G4bool wasQE = ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE();
            if (wasQE) {
              creatorProc += ":QEBooNE";
            }
            else{
              creatorProc += ":" + multStr;
            }
          }
          else{
            creatorProc += ":" + multStr;
          }
        }

      }

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
      trajectory->AddTrajectoryPoint(aTrack, creatorProc, G4ThreeVector(-9999,-9999,-9999));
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

    NuBeamTrackInformation * aTrackInfo = 
    dynamic_cast<NuBeamTrackInformation *>(aTrack->GetUserInformation());
    
    // dynamic final track information
    trajectory->SetFinalEnergy( aTrack->GetTotalEnergy() );
    trajectory->SetFinalMomentum( aTrack->GetMomentum() );
    trajectory->SetFinalPosition( aTrack->GetPosition() );
    trajectory->SetFinalPolarization(aTrack->GetPolarization());
    trajectory->SetFinalTime( aTrack->GetGlobalTime() );
    trajectory->SetFinalStepNumber( aTrack->GetCurrentStepNumber() );
    trajectory->AddTrajectoryPoint(aTrack, "Final",aTrackInfo->GetAuxStoppingMomentum());
  }

  // Get the secondaries produced in this track's processes
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  size_t nSeco = secondaries->size();
  
  // Check if the interaction was inelastic hadronic
  // and add multiplicities to the track information
  bool is_ine = false;   
  for(size_t i=0;i<nSeco;i++){
    G4Track* secondary=(*secondaries)[i];
    int tpdg = secondary->GetParticleDefinition()->GetPDGEncoding();
    if( abs(tpdg)<100 )continue;
    G4String tproc = secondary->GetCreatorProcess()->GetProcessName();
    if(tproc == "BooNEHadronInelastic"){is_ine = true;}
  }
  
  // If inelastic, check secondary particles in the final state
  if( is_ine ){

    // Do the counting of final state particles
    int nProton = 0;
    int nNeutron = 0;
    int nPiMinus = 0;
    int nPiPlus = 0;
    int nKPlus = 0;
    int nKMinus = 0;
    int nOthers = 0;
    for(size_t i=0; i < nSeco; i++){	
      G4Track* secondary=(*secondaries)[i];
      int tpropdg = secondary->GetParticleDefinition()->GetPDGEncoding();
      if( abs(tpropdg)<100 || abs(tpropdg)>1000000 ) continue;
      else if( tpropdg == 2212 ) nProton++;
      else if( tpropdg == 2112 ) nNeutron++;
      else if( tpropdg == 211 ) nPiPlus++;
      else if( tpropdg == -211 ) nPiMinus++;
      else if( tpropdg == 321 ) nKPlus++;
      else if( tpropdg == -321 ) nKMinus++;
      else nOthers++;
    }

    // Now set the multiplicities in the track information of all hadronic secondaries
    for(size_t i=0; i < nSeco; i++){	
      G4Track* secondary=(*secondaries)[i];

      int tpropdg = secondary->GetParticleDefinition()->GetPDGEncoding();
      if( abs(tpropdg)<100 || abs(tpropdg)>1000000 ) continue;

      NuBeamTrackInformation* mainInfo=dynamic_cast<NuBeamTrackInformation*>(secondary->GetUserInformation());
      if (mainInfo) {
        mainInfo->SetCreatorNProtons(nProton);
        mainInfo->SetCreatorNNeutrons(nNeutron);
        mainInfo->SetCreatorNPiPlus(nPiPlus);
        mainInfo->SetCreatorNPiMinus(nPiMinus);
        mainInfo->SetCreatorNKPlus(nKPlus);
        mainInfo->SetCreatorNKMinus(nKMinus);
        mainInfo->SetCreatorNOthers(nOthers);
      }
    }

  }


  
  return;

}





