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

  // Printout new track if <5
  if( aTrack->GetTrackID() < 5 ) {
    //int qestatus = ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE();
    bool hasInfo = (aTrack->GetUserInformation() != 0);
    G4cout << " New track " << aTrack->GetTrackID() << " pdg " << particleType->GetPDGEncoding() << " E " << aTrack->GetKineticEnergy() << " hasInfo " << hasInfo << G4endl;
  }

  if (interestingTrack) {
    if (aTrackInfo == 0) {
      NuBeamTrackInformation* info = new NuBeamTrackInformation(aTrack);
      G4Track *aNCTrack = (G4Track *) aTrack; // defeat constness..
      aNCTrack->SetUserInformation(info);
    }

    if ( aTrack->GetTrackID() < 5 ) {
      // QE Status
      int qestatus = ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE();
      G4cout << "     --> QE status: " << qestatus << G4endl;
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

        // FRAN -- if pBe interaction, check QE
        if ( creatorProc.find("BooNEHadronInelastic")!=G4String::npos ){

          // string with particle multiplicites
          G4String multStr = "";
          NuBeamTrackInformation* tInfo = (NuBeamTrackInformation*)aTrack->GetUserInformation();
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
      trajectory->AddTrajectoryPoint(aTrack, creatorProc);

    }
  }
  else {
    fpTrackingManager->SetStoreTrajectory(false);
  }
  

  // Fran - QE flag
  G4String tproc = "";
  if( aTrack->GetTrackID() == 1 ) 
    tproc = "Primary";
  else if(aTrack->GetCreatorProcess() != 0) 
    tproc = aTrack->GetCreatorProcess()->GetProcessName();
  
  //std::cout << " New tracks " << aTrack->GetTrackID() << " pdg " << particleType->GetPDGEncoding() << " E " << aTrack->GetKineticEnergy()<< " proc " << tproc << ":" << ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorModelName() << " QE: " << ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE() << std::endl;
            
  if(tproc == "Primary") fRecords->SetTrackIdToStartProcMap( aTrack->GetTrackID(), 0); // primary
  else if(tproc == "BooNEHadronElastic"){ fRecords->SetTrackIdToStartProcMap( aTrack->GetTrackID(), 1); } // elastic
  else if(tproc == "BooNEHadronInelastic"){
    if( ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE() )
      fRecords->SetTrackIdToStartProcMap( aTrack->GetTrackID(), 2); // QE
    else
      fRecords->SetTrackIdToStartProcMap( aTrack->GetTrackID(), 3); // production
  }
  else{ fRecords->SetTrackIdToStartProcMap( aTrack->GetTrackID(), -1); }


  // Perform any record-keeping in NuBeamOutput
  const int pdg = particleType->GetPDGEncoding();
  // For dk2nu output..
  if ((fRecords != NULL) &&
      ((std::abs(pdg) == 12) ||  (std::abs(pdg) == 14) || (std::abs(pdg) == 16))) fRecords->RecordNeutrino(aTrack);
}

void NuBeamTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

  

  // Fran - QE flag
  // Get secondaries
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  size_t nSeco = secondaries->size();
  // Process ID
  G4String tprocess = "";
  if(aTrack->GetCreatorProcess() != 0) 
    tprocess = aTrack->GetCreatorProcess()->GetProcessName();

  // Printout secondary if track ID is 1
  if (aTrack->GetTrackID() == 1 ) {
  
    G4cout << "NuBeamTrackingAction::PostUserTrackingAction: TrackID " << aTrack->GetTrackID() 
     << " pdg " << aTrack->GetParticleDefinition()->GetPDGEncoding()
     << " E " << aTrack->GetKineticEnergy()
     << " QE: " << ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE()
     << " produced " << nSeco << " secondaries"
     << " momentum " << std::hypot(aTrack->GetMomentum().x(), aTrack->GetMomentum().y(), aTrack->GetMomentum().z())/CLHEP::GeV
     << G4endl;

    for( size_t i=0; i < nSeco; i++ ) {
        G4Track* secondary=(*secondaries)[i];
        G4String tproc = "";
        if(secondary->GetCreatorProcess() != 0)
          tproc = secondary->GetCreatorProcess()->GetProcessName();
        int tpdg = secondary->GetParticleDefinition()->GetPDGEncoding();
        
        if( abs(tpdg)<100 )continue;
        int qestatus = -1;//((NuBeamTrackInformation*)secondary->GetUserInformation())->GetCreatorWasQE() ? 1 : 0;
        G4cout << "   => secondary pdg " << secondary->GetParticleDefinition()->GetPDGEncoding() 
        << " trackID " << secondary->GetTrackID() 
        << " parentID " << secondary->GetParentID() 
        << " Energy " << secondary->GetKineticEnergy()
        << " momentum " << std::hypot(secondary->GetMomentum().x(), secondary->GetMomentum().y(), secondary->GetMomentum().z())/CLHEP::GeV
        << " creator proc " << tproc
        << " weight " << secondary->GetWeight() << " QE: " << qestatus
        << G4endl;

    }

  }

  
  // First check if it was an inelastic interaction
  bool is_ine = false;   
  for(size_t i=0;i<nSeco;i++){
    G4Track* secondary=(*secondaries)[i];
    int tpdg = secondary->GetParticleDefinition()->GetPDGEncoding();
    if( abs(tpdg)<100 )continue;
    G4String tproc = secondary->GetCreatorProcess()->GetProcessName();
    if(tproc == "BooNEHadronInelastic"){is_ine = true;}
    //if(tproc == "BooNEHadronElastic"){is_had_ela = true;}
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

    bool isQEMultiplicities = false;
    // Nucleon case - there cannot be other than nucleons
    if( (std::abs(aTrack->GetParticleDefinition()->GetPDGEncoding()) == 2212) ||
        (std::abs(aTrack->GetParticleDefinition()->GetPDGEncoding()) == 2112) )
    {
      if( nPiPlus == 0 && nPiMinus == 0 && nKPlus == 0 && nKMinus == 0 && nOthers == 0 ) isQEMultiplicities = true;
    }
    // Pi minus case - only one pi minus in final state !!!!WARNING: confirm this
    else if( std::abs(aTrack->GetParticleDefinition()->GetPDGEncoding()) == 211 )
    {
      if( nPiMinus == 1 && nPiPlus == 0 && nKMinus == 0 && nKPlus == 0 && nOthers == 0 ) isQEMultiplicities = true;
    }
    else if( std::abs(aTrack->GetParticleDefinition()->GetPDGEncoding()) == -211 )
    {
      if( nPiPlus == 1 && nPiMinus == 0 && nKMinus == 0 && nKPlus == 0 && nOthers == 0 ) isQEMultiplicities = true;
    }

    if(isQEMultiplicities)
      fRecords->SetTrackIdToQEMultiplicitiesMap( aTrack->GetTrackID(), 1 );
    else
      fRecords->SetTrackIdToQEMultiplicitiesMap( aTrack->GetTrackID(), 0 );

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
  else{
    fRecords->SetTrackIdToQEMultiplicitiesMap( aTrack->GetTrackID(), 0 );
  }
  

  bool isQEpBEModel =  ((NuBeamTrackInformation*)aTrack->GetUserInformation())->GetCreatorWasQE();
  fRecords->SetTrackIdToQEpBEModelMap( aTrack->GetTrackID(), isQEpBEModel );

  /*// Fran - QE 2
  for(size_t i=0;i<nSeco;i++){
    G4Track* secondary=(*secondaries)[i];
    std::cout << " AAA " << secondary->GetTrackID() << "\n";
    int tpdg = secondary->GetParticleDefinition()->GetPDGEncoding();
    if( abs(tpdg)<100 )continue;
    G4String tproc = secondary->GetCreatorProcess()->GetProcessName();
    if(tproc == "Primary") fRecords->SetTrackIdToStartProcMap( secondary->GetTrackID(), 0); // primary
    else if(tproc == "BooNEHadronElastic"){ fRecords->SetTrackIdToStartProcMap( secondary->GetTrackID(), 1); } // elastic
    else if(tproc == "BooNEHadronInelastic"){
      //std::cout << " FRAN PostUserTrackingAction: TrackID " << aTrack->GetTrackID() << " QE status of creator: "  < ((NuBeamTrackInformation*)secondary->GetUserInformation())->GetCreatorWasQE() << std::endl;
      if( ((NuBeamTrackInformation*)secondary->GetUserInformation())->GetCreatorWasQE() )
        fRecords->SetTrackIdToStartProcMap( secondary->GetTrackID(), 2); // QE
      else
        fRecords->SetTrackIdToStartProcMap( secondary->GetTrackID(), 3); // production
    }
    else{ fRecords->SetTrackIdToStartProcMap( secondary->GetTrackID(), -1); }
  }*/

  

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





