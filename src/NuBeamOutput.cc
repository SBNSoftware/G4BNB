#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "NuBeamOutput.hh"
#include "NuBeamOutputMessenger.hh"
#include "NuBeamTrackInformation.hh"
#include "NuBeamSkinDepthField.hh"
#include "NuBeamLocalField.hh"
#include "NuBeamRunManager.hh"
#include "NuBeamSteppingAction.hh"
#include "NuBeamGeometryConstruction.hh"
#include "NuBeamPrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4TrajectoryContainer.hh"
#include "G4TrajectoryPoint.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4UImessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "G4Version.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4FieldManager.hh"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/readWeightLocations.h"
#include "dk2nu/tree/calcLocationWeights.h"

//Root includes
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TH1D.h"

NuBeamOutput::NuBeamOutput() :
  fMessenger(0)
{
  fMessenger = new NuBeamOutputMessenger(this);
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  UI->ApplyCommand("/boone/output/filename");
  UI->ApplyCommand("/boone/output/saveProductionNtuple");
  UI->ApplyCommand("/boone/output/saveAfterHornNtuple");
  UI->ApplyCommand("/boone/output/nuEnergyThr");
  UI->ApplyCommand("/boone/output/pionMomentumThr");
  UI->ApplyCommand("/boone/output/muonMomentumThr");
  UI->ApplyCommand("/boone/output/kaonMomentumThr");
}

NuBeamOutput::~NuBeamOutput()
{  
  delete fMessenger;
  fMessenger = 0;
}


void NuBeamOutput::RecordBeginOfRun(const G4Run*)
{
  fDk2Nu  = new bsim::Dk2Nu;
  fDkMeta = new bsim::DkMeta;
  if (fOutputFileName.find(".root") == std::string::npos) fOutputFileName += std::string(".dk2nu.root");
  fOutFileDk2Nu = new TFile(fOutputFileName.c_str(),"RECREATE","root ntuple");  
  fOutTreeDk2Nu = new TTree("dk2nuTree", "Neutrino ntuple, dk2Nu format");
  fOutTreeDk2Nu->Branch("dk2nu", "bsim::Dk2Nu", &fDk2Nu, 32000, 99);
  fOutTreeDk2NuMeta  = new TTree("dkmetaTree", "Neutrino ntuple, dk2Nu format, metadata");
  fOutTreeDk2NuMeta->Branch("dkmeta", "bsim::DkMeta", &fDkMeta, 32000,99);
  this->fillDkMeta();   

  if (fSaveProductionNtuple) {
    G4cout<<"Creating aux tree: production"<<G4endl;
    G4cout << "\tRecord pBe inelastic interactions (excluding quasi-elastic)"<<G4endl;
    fProductionTree = new TTree("production"        , "Production Ntuple");
    fProductionTree->Branch("nmult", &fNtuple102.nmult  ,"nmult/I");
    fProductionTree->Branch("id"   , fNtuple102.id      ,"id[nmult]/I");
    fProductionTree->Branch("mom"  , fNtuple102.mom     ,"mom[nmult][3]/F");
    fProductionTree->Branch("wgt"  , fNtuple102.wgt     ,"wgt[nmult]/F");
  }

  for (unsigned int i=0; i<fBoundaryNtp.size();i++) {
    G4cout << "Creating aux tree: "<<fBoundaryNtp[i].fName<<G4endl;
    G4cout << "\tRecord particles on boundary of "<<fBoundaryNtp[i].fExitVol<<" and "<<fBoundaryNtp[i].fEnterVol<<G4endl;
    G4cout << "\t(note that no checks were done to ensure these volumes really exist)"<<G4endl;

    fBoundaryNtp[i].fTree= new TTree(fBoundaryNtp[i].fName, Form("Particles crossing %s and %s boundary",fBoundaryNtp[i].fExitVol.c_str(),fBoundaryNtp[i].fEnterVol.c_str()));
    fBoundaryNtp[i].fTree->Branch("pdg"    , &fBoundaryNtp[i].fNtp.pdg    ,"pdg/I");
    fBoundaryNtp[i].fTree->Branch("beamwgt", &fBoundaryNtp[i].fNtp.beamwgt,"beamwgt/F");
    fBoundaryNtp[i].fTree->Branch("pos"    ,  fBoundaryNtp[i].fNtp.pos    ,"pos[3]/F");
    fBoundaryNtp[i].fTree->Branch("mom"    ,  fBoundaryNtp[i].fNtp.mom    ,"mom[3]/F");
    fBoundaryNtp[i].fTree->Branch("eng"    , &fBoundaryNtp[i].fNtp.eng    ,"eng/F");
  }
}

void NuBeamOutput::fillDkMeta() 
{
   const NuBeamRunManager* theRunManager = dynamic_cast<const NuBeamRunManager*>(G4RunManager::GetRunManager());
   
   fDkMeta->job = theRunManager->GetCurrentRun()->GetRunID();
   fDkMeta->pots = theRunManager->GetCurrentRun()->GetNumberOfEventToBeProcessed();
   fDkMeta->beamsim = std::string("NuBeam/v1");
   fDkMeta->physics = G4Version;
   //   fDkMeta->physcuts = theRunManager->GetPhysicsListName();
   fDkMeta->tgtcfg = std::string("BNB");
   fDkMeta->horncfg = std::string("FHC");
   fDkMeta->dkvolcfg = std::string("Air");
   //
   // Beam parameters 
   //
   const NuBeamPrimaryGeneratorAction* NPGA=
      static_cast<const NuBeamPrimaryGeneratorAction*> (theRunManager->GetUserPrimaryGeneratorAction());
   fDkMeta->beam0x = NPGA->GetPrimaryMeanX();
   fDkMeta->beam0y = NPGA->GetPrimaryMeanY();
   fDkMeta->beam0z = NPGA->GetPrimaryZPosition();
   fDkMeta->beamhwidth = NPGA->GetPrimarySigmaX();
   fDkMeta->beamvwidth = NPGA->GetPrimarySigmaY();
   fDkMeta->beamdxdz = 0.;
   fDkMeta->beamdydz = 0.;
   // Detector positions 
   fDkMeta->location.clear();
   bsim::Location alocationM(0., 189.614, 54134.0, std::string("MiniBooNE")); // PRD says 541 meters..Refined by Zarko, e-mail April 3rd 
   fDkMeta->location.push_back(alocationM);
   bsim::Location alocationMi(0., 0., 47000., std::string("MicroBooNE")); 
   fDkMeta->location.push_back(alocationMi);
   bsim::Location alocationSc(0., 0. , 10000., std::string("SciBooNE")); // Approximate. PRD
   //  http://arxiv.org/ct?url=http%3A%2F%2Fdx.doi.org%2F10%252E1103%2FPhysRevD%252E86%252E052009&v=6c8af350
   fDkMeta->location.push_back(alocationSc);
   bsim::Location alocationT600(0., 0.  , 60000., std::string("T600"));
   fDkMeta->location.push_back(alocationT600);
   bsim::Location alocationSBND(73.78, 0.  , 11000., std::string("SBND"));
   fDkMeta->location.push_back(alocationSBND);
   fOutTreeDk2NuMeta->Fill();
}

void NuBeamOutput::RecordEndOfRun(const G4Run* )
{
  fOutFileDk2Nu->cd();
  if (fSaveProductionNtuple) fProductionTree->Write(); 
  for (unsigned int i=0; i<fBoundaryNtp.size();i++) 
      fBoundaryNtp[i].fTree->Write();

  fOutTreeDk2Nu->Write();
  fOutTreeDk2NuMeta->Write();
  fOutFileDk2Nu->Close();
  
  delete fOutFileDk2Nu;

  
}
void NuBeamOutput::RecordpBeInteraction(G4HadFinalState* aParticleChange)
{
  if (!fSaveProductionNtuple) return;

  fNtuple102.nmult = aParticleChange->GetNumberOfSecondaries();
  for (G4int jsec=0; jsec < aParticleChange->GetNumberOfSecondaries();
       jsec++){
    G4HadSecondary* SecondaryTrack = aParticleChange->GetSecondary(jsec);
    fNtuple102.id[jsec] = 
      SecondaryTrack->GetParticle()->GetDefinition()
      ->GetPDGEncoding();
    fNtuple102.mom[jsec][0] = 
      (SecondaryTrack->GetParticle()->GetMomentum().x())/CLHEP::GeV;
    fNtuple102.mom[jsec][1] = 
      (SecondaryTrack->GetParticle()->GetMomentum().y())/CLHEP::GeV;
    fNtuple102.mom[jsec][2] = 
      (SecondaryTrack->GetParticle()->GetMomentum().z())/CLHEP::GeV;
    fNtuple102.wgt[jsec]    = SecondaryTrack->GetWeight();
  }
  fProductionTree->Fill();
}
void NuBeamOutput::RecordBeginOfEvent(const G4Event*)
{
}

void NuBeamOutput::RecordEndOfEvent(const G4Event*)
{
  /*
  G4TrajectoryContainer* trajectoryContainer =
    a_event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer)
    n_trajectories = trajectoryContainer->entries();

  // loop over interesting trajectories for this event, and write
  // to file the trajectories with a neutrino as the daughter track
  // (PDGEncoding: +/-12, +/- 14)

  G4int NumberOfStoredTrj = 0;
  for (G4int i=0; i<n_trajectories; i++)
    {
      //current trajectory
      trj = (NuBeamTrajectory*)
	((*(a_event->GetTrajectoryContainer()))[i]);

      parentTrj = trj->GetParentTrajectory();

      if(parentTrj != 0) {
	bool writeTrack =
	  (trj->GetPDGEncoding() != parentTrj->GetPDGEncoding() &&
	   (abs(trj->GetPDGEncoding()) == 12 ||
	    abs(trj->GetPDGEncoding()) == 14));
	
	if  (writeTrack) {
	  NumberOfStoredTrj = NumberOfStoredTrj + 1;
	  // replace final muon polarization from final pion momentum
	  // and initial muon momentum. Pion is grandparent, in this case.
	  
	  if (trj->GetParentID()!=0) {
	    grandParentTrj = parentTrj->GetParentTrajectory();
	    if(grandParentTrj != 0) {	    
	      G4bool ParentIsMuon = abs(parentTrj->GetPDGEncoding())==13;
	      G4bool GrandParentIsPion = abs(grandParentTrj->GetPDGEncoding())==211;
	      
	      if (ParentIsMuon && GrandParentIsPion) {		  
		double chargeG4 = parentTrj->GetCharge();
		double epiG4 = sqrt(
				    pow(grandParentTrj->GetFinalMomentum().x(),2)+
				    pow(grandParentTrj->GetFinalMomentum().y(),2)+
				    pow(grandParentTrj->GetFinalMomentum().z(),2)+
				    pow(G4PionPlus::PionPlusDefinition()->GetPDGMass(),2));
		G4ThreeVector ppiG4 = grandParentTrj->GetFinalMomentum();
		double emuG4 = sqrt(
				    pow(parentTrj->GetInitialMomentum().x(),2)+
				    pow(parentTrj->GetInitialMomentum().y(),2)+
				    pow(parentTrj->GetInitialMomentum().z(),2)+
				    pow(G4MuonPlus::MuonPlusDefinition()->GetPDGMass(),2));
		G4ThreeVector  pmuG4 = parentTrj->GetInitialMomentum();
		
		G4ThreeVector muPolarization = GetMuonPolarization(
						       chargeG4,
						       epiG4,
						       ppiG4,
						       emuG4,
						       pmuG4);
	      }
	    }
	  }
	} 
      }
    }
  */
}
void NuBeamOutput::RecordBeginOfTrack(const G4Track*)
{
}

void NuBeamOutput::RecordNeutrino(const G4Track* track)
{
  ///
  // This must be a neutrino!... See clause Tracking action call to this. 
  //
  NuBeamRunManager *pRunManager=
    reinterpret_cast<NuBeamRunManager*>(G4RunManager::GetRunManager());
  
  G4ThreeVector pos = track->GetPosition()/CLHEP::mm; 
  const double x = pos.x();
  const double y = pos.y();
  const double z = pos.z();
  G4ThreeVector NuMomentum = track->GetMomentum();
  G4int parentID           = track->GetParentID();
  
  NuBeamTrajectory* NuParentTrack   = GetTrajectory(parentID);
  G4ThreeVector ParentMomentumFinal = NuParentTrack->GetFinalMomentum();
  G4ThreeVector vertex_r            = NuParentTrack->GetFinalPosition(); //Should be the same as Neutrino vertex
  G4double Parent_mass              = NuParentTrack->GetMass();
  G4double gamma                    = sqrt(ParentMomentumFinal*ParentMomentumFinal+Parent_mass*Parent_mass)/Parent_mass; 
  G4double Parent_energy            = gamma*Parent_mass;
  G4ThreeVector beta_vec            = ParentMomentumFinal/Parent_energy;
  G4double partial                  = gamma*(beta_vec*NuMomentum);
  G4double enuzr                    = gamma*(track->GetTotalEnergy())-partial; //neutrino energy in parent rest frame
  G4double enuzrInGeV               = enuzr/CLHEP::GeV; //neutrino energy in parent rest frame, in GeV 
  G4String parent_name              = NuParentTrack->GetParticleName();

  //check thresholds
  if (NuMomentum.mag()<fNuEnergyThr) return;
  if ((NuParentTrack->GetParticleDefinition()==G4PionPlus::PionPlus() ||
       NuParentTrack->GetParticleDefinition()==G4PionMinus::PionMinus()) &&
      ParentMomentumFinal.mag()<fPionMomentumThr) return;
  if ((NuParentTrack->GetParticleDefinition()==G4MuonPlus::MuonPlus() ||
       NuParentTrack->GetParticleDefinition()==G4MuonMinus::MuonMinus()) &&
      ParentMomentumFinal.mag()<fMuonMomentumThr) return;
  if ((NuParentTrack->GetParticleDefinition()==G4KaonPlus::KaonPlus() ||
       NuParentTrack->GetParticleDefinition()==G4KaonMinus::KaonMinus() ||
       NuParentTrack->GetParticleDefinition()==G4KaonZeroLong::KaonZeroLong()) &&
      ParentMomentumFinal.mag()<fKaonMomentumThr) return;


  G4int Norig = 2;
  if ((parent_name=="mu+") || (parent_name=="mu-")) Norig = 3;
  G4String firstvolname = NuParentTrack->GetInitialVolumeName();
  if (firstvolname.contains("TARG")) Norig = 1;
  //
  // Dk2Nu filling per se 
  // 
  fDk2Nu->nuray.clear();
  fDk2Nu->ancestor.clear();
  fDk2Nu->job = pRunManager->GetCurrentRun()->GetRunID();
  fDk2Nu->potnum = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  
  bsim::NuRay myNu(NuMomentum[0]/CLHEP::GeV, NuMomentum[1]/CLHEP::GeV, NuMomentum[2]/CLHEP::GeV, track->GetTotalEnergy()/CLHEP::GeV, 1); // Sept 24, we do not place the importance weight here.. 
  fDk2Nu->nuray.push_back(myNu);
  fDk2Nu->decay.norig=Norig; // Same convention as in G4BNB old Ntuple, i..e, as above
  fDk2Nu->decay.ntype=track->GetDefinition()->GetPDGEncoding();
  fDk2Nu->decay.ndecay = GetDecayCode(track);
  fDk2Nu->decay.vx = x/CLHEP::cm;
  fDk2Nu->decay.vy = y/CLHEP::cm;
  fDk2Nu->decay.vz = z/CLHEP::cm;
  fDk2Nu->decay.pdpx = ParentMomentumFinal[0]/CLHEP::GeV;
  fDk2Nu->decay.pdpy = ParentMomentumFinal[1]/CLHEP::GeV;
  fDk2Nu->decay.pdpz = ParentMomentumFinal[2]/CLHEP::GeV;
  G4ThreeVector ParentMomentumProduction = NuParentTrack->GetInitialMomentum();
  fDk2Nu->decay.ppdxdz = ParentMomentumProduction[0]/ParentMomentumProduction[2];
  fDk2Nu->decay.ppdydz = ParentMomentumProduction[1]/ParentMomentumProduction[2];
  fDk2Nu->decay.pppz   = ParentMomentumProduction[2]/CLHEP::GeV; 
  G4double parentp = sqrt(ParentMomentumProduction*ParentMomentumProduction);   
  fDk2Nu->decay.ppenergy = sqrt((parentp*parentp+Parent_mass*Parent_mass))/CLHEP::GeV;
  fDk2Nu->decay.ppmedium = NuParentTrack->GetInitialMaterialNumber();
  fDk2Nu->decay.ptype = NuParentTrack->GetPDGEncoding();
   
  // We now look for the grand-parent of the neutrino.
  // Use the trajectory info this time. 
  G4int parID = NuParentTrack->GetParentID();
  NuBeamTrajectory *trajGrandParent = GetTrajectory(parID);
  if (trajGrandParent != 0) {
    G4ThreeVector grandParentMomentum = trajGrandParent->GetFinalMomentum();
    G4ParticleDefinition* grandParDef = trajGrandParent->GetParticleDefinition();
    double grandParMass = grandParDef->GetPDGMass();
    fDk2Nu->decay.muparpx = grandParentMomentum[0]/CLHEP::GeV; 
    fDk2Nu->decay.muparpy = grandParentMomentum[1]/CLHEP::GeV; 
    fDk2Nu->decay.muparpz = grandParentMomentum[2]/CLHEP::GeV;
    fDk2Nu->decay.mupare =  
      sqrt(grandParentMomentum*grandParentMomentum + grandParMass*grandParMass)/CLHEP::GeV; 
   } else {
    fDk2Nu->decay.muparpx = -9999999.; 
    fDk2Nu->decay.muparpy = -9999999.; 
    fDk2Nu->decay.muparpz = -9999999.;
    fDk2Nu->decay.mupare =  -9999999.; 
  }  
  fDk2Nu->decay.necm = enuzrInGeV; // Now in GeV... 
  fDk2Nu->decay.nimpwt = track->GetWeight();
  
  std::vector<NuBeamTrajectory *> trajs;
  G4int trackIDTmp = track->GetParentID();
  while (trackIDTmp > 0) {
    NuBeamTrajectory *tmpTraj = GetTrajectory(trackIDTmp);    
    trajs.push_back(tmpTraj);
    trackIDTmp = tmpTraj->GetParentID();
    if(trackIDTmp > 0) tmpTraj = GetTrajectory(trackIDTmp);  
  }
  std::reverse(trajs.begin(), trajs.end());

  //add neutrino track info to ancestor since it is not in trajectory container yet
  G4String creatorProc=track->GetCreatorProcess()->GetProcessName()+":"+
    ((NuBeamTrackInformation*)track->GetUserInformation())->GetCreatorModelName();
  NuBeamTrajectory* nutraj=new NuBeamTrajectory(track);
  nutraj->AddTrajectoryPoint(track,creatorProc);
  nutraj->AddTrajectoryPoint(track,creatorProc);
  trajs.push_back(nutraj);

  fDk2Nu->ancestor.clear();
  fDk2Nu->vint.clear();

  // Now fill ancestry info. 
  for (auto t: trajs) {
    fDk2Nu->vint.push_back(t->GetTrackID());
    std::vector<NuBeamTrajectory::trajPoint_t> trajPoints=t->GetTrajectoryPoints();
    //hadron elastic scatterings are added as additional points in trajectory
    for (size_t iTP=0; iTP<trajPoints.size();iTP+=2) {
      bsim::Ancestor a;
      a.pdg     = t->GetPDGEncoding();
      a.startx  = trajPoints[iTP].fPosition[0]/CLHEP::cm;
      a.starty  = trajPoints[iTP].fPosition[1]/CLHEP::cm;
      a.startz  = trajPoints[iTP].fPosition[2]/CLHEP::cm;
      a.startt  = trajPoints[iTP].fTime;
      a.startpx = trajPoints[iTP].fMomentum[0]/CLHEP::GeV;
      a.startpy = trajPoints[iTP].fMomentum[1]/CLHEP::GeV;
      a.startpz = trajPoints[iTP].fMomentum[2]/CLHEP::GeV;
      a.stoppx  = trajPoints[iTP+1].fMomentum[0]/CLHEP::GeV;
      a.stoppy  = trajPoints[iTP+1].fMomentum[1]/CLHEP::GeV;
      a.stoppz  = trajPoints[iTP+1].fMomentum[2]/CLHEP::GeV;
      a.polx    = trajPoints[iTP].fPolarization[0];
      a.poly    = trajPoints[iTP].fPolarization[1];
      a.polz    = trajPoints[iTP].fPolarization[2];
      a.pprodpx = 0; //not filled, same as startp* for previous anc
      a.pprodpy = 0;
      a.pprodpz = 0;
      a.nucleus = 0; //need to add this
      a.proc    = trajPoints[iTP].fCreatorProcessName;
      a.ivol    = trajPoints[iTP].fVolumeName;
      a.imat    = trajPoints[iTP].fMaterialName;
      fDk2Nu->ancestor.push_back(a);
    }
  }
  
  if (fDk2Nu->ancestor.size() == 1) { 
    std::cerr << " Incorrect ancestry...  Final number of ancestor at evt " << fDk2Nu->potnum 
	      << " is " << fDk2Nu->ancestor.size() << " num Inel " << trajs.size() << std::endl;
    std::cerr << " Single PDG ancestor " << fDk2Nu->ancestor[0].pdg << "  Start position " 
	      << fDk2Nu->ancestor[0].startx << " / " << fDk2Nu->ancestor[0].starty 
	      << " / " << fDk2Nu->ancestor[0].startz  << std::endl;
    std::cerr << " .. start momentum " << 	fDk2Nu->ancestor[0].startpx << " / " 
	      << fDk2Nu->ancestor[0].startpy << " / " << fDk2Nu->ancestor[0].startpz << std::endl;   
  }       
  fDk2Nu->traj.clear();
  // 
  // Last step : compute the so-called "location weights" for 3 detectors. 
  //
  bsim::calcLocationWeights(fDkMeta, fDk2Nu); 
  fOutTreeDk2Nu->Fill(); 
}


void NuBeamOutput::RecordEndOfTrack(const G4Track*) 
{
} 


void NuBeamOutput::RecordStep(const G4Step* aStep)
{
  for (unsigned int i=0; i<fBoundaryNtp.size();i++) {
    G4String preStepVolName;
    G4String postStepVolName;
    if (aStep->GetPreStepPoint()->GetPhysicalVolume()) {
      preStepVolName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    } else {
      preStepVolName = "NOPOINTER";
    }
    if (aStep->GetPostStepPoint()->GetPhysicalVolume()) {
      postStepVolName = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
    } else {
      postStepVolName = "NOPOINTER";
    }

    if (preStepVolName == fBoundaryNtp[i].fExitVol && 
	postStepVolName == fBoundaryNtp[i].fEnterVol) {
      fBoundaryNtp[i].fNtp.pdg = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
      fBoundaryNtp[i].fNtp.beamwgt = aStep->GetTrack()->GetWeight();
      fBoundaryNtp[i].fNtp.pos[0] = (aStep->GetTrack()->GetPosition().x())/CLHEP::cm;
      fBoundaryNtp[i].fNtp.pos[1] = (aStep->GetTrack()->GetPosition().y())/CLHEP::cm;
      fBoundaryNtp[i].fNtp.pos[2] = (aStep->GetTrack()->GetPosition().z())/CLHEP::cm;
      fBoundaryNtp[i].fNtp.mom[0] = (aStep->GetTrack()->GetMomentum().x())/CLHEP::GeV;
      fBoundaryNtp[i].fNtp.mom[1] = (aStep->GetTrack()->GetMomentum().y())/CLHEP::GeV;
      fBoundaryNtp[i].fNtp.mom[2] = (aStep->GetTrack()->GetMomentum().z())/CLHEP::GeV;
      fBoundaryNtp[i].fNtp.eng = (aStep->GetTrack()->GetTotalEnergy())/CLHEP::GeV;
      
      fBoundaryNtp[i].fTree->Fill();
    }
  }

}

// -------------------------------------------------------

G4ThreeVector NuBeamOutput::GetMuonPolarization(
					       G4double chargeG4,
					       G4double epiG4,
					       G4ThreeVector &ppiG4,
					       G4double emuG4,
					       G4ThreeVector &pmuG4)
{

  //
  // Calculates the polarization vector of the muon in its rest frame
  // with the z axis along the direction of motion of the muon in the
  // lab frame. This is fixed by the V-A interaction since the helicity
  // of the muon in the pion frame is known.
  //
  //  Input:
  //    chargeG4     :    charge (+/- 1.0)
  //    epiG4        :    pion energy in lab frame
  //    ppiG4(3)     :    pion momentum in lab frame
  //    emuG4        :    muon energy in lab frame
  //    pmuG4(3)     :    muon momentum in lab frame
  //   
  //  Output:
  //    polmu(3)   :    muon polarization vector in rest frame of muon
  //                     with z axis pointing along its direction of motion
  //                     the three components are the spin direction.
  //
  //  Note: The program uses CERNLIB routines LORENF and LORENB in KERNLIB
  //
  //  G.B. Mills
  //  5/18/95
  //  Adapted to C by M.Sorel        , 07/01/03
  //  Adapted to Root by E.Gramellini, 12/11/14
  //
  //  Arguments:
  
  bool debug = false;
  
  TLorentzVector muLab(pmuG4.x(),pmuG4.y(),pmuG4.z(),emuG4);
  TLorentzVector piLab(ppiG4.x(),ppiG4.y(),ppiG4.z(),epiG4);
  //TLorentzVector piLab(100.,100.,200.,281.);
  TLorentzVector vMuInPiRestFrame=muLab;
  TVector3 beta (piLab.X()/(piLab.Gamma()*piLab.M()), piLab.Y()/(piLab.Gamma()*piLab.M()), piLab.Z()/(piLab.Gamma()*piLab.M()));
  TVector3 betaMu (muLab.X()/(muLab.Gamma()*muLab.M()), muLab.Y()/(muLab.Gamma()*muLab.M()), muLab.Z()/(muLab.Gamma()*muLab.M()));
  double charge = chargeG4; // why copy ? 
  
  if (debug) 
    {
      G4cout << "  charge :" << charge << G4endl;
      G4cout << "  piLab v,mass:" 
	<<piLab.X() << ", " <<piLab.Y() << ", " 
	<<piLab.Z() << ", " << piLab.T()<< ", " 
	<<piLab.M() << G4endl;

      G4cout << "  muLab v,mass:" 
	<<muLab.X() << ", " <<muLab.Y() << ", " 
	<<muLab.Z() << ", " << muLab.T()<< ", " 
	<<muLab.M() << G4endl;


      G4cout << " vMuInPiRestFrame v,mass:"
	<<vMuInPiRestFrame.X() << ", " <<vMuInPiRestFrame.Y() << ", " 
	<<vMuInPiRestFrame.Z() << ", " <<vMuInPiRestFrame.T()<< ", " 
	<<vMuInPiRestFrame.M() << G4endl;
    } 
  


  //  Boost muon vector to pion frame:
  vMuInPiRestFrame.Boost(beta);
  if (debug) {
    G4cout << "vmup,mass:" 
	   <<vMuInPiRestFrame.X()<<" "
	   <<vMuInPiRestFrame.Y()<<" "
	   <<vMuInPiRestFrame.Z()<<" "
	   <<vMuInPiRestFrame.M()<<"\n";
  }

  //  Now calculate the muon polarization vector in the pion frame
  double pmup = (vMuInPiRestFrame.Vect()).Mag();
  double bmu = pmup/vMuInPiRestFrame.T();
  double gmu = vMuInPiRestFrame.T()/vMuInPiRestFrame.M();

  TLorentzVector vPolMu(-charge*gmu*vMuInPiRestFrame.X()/pmup,
			-charge*gmu*vMuInPiRestFrame.Y()/pmup,
			-charge*gmu*vMuInPiRestFrame.Z()/pmup,
			-charge*bmu*gmu);

  TLorentzVector PolMuLab = vPolMu;  
  if (debug) 
    {
      double spin = vPolMu.M();
      G4cout << "  vpolmu,spin:" 	   
	     <<PolMuLab.X()<<" "
	     <<PolMuLab.Y()<<" "
	     <<PolMuLab.Z()<<" "
	     <<spin<<"\n";
    }
  

      
  //  Boost back to the lab  
  PolMuLab.Boost(-beta);

  if (debug) 
    {
      double spin = PolMuLab.M();
      G4cout << "  pollab,spin:" 
	     <<PolMuLab.X()<<" "
	     <<PolMuLab.Y()<<" "
	     <<PolMuLab.Z()<<" "
	     <<spin<<"\n";
    }


  TLorentzVector PolInMuRestFrame = PolMuLab;
  //  Then finally, boost to the muon rest frame
  PolInMuRestFrame.Boost(betaMu);
      
  if (debug) 
    {
      double spin = PolInMuRestFrame.M();
      G4cout << "  polmu,spin: " 
	     <<PolInMuRestFrame.X()<<" "
	     <<PolInMuRestFrame.Y()<<" "
	     <<PolInMuRestFrame.Z()<<" "
	     <<spin<<"\n";
    }

  
  return G4ThreeVector(PolInMuRestFrame.X(),PolInMuRestFrame.Y(),PolInMuRestFrame.Z());
}

NuBeamTrajectory* NuBeamOutput::GetTrajectory(G4int trackID)
{
   G4TrajectoryContainer* container = G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
   
   if(container==0) 
   {
      G4cout << "NuBeamOutput::GetTrajectory - PROBLEM: No Trajectory Container for track ID = " << trackID << std::endl;
      return 0;
   }   

   TrajectoryVector* vect = container->GetVector();
   G4VTrajectory* tr;
   G4int ii = 0; 
   while (ii < G4int(vect->size()))
   {  
      tr = (*vect)[ii]; 
      NuBeamTrajectory* tr1 = (NuBeamTrajectory*)(tr);  
      if(tr1->GetTrackID() == trackID) return tr1; 
      ++ii; 
   }

   G4cout << "NuBeamAnalysis::GetTrajectory - PROBLEM: Failed to find track with ID = " << trackID << std::endl;
   
   return 0;
}

G4int NuBeamOutput::GetDecayCode(const G4Track* nuTrack)
{
  if (!(nuTrack->GetParticleDefinition()==G4NeutrinoMu::NeutrinoMu() ||
	nuTrack->GetParticleDefinition()==G4AntiNeutrinoMu::AntiNeutrinoMu() ||
	nuTrack->GetParticleDefinition()==G4NeutrinoE::NeutrinoE() ||
	nuTrack->GetParticleDefinition()==G4AntiNeutrinoE::AntiNeutrinoE() ||
	nuTrack->GetParticleDefinition()==G4NeutrinoTau::NeutrinoTau() ||
	nuTrack->GetParticleDefinition()==G4AntiNeutrinoTau::AntiNeutrinoTau()))
    return bsim::dkp_unknown;
  
  NuBeamTrajectory* parentTraj=GetTrajectory(nuTrack->GetParentID());
  G4double parent_mass=parentTraj->GetMass();
  G4double mu_mass=G4MuonPlus::MuonPlus()->GetPDGMass();
  
  //  G4double parent_e=sqrt(parentTraj->GetFinalMomentum().mag2()+parent_mass*parent_mass);
  TLorentzVector parent_p4(parentTraj->GetFinalMomentum().x(),
			   parentTraj->GetFinalMomentum().y(),
			   parentTraj->GetFinalMomentum().z(),
			   parentTraj->GetFinalEnergy());
  TLorentzVector nu_p4(nuTrack->GetMomentum().x(),
		       nuTrack->GetMomentum().y(),
		       nuTrack->GetMomentum().z(),
		       nuTrack->GetTotalEnergy());

  switch (parentTraj->GetPDGEncoding()) {
  case 130:
  case 311:
    switch (nuTrack->GetParticleDefinition()->GetPDGEncoding()) {
    case 12:
      return bsim::dkp_k0l_nuepimep;
    case -12:
      return bsim::dkp_k0l_nuebpipem;
    case 14:
      return bsim::dkp_k0l_numupimmup;
    case -14:
      return bsim::dkp_k0l_numubpipmum;
    default:
      return bsim::dkp_unknown;
    }
    break;
  case 321: 
    switch (nuTrack->GetParticleDefinition()->GetPDGEncoding()) {
    case 14:
      //check if 3-body or 2 body decay  
      nu_p4.Boost(-parent_p4.BoostVector());
      if (fabs((parent_mass*parent_mass-mu_mass*mu_mass)/(2*parent_mass)-nu_p4.Energy())<0.001)
	return bsim::dkp_kp_numumup;
      else 
	return bsim::dkp_kp_numupi0mup;
      break;
    case 12:
      return bsim::dkp_kp_nuepi0ep;
      break;
    default:
      return bsim::dkp_unknown;
    }
    break;
  case -321: 
    switch (nuTrack->GetParticleDefinition()->GetPDGEncoding()) {
    case -14:
      //check if 3-body or 2 body decay  
      nu_p4.Boost(-parent_p4.BoostVector());
      if (fabs((parent_mass*parent_mass-mu_mass*mu_mass)/(2*parent_mass)-nu_p4.Energy())<0.001)
	return bsim::dkp_kp_numubmum;
      else 
	return bsim::dkp_kp_numubpi0mum;
    case -12:
      return bsim::dkp_kp_nuebpi0em;
    default:
      return bsim::dkp_unknown;
    }
    break;
  case 13: 
    return bsim::dkp_mup_nusep;
    break;
  case -13: 
    return bsim::dkp_mum_nusep;
    break;
  case 211:
    switch (nuTrack->GetParticleDefinition()->GetPDGEncoding()) {
    case 14:
      return bsim::dk_pip_numumup;
    case 12:
      return bsim::dkp_other;
    }
    break;
  case -211:
    switch (nuTrack->GetParticleDefinition()->GetPDGEncoding()) {
    case -14:
      return bsim::dk_pim_numubmum;
    case -12:
      return bsim::dkp_other;
    }
    break;
  default:
    return bsim::dkp_unknown;
  }
  return bsim::dkp_unknown;	
}
