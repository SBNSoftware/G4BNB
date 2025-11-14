#include "BooNEHadronElasticProcess.hh"
#include "NuBeamTrackInformation.hh"
#include "NuBeamTrajectory.hh"
#include "NuBeamTrajectoryContainer.hh"

BooNEHadronElasticProcess::
BooNEHadronElasticProcess(const G4String& processName)
  :  G4HadronElasticProcess(processName)
{
}

BooNEHadronElasticProcess::~BooNEHadronElasticProcess()
{
}

G4VParticleChange* BooNEHadronElasticProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  G4ParticleChange* theParticleChange=dynamic_cast<G4ParticleChange*> (G4HadronicProcess::PostStepDoIt(aTrack, aStep));
 
  //Break tracks whenever there is Hadron elastic process to get elastic scatterings into ancestry tree
  G4Track* newTrack=new G4Track(aTrack); 
  newTrack->SetPosition(*theParticleChange->GetPosition());
  newTrack->SetMomentumDirection(*theParticleChange->GetMomentumDirection());
  newTrack->SetKineticEnergy(theParticleChange->GetEnergy());
  newTrack->SetPolarization(*theParticleChange->GetPolarization());
  newTrack->SetVelocity(theParticleChange->GetVelocity());
  
  std::vector<G4Track*> secVec;
  for (int i=0;i<theParticleChange->GetNumberOfSecondaries();i++) {
    G4Track* secTrack=theParticleChange->GetSecondary(i);
    secVec.push_back(secTrack);
  }
  G4int nsec=theParticleChange->GetNumberOfSecondaries();
  theParticleChange->Clear();
  secVec.push_back(newTrack);

  NuBeamTrajectoryContainer & tinst = NuBeamTrajectoryContainer::Instance();
  int nprev = tinst.GetNTrajectories();
  NuBeamTrajectory traj = NuBeamTrajectory(&aTrack);
  tinst.AddTrajectory( traj );
  int naft = tinst.GetNTrajectories();
  if(nprev != naft)
    G4cout << "Adding trajectory with ID " << aTrack.GetTrackID() << "..." << G4endl;

  /*
  //now revert the old track info to get it stored properly in the last step
  //of ancestry tree
  theParticleChange->ProposeMomentumDirection(aTrack.GetMomentumDirection());
  theParticleChange->ProposeEnergy(aTrack.GetKineticEnergy());
  theParticleChange->ProposePolarization(aTrack.GetPolarization());
  theParticleChange->ProposeVelocity(aTrack.GetVelocity());
  theParticleChange->ProposeTrackStatus(fStopAndKill);
  theParticleChange->SetNumberOfSecondaries(nsec+1); //what we had plus the original track
  for (size_t i=0;i<secVec.size();i++) {
    NuBeamTrackInformation* tInfo=dynamic_cast<NuBeamTrackInformation*>
      (secVec[i]->GetUserInformation());
      if (tInfo==NULL) {
	tInfo=new NuBeamTrackInformation();
	secVec[i]->SetUserInformation(tInfo);
      }
      tInfo->SetCreatorModelName(GetHadronicInteraction()->GetModelName());   
      theParticleChange->AddSecondary(secVec[i]);
  }
  */
  
  return theParticleChange;
}
