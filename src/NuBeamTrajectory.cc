#include "NuBeamTrajectory.hh"
#include "NuBeamTrackInformation.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#endif

#include "NuBeamTrackingAction.hh"

G4Allocator<NuBeamTrajectory> aTrajectoryAlloc;

NuBeamTrajectory::NuBeamTrajectory()
{
  fPositionRecord = 0;
  fpParticleDefinition = 0;
  fTrackID = 0;
  fParentID = 0;
  fWeight = 0.0;
  fParticleName = "";
  fPDGCharge = 0.0;
  fPDGEncoding = 0;
  fMass = 0.0;

  fInitialEnergy = 0.0;
  fInitialMomentum = G4ThreeVector(0.0,0.0,0.0);
  fInitialPosition = G4ThreeVector(0.0,0.0,0.0);
  fInitialPolarization = G4ThreeVector(0.0,0.0,0.0);
  fInitialTime = 0.0;
  fInitialStepNumber = 0;

  fFinalEnergy = 0.0;
  fFinalMomentum = G4ThreeVector(0.0,0.0,0.0);
  fFinalPosition = G4ThreeVector(0.0,0.0,0.0);
  fFinalPolarization = G4ThreeVector(0.0,0.0,0.0);
  fFinalTime = 0.0;
  fFinalStepNumber = 0;

  fTrajectoryPoints.clear();
}

NuBeamTrajectory::NuBeamTrajectory(const G4Track* aTrack)
{
  // static quantities
  fPositionRecord = new NuBeamTrajectoryPointContainer();
  fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
  fpParticleDefinition = aTrack->GetDefinition();
  fTrackID = aTrack->GetTrackID();
  fParentID = aTrack->GetParentID();
  fWeight = aTrack->GetWeight();
  fParticleName = fpParticleDefinition->GetParticleName();
  fPDGCharge = fpParticleDefinition->GetPDGCharge();
  fPDGEncoding = fpParticleDefinition->GetPDGEncoding();
  fMass = fpParticleDefinition->GetPDGMass();
}

NuBeamTrajectory::NuBeamTrajectory(NuBeamTrajectory & right)
  : G4VTrajectory(right)
{
  // static quantities
  fPositionRecord = new NuBeamTrajectoryPointContainer();
  for(size_t i=0;i<right.fPositionRecord->size();i++) {
    G4TrajectoryPoint* rightPoint =
      (G4TrajectoryPoint*)((*(right.fPositionRecord))[i]);
    fPositionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
  }
  fpParticleDefinition = right.fpParticleDefinition;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  fWeight = right.fWeight;
  fParticleName = right.fParticleName;
  fPDGCharge = right.fPDGCharge;
  fPDGEncoding = right.fPDGEncoding;
  fMass = right.fMass;
}

NuBeamTrajectory::~NuBeamTrajectory()
{
  if (fPositionRecord == 0) {
    return;
  }
  size_t i;
  for(i=0;i<fPositionRecord->size();i++) {
    delete  (*fPositionRecord)[i];
  }
  fPositionRecord->clear();
  
  delete fPositionRecord;
}

void NuBeamTrajectory::AddTrajectoryPoint(const G4Track* aTrack, G4String creatorProc) 
{
  trajPoint_t newPoint;
  newPoint.fCreatorProcessName=creatorProc;
  newPoint.fEnergy=aTrack->GetTotalEnergy();
  newPoint.fMomentum=aTrack->GetMomentum();
  newPoint.fPosition=aTrack->GetPosition();
  newPoint.fPolarization=aTrack->GetPolarization();
  newPoint.fTime=aTrack->GetGlobalTime();
  newPoint.fStepNumber=aTrack->GetCurrentStepNumber();
  const G4Material *aMat = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  newPoint.fMaterialNumber=aMat->GetIndex();
  newPoint.fMaterialName=aMat->GetName();
  newPoint.fVolumeName=aTrack->GetVolume()->GetName();

  fTrajectoryPoints.push_back(newPoint);
}

void NuBeamTrajectory::ShowTrajectory(std::ostream& os) const
{
  os << G4endl << "TrackID =" << fTrackID 
     << ":ParentID=" << fParentID << G4endl;
  os << "Particle name : " << fParticleName 
     << "  Charge : " << fPDGCharge << G4endl;
  os << "  Current trajectory has " << fPositionRecord->size() 
     << " points." << G4endl;
  
  for( size_t i=0 ; i < fPositionRecord->size() ; i++) {

    G4TrajectoryPoint* aTrajectoryPoint =
      (G4TrajectoryPoint*)((*fPositionRecord)[i]);
    os << "Point[" << i << "]" 
      << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
  }
}

/*
void NuBeamTrajectory::DrawTrajectory() const
{
  const int i_mode = -1;
  // Obsolete, probably.. P.L., March 2016. 
#ifdef G4VIS_USE
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  G4ThreeVector pos;

  if(i_mode>=0)
    {
      G4Polyline pPolyline;
      for (size_t i = 0; i < fPositionRecord->size() ; i++) {
	G4TrajectoryPoint* aTrajectoryPoint =
	  (G4TrajectoryPoint*)((*fPositionRecord)[i]);
	pos = aTrajectoryPoint->GetPosition();
	pPolyline.push_back( pos );
      }

      G4Colour colour;
      if(fPDGCharge<0.) 
	// negatives are blue
	//    if (PDGEncoding == 2212)
	// protons are red
        colour = G4Colour(0.,0.,1.);
      else if(fPDGCharge>0.) 
	// positives are red
	//     else if (PDGEncoding == 130)
	// K0L are blue
        colour = G4Colour(1.,0.,0.);
      //     else if (PDGEncoding == -211)
      // pi- are green 
      else 
	// neutrals are green
        colour = G4Colour(0.,1.,0.);
      //     else if (PDGEncoding == 12)
      // nue are cyan
      //       colour = G4Colour(0.,1.,1.);
      //     else
      //       colour = G4Colour(1.,1.,1.);

      G4VisAttributes attribs(colour);
      pPolyline.SetVisAttributes(attribs);
      //      if(pVVisManager && fPDGCharge>0.) pVVisManager->Draw(pPolyline);
      if(pVVisManager) pVVisManager->Draw(pPolyline);
    }

  if(i_mode!=0) {
    for(size_t j=0; j<fPositionRecord->size(); j++) {
      G4TrajectoryPoint* aTrajectoryPoint =
	(G4TrajectoryPoint*)((*fPositionRecord)[j]);
      pos = aTrajectoryPoint->GetPosition();
      G4Circle circle( pos );
      circle.SetScreenSize(0.001*i_mode);
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colSpot(0.,0.,0.);
      G4VisAttributes attSpot(colSpot);
      circle.SetVisAttributes(attSpot);
      if(pVVisManager) pVVisManager->Draw(circle);
    }
  }
#endif
}
*/
void NuBeamTrajectory::AppendStep(const G4Step* aStep)
{
  fPositionRecord->
    push_back(new
	      G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition() ));
}

G4ParticleDefinition* NuBeamTrajectory::GetParticleDefinition()
{
  return (G4ParticleTable::GetParticleTable()->FindParticle(fParticleName));
}

void NuBeamTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  NuBeamTrajectory* seco = (NuBeamTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(G4int i=1;i<ent;i++) {
    // initial point of the second trajectory should not be merged 
    fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
    //    positionRecord->push_back(seco->positionRecord->removeAt(1));
  }
  delete (*seco->fPositionRecord)[0];
  seco->fPositionRecord->clear();
}
  
NuBeamTrajectory* NuBeamTrajectory::GetParentTrajectory()
{
  G4TrajectoryContainer* container = 
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(container==0) return 0;
  TrajectoryVector* vect = container->GetVector();
  TrajectoryVector::iterator tr = vect->begin();
  for (;tr!=vect->end();++tr) { 
    if((*tr)->GetTrackID()==fParentID) return dynamic_cast<NuBeamTrajectory*>(*tr);
  }
  return 0;
}

std::vector<NuBeamTrajectory*> NuBeamTrajectory::GetDaughterTrajectories()
{
  std::vector<NuBeamTrajectory*> daughters;
  G4TrajectoryContainer* container = 
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(container==0) return daughters;
  TrajectoryVector* vect = container->GetVector();
  TrajectoryVector::iterator tr = vect->begin();

  //G4cout << "Considering " << vect->size() << " trajectories" << G4endl;
  for( ; tr != vect->end() ; ++tr ) {
    /*
    G4cout << "Considering track with ID " << (*tr)->GetTrackID()
	   << ", parent ID " << (*tr)->GetParentID()
	   << " (there are " << container->size() << " trajectories to consider)" << G4endl;
    */
    if( (*tr)->GetParentID() == fTrackID ) 
      daughters.emplace_back( dynamic_cast<NuBeamTrajectory*>( *tr ) );
  }

  return daughters;
}
