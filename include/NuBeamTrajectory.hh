#ifndef BooNETrajectory_h
#define BooNETrajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"               // Include from 'system'
#include "globals.hh"               // Include from 'global'
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryPoint.hh"     // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

class G4Polyline;                   // Forward declaration.

typedef std::vector<G4VTrajectoryPoint*>  NuBeamTrajectoryPointContainer;

class NuBeamTrajectory : public G4VTrajectory
{
public:
  NuBeamTrajectory();
  
  NuBeamTrajectory(const G4Track* aTrack);
  NuBeamTrajectory(NuBeamTrajectory &);
  virtual ~NuBeamTrajectory();

  inline void* operator new(size_t);
  inline void  operator delete(void*);
  inline int operator == (const NuBeamTrajectory& right) const
  {return (this==&right);} 
  
  inline G4int GetTrackID() const
  { return fTrackID; }
  inline void SetTrackID(const G4int aValue)
  { fTrackID = aValue; };

  inline G4int GetParentID() const
  { return fParentID; }
  inline void SetParentID(const G4int aValue)
  { fParentID = aValue; }

  inline G4double GetWeight() const
  { return fWeight; }
  inline void SetWeight(const G4double aValue)
  { fWeight = aValue; }

  inline G4String GetParticleName() const
  { return fParticleName; }
  inline void SetParticleName(const G4String aValue)
  { fParticleName = aValue; }

  inline G4String GetCreatorProcessName() const
  { return fCreatorProcessName; }
  inline void SetCreatorProcessName(const G4String aValue)
  { fCreatorProcessName = aValue; }

  inline G4double GetCharge() const
  { return fPDGCharge; }
  inline void SetCharge(const G4double aValue)
  { fPDGCharge = aValue; }

  inline G4int GetPDGEncoding() const
  { return fPDGEncoding; }
  inline void SetPDGEncoding(const G4int aValue)
  { fPDGEncoding = aValue; }

  inline G4double GetMass() const
  {return fMass;}
  inline void SetMass(const G4double aValue)
  { fMass = aValue; }

  // dynamic initial track information
  inline G4double GetInitialEnergy() const
  {return fInitialEnergy;}
  inline void SetInitialEnergy(const G4double aValue)
  { fInitialEnergy = aValue; }

  inline G4ThreeVector GetInitialMomentum() const
  {return fInitialMomentum;}
  inline void SetInitialMomentum(const G4ThreeVector aValue)
  { fInitialMomentum = aValue; }

  inline G4ThreeVector GetInitialPosition() const
  {return fInitialPosition;}
  inline void SetInitialPosition(const G4ThreeVector aValue)
  { fInitialPosition = aValue; }

  inline const G4ThreeVector GetInitialPolarization() const
  {return fInitialPolarization;}
  inline void SetInitialPolarization(const G4ThreeVector aValue)
  { fInitialPolarization = aValue; }

  inline G4double GetInitialTime() const
  {return fInitialTime;}
  inline void SetInitialTime(const G4double aValue)
  { fInitialTime = aValue; }

  inline G4int GetInitialStepNumber() const
  {return fInitialStepNumber;}
  inline void SetInitialStepNumber(const G4int aValue)
  {fInitialStepNumber = aValue;}

  inline G4String GetInitialVolumeName() const {return fInitialVolumeName; }
  inline void SetInitialVolumeName(const G4String aValue)
  {fInitialVolumeName = aValue;}

  inline G4String GetInitialMaterialName() const {return fInitialMaterialName; }
  inline void SetInitialMaterialName(const G4String aValue)
  {fInitialMaterialName = aValue;}

  inline G4int GetInitialMaterialNumber() const {return fInitialMaterialNumber; }
  inline void SetInitialMaterialNumber(const G4int aValue)
  {fInitialMaterialNumber = aValue;}

  // dynamic final track information
  inline G4double GetFinalEnergy() const
  {return fFinalEnergy;}
  inline void SetFinalEnergy(const G4double aValue)
  { fFinalEnergy = aValue; }

  inline G4ThreeVector GetFinalMomentum() const
  {return fFinalMomentum;}
  inline void SetFinalMomentum(const G4ThreeVector aValue)
  { fFinalMomentum = aValue; }

  inline G4ThreeVector GetFinalPosition() const
  {return fFinalPosition;}
  inline void SetFinalPosition(const G4ThreeVector aValue)
  { fFinalPosition = aValue; }

  inline const G4ThreeVector GetFinalPolarization() const
  {return fFinalPolarization;}
  inline void SetFinalPolarization(const G4ThreeVector aValue)
  { fFinalPolarization = aValue; }

  inline G4double GetFinalTime() const
  {return fFinalTime;}
  inline void SetFinalTime(const G4double aValue)
  { fFinalTime = aValue; }

  inline G4int GetFinalStepNumber() const
  {return fFinalStepNumber;}
  inline void SetFinalStepNumber(const G4int aValue)
  {fFinalStepNumber = aValue;}


  struct trajPoint_t {
    G4String      fCreatorProcessName;
    G4double      fEnergy;
    G4ThreeVector fMomentum;
    G4ThreeVector fPosition;
    G4ThreeVector fPolarization;
    G4double      fTime;
    G4int         fStepNumber;
    G4int         fMaterialNumber;
    G4String      fMaterialName;
    G4String      fVolumeName;
  };


  void AddTrajectoryPoint(const G4Track* aTrack, G4String creatorProc);
  std::vector<trajPoint_t> GetTrajectoryPoints(){return fTrajectoryPoints;};

  // method to retrieve the parent trajectory for any trajectory
  NuBeamTrajectory* GetParentTrajectory();

  // method to retrieve daughter trajectories
  std::vector<NuBeamTrajectory*> GetDaughterTrajectories();

  virtual void ShowTrajectory(std::ostream& os=G4cout) const;
  //  virtual void DrawTrajectory() const;
  virtual void AppendStep(const G4Step* aStep);
  virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);
  
  G4ParticleDefinition* GetParticleDefinition();


  
private:
  
  NuBeamTrajectoryPointContainer* fPositionRecord;

  std::vector<trajPoint_t> fTrajectoryPoints;

  // static track information
  G4ParticleDefinition* fpParticleDefinition;
  G4int    fTrackID;
  G4int    fParentID; 
  G4double fWeight;
  G4String fParticleName;
  G4String fCreatorProcessName;
  G4double fPDGCharge;
  G4int    fPDGEncoding; 
  G4double fMass;

  // dynamic initial track information
  G4double      fInitialEnergy;
  G4ThreeVector fInitialMomentum;
  G4ThreeVector fInitialPosition;
  G4ThreeVector fInitialPolarization;
  G4double      fInitialTime;
  G4int         fInitialStepNumber;
  G4int         fInitialMaterialNumber;
  G4String      fInitialMaterialName;
  G4String      fInitialVolumeName;

  // dynamic final track information
  G4double      fFinalEnergy;
  G4ThreeVector fFinalMomentum;
  G4ThreeVector fFinalPosition;
  G4ThreeVector fFinalPolarization;
  G4double      fFinalTime;
  G4int         fFinalStepNumber;

public:

  virtual int GetPointEntries() const { return fPositionRecord->size(); }
  virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
  { return (*fPositionRecord)[i]; }
};


extern G4Allocator<NuBeamTrajectory> aTrajectoryAlloc;

inline void* NuBeamTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)aTrajectoryAlloc.MallocSingle();
  return aTrajectory;
}

inline void NuBeamTrajectory::operator delete(void* aTrajectory)
{
  aTrajectoryAlloc.FreeSingle((NuBeamTrajectory*)aTrajectory);
}

#endif










