#ifndef NuBeamTrackInformation_h
#define NuBeamTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class NuBeamTrackInformation : public G4VUserTrackInformation 
{
public:
  NuBeamTrackInformation();
  NuBeamTrackInformation(const G4Track* aTrack);
  NuBeamTrackInformation(const NuBeamTrackInformation* aTrackInfo);
  virtual ~NuBeamTrackInformation();
   
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator ==(const NuBeamTrackInformation& right) const
  {return (this==&right);}

  void Print() const;

  G4bool IsBiassed() const { return fIsBiassed; }

  G4bool SetBiassedFlag(G4bool biassed = true) {
    G4bool tmp = fIsBiassed;
    fIsBiassed = biassed;
    return tmp;
  }

  void SetCreatorModelName(G4String val) {fCreatorModelName=val;};
  G4String GetCreatorModelName() {return fCreatorModelName;};

  void SetAuxStoppingMomentum(G4ThreeVector val) {fAuxStoppingMomentum=val;};
  G4ThreeVector GetAuxStoppingMomentum() {return fAuxStoppingMomentum;};
  
private:
  G4String fCreatorModelName;
  G4int fDecayCodeDk2nu; // filled in BooNEStepping, at decay detection. 
  G4int fIsBiassed;
  G4ThreeVector fAuxStoppingMomentum;
  
public:
  inline G4int GetDecayCodeForDk2nu() const {return fDecayCodeDk2nu;}
  inline void  SetDecayCodeForDk2nu(G4int d)  {fDecayCodeDk2nu = d;}
};

extern G4Allocator<NuBeamTrackInformation> aTrackInformationAllocator;

inline void* NuBeamTrackInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void NuBeamTrackInformation::operator delete(void *aTrackInfo)
{ 
  aTrackInformationAllocator.FreeSingle((NuBeamTrackInformation*)aTrackInfo);
}

#endif
