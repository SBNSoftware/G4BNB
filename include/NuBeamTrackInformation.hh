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

  // FRAN QE - keep track if the parent led to a QE interaction
  void SetCreatorWasQE(G4bool val) {fCreatorWasQE=val;};
  G4bool GetCreatorWasQE() const {return fCreatorWasQE;};

private:
  G4String fCreatorModelName;
  G4int fDecayCodeDk2nu; // filled in BooNEStepping, at decay detection. 
  G4int fIsBiassed;

  // FRAN QE - keep track if the parent led to a QE interaction
  G4bool fCreatorWasQE;
  G4int fCreatorNProtons;
  G4int fCreatorNNeutrons;
  G4int fCreatorNPiPlus;
  G4int fCreatorNPiMinus;
  G4int fCreatorNKPlus;
  G4int fCreatorNKMinus;
  G4int fCreatorNOthers;


public:
  inline G4int GetDecayCodeForDk2nu() const {return fDecayCodeDk2nu;}
  inline void  SetDecayCodeForDk2nu(G4int d)  {fDecayCodeDk2nu = d;}

  inline G4int GetCreatorNProtons() const {return fCreatorNProtons;}
  inline void  SetCreatorNProtons(G4int n)  {fCreatorNProtons = n;}
  inline G4int GetCreatorNNeutrons() const {return fCreatorNNeutrons;}
  inline void  SetCreatorNNeutrons(G4int n)  {fCreatorNNeutrons = n;}
  inline G4int GetCreatorNPiPlus() const {return fCreatorNPiPlus;}
  inline void  SetCreatorNPiPlus(G4int n)  {fCreatorNPiPlus = n;}
  inline G4int GetCreatorNPiMinus() const {return fCreatorNPiMinus;}
  inline void  SetCreatorNPiMinus(G4int n)  {fCreatorNPiMinus = n;}
  inline G4int GetCreatorNKPlus() const {return fCreatorNKPlus;}
  inline void  SetCreatorNKPlus(G4int n)  {fCreatorNKPlus = n;}
  inline G4int GetCreatorNKMinus() const {return fCreatorNKMinus;}
  inline void  SetCreatorNKMinus(G4int n)  {fCreatorNKMinus = n;}
  inline G4int GetCreatorNOthers() const {return fCreatorNOthers;}
  inline void  SetCreatorNOthers(G4int n)  {fCreatorNOthers = n;}
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
