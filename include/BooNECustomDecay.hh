#ifndef BooNECustomDecay_h
#define BooNECustomDecay_h 1

#include <memory> // for unique_ptr

#include "BooNEDecayPhysics.hh"
#include "DARTree.hh"
#include "G4Decay.hh"
#include "G4VParticleChange.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4VExtDecayer.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"

const int max_multiplicity = 10; // max N of decay products from a single channel
struct NtupleDAR_def {
  int nmult;
  int ppdg;
  int pdg[max_multiplicity];
  float p4[max_multiplicity][4];
  float wgt[max_multiplicity];
};

class BooNECustomDecay: public G4Decay
{
public:
  BooNECustomDecay(const G4String & processName = "BooNECustomDecay");

  G4VParticleChange * PostStepDoIt(const G4Track & track, const G4Step & step);
  G4VParticleChange * AtRestDoIt(const G4Track & track, const G4Step & step);

private:
  NtupleDAR_def fNtupleDAR;
};

#endif // #ifndef BooNECustomDecay
