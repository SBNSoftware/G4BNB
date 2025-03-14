/*
 * \class BooNEpXInteraction
 *
 * \brief Container for any proton interaction on anything. Implements G4HadronicInteraction and 
 *        delegates to either BooNEpBeInteraction, BertiniCascade (LE model), or FTFP (HE model).
 *        
 *        Also calculates the weight for the secondary tracks originating from the proton.
 *        It takes the cross section arrays from the XSecArrayHolder singleton
 *        (which are themselves calculated from the user-configured physics model
 *        through logic of BooNEpBeInteraction)
 *
 * \date  Nov 26th, 2024
 */

#ifndef BooNEpXInteraction_h
#define BooNEpXInteraction_h 1

#include "G4HadronicInteraction.hh"
#include "G4Material.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4Geantino.hh"
#include "G4ReactionProduct.hh"
#include "G4TrackStatus.hh"
//
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ThreeVector.hh"

#include "NuBeamOutput.hh"
#include "NuBeamPrimaryGeneratorAction.hh"
#include "BooNEHadronInelasticDataSet.hh"
#include "BooNEHadronQuasiElasticModel.hh"
#include "BooNEpBeInteraction.hh"
#include "XSecArrayHolder.hh"

#include "G4HadronicProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4HadFinalState.hh"

#include "G4RunManager.hh"

class BooNEHadronPhysics;
class BooNEpBeInteraction;
class XSecArrayHolder;

class BooNEpXInteraction : public G4HadronicInteraction
{
public:
  BooNEpXInteraction();
  BooNEpXInteraction(G4HadronicInteraction * hadInt, 
		     G4bool proton=true); // to gracefully call ApplyYourself()
  ~BooNEpXInteraction();

private:
  static const G4int kNPzBins=50;
  static const G4int kNPtBins=100;
  static const G4int kNProtonMomentumBins=100;

public:

  // implement the G4HadronicInteraction interface
  
  G4HadFinalState * ApplyYourself( const G4HadProjectile & aTrack,
				   G4Nucleus & targetNucleus );

  G4bool IsApplicable( const G4HadProjectile & /*aTrack*/, 
		       G4Nucleus & /*targetNucleus*/ ) { return true; }
  const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const
  { return std::pair<G4double, G4double>( 100.0 * CLHEP::perCent, 1000.0 * CLHEP::GeV ); }

  ///! Workhorses for calculating the weight, copied from BooNEpBeInteraction
  G4double GetInverseRwgtFactor( G4double protonMomentum, G4double daughterPz, 
				 G4double daughterPt, G4int G3PartID );

  G4double GetInterpolatedXSec( G4double XSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins],
				G4double protonMomentum, G4double daughterPz, G4double daughterPt,
				G4bool debug = false);

  ///! Diagnostics for interaction
  G4bool   HasInteraction()  { return fInteraction != 0; }
  G4String GetIntModelName() { return ( this->HasInteraction() ) ? fInteraction->GetModelName() : this->GetModelName(); }
  // overload GetModelName()
  G4String GetModelName() { return this->GetIntModelName(); }

private:
  G4HadronicInteraction * fInteraction;
  G4bool fIsProton = true;

  ///! Bins
  // total momentum (GeV) for incident protons.
  const G4double fProtonMomentumBins[kNProtonMomentumBins] = {
    0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9,
    3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.1, 5.3, 5.5, 5.7, 5.9,
    6.1, 6.3, 6.5, 6.7, 6.9, 7.1, 7.3, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 8.7, 8.9,
    9.1, 9.3, 9.5, 9.7, 9.9,
    10.1, 10.3, 10.5, 10.7, 10.9, 11.1, 11.3, 11.5, 11.7, 11.9,
    12.1, 12.3, 12.5, 12.7, 12.9, 13.1, 13.3, 13.5, 13.7, 13.9,
    14.1, 14.3, 14.5, 14.7, 14.9, 15.1, 15.3, 15.5, 15.7, 15.9,
    16.1, 16.3, 16.5, 16.7, 16.9, 17.1, 17.3, 17.5, 17.7, 17.9,
    18.1, 18.3, 18.5, 18.7, 18.9, 19.1, 19.3, 19.5, 19.7, 19.9
  };
  
  // longitudinal momentum (GeV/c) array for secondary tracks
  const G4double fPzVec[kNPzBins] = {
    0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9,
    3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.1, 5.3, 5.5, 5.7, 5.9,
    6.1, 6.3, 6.5, 6.7, 6.9, 7.1, 7.3, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 8.7, 8.9,
    9.1, 9.3, 9.5, 9.7, 9.9
  };
  
  // transverse momentum (GeV/c) array for secondary tracks
  const G4double fPtVec[kNPtBins] = {
    .005, .015, .025, .035, .045, .055, .065, .075, .085, .095,
    .105, .115, .125, .135, .145, .155, .165, .175, .185, .195,
    .205, .215, .225, .235, .245, .255, .265, .275, .285, .295,
    .305, .315, .325, .335, .345, .355, .365, .375, .385, .395,
    .405, .415, .425, .435, .445, .455, .465, .475, .485, .495,
    .505, .515, .525, .535, .545, .555, .565, .575, .585, .595,
    .605, .615, .625, .635, .645, .655, .665, .675, .685, .695,
    .705, .715, .725, .735, .745, .755, .765, .775, .785, .795,
    .805, .815, .825, .835, .845, .855, .865, .875, .885, .895,
    .905, .915, .925, .935, .945, .955, .965, .975, .985, .995
  };
};

#endif // #ifndef BooNEpXInteraction_h
