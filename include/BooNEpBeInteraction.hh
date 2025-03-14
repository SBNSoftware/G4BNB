#ifndef BooNEpBeInteraction_h
#define BooNEpBeInteraction_h 1

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

#include "G4HadronicProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4HadFinalState.hh"

#include "G4RunManager.hh"

class BooNEHadronPhysics;
class BooNEpBeInteractionMessenger;

class XSecArrayHolder;

class BooNEpBeInteraction : public G4HadronicInteraction
{
public:
  BooNEpBeInteraction();
  ~BooNEpBeInteraction();
  
private:

  static const G4int kNEBins=20;
  static const G4int kNPzBins=50;
  static const G4int kNPtBins=100;
  static const G4int kNProtonMomentumBins=100;

public:
  bool fIsQE;
  //we don't conserve momentum in this model at all, so set levels very high
  const std::pair< G4double, G4double > GetFatalEnergyCheckLevels() const {return std::pair<G4double, G4double>(100.*CLHEP::perCent, 1000. * CLHEP::GeV);};

  bool IsApplicable(const G4HadProjectile& hadProj, G4Nucleus& nucl) {
    G4cout<<"Nucleus N="<<nucl.GetN_asInt()<<" Z="<<nucl.GetZ_asInt()<<" A="<<nucl.GetA_asInt()<<G4endl;
    G4Element el("tgt","",double(nucl.GetZ_asInt()),double(nucl.GetA_asInt())*(CLHEP::g/CLHEP::mole));
    G4cout<<"min = "<<GetMinEnergy(NULL, &el)
	  <<"max = "<<GetMaxEnergy(NULL, &el)<<G4endl;
    if (hadProj.GetKineticEnergy()>GetMinEnergy(NULL, &el) && hadProj.GetKineticEnergy()<GetMaxEnergy(NULL,&el)) 
      return true;
    else 
      return false;    
  };

    // modifiers for the quasi elastic scattering parameters
  void SetProtonQuasiElasticParameters(G4double* par, G4int npar)
                                      { fQuasiElasticModel.SetProtonBetaParameters(par , npar);    }
  void SetNeutronQuasiElasticParameters(G4double* par, G4int npar)  
                                      { fQuasiElasticModel.SetNeutronBetaParameters(par , npar);   }
  void SetPionPlusQuasiElasticParameters(G4double* par, G4int npar) 
                                      { fQuasiElasticModel.SetPionPlusBetaParameters(par , npar);  }
  void SetPionMinusQuasiElasticParameters(G4double* par, G4int npar)
                                      { fQuasiElasticModel.SetPionMinusBetaParameters(par , npar); }

  void SetPiPlusPhysicsModel(G4String val);
  void SetPiMinusPhysicsModel(G4String val);
  void SetKaonPlusPhysicsModel(G4String val);
  void SetKaonMinusPhysicsModel(G4String val);
  void SetKaonZeroLongPhysicsModel(G4String val);
  void SetProtonPhysicsModel(G4String val);
  void SetNeutronPhysicsModel(G4String val);

  void SetUseBeToAScaling(G4bool valb);
  void SetBeToAScalingPar(G4double c0, G4double c1, G4double c2);
  G4double GetBeToAScalingPar(G4int iParNo);

  void SetPhysicsVerbose(G4bool valb);
  void SetNoBeamPions(G4bool valb);

  G4double GetPrimaryEnergy() { return fPrimaryEnergy; }
  G4Material* GetTargetMaterial();

  G4String GetPiPlusPhysicsModel() { return PiPlusPhysicsModel;};
  G4String GetPiMinusPhysicsModel() { return PiMinusPhysicsModel;};
  G4String GetKaonPlusPhysicsModel() { return KaonPlusPhysicsModel;};
  G4String GetKaonMinusPhysicsModel() { return KaonMinusPhysicsModel;};
  G4String GetKaonZeroLongPhysicsModel() { return KaonZeroLongPhysicsModel;};
  G4String GetProtonPhysicsModel() { return ProtonPhysicsModel;};
  G4String GetNeutronPhysicsModel() { return NeutronPhysicsModel;};

  void SetSWPiPlusPar(
		      G4double val1, G4double val2, G4double val3,
		      G4double val4,G4double val5, G4double val6,
		      G4double val7, G4double val8);
  G4double GetSWPiPlusPar(G4int iParNo);

  void SetSWPiMinusPar(
		      G4double val1, G4double val2, G4double val3,
		      G4double val4,G4double val5, G4double val6,
		      G4double val7, G4double val8);
  G4double GetSWPiMinusPar(G4int iParNo);

  void SetSWKaonPlusPar(
		      G4double val1, G4double val2, G4double val3,
		      G4double val4,G4double val5, G4double val6,
		      G4double val7, G4double val8, G4double val9);
  G4double GetSWKaonPlusPar(G4int iParNo);

  void SetSWKaonZeroLongPar(
		      G4double val1, G4double val2, G4double val3,
		      G4double val4,G4double val5, G4double val6,
		      G4double val7, G4double val8, G4double val9);
  G4double GetSWKaonZeroLongPar(G4int iParNo);

  void SetFSKaonPlusPar(
		      G4double val1, G4double val2, G4double val3,
		      G4double val4,G4double val5, G4double val6,
		      G4double val7, G4double val8);
  G4double GetFSKaonPlusPar(G4int iParNo);

  void SetFSKaonZeroLongPar(
		      G4double val1, G4double val2, G4double val3,
		      G4double val4,G4double val5, G4double val6,
		      G4double val7, G4double val8);
  G4double GetFSKaonZeroLongPar(G4int iParNo);

  G4bool GetUseBeToAScaling() { return UseBeToAScaling; };

  G4bool GetPhysicsVerbose() { return PhysicsVerbose;};

  G4double GetpBeInelasticCrossSection( G4double totalMomentum);
  void FillOtherCrossSections(G4double kineticEnergy);


  G4HadFinalState *ApplyYourself( const G4HadProjectile &aTrack,
				    G4Nucleus &targetNucleus );
  G4double GetMax(G4double Array[kNProtonMomentumBins][kNPzBins][kNPtBins]);
  G4double GetTotalProductionXSect(G4double Array[kNProtonMomentumBins][kNPzBins][kNPtBins],
				   G4int iprotonp);
  void GetpBeSecondaryModel(BooNEHadronPhysics* thePhysics);
  void GetpBeCrossSectionModel(BooNEHadronPhysics* thePhysics,
			       G4double primaryEnergy);

  G4long GetNumberOfProtons(const G4HadProjectile &aTrack);
  G4long GetNumberOfNeutrons(const G4HadProjectile &aTrack);
  G4long GetNumberOfPiPluses(const G4HadProjectile &aTrack);
  G4long GetNumberOfPiMinuses(const G4HadProjectile &aTrack);
  G4long GetNumberOfKPluses(const G4HadProjectile &aTrack);
  G4long GetNumberOfKMinuses(const G4HadProjectile &aTrack);
  G4long GetNumberOfKZeroLongs(const G4HadProjectile &aTrack);

  G4ThreeVector GetMomentumOfSecondary(G4double Array[kNProtonMomentumBins][kNPzBins][kNPtBins],
				       G4double MaxXSec,
				       const G4HadProjectile &aTrack);

  G4double GetInterpolatedXSec(G4double XSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins],
			       G4double protonMomentum,
			       G4double daughterPz,
			       G4double daughterPt);

  G4double GetRwgtXSec(G4double NoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       G4int iprotonp, G4int ipz,
		       G4int ipt, G4int G3PartID);

  G4double GetInverseRwgtFactor(G4double protonMomentum,
				G4double daughterPz,
				G4double daughterPt,
				G4int G3PartID);

  void SetProtonRwgtFunc(G4String newValue);
  void SetNeutronRwgtFunc(G4String newValue);
  void SetPionPlusRwgtFunc(G4String newValue);
  void SetPionMinusRwgtFunc(G4String newValue);
  void SetKaonPlusRwgtFunc(G4String newValue);
  void SetKaonMinusRwgtFunc(G4String newValue);
  void SetKaonZeroLongRwgtFunc(G4String newValue);

  void SetProtonRwgtParams(G4double val0, G4double val1,
			   G4double val2, G4double val3,
			   G4double val4, G4double val5,
			   G4double val6, G4double val7,
			   G4double val8, G4double val9);
  void SetNeutronRwgtParams(G4double val0, G4double val1,
			    G4double val2, G4double val3,
			    G4double val4, G4double val5,
			    G4double val6, G4double val7,
			    G4double val8, G4double val9);
  void SetPionPlusRwgtParams(G4double val0, G4double val1,
			     G4double val2, G4double val3,
			     G4double val4, G4double val5,
			     G4double val6, G4double val7,
			     G4double val8, G4double val9);
  void SetPionMinusRwgtParams(G4double val0, G4double val1,
			      G4double val2, G4double val3,
			      G4double val4, G4double val5,
			      G4double val6, G4double val7,
			      G4double val8, G4double val9);
  void SetKaonPlusRwgtParams(G4double val0, G4double val1,
			     G4double val2, G4double val3,
			     G4double val4, G4double val5,
			     G4double val6, G4double val7,
			     G4double val8, G4double val9);
  void SetKaonMinusRwgtParams(G4double val0, G4double val1,
			      G4double val2, G4double val3,
			      G4double val4, G4double val5,
			      G4double val6, G4double val7,
			      G4double val8, G4double val9);
  void SetKaonZeroLongRwgtParams(G4double val0, G4double val1,
				 G4double val2, G4double val3,
				 G4double val4, G4double val5,
				 G4double val6, G4double val7,
				 G4double val8, G4double val9);

private:

  BooNEpBeInteractionMessenger* PhysicsMessenger;

  //quasi-elastic scattering was added as a separate process in old 
  //MiniBooNE mc, here we use quasielastic modeling included in g4 inelastic
  //models outside 7.5<E<8.5GeV region
  //In energy region where BooNEpBe model is applicable we still have to do 
  //quasi elastic scattering
  BooNEHadronQuasiElasticModel fQuasiElasticModel;

  NuBeamOutput * fRecords;

  G4HadFinalState * theFinalState;

  G4double GetFeynmanX( G4double p, G4double th, G4double m, G4double pbeam, G4double mbeam );

  G4double ProtonXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double NeutronXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiPlusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiMinusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KPlusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KMinusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KZeroLongXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];

  G4double ProtonXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double NeutronXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiPlusXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiMinusXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KPlusXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KMinusXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KZeroLongXSecNoWgtArray[kNProtonMomentumBins][kNPzBins][kNPtBins];

  G4double fPrimaryEnergy;

  //  G4int protonpBins;
  //  const G4double* protonp;

  // The following arrays are declared static to allow the use of initializers.
  // They are initialized in the "physics" subdirectory

  // MARS
  G4double xSectpBeInel_MARS;
  static G4double totalxSectpBeInel_MARS[kNEBins];
  static G4double totalxSectpAlInel_MARS[kNEBins];
  static G4double totalxSectpFeInel_MARS[kNEBins];

  static G4double ProtonXSecArray_MARS[kNPzBins][kNPtBins];
  static G4double NeutronXSecArray_MARS[kNPzBins][kNPtBins];
  static G4double PiPlusXSecArray_MARS[kNPzBins][kNPtBins];
  static G4double PiMinusXSecArray_MARS[kNPzBins][kNPtBins];
  static G4double KPlusXSecArray_MARS[kNPzBins][kNPtBins];
  static G4double KMinusXSecArray_MARS[kNPzBins][kNPtBins];
  // GFLUKA
  G4double xSectpBeInel_GFLUKA;
  static G4double totalxSectpBeInel_GFLUKA[kNEBins];
  static G4double totalxSectpAlInel_GFLUKA[kNEBins];
  static G4double totalxSectpFeInel_GFLUKA[kNEBins];

  static G4double ProtonXSecArray_GFLUKA[kNPzBins][kNPtBins];
  static G4double NeutronXSecArray_GFLUKA[kNPzBins][kNPtBins];
  static G4double PiPlusXSecArray_GFLUKA[kNPzBins][kNPtBins];
  static G4double PiMinusXSecArray_GFLUKA[kNPzBins][kNPtBins];
  static G4double KPlusXSecArray_GFLUKA[kNPzBins][kNPtBins];
  static G4double KMinusXSecArray_GFLUKA[kNPzBins][kNPtBins];
  static G4double KZeroLongXSecArray_GFLUKA[kNPzBins][kNPtBins];
  // ZGS
  G4double xSectpBeInel_ZGS;
  static G4double PiPlusXSecArray_ZGS[kNPzBins][kNPtBins];
  static G4double PiMinusXSecArray_ZGS[kNPzBins][kNPtBins];
  // SWPar
  G4double xSectpBeInel_SWPar;
  G4double PiPlusXSecArray_SWPar[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiMinusXSecArray_SWPar[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KPlusXSecArray_SWPar[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KZeroLongXSecArray_SWPar[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KPlusXSecArray_FSPar[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KZeroLongXSecArray_FSPar[kNProtonMomentumBins][kNPzBins][kNPtBins];

  // G4DEFAULT
  G4double xSectpBeInel_G4DEFAULT;
  static G4double ProtonXSecArray_G4DEFAULT[kNPzBins][kNPtBins];
  static G4double NeutronXSecArray_G4DEFAULT[kNPzBins][kNPtBins];
  static G4double PiPlusXSecArray_G4DEFAULT[kNPzBins][kNPtBins];
  static G4double PiMinusXSecArray_G4DEFAULT[kNPzBins][kNPtBins];
  static G4double KPlusXSecArray_G4DEFAULT[kNPzBins][kNPtBins];
  static G4double KMinusXSecArray_G4DEFAULT[kNPzBins][kNPtBins];
  static G4double KZeroLongXSecArray_G4DEFAULT[kNPzBins][kNPtBins];


  G4ThreeVector secondaryMomentum;
  G4bool acceptSecondaryMomentum;
  G4double XSecValue;  G4double trialPx, trialPy, trialPz;
  G4ThreeVector trialP, p0;
  G4DynamicParticle* aProton;
  G4DynamicParticle* aNeutron;
  G4DynamicParticle* aPiPlus;
  G4DynamicParticle* aPiMinus;
  G4DynamicParticle* aKPlus;
  G4DynamicParticle* aKMinus;
  G4DynamicParticle* aKZeroLong;
  //

  G4String PiPlusPhysicsModel, PiMinusPhysicsModel;
  G4String KaonPlusPhysicsModel, KaonMinusPhysicsModel;
  G4String KaonZeroLongPhysicsModel, ProtonPhysicsModel;
  G4String NeutronPhysicsModel;

  G4bool NoBeamPions;

  G4double SWPiPlusPar1; G4double SWPiPlusPar2; G4double SWPiPlusPar3;
  G4double SWPiPlusPar4; G4double SWPiPlusPar5; G4double SWPiPlusPar6;
  G4double SWPiPlusPar7; G4double SWPiPlusPar8;
  G4double SWPiMinusPar1; G4double SWPiMinusPar2; G4double SWPiMinusPar3;
  G4double SWPiMinusPar4; G4double SWPiMinusPar5; G4double SWPiMinusPar6;
  G4double SWPiMinusPar7; G4double SWPiMinusPar8;
  G4double SWKaonPlusPar1; G4double SWKaonPlusPar2; G4double SWKaonPlusPar3;
  G4double SWKaonPlusPar4; G4double SWKaonPlusPar5; G4double SWKaonPlusPar6;
  G4double SWKaonPlusPar7; G4double SWKaonPlusPar8; G4double SWKaonPlusPar9;
  G4double SWKaonZeroLongPar1; G4double SWKaonZeroLongPar2; G4double SWKaonZeroLongPar3;
  G4double SWKaonZeroLongPar4; G4double SWKaonZeroLongPar5; G4double SWKaonZeroLongPar6;
  G4double SWKaonZeroLongPar7; G4double SWKaonZeroLongPar8; G4double SWKaonZeroLongPar9;
  G4double FSKaonPlusPar1; G4double FSKaonPlusPar2; G4double FSKaonPlusPar3;
  G4double FSKaonPlusPar4; G4double FSKaonPlusPar5; G4double FSKaonPlusPar6;
  G4double FSKaonPlusPar7; G4double FSKaonPlusPar8;
  G4double FSKaonZeroLongPar1; G4double FSKaonZeroLongPar2; G4double FSKaonZeroLongPar3;
  G4double FSKaonZeroLongPar4; G4double FSKaonZeroLongPar5; G4double FSKaonZeroLongPar6;
  G4double FSKaonZeroLongPar7; G4double FSKaonZeroLongPar8;

  G4bool UseBeToAScaling;
  G4double BeToAScalingPar0; G4double BeToAScalingPar1; G4double BeToAScalingPar2;

  G4bool PhysicsVerbose;
  G4double ZGSpBeInelasticXSec, SWParpBeInelasticXSec;

  G4int funcID[15];
  G4double rwgtParams[15][10];
  
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
  
  // kinetic energy (GeV) for incident protons.
  const G4double fProtonKE[kNEBins] = { 
    1., 2., 3., 4., 5., 6., 7., 
    8., 9., 10., 11., 12., 13., 14., 
    15., 16., 17., 18., 19., 20. 
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

#endif
