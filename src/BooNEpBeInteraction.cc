#include "BooNEpBeInteractionMessenger.hh"
#include "BooNEpBeInteraction.hh"
#include "XSecArrayHolder.hh"
#include "BooNEHadronPhysics.hh"
#include "NuBeamRunManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"
#include "G4ProcessManager.hh"
#include "G4UImanager.hh"

//Initialization of static data arrays:
#include "MARSData.hh"
#include "GFLUKAData.hh"
#include "ZGSData.hh"
#include "SWParData.hh"

// include math package
#include <math.h>

BooNEpBeInteraction::BooNEpBeInteraction()
  :G4HadronicInteraction("BooNEpBeInteraction"), fIsQE(false)
{
  PhysicsMessenger = new BooNEpBeInteractionMessenger(this);
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  UI->ApplyCommand("/boone/physics/BooNEpBeModelEnergyRange");
  UI->ApplyCommand("/boone/physics/SWPiPlusPar");
  UI->ApplyCommand("/boone/physics/SWPiMinusPar");
  UI->ApplyCommand("/boone/physics/SWKaonPlusPar");
  UI->ApplyCommand("/boone/physics/SWKaonZeroLongPar");
  UI->ApplyCommand("/boone/physics/FSKaonPlusPar");
  UI->ApplyCommand("/boone/physics/FSKaonZeroLongPar");

  UI->ApplyCommand("boone/physics/quasielastic/proQuasiElasticPar");
  UI->ApplyCommand("boone/physics/quasielastic/neuQuasiElasticPar");
  UI->ApplyCommand("boone/physics/quasielastic/pipQuasiElasticPar");
  UI->ApplyCommand("boone/physics/quasielastic/pimQuasiElasticPar");

  const NuBeamPrimaryGeneratorAction* nbpga=dynamic_cast<const NuBeamPrimaryGeneratorAction*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  fPrimaryEnergy = nbpga->GetParticleGun()->GetParticleEnergy();  
  FillOtherCrossSections(fPrimaryEnergy);
}

BooNEpBeInteraction::~BooNEpBeInteraction()
{
}

void BooNEpBeInteraction::SetSWPiPlusPar(
					 G4double val1, G4double val2,
					 G4double val3, G4double val4,
					 G4double val5, G4double val6,
					 G4double val7, G4double val8)
{
  SWPiPlusPar1 = val1;
  SWPiPlusPar2 = val2;
  SWPiPlusPar3 = val3;
  SWPiPlusPar4 = val4;
  SWPiPlusPar5 = val5;
  SWPiPlusPar6 = val6;
  SWPiPlusPar7 = val7;
  SWPiPlusPar8 = val8;
}

G4double BooNEpBeInteraction::GetSWPiPlusPar(G4int parNo)
{
  if (parNo == 1) return SWPiPlusPar1;
  else if (parNo == 2) return SWPiPlusPar2;
  else if (parNo == 3) return SWPiPlusPar3;
  else if (parNo == 4) return SWPiPlusPar4;
  else if (parNo == 5) return SWPiPlusPar5;
  else if (parNo == 6) return SWPiPlusPar6;
  else if (parNo == 7) return SWPiPlusPar7;
  else if (parNo == 8) return SWPiPlusPar8;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}

// ----------------------------------------


void BooNEpBeInteraction::SetSWPiMinusPar(
					 G4double val1, G4double val2,
					 G4double val3, G4double val4,
					 G4double val5, G4double val6,
					 G4double val7, G4double val8)
{
  SWPiMinusPar1 = val1;
  SWPiMinusPar2 = val2;
  SWPiMinusPar3 = val3;
  SWPiMinusPar4 = val4;
  SWPiMinusPar5 = val5;
  SWPiMinusPar6 = val6;
  SWPiMinusPar7 = val7;
  SWPiMinusPar8 = val8;
}

G4double BooNEpBeInteraction::GetSWPiMinusPar(G4int parNo)
{
  if (parNo == 1) return SWPiMinusPar1;
  else if (parNo == 2) return SWPiMinusPar2;
  else if (parNo == 3) return SWPiMinusPar3;
  else if (parNo == 4) return SWPiMinusPar4;
  else if (parNo == 5) return SWPiMinusPar5;
  else if (parNo == 6) return SWPiMinusPar6;
  else if (parNo == 7) return SWPiMinusPar7;
  else if (parNo == 8) return SWPiMinusPar8;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}

// ----------------------------------------

void BooNEpBeInteraction::SetSWKaonPlusPar(
					 G4double val1, G4double val2,
					 G4double val3, G4double val4,
					 G4double val5, G4double val6,
					 G4double val7, G4double val8,
					 G4double val9)
{
  SWKaonPlusPar1 = val1;
  SWKaonPlusPar2 = val2;
  SWKaonPlusPar3 = val3;
  SWKaonPlusPar4 = val4;
  SWKaonPlusPar5 = val5;
  SWKaonPlusPar6 = val6;
  SWKaonPlusPar7 = val7;
  SWKaonPlusPar8 = val8;
  SWKaonPlusPar9 = val9;
}

G4double BooNEpBeInteraction::GetSWKaonPlusPar(G4int parNo)
{
  if (parNo == 1) return SWKaonPlusPar1;
  else if (parNo == 2) return SWKaonPlusPar2;
  else if (parNo == 3) return SWKaonPlusPar3;
  else if (parNo == 4) return SWKaonPlusPar4;
  else if (parNo == 5) return SWKaonPlusPar5;
  else if (parNo == 6) return SWKaonPlusPar6;
  else if (parNo == 7) return SWKaonPlusPar7;
  else if (parNo == 8) return SWKaonPlusPar8;
  else if (parNo == 9) return SWKaonPlusPar9;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}

// ----------------------------------------

void BooNEpBeInteraction::SetSWKaonZeroLongPar(
					 G4double val1, G4double val2,
					 G4double val3, G4double val4,
					 G4double val5, G4double val6,
					 G4double val7, G4double val8,
					 G4double val9)
{
  SWKaonZeroLongPar1 = val1;
  SWKaonZeroLongPar2 = val2;
  SWKaonZeroLongPar3 = val3;
  SWKaonZeroLongPar4 = val4;
  SWKaonZeroLongPar5 = val5;
  SWKaonZeroLongPar6 = val6;
  SWKaonZeroLongPar7 = val7;
  SWKaonZeroLongPar8 = val8;
  SWKaonZeroLongPar9 = val9;
}

G4double BooNEpBeInteraction::GetSWKaonZeroLongPar(G4int parNo)
{
  if (parNo == 1) return SWKaonZeroLongPar1;
  else if (parNo == 2) return SWKaonZeroLongPar2;
  else if (parNo == 3) return SWKaonZeroLongPar3;
  else if (parNo == 4) return SWKaonZeroLongPar4;
  else if (parNo == 5) return SWKaonZeroLongPar5;
  else if (parNo == 6) return SWKaonZeroLongPar6;
  else if (parNo == 7) return SWKaonZeroLongPar7;
  else if (parNo == 8) return SWKaonZeroLongPar8;
  else if (parNo == 9) return SWKaonZeroLongPar9;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}

// ----------------------------------------

void BooNEpBeInteraction::SetFSKaonPlusPar(
					 G4double val1, G4double val2,
					 G4double val3, G4double val4,
					 G4double val5, G4double val6,
					 G4double val7, G4double val8)
{
  FSKaonPlusPar1 = val1;
  FSKaonPlusPar2 = val2;
  FSKaonPlusPar3 = val3;
  FSKaonPlusPar4 = val4;
  FSKaonPlusPar5 = val5;
  FSKaonPlusPar6 = val6;
  FSKaonPlusPar7 = val7;
  FSKaonPlusPar8 = val8;
}

G4double BooNEpBeInteraction::GetFSKaonPlusPar(G4int parNo)
{
  if (parNo == 1) return FSKaonPlusPar1;
  else if (parNo == 2) return FSKaonPlusPar2;
  else if (parNo == 3) return FSKaonPlusPar3;
  else if (parNo == 4) return FSKaonPlusPar4;
  else if (parNo == 5) return FSKaonPlusPar5;
  else if (parNo == 6) return FSKaonPlusPar6;
  else if (parNo == 7) return FSKaonPlusPar7;
  else if (parNo == 8) return FSKaonPlusPar8;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}

// ----------------------------------------

void BooNEpBeInteraction::SetFSKaonZeroLongPar(
					 G4double val1, G4double val2,
					 G4double val3, G4double val4,
					 G4double val5, G4double val6,
					 G4double val7, G4double val8)
{
  FSKaonZeroLongPar1 = val1;
  FSKaonZeroLongPar2 = val2;
  FSKaonZeroLongPar3 = val3;
  FSKaonZeroLongPar4 = val4;
  FSKaonZeroLongPar5 = val5;
  FSKaonZeroLongPar6 = val6;
  FSKaonZeroLongPar7 = val7;
  FSKaonZeroLongPar8 = val8;
}

G4double BooNEpBeInteraction::GetFSKaonZeroLongPar(G4int parNo)
{
  if (parNo == 1) return FSKaonZeroLongPar1;
  else if (parNo == 2) return FSKaonZeroLongPar2;
  else if (parNo == 3) return FSKaonZeroLongPar3;
  else if (parNo == 4) return FSKaonZeroLongPar4;
  else if (parNo == 5) return FSKaonZeroLongPar5;
  else if (parNo == 6) return FSKaonZeroLongPar6;
  else if (parNo == 7) return FSKaonZeroLongPar7;
  else if (parNo == 8) return FSKaonZeroLongPar8;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}

// -------------------------------------------------------------------------


// set pi+ production physics model

void
BooNEpBeInteraction::SetPiPlusPhysicsModel(G4String val)
{

  PiPlusPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;


  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;
  G4double ProtonMass2 = ProtonMass * ProtonMass;
  G4double PionMass = G4PionPlus::PionPlusDefinition()->GetPDGMass()/CLHEP::GeV;
  
  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );

  // MARS physics
  if (PiPlusPhysicsModel == "MARS") 
  {
    G4cout << "PiPlus Physics Model is MARS" << G4endl;
    // Note: to get the correct MARS multiplicity for pi+, need to scale by the
    // ratio of specfied inelastic cross-section divided by the MARS inelastic
    // cross-section
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  PiPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiPlusXSecArray_MARS[ipz][ipt];
	  PiPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiPlusXSecNoWgtArray,iprotonp,ipz,ipt,8);
	}
      }
    }
  }

  // GFLUKA physics
  else if (PiPlusPhysicsModel == "GFLUKA") 
  {
    G4cout << "PiPlusPhysics Model is GFLUKA" << G4endl;
    // Note: to get the correct GFLUKA multiplicity for pi+, need to scale by the
    // ratio of specfied inelastic cross-section divided by the GFLUKA inelastic
    // cross-section
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  PiPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiPlusXSecArray_GFLUKA[ipz][ipt];
	  PiPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiPlusXSecNoWgtArray,iprotonp,ipz,ipt,8);
	}
      }
    }
  }

  // ZGS physics 
  else if (PiPlusPhysicsModel == "ZGS") 
  {
    G4cout << "PiPlus Physics Model is ZGS" << G4endl;
    // Note: to get the correct ZGS multiplicity for pi+, need to scale by the
    // ratio of specfied inelastic cross-section divided by the ZGS inelastic
    // cross-section
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    { 
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  PiPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiPlusXSecArray_ZGS[ipz][ipt];
	  PiPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiPlusXSecNoWgtArray,iprotonp,ipz,ipt,8);
	}
      }
    }
  }

  // SWPar model 
  else if (PiPlusPhysicsModel == "SWPar") 
  {    
    G4cout << "PiPlus Physics Model is SWPar, with the following choice of parameters:" << G4endl;
    G4cout << "SWPiPlusPar[1-8] = " << 
      SWPiPlusPar1 << ", " << SWPiPlusPar2 << ", " <<
      SWPiPlusPar3 << ", " << SWPiPlusPar4 << ", " <<
      SWPiPlusPar5 << ", " << SWPiPlusPar6 << ", " <<
      SWPiPlusPar7 << ", " << SWPiPlusPar8 << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  EKinetic = sqrt (fProtonMomentumBins[iprotonp]*fProtonMomentumBins[iprotonp]
			   + ProtonMass2) - ProtonMass;
	  pbeam = fProtonMomentumBins[iprotonp];
	  if( scale )
	  {
	    xFey = GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
				     + fPtVec[ipt]*fPtVec[ipt]),
				atan(fPtVec[ipt]/fPzVec[ipz]),
				PionMass, pbeam, ProtonMass );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }

	  G4double p = sqrt(fPzVec[ipz]*fPzVec[ipz]+fPtVec[ipt]*fPtVec[ipt]);
	  G4double theta = atan(fPtVec[ipt]/fPzVec[ipz]);
	  G4double first = SWPiPlusPar1*(pow(p,SWPiPlusPar2))*
	    (1.-(p/(pbeam-1.)));
	  G4double arg1 = -1.*SWPiPlusPar3*(pow(p,SWPiPlusPar4))/
	    pow(pbeam,SWPiPlusPar5);
	  G4double arg2 = -1.*SWPiPlusPar6*theta*
	    (p-SWPiPlusPar7*pbeam*(pow(cos(theta),SWPiPlusPar8)));
	  G4double arg = arg1+arg2;
	  G4double second = exp(arg);
	  G4double jacobian = 2.*M_PI*fPtVec[ipt]/(p*p);
	  G4double valXSec = jacobian * (first*second);
	  if (valXSec >= 0.) {
	    PiPlusXSecArray_SWPar[iprotonp][ipz][ipt] = valXSec;
	  } else {
	    PiPlusXSecArray_SWPar[iprotonp][ipz][ipt] = 0.;
	  }
	  PiPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiPlusXSecArray_SWPar[iprotonp][ipz][ipt];
	  PiPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiPlusXSecNoWgtArray,iprotonp,ipz,ipt,8);
	}
      }
    }
  }


  // no PiPlus Physics model
  else {
    G4cout << "PiPlus Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------


// set pi- production physics model

void
BooNEpBeInteraction::SetPiMinusPhysicsModel(G4String val)
{

  PiMinusPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;


  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;
  G4double ProtonMass2 = ProtonMass * ProtonMass;
  G4double PionMass = G4PionMinus::PionMinusDefinition()->GetPDGMass()/CLHEP::GeV;

  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );


  // MARS physics 
  if (PiMinusPhysicsModel == "MARS") 
  {
    G4cout << "PiMinus Physics Model is MARS" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  PiMinusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiMinusXSecArray_MARS[ipz][ipt];
	  PiMinusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiMinusXSecNoWgtArray,iprotonp,ipz,ipt,9);
	}
      }
    }
  }


  // GFLUKA physics
  else if (PiMinusPhysicsModel == "GFLUKA") 
  {
    G4cout << "PiMinusPhysics Model is GFLUKA" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  PiMinusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiMinusXSecArray_GFLUKA[ipz][ipt];
	  PiMinusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiMinusXSecNoWgtArray,iprotonp,ipz,ipt,9);
	}
      }
    }
  }


  // ZGS physics 
  else if (PiMinusPhysicsModel == "ZGS") 
  {
    G4cout << "PiMinus Physics Model is ZGS" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  PiMinusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiMinusXSecArray_ZGS[ipz][ipt];
	  PiMinusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiMinusXSecNoWgtArray,iprotonp,ipz,ipt,9);
	}
      }
    }
  }


  // SWPar model 
  else if (PiMinusPhysicsModel == "SWPar") {
    G4cout << "PiMinus Physics Model is SWPar, with the following choice of parameters:" << G4endl;
    G4cout << "SWPiMinusPar[1-8] = " << 
      SWPiMinusPar1 << ", " << SWPiMinusPar2 << ", " <<
      SWPiMinusPar3 << ", " << SWPiMinusPar4 << ", " <<
      SWPiMinusPar5 << ", " << SWPiMinusPar6 << ", " <<
      SWPiMinusPar7 << ", " << SWPiMinusPar8 << G4endl;

    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  EKinetic = sqrt (fProtonMomentumBins[iprotonp]*fProtonMomentumBins[iprotonp]
			   + ProtonMass2) - ProtonMass;
	  pbeam = fProtonMomentumBins[iprotonp];
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      PionMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }

	  G4double p = sqrt(fPzVec[ipz]*fPzVec[ipz]+fPtVec[ipt]*fPtVec[ipt]);
	  G4double theta = atan(fPtVec[ipt]/fPzVec[ipz]);
	  G4double first = SWPiMinusPar1*(pow(p,SWPiMinusPar2))*
	    (1.-(p/(pbeam-1.)));
	  G4double arg1 = -1.*SWPiMinusPar3*(pow(p,SWPiMinusPar4))/
	    pow(pbeam,SWPiMinusPar5);
	  G4double arg2 = -1.*SWPiMinusPar6*theta*
	    (p-SWPiMinusPar7*pbeam*(pow(cos(theta),SWPiMinusPar8)));
	  G4double arg = arg1+arg2;
	  G4double second = exp(arg);
	  G4double jacobian = 2.*M_PI*fPtVec[ipt]/(p*p);
	  G4double valXSec = jacobian * (first*second);
	  if (valXSec >= 0.) {
	    PiMinusXSecArray_SWPar[iprotonp][ipz][ipt] = valXSec;
	  } else {
	    PiMinusXSecArray_SWPar[iprotonp][ipz][ipt] = 0.;
	  }
	  PiMinusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*PiMinusXSecArray_SWPar[iprotonp][ipz][ipt];
	  PiMinusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(PiMinusXSecNoWgtArray,iprotonp,ipz,ipt,9);
	}
      }
    }
  }
  
  // no PiMinus Physics model
  else {
    G4cout << "PiMinus Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------


// set k+ production physics model

void
BooNEpBeInteraction::SetKaonPlusPhysicsModel(G4String val)
{

  KaonPlusPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;


  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;
  G4double ProtonMass2 = ProtonMass * ProtonMass;
  G4double KaonMass = G4KaonPlus::KaonPlusDefinition()->GetPDGMass()/CLHEP::GeV;
  G4double LambdaMass = G4Lambda::LambdaDefinition()->GetPDGMass()/CLHEP::GeV;
  
  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );


  // MARS physics 
  if (KaonPlusPhysicsModel == "MARS") 
  {
    G4cout << "KaonPlus Physics Model is MARS" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  KPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KPlusXSecArray_MARS[ipz][ipt];
	  KPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KPlusXSecNoWgtArray,iprotonp,ipz,ipt,11);
	}
      }
    }
  }


  // GFLUKA physics
  else if (KaonPlusPhysicsModel == "GFLUKA") 
  {
    G4cout << "KaonPlusPhysics Model is GFLUKA" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  KPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KPlusXSecArray_GFLUKA[ipz][ipt];
	  KPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KPlusXSecNoWgtArray,iprotonp,ipz,ipt,11);
	}
      }
    }
  }


  // SWPar model 
  else if (KaonPlusPhysicsModel == "SWPar") 
  {
    G4cout << 
      "KaonPlus Physics Model is SWPar, with the following choice of parameters:"
	   << G4endl;
    G4cout << "SWKaonPlusPar[1-9] = " << 
      SWKaonPlusPar1 << ", " << SWKaonPlusPar2 << ", " <<
      SWKaonPlusPar3 << ", " << SWKaonPlusPar4 << ", " <<
      SWKaonPlusPar5 << ", " << SWKaonPlusPar6 << ", " <<
      SWKaonPlusPar7 << ", " << SWKaonPlusPar8 << ", " <<
      SWKaonPlusPar9 << G4endl;

    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  EKinetic = sqrt (fProtonMomentumBins[iprotonp]*fProtonMomentumBins[iprotonp]
			   + ProtonMass2) - ProtonMass;
	  pbeam = fProtonMomentumBins[iprotonp];
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				    KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  G4double p = sqrt(fPzVec[ipz]*fPzVec[ipz]+fPtVec[ipt]*fPtVec[ipt]);
	  G4double theta = atan(fPtVec[ipt]/fPzVec[ipz]);
	  G4double first = SWKaonPlusPar1*(pow(p,SWKaonPlusPar2))*
	    (1.-(p/(pbeam-SWKaonPlusPar9)));
	  G4double arg1 = -1.*SWKaonPlusPar3*(pow(p,SWKaonPlusPar4))/
	    pow(pbeam,SWKaonPlusPar5);
	  G4double arg2 = -1.*SWKaonPlusPar6*theta*
	    (p-SWKaonPlusPar7*pbeam*(pow(cos(theta),SWKaonPlusPar8)));
	  G4double arg = arg1+arg2;
	  G4double second = exp(arg);
	  G4double jacobian = 2.*M_PI*fPtVec[ipt]/(p*p);
	  G4double valXSec = jacobian * (first*second);
	  if (valXSec >= 0.) {
	    KPlusXSecArray_SWPar[iprotonp][ipz][ipt] = valXSec;
	  } else {
	    KPlusXSecArray_SWPar[iprotonp][ipz][ipt] = 0.;
	  }
	  KPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KPlusXSecArray_SWPar[iprotonp][ipz][ipt];
	  KPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KPlusXSecNoWgtArray,iprotonp,ipz,ipt,11);
	}
      }
    }
  }
  
  // FScal model 
  else if (KaonPlusPhysicsModel == "FScal") 
  {
    G4cout << 
      "KaonPlus Physics Model is FScal, with the following choice of parameters:"
	   << G4endl;
    G4cout << "FSKaonPlusPar[1-8] = " << 
      FSKaonPlusPar1 << ", " << FSKaonPlusPar2 << ", " <<
      FSKaonPlusPar3 << ", " << FSKaonPlusPar4 << ", " <<
      FSKaonPlusPar5 << ", " << FSKaonPlusPar6 << ", " <<
      FSKaonPlusPar7 << ", " << FSKaonPlusPar8 << ", " << G4endl;

    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  pbeam = fProtonMomentumBins[iprotonp];
	  G4double ebeam = sqrt(pbeam*pbeam + ProtonMass2);
	  EKinetic = ebeam - ProtonMass;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  G4double pKaonLab = sqrt(fPzVec[ipz]*fPzVec[ipz]+fPtVec[ipt]*fPtVec[ipt]);
	  G4double eKaonLab = sqrt(pKaonLab*pKaonLab + KaonMass*KaonMass);
	  G4double sMand = 2.*ProtonMass*(ProtonMass + ebeam);
	  G4double eKaonCMMax = (sMand + KaonMass*KaonMass - pow(LambdaMass + ProtonMass,2)) / (2. * sqrt(sMand));
	  G4double pKaonCMMax = sqrt(eKaonCMMax*eKaonCMMax - KaonMass*KaonMass);
	  G4double pGamma = sqrt(sMand) / (2. * ProtonMass);
	  G4double pBeta = sqrt(1. - 1./(pGamma * pGamma));
	  G4double pzKaonCM = pGamma*(fPzVec[ipz]-pBeta*eKaonLab);
	  xFey = pzKaonCM / pKaonCMMax;

	  G4double first = pKaonLab*pKaonLab/eKaonLab*FSKaonPlusPar1*pow(1.-fabs(xFey),FSKaonPlusPar8);
	  G4double second = exp(-FSKaonPlusPar3*pow(fabs(xFey),FSKaonPlusPar4) - FSKaonPlusPar7*pow(fabs(fPtVec[ipt]*xFey),FSKaonPlusPar6) - FSKaonPlusPar2*fPtVec[ipt] - FSKaonPlusPar5*fPtVec[ipt]*fPtVec[ipt]);

	  G4double jacobian = 2.*M_PI*fPtVec[ipt]/(pKaonLab*pKaonLab);
	  G4double valXSec = jacobian * (first*second);
	  if ((valXSec >= 0.) && (fabs(xFey) <= 1.) && (eKaonCMMax >= 0.)) {
	    KPlusXSecArray_FSPar[iprotonp][ipz][ipt] = valXSec;
	  } else {
	    KPlusXSecArray_FSPar[iprotonp][ipz][ipt] = 0.;
	  }
	  KPlusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KPlusXSecArray_FSPar[iprotonp][ipz][ipt];
	  KPlusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KPlusXSecNoWgtArray,iprotonp,ipz,ipt,11);
	}
      }
    }
  }
  

  // no KaonPlus Physics model
  else {
    G4cout << "KaonPlus Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------


// set k- production physics model

void
BooNEpBeInteraction::SetKaonMinusPhysicsModel(G4String val)
{

  KaonMinusPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;

  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;
  G4double KaonMass = G4KaonMinus::KaonMinusDefinition()->GetPDGMass()/CLHEP::GeV;

  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );


  // MARS physics 
  if (KaonMinusPhysicsModel == "MARS") 
  {
    G4cout << "KaonMinus Physics Model is MARS" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  KMinusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KMinusXSecArray_MARS[ipz][ipt];
	  KMinusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KMinusXSecNoWgtArray,iprotonp,ipz,ipt,12);
	}
      }
    }
  }


  // GFLUKA physics
  else if (KaonMinusPhysicsModel == "GFLUKA") 
  {
    G4cout << "KaonMinusPhysics Model is GFLUKA" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  KMinusXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KMinusXSecArray_GFLUKA[ipz][ipt];
	  KMinusXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KMinusXSecNoWgtArray,iprotonp,ipz,ipt,12);
	}
      }
    }
  }


  // no KaonMinus Physics model
  else {
    G4cout << "KaonMinus Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------


// set k0L production physics model

void
BooNEpBeInteraction::SetKaonZeroLongPhysicsModel(G4String val)
{

  KaonZeroLongPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;


  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;
  G4double ProtonMass2 = ProtonMass * ProtonMass;
  G4double KaonMass = G4KaonZeroLong::KaonZeroLongDefinition()->GetPDGMass()/CLHEP::GeV;
  G4double SigmaMass = G4SigmaPlus::SigmaPlusDefinition()->GetPDGMass()/CLHEP::GeV;

  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );


  // MARS physics 
  if (KaonZeroLongPhysicsModel == "MARS") 
  {
    G4cout << "KaonZeroLong Physics Model is MARS" << G4endl;
    // MARS does not store KZeroLong information.
    // Temporary kludge: for KZeroLong, use MARS KPlus cross-section,
    // times the ratio of KZeroLong to KPlus in GFLUKA, bin by bin
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  if (KPlusXSecArray_GFLUKA[ipz][ipt]>0.) {
	    KZeroLongXSecNoWgtArray[iprotonp][ipz][ipt] =
	      xSectRatio*KPlusXSecArray_MARS[ipz][ipt]*
	      (KZeroLongXSecArray_GFLUKA[ipz][ipt]/
	       KPlusXSecArray_GFLUKA[ipz][ipt]);
	    KZeroLongXSecArray[iprotonp][ipz][ipt]
	      = GetRwgtXSec(KZeroLongXSecNoWgtArray,iprotonp,ipz,ipt,10);
	  } else {
	    KZeroLongXSecNoWgtArray[iprotonp][ipz][ipt] = 0.;
	    KZeroLongXSecArray[iprotonp][ipz][ipt] = 0.;
	  }
	}
      }
    }
  }


  // GFLUKA physics
  else if (KaonZeroLongPhysicsModel == "GFLUKA") 
  {
    G4cout << "KaonZeroLongPhysics Model is GFLUKA" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  KZeroLongXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KZeroLongXSecArray_GFLUKA[ipz][ipt];
	  KZeroLongXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KZeroLongXSecNoWgtArray,iprotonp,ipz,ipt,10);
	}
      }
    }
  }

  // SWPar model 
  else if (KaonZeroLongPhysicsModel == "SWPar") 
  {
    G4cout << 
      "KaonZeroLong Physics Model is SWPar, with the following choice of parameters:"
	   << G4endl;
    G4cout << "SWKaonZeroLongPar[1-9] = " << 
      SWKaonZeroLongPar1 << ", " << SWKaonZeroLongPar2 << ", " <<
      SWKaonZeroLongPar3 << ", " << SWKaonZeroLongPar4 << ", " <<
      SWKaonZeroLongPar5 << ", " << SWKaonZeroLongPar6 << ", " <<
      SWKaonZeroLongPar7 << ", " << SWKaonZeroLongPar8 << ", " <<
      SWKaonZeroLongPar9 << G4endl;

    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  EKinetic = sqrt (fProtonMomentumBins[iprotonp]*fProtonMomentumBins[iprotonp]
			   + ProtonMass2) - ProtonMass;
	  pbeam = fProtonMomentumBins[iprotonp];
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				    KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  G4double p = sqrt(fPzVec[ipz]*fPzVec[ipz]+fPtVec[ipt]*fPtVec[ipt]);
	  G4double theta = atan(fPtVec[ipt]/fPzVec[ipz]);
	  G4double first = SWKaonZeroLongPar1*(pow(p,SWKaonZeroLongPar2))*
	    (1.-(p/(pbeam-SWKaonZeroLongPar9)));
	  G4double arg1 = -1.*SWKaonZeroLongPar3*(pow(p,SWKaonZeroLongPar4))/
	    pow(pbeam,SWKaonZeroLongPar5);
	  G4double arg2 = -1.*SWKaonZeroLongPar6*theta*
	    (p-SWKaonZeroLongPar7*pbeam*(pow(cos(theta),SWKaonZeroLongPar8)));
	  G4double arg = arg1+arg2;
	  G4double second = exp(arg);
	  G4double jacobian = 2.*M_PI*fPtVec[ipt]/(p*p);
	  G4double valXSec = jacobian * (first*second);
	  if (valXSec >= 0.) {
	    KZeroLongXSecArray_SWPar[iprotonp][ipz][ipt] = valXSec;
	  } else {
	    KZeroLongXSecArray_SWPar[iprotonp][ipz][ipt] = 0.;
	  }
	  KZeroLongXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KZeroLongXSecArray_SWPar[iprotonp][ipz][ipt];
	  KZeroLongXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KZeroLongXSecNoWgtArray,iprotonp,ipz,ipt,10);
	}
      }
    }
  }

  // FScal model 
  else if (KaonZeroLongPhysicsModel == "FScal") 
  {
    G4cout << 
      "KaonZeroLong Physics Model is FScal, with the following choice of parameters:"
	   << G4endl;
    G4cout << "FSKaonZeroLongPar[1-8] = " << 
      FSKaonZeroLongPar1 << ", " << FSKaonZeroLongPar2 << ", " <<
      FSKaonZeroLongPar3 << ", " << FSKaonZeroLongPar4 << ", " <<
      FSKaonZeroLongPar5 << ", " << FSKaonZeroLongPar6 << ", " <<
      FSKaonZeroLongPar7 << ", " << FSKaonZeroLongPar8 << ", " << G4endl;

    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  pbeam = fProtonMomentumBins[iprotonp];
	  G4double ebeam = sqrt(pbeam*pbeam + ProtonMass2);
	  EKinetic = ebeam - ProtonMass;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      KaonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  G4double pKaonLab = sqrt(fPzVec[ipz]*fPzVec[ipz]+fPtVec[ipt]*fPtVec[ipt]);
	  G4double eKaonLab = sqrt(pKaonLab*pKaonLab + KaonMass*KaonMass);
	  G4double sMand = 2.*ProtonMass*(ProtonMass + ebeam);
	  G4double eKaonCMMax = (sMand + KaonMass*KaonMass - pow(SigmaMass + ProtonMass,2)) / (2. * sqrt(sMand));
	  G4double pKaonCMMax = sqrt(eKaonCMMax*eKaonCMMax - KaonMass*KaonMass);
	  G4double pGamma = sqrt(sMand) / (2. * ProtonMass);
	  G4double pBeta = sqrt(1. - 1./(pGamma * pGamma));
	  G4double pzKaonCM = pGamma*(fPzVec[ipz]-pBeta*eKaonLab);
	  xFey = pzKaonCM / pKaonCMMax;

	  G4double first = pKaonLab*pKaonLab/eKaonLab*FSKaonZeroLongPar1*pow(1.-fabs(xFey),FSKaonZeroLongPar8);
	  G4double second = exp(-FSKaonZeroLongPar3*pow(fabs(xFey),FSKaonZeroLongPar4) - FSKaonZeroLongPar7*pow(fabs(fPtVec[ipt]*xFey),FSKaonZeroLongPar6) - FSKaonZeroLongPar2*fPtVec[ipt] - FSKaonZeroLongPar5*fPtVec[ipt]*fPtVec[ipt]);

	  G4double jacobian = 2.*M_PI*fPtVec[ipt]/(pKaonLab*pKaonLab);
	  G4double valXSec = jacobian * (first*second);
	  if ((valXSec >= 0.) && (fabs(xFey) <= 1.) && (eKaonCMMax >= 0.)) {
	    KZeroLongXSecArray_FSPar[iprotonp][ipz][ipt] = valXSec;
	  } else {
	    KZeroLongXSecArray_FSPar[iprotonp][ipz][ipt] = 0.;
	  }
	  KZeroLongXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*KZeroLongXSecArray_FSPar[iprotonp][ipz][ipt];
	  KZeroLongXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(KZeroLongXSecNoWgtArray,iprotonp,ipz,ipt,11);
	}
      }
    }
  }

  // no KaonZeroLong Physics model
  else {
    G4cout << "KaonZeroLong Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------


// set proton production physics model

void
BooNEpBeInteraction::SetProtonPhysicsModel(G4String val)
{

  ProtonPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;


  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;

  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );


  // MARS physics 
  if (ProtonPhysicsModel == "MARS") 
  {
    G4cout << "Proton Physics Model is MARS" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      ProtonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  ProtonXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*ProtonXSecArray_MARS[ipz][ipt];
	  ProtonXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(ProtonXSecNoWgtArray,iprotonp,ipz,ipt,14);
	}
      }
    }
  }


  // GFLUKA physics
  else if (ProtonPhysicsModel == "GFLUKA") 
  {
    G4cout << "ProtonPhysics Model is GFLUKA" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      ProtonMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  ProtonXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*ProtonXSecArray_GFLUKA[ipz][ipt];
	  ProtonXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(ProtonXSecNoWgtArray,iprotonp,ipz,ipt,14);
	}
      }
    }
  }


  // no Proton Physics model
  else {
    G4cout << "Proton Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------


// set neutron production physics model

void
BooNEpBeInteraction::SetNeutronPhysicsModel(G4String val)
{

  NeutronPhysicsModel = val;

  G4bool scale = GetUseBeToAScaling();
  G4double c0=0, c1=0, c2=0, xFey=0;
  if( scale )
  {
    c0 = GetBeToAScalingPar(0);
    c1 = GetBeToAScalingPar(1);
    c2 = GetBeToAScalingPar(2);
    G4cout << "A scaling is active !!  ( c0, c1, c2 ) = ( " 
	   << c0 << ", " 
	   << c1 << ", "
	   << c2 << " )" << G4endl;
  }

  G4double xSectRatio = 1.;
  G4Material* targetMat = GetTargetMaterial();
  const G4Element* ele = targetMat->GetElement(0);
  G4double targetA = ele->GetN();
  G4double BeA = 9.01218;
  G4double ratioA = targetA / BeA;
  if( scale )
    G4cout << "                        ( A_targ / A_Be ) = " << targetA << " / " << BeA << " = " << ratioA << G4endl;


  //--Retrieve the Primary Beam Energy and determine the momentum
  G4double EKinetic = fPrimaryEnergy/CLHEP::GeV;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName = "proton" );
  G4double ProtonMass = particle->GetPDGMass()/CLHEP::GeV;
  G4double NeutronMass = G4Neutron::NeutronDefinition()->GetPDGMass()/CLHEP::GeV;

  G4double pbeam = sqrt( (EKinetic + ProtonMass)*(EKinetic + ProtonMass) - ProtonMass*ProtonMass );


  // MARS physics 
  if (NeutronPhysicsModel == "MARS") 
  {
    G4cout << "Neutron Physics Model is MARS" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
          {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      NeutronMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_MARS;
	  NeutronXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*NeutronXSecArray_MARS[ipz][ipt];
	  NeutronXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(NeutronXSecNoWgtArray,iprotonp,ipz,ipt,13);
	}
      }
    }
  }


  // GFLUKA physics
  else if (NeutronPhysicsModel == "GFLUKA") 
  {
    G4cout << "NeutronPhysics Model is GFLUKA" << G4endl;
    for (G4int ipz=0; ipz<kNPzBins; ipz++)
    {
      for (G4int ipt=0; ipt<kNPtBins; ipt++)
      {
	for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
	  xSectRatio = 1.;
	  if( scale )
	  {
	    xFey = fabs( GetFeynmanX( sqrt(fPzVec[ipz]*fPzVec[ipz]
					   + fPtVec[ipt]*fPtVec[ipt]),
				      atan(fPtVec[ipt]/fPzVec[ipz]),
				      NeutronMass, pbeam, ProtonMass ) );
	    xSectRatio *= pow( ratioA, c0 + c1*xFey + c2*xFey*xFey );
	  }
	  xSectRatio *= GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)
	    / xSectpBeInel_GFLUKA;
	  NeutronXSecNoWgtArray[iprotonp][ipz][ipt]
	    = xSectRatio*NeutronXSecArray_GFLUKA[ipz][ipt];
	  NeutronXSecArray[iprotonp][ipz][ipt]
	    = GetRwgtXSec(NeutronXSecNoWgtArray,iprotonp,ipz,ipt,13);
	}
      }
    }
  }


  // no Neutron Physics model
  else {
    G4cout << "Neutron Physics Model is not valid! Results are nonsense" <<
      G4endl;
  }
}


// ------------------------------------------------------------------------

void BooNEpBeInteraction::SetUseBeToAScaling(G4bool valb)
{
  UseBeToAScaling = valb;
}

void BooNEpBeInteraction::SetBeToAScalingPar(G4double c0, G4double c1, G4double c2)
{
  BeToAScalingPar0 = c0;
  BeToAScalingPar1 = c1;
  BeToAScalingPar2 = c2;
}


G4double BooNEpBeInteraction::GetBeToAScalingPar(G4int parNo)
{
  if (parNo == 0) return BeToAScalingPar0;
  else if (parNo == 1) return BeToAScalingPar1;
  else if (parNo == 2) return BeToAScalingPar2;
  else {
    G4cout << "Problems" << G4endl;
    return 0.;
  }
}


// ------------------------------------------------------------------------

void BooNEpBeInteraction::SetPhysicsVerbose(G4bool valb)
{
  PhysicsVerbose = valb;
}

// -------------------------------------------

void BooNEpBeInteraction::SetNoBeamPions(G4bool valb)
{
  G4cout << "/boone/physics/noBeamPions " << valb << G4endl;
  NoBeamPions = valb;
}


G4HadFinalState *
BooNEpBeInteraction::ApplyYourself( const G4HadProjectile &aTrack,
				    G4Nucleus& aNucl )
{
  fIsQE=false;
  static G4bool FirstProtonBeInteraction = true;

  static G4double ProtonMaxXSec;
  static G4double NeutronMaxXSec;
  static G4double PiPlusMaxXSec;
  static G4double PiMinusMaxXSec;
  static G4double KPlusMaxXSec;
  static G4double KMinusMaxXSec;
  static G4double KZeroLongMaxXSec;

  theParticleChange.Clear();

  if (FirstProtonBeInteraction) {
    FirstProtonBeInteraction = false;
    ProtonMaxXSec = GetMax( ProtonXSecArray );
    NeutronMaxXSec = GetMax( NeutronXSecArray );
    PiPlusMaxXSec = GetMax( PiPlusXSecArray );
    PiMinusMaxXSec = GetMax( PiMinusXSecArray );
    KPlusMaxXSec = GetMax( KPlusXSecArray );
    KMinusMaxXSec = GetMax( KMinusXSecArray );
    KZeroLongMaxXSec = GetMax( KZeroLongXSecArray );

    // Also at this point, we need to let the singleton know what xsec it can expect.
    XSecArrayHolder& xsec_holder = XSecArrayHolder::Instance();
    if( !xsec_holder.IsInitialised() ) {
      xsec_holder.LoadXSecData( ProtonXSecArray,
				NeutronXSecArray,
				PiPlusXSecArray,
				PiMinusXSecArray,
				KPlusXSecArray,
				KMinusXSecArray,
				KZeroLongXSecArray,
				ProtonXSecNoWgtArray,
				NeutronXSecNoWgtArray,
				PiPlusXSecNoWgtArray,
				PiMinusXSecNoWgtArray,
				KPlusXSecNoWgtArray,
				KMinusXSecNoWgtArray,
				KZeroLongXSecNoWgtArray );
    }
  }

  if (PhysicsVerbose) {
    G4cout << "Enter BooNEpBeInteraction::ApplyYourself" << G4endl;
    G4cout << "incident proton kinetic energy = " <<
      aTrack.GetKineticEnergy()/CLHEP::GeV << " GeV" << G4endl;
  }

  G4double trackMom = aTrack.GetTotalMomentum();

  theFinalState = NULL;

  // first decide if this is QE or inelastic interaction
  BooNEHadronCrossSections* boonexsec= BooNEHadronCrossSections::Instance();
  G4DynamicParticle* dynPart=new G4DynamicParticle(aTrack.GetDefinition(),
						   aTrack.Get4Momentum());
  G4double qexsec=boonexsec->GetQuasiElasticCrossSection(dynPart,aNucl.GetZ_asInt(),aNucl.GetA_asInt());
  G4double inelasticxsec=boonexsec->GetInelasticCrossSection(dynPart,aNucl.GetZ_asInt(),aNucl.GetA_asInt());
  delete dynPart;
  G4double qefrac=qexsec/inelasticxsec;
  
  if (G4UniformRand()<qefrac) {
    fIsQE=true;
    G4HadFinalState* part=fQuasiElasticModel.ApplyYourself(aTrack,aNucl);
    return part;
  }
  // if not QE, then generate the secondories from inelastic pBe
  // get number of secondaries, according to average multiplicties
  // and Poisson distributions, and
  // set momentum of secondaries, according to tabulated dataset
  // Note: no need to rotate secondaries' momentum in lab frame
  // using incident proton beam direction. This is now handled 
  // internally by G4

  G4long nProton = 0;
  G4long nNeutron = 0;
  G4long nPiPlus = 0;
  G4long nPiMinus = 0;
  G4long nKPlus = 0;
  G4long nKMinus = 0;
  G4long nKZeroLong = 0;

  // secondary Protons
  nProton = GetNumberOfProtons(aTrack);
  for(G4int iProton=0; iProton <  nProton; iProton++){
    secondaryMomentum =
      (GetMomentumOfSecondary(ProtonXSecArray, ProtonMaxXSec, aTrack));
    // Add this Proton to the list of secondaries
    aProton = new G4DynamicParticle;
    aProton->SetDefinition(G4Proton::ProtonDefinition());
    aProton->SetMomentum(secondaryMomentum);
    G4double invRwgtFactor
      = GetInverseRwgtFactor(trackMom, (aProton->GetMomentum()).z(),
			     (aProton->GetMomentum()).perp(), 14);
    G4HadSecondary *aProtonSec
      = new G4HadSecondary(aProton, invRwgtFactor);
    theParticleChange.AddSecondary(aProtonSec->GetParticle());
    //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
    //
    if (PhysicsVerbose) {
      G4cout << "Proton Number: " << iProton <<
	", (px,py,pz) in projectile frame = "
	     << secondaryMomentum << G4endl;    
    }
  }
  
  // secondary Neutrons
  nNeutron = GetNumberOfNeutrons(aTrack);
  for(G4int iNeutron=0; iNeutron <  nNeutron; iNeutron++){
    secondaryMomentum =
      (GetMomentumOfSecondary(NeutronXSecArray, NeutronMaxXSec, aTrack));
    // Add this Neutron to the list of secondaries
    aNeutron = new G4DynamicParticle;
    aNeutron->SetDefinition(G4Neutron::NeutronDefinition());
    aNeutron->SetMomentum(secondaryMomentum);
    G4double invRwgtFactor
      = GetInverseRwgtFactor(trackMom, (aNeutron->GetMomentum()).z(),
			     (aNeutron->GetMomentum()).perp(), 13);
    G4HadSecondary *aNeutronSec
      = new G4HadSecondary(aNeutron, invRwgtFactor);
    theParticleChange.AddSecondary(aNeutronSec->GetParticle());
    //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
    //
    if (PhysicsVerbose) {
      G4cout << "Neutron Number: " << iNeutron <<
	", (px,py,pz) in projectile frame = "
	     << secondaryMomentum << G4endl;    
    }
  }
  
  if (!NoBeamPions) {
    nPiPlus = GetNumberOfPiPluses(aTrack);
    for(G4int iPiPlus=0; iPiPlus <  nPiPlus; iPiPlus++){
      secondaryMomentum =
	(GetMomentumOfSecondary(PiPlusXSecArray, PiPlusMaxXSec, aTrack));
      // Add this PiPlus to the list of secondaries
      aPiPlus = new G4DynamicParticle;
      aPiPlus->SetDefinition(G4PionPlus::PionPlusDefinition());
      aPiPlus->SetMomentum(secondaryMomentum);
      G4double invRwgtFactor
	= GetInverseRwgtFactor(trackMom, (aPiPlus->GetMomentum()).z(),
			       (aPiPlus->GetMomentum()).perp(), 8);
      G4HadSecondary *aPiPlusSec
	  = new G4HadSecondary(aPiPlus, invRwgtFactor);
      theParticleChange.AddSecondary(aPiPlusSec->GetParticle());
      //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
      //
      if (PhysicsVerbose) {
	G4cout << "PiPlus Number: " << iPiPlus <<
	  ", (px,py,pz) in projectile frame = "
	       << secondaryMomentum << G4endl;    
      }
    }
    
    nPiMinus = GetNumberOfPiMinuses(aTrack);
    for(G4int iPiMinus=0; iPiMinus <  nPiMinus; iPiMinus++){
      secondaryMomentum =
	(GetMomentumOfSecondary(PiMinusXSecArray, PiMinusMaxXSec, aTrack));
      // Add this PiMinus to the list of secondaries
      aPiMinus = new G4DynamicParticle;
      aPiMinus->SetDefinition(G4PionMinus::PionMinusDefinition());
      aPiMinus->SetMomentum(secondaryMomentum);
      G4double invRwgtFactor
	= GetInverseRwgtFactor(trackMom, (aPiMinus->GetMomentum()).z(),
			       (aPiMinus->GetMomentum()).perp(), 9);
      G4HadSecondary *aPiMinusSec
	= new G4HadSecondary(aPiMinus, invRwgtFactor);
      theParticleChange.AddSecondary(aPiMinusSec->GetParticle());
      //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
      //
      if (PhysicsVerbose) {
	G4cout << "PiMinus Number: " << iPiMinus <<
	  ", (px,py,pz) in projectile frame = " 
	       << secondaryMomentum << G4endl;    
      }
    }
  }

  // secondary KPluses
  nKPlus = GetNumberOfKPluses(aTrack);
  for(G4int iKPlus=0; iKPlus <  nKPlus; iKPlus++){
    secondaryMomentum =
      (GetMomentumOfSecondary(KPlusXSecArray, KPlusMaxXSec, aTrack));
    // Add this KPlus to the list of secondaries
    aKPlus = new G4DynamicParticle;
    aKPlus->SetDefinition(G4KaonPlus::KaonPlusDefinition());
    aKPlus->SetMomentum(secondaryMomentum);
    G4double invRwgtFactor
      = GetInverseRwgtFactor(trackMom, (aKPlus->GetMomentum()).z(),
			     (aKPlus->GetMomentum()).perp(), 11);
    G4HadSecondary *aKPlusSec
      = new G4HadSecondary(aKPlus, invRwgtFactor);
    theParticleChange.AddSecondary(aKPlusSec->GetParticle());
    //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
    //
    if (PhysicsVerbose) {
      G4cout << "KPlus Number: " << iKPlus <<
	", (px,py,pz) in projectile frame = "
	     << secondaryMomentum << G4endl;    
    }
  }

  // secondary KMinuses
  nKMinus = GetNumberOfKMinuses(aTrack);
  for(G4int iKMinus=0; iKMinus <  nKMinus; iKMinus++){
    secondaryMomentum =
      (GetMomentumOfSecondary(KMinusXSecArray, KMinusMaxXSec, aTrack));
    // Add this KMinus to the list of secondaries
    aKMinus = new G4DynamicParticle;
    aKMinus->SetDefinition(G4KaonMinus::KaonMinusDefinition());
    aKMinus->SetMomentum(secondaryMomentum);
    G4double invRwgtFactor
      = GetInverseRwgtFactor(trackMom, (aKMinus->GetMomentum()).z(),
			     (aKMinus->GetMomentum()).perp(), 12);
    G4HadSecondary *aKMinusSec
      = new G4HadSecondary(aKMinus, invRwgtFactor);
    theParticleChange.AddSecondary(aKMinusSec->GetParticle());
    //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
    //
    if (PhysicsVerbose) {
      G4cout << "KMinus Number: " << iKMinus <<
	", (px,py,pz) in projectile frame = "
	     << secondaryMomentum << G4endl;    
    }
  }

  // secondary KZeroLongs
  nKZeroLong = GetNumberOfKZeroLongs(aTrack);
  for(G4int iKZeroLong=0; iKZeroLong <  nKZeroLong; iKZeroLong++){
    secondaryMomentum =
      (GetMomentumOfSecondary(KZeroLongXSecArray, KZeroLongMaxXSec, aTrack));
    // Add this KZeroLong to the list of secondaries
    aKZeroLong = new G4DynamicParticle;
    aKZeroLong->SetDefinition(G4KaonZeroLong::KaonZeroLongDefinition());
    aKZeroLong->SetMomentum(secondaryMomentum);
    G4double invRwgtFactor
      = GetInverseRwgtFactor(trackMom, (aKZeroLong->GetMomentum()).z(),
			     (aKZeroLong->GetMomentum()).perp(), 10);
    G4HadSecondary *aKZeroLongSec
      = new G4HadSecondary(aKZeroLong, invRwgtFactor);
    theParticleChange.AddSecondary(aKZeroLongSec->GetParticle());
    //theParticleChange.GetSecondary(theParticleChange.GetNumberOfSecondaries()-1)->SetWeight(invRwgtFactor);
    //
    if (PhysicsVerbose) {
      G4cout << "KZeroLong Number: " << iKZeroLong <<
	", (px,py,pz) in projectile frame = "
	     << secondaryMomentum << G4endl;    
    }
  }
  
  theParticleChange.SetStatusChange( stopAndKill );
  
  const NuBeamRunManager *pRunManager=
    reinterpret_cast<const NuBeamRunManager*>(G4RunManager::GetRunManager());
  NuBeamOutput* records = pRunManager->GetRecordPtr();
  if (records != NULL) records->RecordpBeInteraction(&theParticleChange);

  return &theParticleChange;
}

G4double BooNEpBeInteraction::GetMax(G4double Array[kNProtonMomentumBins][kNPzBins][kNPtBins]){
  G4double current, maxVal = 0.;
  for(G4int iipz=0; iipz<kNPzBins; iipz++) {
    for (G4int iipt=0; iipt<kNPtBins; iipt++) {
      for (G4int iiprotonp=0; iiprotonp<kNProtonMomentumBins; iiprotonp++) {
	current = Array[iiprotonp][iipz][iipt];
	if (current > maxVal) {maxVal = current;}
      }
    }
  }
  return maxVal;
}
    
G4double BooNEpBeInteraction::GetTotalProductionXSect(G4double Array[kNProtonMomentumBins][kNPzBins][kNPtBins], G4int iprotonp){
  G4double sum = 0.;
  for(G4int iipz=0; iipz < kNPzBins; iipz++){
    for(G4int iipt=0; iipt < kNPtBins; iipt++){
      sum = sum + Array[iprotonp][iipz][iipt]*(CLHEP::millibarn/(CLHEP::GeV*CLHEP::GeV));
    }
  }
  return sum;
}

// -----------------------------------------------------------------------------------

G4double BooNEpBeInteraction::GetInterpolatedXSec(G4double XSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins],
						  G4double protonMomentum,
						  G4double daughterPz,
						  G4double daughterPt)
{
  G4double ZerothOrderValue;
  G4double FirstOrderCorrection_x = 0;
  G4double FirstOrderCorrection_y = 0;
  G4double FirstOrderCorrection_z = 0;
  G4double FirstOrderCorrection;
  G4int jpz1, jpz2, jpt1, jpt2, midBin;

  // Verify that the input values are valid
  if (daughterPz < 0.) {
    G4cout << "ERROR in BooNEpBeInteraction::GetInterpolatedXSec: "
	   << "daughterPz (=" << daughterPz << ") < 0." << G4endl;
    G4cout << "     setting daughterPz to 0." << G4endl;
    daughterPz = 0.;
  }

  // Find protonp bin (jprotonp1)
  G4int jprotonp1 = 0;
  G4int jprotonp2 = kNProtonMomentumBins - 1;
  do {
    midBin = (jprotonp1 + jprotonp2)/2;
    if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
      jprotonp2 = midBin;
    else
      jprotonp1 = midBin;
  } while (jprotonp2 - jprotonp1 > 1);

  // Find pz bin (jpz1)
  jpz1 = 0;
  jpz2 = kNPzBins - 1;
  do {
    midBin = (jpz1 + jpz2)/2;
    if (daughterPz/CLHEP::GeV < fPzVec[midBin] )
      jpz2 = midBin;
    else
      jpz1 = midBin;
  } while (jpz2 - jpz1 > 1);

  // Find pt bin (jpt1)
  jpt1 = 0;
  jpt2 = kNPtBins - 1;
  do {
    midBin = (jpt1 + jpt2)/2;
    if (daughterPt/CLHEP::GeV < fPtVec[midBin])
      jpt2 = midBin;
    else
      jpt1 = midBin;
  } while (jpt2 - jpt1 > 1);

  // interpolated value (XSecValue)
  // 0th order
  ZerothOrderValue = XSecArray[jprotonp1][jpz1][jpt1]
    *(CLHEP::millibarn/(CLHEP::GeV*CLHEP::GeV));
  // 1st order
  FirstOrderCorrection_x = ((XSecArray[jprotonp1][jpz2][jpt1]
			     - XSecArray[jprotonp1][jpz1][jpt1])
			    * (daughterPz/CLHEP::GeV-fPzVec[jpz1])
			    / (fPzVec[jpz2]-fPzVec[jpz1]))
    *(CLHEP::millibarn/(CLHEP::GeV*CLHEP::GeV));
  FirstOrderCorrection_y = ((XSecArray[jprotonp1][jpz1][jpt2]
			     - XSecArray[jprotonp1][jpz1][jpt1])
			    * (daughterPt/CLHEP::GeV-fPtVec[jpt1])
			    / (fPtVec[jpt2]-fPtVec[jpt1]))
    *(CLHEP::millibarn/(CLHEP::GeV*CLHEP::GeV));
  if (kNProtonMomentumBins > 1) {
    FirstOrderCorrection_z = ((XSecArray[jprotonp2][jpz1][jpt1]
			       - XSecArray[jprotonp1][jpz1][jpt1])
			      * (protonMomentum/CLHEP::GeV-fProtonMomentumBins[jprotonp1])
			      / (fProtonMomentumBins[jprotonp2]-fProtonMomentumBins[jprotonp1]))
      *(CLHEP::millibarn/(CLHEP::GeV*CLHEP::GeV));
  } else {
    FirstOrderCorrection_z = 0.;
  }

  FirstOrderCorrection = (FirstOrderCorrection_x
			  + FirstOrderCorrection_y
			  + FirstOrderCorrection_z);

  return ZerothOrderValue + FirstOrderCorrection;
}

G4double BooNEpBeInteraction::GetInverseRwgtFactor(G4double protonMomentum,
						   G4double daughterPz,
						   G4double daughterPt,
						   G4int G3PartID)
{
  // Determine the weighted and unweighted interpolated cross-section value
  G4double noWgtXSec;
  G4double rwgtXSec;
  G4double rwgtFactor;

  if (G3PartID == 8) {
    noWgtXSec = GetInterpolatedXSec(PiPlusXSecNoWgtArray, protonMomentum,
				    daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(PiPlusXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else if (G3PartID == 9) {
    noWgtXSec = GetInterpolatedXSec(PiMinusXSecNoWgtArray, protonMomentum,
				    daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(PiMinusXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else if (G3PartID == 10) {
    noWgtXSec = GetInterpolatedXSec(KZeroLongXSecNoWgtArray, protonMomentum,
				    daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(KZeroLongXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else if (G3PartID == 11) {
    noWgtXSec = GetInterpolatedXSec(KPlusXSecNoWgtArray, protonMomentum,
				    daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(KPlusXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else if (G3PartID == 12) {
    noWgtXSec =  GetInterpolatedXSec(KMinusXSecNoWgtArray, protonMomentum,
				     daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(KMinusXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else if (G3PartID == 13) {
    noWgtXSec =  GetInterpolatedXSec(NeutronXSecNoWgtArray, protonMomentum,
				     daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(NeutronXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else if (G3PartID == 14) {
    noWgtXSec =  GetInterpolatedXSec(ProtonXSecNoWgtArray, protonMomentum,
				     daughterPz, daughterPt);
    rwgtXSec = GetInterpolatedXSec(ProtonXSecArray, protonMomentum,
				   daughterPz, daughterPt);
  } else {
    noWgtXSec = 1.;
    rwgtXSec = 1.;
  }

  // calculate and return the reweighting factor
  if (rwgtXSec > 0.) {
    rwgtFactor = noWgtXSec / rwgtXSec;
  } else {
    G4cout << "ERROR:  rwgtXSec (=" << rwgtXSec
	   << ") <= 0, noWgtXSec =" << noWgtXSec << G4endl;
    G4cout << "        protonMomentum, pz, pt = " << protonMomentum << ", "
	   << daughterPz << ", " << daughterPt << G4endl;
    rwgtFactor = 1.;
  }

  return rwgtFactor;

}

G4ThreeVector 
BooNEpBeInteraction::GetMomentumOfSecondary(G4double XSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins], G4double MaxXSec, const G4HadProjectile &aTrack){

  G4double trialPt, phi;
  G4double trialXSec;

  G4double protonMomentum = aTrack.GetTotalMomentum();

  do {
    // generate trial set of momentum components for the secondary. this
    // direction is still with respect to the coordinate system where the
    // primary particle direction is aligned with z axis
    trialPz = (G4UniformRand() * fPzVec[kNPzBins-1])*CLHEP::GeV;
    trialPt = G4UniformRand()*CLHEP::GeV;   
    phi = G4UniformRand() * CLHEP::twopi;
    trialPx = (trialPt * cos(phi));
    trialPy = (trialPt * sin(phi));
    trialP = G4ThreeVector(trialPx,trialPy,trialPz);

    XSecValue = GetInterpolatedXSec(XSecArray, protonMomentum, trialPz, trialPt);

    // Check whether to accept this set of (px,py,pz) values for this secondary,
    // based on the tabulated production cross-section 

    trialXSec = G4UniformRand()*MaxXSec*(CLHEP::millibarn/(CLHEP::GeV*CLHEP::GeV));

    if ( trialXSec < XSecValue )
      acceptSecondaryMomentum = true;
    else
      acceptSecondaryMomentum = false;
  } while ( !acceptSecondaryMomentum );
  // rotate momentum direction in global reference system is no longer necessary! (MS)

  // TRICK FOR COMPARING TO K2K !!!!!!!
  //if( atan2( sqrt(trialP.x()*trialP.x() + trialP.y()*trialP.y()), trialP.z()) > 0.4 )
  //{
  //  G4ThreeVector tmp = G4ThreeVector(0,0,0);
  //  return tmp;
  //}
  return trialP;
}

// -------------------------------------------------------------------------------

G4long BooNEpBeInteraction::GetNumberOfProtons(const G4HadProjectile &aTrack){

  static G4double ProtonAvMult[kNProtonMomentumBins];
  static G4bool FirstTimeProton = true;

  if (FirstTimeProton) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // Proton average multiplicity per inelastic collision
    // If sigma_proton is the inclusive proton production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average proton multiplicity per inelastic collision is
    // sigma_proton/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      ProtonAvMult[iprotonp] = (GetTotalProductionXSect(ProtonXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Average Proton Multiplicity at proton momentum "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << ProtonAvMult[iprotonp] << G4endl;
      }

    }

    FirstTimeProton = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = ProtonAvMult[jprotonp1]
      + (ProtonAvMult[jprotonp2] - ProtonAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = ProtonAvMult[0];
  }

  // The number of protons produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by avMult
  G4long nProton = G4Poisson( avMult );
  if (PhysicsVerbose) {
    G4cout << "Proton Multiplicity in this interaction = " << nProton
	   << G4endl;
  }
  return nProton;
}


// -------------------------------------------------------------------------------

G4long BooNEpBeInteraction::GetNumberOfNeutrons(const G4HadProjectile &aTrack){

  static G4double NeutronAvMult[kNProtonMomentumBins];
  static G4bool FirstTimeNeutron = true;

  if (FirstTimeNeutron) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // Neutron average multiplicity per inelastic collision
    // If sigma_neutron is the inclusive neutron production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average proton multiplicity per inelastic collision is
    // sigma_neutron/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      NeutronAvMult[iprotonp] = (GetTotalProductionXSect(NeutronXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Average Neutron Multiplicity at proton momentum "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << NeutronAvMult[iprotonp] << G4endl;
      }

    }

    FirstTimeNeutron = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = NeutronAvMult[jprotonp1]
      + (NeutronAvMult[jprotonp2] - NeutronAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = NeutronAvMult[0];
  }

  // The number of neutrons produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by NeutronAvMult
  G4long nNeutron = G4Poisson( avMult );
  if (PhysicsVerbose) {
    G4cout << "Neutron Multiplicity in this interaction = " << nNeutron
	   << G4endl;
  }
  return nNeutron;
}

// -------------------------------------------------------------------------


G4long BooNEpBeInteraction::GetNumberOfPiPluses(const G4HadProjectile &aTrack){

  static G4double PiPlusAvMult[kNProtonMomentumBins];
  static G4bool FirstTimePiPlus = true;

  if (FirstTimePiPlus) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // PiPlus average multiplicity per inelastic collision
    // If sigma_piplus is the inclusive pion production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average piplus multiplicity per inelastic collision is
    // sigma_piplus/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      PiPlusAvMult[iprotonp] = (GetTotalProductionXSect(PiPlusXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Avg PiPlus Multiplicity at p(p) "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << PiPlusAvMult[iprotonp] 
	       << "\t" <<iprotonp<<"\t"<<GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV)<<", "<<GetTotalProductionXSect(PiPlusXSecArray,iprotonp)<<
	  G4endl;
      }

    }
    FirstTimePiPlus = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = PiPlusAvMult[jprotonp1]
      + (PiPlusAvMult[jprotonp2] - PiPlusAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = PiPlusAvMult[0];
  }

  // The number of pipluses produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by PiPlusAvMult
  G4long nPiPlus = G4Poisson( avMult );

  if (PhysicsVerbose) {
    G4cout <<
      "PiPlus Multiplicity in this interaction = " << nPiPlus << G4endl;
  }
  return nPiPlus;
}

// -------------------------------------------------------------------------

G4long BooNEpBeInteraction::GetNumberOfPiMinuses(const G4HadProjectile &aTrack){

  static G4double PiMinusAvMult[kNProtonMomentumBins];
  static G4bool FirstTimePiMinus = true;

  if (FirstTimePiMinus) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // PiMinus average multiplicity per inelastic collision
    // If sigma_piminus is the inclusive pion production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average piminus multiplicity per inelastic collision is
    // sigma_piminus/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      PiMinusAvMult[iprotonp] = (GetTotalProductionXSect(PiMinusXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Average PiMinus Multiplicity at proton momentum "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << PiMinusAvMult[iprotonp] << G4endl;
      }

    }

    FirstTimePiMinus = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = PiMinusAvMult[jprotonp1]
      + (PiMinusAvMult[jprotonp2] - PiMinusAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = PiMinusAvMult[0];
  }

  // The number of piminuses produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by PiMinusAvMult
  G4long nPiMinus = G4Poisson( avMult );
  if (PhysicsVerbose) {
    G4cout <<
      "PiMinus Multiplicity in this interaction = " << nPiMinus << G4endl;
  }
  return nPiMinus;
}

// -------------------------------------------------------------------------

G4long BooNEpBeInteraction::GetNumberOfKPluses(const G4HadProjectile &aTrack){

  static G4double KPlusAvMult[kNProtonMomentumBins];
  static G4bool FirstTimeKPlus = true;

  if (FirstTimeKPlus) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // KPlus average multiplicity per inelastic collision
    // If sigma_kplus is the inclusive k production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average kplus multiplicity per inelastic collision is
    // sigma_kplus/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      KPlusAvMult[iprotonp] = (GetTotalProductionXSect(KPlusXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Average KPlus Multiplicity at proton momentum "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << KPlusAvMult[iprotonp] << G4endl;
      }

    }

    FirstTimeKPlus = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = KPlusAvMult[jprotonp1]
      + (KPlusAvMult[jprotonp2] - KPlusAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = KPlusAvMult[0];
  }

  // The number of kpluses produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by KPlusAvMult
  G4long nKPlus = G4Poisson( avMult );
  if (PhysicsVerbose) {
    G4cout <<
      "KPlus Multiplicity in this interaction = " << nKPlus << G4endl;
  }
  return nKPlus;
}

// -------------------------------------------------------------------------

G4long BooNEpBeInteraction::GetNumberOfKMinuses(const G4HadProjectile &aTrack){

  static G4double KMinusAvMult[kNProtonMomentumBins];
  static G4bool FirstTimeKMinus = true;

  if (FirstTimeKMinus) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // KMinus average multiplicity per inelastic collision
    // If sigma_kminus is the inclusive k production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average kminus multiplicity per inelastic collision is
    // sigma_kminus/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      KMinusAvMult[iprotonp] = (GetTotalProductionXSect(KMinusXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Average KMinus Multiplicity at proton momentum "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << KMinusAvMult[iprotonp] << G4endl;
      }

    }

    FirstTimeKMinus = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = KMinusAvMult[jprotonp1]
      + (KMinusAvMult[jprotonp2] - KMinusAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = KMinusAvMult[0];
  }

  // The number of kminuses produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by KMinusAvMult
  G4long nKMinus = G4Poisson( avMult );
  if (PhysicsVerbose) {
    G4cout <<
      "KMinus Multiplicity in this interaction = " << nKMinus << G4endl;
  }
  return nKMinus;
}


// -------------------------------------------------------------------------

G4long BooNEpBeInteraction::GetNumberOfKZeroLongs(const G4HadProjectile &aTrack){

  static G4double KZeroLongAvMult[kNProtonMomentumBins];
  static G4bool FirstTimeKZeroLong = true;

  if (FirstTimeKZeroLong) {

    // longitudinal and transverse momentum bin widths, from tabulated data set
    G4double Deltapz = (fPzVec[1]-fPzVec[0])*CLHEP::GeV;
    G4double Deltapt = (fPtVec[1]-fPtVec[0])*CLHEP::GeV;

    // KZeroLong average multiplicity per inelastic collision
    // If sigma_kminus is the inclusive k production xsect in p-Be
    // interactions,
    // and sigma_inel is the total p-Be inelastic cross-section,
    // the average kminus multiplicity per inelastic collision is
    // sigma_kminus/sigma_inel
    for (G4int iprotonp=0; iprotonp<kNProtonMomentumBins; iprotonp++) {
      KZeroLongAvMult[iprotonp] = (GetTotalProductionXSect(KZeroLongXSecArray,iprotonp))*
	Deltapz*Deltapt / GetpBeInelasticCrossSection(fProtonMomentumBins[iprotonp]*CLHEP::GeV);

      if (PhysicsVerbose) {
	G4cout << "Average KZeroLong Multiplicity at proton momentum "
	       << fProtonMomentumBins[iprotonp] << " GeV = "
	       << KZeroLongAvMult[iprotonp] << G4endl;
      }

    }

    FirstTimeKZeroLong = false;
  }

  // Find the proton momentum bin
  G4double avMult;
  if (kNProtonMomentumBins > 1) {
    G4double protonMomentum = aTrack.GetTotalMomentum();
    G4int jprotonp1 = 0;
    G4int jprotonp2 = kNProtonMomentumBins - 1;
    G4int midBin;
    do {
      midBin = (jprotonp1 + jprotonp2)/2;
      if (protonMomentum/CLHEP::GeV < fProtonMomentumBins[midBin] )
	jprotonp2 = midBin;
      else
	jprotonp1 = midBin;
    } while (jprotonp2 - jprotonp1 > 1);

    avMult = KZeroLongAvMult[jprotonp1]
      + (KZeroLongAvMult[jprotonp2] - KZeroLongAvMult[jprotonp1])
      * (protonMomentum/CLHEP::GeV - fProtonMomentumBins[jprotonp1])
      / (fProtonMomentumBins[jprotonp2] - fProtonMomentumBins[jprotonp1]);
  } else {
    avMult = KZeroLongAvMult[0];
  }

  // The number of kminuses produced in this particular inelastic collision
  // is determined by a Poisson distr., with mean given by KZeroLongAvMult
  G4long nKZeroLong = G4Poisson( avMult );
  if (PhysicsVerbose) {
    G4cout <<
      "KZeroLong Multiplicity in this interaction = " << nKZeroLong << G4endl;
  }
  return nKZeroLong;
}

// -------------------------------------------------------------------------

G4double BooNEpBeInteraction::GetpBeInelasticCrossSection( G4double totalMomentum)
{
  G4double xsec=0;
  //First check if BooNExsec is applicable,
  BooNEHadronCrossSections* boonexsec= BooNEHadronCrossSections::Instance();
  // Create a proton with momentum = totalMomentum.
  // The direction doesn't matter since we are only interested in the
  // value of the p-Be cross section
  G4ThreeVector momentumVector;
  momentumVector.set(0.,0.,totalMomentum);
  G4DynamicParticle* dummyProton = new G4DynamicParticle(G4Proton::ProtonDefinition(), momentumVector); 
  if (boonexsec->IsApplicable(dummyProton, 4, 9)) {
    //old g4 mc was subtracting quasielastic cross section from inelastic, but we don't do that any more since quasi elastic process is included in inelastic. to get the same multiplicities as before we need to subtract the quasi elastic here
    G4double qexsec=boonexsec->GetQuasiElasticCrossSection(dummyProton,4,9);
    xsec=boonexsec->GetInelasticCrossSection(dummyProton,4,9)-qexsec;
  } else {
    //find hadronic process and get cross section
    G4ProcessManager* pManager=G4Proton::Proton()->GetProcessManager();
    G4ProcessVector* procVec=pManager->GetProcessList();

    for (int iproc=0; iproc<pManager->GetProcessListLength();iproc++) {
      // Get the total cross section at this momentum from hadron inelastic process
      if ((*procVec)[iproc]->GetProcessType()==G4ProcessType::fHadronic &&
	  (*procVec)[iproc]->GetProcessSubType()==G4HadronicProcessType::fHadronInelastic) {
	G4HadronicProcess* hadProc=dynamic_cast<G4HadronicProcess*> ((*procVec)[iproc]);

	G4Element * myElement = G4Element::GetElement("Be");
	G4Material * myMat = G4Material::GetMaterial("Be");
	if( ! myMat ) {
	  myMat = new G4Material("Be", 1.85*CLHEP::g/CLHEP::cm3, 1);
	  myMat->AddElement(myElement, 1.0);
	}

	//xsec = hadProc->GetElementCrossSection(dummyProton, G4Element::GetElement("Be"));
	xsec = hadProc->GetElementCrossSection(dummyProton, myElement, myMat);
      }
    }
  }
  delete dummyProton;

  return xsec;

}

// -------------------------------------------------------------------------

G4Material* BooNEpBeInteraction::GetTargetMaterial( ){

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();

  G4Material* targetMat;
  targetMat = NULL;

  for( size_t i = 0; i < pLVStore->size(); i++ ) 
  {
    G4LogicalVolume* lv = (*pLVStore)[i];

    if ( lv->GetName() == "TARG" )
    {
      targetMat = lv->GetMaterial();
    }
  }
  return targetMat;;

}

G4double BooNEpBeInteraction::GetFeynmanX( double p, double th, double m, double pbeam, double mbeam ){

  double energy = sqrt( mbeam*mbeam + pbeam*pbeam );
  double Ecm = sqrt(2.*mbeam*mbeam + 2.*mbeam*energy);
  double pcmmax = Ecm/2.;
  double gamma = pbeam/Ecm;
  double beta = sqrt(1. - 1./(gamma*gamma));
  double e = sqrt(p*p + m*m);
  double pcmz = gamma*(p*cos(th) - beta*e);
  double x = pcmz/pcmmax;
  
  return x;

}

void BooNEpBeInteraction::FillOtherCrossSections(G4double kineticEnergy)
{
  G4int bin = 0;
  
  kineticEnergy /=CLHEP::GeV;

  //  MARS
  bin = (G4int)((kineticEnergy-fProtonKE[0])/((fProtonKE[kNEBins-1]-fProtonKE[0])/double(kNEBins)));
  G4cout<<"Filling cross section from bin "<<bin<<"\t"<<fProtonKE[kNEBins-1]<<"\t"<<fProtonKE[0]<<"\t"<<kNEBins<<G4endl;
  if( bin >= kNEBins-1 ) {
    bin = kNEBins - 1;
    xSectpBeInel_MARS = totalxSectpBeInel_MARS[kNEBins-1]*CLHEP::millibarn;
  } else {
    G4double m = (totalxSectpBeInel_MARS[bin] - totalxSectpBeInel_MARS[bin+1]) / (fProtonKE[bin] - fProtonKE[bin+1]);
    G4double b =  totalxSectpBeInel_MARS[bin] - m*fProtonKE[bin];
    xSectpBeInel_MARS = (m*kineticEnergy + b)*CLHEP::millibarn;
  }

  //  GFLUKA
  bin = (G4int)(( kineticEnergy-fProtonKE[0] ) / ((fProtonKE[kNEBins-1] - fProtonKE[0])/double(kNEBins)) );
  if( bin >= kNEBins-1 ) {
    bin = kNEBins - 1;
    xSectpBeInel_GFLUKA = totalxSectpBeInel_GFLUKA[kNEBins-1]*CLHEP::millibarn;
  } else {  
    G4double m = (totalxSectpBeInel_GFLUKA[bin] - totalxSectpBeInel_GFLUKA[bin+1]) / (fProtonKE[bin] - fProtonKE[bin+1]);
    G4double b =  totalxSectpBeInel_GFLUKA[bin] - m*fProtonKE[bin];
    xSectpBeInel_GFLUKA = (m*kineticEnergy + b)*CLHEP::millibarn;
  }

}

G4double BooNEpBeInteraction::GetRwgtXSec(G4double NoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
					  G4int iprotonp, G4int ipz,
					  G4int ipt, G4int G3PartID)
{
  G4double rwgtXSec;

  if (funcID[G3PartID] == 0) {

    rwgtXSec = NoWgtXSec[iprotonp][ipz][ipt];

  }  else if (funcID[G3PartID] == 1) {

    rwgtXSec = NoWgtXSec[iprotonp][ipz][ipt]
      *(rwgtParams[G3PartID][0]
	+ rwgtParams[G3PartID][1] * fPzVec[ipz]
	+ rwgtParams[G3PartID][2] * pow(fPzVec[ipz],2)
	+ rwgtParams[G3PartID][3] * pow(fPzVec[ipz],3)
	+ rwgtParams[G3PartID][4] * pow(fPzVec[ipz],4)
	+ rwgtParams[G3PartID][5] * pow(fPzVec[ipz],5)
	+ rwgtParams[G3PartID][6] * pow(fPzVec[ipz],6)
	+ rwgtParams[G3PartID][7] * pow(fPzVec[ipz],7)
	+ rwgtParams[G3PartID][8] * pow(fPzVec[ipz],8)
	+ rwgtParams[G3PartID][9] * pow(fPzVec[ipz],9));

  }  else if (funcID[G3PartID] == 2) {

    rwgtXSec =  NoWgtXSec[iprotonp][ipz][ipt]
      *(rwgtParams[G3PartID][0]
	+ rwgtParams[G3PartID][1]
	* exp(rwgtParams[G3PartID][2] * (fPzVec[ipz]-rwgtParams[G3PartID][3])));

  }  else if (funcID[G3PartID] == 3) {

    G4double p = sqrt(fPzVec[ipz]*fPzVec[ipz] + fPtVec[ipt]*fPtVec[ipt]);
    if (p < fProtonMomentumBins[iprotonp]) {
      rwgtXSec =  rwgtParams[G3PartID][0] * 2.*M_PI*fPtVec[ipt]/(p*p);
    } else {
      rwgtXSec =  0.;
    }

  }  else if (funcID[G3PartID] == 4) {

    G4double p = sqrt(fPzVec[ipz]*fPzVec[ipz] + fPtVec[ipt]*fPtVec[ipt]);
    if (p < fProtonMomentumBins[iprotonp]) {
      rwgtXSec =  rwgtParams[G3PartID][0];
    } else {
      rwgtXSec =  0.;
    }

  } else {

    G4cerr << "ERROR in BooNEpBeInteraction::GetRwgtXSec:  funcID["
	   << G3PartID << "] = " << funcID[G3PartID]
	   << ", which is outside the range [0,4]" << G4endl;
    rwgtXSec = 1.;

  }

  return rwgtXSec;

}

void BooNEpBeInteraction::SetProtonRwgtFunc(G4String newValue) 
{

  G4cout << "Proton reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[14] = 0;
  } else if (newValue == "POLY") {
    funcID[14] = 1;
  } else if (newValue == "EXP") {
    funcID[14] = 2;
  } else if (newValue == "FLAT1") {
    funcID[14] = 3;
  } else if (newValue == "FLAT2") {
    funcID[14] = 4;
  } else {
    funcID[14] = 0;
    G4cout << "ERROR:  unknown proton reweighting function option: "
	   << newValue << G4endl;
  }

}

void BooNEpBeInteraction::SetNeutronRwgtFunc(G4String newValue) 
{

  G4cout << "Neutron reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[13] = 0;
  } else if (newValue == "POLY") {
    funcID[13] = 1;
  } else if (newValue == "EXP") {
    funcID[13] = 2;
  } else if (newValue == "FLAT1") {
    funcID[13] = 3;
  } else if (newValue == "FLAT2") {
    funcID[13] = 4;
  } else {
    funcID[13] = 0;
    G4cout << "ERROR:  unknown neutron reweighting function option: "
	   << newValue << G4endl;
  }

}

void BooNEpBeInteraction::SetPionPlusRwgtFunc(G4String newValue) 
{

  G4cout << "PionPlus reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[8] = 0;
  } else if (newValue == "POLY") {
    funcID[8] = 1;
  } else if (newValue == "EXP") {
    funcID[8] = 2;
  } else if (newValue == "FLAT1") {
    funcID[8] = 3;
  } else if (newValue == "FLAT2") {
    funcID[8] = 4;
  } else {
    funcID[8] = 0;
    G4cout << "ERROR:  unknown pionPlus reweighting function option: "
	   << newValue << G4endl;
  }
}

void BooNEpBeInteraction::SetPionMinusRwgtFunc(G4String newValue) 
{

  G4cout << "PionMinus reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[9] = 0;
  } else if (newValue == "POLY") {
    funcID[9] = 1;
  } else if (newValue == "EXP") {
    funcID[9] = 2;
  } else if (newValue == "FLAT1") {
    funcID[9] = 3;
  } else if (newValue == "FLAT2") {
    funcID[9] = 4;
  } else {
    funcID[9] = 0;
    G4cout << "ERROR:  unknown pionMinus reweighting function option: "
	   << newValue << G4endl;
  }

}

void BooNEpBeInteraction::SetKaonPlusRwgtFunc(G4String newValue) 
{

  G4cout << "KaonPlus reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[11] = 0;
  } else if (newValue == "POLY") {
    funcID[11] = 1;
  } else if (newValue == "EXP") {
    funcID[11] = 2;
  } else if (newValue == "FLAT1") {
    funcID[11] = 3;
  } else if (newValue == "FLAT2") {
    funcID[11] = 4;
  } else {
    funcID[11] = 0;
    G4cout << "ERROR:  unknown kaonPlus reweighting function option: "
	   << newValue << G4endl;
  }

}

void BooNEpBeInteraction::SetKaonMinusRwgtFunc(G4String newValue) 
{

  G4cout << "KaonMinus reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[12] = 0;
  } else if (newValue == "POLY") {
    funcID[12] = 1;
  } else if (newValue == "EXP") {
    funcID[12] = 2;
  } else if (newValue == "FLAT1") {
    funcID[12] = 3;
  } else if (newValue == "FLAT2") {
    funcID[12] = 4;
  } else {
    funcID[12] = 0;
    G4cout << "ERROR:  unknown kaonMinus reweighting function option: "
	   << newValue << G4endl;
  }

}

void BooNEpBeInteraction::SetKaonZeroLongRwgtFunc(G4String newValue) 
{

  G4cout << "KaonZeroLong reweighting function set to " << newValue << G4endl;
  if (newValue == "NONE") {
    funcID[10] = 0;
  } else if (newValue == "POLY") {
    funcID[10] = 1;
  } else if (newValue == "EXP") {
    funcID[10] = 2;
  } else if (newValue == "FLAT1") {
    funcID[10] = 3;
  } else if (newValue == "FLAT2") {
    funcID[10] = 4;
  } else {
    funcID[10] = 0;
    G4cout << "ERROR:  unknown kaonZeroLong reweighting function option: "
	   << newValue << G4endl;
  }

}

void BooNEpBeInteraction::SetProtonRwgtParams(G4double val0, G4double val1,
					      G4double val2, G4double val3,
					      G4double val4, G4double val5,
					      G4double val6, G4double val7,
					      G4double val8, G4double val9)
{

  rwgtParams[14][0] = val0;
  rwgtParams[14][1] = val1;
  rwgtParams[14][2] = val2;
  rwgtParams[14][3] = val3;
  rwgtParams[14][4] = val4;
  rwgtParams[14][5] = val5;
  rwgtParams[14][6] = val6;
  rwgtParams[14][7] = val7;
  rwgtParams[14][8] = val8;
  rwgtParams[14][9] = val9;

  G4cout << "Proton reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}

void BooNEpBeInteraction::SetNeutronRwgtParams(G4double val0, G4double val1,
					       G4double val2, G4double val3,
					       G4double val4, G4double val5,
					       G4double val6, G4double val7,
					       G4double val8, G4double val9)
{

  rwgtParams[13][0] = val0;
  rwgtParams[13][1] = val1;
  rwgtParams[13][2] = val2;
  rwgtParams[13][3] = val3;
  rwgtParams[13][4] = val4;
  rwgtParams[13][5] = val5;
  rwgtParams[13][6] = val6;
  rwgtParams[13][7] = val7;
  rwgtParams[13][8] = val8;
  rwgtParams[13][9] = val9;

  G4cout << "Neutron reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}

void BooNEpBeInteraction::SetPionPlusRwgtParams(G4double val0, G4double val1,
						G4double val2, G4double val3,
						G4double val4, G4double val5,
						G4double val6, G4double val7,
						G4double val8, G4double val9)
{

  rwgtParams[8][0] = val0;
  rwgtParams[8][1] = val1;
  rwgtParams[8][2] = val2;
  rwgtParams[8][3] = val3;
  rwgtParams[8][4] = val4;
  rwgtParams[8][5] = val5;
  rwgtParams[8][6] = val6;
  rwgtParams[8][7] = val7;
  rwgtParams[8][8] = val8;
  rwgtParams[8][9] = val9;

  G4cout << "PionPlus reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}

void BooNEpBeInteraction::SetPionMinusRwgtParams(G4double val0, G4double val1,
						 G4double val2, G4double val3,
						 G4double val4, G4double val5,
						 G4double val6, G4double val7,
						 G4double val8, G4double val9)
{

  rwgtParams[9][0] = val0;
  rwgtParams[9][1] = val1;
  rwgtParams[9][2] = val2;
  rwgtParams[9][3] = val3;
  rwgtParams[9][4] = val4;
  rwgtParams[9][5] = val5;
  rwgtParams[9][6] = val6;
  rwgtParams[9][7] = val7;
  rwgtParams[9][8] = val8;
  rwgtParams[9][9] = val9;

  G4cout << "PionMinus reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}

void BooNEpBeInteraction::SetKaonPlusRwgtParams(G4double val0, G4double val1,
						G4double val2, G4double val3,
						G4double val4, G4double val5,
						G4double val6, G4double val7,
						G4double val8, G4double val9)
{

  rwgtParams[11][0] = val0;
  rwgtParams[11][1] = val1;
  rwgtParams[11][2] = val2;
  rwgtParams[11][3] = val3;
  rwgtParams[11][4] = val4;
  rwgtParams[11][5] = val5;
  rwgtParams[11][6] = val6;
  rwgtParams[11][7] = val7;
  rwgtParams[11][8] = val8;
  rwgtParams[11][9] = val9;

  G4cout << "KaonPlus reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}

void BooNEpBeInteraction::SetKaonMinusRwgtParams(G4double val0, G4double val1,
						 G4double val2, G4double val3,
						 G4double val4, G4double val5,
						 G4double val6, G4double val7,
						 G4double val8, G4double val9)
{

  rwgtParams[12][0] = val0;
  rwgtParams[12][1] = val1;
  rwgtParams[12][2] = val2;
  rwgtParams[12][3] = val3;
  rwgtParams[12][4] = val4;
  rwgtParams[12][5] = val5;
  rwgtParams[12][6] = val6;
  rwgtParams[12][7] = val7;
  rwgtParams[12][8] = val8;
  rwgtParams[12][9] = val9;

  G4cout << "KaonMinus reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}

void BooNEpBeInteraction::SetKaonZeroLongRwgtParams(G4double val0, G4double val1,
						    G4double val2, G4double val3,
						    G4double val4, G4double val5,
						    G4double val6, G4double val7,
						    G4double val8, G4double val9)
{

  rwgtParams[10][0] = val0;
  rwgtParams[10][1] = val1;
  rwgtParams[10][2] = val2;
  rwgtParams[10][3] = val3;
  rwgtParams[10][4] = val4;
  rwgtParams[10][5] = val5;
  rwgtParams[10][6] = val6;
  rwgtParams[10][7] = val7;
  rwgtParams[10][8] = val8;
  rwgtParams[10][9] = val9;

  G4cout << "KaonZeroLong reweighting parameters set to: "
	 << val0 << " "
	 << val1 << " "
	 << val2 << " "
	 << val3 << " "
	 << val4 << " "
	 << val5 << " "
	 << val6 << " "
	 << val7 << " "
	 << val8 << " "
	 << val9 << " "
	 << G4endl;

}
