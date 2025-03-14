#include "BooNEpXInteraction.hh"
#include "BooNEpBeInteraction.hh"
#include "BooNEHadronPhysics.hh"
#include "XSecArrayHolder.hh"
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

// include math package
#include <math.h>

// -----------------------------------------------------------------------------------
BooNEpXInteraction::BooNEpXInteraction()
  :G4HadronicInteraction("BooNEpXInteraction"), fInteraction(0), fIsProton(false)
{
}
// -----------------------------------------------------------------------------------
BooNEpXInteraction::BooNEpXInteraction(G4HadronicInteraction * hadInt, G4bool proton)
  :G4HadronicInteraction("BooNEpXInteraction"), fInteraction(hadInt), fIsProton(proton)
{   
}
// -----------------------------------------------------------------------------------
BooNEpXInteraction::~BooNEpXInteraction()
{
}
// -----------------------------------------------------------------------------------
G4HadFinalState * BooNEpXInteraction::ApplyYourself( const G4HadProjectile & aTrack,
						     G4Nucleus & targetNucleus )
{
  G4HadFinalState * hadFS = 0;

  hadFS = fInteraction->ApplyYourself( aTrack, targetNucleus );

  /*
  G4cout << "\nFor interaction of track with particle name "
	 << aTrack.GetDefinition()->GetParticleName() << " and model name " 
	 << this->GetIntModelName()
	 << " and KE = "
	 << aTrack.GetKineticEnergy() << ":" << G4endl;

  G4cout << "\n";
  for( G4int i = 0; i < hadFS->GetNumberOfSecondaries(); i++ ) {
    G4HadSecondary * sec = hadFS->GetSecondary( static_cast<size_t>(i) );
    G4cout << "sec " << i << " name " << sec->GetParticle()->GetParticleDefinition()->GetParticleName()
	   << " weight " << sec->GetWeight()
	   << " KE " << sec->GetParticle()->GetKineticEnergy() << "\n"; 
  }
  G4cout << "\n" << G4endl;
  */

  if( !fIsProton ) return hadFS; // no further modifications if not a proton interaction

  // now we need to calculate the weight of each of these secondaries!
  for( G4int i = 0; i < hadFS->GetNumberOfSecondaries(); i++ ) {
    G4HadSecondary * sec = hadFS->GetSecondary( static_cast<size_t>(i) );

    G4double trackMom = aTrack.GetTotalMomentum();
    // RETHERE: If need to consider additional multiplicities, can do here.

    G4double invRwgtFactor = 1.0;
    G4DynamicParticle * particle = sec->GetParticle();
    G4double pz = particle->Get4Momentum().z();
    G4double pT = particle->Get4Momentum().perp();

    G4double part_id = 0;
    if( particle->GetParticleDefinition() == G4Proton::ProtonDefinition() ) part_id = 14;
    else if( particle->GetParticleDefinition() == G4Neutron::NeutronDefinition() ) part_id = 13;
    else if( particle->GetParticleDefinition() == G4PionPlus::PionPlusDefinition() ) part_id = 8;
    else if( particle->GetParticleDefinition() == G4PionMinus::PionMinusDefinition() ) part_id = 9;
    else if( particle->GetParticleDefinition() == G4KaonPlus::KaonPlusDefinition() ) part_id = 11;
    else if( particle->GetParticleDefinition() == G4KaonMinus::KaonMinusDefinition() ) part_id = 12;
    else if( particle->GetParticleDefinition() == G4KaonZeroLong::KaonZeroLongDefinition() ) part_id = 10;

    // part_id == 0 is things like gammas, alphas, deuterons, rogue elements such as Cr51/Cr52...
    if( aTrack.GetDefinition() == G4Proton::ProtonDefinition() && part_id != 0
	&& strcmp( this->GetIntModelName(), "BooNEpBeInteraction" ) == 0 )  {
      //G4cout << sec->GetParticle()->GetParticleDefinition()->GetParticleName() << " " <<
      //	this->GetIntModelName() << G4endl;
      invRwgtFactor = sec->GetWeight() * GetInverseRwgtFactor( trackMom, pz, pT, part_id ); // if BooNEpBe didn't already do it
    } else invRwgtFactor = sec->GetWeight();

    sec->SetWeight( invRwgtFactor );

  } // for secondary

  /*
  G4cout << "\nI did a proton! Showing updated secondaries with weights " << G4endl;
  G4cout << "\n";
  for( G4int i = 0; i < hadFS->GetNumberOfSecondaries(); i++ ) {
    G4HadSecondary * sec = hadFS->GetSecondary( static_cast<size_t>(i) );
    G4cout << "sec " << i << " name " << sec->GetParticle()->GetParticleDefinition()->GetParticleName()
	   << " weight " << sec->GetWeight()
	   << " KE " << sec->GetParticle()->GetKineticEnergy() << "\n"; 
  }
  G4cout << "\n" << G4endl;
  */
  
  return hadFS;
}
// -----------------------------------------------------------------------------------
G4double BooNEpXInteraction::GetInterpolatedXSec( G4double XSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins],
						  G4double protonMomentum, 
						  G4double daughterPz, G4double daughterPt,
						  G4bool debug)
{
  G4double ZerothOrderValue;
  G4double FirstOrderCorrection_x = 0;
  G4double FirstOrderCorrection_y = 0;
  G4double FirstOrderCorrection_z = 0;
  G4double FirstOrderCorrection;
  G4int jpz1, jpz2, jpt1, jpt2, midBin;

  // Verify that the input values are valid
  if (daughterPz < 0.) {
    G4cout << "ERROR in BooNEpXInteraction::GetInterpolatedXSec: "
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

  // Output the value of the thing
  if( debug )
    {
      G4cout << "Getting array values" << G4endl;
      G4cout << "ranges: proton p: [ " << fProtonMomentumBins[jprotonp1]
	     << ", " << fProtonMomentumBins[jprotonp2] << " ] , daughter pz : [ "
	     << fPzVec[jpz1] << ", " << fPzVec[jpz2] << " ] , daughter pt: [ "
	     << fPtVec[jpt1] << ", " << fPtVec[jpt2] << " ]"
	     << "\n\tarray[" << jprotonp1 << "][" << jpz1 << "][" << jpt1 << "] = "
	     << XSecArray[jprotonp1][jpz1][jpt1]
	     << "\n\tarray_x[" << jprotonp1 << "][" << jpz2 << "][" << jpt1 << "] = "
	     << XSecArray[jprotonp1][jpz2][jpt1]
	     << "\n\tarray_y[" << jprotonp1 << "][" << jpz1 << "][" << jpt2 << "] = "
	     << XSecArray[jprotonp1][jpz1][jpt2]
	     << "\n\tarray_z[" << jprotonp2 << "][" << jpz1 << "][" << jpt1 << "] = "
	     << XSecArray[jprotonp2][jpz1][jpt1] << G4endl;
    }

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
// -----------------------------------------------------------------------------------
G4double BooNEpXInteraction::GetInverseRwgtFactor( G4double protonMomentum, G4double daughterPz,
						   G4double daughterPt, G4int G3PartID )
{
  // Determine the weighted and unweighted interpolated cross-section value
  G4double noWgtXSec;
  G4double rwgtXSec;
  G4double rwgtFactor;

  if( daughterPz < 0.0 ) {
    //G4cout << "Negative daughter pz = " << daughterPz << " for particle with G3PartID = "
    //	   << G3PartID << ", returning 1..." << G4endl;
    return 1.0;
  }

  XSecArrayHolder& xsec_holder = XSecArrayHolder::Instance();
  switch( G3PartID ) {
  case 8: 
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetPiPlusNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetPiPlusXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  case 9:
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetPiMinusNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetPiMinusXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  case 10:
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetKZeroLongNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetKZeroLongXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  case 11:
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetKPlusNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetKPlusXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  case 12:
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetKMinusNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetKMinusXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  case 13:
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetNeutronNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetNeutronXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  case 14:
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetProtonNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetProtonXSec(),      protonMomentum, daughterPz, daughterPt);
    break;
  default:
    // return 1.0 reweight
    noWgtXSec = GetInterpolatedXSec(xsec_holder.GetPiPlusNoWgtXSec(), protonMomentum, daughterPz, daughterPt);
    rwgtXSec  = GetInterpolatedXSec(xsec_holder.GetPiPlusNoWgtXSec(), protonMomentum, daughterPz, daughterPt); 
    break;
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
