#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
//#include "G4KL3DecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

#include "BooNEKaonDecayChannel.hh"
#include <memory>

// -----------------------------------------------------------------------------------
BooNEKaonDecayChannel::BooNEKaonDecayChannel() : G4VDecayChannel(), pLambda( 0.0 ), pXi0( 0.0 ) {}
// -----------------------------------------------------------------------------------
BooNEKaonDecayChannel::BooNEKaonDecayChannel( const G4String & theParentName,
					      G4double         theBR,
					      const G4String & thePionName,
					      const G4String & theLeptonName,
					      const G4String & theNeutrinoName ) :
  G4VDecayChannel( "KL3 Decay", theParentName, theBR, 3,
		   thePionName, theLeptonName, theNeutrinoName )
{
  static const G4String K_plus("kaon+");
  static const G4String K_minus("kaon-");
  static const G4String K_L("kaon0L");
  static const G4String Mu_plus("mu+");
  static const G4String Mu_minus("mu-");
  static const G4String E_plus("e+");
  static const G4String E_minus("e-");

  // Here the constants are set. G4KL3DecayChannel uses up-to-date PDG values.
  // For now, I make the decision to ignore these and set the mBooNE values

  pLambda = 0.0282;
  pXi0    = -0.19;

  /*
  // check modes
  if ( ((theParentName == K_plus)&&(theLeptonName == E_plus)) ||
       ((theParentName == K_minus)&&(theLeptonName == E_minus))   ) {
    // K+- (Ke3)
    pLambda = 0.0286;
    pXi0    = -0.35;
   } else if ( ((theParentName == K_plus)&&(theLeptonName == Mu_plus)) ||
       ((theParentName == K_minus)&&(theLeptonName == Mu_minus))   ) {
    // K+- (Kmu3)
    pLambda = 0.033;
    pXi0    = -0.35;
  } else if ( (theParentName == K_L) && 
              ((theLeptonName == E_plus) ||(theLeptonName == E_minus))  ){
    // K0L (Ke3)
    pLambda = 0.0300;
    pXi0    = -0.11;
  } else if ( (theParentName == K_L) && 
              ((theLeptonName == Mu_plus) ||(theLeptonName == Mu_minus))  ){
    // K0L (Kmu3)
    pLambda = 0.034;
    pXi0    = -0.11;
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4KL3DecayChannel:: constructor :";
      G4cout << "illegal arguments " << G4endl;;
      DumpInfo();
    }
#endif
    // set values for K0L (Ke3) temporarily
    pLambda = 0.0300;
    pXi0    = -0.11;
  }
  */
}
BooNEKaonDecayChannel::BooNEKaonDecayChannel( const G4String & theParentName,
					      G4double         theBR,
					      const G4String & thePionName,
					      const G4String & theLeptonName,
					      const G4String & theNeutrinoName,
					      const G4double theLambda, const G4double theXi0 ) :
  G4VDecayChannel( "KL3 Decay", theParentName, theBR, 3,
		   thePionName, theLeptonName, theNeutrinoName ), 
  pLambda(theLambda), pXi0(theXi0)
{
  G4cout << "Constructing decay channel with lambda = " << pLambda 
	 << " and xi0 = " << pXi0 << G4endl;
}
// -----------------------------------------------------------------------------------
BooNEKaonDecayChannel::~BooNEKaonDecayChannel() {}
// -----------------------------------------------------------------------------------
BooNEKaonDecayChannel::BooNEKaonDecayChannel(const BooNEKaonDecayChannel & right) :
  G4VDecayChannel(right), pLambda(right.pLambda), pXi0(right.pXi0) {}
// -----------------------------------------------------------------------------------
BooNEKaonDecayChannel & BooNEKaonDecayChannel::operator = (const BooNEKaonDecayChannel & right)
{
  if( this != &right ) { // only update if diff pointer
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    parent_name = new G4String( *right.parent_name );
    // Clear array of daughter names and recreate it
    ClearDaughtersName();
    numberOfDaughters = right.numberOfDaughters;
    if( numberOfDaughters > 0 ) {
      if( daughters_name != 0 ) ClearDaughtersName();
      daughters_name = new G4String * [numberOfDaughters];
      for( G4int index = 0; index < numberOfDaughters; index++ )
	daughters_name[index] = new G4String( *right.daughters_name[index] );
    }

    pLambda = right.pLambda;
    pXi0 = right.pXi0;
  }
  return *this;
}
// -----------------------------------------------------------------------------------
void BooNEKaonDecayChannel::PhaseSpace( G4double parentM, const G4double * Mdaughter,
					G4double * Edaughter, G4double * Pdaughter )
{
  /*
  // first draw neutrino energy
  const G4double Q0 = 0.5 * ( parentM*parentM - Mdaughter[idPi]*Mdaughter[idPi] ) / parentM;
  G4double Enu = Q0 * G4UniformRand();

  // Then draw lepton energy
  G4double Elep_minus = Q0 - Enu;
  G4double Elep_plus  = 0.5 * ( parentM - Mdaughter[idPi]*Mdaughter[idPi] / (parentM - 2.0 * Enu) );
  G4double Elep = Elep_minus + (Elep_plus - Elep_minus) * G4UniformRand();

  // Rest goes to the pion
  G4double Epi = parentM - (Enu + Elep);

  Edaughter[idNeutrino] = Enu; Edaughter[idLepton] = Elep; Edaughter[idPi] = Epi;
  Pdaughter[idNeutrino] = std::sqrt( Enu*Enu - Mdaughter[idNeutrino]*Mdaughter[idNeutrino] );
  Pdaughter[idLepton] = std::sqrt( Enu*Enu - Mdaughter[idLepton]*Mdaughter[idLepton] );
  Pdaughter[idPi] = std::sqrt( Enu*Enu - Mdaughter[idPi]*Mdaughter[idPi] );
  */

  // Incorrect. Use the routine from G4KL3

  G4double sum_masses = 0.0;
  G4int index;
  const G4int N_DAUGHTER = 3;

  for( index = 0; index < N_DAUGHTER ; index++ ) sum_masses += Mdaughter[index];

  // Generate two daughters, add the third in later
  G4double rd1, rd2, rd;
  G4double pmax = 0.0, psum = 0.0;
  G4double energy;
  const size_t MAX_LOOP=10000;
  for( size_t loop_counter = 0; loop_counter < MAX_LOOP ; ++loop_counter ) {
    rd1 = G4UniformRand();
    rd2 = G4UniformRand();
    if (rd2 > rd1) {
      rd  = rd1;
      rd1 = rd2;
      rd2 = rd;
    } 
    pmax = 0.0;
    psum = 0.0;
    // daughter 0
    energy = rd2*(parentM - sum_masses);
    Pdaughter[0] = std::sqrt(energy*energy + 2.0*energy*Mdaughter[0]);
    Edaughter[0] = energy;
    if ( Pdaughter[0] >pmax )pmax =  Pdaughter[0];
    psum  +=  Pdaughter[0];
    // daughter 1
    energy = (1.-rd1)*(parentM - sum_masses);
    Pdaughter[1] = std::sqrt(energy*energy + 2.0*energy*Mdaughter[1]);
    Edaughter[1] = energy;
    if ( Pdaughter[1] >pmax )pmax =  Pdaughter[1];
    psum  +=  Pdaughter[1];
    // daughter 2
    energy = (rd1-rd2)*(parentM - sum_masses);
    Pdaughter[2] = std::sqrt(energy*energy + 2.0*energy*Mdaughter[2]);
    Edaughter[2] = energy;
    if ( Pdaughter[2] >pmax )pmax =  Pdaughter[2];
    psum  +=  Pdaughter[2];
    if (pmax <=  psum - pmax ) break;
  } 
}
// -----------------------------------------------------------------------------------
G4double BooNEKaonDecayChannel::DalitzDensity( G4double Mparent, 
					       G4double Epi, G4double El, G4double Enu,
					       G4double Mpi, G4double Ml )
{
  G4double Epi_prime = 0.5 * (Mparent*Mparent + Mpi*Mpi - Ml*Ml) / Mparent - Epi;
  G4double t = Mparent*Mparent + Mpi*Mpi - 2.0*Mparent*Epi;
  G4double form_fac = 1 + pLambda * t*t / (Mpi*Mpi);

  //G4cout << "ARGH. Epi, El, Enu = " << Epi << ", " << El << ", " << Enu << G4endl;

  // Now we need the right integrand, which depends if we're decaying to an electron or a muon
  if( Ml == fElectronMass ) {
    return (2.0 * El * Enu - Mparent * Epi_prime) * form_fac * form_fac;
  } else {
    // three terms in the quadratic expansion
    G4double A_fac = Mparent * (2.0 * El * Enu - Mparent * Epi_prime) + Ml*Ml * (0.25 * Epi_prime - Enu);
    G4double B_fac = Ml*Ml * (Enu - 0.5 * Epi_prime);
    G4double C_fac = Ml*Ml * 0.25 * Epi_prime;

    return ( A_fac + B_fac * pXi0 + C_fac * pXi0*pXi0 ) * form_fac * form_fac;
  }
}
// -----------------------------------------------------------------------------------
G4DecayProducts * BooNEKaonDecayChannel::DecayIt(G4double)
{

  // For accounting, you need to run CheckAndFillParent(), Daughters();
  CheckAndFillParent(); // fills G4MT_parent ParticleDefinition*
  CheckAndFillDaughters(); // fills G4MT_daughters ParticleDefinition*[]
  
  // Good old accept-reject.
  // Throw a random number, compare with maximum.
  // The maxima have been computed with Mathematica at the PDG values.. Analytically too expensive

  G4double MK = 0.0, Ml = 0.0, Mpi = 0.0;
  if( G4MT_parent == G4KaonPlus::KaonPlusDefinition() || 
      G4MT_parent == G4KaonMinus::KaonMinusDefinition() )
    { MK = fKPlusMass; Mpi = fPiZeroMass;  }
  else if( G4MT_parent == G4KaonZeroLong::KaonZeroLongDefinition() )
    { MK = fKZeroMass; Mpi = fPiPlusMass; }

  if( G4MT_daughters[idLepton] == G4Electron::ElectronDefinition() ||
      G4MT_daughters[idLepton] == G4Positron::PositronDefinition() )
    Ml = fElectronMass;
  else if( G4MT_daughters[idLepton] == G4MuonMinus::MuonMinusDefinition() ||
	   G4MT_daughters[idLepton] == G4MuonPlus::MuonPlusDefinition() )
    Ml = fMuonMass;

  G4double max_dens = 0.0;

  const G4int max_idx = 10000;
  const G4double Mdau[3] = { Mpi, Ml, 0.0 };
  G4double Edau[3] = { 0.0, 0.0, 0.0 }, Pdau[3] = { 0.0, 0.0, 0.0 };
  if( MK == fKPlusMass && Ml == fElectronMass ) max_dens = 0.026079877579886;
  else if( MK == fKPlusMass && Ml == fMuonMass ) max_dens = 0.011255195264726;
  else if( MK == fKZeroMass && Ml == fElectronMass ) max_dens = 0.026273645460311;
  else if( MK == fKZeroMass && Ml == fMuonMass ) max_dens = 0.011434838735247;
  else {
    G4cerr << "ERROR in BooNEKaonDecayChannel::DalitzDensity: Unrecognised parent with mass "
	   << MK << " and lepton with mass " << Ml << G4endl;
    G4cerr << "Returning now." << G4endl;
    return 0;
  }

  G4double w = -1.0; G4double r = G4UniformRand() * max_dens;
  G4int idx = -1; 
  
  // Loop to find a suitable point
  while( ! ( w > 0 && w > r ) && idx < max_idx ) {
    idx++;
    
    // First get the phase space density
    this->PhaseSpace( MK, &Mdau[0], &Edau[0], &Pdau[0] );
    // And add the masses of the daughters!
    Edau[0] += Mdau[0];
    Edau[1] += Mdau[1];
    Edau[2] += Mdau[2];

    // Next get the Dalitz density (assuming f+(0) = 1 as it'll cancel out)
    w = this->DalitzDensity( MK, Edau[idPi], Edau[idLepton], Edau[idNeutrino],
			     Mdau[idPi], Mdau[idLepton] );
    r = G4UniformRand() * max_dens;

    //G4cout << "w = " << w << " r = " << r << G4endl;
  }

  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << *daughters_name[0] << ":" << Pdau[0] << "[GeV/c]" <<G4endl;
    G4cout << *daughters_name[1] << ":" << Pdau[1] << "[GeV/c]" <<G4endl;
    G4cout << *daughters_name[2] << ":" << Pdau[2] << "[GeV/c]" <<G4endl;
  }
#endif

  // Turn into MeV - apparently that's what G4 wants...
  for( G4int midx = 0; midx < 3; midx++ ) {
    Pdau[midx] *= 1000.0;
    Edau[midx] *= 1000.0;
  }

  // Deal with kinematics.
  G4ThreeVector* direction = new G4ThreeVector(1.0,0.0,0.0);
  G4DynamicParticle * parent_particle = new G4DynamicParticle( G4MT_parent, *direction, 0.0 );
  //delete direction; // ok, bye now

  G4DecayProducts * products = new G4DecayProducts(*parent_particle);
  //delete parent_particle; // you too? ok, bye

  G4double cthPi, sthPi, phiPi, cphiPi, sphiPi;
  G4double cthNu, sthNu, phiNu, cphiNu, sphiNu;

  // rehashed from G4KL3DecayChannel

  // Pion kinematics
  
  cthPi  = 2.0 * G4UniformRand() - 1.0;
  sthPi  = std::sqrt( 1.0 - cthPi * cthPi );
  phiPi  = twopi * G4UniformRand();
  sphiPi = std::sin(phiPi);
  cphiPi = std::cos(phiPi);

  direction = new G4ThreeVector(sthPi * cphiPi, sthPi * sphiPi, cthPi);
  G4ThreeVector p3Pi = (*direction) * Pdau[idPi];
  std::unique_ptr<G4DynamicParticle> daughter_particle = std::make_unique<G4DynamicParticle>( G4DynamicParticle( G4MT_daughters[idPi], p3Pi ) );
  products->PushProducts( new G4DynamicParticle( *daughter_particle ) );
  // this was not in original code but will leak memory otherwise
  //delete daughter_particle; daughter_particle = 0;
  //G4cout << "AXA - p3 pi = " << daughter_particle->GetTotalMomentum() << " MeV" << G4endl;
  //products->DumpInfo();

  // Neutrino kinematics

  cthNu  = (Pdau[idLepton]*Pdau[idLepton] -
	    Pdau[idNeutrino]*Pdau[idNeutrino] -
	    Pdau[idPi]*Pdau[idPi]) / (2.0 * Pdau[idNeutrino] * Pdau[idPi]);
  sthNu  = std::sqrt(1.0 - cthNu*cthNu);
  phiNu  = twopi * G4UniformRand();
  sphiNu = std::sin(phiNu);
  cphiNu = std::cos(phiNu);

  // I trust ye...
  direction->setX( sthNu * cphiNu * cthPi * cphiPi -
		   sthNu * sphiNu * sphiPi +
		   cthNu * sthPi * cphiPi );
  direction->setY( sthNu * cphiNu * cthPi * sphiPi +
		   sthNu * sphiNu * cphiPi +
		   cthNu * sthPi * sphiPi );
  direction->setZ( -sthNu * cphiNu * sthPi + cthNu * cthPi );

  G4ThreeVector p3Nu = (*direction) * Pdau[idNeutrino];
  daughter_particle = std::make_unique<G4DynamicParticle>( G4DynamicParticle( G4MT_daughters[idNeutrino], p3Nu ) );
  products->PushProducts( new G4DynamicParticle( *daughter_particle ) );
  //delete daughter_particle; daughter_particle = 0;

  // Lepton kinematics
  G4ThreeVector p3Lep = (p3Pi + p3Nu) * (-1.0);
  daughter_particle = std::make_unique<G4DynamicParticle>( G4DynamicParticle( G4MT_daughters[idLepton], p3Lep ) );
  products->PushProducts( new G4DynamicParticle( *daughter_particle ) );
  //delete daughter_particle; daughter_particle = 0;

  #ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4KL3DecayChannel::DecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     G4cout << "  decay products address=" << products << G4endl;
     products->DumpInfo();
  }
#endif
  delete direction; direction = 0;
  return products;
}
// -----------------------------------------------------------------------------------
