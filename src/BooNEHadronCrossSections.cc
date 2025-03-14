#include "BooNEHadronCrossSections.hh"
#include "G4ios.hh"

// Initialize static pointer for singleton instance
BooNEHadronCrossSections* 
BooNEHadronCrossSections::theInstance = 0;

BooNEHadronCrossSections::BooNEHadronCrossSections()
{
  theCrossSectionsMessenger = new BooNEHadronCrossSectionsMessenger(this);
}

G4bool 
BooNEHadronCrossSections::IsApplicable(const G4DynamicParticle* aParticle,  
				       G4int Z, G4int A)
{
  // currently only Beryllium and Aluminum
  bool isBe =  ( Z ==  4 ) && ( A ==  9) ;
  bool isAl =  ( Z == 13 ) && ( A == 27) ;

  // return if neither
  if( !(isBe || isAl) ) return false;

  // get the particle momentum to ensure that we are in the valide momentum range
  G4double momentum = aParticle->GetTotalMomentum()/CLHEP::GeV;

  // currently only proton, neutron, piplus and piminus
  G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();

  // Beryllium
  if(aParticleDef == G4Proton::Proton()       && isBe)return (momentum > pBeMin   && momentum < pBeMax  );
  if(aParticleDef == G4Neutron::Neutron()     && isBe)return (momentum > nBeMin   && momentum < nBeMax  );
  if(aParticleDef == G4PionPlus::PionPlus()   && isBe)return (momentum > pipBeMin && momentum < pipBeMax);
  if(aParticleDef == G4PionMinus::PionMinus() && isBe)return (momentum > pimBeMin && momentum < pimBeMax);

  // Aluminum
  if(aParticleDef == G4Proton::Proton()       && isAl)return (momentum > pAlMin   && momentum < pAlMax  );
  if(aParticleDef == G4Neutron::Neutron()     && isAl)return (momentum > nAlMin   && momentum < nAlMax  );
  if(aParticleDef == G4PionPlus::PionPlus()   && isAl)return (momentum > pipAlMin && momentum < pipAlMax);
  if(aParticleDef == G4PionMinus::PionMinus() && isAl)return (momentum > pimAlMin && momentum < pimAlMax);

  // if we're here, we ain't got it
  return false;

}

G4double
BooNEHadronCrossSections::GetTotalCrossSection(const G4DynamicParticle* aParticle,
					       G4int ZNucleus, G4int ANucleus)
{

  G4double theCrossSection = 0;
  G4double theMomentum     = aParticle->GetTotalMomentum()/CLHEP::GeV;

  // go through the various cases:
  if(aParticle->GetDefinition()      == G4Proton::Proton()){

    // Proton total cross section is based on Regge parametrization
    G4double A, B, n, C, D;

    if( ZNucleus == 4 && ANucleus == 9){
      A = pBeTotParameters[0];
      B = pBeTotParameters[1];
      n = pBeTotParameters[2];
      C = pBeTotParameters[3];
      D = pBeTotParameters[4];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      A = pAlTotParameters[0];
      B = pAlTotParameters[1];
      n = pAlTotParameters[2];
      C = pAlTotParameters[3];
      D = pAlTotParameters[4];
    }
    theCrossSection =  Regge(theMomentum, A, B, n, C, D);

  }

  else if(aParticle->GetDefinition() == G4Neutron::Neutron()){

    // Neutron total cross section is based on Regge parametrization
    G4double A, B, n, C, D;

    if( ZNucleus == 4 && ANucleus == 9){
      A = nBeTotParameters[0];
      B = nBeTotParameters[1];
      n = nBeTotParameters[2];
      C = nBeTotParameters[3];
      D = nBeTotParameters[4];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      A = nAlTotParameters[0];
      B = nAlTotParameters[1];
      n = nAlTotParameters[2];
      C = nAlTotParameters[3];
      D = nAlTotParameters[4];
    }
    theCrossSection =  Regge(theMomentum, A, B, n, C, D);
  }

  else if(aParticle->GetDefinition() == G4PionPlus::PionPlus()){

    // piplus total invovles a Caroll et al.Breit Wigner and a Regge parametrization with threshold
    // however, instead of summing the two parts, we have a transition parameter that switches
    // between the two functional forms. We'll call it tr
    G4double tr, p0, s, A, B, C, n, D;
    
    if( ZNucleus == 4 && ANucleus == 9){
      tr = pipBeTotParameters[0];
      p0 = pipBeTotParameters[1];
      s  = pipBeTotParameters[2];
      A  = pipBeTotParameters[3];
      B  = pipBeTotParameters[4];
      C  = pipBeTotParameters[5];
      n  = pipBeTotParameters[6];
      D  = pipBeTotParameters[7];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      tr    = pipAlTotParameters[0];
      p0    = pipAlTotParameters[1];
      s     = pipAlTotParameters[2];
      A     = pipAlTotParameters[3];
      B     = pipAlTotParameters[4];
      C     = pipAlTotParameters[5];
      n     = pipAlTotParameters[6];
      D     = pipAlTotParameters[7];
    }

    if( theMomentum < tr)
      theCrossSection =  CarrollBreitWigner(theMomentum, ANucleus, G4PionPlus::PionPlus()->GetPDGMass()/CLHEP::GeV);
    else 
      theCrossSection = ReggeWithThreshold(theMomentum, p0, s, A, B, C, n, D);


  }

  else if(aParticle->GetDefinition() == G4PionMinus::PionMinus()){

    // piminus total invovles a Caroll et al.Breit Wigner and a Regge parametrization with threshold
    // however, instead of summing the two parts, we have a transition parameter that switches
    // between the two functional forms. We'll call it tr
    G4double tr, p0, s, A, B, C, n, D; 

    if( ZNucleus == 4 && ANucleus == 9){
      tr = pimBeTotParameters[0];
      p0 = pimBeTotParameters[1];
      s  = pimBeTotParameters[2];
      A  = pimBeTotParameters[3];
      B  = pimBeTotParameters[4];
      C  = pimBeTotParameters[5];
      n  = pimBeTotParameters[6];
      D  = pimBeTotParameters[7];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      tr    = pimAlTotParameters[0];
      p0    = pimAlTotParameters[1];
      s     = pimAlTotParameters[2];
      A     = pimAlTotParameters[3];
      B     = pimAlTotParameters[4];
      C     = pimAlTotParameters[5];
      n     = pimAlTotParameters[6];
      D     = pimAlTotParameters[7];
    }


    if( theMomentum < tr)
      theCrossSection =  CarrollBreitWigner(theMomentum, ANucleus, G4PionMinus::PionMinus()->GetPDGMass()/CLHEP::GeV);
    else 
      theCrossSection = ReggeWithThreshold(theMomentum, p0, s, A, B, C, n, D);
  }

  // add units
  return theCrossSection * CLHEP::millibarn;

}


G4double
BooNEHadronCrossSections::GetElasticCrossSection(const G4DynamicParticle* aParticle,
						 G4int Z, G4int A)
{

  // we define the elastic cross section 
  // as the difference between total and inelastic
  G4double theCrossSection = GetTotalCrossSection(aParticle, Z, A) - 
    GetInelasticCrossSection(aParticle, Z, A);
  //if we are not subtracting QE from inelastic then we should not do it here either 
  //see comment under GetInelasticCrossSection
  //  -GetQuasiElasticCrossSection(aParticle, Z, A);

  return theCrossSection;
}



G4double
BooNEHadronCrossSections::GetQuasiElasticCrossSection(const G4DynamicParticle* aParticle,
						      G4int ZNucleus, G4int ANucleus)
{

  G4double theCrossSection = 0;
  G4double theMomentum     = aParticle->GetTotalMomentum()/CLHEP::GeV;

  // go through the various cases:
  if(aParticle->GetDefinition()      == G4Proton::Proton()){

    // Proton quasi-elastic cross section is based on Regge parametrization
    G4double A, B, n, C, D;

    if( ZNucleus == 4 && ANucleus == 9){
      A = pBeQelParameters[0];
      B = pBeQelParameters[1];
      n = pBeQelParameters[2];
      C = pBeQelParameters[3];
      D = pBeQelParameters[4];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      A = pAlQelParameters[0];
      B = pAlQelParameters[1];
      n = pAlQelParameters[2];
      C = pAlQelParameters[3];
      D = pAlQelParameters[4];
    }
    theCrossSection =  Regge(theMomentum, A, B, n, C, D);

  }

  else if(aParticle->GetDefinition() == G4Neutron::Neutron()){

    // Neutron quasi-elastic cross section is based on Regge parametrization
    G4double A, B, n, C, D;

    if( ZNucleus == 4 && ANucleus == 9){
      A = nBeQelParameters[0];
      B = nBeQelParameters[1];
      n = nBeQelParameters[2];
      C = nBeQelParameters[3];
      D = nBeQelParameters[4];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      A = nAlQelParameters[0];
      B = nAlQelParameters[1];
      n = nAlQelParameters[2];
      C = nAlQelParameters[3];
      D = nAlQelParameters[4];
    }
    theCrossSection =  Regge(theMomentum, A, B, n, C, D);
  }

  else if(aParticle->GetDefinition() == G4PionPlus::PionPlus()){

    // piplus quasi-elastic invovles a Breit Wigner and a Regge parametrization with threshold

    G4double N, M0, Gamma, p0, s, A, B, C, n, D;
    
    if( ZNucleus == 4 && ANucleus == 9){
      N     = pipBeQelParameters[0];
      M0    = pipBeQelParameters[1];
      Gamma = pipBeQelParameters[2];
      p0    = pipBeQelParameters[3];
      s     = pipBeQelParameters[4];
      A     = pipBeQelParameters[5];
      B     = pipBeQelParameters[6];
      C     = pipBeQelParameters[7];
      n     = pipBeQelParameters[8];
      D     = pipBeQelParameters[9];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      N     = pipAlQelParameters[0];
      M0    = pipAlQelParameters[1];
      Gamma = pipAlQelParameters[2];
      p0    = pipAlQelParameters[3];
      s     = pipAlQelParameters[4];
      A     = pipAlQelParameters[5];
      B     = pipAlQelParameters[6];
      C     = pipAlQelParameters[7];
      n     = pipAlQelParameters[8];
      D     = pipAlQelParameters[9];
    }

    theCrossSection =  BreitWigner(theMomentum, 
				   G4PionPlus::PionPlus()->GetPDGMass()/CLHEP::GeV, 
				   G4Proton::Proton()->GetPDGMass()/CLHEP::GeV, 
				   N, M0, Gamma) + 
                       ReggeWithThreshold(theMomentum, p0, s, A, B, C, n, D);
  }

  else if(aParticle->GetDefinition() == G4PionMinus::PionMinus()){
    // piminus quasi-elastic invovles a Breit Wigner and a Regge parametrization with threshold
    G4double N, M0, Gamma, p0, s, A, B, C, n, D;
    
    if( ZNucleus == 4 && ANucleus == 9){
      N     = pimBeQelParameters[0];
      M0    = pimBeQelParameters[1];
      Gamma = pimBeQelParameters[2];
      p0    = pimBeQelParameters[3];
      s     = pimBeQelParameters[4];
      A     = pimBeQelParameters[5];
      B     = pimBeQelParameters[6];
      C     = pimBeQelParameters[7];
      n     = pimBeQelParameters[8];
      D     = pimBeQelParameters[9];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      N     = pimAlQelParameters[0];
      M0    = pimAlQelParameters[1];
      Gamma = pimAlQelParameters[2];
      p0    = pimAlQelParameters[3];
      s     = pimAlQelParameters[4];
      A     = pimAlQelParameters[5];
      B     = pimAlQelParameters[6];
      C     = pimAlQelParameters[7];
      n     = pimAlQelParameters[8];
      D     = pimAlQelParameters[9];
    }
    theCrossSection =  BreitWigner(theMomentum, 
				   G4PionMinus::PionMinus()->GetPDGMass()/CLHEP::GeV, 
				   G4Proton::Proton()->GetPDGMass()/CLHEP::GeV, 
				   N, M0, Gamma) +       
                       ReggeWithThreshold(theMomentum, p0, s, A, B, C, n, D);
 
  }

  // add units
  return theCrossSection * CLHEP::millibarn;
}

G4double
BooNEHadronCrossSections::GetInelasticCrossSection(const G4DynamicParticle* aParticle,
						   G4int ZNucleus, G4int ANucleus)
{

  // note strange situation here:
  // the G4 inelastic models typically do not include quasi elastic interactions so
  // we will subtract the quasi elastic cross section from whatever we get here
  // under the assumption that the quasi-elastic process gets assigned to its own model

  // Actually with geant4.10 we will let geant handle quasielastic within inelastic model, 
  // so we will not subtract quasielastic part
  G4double theCrossSection             = 0;
  //  G4double theQuasiElasticCrossSection = GetQuasiElasticCrossSection(aParticle,ZNucleus, ANucleus);

  G4double theMomentum     = aParticle->GetTotalMomentum()/CLHEP::GeV;

  // go through the various cases:
  if(aParticle->GetDefinition()      == G4Proton::Proton()){

    // Proton inelastic cross section is based on Regge parametrization
    G4double A, B, n, C, D;

    if( ZNucleus == 4 && ANucleus == 9){
      A = pBeIneParameters[0];
      B = pBeIneParameters[1];
      n = pBeIneParameters[2];
      C = pBeIneParameters[3];
      D = pBeIneParameters[4];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      A = pAlIneParameters[0];
      B = pAlIneParameters[1];
      n = pAlIneParameters[2];
      C = pAlIneParameters[3];
      D = pAlIneParameters[4];
    }
    theCrossSection =  Regge(theMomentum, A, B, n, C, D);

  }

  else if(aParticle->GetDefinition() == G4Neutron::Neutron()){

    // Neutron inelastic cross section is based on Regge parametrization
    G4double A, B, n, C, D;

    if( ZNucleus == 4 && ANucleus == 9){
      A = nBeIneParameters[0];
      B = nBeIneParameters[1];
      n = nBeIneParameters[2];
      C = nBeIneParameters[3];
      D = nBeIneParameters[4];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      A = nAlIneParameters[0];
      B = nAlIneParameters[1];
      n = nAlIneParameters[2];
      C = nAlIneParameters[3];
      D = nAlIneParameters[4];
    }
    theCrossSection =  Regge(theMomentum, A, B, n, C, D);
  }

  else if(aParticle->GetDefinition() == G4PionPlus::PionPlus()){

    // piplus inelastic invovles a Breit Wigner and a Regge parametrization with threshold
    G4double N, M0, Gamma, p0, s, A, B, C, n, D;
    
    if( ZNucleus == 4 && ANucleus == 9){
      N     = pipBeIneParameters[0];
      M0    = pipBeIneParameters[1];
      Gamma = pipBeIneParameters[2];
      p0    = pipBeIneParameters[3];
      s     = pipBeIneParameters[4];
      A     = pipBeIneParameters[5];
      B     = pipBeIneParameters[6];
      C     = pipBeIneParameters[7];
      n     = pipBeIneParameters[8];
      D     = pipBeIneParameters[9];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      N     = pipAlIneParameters[0];
      M0    = pipAlIneParameters[1];
      Gamma = pipAlIneParameters[2];
      p0    = pipAlIneParameters[3];
      s     = pipAlIneParameters[4];
      A     = pipAlIneParameters[5];
      B     = pipAlIneParameters[6];
      C     = pipAlIneParameters[7];
      n     = pipAlIneParameters[8];
      D     = pipAlIneParameters[9];
    }
    theCrossSection =  BreitWigner(theMomentum, 
				   G4PionPlus::PionPlus()->GetPDGMass()/CLHEP::GeV, 
				   G4Proton::Proton()->GetPDGMass()/CLHEP::GeV, 
				   N, M0, Gamma) + 
                       ReggeWithThreshold(theMomentum, p0, s, A, B, C, n, D);

  }

  else if(aParticle->GetDefinition() == G4PionMinus::PionMinus()){
  // piminus inelastic invovles a Breit Wigner and a Regge parametrization with threshold
    G4double N, M0, Gamma, p0, s, A, B, C, n, D;
    
    if( ZNucleus == 4 && ANucleus == 9){
      N     = pimBeIneParameters[0];
      M0    = pimBeIneParameters[1];
      Gamma = pimBeIneParameters[2];
      p0    = pimBeIneParameters[3];
      s     = pimBeIneParameters[4];
      A     = pimBeIneParameters[5];
      B     = pimBeIneParameters[6];
      C     = pimBeIneParameters[7];
      n     = pimBeIneParameters[8];
      D     = pimBeIneParameters[9];
    } 
    else if( ZNucleus == 13 && ANucleus == 27){
      N     = pimAlIneParameters[0];
      M0    = pimAlIneParameters[1];
      Gamma = pimAlIneParameters[2];
      p0    = pimAlIneParameters[3];
      s     = pimAlIneParameters[4];
      A     = pimAlIneParameters[5];
      B     = pimAlIneParameters[6];
      C     = pimAlIneParameters[7];
      n     = pimAlIneParameters[8];
      D     = pimAlIneParameters[9];
    }


    theCrossSection =  BreitWigner(theMomentum, 
				   G4PionMinus::PionMinus()->GetPDGMass()/CLHEP::GeV, 
				   G4Proton::Proton()->GetPDGMass()/CLHEP::GeV, 
				   N, M0, Gamma) + 
                       ReggeWithThreshold(theMomentum, p0, s, A, B, C, n, D);
  }

  // add units
  return theCrossSection * CLHEP::millibarn;// - theQuasiElasticCrossSection;

}


// Routines: all assume p is in GeV/c (c=1) and return cross sections in mb.   

// Breit Wigner functions
G4double BooNEHadronCrossSections::
BreitWigner(G4double p, G4double mHadron, G4double mTarget, G4double N, G4double M0, G4double Gamma)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  G4double EHadron     = sqrt(mHadron*mHadron + p*p);
  G4double m           = sqrt( pow(EHadron+mTarget,2) - pow(p,2) );

  G4double numerator   = N * pow(m*Gamma, 2);
  G4double denominator = pow(M0*M0-m*m, 2) + pow(m*Gamma, 2);
  G4double bw          = numerator/denominator;

  return bw;
#pragma GCC diagnostic pop
}

// note this function works for Be and Al only!
G4double BooNEHadronCrossSections::
CarrollBreitWigner(G4double p,  G4double A,  G4double mHadron)
{
  // this nasty function reproduces equation 2 of Carroll et al.

  // hadron kinetic energy in GeV
  G4double hadronKE = sqrt( pow(p, 2) + pow(mHadron, 2) ) - mHadron;

  // hbarc in GeV fm
  G4double myhbarc  = 0.197327;

  // All masses in GeV (c=1)
  G4double mAMU     = 0.931494;
  G4double mNucleus = A*mAMU;
  G4double mNucleon = 0.938;

  // radius of nucleus
  G4double rNucleus = 0;
  if( A == 9){
    rNucleus = 3.12;
  } else if (A == 27){
    rNucleus = 3.76;
  }

  // Get the pion momentum in the Center of Mass of the pi-Nucleus system
  G4double eHadron = hadronKE + mHadron; 
  G4double s       = pow(mHadron, 2) + pow(mNucleus, 2) + (2 * mNucleus * eHadron);
  G4double pCM     = p * mNucleus / sqrt(s);
  G4double kCM     = pCM / myhbarc;

  // lbar in fermi
  G4double lbarCM = 1.0/kCM;
  
  G4double Lmax = (kCM * rNucleus) + 1.0;


  // A-dependent parameters from Carroll et al.
  G4double E0     = 0.001 * ( 1227.0 -16.0 * pow(A,  0.33333) );
  G4double E1     = 0.001 * ( 15.0         * pow(A, -0.66667) );

  // override parametrization for Be; doesn't seem to work well here
  if(A== 9) E1 = 0.001 * 2.5;

  G4double G10    = 0.001 * (-67.0 + 54.0 * pow(A, 0.33333) );
  G4double Gamma2 = 0.001 * (- 9.0 + 55.0 * pow(A, 0.33333) );

  // CM momentum corresponding to E0
  G4double p0     =  sqrt( pow(E0 - mNucleon, 2) - pow(mHadron,2) );

  // construct the arguments of the atan function 
  G4double E = sqrt(pCM*pCM + mHadron *mHadron )  + 
               sqrt(pCM*pCM + mNucleus*mNucleus)  - mNucleus + mNucleon;


  G4double Gamma1 = pCM  * (G10/p0);
  G4double Gamma  = Gamma1 + Gamma2;

  G4double coeff = (2 * M_PI * lbarCM*lbarCM * Gamma1)/E1;
  G4double arg1  = (E0 + (E1 * Lmax * Lmax) - E)/ (0.5 * Gamma);
  G4double arg2  = (E0 - E)/ (0.5 * Gamma);

  // the factor of ten is to translate fm^2 to mb
  return 10.0 * coeff * (atan(arg1) - atan(arg2));

}


// Regge parametrization of high energy behavior
G4double BooNEHadronCrossSections::
Regge(G4double p, G4double A, G4double B, G4double n, G4double C, G4double D)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  G4double regge = A + B * pow(p,n) + C * pow(log(p), 2) + D * log(p);
  return regge;
#pragma GCC diagnostic pop
}

// Regge parameterization with threshold
G4double BooNEHadronCrossSections::
ReggeWithThreshold(G4double p, G4double p0, G4double s, 
		   G4double A, G4double B, G4double n, G4double C, G4double D)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  G4double threshold = 0.5 * (1 + tanh( s*(p - p0) ) );
  G4double regge     = Regge(p, A, B, n, C, D);
  return threshold * regge;
#pragma GCC diagnostic pop
}

