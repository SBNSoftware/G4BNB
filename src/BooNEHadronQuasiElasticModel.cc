#include  "BooNEHadronQuasiElasticModel.hh"
#include  <sstream>
#include  <iostream>
#include  <fstream>
#include  "G4ios.hh"

using namespace std;

BooNEHadronQuasiElasticModel::BooNEHadronQuasiElasticModel()       
{     }



G4HadFinalState*  
BooNEHadronQuasiElasticModel::ApplyYourself(const  G4HadProjectile  &aHadron,
 					           G4Nucleus        &aNucleus)
{

  // Everything here is done in GeV
  theParticleChange.Clear();

  // Get information on the initial statenucleus
  G4Nucleus * theNucleus      = &aNucleus;
  G4int nN                    = (int) theNucleus->GetN_asInt();
  G4int nZ                    = (int) theNucleus->GetZ_asInt();  

  // Get information on the intial state hadron
  const G4HadProjectile * theHadron   = &aHadron;

  G4double    hadrMass    = theHadron->GetDefinition()->GetPDGMass()/CLHEP::GeV;
  G4String    hadrName    = theHadron->GetDefinition()->GetParticleName();
  G4double    nucleonMass = G4Proton::Proton()->GetPDGMass()/CLHEP::GeV;
  
  G4ThreeVector    hadrMomentum = theHadron->Get4Momentum().vect()/CLHEP::GeV;
  G4LorentzVector  oldHadrMomentum = theHadron->Get4Momentum()/CLHEP::GeV;
  
   
  G4double   inLabMom  = theHadron->GetTotalMomentum()/CLHEP::GeV; 
  G4double   inEneHadr = theHadron->GetTotalEnergy()/CLHEP::GeV;   

  // Get kinematics of outgoing hadron: 
  // calculations are done with collision axis along z in the CM frame
  // first in CM frame
  G4double   s         = hadrMass * hadrMass +  nucleonMass*nucleonMass + 2*inEneHadr*nucleonMass;
  G4double   outCMMom  = inLabMom * nucleonMass / sqrt(s);

  // draw a q^2 (GeV^2), phi based on the beta
   G4double q2          = 0;

  // draw q2 until it makes sense (i.e. q2 is within kinematic limits)
  // try 100 times and give up, in which case proceed with q2 = 0
  G4double test         = 2;
  G4int    itest        = 0;
  G4double ran1         = 0;

  while(test > 1.0){
    ran1  = G4UniformRand();
    q2    = -std::log(ran1)/Beta(theHadron);
    test  = sqrt(q2/(4*outCMMom*outCMMom));
    itest++;
    if(itest == 100){
      test = 0;
      G4cout << "BooNEHadronQuasiElasticModel: "
	     << " Drew 100 times but couldn't get acceptable q^2."
	     << " Do not Scatter:set q^2 = 0" << G4endl;
     }
  }


  G4double ran2        = G4UniformRand();
  G4double cmsAngle    = 2 * asin( test );
  G4double cmsPhi      = CLHEP::twopi*ran2;
 

  // create a lorentz vector and boost to lab frame
  // outgoing hadron
  G4LorentzVector hadronOutMomentum(outCMMom*sin(cmsAngle)*cos(cmsPhi), 
				    outCMMom*sin(cmsAngle)*sin(cmsPhi), 
				    outCMMom*cos(cmsAngle),
				    sqrt(outCMMom*outCMMom + hadrMass*hadrMass));
  // outgoing nucleon
  G4LorentzVector nucleonOutMomentum(-outCMMom*sin(cmsAngle)*cos(cmsPhi), 
				     -outCMMom*sin(cmsAngle)*sin(cmsPhi), 
				     -outCMMom*cos(cmsAngle),
				     sqrt(outCMMom*outCMMom + nucleonMass*nucleonMass));
  
  //get boost vector
  G4double      boostHadNucleon       = inLabMom/(inEneHadr + nucleonMass);
  G4ThreeVector boostHadNucleonVector = hadrMomentum.unit();

  // boost the momentum
  hadronOutMomentum.boost(boostHadNucleonVector,  boostHadNucleon);
  nucleonOutMomentum.boost(boostHadNucleonVector, boostHadNucleon);

  //   GEANT4 hands you the incident particle with its momentum along the z-axis.
  //   All final state kinematics are performed according to this frame.
  //   GEANT4 then automatically rotates the particles back to the lab frame
  //   so no explicit rotation back to the lab frame is necessary.
  //   Even if you tried, you don't have the direction of the incident particle
  //   in the lab frame anyways!

  // !!!! Removed by FRAN  
  // Define the secondaries necessary for the final state
  // first set final state of hadron, now back to MeV 
  // Kinetic energy
  // ///theParticleChange.SetEnergyChange( (hadronOutMomentum.e() - hadrMass) * CLHEP::GeV);

  // Direction
  // //// theParticleChange.SetMomentumChange(hadronOutMomentum.vect().unit());

  // !!! FRAN Add secondary for QE-scattered hadron
  G4DynamicParticle * outHadron = 
    new G4DynamicParticle(theHadron->GetDefinition(), hadronOutMomentum.vect() * CLHEP::GeV);
  theParticleChange.AddSecondary(outHadron);
  
  // decide if we scattered off at proton or a neutron
  // we'll use a simple scaling based on Z protons and (A-Z) neutrons;
  G4double protonRatio = (G4double)nZ/(G4double)(nN);
  G4double random3 = G4UniformRand();

  if(random3 < protonRatio){
    // proton
    G4DynamicParticle * outNucleon = 
      new G4DynamicParticle(G4Proton::Proton(), nucleonOutMomentum.vect() * CLHEP::GeV);
    theParticleChange.AddSecondary(outNucleon);
    
    // Do not add recoil nucleus
    // G4ParticleDefinition * outNucleusDef = G4ParticleTable::GetParticleTable()->FindIon(nZ-1,nN-1,0); 
    // G4ThreeVector outNucleusMomentum(0,0,0);
    // G4DynamicParticle * outNucleus =   new G4DynamicParticle(outNucleusDef, outNucleusMomentum * GeV);
    // theParticleChange.AddSecondary(outNucleus);

  } else {
    // neutron
    G4DynamicParticle * outNucleon =   
      new G4DynamicParticle(G4Neutron::Neutron(), nucleonOutMomentum.vect() * CLHEP::GeV);
    theParticleChange.AddSecondary(outNucleon);

    // Do not add recoil nucleus
    // G4ParticleDefinition * outNucleusDef = G4ParticleTable::GetParticleTable()->FindIon(nZ,nN-1,0); 
    // G4ThreeVector outNucleusMomentum(0,0,0);
    // G4DynamicParticle * outNucleus =   new G4DynamicParticle(outNucleusDef, outNucleusMomentum);

  }

  // !!! FRAN Status set to stop the incident hadron
  theParticleChange.SetStatusChange( stopAndKill );


  return &theParticleChange;

 }


G4double BooNEHadronQuasiElasticModel::Beta(const G4HadProjectile* theIncidentHadron)
{

  // Return slope of h-(p/n) elastic scattering (Beta) 
  // according to incident hadron type and momentum in units of (c/GeV)^2
  G4double betaParameter = 0;
  G4double p = theIncidentHadron->GetTotalMomentum()/CLHEP::GeV; 

  // get a random number
  G4double ran = G4UniformRand();  

  G4double tNorm, tOffset, tSlope, gNorm, gMean, gWidth;
  G4double p0, p1, p2, p3, p4, p5;

  // assuming isoscalar target (not quite true) and isospin symmetry, 
  // we'll randomly use proton/neutron or pi+/pi- parameters

  if( theIncidentHadron->GetDefinition() == G4Proton::Proton() ||
      theIncidentHadron->GetDefinition() == G4Neutron::Neutron() ){
    if(ran < 0.5){
      p0  = protonBetaParameters[0];
      p1  = protonBetaParameters[1];
      p2  = protonBetaParameters[2];
      p3  = protonBetaParameters[3];
      p4  = protonBetaParameters[4];
      p5  = protonBetaParameters[5];
    } else {
      p0  = neutronBetaParameters[0];
      p1  = neutronBetaParameters[1];
      p2  = neutronBetaParameters[2];
      p3  = neutronBetaParameters[3];
      p4  = neutronBetaParameters[4];
      p5  = neutronBetaParameters[5];
    }

    betaParameter = p0 + p1*p + p2*pow(p,2) + p3* pow(p,3) + p4 *pow(p,4) + p5 * pow(p,5);
  }
  else if( theIncidentHadron->GetDefinition() == G4PionPlus::PionPlus() ||
	   theIncidentHadron->GetDefinition() == G4PionMinus::PionMinus() ){
    if(ran < 0.5){
      tNorm   = pionPlusBetaParameters[0];
      tOffset = pionPlusBetaParameters[1];
      tSlope  = pionPlusBetaParameters[2];
      gNorm   = pionPlusBetaParameters[3];
      gMean   = pionPlusBetaParameters[4];
      gWidth  = pionPlusBetaParameters[5];
    } else {
      tNorm   = pionMinusBetaParameters[0];
      tOffset = pionMinusBetaParameters[1];
      tSlope  = pionMinusBetaParameters[2];
      gNorm   = pionMinusBetaParameters[3];
      gMean   = pionMinusBetaParameters[4];
      gWidth  = pionMinusBetaParameters[5];
    }

    betaParameter = TanhWithGaussian( p, 
				      tNorm,  tOffset,  tSlope, 
				      gNorm,  gMean,    gWidth);

  }


  return betaParameter;

}


G4double BooNEHadronQuasiElasticModel::TanhWithGaussian(G4double p, 
							G4double tNorm, G4double tOffset, G4double tSlope, 
							G4double gNorm, G4double gMean,   G4double gWidth){

  G4double Tanh  = tNorm * ( 1 + tanh( tSlope*(p - tOffset) ) ); 
  G4double Gauss = ( gNorm/sqrt(2 * M_PI * gWidth * gWidth) ) * exp(-(p-gMean)*(p-gMean)/(2*gWidth*gWidth) );

  return Tanh+Gauss;
}


/*  End of file  */
