#include "BooNEHadronElasticPhysics.hh"
#include "BooNEHadronElasticProcess.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"


BooNEHadronElasticPhysics::BooNEHadronElasticPhysics(G4int)
  :  G4VPhysicsConstructor("BooNEHadronElasticPhysics")
{}

BooNEHadronElasticPhysics::~BooNEHadronElasticPhysics()
{
}

void BooNEHadronElasticPhysics::ConstructParticle()
{

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  return; 
}


void BooNEHadronElasticPhysics::ConstructProcess()
{  
  G4ProcessManager* pManager = 0;
  auto aParticleIterator = GetParticleIterator();
  aParticleIterator->reset();
  while( (*aParticleIterator)() ) {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String pname=particle->GetParticleName();
    if(pname == "proton"   ||
       pname == "anti_proton"   ||
       pname == "neutron"   ||
       pname == "anti_neutron"   ||
       pname == "pi+"   ||
       pname == "pi-"   ||
       pname == "pi0"   ||
       pname == "kaon+"   ||
       pname == "kaon-"   ||
       pname == "kaon0"   ||
       pname == "kaon0L"   ||
       pname == "kaon0S"   ||
       pname == "anti_lambda"  ||
       pname == "anti_neutron" ||
       pname == "anti_omega-"  || 
       pname == "anti_sigma-"  || 
       pname == "anti_sigma+"  || 
       pname == "anti_xi-"  || 
       pname == "anti_xi0"  || 
       pname == "lambda"    || 
       pname == "omega-"    || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton"   
       ) {
      pManager = particle->GetProcessManager();
      pManager->AddDiscreteProcess(&fBooNEHadronElasticProcess);
    }
  }
  fBooNEHadronElasticProcess.RegisterMe(&fBooNEHadronElasticModel);
  fBooNEHadronElasticProcess.AddDataSet(&fBooNEHadronElasticData);

}
