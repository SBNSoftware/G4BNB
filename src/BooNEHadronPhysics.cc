#include "BooNEHadronPhysics.hh"
#include "BooNEpBeInteraction.hh"
#include "BooNEHadronZeroModel.hh"
#include "BooNEHadronZeroCrossSection.hh"
#include "BooNEHadronInelasticProcess.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"
#include "G4HadronElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4Material.hh"

#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4QuasiElasticChannel.hh"

#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"
#include "G4HyperonFTFPBuilder.hh"

BooNEHadronPhysics::BooNEHadronPhysics(G4int)
  :  G4VPhysicsConstructor("BooNEPhysics")
{}

BooNEHadronPhysics::~BooNEHadronPhysics()
{
  delete fHyperonBuilder;
  delete fAntiBarionBuilder;
  delete fFTFPAntiBarion;
  delete fBooNEHadronInelastic;
}

void BooNEHadronPhysics::ConstructParticle()
{

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle(); 
  
  return; 

}

#include "G4ProcessManager.hh"
void BooNEHadronPhysics::ConstructProcess()
{  
  G4ProcessManager* pManager = 0;
  auto aParticleIterator = GetParticleIterator();
  aParticleIterator->reset();
  while( (*aParticleIterator)() ) {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String pname=particle->GetParticleName();
    if(pname == "proton"   ||
       pname == "neutron"   ||
       pname == "pi+"   ||
       pname == "pi-" ||
       pname == "kaon+"   ||
       pname == "kaon-"   ||
       pname == "kaon0L"   ||
       pname == "kaon0S"   ) {
      pManager = particle->GetProcessManager();
      fBooNEHadronInelastic=new BooNEHadronInelasticProcess(*particle);
      pManager->AddDiscreteProcess(fBooNEHadronInelastic);
    }
  }
  
  //Anti barions
  fAntiBarionBuilder=new G4AntiBarionBuilder;
  fAntiBarionBuilder->RegisterMe(fFTFPAntiBarion=new G4FTFPAntiBarionBuilder(true)); //QuasiElastic = true/false
  fAntiBarionBuilder->Build();

  //Hyperons
  fHyperonBuilder=new G4HyperonFTFPBuilder; 
  fHyperonBuilder->Build();

  return;
}
