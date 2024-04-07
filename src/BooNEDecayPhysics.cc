#include "BooNEDecayPhysics.hh"

#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4KL3DecayChannel.hh"

BooNEDecayPhysics::BooNEDecayPhysics(G4int ver)
  :  G4VPhysicsConstructor("BooNEDecay"), verbose(ver), wasActivated(false)
{
}

BooNEDecayPhysics::BooNEDecayPhysics(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
}

BooNEDecayPhysics::~BooNEDecayPhysics()
{
}

void BooNEDecayPhysics::ConstructParticle()
{

  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

  this->CorrectDecayTables();
}

void BooNEDecayPhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // Add Decay Process
  G4Decay* decayProcess = new G4Decay();
  auto aParticleIterator = GetParticleIterator();
  aParticleIterator->reset();
  G4ParticleDefinition* particle=0;

  while( (*aParticleIterator)() )
  {
    particle = aParticleIterator->value();
    if( decayProcess->IsApplicable(*particle) && particle->GetParticleName()!="neutron" ) 
    { 
      if(verbose > 1) {
        G4cout << "### Decays for " << particle->GetParticleName() << G4endl;
      }
      ph->RegisterProcess(decayProcess, particle);
    }
  }
}

void BooNEDecayPhysics::CorrectDecayTables()
{
  G4cout<<"BooNEDecayPhysics::Correcting decay tables for pi+-, K+-, K0L"<<G4endl;

  // add pi->e nue to decay table  
  G4DecayTable* pPionPlusDecayTable = G4PionPlus::PionPlus()->GetDecayTable();
  G4VDecayChannel* munuMode;
  G4VDecayChannel* enuMode;
  G4double enuModeBR = 1.230e-4;
  munuMode = pPionPlusDecayTable->GetDecayChannel(0);
  munuMode->SetBR(1.-enuModeBR);
  enuMode = new G4PhaseSpaceDecayChannel("pi+",enuModeBR,2,"e+","nu_e");
  pPionPlusDecayTable->Insert(enuMode);
  pPionPlusDecayTable->DumpInfo();

  // add pi->e nue to decay table
  G4DecayTable* pPionMinusDecayTable = G4PionMinus::PionMinus()->GetDecayTable();
  munuMode = pPionMinusDecayTable->GetDecayChannel(0);
  munuMode->SetBR(1.-enuModeBR);
  enuMode = new G4PhaseSpaceDecayChannel("pi-",enuModeBR,2,"e-","anti_nu_e");
  pPionMinusDecayTable->Insert(enuMode);
  pPionMinusDecayTable->DumpInfo();

  // Redefinition of KaonPlus branching ratios (numbers from PDG 2016).
  G4DecayTable* pKaonPlusDecayTable = new G4DecayTable();
  // create decay channels
  G4int nmodes = 6;
  G4VDecayChannel** mode = new G4VDecayChannel*[nmodes];
  G4double kp_mup_numu = 0.6356;     // 0.6356(11)
  G4double kp_pip_pi0 = 0.2067;      // 0.2067(8)
  G4double kp_pip_pip_pim = 0.05583; // 0.05583(24)
  G4double kp_pip_pi0_pi0 = 0;       // set value below
  G4double kp_pi0_ep_nue = 0.0507;   // 0.0507(4)
  G4double kp_pi0_mup_numu = 0.03352; // 0.03352(33)
    // put the remainder in K+ -> pi+ pi0 pi0 to make the sum of all modes unity
  G4double brsum = kp_mup_numu + kp_pip_pi0 + kp_pip_pip_pim
    + kp_pi0_ep_nue + kp_pi0_mup_numu;
  if (brsum < 1.) {
    kp_pip_pi0_pi0 = 1. - brsum; // actual value = 0.01757(24)
  } else {
    kp_pip_pi0_pi0 = 0.00001;
    G4cout << "Error in charged K branching ratios:  sum of BRs > 1   (sum = "
	   << brsum << ")" << G4endl;
  }

  // kaon+ -> mu+ + nu_mu
  mode[0] = new G4PhaseSpaceDecayChannel("kaon+",kp_mup_numu,2,"mu+","nu_mu");
  // kaon+ -> pi+ + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("kaon+",kp_pip_pi0,2,"pi+","pi0");
  // kaon+ -> pi+ + pi+ + pi-
  mode[2] = new G4PhaseSpaceDecayChannel("kaon+",kp_pip_pip_pim,3,"pi+","pi+","pi-");
  // kaon+ -> pi+ + pi0 + pi0
  mode[3] = new G4PhaseSpaceDecayChannel("kaon+",kp_pip_pi0_pi0,3,"pi+","pi0","pi0");
  // kaon+ -> pi0 + e+ + nu_e (Ke3)
  mode[4] = new G4KL3DecayChannel("kaon+",kp_pi0_ep_nue,"pi0","e+","nu_e");
  // kaon+ -> pi0 + mu+ + nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon+",kp_pi0_mup_numu,"pi0","mu+","nu_mu");

  // load decay modes into the decay table
  for (G4int index=0; index<nmodes; index++ ) {
    pKaonPlusDecayTable->Insert(mode[index]);
  }
  delete [] mode;

  // delete the old decay table and replace with our new table
  G4DecayTable* oldtable = G4KaonPlus::KaonPlus()->GetDecayTable(); 
  delete oldtable;
  G4KaonPlus::KaonPlus()->SetDecayTable(pKaonPlusDecayTable);
  pKaonPlusDecayTable->DumpInfo();

  // Redefinition of KaonMinus branching ratios (numbers from PDG 2006).
  G4DecayTable* pKaonMinusDecayTable = new G4DecayTable();

  // create decay channels with same branching ratios as the K+ channels
  nmodes = 6;
  mode = new G4VDecayChannel*[nmodes];
  G4double km_mum_numubar = kp_mup_numu;
  G4double km_pim_pi0 = kp_pip_pi0;
  G4double km_pim_pip_pim = kp_pip_pip_pim;
  G4double km_pim_pi0_pi0 = kp_pip_pi0_pi0;
  G4double km_pi0_em_nuebar = kp_pi0_ep_nue;
  G4double km_pi0_mum_numubar = kp_pi0_mup_numu;

  // kaon- -> mu- + anti_nu_mu
  mode[0] = new G4PhaseSpaceDecayChannel("kaon-",km_mum_numubar,2,"mu-","anti_nu_mu");
  // kaon- -> pi- + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("kaon-",km_pim_pi0,2,"pi-","pi0");
  // kaon- -> pi- + pi+ + pi-
  mode[2] = new G4PhaseSpaceDecayChannel("kaon-",km_pim_pip_pim,3,"pi-","pi+","pi-");
  // kaon- -> pi- + pi0 + pi0
  mode[3] = new G4PhaseSpaceDecayChannel("kaon-",km_pim_pi0_pi0,3,"pi-","pi0","pi0");
  // kaon- -> pi0 + e- + anti_nu_e (Ke3)
  mode[4] = new G4KL3DecayChannel("kaon-",km_pi0_em_nuebar,"pi0","e-","anti_nu_e");
  // kaon- -> pi0 + mu- + anti_nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon-",km_pi0_mum_numubar,"pi0","mu-","anti_nu_mu");

  for (G4int index=0; index<nmodes; index++ ) {
    pKaonMinusDecayTable->Insert(mode[index]);
  }
  delete [] mode;

  // delete the old decay table and replace with our new table
  oldtable = G4KaonMinus::KaonMinus()->GetDecayTable(); 
  delete oldtable;
  G4KaonMinus::KaonMinus()->SetDecayTable(pKaonMinusDecayTable);
  pKaonMinusDecayTable->DumpInfo();

  // Redefinition of KaonZeroLong branching ratios (numbers from PDG 2016).
  G4DecayTable* pKaonZeroLDecayTable = new G4DecayTable();

  // create decay channels
  nmodes = 9;
  mode = new G4VDecayChannel*[nmodes];
  G4double kl_pi0_pi0_pi0 = 0.1952;      // 0.1952(12)
  G4double kl_pip_pim_pi0 = 0.1254;      // 0.1254(5)
  G4double kl_pim_ep_nue = 0.20275;      // 0.4055(11)
  G4double kl_pip_em_nuebar = 0.20275;   // 
  G4double kl_pim_mup_numu = 0.1352;    // 0.2704(7)
  G4double kl_pip_mum_numubar = 0.1352; // 
  G4double kl_pip_pim = 0.001967;        // 0.001967(10)
  G4double kl_pi0_pi0 = 0.000864;        // 0.000864(6)
  G4double kl_gamma_gamma;               // set value below

  // put the remainder in KL -> gg to make the sum of all modes unity
  brsum = kl_pi0_pi0_pi0 + kl_pip_pim_pi0 + kl_pim_ep_nue + kl_pip_em_nuebar
    + kl_pim_mup_numu + kl_pip_mum_numubar + kl_pip_pim + kl_pi0_pi0;
  if (brsum < 1.) {
    kl_gamma_gamma = 1. - brsum;  // actual value = 0.000548(5)
  } else {
    kl_gamma_gamma = 0.00001;
    G4cout << "Error in KL branching ratios:  sum of BRs > 1   (sum = "
	   << brsum << ")" << G4endl;
  }

  // kaon0L -> pi0 + pi0 + pi0
  mode[0] = new G4PhaseSpaceDecayChannel("kaon0L",kl_pi0_pi0_pi0,3,"pi0","pi0","pi0");
  // kaon0L -> pi0 + pi+ + pi-
  mode[1] = new G4PhaseSpaceDecayChannel("kaon0L",kl_pip_pim_pi0,3,"pi0","pi+","pi-");
  // kaon0L -> pi- + e+ + nu_e (Ke3)
  mode[2] = new G4KL3DecayChannel("kaon0L",kl_pim_ep_nue,"pi-","e+","nu_e");
  // kaon0L -> pi+ + e- + anti_nu_e (Ke3)
  mode[3] = new G4KL3DecayChannel("kaon0L",kl_pip_em_nuebar,"pi+","e-","anti_nu_e");
  // kaon0L -> pi- + mu+ + nu_mu (Kmu3)
  mode[4] = new G4KL3DecayChannel("kaon0L",kl_pim_mup_numu,"pi-","mu+","nu_mu");
  // kaon0L -> pi+ + mu- + anti_nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon0L",kl_pip_mum_numubar,"pi+","mu-","anti_nu_mu");
  // kaon0L -> pi+ + pi-
  mode[6] = new G4PhaseSpaceDecayChannel("kaon0L",kl_pip_pim,2,"pi+","pi-");
  // kaon0L -> pi0 + pi0
  mode[7] = new G4PhaseSpaceDecayChannel("kaon0L",kl_pi0_pi0,2,"pi0","pi0");
  // kaon0L -> gamma + gamma
  mode[8] = new G4PhaseSpaceDecayChannel("kaon0L",kl_gamma_gamma,2,"gamma","gamma");

  // load decay modes into the decay table
  for (G4int index=0; index<nmodes; index++ ) {
    pKaonZeroLDecayTable->Insert(mode[index]);
  }
  delete [] mode;

  // delete the old decay table and replace with our new table
  oldtable = G4KaonZeroLong::KaonZeroLong()->GetDecayTable(); 
  delete oldtable;
  G4KaonZeroLong::KaonZeroLong()->SetDecayTable(pKaonZeroLDecayTable);
  pKaonZeroLDecayTable->DumpInfo();

}
