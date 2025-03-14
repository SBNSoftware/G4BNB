#ifndef BooNEHadronInelasticProcess_h
#define BooNEHadronInelasticProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"
#include "BooNEpXInteraction.hh"

class BooNEpBeInteraction;
class BooNEpXInteraction;
class G4InelasticInteraction;

class XSecArrayHolder;

class BooNEHadronInelasticProcess : public G4HadronicProcess
{
public:
  BooNEHadronInelasticProcess(const G4ParticleDefinition& aParticleType);
  ~BooNEHadronInelasticProcess();
  
  //copy PostStepDoIt from G4HadronicProcess and change to call ChooseHadronicInteraction from this class (not base)
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step&) override;  
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;

protected:

  //reimplement ChooseHadronicInteraction to allow BooNEpBe model to completely
  //overlap in energy with otherwise default model and choose BooNEpBe model
  //exclusively in applicable region
  G4HadronicInteraction* ChooseHadronicInteraction(G4HadProjectile&, G4Nucleus&, G4Material*&, G4Element*&);
  
  // access to the chosen generator
  BooNEpXInteraction* GetHadronicInteraction() const
  { return fInteraction; }

private:

  const G4ParticleDefinition* fParticle;
  BooNEpBeInteraction*        fBooNEpBeModel;
  G4HadronicInteraction*      fLEProtonModel;
  G4HadronicInteraction*      fHEProtonModel;
  // Wrappers around the above models with abstracted ApplyYourself()
  BooNEpXInteraction*         fBooNEModel_pBe;
  BooNEpXInteraction*         fBooNEModel_LEp;
  BooNEpXInteraction*         fBooNEModel_HEp;
  //G4HadronicInteraction*      fInteraction;
  BooNEpXInteraction*         fInteraction;

};
#endif
