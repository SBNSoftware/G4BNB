#ifndef BooNEHadronElasticProcess_h
#define BooNEHadronElasticProcess_h 1
 
#include "globals.hh"
#include "G4HadronElasticProcess.hh"
#include "G4CrossSectionDataStore.hh"
//#include "G4HadronElasticDataSet.hh"
#include "G4ChipsHyperonElasticXS.hh"
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonZeroElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsPionMinusElasticXS.hh"
#include "G4ChipsPionPlusElasticXS.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4VDiscreteProcess.hh"
//#include "G4LightMedia.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"


class BooNEHadronElasticProcess : public G4HadronElasticProcess
{
public:

  BooNEHadronElasticProcess(const G4String& processName = "BooNEHadronElastic");
  ~BooNEHadronElasticProcess();
  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  /*
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
  
  void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);
  */
private:
  /*  
  G4double GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
				      const G4Element* anElement,
				      const G4Material* aMaterial);
  */
};
#endif
