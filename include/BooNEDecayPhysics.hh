#ifndef BooNEDecayPhysics_h
#define BooNEDecayPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4Decay.hh"
#include "BooNECustomDecay.hh"

class BooNEDecayPhysics : public G4VPhysicsConstructor
{
  public: 
    BooNEDecayPhysics(G4int ver = 1);
    BooNEDecayPhysics(const G4String& name, G4int ver = 1);
    virtual ~BooNEDecayPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
  virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
  virtual void ConstructProcess();

  void CorrectDecayTables();

private:
  G4int    verbose;
  G4bool   wasActivated;
};

#endif
