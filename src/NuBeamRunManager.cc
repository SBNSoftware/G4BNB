#include "NuBeamRunManager.hh"
#include "G4ios.hh"

NuBeamRunManager::NuBeamRunManager(NuBeamOutput *rec):
fGeometryInitializedHere(false),
fPhysicsInitializedHere(false),
fRecord(rec)
{
}

NuBeamRunManager::~NuBeamRunManager()
{
}

void NuBeamRunManager::InitializeGeometry() 
{
  if (fGeometryInitializedHere) {
    G4cout<<"NuBeamRunManager::InitializeGeometry() already called!"<<G4endl;
    return;
  }
  G4cout<<"NuBeamRunManager::Initialize geometry."<<G4endl;
  G4RunManager::InitializeGeometry();
  fGeometryInitializedHere = true;
}

void NuBeamRunManager::InitializePhysics() 
{
  if (fPhysicsInitializedHere) {
    G4cout<<"NuBeamRunManager::InitializePhysics() already called!"<<G4endl;
    return;
  }
  G4cout<<"NuBeamRunManager::Initialize physics."<<G4endl;
  G4RunManager::InitializePhysics();
  fPhysicsInitializedHere = true;
  G4cout<<"NuBeamRunManager::Physics initialized successfully."<<G4endl;
 }

