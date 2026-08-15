#ifndef NuBeamGeometryConstruction_h
#define NuBeamGeometryConstruction_h 1

#include <stdlib.h>
#include <cstdlib>
#include <cstdio>

#include "G4VUserDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include "G4PhysicalVolumeStore.hh"

class NuBeamGeometryConstructionMessenger;
class NuBeamLocalField;
class NuBeamSkinDepthField;
class G4Element;

class NuBeamGeometryConstruction : public G4VUserDetectorConstruction
{
public:
  NuBeamGeometryConstruction();
  ~NuBeamGeometryConstruction();
  G4VPhysicalVolume* Construct();

  void SetGeometryFile (G4String val) { 
    if( std::getenv("G4BNB_DIR") == 0 ) fGeometryFile = val;
    else fGeometryFile = G4String(G4String(std::getenv("G4BNB_DIR")) + G4String("/") + G4String(val.c_str()));
  };
  G4String GetGeometryFile() {return fGeometryFile;};
  void SetHornCurrent(G4double);

  inline NuBeamLocalField* GetLocalField() const
    { return pMyField; }
  inline void SetLocalField(NuBeamLocalField* aValue)
    { pMyField = aValue; }
  G4ThreeVector GetAbsolutePosition(G4String volString);

private:

  G4String _inFile;
  G4VPhysicalVolume* _pv;
  G4LogicalVolume*   _lv;

  G4LogicalVolume* nuspek_log;
  G4VPhysicalVolume* nuspek_phys;  

  G4String fGeometryFile;
  NuBeamGeometryConstructionMessenger* fMessenger;

  NuBeamLocalField*     pMyField;
  NuBeamSkinDepthField* pMySkinDepthField;

};

#endif
