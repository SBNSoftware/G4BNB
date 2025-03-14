#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4UImanager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Polycone.hh"
#include "G4UnitsTable.hh" 

#include "NuBeamSkinDepthField.hh"
#include "NuBeamSteppingAction.hh"

NuBeamSkinDepthField::NuBeamSkinDepthField()
 :fMessenger(0),
 fSkinDepthHornCurrent(174000.*CLHEP::ampere), 
 fSkinDepth(1.4*CLHEP::mm),
 fSkinDepthModel("EXPI")
{ 
  fMessenger = new NuBeamSkinDepthFieldMessenger(this);
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/boone/field/skin/SkinDepthModel");
  UI->ApplyCommand("/boone/field/skin/SkinDepth");
  UI->ApplyCommand("/boone/field/skin/SkinDepthHornCurrent");

  fEquation = new G4Mag_UsualEqRhs(this); 
  fStepperType = 4 ;      // ClassicalRK4 is default stepper 

  G4cout << "Constructing Skin Depth Model" << G4endl;

  mu0 = 4.*M_PI*0.0000001;
  
  // Find the pointers to the HWL1 and HRN1 physical volumes
  G4PhysicalVolumeStore* pPVStore = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* HWL1Vol = NULL;
  G4VPhysicalVolume* HRN1Vol = NULL;
  for( size_t i = 0; i < pPVStore->size(); i++ ) {
    G4VPhysicalVolume* pv = (*pPVStore)[i];
    if ( pv->GetName() == "HWL1" ) {
      HWL1Vol = pv;
    } else if ( pv->GetName() == "HRN1" ) {
      HRN1Vol = pv;
    }
  }

  // Get the polycone details of the inner conductor from the geometry specification
  G4Polycone* HRN1Poly = (G4Polycone*)HRN1Vol->GetLogicalVolume()->GetSolid();
  Npoints = HRN1Poly->GetOriginalParameters()->Num_z_planes;
  G4double *zoutertemp =  HRN1Poly->GetOriginalParameters()->Z_values;
  G4double *routertemp =  HRN1Poly->GetOriginalParameters()->Rmin;

  G4Polycone* HWL1Poly = (G4Polycone*)HWL1Vol->GetLogicalVolume()->GetSolid();
  G4double *zinnertemp =  HWL1Poly->GetOriginalParameters()->Z_values;
  G4double *rinnertemp =  HWL1Poly->GetOriginalParameters()->Rmin;

  for(int i = 0; i < Npoints; i++){
    innerZ[i] = zinnertemp[i]/CLHEP::cm;
    innerR[i] = rinnertemp[i]/CLHEP::cm;
    outerZ[i] = zoutertemp[i]/CLHEP::cm;
    outerR[i] = routertemp[i]/CLHEP::cm;
  }

  G4cout << Npoints << " points in the Polycone" << G4endl;

  G4cout << "outerZ = ";
  for(int i = 0; i < Npoints; i++){ G4cout << outerZ[i] << "  "; }
  G4cout << G4endl;
  G4cout << "outerR = ";
  for(int i = 0; i < Npoints; i++){ G4cout << outerR[i] << "  "; }
  G4cout << G4endl;
  G4cout << "innerZ = ";
  for(int i = 0; i < Npoints; i++){ G4cout << innerZ[i] << "  "; }
  G4cout << G4endl;
  G4cout << "innerR = ";
  for(int i = 0; i < Npoints; i++){ G4cout << innerR[i] << "  "; }
  G4cout << G4endl;

}


NuBeamSkinDepthField::~NuBeamSkinDepthField()
{
  if (fStepper) delete fStepper;
}

void NuBeamSkinDepthField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  // Check to see if we are in the appropriate  z range. Otherwise return zero field
  if(Point[2]/CLHEP::cm < outerZ[0] || Point[2]/CLHEP::cm > outerZ[Npoints-1]){
    Bfield[0] = 0.;
    Bfield[1] = 0.;
    Bfield[2] = 0.;
    return;
  }

  // If the Skin Depth  Model is "NONE", we want zero field in all cases.
  if(this->GetSkinDepthModel() == "NONE"){
    Bfield[0] = 0.;
    Bfield[1] = 0.;
    Bfield[2] = 0.;
    return;
  }
  
  // Now we will deal with real models: some preliminary geometry
  // we are using MKS units; all lengths in meters, currents in amperes, fields in tesla
  // find how far in to the alumininum we are
  G4double r   = sqrt(Point[0]*Point[0] + Point[1]*Point[1])/CLHEP::m;
  // okay, here's the exception: we pulled info from the Geometry file which is in cm
  // we need to do this part in cm, but then rinner/router will be in meters.
  G4double zcm  = Point[2]/CLHEP::cm; 
  G4double rinnercm = GetRInnercm(zcm);
  G4double routercm = GetROutercm(zcm);
  double rinner   = rinnercm/100.0;
  double router   = routercm/100.0;
  
  // eject with zero field if we are outside the appropriate transverse radius.
  if(r < rinner || r > router){
    Bfield[0] = 0.;
    Bfield[1] = 0.;
    Bfield[2] = 0.;
    return;
  }
  
  // Calculate Magnitude of Field in non-trivial cases
  G4double B = 0;
  if(this->GetSkinDepthModel() == "IN"){
    //The skin depth model is IN, all current lies on inner surface, field in conductor is 1/r
    B = ( (mu0/(2.* M_PI) ) * (fSkinDepthHornCurrent/CLHEP::ampere)/r ) * CLHEP::tesla;
  } else if(this->GetSkinDepthModel() == "FLAT"){
    // The skin depth model is FLAT, current distribution is flat within conductor: use Ampere's law
    //G4double IFlat = (fSkinDepthHornCurrent/CLHEP::ampere)*(r-rinner)/(router-rinner);
    //my flat
    G4double IFlat = (fSkinDepthHornCurrent/CLHEP::ampere)*(r*r-rinner*rinner)/(router*router-rinner*rinner);
    B = ( (mu0/(2.* M_PI) ) * IFlat/r ) * CLHEP::tesla;
  } else if(this->GetSkinDepthModel() == "2D"){
    //  2 Dimensional field penetration model from Jackson 
    G4double B0 = ( (mu0/(2.* M_PI) ) * (fSkinDepthHornCurrent/CLHEP::ampere)/router ) * CLHEP::tesla;
    B  = B0 * exp(-(router-r) / (fSkinDepth/CLHEP::m) ) * cos( (router-r) / (fSkinDepth/CLHEP::m) );
  } else if(this->GetSkinDepthModel() == "EXPI"){
    //  Exponentially damped current from the outer surface
    G4double I0 = (fSkinDepthHornCurrent/CLHEP::ampere) / 
      ( (fSkinDepth/CLHEP::m) * (1.0 - exp( -(router-rinner)/ (fSkinDepth/CLHEP::m)) ) );
    G4double Ienclosed = I0 * (fSkinDepth/CLHEP::m) * 
      ( exp( -(router-r)/(fSkinDepth/CLHEP::m) ) - exp( -(router-rinner)/(fSkinDepth/CLHEP::m) ) );
    B = ( (mu0/(2.* M_PI) ) * (Ienclosed)/r ) * CLHEP::tesla;
  }
  
  Bfield[0] = - B * ( (Point[1]/CLHEP::m)/r);
  Bfield[1] =   B * ( (Point[0]/CLHEP::m)/r);
  Bfield[2] =   0.;
  /*
  G4cout<<"Skin depth "
	    <<Point[0]/CLHEP::mm<<"\t"<<Point[1]/CLHEP::mm<<"\t"
	<<sqrt(Point[0]/CLHEP::mm*Point[0]/CLHEP::mm+Point[1]/CLHEP::mm*Point[1]/CLHEP::mm)<<"\t"
	    <<Bfield[0]<<"\t"<<Bfield[1]<<G4endl;
  */
}


G4double NuBeamSkinDepthField::GetRInnercm(const G4double zcm) const
{
  // z here is in cm. Radius also in cm;
  if(zcm < innerZ[0] || zcm > innerZ[Npoints-1]){
     G4Exception("NuBeamSkinDepthField", "111", JustWarning, 
		 "NuBeamSkinDepthField: z is out range"); 
     return 0;
  }

  G4double rcm = 0;
  for(int i = 0; i < Npoints; i++){
    if(zcm > innerZ[i] && innerZ[i+1]){
      rcm = innerR[i] + (innerR[i+1] - innerR[i]) * (zcm-innerZ[i])/(innerZ[i+1]-innerZ[i]);
    }
  }

  if(rcm == 0){
    G4Exception("NuBeamSkinDepthField", "111", JustWarning, 
		"NuBeamSkinDepthField: z is out range, inner radius undetermined"); 
  }

  return rcm;
}

G4double NuBeamSkinDepthField::GetROutercm(const G4double zcm) const
{
  // z here is in cm. Radius also in cm;
  if(zcm < outerZ[0] || zcm > outerZ[Npoints-1]){
     G4Exception("NuBeamSkinDepthField", "111", JustWarning, 
		 "NuBeamSkinDepthField: z is out range"); 
     return 0;
  }

  G4double rcm = 0;
  for(int i = 0; i < Npoints; i++){
    if(zcm > outerZ[i] && outerZ[i+1]){
      rcm = outerR[i] + (outerR[i+1] - outerR[i]) * (zcm-outerZ[i])/(outerZ[i+1]-outerZ[i]);
    }
  }

  if(rcm == 0){
    G4Exception("NuBeamSkinDepthField", "111", JustWarning, 
		"NuBeamSkinDepthField: z is out range, outer radius undetermined"); 
  }

  return rcm;
}
