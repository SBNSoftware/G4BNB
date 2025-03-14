#include "G4ios.hh"
#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Polycone.hh"

#include "NuBeamSkinDepthFieldMessenger.hh"
#include "NuBeamSkinDepthField.hh"

NuBeamSkinDepthFieldMessenger::NuBeamSkinDepthFieldMessenger(NuBeamSkinDepthField* SDF)
  :fSkinDepthField(SDF)
{ 
  fDirectory = new G4UIdirectory("/boone/field/skin/");
  fDirectory->SetGuidance( "BooNE magnetic field commands.");

  // Set Skin depth model
  fSkinDepthModelCmd = new G4UIcmdWithAString("/boone/field/skin/SkinDepthModel",this);
  fSkinDepthModelCmd->SetGuidance("Model for Skin Depth Field");
  fSkinDepthModelCmd->SetGuidance("Description:");
  fSkinDepthModelCmd->SetGuidance("NONE:    Current entirely on outer radius of inner conductor");
  fSkinDepthModelCmd->SetGuidance("FLAT:    Current uniformly distrbuted in inner conductor");
  fSkinDepthModelCmd->SetGuidance("IN:      Current entirely on inner radius of inner conductor");
  fSkinDepthModelCmd->SetGuidance("2D:      Two-dimensional field penetration from Jackson");
  fSkinDepthModelCmd->SetGuidance("EXPI:    Current exponential damped from outer surface of inner conductor");
  fSkinDepthModelCmd->SetParameterName("SkinDepthModel", true, false);
  fSkinDepthModelCmd->SetDefaultValue(fSkinDepthField->GetSkinDepthModel());
  fSkinDepthModelCmd->SetCandidates("NONE FLAT IN 2D EXPI");

  //Set the Skin depth 
  fSkinDepthCmd = new G4UIcmdWithADoubleAndUnit("/boone/field/skin/SkinDepth",this);
  fSkinDepthCmd->SetGuidance("Skin depth");
  fSkinDepthCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit for the length (e.g. mm, cm, etc.)");
  fSkinDepthCmd->SetGuidance("Will be used only when skinDepthModel is set to 2d");
  fSkinDepthCmd->SetParameterName("SkinDepth", true, false);
  fSkinDepthCmd->SetDefaultUnit("mm");
  fSkinDepthCmd->SetDefaultValue(fSkinDepthField->GetSkinDepth());

  fSkinDepthHornCurrentCmd = new G4UIcmdWithADoubleAndUnit("/boone/field/skin/SkinDepthHornCurrent",this);
  fSkinDepthHornCurrentCmd->SetGuidance("Skin horn current");
  fSkinDepthHornCurrentCmd->SetGuidance("Available options: any positive or negative floating-point number, with a valid unit for the current (e.g. kA, A, etc.)");
  fSkinDepthHornCurrentCmd->SetParameterName("hornCur", true, false);
  fSkinDepthHornCurrentCmd->SetDefaultUnit("ampere");
  fSkinDepthHornCurrentCmd->SetDefaultValue(fSkinDepthField->GetSkinDepthHornCurrent()/CLHEP::ampere);

}

NuBeamSkinDepthFieldMessenger::~NuBeamSkinDepthFieldMessenger()
{

  delete fSkinDepthModelCmd;
  delete fSkinDepthHornCurrentCmd;
  delete fSkinDepthCmd;

  delete fDirectory;

}

void NuBeamSkinDepthFieldMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 

  if(command == fSkinDepthModelCmd){
    fSkinDepthField->SetSkinDepthModel(newValue);
  }

  if(command == fSkinDepthCmd){
    fSkinDepthField->SetSkinDepth(fSkinDepthCmd->GetNewDoubleValue(newValue));
  }
  
  if(command == fSkinDepthHornCurrentCmd){
    fSkinDepthField->SetSkinDepthHornCurrent(fSkinDepthHornCurrentCmd->GetNewDoubleValue(newValue));
  }

}

