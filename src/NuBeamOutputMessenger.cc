#include "NuBeamOutputMessenger.hh"
#include "NuBeamOutput.hh"

#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UIdirectory.hh"

#include <sstream>
using namespace std;

NuBeamOutputMessenger::NuBeamOutputMessenger(NuBeamOutput* O)
:fOutput(O)
{ 
  fDirectory = new G4UIdirectory("/boone/output/");
  fDirectory->SetGuidance( "BooNE output control commands.");

  fOutputFileName = new G4UIcmdWithAString("/boone/output/filename",this);
  fOutputFileName->SetGuidance("Set output ntuple file name");
  fOutputFileName->SetParameterName("filename",true);
  fOutputFileName->SetDefaultValue (G4String(""));
  fOutputFileName->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSaveProductionNtuple=new G4UIcmdWithABool("/boone/output/saveProductionNtuple",this);
  fSaveProductionNtuple->SetGuidance("Save production ntuple (true/false)");
  fSaveProductionNtuple->SetParameterName("saveProductionNtuple",true,false);
  fSaveProductionNtuple->SetDefaultValue(false);

  G4UIparameter* param;
  fBoundaryNtuple=new G4UIcommand("/boone/output/boundaryNtuple",this);
  fBoundaryNtuple->SetGuidance("Save particles at boundary of two volumes (exitVolume and enterVolume)");
  param=new G4UIparameter("ntupleName",'s',true);
  fBoundaryNtuple->SetParameter(param);
  param=new G4UIparameter("exitVolumeName",'s',true);
  fBoundaryNtuple->SetParameter(param);
  param=new G4UIparameter("enterVolumeName",'s',true);
  fBoundaryNtuple->SetParameter(param);

  fNuEnergyThrCmd = new G4UIcmdWithADoubleAndUnit("/boone/output/nuEnergyThr",this);
  fNuEnergyThrCmd->SetGuidance("Energy threshold for neutrinos in neutrino ntuple");
  fNuEnergyThrCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit (e.g. MeV, GeV, etc.)");
  fNuEnergyThrCmd->SetParameterName("nuEnergyThr", true, false);
  fNuEnergyThrCmd->SetDefaultUnit("MeV");
  fNuEnergyThrCmd->SetDefaultValue(0.);

  fNuEnergyUpperThrCmd = new G4UIcmdWithADoubleAndUnit("/boone/output/nuEnergyUpperThr",this);
  fNuEnergyUpperThrCmd->SetGuidance("Upper energy threshold for neutrinos in neutrino ntuple");
  fNuEnergyUpperThrCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit (e.g. MeV, GeV, etc.)");
  fNuEnergyUpperThrCmd->SetParameterName("nuEnergyUpperThr", true, false);
  fNuEnergyUpperThrCmd->SetDefaultUnit("MeV");
  fNuEnergyUpperThrCmd->SetDefaultValue(0.);
  
  fPionMomentumThrCmd = new G4UIcmdWithADoubleAndUnit("/boone/output/pionMomentumThr",this);
  fPionMomentumThrCmd->SetGuidance("Momentum threshold for parent pions in neutrino ntuple");
  fPionMomentumThrCmd->SetGuidance("To get neutrinos from pion decays at rest, choose any negative number");
  fPionMomentumThrCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit (e.g. MeV, GeV, etc.)");
  fPionMomentumThrCmd->SetParameterName("pionMomentumThr", true, false);
  fPionMomentumThrCmd->SetDefaultUnit("MeV");
  fPionMomentumThrCmd->SetDefaultValue(1.);
  
  fMuonMomentumThrCmd = new G4UIcmdWithADoubleAndUnit("/boone/output/muonMomentumThr",this);
  fMuonMomentumThrCmd->SetGuidance("Momentum threshold for parent muons in neutrino ntuple");
  fMuonMomentumThrCmd->SetGuidance("To get neutrinos from muon decays at rest, choose any negative number");
  fMuonMomentumThrCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit (e.g. MeV, GeV, etc.)");
  fMuonMomentumThrCmd->SetParameterName("muonMomentumThr", true, false);
  fMuonMomentumThrCmd->SetDefaultUnit("MeV");
  fMuonMomentumThrCmd->SetDefaultValue(1.);
  
  fKaonMomentumThrCmd = new G4UIcmdWithADoubleAndUnit("/boone/output/kaonMomentumThr",this);
  fKaonMomentumThrCmd->SetGuidance("Momentum threshold for parent kaons in neutrino ntuple");
  fKaonMomentumThrCmd->SetGuidance("To get neutrinos from kaon decays at rest, choose any negative number");
  fKaonMomentumThrCmd->SetGuidance("Available options: any positive floating-point number, with a valid unit (e.g. MeV, GeV, etc.)");
  fKaonMomentumThrCmd->SetParameterName("kaonMomentumThr", true, false);
  fKaonMomentumThrCmd->SetDefaultUnit("MeV");
  fKaonMomentumThrCmd->SetDefaultValue(0.);
}

NuBeamOutputMessenger::~NuBeamOutputMessenger()
{
  delete fOutputFileName;
  delete fNuEnergyThrCmd;
  delete fPionMomentumThrCmd;
  delete fMuonMomentumThrCmd;
  delete fKaonMomentumThrCmd;
  delete fDirectory;
}

void NuBeamOutputMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == fOutputFileName)
    fOutput->SetOutputFileName(newValue);

  if(command == fSaveProductionNtuple) 
    fOutput->SaveProductionNtuple(StoB(newValue));

  if(command == fBoundaryNtuple) {
    G4Tokenizer next(newValue);
    G4String ntpName=next();
    G4String exitVol=next();
    G4String enterVol=next();
    fOutput->AddBoundaryNtuple(ntpName, exitVol, enterVol);
  }

  if(command == fNuEnergyThrCmd)
    fOutput->SetNuEnergyThr(fNuEnergyThrCmd->GetNewDoubleValue(newValue));

  if(command == fNuEnergyUpperThrCmd)
    fOutput->SetNuEnergyUpperThr(fNuEnergyUpperThrCmd->GetNewDoubleValue(newValue));

  if(command == fPionMomentumThrCmd)
    fOutput->SetPionMomentumThr(fPionMomentumThrCmd->GetNewDoubleValue(newValue));

  if(command == fMuonMomentumThrCmd)
    fOutput->SetMuonMomentumThr(fMuonMomentumThrCmd->GetNewDoubleValue(newValue));

  if(command == fKaonMomentumThrCmd)
    fOutput->SetKaonMomentumThr(fKaonMomentumThrCmd->GetNewDoubleValue(newValue));
}

