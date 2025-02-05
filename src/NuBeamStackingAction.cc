#include "NuBeamStackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "NuBeamTrackInformation.hh"
#include "NuBeamTrajectory.hh"
#include "NuBeamRunManager.hh"

NuBeamStackingAction::NuBeamStackingAction(): fMessenger(0)
{ 
  fMessenger = new NuBeamStackingActionMessenger(this);
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/boone/stacking/lowKEThreshold");
}

NuBeamStackingAction::~NuBeamStackingAction()
{
  delete fMessenger;
  fMessenger = 0;
}

G4ClassificationOfNewTrack NuBeamStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;

  KillEMParticles(classification, aTrack);
  KillThresholdParticles(classification, aTrack);

  return classification;
}

void NuBeamStackingAction::NewStage() 
{
}

void NuBeamStackingAction::PrepareNewEvent()
{
}

void NuBeamStackingAction::KillEMParticles(G4ClassificationOfNewTrack &classification, 
					 const G4Track * aTrack)
{
   G4ParticleDefinition * particleType = aTrack->GetDefinition();
   bool isEM =   (particleType==G4Gamma::GammaDefinition()       ||
		  particleType==G4Electron::ElectronDefinition()  ||
		  particleType==G4Positron::PositronDefinition());

   if (isEM)
     classification=fKill;
}

void NuBeamStackingAction::KillThresholdParticles(G4ClassificationOfNewTrack& classification, 
						  const G4Track * aTrack)
{
  G4ParticleDefinition * particleType = aTrack->GetDefinition();
   //Discard particles with kinetic energy < KillTracking Threshold.GeV (that are not neutrinos, muons, pions, or kaons )
  if ((std::abs(particleType->GetPDGEncoding()) != G4MuonPlus::MuonPlus()->GetPDGEncoding()) &&
      (std::abs(particleType->GetPDGEncoding()) != G4PionPlus::PionPlus()->GetPDGEncoding()) &&
      (std::abs(particleType->GetPDGEncoding()) != G4KaonPlus::KaonPlus()->GetPDGEncoding()) &&
      (particleType!=G4Proton::Proton())&&
      (particleType!=G4KaonZeroLong::KaonZeroLong())&&
      (particleType!=G4KaonZero::KaonZero())&&
      (particleType!=G4NeutrinoE::NeutrinoE())&&
      (particleType!=G4NeutrinoMu::NeutrinoMu())&&
      (particleType!=G4NeutrinoTau::NeutrinoTau())&&
      (particleType!=G4AntiNeutrinoE::AntiNeutrinoE())&&
      (particleType!=G4AntiNeutrinoMu::AntiNeutrinoMu())&&
      (particleType!=G4AntiNeutrinoTau::AntiNeutrinoTau())) {
    G4double energy = aTrack->GetKineticEnergy();
    
    if ((energy < fLowKEThreshold * CLHEP::MeV/CLHEP::GeV ) && (classification != fKill))
      {classification = fKill;}
  }  
}

