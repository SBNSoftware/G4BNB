#ifndef NuBeamOutput_h
#define NuBeamOutput_h 1 
#include <iostream>
#include "globals.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4HadFinalState.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4DecayProducts.hh"

#include "NuBeamTrajectory.hh"
#include "NuBeamLocalField.hh"
#include "NuBeamSkinDepthField.hh"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

//ROOT includes
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

const G4int max_ancestors = 20;
const G4int max_multiplicity = 200;

class NuBeamOutputMessenger;
struct Ntuple101_def {
  float beamwgt;
  int ntp;
  int npart;
  int id[max_ancestors];
  float ini_pos[max_ancestors][3];
  float ini_mom[max_ancestors][3];
  float ini_eng[max_ancestors];
  float ini_t[max_ancestors];
  float fin_mom[max_ancestors][3];
  float fin_pol[max_ancestors][3];
};
struct Ntuple102_def {
  int nmult;
  int id[max_multiplicity];
  float mom[max_multiplicity][3];
  float wgt[max_multiplicity];
};
struct Ntuple103_def {
  int pdg;
  float beamwgt;
  float pos[3];
  float mom[3];
  float eng;
};

struct boundaryNtp_t {
  G4String fName;
  G4String fExitVol;
  G4String fEnterVol;
  Ntuple103_def fNtp;
  TTree* fTree;
};

class NuBeamOutput 
{
public:
  NuBeamOutput();
  virtual ~NuBeamOutput();

  void RecordBeginOfRun(const G4Run*);
  void RecordEndOfRun(const G4Run*);
  void RecordBeginOfEvent(const G4Event*);
  void RecordEndOfEvent(const G4Event*);
  void RecordBeginOfTrack(const G4Track*);
  void RecordEndOfTrack(const G4Track*);
  void RecordStep(const G4Step*);
  void RecordNeutrino(const G4Track*);

  void RecordpBeInteraction(G4HadFinalState* aParticleChange);
  void RecordNuParentDecay(const G4Track* aTrack);
    
  void SetOutputFileName (G4String val) { fOutputFileName = val;};
  G4String GetOutputFileName() { return fOutputFileName;}

  void SaveProductionNtuple(bool val) {fSaveProductionNtuple = val;};
  void SaveRestFrameDecay(bool val) {fSaveRestFrameDecay = val;}

  void SetNuEnergyThr(G4double val){ fNuEnergyThr = val;};
  G4double GetNuEnergyThr() const { return fNuEnergyThr;};
  void SetPionMomentumThr(G4double val){ fPionMomentumThr = val;};
  G4double GetPionMomentumThr() const { return fPionMomentumThr;};
  void SetMuonMomentumThr(G4double val){ fMuonMomentumThr = val;};
  G4double GetMuonMomentumThr() const { return fMuonMomentumThr;};
  void SetKaonMomentumThr(G4double val){ fKaonMomentumThr = val;};
  G4double GetKaonMomentumThr() const { return fKaonMomentumThr;};
  
  // gets the muon polarization from the pion and muon momenta 
  virtual G4ThreeVector GetMuonPolarization(
					    G4double chargeG4,
					    G4double epiG4,
					    G4ThreeVector &ppiG4,
					    G4double emuG4,
					    G4ThreeVector &pmuG4);

  void AddBoundaryNtuple(G4String ntpName, G4String exitVol, G4String enterVol) {
    boundaryNtp_t ntp;
    ntp.fName=ntpName;
    ntp.fExitVol=exitVol;
    ntp.fEnterVol=enterVol;
    fBoundaryNtp.push_back(ntp);
  }

private:
  void fillDkMeta();
  NuBeamTrajectory* GetTrajectory(G4int parentID);

  G4int GetDecayCode(const G4Track* nuTrack);

  G4String fGeometryFile;
  G4String fOutputFileName;

  NuBeamOutputMessenger* fMessenger;

  // trajectory objects, used when tracing back in neutrinos' ancestor tree
  NuBeamTrajectory* trj;
  NuBeamTrajectory* parentTrj;
  NuBeamTrajectory* grandParentTrj;
  NuBeamTrajectory* ancestorTrj;

  TFile*        fOutFileDk2Nu; 
  TTree*        fOutTreeDk2Nu;
  TTree*        fOutTreeDk2NuMeta;
  bsim::Dk2Nu*  fDk2Nu;
  bsim::DkMeta* fDkMeta;

  bool   fSaveProductionNtuple;
  TTree* fProductionTree;

  bool   fSaveRestFrameDecay;
  TTree* fDARTree;

  G4double fNuEnergyThr;
  G4double fPionMomentumThr;
  G4double fMuonMomentumThr;
  G4double fKaonMomentumThr;

  // ntuple structs
  Ntuple101_def fNtuple101;
  Ntuple102_def fNtuple102;

  //aux ntuples
  std::vector<boundaryNtp_t> fBoundaryNtp;
};

#endif
