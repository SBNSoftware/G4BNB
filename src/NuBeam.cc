// header files includes
#include "NuBeamRunManager.hh"
#include "NuBeamGeometryConstruction.hh"
#include "NuBeamPrimaryGeneratorAction.hh"
#include "NuBeamRunAction.hh"
#include "NuBeamEventAction.hh"
#include "NuBeamStackingAction.hh"
#include "NuBeamTrackingAction.hh"
#include "NuBeamSteppingAction.hh"
#include "NuBeamOutput.hh"
#include "NuBeamPhysicsList.hh"

// geant4 includes
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc,char** argv)
{

  G4cout << "Running NuBeam." << G4endl;

  G4String macroFile = "";
  if (argc >= 2) {
    macroFile = argv[1];
    std::ifstream mac(macroFile);
    if (!mac) {
      G4cout << "Cannot open macro file """ << macroFile << """" << G4endl;
      return 2;
    }
  }
  if (argc >= 3) {
    G4cerr << "Too many command line arguments (" << argc <<")" << G4endl;
    return EXIT_FAILURE;
  }
  
  NuBeamOutput* myRecords = new NuBeamOutput();
  NuBeamRunManager *runManager = new NuBeamRunManager(myRecords); 
  // set mandatory initialization classes
  runManager->SetUserInitialization(new NuBeamGeometryConstruction());
  runManager->SetUserInitialization(new NuBeamPhysicsList());
  runManager->SetUserAction(new NuBeamRunAction());
  runManager->SetUserAction(new NuBeamPrimaryGeneratorAction());
  runManager->SetUserAction(new NuBeamEventAction());
  runManager->SetUserAction(new NuBeamStackingAction());
  runManager->SetUserAction(new NuBeamTrackingAction());
  runManager->SetUserAction(new NuBeamSteppingAction());
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif
  
  G4bool batch_mode = macroFile != "";
  
  G4UImanager * UI = G4UImanager::GetUIpointer();    
  if(!batch_mode) {
    G4UIsession * session = new G4UIterminal;
    if (session != 0) {
      session->SessionStart();
      delete session;
    }
  } else {
    // Batch mode
    G4String command = "/control/execute ";
    UI->ApplyCommand(command+macroFile);
  }
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  if (myRecords) delete myRecords;
  if (runManager) delete runManager;
  return EXIT_SUCCESS;
}
