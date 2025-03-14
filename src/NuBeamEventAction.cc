#include "NuBeamOutput.hh"
#include "NuBeamRunManager.hh"
#include "NuBeamEventAction.hh"
#include "NuBeamEventActionMessenger.hh"

#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

NuBeamEventAction::NuBeamEventAction()
  :fPrintModulo(10),fMessenger(0)
{
  const NuBeamRunManager *pRunManager=
    reinterpret_cast<const NuBeamRunManager*>(G4RunManager::GetRunManager());
  fRecords = pRunManager->GetRecordPtr();
  fMessenger = new NuBeamEventActionMessenger(this);
}

// ---------------

NuBeamEventAction::~NuBeamEventAction()
{
  delete fMessenger;
}

// ----------------

void NuBeamEventAction::BeginOfEventAction(const G4Event* anEvent)
{
  // Initiate record-keeping for this event.
  if (fRecords != NULL) fRecords->RecordBeginOfEvent(anEvent);

  G4int evtNb = anEvent->GetEventID();
  //printing survey
  if (evtNb > 100)  fPrintModulo = 100;
  if (evtNb > 1000)  fPrintModulo = 1000;
  if (evtNb > 10000 && false)  fPrintModulo = 10000;
  if (evtNb%fPrintModulo == 0) {
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  }
}

void NuBeamEventAction::EndOfEventAction(const G4Event* anEvent)
{
  // Terminate record-keeping for this event.
  if (fRecords != NULL)
    fRecords->RecordEndOfEvent(anEvent);	
}

