#include "NuBeamTrajectoryContainer.hh"
#include "G4ios.hh"
#include <exception>

NuBeamTrajectoryContainer & NuBeamTrajectoryContainer::Instance() { 
  static NuBeamTrajectoryContainer instance;
  return instance;
}

NuBeamTrajectoryContainer::~NuBeamTrajectoryContainer() {
  fTrajectories.clear();
}

NuBeamTrajectory NuBeamTrajectoryContainer::GetTrajectory(int id) {
  if( fTrajectories.find(id) == fTrajectories.end() ) {
    G4cout << "[NuBeamTrajectoryContainer]: " <<
      "Trajectory with id = " << id << " not found in trajectory map!!" << G4endl;
    throw std::runtime_error("Trajectory not handled");
  } else
    return fTrajectories[id]; // use of [id] means this can't be a const method
}

void NuBeamTrajectoryContainer::AddTrajectory(NuBeamTrajectory traj) {
  int par_id = traj.GetParentID();
  NuBeamTrajectory ins = traj;
  fTrajectories[par_id] = ins;
}
