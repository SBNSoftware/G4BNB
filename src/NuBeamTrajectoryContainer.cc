#include "NuBeamTrajectoryContainer.hh"
#include "G4ios.hh"

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
    NuBeamTrajectory traj = NuBeamTrajectory();
    return traj;
  } else
    return fTrajectories[id]; // use of [id] means this can't be a const method
}

void NuBeamTrajectoryContainer::AddTrajectory(NuBeamTrajectory traj) {
  int par_id = traj.GetParentID();
  /*
  if( fTrajectories.find(par_id) != fTrajectories.end() ) {
    G4cout << "[NuBeamTrajectoryContainer]: " <<
      "WARNING: Trajectory with id = " << par_id << 
      " already in trajectory map!!! Not overriding" << G4endl;
    return;
  } else
  */
  NuBeamTrajectory ins = traj;
  fTrajectories[par_id] = ins;
}
