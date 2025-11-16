#ifndef BooNETrajectoryContainer_h
#define BooNETrajectoryContainer_h 1

#include "NuBeamTrajectory.hh"
#include <vector>
#include <unordered_map> // make a hashmap of NuBeamTrajectories

// This class is a singleton whose only goal is to store the vector of NuBeamTrajectory objects
// later used in NuBeamOutput

class NuBeamTrajectoryContainer {
public:

  static NuBeamTrajectoryContainer & Instance();
  NuBeamTrajectoryContainer() = default;
  ~NuBeamTrajectoryContainer();
  
  NuBeamTrajectory & GetTrajectory(int id);
  inline size_t GetNTrajectories() { return fTrajectories.size(); }
  inline bool ContainsTrajectory(int id) { return fTrajectories.count(id); }

  // Add a trajectory
  void AddTrajectory(std::unique_ptr<NuBeamTrajectory> traj);

  // Which event are we on?
  inline G4int GetEvent() { return fEvtId; }
  void SetEvent(G4int id) { fEvtId = id; }

  // Clear in between events
  inline void Clear() { fTrajectories.clear(); }

private:
  
  // Singletons don't like the Rule of Three. So delete the copy and move operations
  NuBeamTrajectoryContainer(const NuBeamTrajectoryContainer&) = delete;
  NuBeamTrajectoryContainer & operator=(const NuBeamTrajectoryContainer &) = delete;

  std::unordered_map<int, std::unique_ptr<NuBeamTrajectory>> fTrajectories;
  // Keep track of the event to use
  G4int fEvtId;
};

#endif // # ifndef BooNETrajectoryContainer_h
