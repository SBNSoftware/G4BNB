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
  
  NuBeamTrajectory GetTrajectory(int id);
  inline size_t GetNTrajectories() { return fTrajectories.size(); }

  // Add a trajectory
  // Don't add by reference, as this shares a pointer to the PositionRecord of the NBTrajectory...
  void AddTrajectory(NuBeamTrajectory traj);

  // Clear in between events
  inline void Clear() { fTrajectories.clear(); }

private:
  
  // Singletons don't like the Rule of Three. So delete the copy and move operations
  NuBeamTrajectoryContainer(const NuBeamTrajectoryContainer&) = delete;
  NuBeamTrajectoryContainer & operator=(const NuBeamTrajectoryContainer &) = delete;

  std::unordered_map<int, NuBeamTrajectory> fTrajectories;
};

#endif // # ifndef BooNETrajectoryContainer_h
