#ifndef DARTree_h
#define DARTree_h 1

#include "TTree.h"

class DARTree {
public:
  // Retrieves the tree singleton
  static TTree * RetrieveTree();

private:
  // Private c'tor prevents instantiation
  DARTree();
  static TTree * fDARTree;
};

#endif // #ifndef DARTree_h
