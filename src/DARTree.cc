#include "DARTree.hh"

// Initialisation of static member
TTree * DARTree::fDARTree = 0;

TTree * DARTree::RetrieveTree() {
  if( ! fDARTree ) {
    fDARTree = new TTree("rest_frame", "Rest frame decay ntuple");
  }
  return fDARTree;
}
