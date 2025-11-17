#include "NuBeamTrackInformation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"

G4Allocator<NuBeamTrackInformation> aTrackInformationAllocator;

NuBeamTrackInformation::NuBeamTrackInformation()
  :
  fCreatorModelName(""),
  fDecayCodeDk2nu(0),
  fIsBiassed(false),
  fCreatorWasQE(false),
  fCreatorNProtons(0),
  fCreatorNNeutrons(0),
  fCreatorNPiPlus(0),
  fCreatorNPiMinus(0),
  fCreatorNKPlus(0),
  fCreatorNKMinus(0),
  fCreatorNOthers(0)
{   
}

NuBeamTrackInformation::NuBeamTrackInformation(const G4Track*)
  :
  fIsBiassed(false),
  fCreatorWasQE(false),
  fCreatorNProtons(0),
  fCreatorNNeutrons(0),
  fCreatorNPiPlus(0),
  fCreatorNPiMinus(0),
  fCreatorNKPlus(0),
  fCreatorNKMinus(0),
  fCreatorNOthers(0)
{
}

NuBeamTrackInformation::NuBeamTrackInformation(const NuBeamTrackInformation* aTrackInfo)
{
  fCreatorModelName = aTrackInfo->fCreatorModelName;
  fDecayCodeDk2nu = aTrackInfo->fDecayCodeDk2nu;
  fIsBiassed = aTrackInfo->fIsBiassed;
  fCreatorWasQE = aTrackInfo->fCreatorWasQE;
  fCreatorNProtons = aTrackInfo->fCreatorNProtons;
  fCreatorNNeutrons = aTrackInfo->fCreatorNNeutrons;
  fCreatorNPiPlus = aTrackInfo->fCreatorNPiPlus;
  fCreatorNPiMinus = aTrackInfo->fCreatorNPiMinus;
  fCreatorNKPlus = aTrackInfo->fCreatorNKPlus;
  fCreatorNKMinus = aTrackInfo->fCreatorNKMinus;
  fCreatorNOthers = aTrackInfo->fCreatorNOthers;
}

NuBeamTrackInformation::~NuBeamTrackInformation()
{
}

void NuBeamTrackInformation::Print() const
{
}
