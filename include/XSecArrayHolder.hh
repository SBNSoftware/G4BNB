#ifndef XSecArrayHolder_h
#define XSecArrayHolder_h 1

#include "globals.hh"

// include math package
#include <math.h>

class XSecArrayHolder
{
public:

  static XSecArrayHolder & Instance(); ///! Accesses instance

private:
  static const G4int kNPzBins=50;
  static const G4int kNPtBins=100;
  static const G4int kNProtonMomentumBins=100;

public:

  XSecArrayHolder();
  ~XSecArrayHolder() = default;

  // Delete copy c'tor and assignment op, keeps rule of 3 happy and this is a singleton
  XSecArrayHolder(const XSecArrayHolder&) = delete;
  XSecArrayHolder & operator= (const XSecArrayHolder &) = delete;

  ///! Gets an array from the instance
  G4double (&GetProtonXSec()    )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetNeutronXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetPiPlusXSec()    )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetPiMinusXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetKPlusXSec()     )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetKMinusXSec()    )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetKZeroLongXSec() )[kNProtonMomentumBins][kNPzBins][kNPtBins];

  G4double (&GetProtonNoWgtXSec()    )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetNeutronNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetPiPlusNoWgtXSec()    )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetPiMinusNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetKPlusNoWgtXSec()     )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetKMinusNoWgtXSec()    )[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double (&GetKZeroLongNoWgtXSec() )[kNProtonMomentumBins][kNPzBins][kNPtBins];

  ///! Loads array data only if not initialised
  G4bool LoadXSecData( const G4double ProtonXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double NeutronXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double PiPlusXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double PiMinusXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double KPlusXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double KMinusXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double KZeroLongXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double ProtonNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double NeutronNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double PiPlusNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double PiMinusNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double KPlusNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double KMinusNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
		       const G4double KZeroLongNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins] );

  G4bool IsInitialised() const;
  
private:

  ///! Members
  
  G4double ProtonXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double NeutronXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiPlusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiMinusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KPlusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KMinusXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KZeroLongXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];

  G4double ProtonNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double NeutronNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiPlusNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double PiMinusNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KPlusNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KMinusNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];
  G4double KZeroLongNoWgtXSecArray[kNProtonMomentumBins][kNPzBins][kNPtBins];

  G4bool fIsInitialised;

};

#endif
