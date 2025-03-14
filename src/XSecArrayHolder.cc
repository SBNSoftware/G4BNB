#include "XSecArrayHolder.hh"

///! Accesses instance
XSecArrayHolder & XSecArrayHolder::Instance() {
  static XSecArrayHolder instance;
  return instance;
}

XSecArrayHolder::XSecArrayHolder() : fIsInitialised(false) 
{
  for( int i = 0; i < kNProtonMomentumBins; i++ ) {
    for( int j = 0; j < kNPzBins; j++ ) {
      for( int k = 0; k < kNPtBins; k++ ) {
	ProtonXSecArray[i][j][k]    = 0.0;
	NeutronXSecArray[i][j][k]   = 0.0;
	PiPlusXSecArray[i][j][k]    = 0.0;
	PiMinusXSecArray[i][j][k]   = 0.0;
	KPlusXSecArray[i][j][k]     = 0.0;
	KMinusXSecArray[i][j][k]    = 0.0;
	KZeroLongXSecArray[i][j][k] = 0.0;

	ProtonNoWgtXSecArray[i][j][k]    = 0.0;
	NeutronNoWgtXSecArray[i][j][k]   = 0.0;
	PiPlusNoWgtXSecArray[i][j][k]    = 0.0;
	PiMinusNoWgtXSecArray[i][j][k]   = 0.0;
	KPlusNoWgtXSecArray[i][j][k]     = 0.0;
	KMinusNoWgtXSecArray[i][j][k]    = 0.0;
	KZeroLongNoWgtXSecArray[i][j][k] = 0.0;
      }
    }
  }
}

G4bool XSecArrayHolder::LoadXSecData( const G4double ProtonXSec[kNProtonMomentumBins][kNPzBins][kNPtBins],
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
				      const G4double KZeroLongNoWgtXSec[kNProtonMomentumBins][kNPzBins][kNPtBins] )
{
  if( fIsInitialised ) return false;

  for( int i = 0 ; i < kNProtonMomentumBins ; i++ ) {
    for( int j = 0 ; j < kNPzBins ; j++ ) { 
      for( int k = 0 ; k < kNPtBins ; k++ ) {
	ProtonXSecArray[i][j][k]    = ProtonXSec[i][j][k];
	NeutronXSecArray[i][j][k]   = NeutronXSec[i][j][k];
	PiPlusXSecArray[i][j][k]    = PiPlusXSec[i][j][k];
	PiMinusXSecArray[i][j][k]   = PiMinusXSec[i][j][k];
	KPlusXSecArray[i][j][k]     = KPlusXSec[i][j][k];
	KMinusXSecArray[i][j][k]    = KMinusXSec[i][j][k];
	KZeroLongXSecArray[i][j][k] = KZeroLongXSec[i][j][k];

	ProtonNoWgtXSecArray[i][j][k]    = ProtonNoWgtXSec[i][j][k];
	NeutronNoWgtXSecArray[i][j][k]   = NeutronNoWgtXSec[i][j][k];
	PiPlusNoWgtXSecArray[i][j][k]    = PiPlusNoWgtXSec[i][j][k];
	PiMinusNoWgtXSecArray[i][j][k]   = PiMinusNoWgtXSec[i][j][k];
	KPlusNoWgtXSecArray[i][j][k]     = KPlusNoWgtXSec[i][j][k];
	KMinusNoWgtXSecArray[i][j][k]    = KMinusNoWgtXSec[i][j][k];
	KZeroLongNoWgtXSecArray[i][j][k] = KZeroLongNoWgtXSec[i][j][k];
      } // loop
    } // over
  } // bins

  fIsInitialised = true;
  return true;
}

// Getters

G4double (&XSecArrayHolder::GetProtonXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return ProtonXSecArray;
}

G4double (&XSecArrayHolder::GetNeutronXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return NeutronXSecArray;
}

G4double (&XSecArrayHolder::GetPiPlusXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return PiPlusXSecArray;
}

G4double (&XSecArrayHolder::GetPiMinusXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return PiMinusXSecArray;
}

G4double (&XSecArrayHolder::GetKPlusXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return KPlusXSecArray;
}

G4double (&XSecArrayHolder::GetKMinusXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return KMinusXSecArray;
}

G4double (&XSecArrayHolder::GetKZeroLongXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return KZeroLongXSecArray;
}

G4double (&XSecArrayHolder::GetProtonNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return ProtonNoWgtXSecArray;
}

G4double (&XSecArrayHolder::GetNeutronNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return NeutronNoWgtXSecArray;
}

G4double (&XSecArrayHolder::GetPiPlusNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return PiPlusNoWgtXSecArray;
}

G4double (&XSecArrayHolder::GetPiMinusNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return PiMinusNoWgtXSecArray;
}

G4double (&XSecArrayHolder::GetKPlusNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return KPlusNoWgtXSecArray;
}

G4double (&XSecArrayHolder::GetKMinusNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return KMinusNoWgtXSecArray;
}

G4double (&XSecArrayHolder::GetKZeroLongNoWgtXSec()   )[kNProtonMomentumBins][kNPzBins][kNPtBins] 
{
  return KZeroLongNoWgtXSecArray;
}

G4bool XSecArrayHolder::IsInitialised() const{ return fIsInitialised; }
