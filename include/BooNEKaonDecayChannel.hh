/*
 * \class BooNEKaonDecayChannel
 *
 * \brief Implements the correct neutrino energy spectrum for K+- and K0L decays to neutrinos.
 *        The method in G4KL3DecayChannel.hh is wrong.
 *         
 *        Essentially reimplements G4KL3DecayChannel to use the miniBooNE custom decay model.
 *        Refer to arXiv:0806.1449 section III.F
 *
 * \date  Jan 13th, 2025
 */

#ifndef BooNEKaonDecayChannel_h
#define BooNEKaonDecayChannel_h 1

#include <cmath>
#include "globals.hh"
#include "G4ios.hh"
#include "G4VDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"

// the particles
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"

class BooNEKaonDecayChannel: public G4VDecayChannel
{
public:
  BooNEKaonDecayChannel(const G4String & theParentName,
			G4double         theBR,
			const G4String & thePionName,
			const G4String & theLeptonName,
			const G4String & theNeutrinoName);

  BooNEKaonDecayChannel(const G4String & theParentName,
			G4double         theBR,
			const G4String & thePionName,
			const G4String & theLeptonName,
			const G4String & theNeutrinoName,
			const G4double theLambda, const G4double theXi0);
  
  virtual ~BooNEKaonDecayChannel();

protected:
  // Copy c'tor and assignment are protected
  BooNEKaonDecayChannel(const BooNEKaonDecayChannel &);
  BooNEKaonDecayChannel & operator=(const BooNEKaonDecayChannel &);

private:
  BooNEKaonDecayChannel();

public:
  virtual G4DecayProducts * DecayIt(G4double);

protected:
  // yes we use arrays of pointers.
  enum{idPi=0, idLepton=1, idNeutrino=2};

protected:
  // This method samples a random point from the allowed phase space.
  // It stores the candidate momentum and energy in E,Pdaughter[3]
  void PhaseSpace( G4double Mparent,
		   const G4double * Mdaughter,
		   G4double *       Edaughter,
		   G4double *       Pdaughter );

  // This method constructs the Dalitz density for a pure V-A decay
  // Modelled after Chounet, Gaillard and Gaillard, Phys Rep 4 (1972) 199-324
  G4double DalitzDensity( G4double Mparent,
			  G4double Epi, G4double El, G4double Enu,
			  G4double Mpi, G4double Ml/*, G4double Mnu */ );

private:
  // Constants for DalitzDensity()
  G4double pLambda, pXi0;

public:
  void SetDalitzParameters( G4double aLambda, G4double aXi );
  G4double GetDalitzLambda() const;
  G4double GetDalitzXi() const;

private:
  // The PDG masses we use, PDG 2002
  const G4double fKPlusMass  = 0.493677;
  const G4double fKZeroMass  = 0.497611;
  const G4double fPiPlusMass = 0.13957;
  const G4double fPiZeroMass = 0.13498;
  const G4double fMuonMass   = 0.10566;
  const G4double fElectronMass = 0.000511;
};

// Keep same imp as G4KL3DecayChannel... ugh
inline void BooNEKaonDecayChannel::SetDalitzParameters(G4double aLambda, G4double aXi)
{
  pLambda = aLambda;
  pXi0    = aXi;
}

inline G4double BooNEKaonDecayChannel::GetDalitzLambda() const { return pLambda; }
inline G4double BooNEKaonDecayChannel::GetDalitzXi() const { return pXi0; }

#endif // #ifndef BooNEKaonDecayChannel_h
