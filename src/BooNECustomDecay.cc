#include "BooNECustomDecay.hh"

BooNECustomDecay::BooNECustomDecay(const G4String & processName) 
  : G4Decay(processName) {
}

G4VParticleChange * BooNECustomDecay::AtRestDoIt(const G4Track & track,
						 const G4Step & step)
{
  return this->PostStepDoIt(track, step);
}

G4VParticleChange * BooNECustomDecay::PostStepDoIt(const G4Track & track,
						   const G4Step & step)
{
  //return G4Decay::PostStepDoIt(track, step);

  // Clear the decay tuple
  fNtupleDAR.nmult = 0;
  fNtupleDAR.ppdg = 0;
  for( G4int im = 0; im < max_multiplicity; im++ ){
    fNtupleDAR.pdg[im] = 0;
    fNtupleDAR.p4[im][0] = 0.0;
    fNtupleDAR.p4[im][1] = 0.0;
    fNtupleDAR.p4[im][2] = 0.0;
    fNtupleDAR.p4[im][3] = 0.0;
    fNtupleDAR.wgt[im] = 0.0;
  }

  // Lift the code from https://apc.u-paris.fr/~franco/g4doxy/html/classG4Decay.html
  // Initialize ParticleChange
  fParticleChangeForDecay.Initialize(track);

  // Get particle
  const G4DynamicParticle * aParticle = track.GetDynamicParticle();
  const G4ParticleDefinition * aParticleDef = aParticle->GetDefinition();

  // Check if it's stable
  if( aParticleDef->GetPDGStable() ) return &fParticleChangeForDecay;

  // Check if thePreAssignedDecayProducts exists
  const G4DecayProducts * o_products = (aParticle->GetPreAssignedDecayProducts());
  G4bool isPreAssigned = (o_products != 0);
  G4DecayProducts * products = 0;
  G4DecayProducts * rf_products = 0;

  // Decay table
  G4DecayTable * decaytable = aParticleDef->GetDecayTable();
  // Dump it
  //if( decaytable ) decaytable->DumpInfo();

  // Check if external decayer exists
  G4bool isExtDecayer = (decaytable == 0) && (pExtDecayer != 0);

  // Error due to NO decay table
  if ( (decaytable == 0) && !isExtDecayer && !isPreAssigned ) {
    if (GetVerboseLevel() > 0) {
      G4cout << "BooNECustomDecay::DoIt : decay table not defined for "
	     << aParticle->GetDefinition()->GetParticleName() << G4endl;
    }
    G4Exception( "BooNECustomDecay::DecayIt", "DECAY101", JustWarning, "Decay table not defined" );
    fParticleChangeForDecay.SetNumberOfSecondaries(0);
    // Kill the parent particle
    fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill );
    fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.0);

    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForDecay;
  }
  
  G4VDecayChannel * decaychannel = 0;
  if( isPreAssigned ) {
    products = new G4DecayProducts(*o_products);
  } else if( isExtDecayer ) {
    products = pExtDecayer->ImportDecayProducts(track);
  } else {
    // Choose a decay channel and tell me what is picked
    decaychannel = decaytable->SelectADecayChannel();
    if( decaychannel == 0 ) {
      G4Exception("BooNECustomDecay::DoIt", "DECAY003", FatalException, "Cannot determine decay channel");
    } else {
      //G4cout << "Picked a decay channel! Dumping it" << G4endl;
      //decaychannel->DumpInfo();
      products = decaychannel->DecayIt(aParticle->GetMass());
      //G4cout << "Dumping the products, too" << G4endl;
      //products->DumpInfo();
    }
  }

  // Get parent particle information
  G4double ParentEnergy = aParticle->GetTotalEnergy();
  G4double ParentMass   = aParticle->GetMass();
  if( ParentEnergy < ParentMass ) {
    G4Exception( "BooNECustomDecay::DoIt", "DECAY102", JustWarning, "Setting energy = mass" );
    ParentEnergy = ParentMass;
  }

  G4ThreeVector ParentDirection(aParticle->GetMomentumDirection());

  // Boost all decay products to lab frame
  rf_products = new G4DecayProducts(*products);
  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = track.GetGlobalTime();
  G4double finalLocalTime  = track.GetLocalTime();
  if( track.GetTrackStatus() == fStopButAlive ) {
    // if at rest
    finalGlobalTime += fRemainderLifeTime;
    finalLocalTime  += fRemainderLifeTime;
    energyDeposit   += aParticle->GetKineticEnergy();
    if( isPreAssigned ) products->Boost( ParentEnergy, ParentDirection );
  } else {
    // PostStep case
    if( !isExtDecayer ) products->Boost( ParentEnergy, ParentDirection );
  }

  // Set polarisation for daughter particles
  DaughterPolarization(track, products);

  // Add products in fParticleChangeForDecay
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForDecay.SetNumberOfSecondaries(numberOfSecondaries);

  // Only dump if there is a neutrino in the decay channel
  G4bool found_a_neutrino = false;
  std::unique_ptr<G4DecayProducts> poproducts = std::make_unique<G4DecayProducts>(*products);
  while( poproducts->entries() > 0 ) {
    G4DynamicParticle * last = poproducts->PopProducts();
    G4int pdg = last->GetPDGcode();
    if( std::abs(pdg) == G4NeutrinoE::NeutrinoE()->GetPDGEncoding() ||
	std::abs(pdg) == G4NeutrinoMu::NeutrinoMu()->GetPDGEncoding() || 
	std::abs(pdg) == G4NeutrinoTau::NeutrinoTau()->GetPDGEncoding() ) found_a_neutrino = true;
  }
  if( found_a_neutrino ) {
    /*
    G4cout << "APAPAPA. CALLING DECAY FOR TRACK WITH ID " << track.GetTrackID()
	   << " AND PARTICLE NAME " << track.GetDefinition()->GetParticleName() << G4endl;
    G4cout << "Dumping decay channel with " << products->entries() << " entries:" << G4endl;
    
    decaychannel->DumpInfo();
    G4cout << "Dumping decay products in lab frame:" << G4endl;
    products->DumpInfo();
    G4cout << "Dumping decay products in rest frame:" << G4endl;
    rf_products->DumpInfo();
    */

    // Save a neutrino, regardless of if it will make RecordNeutrino() !
    if( /* fSaveRestFrameDecay && */ true ) {
      TTree * fDARTree = DARTree::RetrieveTree();
      if( fDARTree->GetListOfBranches()->GetEntries() == 0 ) {
	fDARTree->Branch( "nmult", &fNtupleDAR.nmult, "nmult/I" );
	fDARTree->Branch( "ppdg", &fNtupleDAR.ppdg, "ppdg/I" );
	fDARTree->Branch( "pdg", fNtupleDAR.pdg, "pdg[nmult]/I" );
	fDARTree->Branch( "p4", fNtupleDAR.p4, "p4[nmult][4]/F" );
	fDARTree->Branch( "wgt", fNtupleDAR.wgt, "wgt[nmult]/F" );

	fNtupleDAR.nmult = 0;
	fNtupleDAR.ppdg = 0;
	for( G4int im = 0; im < max_multiplicity; im++ ){
	  fNtupleDAR.pdg[im] = 0;
	  fNtupleDAR.p4[im][0] = 0.0;
	  fNtupleDAR.p4[im][1] = 0.0;
	  fNtupleDAR.p4[im][2] = 0.0;
	  fNtupleDAR.p4[im][3] = 0.0;
	  fNtupleDAR.wgt[im] = 0.0;
	}
      }

      fNtupleDAR.nmult = rf_products->entries();
      fNtupleDAR.ppdg = aParticle->GetDefinition()->GetPDGEncoding();

      poproducts = std::make_unique<G4DecayProducts>(*rf_products); // Reassign memory
      //G4cout << "There are " << poproducts->entries() << " in the reassigned products" << G4endl;
      G4int imult = -1;
      while( poproducts->entries() > 0 ) {
	G4DynamicParticle * last = poproducts->PopProducts();
	G4int pdg = last->GetPDGcode();
	//G4cout << "imult = " << imult+1 << " ppdg = " << aParticle->GetDefinition()->GetPDGEncoding()
	//       << " pdg = " << pdg << G4endl;
	G4LorentzVector p4 = last->Get4Momentum();

	fNtupleDAR.pdg[++imult] = pdg;
	fNtupleDAR.p4[imult][0] = p4.px();
	fNtupleDAR.p4[imult][1] = p4.py();
	fNtupleDAR.p4[imult][2] = p4.pz();
	fNtupleDAR.p4[imult][3] = p4.e();
	fNtupleDAR.wgt[imult] = track.GetWeight();
      }
      fDARTree->Fill();
    }
  }

  G4int index;
  G4ThreeVector currentPosition;
  const G4TouchableHandle thand = track.GetTouchableHandle();
  for( index = 0; index < numberOfSecondaries; index++ ) {
    // Get current position of the track
    currentPosition = track.GetPosition();
    // Create a new track object
    G4Track * secondary = new G4Track( products->PopProducts(),
				       finalGlobalTime, currentPosition );
    // Switch on good-for-tracking flag
    secondary->SetGoodForTrackingFlag();
    secondary->SetTouchableHandle(thand);
    // Add the secondary to the list
    fParticleChangeForDecay.AddSecondary(secondary);
  }
  delete products; delete rf_products;

  // Kill the parent particle
  fParticleChangeForDecay.ProposeTrackStatus(fStopAndKill);
  fParticleChangeForDecay.ProposeLocalEnergyDeposit(energyDeposit);
  fParticleChangeForDecay.ProposeLocalTime(finalLocalTime);

  // Clear the number of interaction lengths
  ClearNumberOfInteractionLengthLeft();

  return &fParticleChangeForDecay;
}
