// === Main Script =====================================================
void hadTau1(const TString &inputEvents = "inputEvents", int nEvts = -1) {
  // --- Declare the Output Histograms ---------------------------------
  TH1* hGenMuPt = new TH1F("hGenMuPt",";p_{T}(#mu^{gen}) [GeV];N(#mu^{gen})",20,0,500);
  hGenMuPt->Sumw2();		// The error per bin will be computed as sqrt(sum of squares of weight) for each bin. Important for errors in ratio plot
  TH1* hGenTauPt = static_cast<TH1*>(hGenMuPt->Clone("hGenTauPt"));
  hGenTauPt->SetTitle(";p_{T}(#tau^{gen}) [GeV];N(#tau^{gen})");


  // --- Declare the Variables Read from the Tree ----------------------
  // Array dimensions
  const int kGenLepColSize = 6;

  // W-decay mode
  int flgW = 0;			// PdgId of lepton the W decays into
  int flgTauHad = 0;		// 1: fully-hadronic decay

  // Gen-level leptons from W decay
  float genLepPx[kGenLepColSize];
  float genLepPy[kGenLepColSize];
  float genLepPz[kGenLepColSize];
  float genLepE[kGenLepColSize];
  float genLepID[kGenLepColSize];



  // --- Set Up the Tree -----------------------------------------------

  // Get the tree from file
  TChain* tr = new TChain("AnaTree");
  tr->Add(inputEvents);

  // Set the branches
  tr->SetBranchAddress("flgW",&flgW);
  tr->SetBranchAddress("flgTauHad",&flgTauHad);
  tr->SetBranchAddress("lepPx",genLepPx);
  tr->SetBranchAddress("lepPy",genLepPy);
  tr->SetBranchAddress("lepPz",genLepPz);
  tr->SetBranchAddress("lepE",genLepE);
  tr->SetBranchAddress("lepID",genLepID);


  // --- Process the Events in the Tree --------------------------------
  int nEvtsToProcess = tr->GetEntries();
  if( nEvts > 0 && nEvts < nEvtsToProcess ) nEvtsToProcess = nEvts;
  std::cout << "Processing " << nEvtsToProcess << " events" << std::endl;

  // Auxiliary variables (initialize outside event loop
  // for performance reasons)
  TLorentzVector aux4Vector;

  // Loop over the tree entries
  for(int evtIdx = 0; evtIdx < nEvtsToProcess; ++evtIdx) {
    if( evtIdx%10000 == 0 ) std::cout<<"  Event: " << evtIdx << std::endl;

    // Get the variables' values for this event
    // from tree
    tr->GetEntry(evtIdx);


    // Select only events where the W decayed either
    // - into a muon (pdgId 13)
    // - into a tau (pdgId 15)
    if( flgW != 13 && flgW != 15 ) continue;
    // Select only events where the tau decayed hadronically
    if( flgW == 15 && flgTauHad != 1) continue;


    // Identify the charged generator-level lepton (muon or tau)
    // from the W decay.
    int chargedGenLepIdx = -1;
    if( !findMuOrTau(chargedGenLepIdx,genLepID,kGenLepColSize) ) continue; // Want exactly one lepton
    
    // Compute some more kinematic quantities of the 
    // generator-level lepton
    aux4Vector.SetPxPyPzE(genLepPx[chargedGenLepIdx],genLepPy[chargedGenLepIdx],genLepPz[chargedGenLepIdx],genLepE[chargedGenLepIdx]);
    float genLepEta = aux4Vector.Eta();
    float genLepPhi = aux4Vector.Phi();
    float genLepPt  = aux4Vector.Pt();


    // Fill generator-level lepton pt into histograms
    if( flgW == 13 ) hGenMuPt->Fill(genLepPt);
    else if( flgW == 15 ) hGenTauPt->Fill(genLepPt);

  } // End of loop over tree entries

  // --- Save the Histograms to File -----------------------------------
  TFile outFile("HadTau_WJetMC.root","RECREATE");
  hGenMuPt->Write();
  hGenTauPt->Write();
  outFile.Close();
}




// === Auxiliary Functions =======================================

// Find index 'lepIdx' of muon or tau in lepton collection
// Returns
// - true  : if exactly one lepton (muon or tau) has been found
// - false : otherwise. In that case, 'chargedLepIdx == -1'
bool findMuOrTau(int &lepIdx, const float* lepID, int nLep) {
  bool unambigiousHit = false;
  lepIdx = -1;
  for(int i = 0; i < nLep; ++i) { // Loop over leptons
    // Look for muon (pdgId 13) or tau (pdgId 15)
    if( TMath::Abs(lepID[i]) == 13 || TMath::Abs(lepID[i]) == 15 ) {
      if( unambigiousHit ) {
	lepIdx = -1;
	unambigiousHit = false;	// We only want exactly one match
	break;
      }
      lepIdx = i;
      unambigiousHit = true;
    }
  } // End of loop over leptons
  
  return unambigiousHit;
}

