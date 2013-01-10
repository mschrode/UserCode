// === Global Variables ================================================
// RA2 selection cuts
const float kHtJetPtMin   = 50.;
const float kHtJetEtaMax  = 2.5;



// === Main Script =====================================================
void hadTau2(const TString &inputEvents = "inputEvents", int nEvts = -1) {
  // --- Declare the Output Histograms ---------------------------------
  const int kNRespPtBins = 4;
  TH1* hTauResp[kNRespPtBins];
  for(int i = 0; i < kNRespPtBins; ++i) {
    TString name = "hTauResp_";
    name += i;
    hTauResp[i] = new TH1F(name,";p_{T}(visible) / p_{T}(generated);N(#tau)",60,-0.5,2.5);
  }


  // --- Declare the Variables Read from the Tree ----------------------
  // Array dimensions
  const int kRecoJetColSize = 15;
  const int kGenLepColSize = 6;

  // Reco-level jets
  int nRecoJets = 0;
  float recoJetPt[kRecoJetColSize];
  float recoJetPhi[kRecoJetColSize];
  float recoJetEta[kRecoJetColSize];

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
  tr->SetBranchAddress("NrecoJet",&nRecoJets);
  tr->SetBranchAddress("recoJetPt",recoJetPt);
  tr->SetBranchAddress("recoJetPhi",recoJetPhi);
  tr->SetBranchAddress("recoJetEta",recoJetEta);
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
    tr->GetEntry(evtIdx);
    if( nRecoJets > kRecoJetColSize ) {
      std::cerr << "ERROR: more than " << kRecoJetColSize << " reco jets in event " << evtIdx << std::endl;
      exit(-1);
    }


    // Select only events where the W decayed either
    // - into a muon (pdgId 13)
    // - into a tau (pdgId 15)
    if( flgW != 13 && flgW != 15 ) continue;
    // Select only events where the tau decayed hadronically
    if( flgW == 15 && flgTauHad != 1) continue;


    // Identify the charged generator-level lepton
    // from the W decay
    int chargedGenLepIdx = -1;
    for(int i = 0; i < kGenLepColSize; ++i) { // Loop over gen leptons
      if( TMath::Abs(genLepID[i]) == 13 || TMath::Abs(genLepID[i]) == 15 ) {
	chargedGenLepIdx = i;
	break;
      }
    } // End of loop over gen leptons
    if( chargedGenLepIdx < 0 ) {
      std::cerr << "ERROR: Could not find charged lepton from W decay in event " << evtIdx << std::endl;
      exit(-1);
    }

    
    // Compute some more kinematic quantities of the 
    // generator-level lepton from the W decay
    aux4Vector.SetPxPyPzE(genLepPx[chargedGenLepIdx],genLepPy[chargedGenLepIdx],genLepPz[chargedGenLepIdx],genLepE[chargedGenLepIdx]);
    float genLepEta = aux4Vector.Eta();
    float genLepPhi = aux4Vector.Phi();
    float genLepPt  = aux4Vector.Pt();



    // "Cross cleaning": find the jet that
    // - is the muon (each muon is also considered a Particle-Flow jet!)
    // - originates in the hadronic-tau decay
    // Associate jet and lepton if they are closer in eta-phi space
    // than deltaRMax
    int leptonJetIdx = -1;
    float deltaRMax = 0.1;
    if( flgW == 15 && genLepPt < 50. ) deltaRMax = 0.2;
    if( !findMatchedJet(leptonJetIdx,recoJetEta,recoJetPhi,nRecoJets,genLepEta,genLepPhi,deltaRMax) ) continue; // Want unambigious matching (exactly one jet within deltaRMax)


    // Calculate RA2 selection-variables from "cleaned" jets
    float selNJet = 0; // Number of jets with pt > 50 GeV and |eta| < 2.5 (HT jets)

    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
      // Skip this jet if it is the lepton
      if( jetIdx == leptonJetIdx ) continue;

      // Calculate NJet 
      if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) selNJet++;
    } // End of loop over reco jets

    
    // Select only events with at least 2 HT jets
    if( selNJet < 2 ) continue;
    

    // Fill histogram with relative visible energy of the tau
    // ("tau response template") for hadronically decaying tau
    if( flgW == 15 && flgTauHad == 1 ) {
      for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
	// Select tau jet
	if( jetIdx == leptonJetIdx ) {
	  int ptBin = tauPtBin(genLepPt);
	  if( ptBin >= 0 && ptBin < kNRespPtBins && TMath::Abs(genLepEta) < 2.1 ) {
	    hTauResp[ptBin]->Fill( recoJetPt[jetIdx]/genLepPt );
	  }
	  break;		// End the jet loop once the tau jet has been found
	}
      }	// End of loop over reco jets
    }
  } // End of loop over tree entries


  // --- Save the Histograms to File -----------------------------------
  TFile outFile("HadTau_WJetMC.root","UPDATE");
  for(int i = 0; i < kNRespPtBins; ++i) {
    hTauResp[i]->Write();
  }
  outFile.Close();
}




// === Auxiliary Functions =======================================

// Find index 'matchedJetIdx' of the jet that is within deltaRMax
// around the lepton. Returns
//  - true  : if exactly one jet has been found.
//  - false : otherwise. In that case, 'matchedJetIdx == -1'
bool findMatchedJet(int &matchedJetIdx, const float* jetEta, const float* jetPhi, int nJets, float lepEta, float lepPhi, float deltaRMax) {
  bool unambigiousMatch = false;
  matchedJetIdx = -1;
  for(int jetIdx = 0; jetIdx < nJets; ++jetIdx) { // Loop over jets
    if( deltaR(jetEta[jetIdx],lepEta,jetPhi[jetIdx],lepPhi) < deltaRMax ) {
      // Check if a jet has already been matched
      if( unambigiousMatch ) {
	matchedJetIdx = -1;
	unambigiousMatch = false; // We only want exactly one matched jet
	break;
      }
      matchedJetIdx = jetIdx;
      unambigiousMatch = true;
    }
  } // End of loop over jets

  return unambigiousMatch;
}


float deltaR(float eta1, float eta2, float phi1, float phi2) {
  float dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  float deta = eta1 - eta2;

  return sqrt( deta*deta + dphi*dphi );
}


int tauPtBin(float pt) {
  int bin = -1;
  if( pt >= 20 )  bin = 0;
  if( pt >= 30 )  bin = 1;
  if( pt >= 50 )  bin = 2;
  if( pt >= 100 ) bin = 3;

  return bin;
}
