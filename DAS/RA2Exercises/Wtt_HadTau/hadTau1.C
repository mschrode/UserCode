// === Global Variables ================================================
// RA2 selection cuts
const float kHtJetPtMin   = 50.;
const float kHtJetEtaMax  = 2.5;



// === Main Script =====================================================
void hadTau1(int nEvts = -1) {
  // --- Declare the Output Histograms ---------------------------------
  TH1* hGenMuPt = new TH1F("hGenMuPt",";p_{T}(#mu^{gen});N(#mu^{gen})",30,0,300);
  TH1* hGenTauPt = new TH1F("hGenTauPt",";p_{T}(#tau^{gen});N(#tau^{gen})",30,0,300);


  // --- Declare the Variables Read from the Tree ----------------------
  // Array dimensions
  const unsigned int kRecoJetColSize = 200;
  const unsigned int kGenLepColSize = 6;

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
  tr->Add("/nfs/dust/test/cmsdas/school61/susy/ntuple/test/WJets.root");

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
  std::set<int> auxIdcsLeptonJets;

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
    // generator-level lepton (electron or the muon)
    // from the W decay
    aux4Vector.SetPxPyPzE(genLepPx[chargedGenLepIdx],genLepPy[chargedGenLepIdx],genLepPz[chargedGenLepIdx],genLepE[chargedGenLepIdx]);
    float genLepEta = aux4Vector.Eta();
    float genLepPhi = aux4Vector.Phi();
    float genLepPt  = aux4Vector.Pt();



    // "Cross cleaning": find the jet(s) that
    // - is the muon (each muon is also considered a Particle-Flow jet!)
    // - originate in the hadronic-tau decay
    // Associate jets and leptons if they are closer in eta-phi space
    // than deltaRMax
    auxIdcsLeptonJets.clear();
    float deltaRMax = 0.1;
    if( flgW == 15 && genLepPt < 50. ) deltaRMax = 0.2;
    findMatchedJets(auxIdcsLeptonJets,recoJetEta,recoJetPhi,nRecoJets,genLepEta,genLepPhi,deltaRMax);


    // Calculate RA2 selection-variables from "cleaned" jets
    float selNJet = 0; // Number of jets with pt > 50 GeV and |eta| < 2.5 (HT jets)

    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
      // Skip this jet if it is the lepton
      if( auxIdcsLeptonJets.find(jetIdx) != auxIdcsLeptonJets.end() ) continue;

      // Calculate NJet 
      if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) selNJet++;
    } // End of loop over reco jets

    
    // Select only events with at least 2 HT jets
    if( selNJet < 2 ) continue;

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

void findMatchedJets(std::set<int> &jetIdcs, const float* jetEta, const float* jetPhi, int nJets, float lepEta, float lepPhi, float deltaRMax) {
  for(int jetIdx = 0; jetIdx < nJets; ++jetIdx) { // Loop over jets
    if( deltaR(jetEta[jetIdx],lepEta,jetPhi[jetIdx],lepPhi) < deltaRMax ) {
      jetIdcs.insert(jetIdx);
    }
  } // End of loop over jets
}


float deltaR(float eta1, float eta2, float phi1, float phi2) {
  float deta, dphi, dr;
  
  dphi = TMath::Abs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
  deta = TMath::Abs(eta1- eta2);

  dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

