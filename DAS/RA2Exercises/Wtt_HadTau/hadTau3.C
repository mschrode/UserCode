// === Global Variables ================================================
const float kHtJetPtMin   = 50.;
const float kHtJetEtaMax  = 2.5;
const float kMhtJetPtMin  = 30.;
const float kMhtJetEtaMax = 5.0;
const float kMuonPtMin    = 20.;



// === Main Script =====================================================
void hadTau3(int nSimSteps,
	     const TString &inputEvents = "inputEvents",
	     const TString &respTempl = "HadTau_WJetMC.root",
	     int nEvts = -1) {

  const int kNSimStepsMax = 500;
  if( nSimSteps > kNSimStepsMax ) {
    std::cerr << "ERROR: Number of simulation steps larger than allowed." << std::endl;
    std::cerr << "       'nSimSteps' must be less than " << kNSimStepsMax << std::endl;
    exit(-1);
  }


  // --- Binning for Prediction and Closure Test -----------------------
  const int kBinsPtN = 30;
  const float kBinsPtMin = 0.;
  const float kBinsPtMax = 300.;
  const int kBinsHtN = 30;
  const float kBinsHtMin = 0.;
  const float kBinsHtMax = 3000.;
  const int kBinMhtN = 20;
  const float kBinsMhtMin = 0.;
  const float kBinsMhtMax = 1000.;

  // Counters: store for each simulation step, how many events
  // are predicted in each pt / HT / MHT bin
  int nPredTauJetPt2D[kNSimStepsMax][kBinsPtN];
  int nPredHT2D[kNSimStepsMax][kBinsHtN];
  int nPredMHT2D[kNSimStepsMax][kBinMhtN];
  for(int i = 0; i < nSimSteps; ++i) {
    for(int j = 0; j < kBinsPtN; ++j) nPredTauJetPt2D[i][j] = 0;
    for(int j = 0; j < kBinsHtN; ++j) nPredHT2D[i][j] = 0;
    for(int j = 0; j < kBinMhtN; ++j) nPredMHT2D[i][j] = 0;
  }



  // --- Declare the Output Histograms ---------------------------------
  TH1* hPredTauJetPt = new TH1F("hPredTauJetPt",";Predicted p_{T}(#tau) [GeV];N(#tau)",kBinsPtN,kBinsPtMin,kBinsPtMax);   
  TH1* hPredHT = new TH1F("hPredHT",";Predicted H_{T} [GeV];N(events)",kBinsHtN,kBinsHtMin,kBinsHtMax);   
  TH1* hPredMHT = new TH1F("hPredMHT",";Predicted #slash{H}_{T} [GeV];N(events)",kBinMhtN,kBinsMhtMin,kBinsMhtMax);   
  TH1* hTrueTauJetPt = new TH1F("hTrueTauJetPt",";True p_{T}(#tau) [GeV];N(#tau)",kBinsPtN,kBinsPtMin,kBinsPtMax);   
  TH1* hTrueHT = new TH1F("hTrueHT",";True H_{T} [GeV];N(events)",kBinsHtN,kBinsHtMin,kBinsHtMax);   
  TH1* hTrueMHT = new TH1F("hTrueMHT",";True #slash{H}_{T} [GeV];N(events)",kBinMhtN,kBinsMhtMin,kBinsMhtMax);   



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



  // --- Get the Tau Templates from File -------------------------------
  TFile fTauResp(respTempl,"READ");
  const int kNRespPtBins = 4;
  TH1* hTauResp[kNRespPtBins];
  for(int i = 0; i < kNRespPtBins; ++i) {
    hTauResp[i] = 0;
    TString name = "hTauResp_";
    name += i;
    fTauResp.GetObject(name,hTauResp[i]);
    if( hTauResp[i] ) {
      hTauResp[i]->SetDirectory(0);
    } else {
      std::cerr << "ERROR: Histogram '" << name << "' not found in file '" << fTauResp.GetName() << "'" << std::endl;
      exit(-1);
    }
  }
  fTauResp.Close();



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
    float selHt   = 0.;
    float selMhtX = 0.;
    float selMhtY = 0.;
    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
      // Skip this jet if it is the lepton
      if( jetIdx == leptonJetIdx ) continue;

      // Calculate NJet and HT
      if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) {
	selNJet++;
	selHt += recoJetPt[jetIdx];
      }
      // Calculate MHT
      if( recoJetPt[jetIdx] > kMhtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kMhtJetEtaMax ) {
	selMhtX -= recoJetPt[jetIdx]*TMath::Cos(recoJetPhi[jetIdx]);
	selMhtX -= recoJetPt[jetIdx]*TMath::Sin(recoJetPhi[jetIdx]);
      }
    } // End of loop over reco jets


    // Select only events with at least 2 HT jets
    if( selNJet < 2 ) continue;

    
    // Now, the lepton
    if( flgW == 15 ) {
      // In case the W decayed into a hadronically decaying tau
      // - store the reco tau-jet pt
      // - add the reco tau-jet pt to HT and MHT
      float tauJetPt  = recoJetPt[leptonJetIdx];
      float tauJetEta = TMath::Abs(recoJetEta[leptonJetIdx]);
      float tauJetPhi = recoJetPhi[leptonJetIdx];

      // If tau jet meets same criteria as other RA2 jets for HT,
      // recompute NJets and check if NJets >= 3
      if( tauJetPt > kHtJetPtMin && tauJetEta < kHtJetEtaMax ) selNJet++;
      if( selNJet < 3 ) continue;

      // Apply same pt cut as for muon: response templates
      // available for pt > 20
      if( tauJetPt > kMuonPtMin ) {
	hTrueTauJetPt->Fill(tauJetPt);

	// If tau jet meets same criteria as other RA2 jets for HT,
	// add tau-jet pt to HT and fill histogram of true HT
	if( tauJetPt > kHtJetPtMin && tauJetEta < kHtJetEtaMax ) {
	  selHt += tauJetPt;
	  hTrueHT->Fill(selHt);
	}
	// If tau jet meets same criteria as other RA2 jets for MHT,
	// add tau-jet pt to MHT and fill histogram of true MHT
	if( tauJetPt > kMhtJetPtMin && tauJetEta < kMhtJetEtaMax ) {
	  selMhtX -= tauJetPt*cos(tauJetPhi);
	  selMhtY -= tauJetPt*sin(tauJetPhi);
	  // Calculate the MHT from x and y component
	  float selMht = sqrt( selMhtX*selMhtX + selMhtY*selMhtY );
	  hTrueMHT->Fill(selMht);
	}
      }
    } else if( flgW == 13 ) {
      // In case the W decayed into a muon
      // - scale its pt by a random factor drawn from the
      //   tau-reponse template to simulate the tau measurement
      // - use simulated tau-pt to predict HT and MHT

      // Response templates available for pt > 20
      if( genLepPt < kMuonPtMin ) continue;

      // Find the pt bin of the tau-response template
      int tauRespPtBin = tauPtBin(genLepPt);
      if( tauRespPtBin < 0 || tauRespPtBin >= kNRespPtBins ) continue;

      // Perform nSimSteps simulations
      for(int it = 0; it < nSimSteps; ++it) {
	// Get random number from tau-response template
	float scale = hTauResp[tauRespPtBin]->GetRandom();
	// Scale muon pt with tau response --> simulate tau jet pt
	float simTauJetPt = scale * genLepPt;

	// If simulted tau-jet meets same criteria as other RA2 jets for HT,
	// recompute NJets and check if NJets >= 3
	float simNJet = selNJet;
	if( simTauJetPt > kHtJetPtMin && TMath::Abs(genLepEta) < kHtJetEtaMax ) simNJet++;
	if( simNJet < 3 ) continue;
	  
	// Tau-jet pt bin
	int binSimTauJetPt = bin(simTauJetPt,kBinsPtN,kBinsPtMin,kBinsPtMax);
	if( binSimTauJetPt > -1 ) nPredTauJetPt2D[it][binSimTauJetPt]++;
	  	
	// If simulated tau-jet meets same criteria as RA2 jets for HT,
	// add simulated tau-jet pt to HT and increment number of events
	// with that HT
	if( simTauJetPt > kHtJetPtMin && TMath::Abs(genLepEta) < kHtJetEtaMax ) {
	  float simHt = selHt + simTauJetPt;
	  int binSimHT = bin(simHt,kBinsHtN,kBinsHtMin,kBinsHtMax);
	  if( binSimHT > -1 ) nPredHT2D[it][binSimHT]++;
	}
	// If simulated tau-jet meets same criteria as RA2 jets for MHT,
	// add simulated tau-jet pt to MHT and increment number of events
	// with that MHT
	if( simTauJetPt > kMhtJetPtMin && TMath::Abs(genLepEta) < kMhtJetEtaMax ) {
	  float simMhtX = selMhtX - simTauJetPt*cos(genLepPhi);
	  float simMhtY = selMhtY - simTauJetPt*sin(genLepPhi);
	  float simMht = sqrt( simMhtX*simMhtX + simMhtY*simMhtY );
	  int binSimMHT = bin(simMht,kBinMhtN,kBinsMhtMin,kBinsMhtMax);
	  if( binSimMHT > -1 ) nPredMHT2D[it][binSimMHT]++;
	}
      }	// End of loop over simulations
    } // End case of W --> tau
  } // End of loop over events
  


  // Compute the predictions for pt, HT, and MHT from
  // hadronically decaying taus

  // Prediction of tau-jet pt
  // Loop over pt bins
  for(int ptBinIdx = 0; ptBinIdx < kBinsPtN; ++ptBinIdx) {
    // Compute mean and standard deviation
    // of number of events in this pt bin
    float mean  = 0.;
    float mean2 = 0.;
    // Loop over number of simulations
    for(int simIdx = 0; simIdx < nSimSteps; ++simIdx) {
      mean  += 1. * nPredTauJetPt2D[simIdx][ptBinIdx];
      mean2 += 1. * nPredTauJetPt2D[simIdx][ptBinIdx]*nPredTauJetPt2D[simIdx][ptBinIdx];
    }
    mean =  mean / nSimSteps;
    mean2 = mean2 / nSimSteps;
    float sig = sqrt( mean2 - mean*mean );
    hPredTauJetPt->SetBinContent(ptBinIdx+1,mean);
    hPredTauJetPt->SetBinError(ptBinIdx+1,sig);
  } // End of loop over pt bins

  // Prediction of tau-jet HT
  // Loop over HT bins
  for(int htBinIdx = 0; htBinIdx < kBinsHtN; ++htBinIdx) {
    // Compute mean and standard deviation
    // of number of events in this HT bin
    float mean  = 0.;
    float mean2 = 0.;
    // Loop over number of simulations
    for(int simIdx = 0; simIdx < nSimSteps; ++simIdx) {
      mean  += 1. * nPredHT2D[simIdx][htBinIdx];
      mean2 += 1. * nPredHT2D[simIdx][htBinIdx]*nPredHT2D[simIdx][htBinIdx];
    }
    mean =  mean / nSimSteps;
    mean2 = mean2 / nSimSteps;
    float sig = sqrt( mean2 - mean*mean );
    hPredHT->SetBinContent(htBinIdx+1,mean);
    hPredHT->SetBinError(htBinIdx+1,sig);
  } // End of loop over HT bins


  // Prediction of tau-jet MHT
  // Loop over MHT bins
  for(int mhtBinIdx = 0; mhtBinIdx < kBinMhtN; ++mhtBinIdx) {
    // Compute mean and standard deviation
    // of number of events in this MHT bin
    float mean  = 0.;
    float mean2 = 0.;
    // Loop over number of simulations
    for(int simIdx = 0; simIdx < nSimSteps; ++simIdx) {
      mean  += 1. * nPredMHT2D[simIdx][mhtBinIdx];
      mean2 += 1. * nPredMHT2D[simIdx][mhtBinIdx]*nPredMHT2D[simIdx][mhtBinIdx];
    }
    mean =  mean / nSimSteps;
    mean2 = mean2 / nSimSteps;
    float sig = sqrt( mean2 - mean*mean );
    hPredMHT->SetBinContent(mhtBinIdx+1,mean);
    hPredMHT->SetBinError(mhtBinIdx+1,sig);
  } // End of loop over MHT bins


  // --- Save the Histograms to File -----------------------------------
  TFile outFile("HadTau_WJetMC_PredGen.root","RECREATE");
  hTrueTauJetPt->Write();
  hTrueHT->Write();
  hTrueMHT->Write();
  hPredTauJetPt->Write();
  hPredHT->Write();
  hPredMHT->Write();
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


// Returns the index between 0 and nBins-1 of that bin
// x falls into given nBins equidistant bins between
// min and max. If x falls out of this range, -1 is
// returned
int bin(float x, int nBins, float min, float max) {
  int b = -1;
  if( x >= min && x <= max ) {
    float dx = (max - min)/nBins;
    for(int i = 0; i < nBins; ++i) {
      if( x <= min + (1+i)*dx ) {
	b = i;
	break;
      }
    }
  }

  return b;
}
