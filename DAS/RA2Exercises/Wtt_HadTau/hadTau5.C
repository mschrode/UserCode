#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"



// === Global Variables ================================================

// Array dimensions in tree
const int kRecoJetColSize = 15;
const int kRecoMuColSize = 15;
const int kGenLepColSize = 6;

// RA2 selection cuts
const float kHtJetPtMin   = 50.;
const float kHtJetEtaMax  = 2.5;
const float kMhtJetPtMin  = 30.;
const float kMhtJetEtaMax = 5.0;

// Tau-template and muon related
const int kNRespPtBins  = 4;
const float kMuonPtMin  = 20.;
const float kMuonEtaMax = 2.1;

// Binning for prediction
const int kBinsPtN      = 50;
const float kBinsPtMin  = 0.;
const float kBinsPtMax  = 500.;
const int kBinsHtN      = 25;
const float kBinsHtMin  = 0.;
const float kBinsHtMax  = 2500.;
const int kBinsMhtN     = 20;
const float kBinsMhtMin = 0.;
const float kBinsMhtMax = 1000.;



// === Declaration of Auxiliary Functions ==============================
bool findMuOrTau(int &lepIdx, const float* lepID, int nLep);
bool findMatchedJet(int &matchedJetIdx, const float* jetEta, const float* jetPhi, int nJets, float lepEta, float lepPhi, float deltaRMax);
float deltaR(float eta1, float eta2, float phi1, float phi2);
int tauPtBin(float pt);
int bin(float x, int nBins, float min, float max);
std::vector<TH1*> initializeHistsForNPred(const TH1* hPredDist, const std::vector< std::vector<int> > &nPred2D, const TString &histName, const TString &varName);
void getPrediction(std::vector<TH1*> &hNPredPerBin, TH1* &hPredDist, const std::vector< std::vector<int> > &nPred2D);
bool passesDeltaPhiCut(float mhtPhi, int nJets, const float* jetPt, const float* jetEta, const float* jetPhi, int tauJetMatchIdx, float tauJetPt, float tauJetEta, float tauJetPhi);



// === Main Function ===================================================
void hadTau5(bool isMC,
	     int nSimSteps,
	     const TString &respTempl = "HadTau_WJetMC.root",
	     int nEvts = -1) {

  // Counters: store for each simulation step, how many events
  // are predicted in each pt / HT / MHT bin
  // Dimension is [nVarBins][nSimSteps]
  std::vector< std::vector<int> > nPredTauJetPt2D(kBinsPtN);
  std::vector< std::vector<int> > nPredHt2D(kBinsHtN);
  std::vector< std::vector<int> > nPredMht2D(kBinsMhtN);
  for(int i = 0; i < kBinsPtN; ++i) {
    nPredTauJetPt2D[i] = std::vector<int>(nSimSteps,0);
  }
  for(int i = 0; i < kBinsHtN; ++i) {
    nPredHt2D[i] = std::vector<int>(nSimSteps,0);
  }
  for(int i = 0; i < kBinsMhtN; ++i) {
    nPredMht2D[i] = std::vector<int>(nSimSteps,0);
  }



  // --- Declare the Output Histograms ---------------------------------
  TH1* hMuonPt = new TH1F("hMuonPt",";p_{T}(#mu) [GeV];N(events)",kBinsPtN,kBinsPtMin,kBinsPtMax);
  hMuonPt->Sumw2();
  TH1* hPredTauJetPt = new TH1F("hPredTauJetPt",";Predicted p_{T}(#tau) [GeV];N(#tau)",kBinsPtN,kBinsPtMin,kBinsPtMax);   
  TH1* hPredHt = new TH1F("hPredHt",";Predicted H_{T} [GeV];N(events)",kBinsHtN,kBinsHtMin,kBinsHtMax);   
  TH1* hPredMht = new TH1F("hPredMht",";Predicted #slash{H}_{T} [GeV];N(events)",kBinsMhtN,kBinsMhtMin,kBinsMhtMax);   
  TH1* hTrueTauJetPt = new TH1F("hTrueTauJetPt",";True p_{T}(#tau) [GeV];N(#tau)",kBinsPtN,kBinsPtMin,kBinsPtMax);   
  TH1* hTrueHt = new TH1F("hTrueHt",";True H_{T} [GeV];N(events)",kBinsHtN,kBinsHtMin,kBinsHtMax);   
  TH1* hTrueMht = new TH1F("hTrueMht",";True #slash{H}_{T} [GeV];N(events)",kBinsMhtN,kBinsMhtMin,kBinsMhtMax);   



  // --- Declare the Variables Read from the Tree ----------------------
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
  
  // Reco-level muons
  int nRecoMus = 0;
  float recoMuPt[kRecoMuColSize];
  float recoMuPhi[kRecoMuColSize];
  float recoMuEta[kRecoMuColSize];

  // Reco-level electrons
  int nRecoEle = 0;



  // --- Set Up the Tree -----------------------------------------------

  // Get the tree from file
  TChain* tr = new TChain("AnaTree");
  tr->Add("/nfs/dust/test/cmsdas/school61/susy/ntuple/2013-v1/WJets_*.root");

  // Set the branches
  tr->SetBranchAddress("NrecoJet",&nRecoJets);
  tr->SetBranchAddress("recoJetPt",recoJetPt);
  tr->SetBranchAddress("recoJetPhi",recoJetPhi);
  tr->SetBranchAddress("recoJetEta",recoJetEta);
  tr->SetBranchAddress("NrecoMu",&nRecoMus);
  tr->SetBranchAddress("recoMuPt",recoMuPt);
  tr->SetBranchAddress("recoMuPhi",recoMuPhi);
  tr->SetBranchAddress("recoMuEta",recoMuEta);
  tr->SetBranchAddress("NrecoEle",&nRecoEle);
  if( isMC ) {			// Generator quantities only in MC...
    tr->SetBranchAddress("flgW",&flgW);
    tr->SetBranchAddress("flgTauHad",&flgTauHad);
    tr->SetBranchAddress("lepPx",genLepPx);
    tr->SetBranchAddress("lepPy",genLepPy);
    tr->SetBranchAddress("lepPz",genLepPz);
    tr->SetBranchAddress("lepE",genLepE);
    tr->SetBranchAddress("lepID",genLepID);
  }



  // --- Get the Tau Templates from File -------------------------------
  TFile fTauResp(respTempl,"READ");
  std::vector<TH1*> hTauResp(kNRespPtBins);
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
    if( evtIdx%100000 == 0 ) std::cout<<"  Event: " << evtIdx << std::endl;

    // Get the variables' values for this event
    tr->GetEntry(evtIdx);
    if( nRecoJets > kRecoJetColSize ) {
      std::cerr << "ERROR: more than " << kRecoJetColSize << " reco jets in event " << evtIdx << std::endl;
      exit(-1);
    }


    // Select events with exactly one well-isolated reco muon
    // - scale its pt by a random factor drawn from the
    //   tau-reponse template to simulate the tau measurement
    // - use simulated tau-pt to predict HT and MHT
    if( nRecoMus == 1 && nRecoEle == 0 ) {
      float muPt = recoMuPt[0];
      float muEta = recoMuEta[0];
      float muPhi = recoMuPhi[0];

      // Use only events where the muon is inside acceptance
      if( muPt < kMuonPtMin ) continue;
      if( muEta > kMuonEtaMax ) continue;


      // "Cross cleaning": find the jet that corresponds to
      // the muon. Associate jet and lepton if they are closer
      // in eta-phi space than deltaRMax
      int muJetIdx = -1;
      float deltaRMax = 0.1;
      if( !findMatchedJet(muJetIdx,recoJetEta,recoJetPhi,nRecoJets,muEta,muPhi,deltaRMax) ) continue; // Want unambigious matching (exactly one jet within deltaRMax)


      // Calculate RA2 selection-variables from "cleaned" jets
      float selNJet = 0; // Number of jets with pt > 50 GeV and |eta| < 2.5 (HT jets)
      float selHt   = 0.;
      float selMhtX = 0.;
      float selMhtY = 0.;
      for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
	// Skip this jet if it is the muon
	if( jetIdx == muJetIdx ) continue;
	
	// Calculate NJet and HT
	if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) {
	  selNJet++;
	  selHt += recoJetPt[jetIdx];
	}
	// Calculate MHT
	if( recoJetPt[jetIdx] > kMhtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kMhtJetEtaMax ) {
	  selMhtX -= recoJetPt[jetIdx]*TMath::Cos(recoJetPhi[jetIdx]);
	  selMhtY -= recoJetPt[jetIdx]*TMath::Sin(recoJetPhi[jetIdx]);
	}
      } // End of loop over reco jets
      

      // Select only events with at least 2 HT jets
      if( selNJet < 2 ) continue;

      // Find the pt bin of the tau-response template
      int tauRespPtBin = tauPtBin(muPt);
      if( tauRespPtBin < 0 || tauRespPtBin >= kNRespPtBins ) continue;

      // Plot the muon pt as control plot
      hMuonPt->Fill(muPt);
      
      // Perform nSimSteps simulations
      for(int it = 0; it < nSimSteps; ++it) {
	// Get random number from tau-response template
	float scale = hTauResp[tauRespPtBin]->GetRandom();
	// Scale muon pt with tau response --> simulate tau jet pt
	float simTauJetPt = scale * muPt;

	// If simulted tau-jet meets same criteria as other RA2 jets for HT,
	// recompute NJets and check if NJets >= 3
	float simNJet = selNJet;
	if( simTauJetPt > kHtJetPtMin && TMath::Abs(muEta) < kHtJetEtaMax ) simNJet++;
	if( simNJet < 3 ) continue;

	float simHt = selHt;
	float simMhtX = selMhtX;
	float simMhtY = selMhtY;
	  	
	// If simulated tau-jet meets same criteria as RA2 jets for HT,
	// add simulated tau-jet pt to HT and increment number of events
	// with that HT
	if( simTauJetPt > kHtJetPtMin && TMath::Abs(muEta) < kHtJetEtaMax ) {
	  simHt += simTauJetPt;
	}
	// If simulated tau-jet meets same criteria as RA2 jets for MHT,
	// add simulated tau-jet pt to MHT
	if( simTauJetPt > kMhtJetPtMin && TMath::Abs(muEta) < kMhtJetEtaMax ) {
	  simMhtX -= simTauJetPt*cos(muPhi);
	  simMhtY -= simTauJetPt*sin(muPhi);
	}
	float simMht = sqrt( simMhtX*simMhtX + simMhtY*simMhtY );

	if( simHt < 500 ) continue;
	if( simMht < 200 ) continue;
	float simMhtPhi = std::atan2(simMhtY,simMhtX);
	if( !passesDeltaPhiCut(simMhtPhi,nRecoJets,recoJetPt,recoJetEta,recoJetPhi,muJetIdx,simTauJetPt,muEta,muPhi) ) continue;

	    
	// Increment number of events with the simulated tau-jet pt, HT, MHT
	int binSimTauJetPt = bin(simTauJetPt,kBinsPtN,kBinsPtMin,kBinsPtMax);
	if( binSimTauJetPt > -1 ) nPredTauJetPt2D[binSimTauJetPt][it]++;
	int binSimHT = bin(simHt,kBinsHtN,kBinsHtMin,kBinsHtMax);
	if( binSimHT > -1 ) nPredHt2D[binSimHT][it]++;
	int binSimMHT = bin(simMht,kBinsMhtN,kBinsMhtMin,kBinsMhtMax);
	if( binSimMHT > -1 ) nPredMht2D[binSimMHT][it]++;

      }	// End of loop over simulations
    } // End if exactly one muon


    // In case of MC, plot the true distributions for comparison
    if( isMC ) {
      // Select only events where the W decayed into a hadronically
      // decaying tau
      if( !(flgW == 15 && flgTauHad == 1) ) continue;


      // Identify the generator-level tau from the W decay
      int genTauIdx = -1;
      if( !findMuOrTau(genTauIdx,genLepID,kGenLepColSize) ) continue; // Want exactly one lepton
    
      // Compute some more kinematic quantities of the tau
      // generator-level lepton from the W decay
      aux4Vector.SetPxPyPzE(genLepPx[genTauIdx],genLepPy[genTauIdx],genLepPz[genTauIdx],genLepE[genTauIdx]);
      float genTauEta = aux4Vector.Eta();
      float genTauPhi = aux4Vector.Phi();
      float genTauPt  = aux4Vector.Pt();

      // Use only events where the lepton is inside the muon acceptance
      if( genTauPt < kMuonPtMin ) continue;
      if( TMath::Abs(genTauEta) > kMuonEtaMax ) continue;


      // "Cross cleaning": find the jet that originates in the
      // hadronic-tau decay. Associate jet and tau if they are
      // closer in eta-phi space than deltaRMax
      int tauJetIdx = -1;
      float deltaRMax = 0.1;
      if( genTauPt < 50. ) deltaRMax = 0.2;
      if( !findMatchedJet(tauJetIdx,recoJetEta,recoJetPhi,nRecoJets,genTauEta,genTauPhi,deltaRMax) ) continue; // Want unambigious matching (exactly one jet within deltaRMax)


      // Calculate RA2 selection-variables from "cleaned" jets
      float selNJet = 0; // Number of jets with pt > 50 GeV and |eta| < 2.5 (HT jets)
      float selHt   = 0.;
      float selMhtX = 0.;
      float selMhtY = 0.;
      for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
	// Skip this jet if it is the tau
	if( jetIdx == tauJetIdx ) continue;
	
	// Calculate NJet and HT
	if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) {
	  selNJet++;
	  selHt += recoJetPt[jetIdx];
	}
	// Calculate MHT
	if( recoJetPt[jetIdx] > kMhtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kMhtJetEtaMax ) {
	  selMhtX -= recoJetPt[jetIdx]*TMath::Cos(recoJetPhi[jetIdx]);
	  selMhtY -= recoJetPt[jetIdx]*TMath::Sin(recoJetPhi[jetIdx]);
	}
      } // End of loop over reco jets

      // Select only events with at least 2 HT jets
      if( selNJet < 2 ) continue;

      float tauJetPt  = recoJetPt[tauJetIdx];

      // If tau jet meets same criteria as other RA2 jets for HT,
      // recompute NJets and check if NJets >= 3
      float trueNJet = selNJet;
      if( tauJetPt > kHtJetPtMin && TMath::Abs(genTauEta) < kHtJetEtaMax ) trueNJet++;
      if( trueNJet < 3 ) continue;

      float trueHt = selHt;
      float trueMhtX = selMhtX;
      float trueMhtY = selMhtY;

      // If tau jet meets same criteria as other RA2 jets for HT,
      // add tau-jet pt to HT
      if( tauJetPt > kHtJetPtMin && TMath::Abs(genTauEta) < kHtJetEtaMax ) {
	trueHt += tauJetPt;
      }
      // If tau jet meets same criteria as other RA2 jets for MHT,
      // add tau-jet pt to MHT
      if( tauJetPt > kMhtJetPtMin && TMath::Abs(genTauEta) < kMhtJetEtaMax ) {
	trueMhtX -= tauJetPt*cos(genTauPhi);
	trueMhtY -= tauJetPt*sin(genTauPhi);
      }
      float trueMht = sqrt( trueMhtX*trueMhtX + trueMhtY*trueMhtY );
      if( trueHt < 500 ) continue;
      if( trueMht < 200 ) continue;

      float trueMhtPhi = std::atan2(trueMhtY,trueMhtX);
      if( !passesDeltaPhiCut(trueMhtPhi,nRecoJets,recoJetPt,recoJetEta,recoJetPhi,tauJetIdx,tauJetPt,genTauEta,genTauPhi) ) continue;

      hTrueTauJetPt->Fill(tauJetPt);
      hTrueHt->Fill(trueHt);
      hTrueMht->Fill(trueMht);
    } // End of isMC
  } // End of loop over events
  


  // --- Compute Prediction --------------------------------------------

  // Prepare plots to store number of predicted events
  // per simulation step
  std::vector<TH1*> hNPredTauJetPt = initializeHistsForNPred(hPredTauJetPt,nPredTauJetPt2D,"hNPredTauJetPt","p_{T}(#tau)");
  std::vector<TH1*> hNPredHt = initializeHistsForNPred(hPredHt,nPredHt2D,"hNPredHt","H_{T}");
  std::vector<TH1*> hNPredMht = initializeHistsForNPred(hPredMht,nPredMht2D,"hNPredMht","#slash{H}_{T}");

  // Compute the predictions for pt, HT, and MHT from
  // hadronically decaying taus
  getPrediction(hNPredTauJetPt,hPredTauJetPt,nPredTauJetPt2D);
  getPrediction(hNPredHt,hPredHt,nPredHt2D);
  getPrediction(hNPredMht,hPredMht,nPredMht2D);


  // --- Save the Histograms to File -----------------------------------
  TString outFileName = "HadTau_";
  if( isMC ) outFileName += "WJetMC_PredReco.root";
  else       outFileName += "Data_Pred.root";
  TFile outFile(outFileName,"RECREATE");
  hMuonPt->Write();
  hTrueTauJetPt->Write();
  hTrueHt->Write();
  hTrueMht->Write();
  hPredTauJetPt->Write();
  hPredHt->Write();
  hPredMht->Write();
  for(std::vector<TH1*>::const_iterator it = hNPredTauJetPt.begin();
      it != hNPredTauJetPt.end(); ++it) {
    (*it)->Write();
  }
  for(std::vector<TH1*>::const_iterator it = hNPredHt.begin();
      it != hNPredHt.end(); ++it) {
    (*it)->Write();
  }
  for(std::vector<TH1*>::const_iterator it = hNPredMht.begin();
      it != hNPredMht.end(); ++it) {
    (*it)->Write();
  }
  outFile.Close();
}



// === Implementation of Auxiliary Functions =====================

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


std::vector<TH1*> initializeHistsForNPred(const TH1* hPredDist, const std::vector< std::vector<int> > &nPred2D, const TString &histName, const TString &varName) {
  std::vector<TH1*> hists;
  
  // Loop over bins e.g. in pt or HT
  for(unsigned int varBin = 0; varBin < nPred2D.size(); ++varBin) {
    // Prepare names and titles
    TString name = histName+"_";
    name += varBin;
    char xTitle[200];
    sprintf(xTitle,"N(pred) with %.0f < %s < %.0f GeV",hPredDist->GetXaxis()->GetBinLowEdge(varBin+1),varName.Data(),hPredDist->GetXaxis()->GetBinUpEdge(varBin+1));
    
    // Find maximum prediction to set histogram range
    int max = *(std::max_element(nPred2D[varBin].begin(),nPred2D[varBin].end()));
		
    // Create histogram
    TH1* h = new TH1F(name,"",max+2,-0.5,max+1.5);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle("Entries");
    h->GetXaxis()->SetNdivisions(505);

    // Store histogram in vector
    hists.push_back(h);
  } // End of loop over bins

  return hists;
}


void getPrediction(std::vector<TH1*> &hNPredPerBin, TH1* &hPredDist, const std::vector< std::vector<int> > &nPred2D) {
  // Loop over bins e.g. in pt or HT
  for(unsigned int varBin = 0; varBin < nPred2D.size(); ++varBin) {
    // Store number of predictions for this bin
    // per simulation step 
    for(unsigned int simIdx = 0; simIdx < nPred2D[varBin].size(); ++simIdx) {
      hNPredPerBin[varBin]->Fill(nPred2D[varBin][simIdx]);
    }

    // Compute mean and standard deviation
    // of number of events in this bin
    float mean  = hNPredPerBin[varBin]->GetMean();
    float sig   = hNPredPerBin[varBin]->GetRMS();

    // Fill prediction into histogram
    hPredDist->SetBinContent(varBin+1,mean);
    hPredDist->SetBinError(varBin+1,sig);
  } // End of loop over pt bins
}



bool passesDeltaPhiCut(float mhtPhi, int nJets, const float* jetPt, const float* jetEta, const float* jetPhi, int tauJetMatchIdx, float tauJetPt, float tauJetEta, float tauJetPhi) {
  // Find leading three MHT jets, i.e. with pt > 30 GeV,
  //that are not matched to the tau jet
  int nMhtJets = 0;
  float mhtJetPt[3] = { 0., 0., 0. };
  float mhtJetPhi[3] = { 0., 0., 0. };
  for(int jetIdx = 0; jetIdx < nJets; ++jetIdx) {
    // Skip this jet if it is the tau
    if( jetIdx == tauJetMatchIdx ) continue;
    
    if( jetPt[jetIdx] > kMhtJetPtMin && TMath::Abs(jetEta[jetIdx]) < kMhtJetEtaMax ) {
      nMhtJets++;
      if( nMhtJets == 1 ) {
	mhtJetPt[0] = jetPt[jetIdx];
	mhtJetPhi[0] = jetPhi[jetIdx];
      } else if( nMhtJets == 2 ) {
	mhtJetPt[1] = jetPt[jetIdx];
	mhtJetPhi[1] = jetPhi[jetIdx];
      } else if( nMhtJets == 3 ) {
	mhtJetPt[2] = jetPt[jetIdx];
	mhtJetPhi[2] = jetPhi[jetIdx];
	break;
      }
    }
  }

  if( nMhtJets < 2 ) {		// One or no MHT jet
    return false;
  }

  // Check if tau jet counts as MHT jet
  // and in that case reorder the jets
  if( tauJetPt > kMhtJetPtMin && TMath::Abs(tauJetEta) < kMhtJetEtaMax ) {
    if( nMhtJets == 2 ) {	// Two MHT jets
      if( tauJetPt > mhtJetPt[0] ) {
	mhtJetPhi[2] = mhtJetPhi[1];
	mhtJetPhi[1] = mhtJetPhi[0];
	mhtJetPhi[0] = tauJetPhi;
      } else if( tauJetPt > mhtJetPt[1] ) {
	mhtJetPhi[2] = mhtJetPhi[1];
	mhtJetPhi[1] = tauJetPhi;
      } else {
	mhtJetPhi[2] = tauJetPhi;
      }
      nMhtJets++;
    } else {			// Three or more MHT jets
      if( tauJetPt > mhtJetPt[0] ) {
	mhtJetPhi[2] = mhtJetPhi[1];
	mhtJetPhi[1] = mhtJetPhi[0];
	mhtJetPhi[0] = tauJetPhi;
      } else if( tauJetPt > mhtJetPt[1] ) {
	mhtJetPhi[2] = mhtJetPhi[1];
	mhtJetPhi[1] = tauJetPhi;
      } else {
	mhtJetPhi[2] = tauJetPhi;
      }
    }
  }

  float dphi[3];
  for(int i = 0; i < nMhtJets; ++i) {
    dphi[i] = std::abs(TVector2::Phi_mpi_pi(mhtJetPhi[i]-mhtPhi));    
  }

  // Apply DeltaPhi cuts
  if( dphi[0] < 0.5 ) return false;
  if( dphi[1] < 0.5 ) return false;
  if( dphi[2] < 0.3 ) return false;
  return true;
}
