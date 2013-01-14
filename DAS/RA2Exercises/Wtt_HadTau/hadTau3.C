#include <algorithm>
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



// === Main Function ===================================================
void hadTau3(int nSimSteps,
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
  


  // --- Set Up the Tree -----------------------------------------------

  // Get the tree from file
  TChain* tr = new TChain("AnaTree");
  tr->Add("/nfs/dust/test/cmsdas/school61/susy/ntuple/2013-v1/WJets_0.root");

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


    // Select only events where the W decayed either
    // - into a muon (pdgId 13)
    // - into a tau (pdgId 15)
    if( flgW != 13 && flgW != 15 ) continue;
    // Select only events where the tau decayed hadronically
    if( flgW == 15 && flgTauHad != 1) continue;


    // Identify the charged generator-level lepton
    // from the W decay
    int chargedGenLepIdx = -1;
    if( !findMuOrTau(chargedGenLepIdx,genLepID,kGenLepColSize) ) continue; // Want exactly one lepton

    // Compute some more kinematic quantities of the 
    // generator-level lepton from the W decay
    aux4Vector.SetPxPyPzE(genLepPx[chargedGenLepIdx],genLepPy[chargedGenLepIdx],genLepPz[chargedGenLepIdx],genLepE[chargedGenLepIdx]);
    float genLepEta = aux4Vector.Eta();
    float genLepPhi = aux4Vector.Phi();
    float genLepPt  = aux4Vector.Pt();

    // Use only events where the lepton is inside the muon acceptance
    if( genLepPt < kMuonPtMin ) continue;
    if( TMath::Abs(genLepEta) > kMuonEtaMax ) continue;



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

      // If tau jet meets same criteria as other RA2 jets for HT,
      // recompute NJets and check if NJets >= 3
      float trueNJet = selNJet;
      if( tauJetPt > kHtJetPtMin && TMath::Abs(genLepEta) < kHtJetEtaMax ) trueNJet++;
      if( trueNJet < 3 ) continue;

      // Fill histogram of true tau-jet pt
      hTrueTauJetPt->Fill(tauJetPt);

      // If tau jet meets same criteria as other RA2 jets for HT,
      // add tau-jet pt to HT and fill histogram of true HT
      if( tauJetPt > kHtJetPtMin && TMath::Abs(genLepEta) < kHtJetEtaMax ) {
	hTrueHt->Fill(selHt+tauJetPt);
      }
      // If tau jet meets same criteria as other RA2 jets for MHT,
      // add tau-jet pt to MHT and fill histogram of true MHT
      if( tauJetPt > kMhtJetPtMin && TMath::Abs(genLepEta) < kMhtJetEtaMax ) {
	float trueMhtX = selMhtX - tauJetPt*cos(genLepPhi);
	float trueMhtY = selMhtY - tauJetPt*sin(genLepPhi);
	// Calculate the MHT from x and y component
	hTrueMht->Fill(sqrt( trueMhtX*trueMhtX + trueMhtY*trueMhtY ));
      }
    } else if( flgW == 13 ) {
      // In case the W decayed into a muon
      // - scale its pt by a random factor drawn from the
      //   tau-reponse template to simulate the tau measurement
      // - use simulated tau-pt to predict HT and MHT

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
	if( binSimTauJetPt > -1 ) nPredTauJetPt2D[binSimTauJetPt][it]++;
	  	
	// If simulated tau-jet meets same criteria as RA2 jets for HT,
	// add simulated tau-jet pt to HT and increment number of events
	// with that HT
	if( simTauJetPt > kHtJetPtMin && TMath::Abs(genLepEta) < kHtJetEtaMax ) {
	  float simHt = selHt + simTauJetPt;
	  int binSimHT = bin(simHt,kBinsHtN,kBinsHtMin,kBinsHtMax);
	  if( binSimHT > -1 ) nPredHt2D[binSimHT][it]++;
	}
	// If simulated tau-jet meets same criteria as RA2 jets for MHT,
	// add simulated tau-jet pt to MHT and increment number of events
	// with that MHT
	if( simTauJetPt > kMhtJetPtMin && TMath::Abs(genLepEta) < kMhtJetEtaMax ) {
	  float simMhtX = selMhtX - simTauJetPt*cos(genLepPhi);
	  float simMhtY = selMhtY - simTauJetPt*sin(genLepPhi);
	  float simMht = sqrt( simMhtX*simMhtX + simMhtY*simMhtY );
	  int binSimMHT = bin(simMht,kBinsMhtN,kBinsMhtMin,kBinsMhtMax);
	  if( binSimMHT > -1 ) nPredMht2D[binSimMHT][it]++;
	}
      }	// End of loop over simulations
    } // End case of W --> tau
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
  TFile outFile("HadTau_WJetMC_PredGen.root","RECREATE");
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
