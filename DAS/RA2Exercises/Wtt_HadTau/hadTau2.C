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



// === Declaration of Auxiliary Functions ==============================
bool findMuOrTau(int &lepIdx, const float* lepID, int nLep);
bool findMatchedJet(int &matchedJetIdx, const float* jetEta, const float* jetPhi, int nJets, float lepEta, float lepPhi, float deltaRMax);
float deltaR(float eta1, float eta2, float phi1, float phi2);
int tauPtBin(float pt);



// === Main Function ===================================================
void hadTau2(const TString &inputEvents = "inputEvents", int nEvts = -1) {
  // --- Declare the Output Histograms ---------------------------------
  std::vector<TH1*> hTauResp(kNRespPtBins);
  for(int i = 0; i < kNRespPtBins; ++i) {
    TString name = "hTauResp_";
    name += i;
    hTauResp[i] = new TH1F(name,";p_{T}(visible) / p_{T}(generated);N(#tau)",50,0.,2.5);
  }


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


    // Select only events where the W decayed into a hadronically
    // decaying tau
    if( !(flgW == 15 && flgTauHad == 1) ) continue;

    // Identify the charged generator-level taus
    // from the W decay
    int genTauIdx = -1;
    if( !findMuOrTau(genTauIdx,genLepID,kGenLepColSize) ) continue; // Want exactly one lepton
    if( TMath::Abs(genLepID[genTauIdx]) != 15 ) continue; // Should be a tau!

    
    // Compute some more kinematic quantities of the tau
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
    if( flgW == 15 && genTauPt < 50. ) deltaRMax = 0.2;
    if( !findMatchedJet(tauJetIdx,recoJetEta,recoJetPhi,nRecoJets,genTauEta,genTauPhi,deltaRMax) ) continue; // Want unambigious matching (exactly one jet within deltaRMax)


    // Calculate RA2 selection-variables from "cleaned" jets
    float selNJet = 0; // Number of jets with pt > 50 GeV and |eta| < 2.5 (HT jets)

    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
      // Skip this jet if it is the lepton
      if( jetIdx == tauJetIdx ) continue;

      // Calculate NJet 
      if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) selNJet++;
    } // End of loop over reco jets

    
    // Select only events with at least 2 HT jets
    if( selNJet < 2 ) continue;
    

    // Fill histogram with relative visible energy of the tau
    // ("tau response template") for hadronically decaying tau
    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {	// Loop over reco jets
      // Select tau jet
      if( jetIdx == tauJetIdx ) {
	int ptBin = tauPtBin(genTauPt);
	if( ptBin >= 0 && ptBin < kNRespPtBins ) {
	  hTauResp[ptBin]->Fill( recoJetPt[jetIdx]/genTauPt );
	}
	break;		// End the jet loop once the tau jet has been found
      }
    }	// End of loop over reco jets
  } // End of loop over tree entries


  // --- Save the Histograms to File -----------------------------------
  TFile outFile("HadTau_WJetMC.root","UPDATE");
  for(int i = 0; i < kNRespPtBins; ++i) {
    hTauResp[i]->Write();
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
