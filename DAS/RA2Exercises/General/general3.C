#include <iostream>
#include <vector>

#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector2.h"



// === Global Variables ================================================

// Array dimensions in tree
const int kRecoJetColSize = 15;

// RA2 selection cuts
const float kHtJetPtMin   = 50.;
const float kHtJetEtaMax  = 2.5;
const float kMhtJetPtMin  = 30.;
const float kMhtJetEtaMax = 5.0;




// === Declaration of Auxiliary Functions ==============================
TString sampleLabel(int sampleId);
TString fileName(int sampleId);




// === Main Function ===================================================
void general3(int sampleId) {
  std::cout << "Analysing the " << sampleLabel(sampleId) << " sample" << std::endl;


  // --- Declare the Output Histograms ---------------------------------
  TH1* hNJets = new TH1F("hNJets",";N(jets);N(events)",12,0,12);
  hNJets->Sumw2();
  TH1* hHt = new TH1F("hHt",";H_{T} [GeV]",30,0,3000);
  hHt->Sumw2();
  hHt->GetXaxis()->SetNdivisions(505);
  TH1* hMht = new TH1F("hMht",";#slash{H}_{T} [GeV]",30,0,1500);
  hMht->Sumw2();
  hMht->GetXaxis()->SetNdivisions(505);
  TH1* hMEff = new TH1F("hMEff",";M_{eff} [GeV]",50,0,5000);
  hMEff->Sumw2();
  hMEff->GetXaxis()->SetNdivisions(505);
  std::vector<TH1*> hJetPt(6);
  std::vector<TH1*> hJetPhi(6);
  std::vector<TH1*> hJetEta(6);
  for(unsigned int i = 0; i < hJetEta.size(); ++i) {
    TString name = "hJetPt_";
    name += i;
    TString title = ";p_{T}(jet ";
    title += i+1;
    title += ") [GeV];N(events)";
    hJetPt.at(i) = new TH1F(name,title,30,0,1500);
    hJetPt.at(i)->Sumw2();

    name = "hJetPhi_";
    name += i;
    title = ";#phi(jet ";
    title += i+1;
    title += ");N(events)";
    hJetPhi.at(i) = new TH1F(name,title,24,-4,4);
    hJetPhi.at(i)->Sumw2();

    name = "hJetEta_";
    name += i;
    title = ";#eta(jet ";
    title += i+1;
    title += ");N(events)";
    hJetEta.at(i) = new TH1F(name,title,25,-5,5);
    hJetEta.at(i)->Sumw2();
  }



  // --- Declare the Variables Read from the Tree ----------------------
  // Reco-level jets
  int nRecoJets = 0;
  float recoJetPt[kRecoJetColSize];
  float recoJetPhi[kRecoJetColSize];
  float recoJetEta[kRecoJetColSize];

  // Number of reco-level muons and electrons
  int nRecoMus = 0;
  int nRecoEle = 0;

  // MC Event weight
  float evtWgt = 1.;



  // --- Set Up the Tree -----------------------------------------------

  // Get the tree from file
  TChain* tr = new TChain("AnaTree");
  tr->Add("/nfs/dust/test/cmsdas/school61/susy/ntuple/2013-v1/"+fileName(sampleId)+"_*.root");

  // Set the branches
  tr->SetBranchAddress("NrecoJet",&nRecoJets);
  tr->SetBranchAddress("recoJetPt",recoJetPt);
  tr->SetBranchAddress("recoJetPhi",recoJetPhi);
  tr->SetBranchAddress("recoJetEta",recoJetEta);
  tr->SetBranchAddress("NrecoMu",&nRecoMus);
  tr->SetBranchAddress("NrecoEle",&nRecoEle);
  if( sampleId > 0 ) tr->SetBranchAddress("EvtWgt",&evtWgt);



  // --- Process the Events in the Tree --------------------------------
  int nEvtsToProcess = tr->GetEntries();
  std::cout << "Processing " << nEvtsToProcess << " events" << std::endl;

  // Loop over the tree entries
  for(int evtIdx = 0; evtIdx < nEvtsToProcess; ++evtIdx) {
    if( evtIdx%100000 == 0 ) std::cout<<"  Event: " << evtIdx << std::endl;

    // Get the variables' values for this event
    tr->GetEntry(evtIdx);
    if( nRecoJets > kRecoJetColSize ) {
      std::cerr << "ERROR: more than " << kRecoJetColSize << " reco jets in event " << evtIdx << std::endl;
      exit(-1);
    }


    // Apply the lepton veto
    if( nRecoEle > 0 ) continue;
    if( nRecoMus > 0 ) continue;


    // Calculate RA2 selection-variables from jets
    float selNJet = 0; // Number of jets with pt > 50 GeV and |eta| < 2.5 (HT jets)
    float selHt   = 0.;
    float selMhtX = 0.;
    float selMhtY = 0.;
    // Loop over reco jets: they are ordered in pt
    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {
      // Calculate NJet and HT
      if( recoJetPt[jetIdx] > kHtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kHtJetEtaMax ) {
	selNJet++;
	selHt += recoJetPt[jetIdx];
      }
      // Calculate MHT components
      if( recoJetPt[jetIdx] > kMhtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kMhtJetEtaMax ) {
	selMhtX -= recoJetPt[jetIdx]*TMath::Cos(recoJetPhi[jetIdx]);
	selMhtY -= recoJetPt[jetIdx]*TMath::Sin(recoJetPhi[jetIdx]);
      }
    } // End of loop over reco jets
    float selMht = sqrt( selMhtX*selMhtX + selMhtY*selMhtY );
    

    // Select only events with at least 2 HT jets
    if( selNJet < 3 ) continue;



    //>>> PLACE OTHER RA2 CUTS HERE
    if( selHt < 500 ) continue;
    if( selMht < 200 ) continue;

    float phiMht = TMath::ATan2(selMhtY,selMhtX);
    int nMhtJets = 0;
    bool passesDeltaPhiCut = true;
    // Loop over reco jets: remember, they are ordered in pt!
    for(int jetIdx = 0; jetIdx < nRecoJets; ++jetIdx) {
      // Select MHT jets
      if( recoJetPt[jetIdx] > kMhtJetPtMin && TMath::Abs(recoJetEta[jetIdx]) < kMhtJetEtaMax ) {
	// Increase counter of MHT jets
	nMhtJets++;		

	// Compute deltaphi between -Pi and Pi
	float deltaPhi = TVector2::Phi_mpi_pi(recoJetPhi[jetIdx]-phiMht); 

	// Check DeltaPhi selection criterion
	if( nMhtJets == 1 || nMhtJets == 2 ) {
	  passesDeltaPhiCut = TMath::Abs(deltaPhi) > 0.5;
	} else {
	  passesDeltaPhiCut = TMath::Abs(deltaPhi) > 0.3;
	}
	if( !passesDeltaPhiCut ) break;
      }
      if( nMhtJets == 3 ) break; // DeltaPhi cut only for first three jets
    }
    if( !passesDeltaPhiCut ) continue;





    
    // Event weight in plots
    float weight = evtWgt;
    if( sampleId == 4 ) weight *= 6.26/5.274; // Correct Z-->inv xs from LO to NNLO
    if( sampleId == 5 ) weight *= 1.3; // 7 TeV k-factor

    
    // Fill histogram
    hNJets->Fill(selNJet,weight);
    hHt->Fill(selHt,weight);
    hMht->Fill(selMht,weight);
    hMEff->Fill(selHt+selMht,weight);
    for(int i = 0; i < static_cast<int>(hJetPt.size()); ++i) {
      if( i == nRecoJets ) break;
      hJetPt.at(i)->Fill(recoJetPt[i],weight);
      hJetPhi.at(i)->Fill(recoJetPhi[i],weight);
      hJetEta.at(i)->Fill(recoJetEta[i],weight);
    }
  } // End of loop over events



  // --- Save the Histograms to File -----------------------------------
  TFile outFile("General_"+fileName(sampleId)+".root","RECREATE");
  hNJets->Write();
  hHt->Write();
  hMht->Write();
  hMEff->Write();
  for(unsigned int i = 0; i < hJetPt.size(); ++i) {
    hJetPt[i]->Write();
    hJetEta[i]->Write();
    hJetPhi[i]->Write();
  }  
}




// === Implementation of Auxiliary Functions =====================

// Return the label for a given sample
TString sampleLabel(int sampleId) {
  TString label = "";
  if( sampleId == 0 )      label += "Data";
  else if( sampleId == 1 ) label += "QCD";
  else if( sampleId == 2 ) label += "t#bar{t}+Jets";
  else if( sampleId == 3 ) label += "W(l#nu)+Jets";
  else if( sampleId == 4 ) label += "Z(#nu#bar{#nu})+Jets";
  else if( sampleId == 5 ) label += "LM6";
  else {
    std::cerr << "ERROR: no sample with id " << sampleId << std::endl;
    exit(-1);
  }

  return label;
}


// Return the file name for a given sample
TString fileName(int sampleId) {
  TString name = "";
  if( sampleId == 0 )      name += "Data";
  else if( sampleId == 1 ) name += "QCD";
  else if( sampleId == 2 ) name += "TTJets";
  else if( sampleId == 3 ) name += "WJets";
  else if( sampleId == 4 ) name += "ZInv";
  else if( sampleId == 5 ) name += "LM6";
  else {
    std::cerr << "ERROR: no sample with id " << sampleId << std::endl;
    exit(-1);
  }

  return name;
}
