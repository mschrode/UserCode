//  Example script showing how to read
//  dijets from Kalibri ntuple and plot some quantities.
//
//  For questions on ROOT see also the manual at
//  http://root.cern.ch/drupal/content/users-guide
// ---------------------------------------------

#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"



// ---- Global variables -----------------------
// ---------------------------------------------
// Auxiliary quantities
const int NJETS  = 100;   // Maximum number of jets per event

// Histograms
TH1* hJet1Pt    = 0;   // Pt of first (highest corrected pt) jet



// ---- Function declarations ------------------
// ---------------------------------------------
void bookHistos();
void draw();
bool readDijets(int nEvents);
void runDijets(int nEvts);



// ---- Function implementations ---------------
// ---------------------------------------------

// Run example script
// This is the actual command to be typed
// used in the ROOT session.
// ---------------------------------------------
void runDijets(int nEvts) {
  bookHistos();
  readDijets(nEvts);
  draw();
}


// Booking histograms
// Old histograms are deleted
// ---------------------------------------------
void bookHistos() {
  std::cout << "Creating histograms" << std::endl;

  // Distribution of pt of both jets
  if( hJet1Pt ) delete hJet1Pt;                          // Delete potential old histogram
  hJet1Pt = new TH1D("hJet1Pt","Jet1 p_{T}",100,0,1000); // Create histogram
  hJet1Pt->SetXTitle("p_{T,1} (GeV)");                   // Set title of x-axis
  hJet1Pt->SetYTitle("Number of events");                // Set title of y-axis
}



// Read dijets from file and fill histograms. 
// Param:
//  input    Name of input file
//  nEvents  Maximum number of events to be read
// Return:
//  true for success and false for failure.
// ---------------------------------------------
bool readDijets(int nEvents) {
  // ---- Reading tree from file ---------------

  // Create TChain
  TChain* chain = new TChain("DiJetTree");
  // Add files with Trees (hard-coded for the time being)
  chain->Add("~/lustre/data/Jet_2010B_Nov4ReReco_148822-149294/Jetjob_34_Dijet-ak5Calo.root");
  chain->Add("~/lustre/data/Jet_2010B_Nov4ReReco_148822-149294/Jetjob_35_Dijet-ak5Calo.root");
  chain->Add("~/lustre/data/Jet_2010B_Nov4ReReco_148822-149294/Jetjob_36_Dijet-ak5Calo.root");


  // ---- Reading variables from tree ----------
  // ---- and filling histograms ---------------

  // See also the chapter 'Trees' in the ROOT Manual

  // Book variables to be read from tree
  int nobjJet;                   // Number of jets in an event
  float jetPt[NJETS];            // Pt of jets in an event
  float jetCorrL2L3[NJETS];	 // Jet energy correction factors
  float jetPhi[NJETS];           // Phi of jets in an event
  float jetEta[NJETS];		 // Jet eta

  // Set branch addresses
  chain->SetBranchAddress("NobjJet",&nobjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetPhi",jetPhi);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);

  std::cout << "Processing events" << std::endl;

  // If tree contains less entries than nEvents
  // use less events
  int nMax = nEvents;
  if( chain->GetEntries() < nEvents ) nMax = chain->GetEntries();

  // Loop over nMax entries, read variables
  // and fill histogram
  for(int n = 0; n < nMax; n++) {
    if( n%1000 == 0 ) std::cout << "Event " << n << std::endl;
    
    // Get this event i.e. entry n
    // and fill values into variables
    chain->GetEntry(n);
    
    // Check if array is large enough
    if( nobjJet > NJETS ) {
      std::cerr << "\nERROR: 'nobjJet = " << nobjJet << " > NJETS = " << NJETS << "'" << std::endl;
      continue;
    }

    // Select events with at least 2 jets
    if( nobjJet < 2 ) continue;
    
    // Select events in central detector region
    if( std::abs(jetEta[0]) > 1.1 || std::abs(jetEta[1]) > 1.1 ) continue;
    
    // Select events with back-to-back jets (dijets)
    double deltaPhi = TVector2::Phi_mpi_pi(jetPhi[0] - jetPhi[1]);
    if( std::abs(deltaPhi) < 2.7 ) continue;
    
    // Apply jet energy corrections
    double ptCorr = jetCorrL2L3[0]*jetPt[0];
    
    // Fill pt into histogram
    hJet1Pt->Fill(ptCorr);
  } // End of loop over entries

  delete chain;
    
  return true;
}



// Draw histograms
// ---------------------------------------------
void draw() {
  std::cout << "Plotting histograms" << std::endl;

  // Create a canvas, on which the histogram of
  // jet pt is drawn
  TCanvas * cJetPt = new TCanvas("cJetPt","Jet Pt",0,0,500,500);
  
  // Go onto that canvas
  cJetPt->cd();
  
  // Draw histogram onto that canvas
  hJet1Pt->Draw();

  // Set logarithmic scale on y-axis
  cJetPt->SetLogy();
}
