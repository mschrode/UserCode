// $Id: computeTriggerEfficiencies.C,v 1.3 2011/05/09 18:37:25 mschrode Exp $

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


// Trigger efficiency plots from Kalibri ntuple
//
// Note: traditionally you would compute the efficiency of a
// trigger A using a reference trigger B as
//
//   eff = N(passing A && passing B) / N(passing B) .
//
// However, due to the sometimes large pre-scales there might
// only be very few events available for this method preventing
// a reliable measurement. Instead, the efficiency is defined
// as
//
//   eff = N(passing A) / N(passing B) .
//
// Of course, now 0 <= eff <= 1 is not true anymore but without
// knowledge of the pre-scales the previous definition is not
// exact either. The pt of the turn-on can still be correctly
// determined as the pt where eff is 99% of the plateau value
// of eff.
//
// TODO: Re-write to run over all triggers at the same time and 
// divide the different histograms.



// --------------------------------------------------
double eff(double *x, double *par) {
  return 0.5*par[2]*(erf(par[0]*(x[0]-par[1]))+1);
}



//! Create TChain from input root files. The root
//! files are expected to contain a TTree "DiJetTree".
//! There are two possible input options:
//!
//! 1) 'fileName' specifies a single root file; it ends
//!    with '.root';
//! 2) 'fileName' contains a list of root file names.
// --------------------------------------------------
TChain *createTChain(const TString &fileName) {
  TChain* chain = new TChain("DiJetTree"); 

  // Option 1: single root file
  if( fileName.EndsWith(".root") ) {
    chain->Add(fileName);
  }
  // Option 2: list of root files
  else {
    std::ifstream filelist;
    filelist.open(fileName);
    int nOpenedFiles = 0;
    if( filelist.is_open() ) {
      TString name = "";
      while( !filelist.eof() ) {
	filelist >> name;
	if( filelist.eof() ) break;
	chain->Add(name);
	nOpenedFiles++;
      }
    } else {
      std::cerr << "ERROR opening file '" << fileName << "'\n";
      exit(1);
    }
    filelist.close();
  }

  return chain;
}



// --------------------------------------------------
void computeTriggerEfficiences(const TString &fileNames, int runStart, int runEnd, double lumi, const TString &tagTrigger, const TString &probeTrigger, unsigned int nBins, double minPtAve, double maxPtAve, int nEvts = -1) {
  util::StyleSettings::presentationNoTitle();

  TString jetType = "";
  TString jetLabel = "";
  if( fileNames.Contains("AK5") || fileNames.Contains("ak5") ) {
    jetType += "AK5";
    jetLabel += "AK5";
  }
  if( fileNames.Contains("PF") || fileNames.Contains("pf") ) {
    jetType += "PF";
    jetLabel += " PF";
  } else if( fileNames.Contains("Calo") || fileNames.Contains("calo") ) {
    jetType += "Calo";
    jetLabel += " Calo";
  }
  jetLabel += " Jets";


  std::cout << "Creating histograms...\n";
  // Event counts
  TH1 *hTotal = util::HistOps::createTH1D("hTotal",nBins,minPtAve,maxPtAve,"p^{ave}_{T}","GeV","events");
  hTotal->Sumw2();
  TH1 *hPass = util::HistOps::createTH1D("hPass",nBins,minPtAve,maxPtAve,"p^{ave}_{T}","GeV","events");
  hPass->Sumw2();


  std::cout << "Adding files to chain...\n";
  // Set up chain
  TChain *chain = createTChain(fileNames);

  // Set branch addresses
  const int maxNJet = 50;
  
  int nTotal = 0;
  int nPass  =0;

  std::vector<double> corrJetPt(20);

  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetCorrL1[maxNJet];
  float jetCorrL2[maxNJet];
  float jetCorrL3[maxNJet];
  bool hltTag = false;
  bool hltProbe = false;

  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress(tagTrigger,&hltTag);
  chain->SetBranchAddress(probeTrigger,&hltProbe);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetCorrL1",jetCorrL1);
  chain->SetBranchAddress("JetCorrL2",jetCorrL2);
  chain->SetBranchAddress("JetCorrL3",jetCorrL3);

  // Loop over tree entries and fill histograms
  if( nEvts < 0 || nEvts > chain->GetEntries() ) nEvts = chain->GetEntries();
  std::cout << "Reading " << nEvts << " entries...\n";
  for(int n = 0; n < nEvts; ++n) {
    if( n%50000 == 0 ) std::cout << "  " << n << std::endl;
    chain->GetEntry(n);

    if( nObjJet > maxNJet ) {
      std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event.\n";
      continue;
    }

    if( nObjJet > 1 && (hltTag || hltProbe) ) {
      ++nTotal;

      // Sort jets by corrected pt
      for(int i = 0; i < static_cast<int>(corrJetPt.size()); ++i) {
	if( i < nObjJet ) corrJetPt[i] = jetCorrL1[i]*jetCorrL2[i]*jetCorrL3[i]*jetPt[i];
	else corrJetPt[i] = 0.;      
      }
      std::sort(corrJetPt.begin(),corrJetPt.end());

      // Compute average corrected pt
      double ptAve = 0.5*(corrJetPt[corrJetPt.size()-1]+corrJetPt[corrJetPt.size()-2]);
      
      // Fill histogram of all events
      if( hltTag ) {
	++nTotal;
	hTotal->Fill(ptAve);
      }

      // If triggered by probe trigger, fill histogram of
      // events passing
      if( hltProbe ) {
	++nPass;
	hPass->Fill(ptAve);
      }
    } // End if( hltTag )
  } // End loop over entries

  std::cout << "\n\n------------------------------------------\n";
  std::cout << "  Entries                    \t" << nEvts << std::endl;
  std::cout << "  N(passing " << tagTrigger << ") \t" << nTotal << std::endl;
  std::cout << "  N(passing " << probeTrigger << ") \t" << nPass << std::endl;
  std::cout << "------------------------------------------\n\n";

  // Efficiency graph
//   TGraphAsymmErrors *gEff = new TGraphAsymmErrors();
//   gEff->BayesDivide(hPass,hTotal);
//   gEff->SetMarkerStyle(20);

  TH1* hEff = static_cast<TH1*>(hPass->Clone("hEff"));
  hEff->Divide(hPass,hTotal,1,1);
  hEff->SetMarkerStyle(20);

  // Fit efficiency
  TF1* fit = new TF1("efficiency",eff,minPtAve,maxPtAve,3);
  fit->SetParameter(0,1.);
  fit->SetParameter(1,0.5*(minPtAve+maxPtAve));
  fit->SetParameter(2,1.);
  fit->SetLineColor(kBlue);
  //gEff->Fit(fit,"0QR");
  hEff->Fit(fit,"0QR");

  // Determine 99% threshold
  double ptAve99 = fit->GetX(0.99*fit->GetParameter(2),minPtAve,maxPtAve);
  TLine *lPtAve99 = new TLine(ptAve99,0.,ptAve99,fit->GetParameter(2));
  lPtAve99->SetLineWidth(2);
  lPtAve99->SetLineColor(kBlue);

  std::cout << "\\texttt{" << probeTrigger << "} & " << ptAve99 << " \\\\ \n";

  // Label efficiency
  TPaveText *label = util::LabelFactory::createPaveText(3);
  label->AddText(jetLabel+", Runs "+util::toTString(runStart)+" - "+util::toTString(runEnd));
  if( lumi > 0. ) label->AddText("L = "+util::toTString(lumi)+"/pb,  Reference: "+tagTrigger);
  else label->AddText("Tag trigger "+tagTrigger);
  label->AddText("p^{ave}_{T}(#epsilon > 99%) = "+util::toTString(ptAve99)+" GeV");
  label->GetLine(2)->SetTextColor(lPtAve99->GetLineColor());

  TH1 *hFrame = util::HistOps::createRatioFrame(minPtAve,maxPtAve,0.,1.8*fit->GetParameter(2),"p^{ave}_{T} (GeV)",probeTrigger+" Efficiency");
  hFrame->SetLineColor(fit->GetLineColor());
  for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
    hFrame->SetBinContent(i,fit->GetParameter(2));
  }

  TString outNamePrefix = "TurnOn_"+jetType+"_"+probeTrigger+"_"+util::toTString(runStart)+"-"+util::toTString(runEnd);

  TCanvas *canEff = new TCanvas("canEff","Eff",500,500);
  canEff->cd();
  hFrame->Draw();
  //gEff->Draw("PE1same");
  hEff->Draw("PE1same");
  fit->Draw("same");
  lPtAve99->Draw("same");
  label->Draw("same");
  canEff->SaveAs(outNamePrefix+".eps","eps");


  // Label events
  TLegend *legEvts = util::LabelFactory::createLegend(2);
  legEvts->AddEntry(hTotal,tagTrigger,"L");
  //legEvts->AddEntry(hPass,tagTrigger+" && "+probeTrigger,"L");
  legEvts->AddEntry(hPass,probeTrigger,"L");

  TCanvas *canEvts = new TCanvas("canEvts","Evts",500,500);
  canEvts->cd();
  util::HistOps::setYRange(hTotal,2);
  hTotal->Draw("HIST");
  hPass->SetLineColor(kRed);
  hPass->Draw("HISTsame");
  legEvts->Draw("same");
}


// --------------------------------------------------
void run(int nEvts = -1) {
  TString input = "~/lustre/data/Jet_Run2011A-PromptReco-v2_163337-163757/Jet_Run2011A-PromptReco-v2_163337-163757_ak5PF.root";
  int runStart = 163337;
  int runEnd = 163757;
  double lumi = 62.;
  int nBins1 = 35;
  int nBins2 = 70;

  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve30","HltDiJetAve60",nBins1,40,100,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve60","HltDiJetAve80",nBins1,65,130,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve80","HltDiJetAve110",nBins1,80,150,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve110","HltDiJetAve150",nBins1,120,200,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve150","HltDiJetAve190",nBins2,160,270,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve190","HltDiJetAve240",nBins2,200,320,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve240","HltDiJetAve300",nBins2,250,380,nEvts);
  computeTriggerEfficiences(input,runStart,runEnd,lumi,"HltDiJetAve300","HltDiJetAve370",nBins2,350,450,nEvts);
}
