// $Id: computeTriggerEfficiencies.C,v 1.1 2010/11/03 15:10:46 mschrode Exp $

#include <algorithm>
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


// --------------------------------------------------
TChain *createTChain(const TString &fileListName) {
  TChain* chain = new TChain("DiJetTree"); 
  std::ifstream filelist;
  filelist.open(fileListName);
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
    std::cerr << "ERROR opening file '" << fileListName << "'\n";
    exit(1);
  }
  filelist.close();
  return chain;
}



// --------------------------------------------------
void computeTriggerEfficiences(const TString &fileNames, int runStart, int runEnd, double lumi, const TString &tagTrigger, const TString &probeTrigger, unsigned int nBins, double minPtAve, double maxPtAve, int nEvts = -1) {
  util::StyleSettings::presentationNoTitle();

  TString jetType = "";
  TString jetLabel = "";
  if( fileNames.Contains("AK5") ) {
    jetType += "AK5";
    jetLabel += "AK5";
  }
  if( fileNames.Contains("PF") ) {
    jetType += "PF";
    jetLabel += " PF";
  } else if( fileNames.Contains("Calo") ) {
    jetType += "Calo";
    jetLabel += " Calo";
  }
  jetLabel += " Jets";


  std::cout << "Creating histograms...\n";
  // Event counts
  TH1 *hFrameU = util::HistOps::createRatioFrame(minPtAve,maxPtAve,0.,1.4,"Uncorrected p^{ave}_{T} (GeV)",probeTrigger+" Efficiency");
  TH1 *hTotalU = util::HistOps::createTH1D("hTotalU",nBins,minPtAve,maxPtAve,"Uncorrected p^{ave}_{T}","GeV","events");
  TH1 *hPassU = util::HistOps::createTH1D("hPassU",nBins,minPtAve,maxPtAve,"Uncorrected p^{ave}_{T}","GeV","events");
  hPassU->SetLineColor(kRed);

  double minPtAveCorr = minPtAve;
  double maxPtAveCorr = maxPtAve;
  if( jetType.Contains("Calo") ) {
    minPtAveCorr *= 1.5;
    maxPtAveCorr *= 1.5;
  } else if( jetType.Contains("PF") ) {
    minPtAveCorr *= 1.2;
    maxPtAveCorr *= 1.2;
  }

  TH1 *hFrame = util::HistOps::createRatioFrame(minPtAveCorr,maxPtAveCorr,0.001,1.59,"p^{ave}_{T} (GeV)",probeTrigger+" Efficiency");
  TH1 *hTotal = util::HistOps::createTH1D("hTotal",nBins,minPtAveCorr,maxPtAveCorr,"p^{ave}_{T}","GeV","events");
  TH1 *hPass = util::HistOps::createTH1D("hPass",nBins,minPtAveCorr,maxPtAveCorr,"p^{ave}_{T}","GeV","events");


  std::cout << "Adding files to chain...\n";
  // Set up chain
  TChain *chain = createTChain(fileNames);

  // Set branch addresses
  const int maxNJet = 50;
  
  int nTotal = 0;
  int nPass  =0;

  std::vector<double> corrJetPt(10);

  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetCorrL2L3[maxNJet];
  bool hltTag = false;
  bool hltProbe = false;

  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress(tagTrigger,&hltTag);
  chain->SetBranchAddress(probeTrigger,&hltProbe);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);

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

    if( nObjJet > 1 && hltTag ) {
      ++nTotal;

      // Uncorrected ptAve
      double ptAveU = 0.5*(jetPt[0]+jetPt[1]);
      
      // Corrected ptAve
      for(int i = 0; i < static_cast<int>(corrJetPt.size()); ++i) {
	if( i < nObjJet ) corrJetPt[i] = jetCorrL2L3[i]*jetPt[i];
	else corrJetPt[i] = 0.;      
      }
      std::sort(corrJetPt.begin(),corrJetPt.end());
      double ptAve = 0.5*(corrJetPt[corrJetPt.size()-1]+corrJetPt[corrJetPt.size()-2]);
      
      hTotalU->Fill(ptAveU);
      hTotal->Fill(ptAve);

      if( hltProbe ) {
	++nPass;
	hPassU->Fill(ptAveU);
	hPass->Fill(ptAve);
      }
    } // End if( hltTag )
  } // End loop over entries

  std::cout << "\n\n--------------------------------------\n";
  std::cout << "  Entries                    \t" << nEvts << std::endl;
  std::cout << "  Total (" << tagTrigger << ") \t" << nTotal << std::endl;
  std::cout << "  Passed (" << probeTrigger << ") \t" << nPass << std::endl;
  std::cout << "--------------------------------------\n\n";


  // Efficiency graph
  TGraphAsymmErrors *gEffU = new TGraphAsymmErrors();
  gEffU->BayesDivide(hPassU,hTotalU);
  gEffU->SetMarkerStyle(20);

  TGraphAsymmErrors *gEff = new TGraphAsymmErrors();
  gEff->BayesDivide(hPass,hTotal);
  gEff->SetMarkerStyle(20);

  // Determine 99% threshold
  double ptAve99U = maxPtAve;
  for(int i = 0; i < gEffU->GetN(); ++i) {
    if( gEffU->GetY()[i] > 0.99 ) {
      ptAve99U = gEffU->GetX()[i];
      break;
    }
  }
  TLine *lPtAve99U = new TLine(ptAve99U,0.,ptAve99U,1.);
  lPtAve99U->SetLineWidth(2);
  lPtAve99U->SetLineColor(kBlue);

  double ptAve99 = maxPtAveCorr;
  for(int i = 0; i < gEff->GetN(); ++i) {
    if( gEff->GetY()[i] > 0.99 ) {
      ptAve99 = gEff->GetX()[i];
      break;
    }
  }
  TLine *lPtAve99 = new TLine(ptAve99,0.,ptAve99,1.);
  lPtAve99->SetLineWidth(2);
  lPtAve99->SetLineColor(kBlue);

  std::cout << "\\texttt{" << probeTrigger << "} & " << ptAve99U << " & " << ptAve99 << " \\\\ \n";

  // Label efficiency
  TPaveText *label = util::LabelFactory::createPaveText(3);
  label->AddText(jetLabel+", Runs "+util::toTString(runStart)+" - "+util::toTString(runEnd));
  if( lumi > 0. ) label->AddText("Tag trigger "+tagTrigger+" ("+util::toTString(lumi)+"/pb)");
  else label->AddText("Tag trigger "+tagTrigger);

  TPaveText *labelU = static_cast<TPaveText*>(label->Clone());
  label->AddText("p^{ave}_{T}(#epsilon > 99%) = "+util::toTString(ptAve99)+" GeV");
  label->GetLine(2)->SetTextColor(lPtAve99->GetLineColor());
  labelU->AddText("p^{ave}_{T}(#epsilon > 99%) = "+util::toTString(ptAve99U)+" GeV");
  labelU->GetLine(2)->SetTextColor(lPtAve99->GetLineColor());

  TString outNamePrefix = "Efficiency_"+jetType+"_"+probeTrigger+"_"+util::toTString(runStart)+"-"+util::toTString(runEnd)+"_";

  TCanvas *canEffU = new TCanvas("canEffU","EffU",500,500);
  canEffU->cd();
  hFrameU->Draw();
  gEffU->Draw("PE1same");
  lPtAve99U->Draw("same");
  labelU->Draw("same");
  canEffU->SaveAs(outNamePrefix+"Uncorrected.eps","eps");

  TCanvas *canEff = new TCanvas("canEff","Eff",500,500);
  canEff->cd();
  hFrame->Draw();
  gEff->Draw("PE1same");
  lPtAve99->Draw("same");
  label->Draw("same");
  canEff->SaveAs(outNamePrefix+"Corrected.eps","eps");


  // Label events
  TLegend *legEvts = util::LabelFactory::createLegend(2);
  legEvts->AddEntry(hTotalU,tagTrigger,"L");
  legEvts->AddEntry(hPassU,tagTrigger+" && "+probeTrigger,"L");

  TCanvas *canEvts = new TCanvas("canEvts","Evts",500,500);
  canEvts->cd();
  util::HistOps::setYRange(hTotalU,2);
  hTotalU->Draw("HIST");
  hPassU->Draw("HISTsame");
  legEvts->Draw("same");
}
