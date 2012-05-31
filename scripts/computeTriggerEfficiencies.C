// $Id: computeTriggerEfficiencies.C,v 1.5 2012/01/24 10:16:38 mschrode Exp $

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


//! Trigger efficiency plots from Kalibri ntuple
//!
//! Note: traditionally you would compute the efficiency of a
//! trigger A using a reference trigger B as
//!
//!   eff = N(passing A && passing B) / N(passing B) .
//!
//! However, due to the sometimes large pre-scales there might
//! only be very few events available for this method preventing
//! a reliable measurement. Instead, the efficiency is defined
//! as
//!
//!   eff = N(passing A) / N(passing B) .
//!
//! Of course, now 0 <= eff <= 1 is not true anymore but without
//! knowledge of the pre-scales the previous definition is not
//! exact either. The pt of the turn-on can still be correctly
//! determined as the pt where eff is 99% of the plateau value
//! of eff.
//!
//! TODO: Encapsulate information per trigger into type


//! Fit function for turn-on
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
    chain->Add(util::absolutePath(fileName));
  }
  // Option 2: list of root files
  else {
    std::ifstream filelist;
    filelist.open(util::absolutePath(fileName));
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
void getTriggerInfo(std::vector<TString> &tagTriggers, std::vector<TString> &probeTriggers, std::vector<double> &probeTriggerVals) {
  std::vector<TString> triggers;
  triggers.push_back("HltDiJetAve30"); 
  triggers.push_back("HltDiJetAve60"); 
  triggers.push_back("HltDiJetAve80"); 
  triggers.push_back("HltDiJetAve110");
  triggers.push_back("HltDiJetAve150");
  triggers.push_back("HltDiJetAve190");
  triggers.push_back("HltDiJetAve240");
  triggers.push_back("HltDiJetAve300");
  triggers.push_back("HltDiJetAve370");

  for(size_t i = 0; i < triggers.size()-1; ++i) {
    tagTriggers.push_back(triggers.at(i));
    probeTriggers.push_back(triggers.at(i+1));
  }

  probeTriggerVals.push_back(60); 
  probeTriggerVals.push_back(80); 
  probeTriggerVals.push_back(110);
  probeTriggerVals.push_back(150);
  probeTriggerVals.push_back(190);
  probeTriggerVals.push_back(240);
  probeTriggerVals.push_back(300);
  probeTriggerVals.push_back(370);

  assert( probeTriggerVals.size() == probeTriggers.size() );
}



// --------------------------------------------------
void computeTriggerEfficiences(const TString &fileNames, double lumi, int nEvts = -1, double etaMin = 0.) {

  std::cout << "Preparing input...\n";

  // Style
  util::StyleSettings::setStyleNoteNoTitle();
  TString jetType = util::LabelFactory::jetAlgo(fileNames)+util::LabelFactory::jetType(fileNames);
  TString jetLabel = util::LabelFactory::labelJet(fileNames);
  int nBins = 250;
  if( etaMin > 0. ) nBins = 100;
  double histMin = 0.;
  double histMax = 500.;

  // The trigger names
  std::vector<TString> tagTriggers;
  std::vector<TString> probeTriggers;
  std::vector<double> probeTriggerVals;
  getTriggerInfo(tagTriggers,probeTriggers,probeTriggerVals);

  // Fit ranges
  std::vector<double> ptMin;
  std::vector<double> ptMax;
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    ptMin.push_back(0.7*probeTriggerVals.at(i));
    ptMax.push_back(std::min(1.5*probeTriggerVals.at(i),histMax));
  }

  // Event counters
  std::vector<TH1*> hNTag;
  std::vector<TH1*> hNProbe;
  std::vector<TH1*> hEff;
  std::vector<TF1*> fEff;
  std::vector<TLine*> lPtAve99;

  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    TH1* h = util::HistOps::createTH1D("hNTag_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin),nBins,histMin,histMax,"p^{ave}_{T}","GeV","events");
    h->Sumw2();
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.1);
    hNTag.push_back(h);

    h = static_cast<TH1*>(hNTag.back()->Clone("hNProbe_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin)));
    hNProbe.push_back(h);

    h = static_cast<TH1*>(hNTag.back()->Clone("hEff_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin)));
    hEff.push_back(h);

    TF1* f = new TF1("efficiency_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin),eff,ptMin.at(i),ptMax.at(i),3);
    f->SetLineColor(kBlue);
    f->SetLineWidth(2);
    fEff.push_back(f);

    TLine* l = new TLine(1.,0.,1.,1.);
    l->SetLineWidth(fEff.back()->GetLineWidth());
    l->SetLineColor(fEff.back()->GetLineColor());
    lPtAve99.push_back(l);
  }
  std::vector<double> plateauVals;
  std::vector<double> ptAve99;


  std::cout << "Preparing chain...\n";
  // Set up chain
  TChain *chain = createTChain(fileNames);

  // Set branch addresses
  const int maxNJet = 50;
  
  std::vector<double> corrJetPt(20);

  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetEta[maxNJet];
  float jetCorrL1[maxNJet];
  float jetCorrL2L3[maxNJet];
  bool hlt30 = false; 
  bool hlt60 = false; 
  bool hlt80 = false; 
  bool hlt110 = false;
  bool hlt150 = false;
  bool hlt190 = false;
  bool hlt240 = false;
  bool hlt300 = false;
  bool hlt370 = false;
  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetCorrL1",jetCorrL1);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  chain->SetBranchAddress("HltDiJetAve30",&hlt30); 
  chain->SetBranchAddress("HltDiJetAve60",&hlt60); 
  chain->SetBranchAddress("HltDiJetAve80",&hlt80); 
  chain->SetBranchAddress("HltDiJetAve110",&hlt110);
  chain->SetBranchAddress("HltDiJetAve150",&hlt150);
  chain->SetBranchAddress("HltDiJetAve190",&hlt190);
  chain->SetBranchAddress("HltDiJetAve240",&hlt240);
  chain->SetBranchAddress("HltDiJetAve300",&hlt300);
  chain->SetBranchAddress("HltDiJetAve370",&hlt370);

  // Association of trigger names and booleans
  std::map<TString,bool> triggerDecisions;
  triggerDecisions["HltDiJetAve30"] = false;
  triggerDecisions["HltDiJetAve60"] = false;
  triggerDecisions["HltDiJetAve80"] = false;
  triggerDecisions["HltDiJetAve110"] = false;
  triggerDecisions["HltDiJetAve150"] = false;
  triggerDecisions["HltDiJetAve190"] = false;
  triggerDecisions["HltDiJetAve240"] = false;
  triggerDecisions["HltDiJetAve300"] = false;
  triggerDecisions["HltDiJetAve370"] = false;

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

    if( nObjJet > 1 ) {
      // Sort jets by corrected pt
      for(int i = 0; i < static_cast<int>(corrJetPt.size()); ++i) {
	if( i < nObjJet ) corrJetPt[i] = jetCorrL1[i]*jetCorrL2L3[i]*jetPt[i];
	else corrJetPt[i] = 0.;      
      }
      std::sort(corrJetPt.begin(),corrJetPt.end());
      
      // Compute average corrected pt
      double ptAve = 0.5*(corrJetPt[corrJetPt.size()-1]+corrJetPt[corrJetPt.size()-2]);

      if( std::abs(jetEta[0]) < etaMin || std::abs(jetEta[1]) < etaMin ) continue;
	
      triggerDecisions["HltDiJetAve30"] = hlt30; 
      triggerDecisions["HltDiJetAve60"] = hlt60; 
      triggerDecisions["HltDiJetAve80"] = hlt80; 
      triggerDecisions["HltDiJetAve110"] = hlt110;
      triggerDecisions["HltDiJetAve150"] = hlt150;
      triggerDecisions["HltDiJetAve190"] = hlt190;
      triggerDecisions["HltDiJetAve240"] = hlt240;
      triggerDecisions["HltDiJetAve300"] = hlt300;
      triggerDecisions["HltDiJetAve370"] = hlt370;

      // Loop over triggers
      for(size_t i = 0; i < probeTriggers.size(); ++i) {
	if( ptAve > ptMin.at(i) && ptAve < ptMax.at(i) ) {
	  // Count events passing reference (tag) trigger
	  if( triggerDecisions[tagTriggers.at(i)] ) hNTag.at(i)->Fill(ptAve);
	  // Count events passing probe trigger
	  if( triggerDecisions[probeTriggers.at(i)] ) hNProbe.at(i)->Fill(ptAve);
	}
      }
    } // End if( nObjJet > 1 )
  } // End loop over entries


  // Compute and fit efficiencies,
  // determine 99% threshold
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    hEff.at(i)->Divide(hNProbe.at(i),hNTag.at(i),1,1);
    fEff.at(i)->SetParameter(0,1.);
    fEff.at(i)->SetParameter(1,0.5*(ptMin.at(i)+ptMax.at(i)));
    fEff.at(i)->SetParameter(2,1.);
    hEff.at(i)->Fit(fEff.at(i),"NQRI","",ptMin.at(i),ptMax.at(i));

    plateauVals.push_back(fEff.at(i)->GetParameter(2));
    ptAve99.push_back(fEff.at(i)->GetX(0.99*plateauVals.at(i),ptMin.at(i),ptMax.at(i)));
    lPtAve99.at(i)->SetX1(ptAve99.back());
    lPtAve99.at(i)->SetX2(ptAve99.back());
    lPtAve99.at(i)->SetY2(0.99*plateauVals.at(i));
  }


  // Plot efficiencies
  TString outNamePrefix = "TurnOn_"+jetType;
  if( etaMin > 0. ) outNamePrefix += "_MinEta"+util::toTStringNoPoint(etaMin,1);
  TFile* outFile = new TFile(outNamePrefix+".root","RECREATE");

  for(size_t i = 0; i < probeTriggers.size(); ++i) {

    // Labels
    TPaveText *label = util::LabelFactory::createPaveText(3);
    if( etaMin > 0. ) label->AddText(jetLabel+",  |#eta_{1,2}| > "+util::toTString(etaMin)+",  L = "+util::StyleSettings::luminosity(lumi));
    else label->AddText(jetLabel+",  L = "+util::StyleSettings::luminosity(lumi));
    label->AddText("Reference: "+tagTriggers.at(i));
    label->AddText("p^{ave}_{T}(#epsilon > 99%) = "+util::toTString(ptAve99.at(i))+" GeV");
    label->GetLine(2)->SetTextColor(lPtAve99.at(i)->GetLineColor());
  
    TH1 *hFrame = util::HistOps::createRatioFrame(ptMin.at(i),ptMax.at(i),0.,1.6*plateauVals.at(i),"p^{ave}_{T} (GeV)",probeTriggers.at(i)+" Efficiency");
    hFrame->SetLineColor(fEff.at(i)->GetLineColor());
    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
      hFrame->SetBinContent(bin,plateauVals.at(i));
    }
    
    TCanvas *canEff = new TCanvas("canEff_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin),probeTriggers.at(i)+" #hat{#epsilon}",500,500);
    canEff->cd();
    hFrame->Draw();
    hEff.at(i)->Draw("PE1same");
    fEff.at(i)->Draw("same");
    lPtAve99.at(i)->Draw("same");
    label->Draw("same");
    util::FileOps::toFiles(canEff,outFile,outNamePrefix+"_"+probeTriggers.at(i));
  }

  // All in one plot
  TH1 *hFrame = util::HistOps::createRatioFrame(histMin,histMax,0.,2.45,"p^{ave}_{T} (GeV)","#hat{#epsilon}");
  for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
    hFrame->SetBinContent(bin,1.);
  }
  TCanvas *canEff = new TCanvas("canEffAll_Eta"+util::toTString(etaMin),"efficiency",500,500);
  canEff->cd();
  hFrame->Draw();
  for(int i = static_cast<int>(probeTriggers.size()-1); i >= 0; --i) {
    if( plateauVals.at(i) > 0. ) hEff.at(i)->Scale(1./plateauVals.at(i));
    fEff.at(i)->SetParameter(2,1.);

    hEff.at(i)->SetMarkerStyle(20+(i%4));
    hEff.at(i)->SetMarkerColor(util::StyleSettings::color(i));
    hEff.at(i)->SetLineColor(hEff.at(i)->GetMarkerColor());
    fEff.at(i)->SetLineColor(hEff.at(i)->GetMarkerColor());

    hEff.at(i)->Draw("PE1same");
    fEff.at(i)->Draw("same");
  }
  TPaveText *label = util::LabelFactory::createPaveText(1);
  if( etaMin > 0. ) label->AddText("|#eta_{1,2}| > "+util::toTString(etaMin)+",  "+util::LabelFactory::labelData(util::StyleSettings::luminosity(lumi)));
  else label->AddText(util::LabelFactory::labelData(util::StyleSettings::luminosity(lumi)));
  label->Draw("same");
  TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(4,-0.5,1);
  TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(4,0.5,1);
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    if( i < 4 ) {
      leg1->AddEntry(hEff.at(i),probeTriggers.at(i),"LP");
    } else {
      leg2->AddEntry(hEff.at(i),probeTriggers.at(i),"LP");
    }
  }
  leg1->Draw("same");
  leg2->Draw("same");
  util::FileOps::toFiles(canEff,outFile,outNamePrefix);


  // Print efficiencies
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    std::cout << "\\texttt{" << probeTriggers.at(i) << "} & " << util::toTString(ptAve99.at(i),1) << " \\\\"  << std::endl;
  }

  std::cout << std::endl << std::endl;
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    std::cout << probeTriggers.at(i) << " : " << util::toTString(ptAve99.at(i),1) << std::endl;
  }

  delete outFile;
}


// --------------------------------------------------
void run(int nEvts = -1) {
  TString input = "~/UserCode/mschrode/resolutionFit/input/Analysis2011/Kalibri_JetRun2011A_V3_163337-167151_L1FastJet_AK5PF";
  double lumi = 855.;

  computeTriggerEfficiences(input,lumi,nEvts);
//   computeTriggerEfficiences(input,lumi,nEvts,1.7);
//   computeTriggerEfficiences(input,lumi,nEvts,2.3);
}
