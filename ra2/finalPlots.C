// $ Id: $

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "THStack.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


class EventProvenance {
public:
  EventProvenance(unsigned int runNumber, unsigned int lumiBlockNumber, unsigned int eventNumber)
    : runNumber_(runNumber), lumiBlockNumber_(lumiBlockNumber), eventNumber_(eventNumber) {};

  bool operator<(const EventProvenance &rhs) const {
    bool lt = false;
    if( runNumber_ < rhs.runNumber_ ) {
      lt = true;
    } else if( runNumber_ > rhs.runNumber_ ) {
      lt = false;
    } else {
      if( lumiBlockNumber_ < rhs.lumiBlockNumber_ ) {
	lt = true;
      } else if( lumiBlockNumber_ > rhs.lumiBlockNumber_ ) {
	lt = false;
      } else {
	if( eventNumber_ < rhs.eventNumber_ ) {
	  lt = true;
	}
      }
    }

    return lt;
  }

  void print() const {
    std::cout << runNumber_ << " \t: " << lumiBlockNumber_ << " \t: " << eventNumber_ << std::endl;
  }

  unsigned int runNumber_;
  unsigned int lumiBlockNumber_;
  unsigned int eventNumber_;
};


void compareDataAndPrediction(const TString &fileNameData, const std::vector<TString> &fileNamesMC, const std::vector<TString> &labelsMC, const TString &histName, const TString &uid, const TString &xTitle, double xMax, const TString &yTitle, double lumi, double minHT, double minMHT) {
  if( fileNamesMC.size() != labelsMC.size() ) {
    std::cerr << "ERROR: Number of MC labels does not match number of files!" << std::endl;
    exit(1);
  }
  
  // The data distribution
  TH1* hData = util::FileOps::readTH1(fileNameData,fileNameData+":///"+histName,uid+"_Data");
  hData->SetTitle("");
  hData->UseCurrentStyle();
  hData->Rebin(4);
  hData->SetNdivisions(505);
  hData->SetMarkerStyle(20);
  hData->GetXaxis()->SetRange(1,hData->FindBin(xMax));
  hData->GetXaxis()->SetTitle(xTitle);
  hData->GetYaxis()->SetTitle(yTitle);
  util::HistOps::setYRange(hData,1,3E-1);

  // The background prediction
  THStack* hMCBgrs = new THStack(uid+"_MCBgrs","");
  std::vector<TH1*> hMCs;
  std::vector<TString> entriesMC;
  double sumEvtsMC = 0;
  for(unsigned int i = 0; i != fileNamesMC.size(); ++i) {
    TH1* hMC = util::FileOps::readTH1(fileNamesMC.at(i),fileNamesMC.at(i)+":///"+histName,uid+"_MC"+util::toTString(i));
    hMC->SetTitle("");
    hMC->UseCurrentStyle();
    hMC->Rebin(4);
    hMC->SetFillColor(1+i);
    hMC->SetLineColor(hMC->GetFillColor());
    hMC->GetXaxis()->SetRange(1,hMC->FindBin(xMax));
    hMCBgrs->Add(hMC);
    hMCs.push_back(hMC);
    //entriesMC.push_back(labelsMC.at(i)+" N = "+util::toTString(static_cast<int>(hMC->Integral())));
    entriesMC.push_back(labelsMC.at(i)+" N = "+util::toTString(hMC->Integral(),0));
    sumEvtsMC += hMC->Integral();
  }

  TPaveText* lab = util::LabelFactory::createPaveText(1);
  lab->AddText("L = "+util::toTString(lumi)+" pb^{-1},  H_{T} > "+util::toTString(minHT)+" GeV,  #slash{H}_{T} > "+util::toTString(minMHT)+" GeV");

  TLegend* leg = util::LabelFactory::createLegendColWithOffset(labelsMC.size()+1,0.5,1);
  leg->AddEntry(hData,"Data N = "+util::toTString(hData->Integral()),"P");
  for(int i = static_cast<int>(entriesMC.size()-1); i >= 0; --i) {
    leg->AddEntry(hMCs.at(i),entriesMC.at(i),"F");
  }

  TCanvas* canData = new TCanvas(uid+"canData",uid,550,550);
  canData->cd();
  hData->Draw("PE");
  lab->Draw("same");
  leg->Draw("same");
  canData->SetLogy();
  canData->SaveAs(uid+"_Data.eps","eps");

  if( fileNamesMC.size() ) {
    TCanvas* canComp = new TCanvas(uid+"canComp",uid,550,550);
    canComp->cd();
    hData->Draw("PE");
    hMCBgrs->Draw("HISTEsame");
    hData->Draw("PEsame");
    lab->Draw("same");
    leg->Draw("same");
    canComp->SetLogy();
    canComp->SaveAs(uid+"_DataMC.eps","eps");
  }

  std::cout << "N(Data) = " << hData->Integral() << ",  N(MC) = " << sumEvtsMC << std::endl;
}


void compareDistributions(const std::vector<TString> &fileNames, const std::vector<TString> &labels, const TString &histName, const TString &uid, const TString &xTitle, double xMax, const TString &yTitle, double lumi, double minHT, double minMHT) {
  if( fileNames.size() != labels.size() ) {
    std::cerr << "ERROR: Number of labels does not match number of files!" << std::endl;
    exit(1);
  }
  
  std::vector<TH1*> hDists;
  std::vector<TString> entries;
  for(unsigned int i = 0; i != fileNames.size(); ++i) {
    TH1* h = util::FileOps::readTH1(fileNames.at(i),fileNames.at(i)+":///"+histName,uid+"_"+util::toTString(i));
    h->SetTitle("");
    h->UseCurrentStyle();
    //h->Rebin(4);
    h->SetNdivisions(505);
    h->GetXaxis()->SetRange(1,h->FindBin(xMax));
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    util::HistOps::setYRange(h,1,3E-1);
    h->SetLineColor(1+i);
    h->SetLineWidth(2);
    hDists.push_back(h);
    entries.push_back(labels.at(i)+": N = "+util::toTString(h->Integral(),0));
  }

  TPaveText* lab = util::LabelFactory::createPaveText(1);
  lab->AddText("L = "+util::toTString(lumi)+" pb^{-1},  H_{T} > "+util::toTString(minHT)+" GeV,  #slash{H}_{T} > "+util::toTString(minMHT)+" GeV");

  TLegend* leg = util::LabelFactory::createLegendColWithOffset(labels.size(),0.8,1);
  for(int i = static_cast<int>(entries.size()-1); i >= 0; --i) {
    leg->AddEntry(hDists.at(i),entries.at(i),"L");
  }

  TCanvas* can = new TCanvas(uid+"can",uid,550,550);
  can->cd();
  hDists.front()->Draw("HIST");
  for(size_t i = 1; i < hDists.size(); ++i) {
    hDists.at(i)->Draw("HISTsame");
  }
  lab->Draw("same");
  leg->Draw("same");
  can->SetLogy();
  gPad->RedrawAxis();
  can->SaveAs(uid+".eps","eps");
}


void printEventProvenance(const TString &fileName, const TString &dir, const TString &outFileName) {
  TChain* chain = util::FileOps::createTChain(fileName,fileName+":///"+dir+"/Provenance");
  unsigned int runNumber = 0;
  unsigned int lumiBlockNumber = 0;
  unsigned int eventNumber = 0;
  chain->SetBranchAddress("RunNumber",&runNumber);
  chain->SetBranchAddress("LumiBlockNumber",&lumiBlockNumber);
  chain->SetBranchAddress("EventNumber",&eventNumber);
  std::vector<EventProvenance> evtPr;
  for(int i = 0; i < chain->GetEntries(); ++i) {
    chain->GetEntry(i);
    evtPr.push_back(EventProvenance(runNumber,lumiBlockNumber,eventNumber));
  }
  std::sort(evtPr.begin(),evtPr.end());
  for(std::vector<EventProvenance>::const_iterator it = evtPr.begin();
      it != evtPr.end(); ++it) {
    it->print();
  }

  delete chain;
}




void finalPlots() {
  util::StyleSettings::presentationNoTitle();

  TString path = "/afs/naf.desy.de/user/m/mschrode/lustre/RA2/FinalPlots";
  TString fileNameData = path+"/Data/FinalPlots_Data_HT_Run2011A_PromptReco_V3TEST_2011-05-18.root";

  std::vector<TString> fileNamesMC;
  std::vector<TString> labelsMC;

  fileNamesMC.push_back(path+"/MC/FinalPlots_MC_DYJetsToLL_TuneD6T_7TeV-madgraph-tauola_V2.root");
  labelsMC.push_back("DY+Jets");
  fileNamesMC.push_back(path+"/MC/FinalPlots_MC_Other_V2.root");
  labelsMC.push_back("Other EWK");
  fileNamesMC.push_back(path+"/MC/FinalPlots_MC_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_V2.root");
  labelsMC.push_back("QCD");
  fileNamesMC.push_back(path+"/MC/FinalPlots_MC_TTJets_TuneZ2_7TeV-madgraph-tauola_V2.root");
  labelsMC.push_back("t#bar{t}+Jets");
  fileNamesMC.push_back(path+"/MC/FinalPlots_MC_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_V2.root");
  labelsMC.push_back("W+Jets");
  fileNamesMC.push_back(path+"/MC/FinalPlots_MC_ZinvisibleJets_7TeV-madgraph_V2.root");
  labelsMC.push_back("Z#nu#bar{#nu}+Jets");


//   // Baseline selection
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"BaselinePlot/HTHT","RA2FinalPlots_V2_BaseLine_HT","H_{T} (GeV)",2150,"Events",191.,350.,150.);
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"BaselinePlot/MHTMHT","RA2FinalPlots_V2_BaseLine_MHT","#slash{H}_{T} (GeV)",950,"Events",191.,350.,150.);
//   // High MHT
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalMHT250Plot/HTHT","RA2FinalPlots_V2_HighMHT_HT","H_{T} (GeV)",2150,"Events",191.,350.,250.);
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalMHT250Plot/MHTMHT","RA2FinalPlots_V2_HighMHT_MHT","#slash{H}_{T} (GeV)",950,"Events",191.,350.,250.);
//   // High HT
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalHT500Plot/HTHT","RA2FinalPlots_V2_HighHT_HT","H_{T} (GeV)",2150,"Events",191.,500.,150.);
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalHT500Plot/MHTMHT","RA2FinalPlots_V2_HighHT_MHT","#slash{H}_{T} (GeV)",950,"Events",191.,500.,150.);


  std::vector<TString> files;
  std::vector<TString> labels;
  files.push_back(path+"/Data/FinalPlots_Data_HT_Run2011A_May10ReReco_V4.root");
  labels.push_back("No TP+BE filter");
  files.push_back(path+"/Data/FinalPlots_Data_HT_Run2011A_May10ReReco_V4_TPBEFilter.root");
  labels.push_back("TP+BE filter");
  compareDistributions(files,labels,"BaselinePlot/HTHT","RA2TPBEFilterImpact_V4_BaseLine_HT","H_{T} (GeV)",2150,"Events",204.2,350.,150.);
  compareDistributions(files,labels,"BaselinePlot/MHTMHT","RA2TPBEFilterImpact_V4_BaseLine_MHT","#slash{H}_{T} (GeV)",950,"Events",204.2,350.,150.);
}
