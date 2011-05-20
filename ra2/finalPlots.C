#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


void compareDistribution(const TString &fileNameData, const std::vector<TString> &fileNamesMC, const std::vector<TString> &labelsMC, const TString &histName, const TString &uid, const TString &xTitle, double xMax, const TString &yTitle, double lumi, double minHT, double minMHT) {
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

  // The MC background prediction
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




void finalPlots() {
  util::StyleSettings::presentationNoTitle();

  TString path = "/afs/naf.desy.de/user/m/mschrode/lustre/RA2/FinalPlots";
  TString fileNameData = path+"/Data/FinalPlots_Data_HTPromptReco_V2_2011-05-13.root";
  //TString fileNameData = path+"/Data/FinalPlots_Data_HT_Run2011A_PromptReco_V3TEST_2011-05-18.root";

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


  // Baseline selection
  compareDistribution(fileNameData,fileNamesMC,labelsMC,"BaselinePlot/HTHT","RA2FinalPlots_V2_BaseLine_HT","H_{T} (GeV)",2150,"Events",191.,350.,150.);
  compareDistribution(fileNameData,fileNamesMC,labelsMC,"BaselinePlot/MHTMHT","RA2FinalPlots_V2_BaseLine_MHT","#slash{H}_{T} (GeV)",950,"Events",191.,350.,150.);
  // High MHT
  compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalMHT250Plot/HTHT","RA2FinalPlots_V2_HighMHT_HT","H_{T} (GeV)",2150,"Events",191.,350.,250.);
  compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalMHT250Plot/MHTMHT","RA2FinalPlots_V2_HighMHT_MHT","#slash{H}_{T} (GeV)",950,"Events",191.,350.,250.);
  // High HT
  compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalHT500Plot/HTHT","RA2FinalPlots_V2_HighHT_HT","H_{T} (GeV)",2150,"Events",191.,500.,150.);
  compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalHT500Plot/MHTMHT","RA2FinalPlots_V2_HighHT_MHT","#slash{H}_{T} (GeV)",950,"Events",191.,500.,150.);
  
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalPlot/HTHT","RA2FinalPlots_V3TEST_BaseLine_HT","H_{T} (GeV)",2150,"Events",191,350.,150.);
//   compareDistribution(fileNameData,fileNamesMC,labelsMC,"FinalPlot/MHTMHT","RA2FinalPlots_V3TEST_BaseLine_MHT","#slash{H}_{T} (GeV)",950,"Events",191,350.,150.);
}
