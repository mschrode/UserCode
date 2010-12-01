#include <algorithm>
#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "globalFunctions.h"
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



bool fitWidth(const util::HistVec &hists, double nSig, TH1* &hWidth, TH1* &hRMS) {
  bool ok = true;

  if( static_cast<int>(hists.size()) != hWidth->GetNbinsX() ||
      static_cast<int>(hists.size()) != hRMS->GetNbinsX() ) {
    std::cerr << "ERROR in fitWidth(): number of asymmetry histograms and bins do not match" << std::endl;
    ok = false;
  } else {
    int bin = 1;
    for(util::HistItConst hIt = hists.begin(); hIt != hists.end(); ++hIt, ++bin) {
      double width = 0.;
      double widthErr = 0.;
      double rms = 0.;
      double rmsErr = 0.;
      if( func::fitCoreWidth(*hIt,nSig,width,widthErr,rms,rmsErr) ) {
	hWidth->SetBinContent(bin,sqrt(2.)*width);
	hWidth->SetBinError(bin,sqrt(2.)*widthErr);

	hRMS->SetBinContent(bin,sqrt(2.)*rms);
	hRMS->SetBinError(bin,sqrt(2.)*rmsErr);
      }
    }
  }

  return ok;
}



void setStyle(util::HistVec &hists) {
  for(util::HistIt it = hists.begin(); it != hists.end(); ++it) {
    (*it)->GetXaxis()->SetRangeUser(-1.,1.);
  }
}


void compareDijetAsymmetry(const TString &binCfg = "BinningAdmin.cfg") {
  
  // Prepare parameters
  std::cout << "Preparing parameters" << std::endl;

  util::StyleSettings::presentationNoTitle();
  sampleTools::BinningAdmin admin(binCfg);
  const TString inData = "Tails_Calo_Data_Pt3Cut_Eta0_PtSoft1.root";
  const TString inMC = "Tails_Calo_MCFall10_Pt3Cut_Eta0_PtSoft1.root";

  

  // Read asymmetry 
  std::cout << "Reading histograms" << std::endl;

  util::HistVec hAsymData;
  util::HistVec hAsymMC;
  for(unsigned int ptBin = 0; ptBin < admin.nPtBins(0); ++ptBin) {
    TString histName = "hPtAsym_Eta0_Pt"+util::toTString(ptBin);
    hAsymData.push_back(util::FileOps::readTH1(inData,histName,histName+"_Data"));
    hAsymMC.push_back(util::FileOps::readTH1(inMC,histName,histName+"_MC"));
  }
  setStyle(hAsymData);
  setStyle(hAsymMC);


  
  // Fitting width
  std::cout << "Fitting width" << std::endl;

  std::vector<double> nSig;
  nSig.push_back(1.5);
  nSig.push_back(2.);
  nSig.push_back(2.5);
  nSig.push_back(3.);

  util::HistVec hWidthData;
  util::HistVec hWidthMC;
  for(unsigned int i = 0; i < nSig.size(); ++i) {
    hWidthData.push_back(util::HistOps::createTH1D("hWidthData"+util::toTString(i),admin.ptBinEdges(0),"p^{ave}_{T}","GeV","#sqrt{2} #upoint Gaussian #sigma(Asymmetry)",true));
    hWidthData.back()->GetXaxis()->SetMoreLogLabels();
    hWidthMC.push_back(static_cast<TH1D*>(hWidthData.back()->Clone("hWidthMC"+util::toTString(i))));
    hWidthData.back()->SetMarkerStyle(20);
  }
  
  TH1* hRMSData = util::HistOps::createTH1D("hRMSData",admin.ptBinEdges(0),"p^{ave}_{T}","GeV","#sqrt{2} #upoint Std Deviation(Asymmetry)",true);
  hRMSData->GetXaxis()->SetMoreLogLabels();
  TH1* hRMSMC = static_cast<TH1D*>(hRMSData->Clone("hRMSMC"));
  hRMSData->SetMarkerStyle(24);

  for(unsigned int i = 0; i < nSig.size(); ++i) {
    fitWidth(hAsymMC,nSig[i],hWidthMC[i],hRMSMC);
    fitWidth(hAsymData,nSig[i],hWidthData[i],hRMSData);
  }

  util::HistVec hWidthRatio;
  for(unsigned int i = 0; i < nSig.size(); ++i) {
    hWidthRatio.push_back(util::HistOps::createRatioPlot(hWidthData[i],hWidthMC[i]));
  }
  TH1* hRMSRatio = util::HistOps::createRatioPlot(hRMSData,hRMSMC);



  // Plotting
  std::cout << "Creating labels" << std::endl;

  TPaveText* label = util::LabelFactory::createPaveText(2,-0.5);
  label->AddText(util::LabelFactory::labelJetAlgo(inData,inMC));
  label->AddText(util::LabelFactory::labelEta(admin.etaMin(0),admin.etaMax(0)));

  std::vector<TLegend*> legWidth;
  TLegend* legWidthRatio = util::LabelFactory::createLegendCol(nSig.size()+1,0.6);
  util::LabelFactory::addExtraLegLine(legWidthRatio,"Fit range");
  for(unsigned int i = 0; i < nSig.size(); ++i) {
    legWidth.push_back(util::LabelFactory::createLegendCol(3,0.6));
    util::LabelFactory::addExtraLegLine(legWidth[i],"Fit range "+util::toTString(nSig[i])+" #sigma");
    legWidth[i]->AddEntry(hWidthData[i],"Data","P");
    legWidth[i]->AddEntry(hWidthMC[i],"MC","L");

    legWidthRatio->AddEntry(hWidthRatio[i],util::toTString(nSig[i])+" #sigma");
  }



  TLegend* legRMS = util::LabelFactory::createLegendCol(2,0.4);
  legRMS->AddEntry(hRMSData,"Data","P");
  legRMS->AddEntry(hRMSMC,"MC","L");



  std::cout << "Plotting histograms" << std::endl;

  for(unsigned int i = 0; i < nSig.size(); ++i) {
    TCanvas* canWidth = util::HistOps::createRatioTopCanvas();
    TPad *bottomPadWidth = util::HistOps::createRatioBottomPad();
    TH1 *topFrameWidth = util::HistOps::createRatioTopHist(hWidthMC[i]);
    TH1 *bottomFrameWidth = util::HistOps::createRatioBottomFrame(hWidthMC[i],"p^{ave}_{T}","GeV",0.91,1.29);
    canWidth->cd();
    topFrameWidth->GetYaxis()->SetRangeUser(0.,0.3);
    topFrameWidth->Draw("HISTE");
    hWidthData[i]->Draw("PE1same");
    label->Draw("same");
    legWidth[i]->Draw("same");
    canWidth->SetLogx();
    bottomPadWidth->Draw();
    bottomPadWidth->cd();
    bottomFrameWidth->GetXaxis()->SetMoreLogLabels();
    bottomFrameWidth->Draw();
    hWidthRatio[i]->Draw("PE1same");
    bottomPadWidth->SetLogx();
    canWidth->SaveAs("AsymmetryWidth_nSig"+util::toTString(nSig[i])+".eps","eps");
  }

  TCanvas* canNSig = new TCanvas("canNSig","",500,500);
  canNSig->cd();
  TH1* hWidthRatioFrame = util::HistOps::createRatioFrame(hWidthRatio[0],hWidthRatio[0]->GetXaxis()->GetTitle(),0.9,1.7);
  hWidthRatioFrame->GetXaxis()->SetMoreLogLabels();
  hWidthRatioFrame->Draw("HIST");
  for(unsigned int i = 0; i < nSig.size(); ++i) {
    util::HistOps::setStyleColor(hWidthRatio[i],i);
    hWidthRatio[i]->SetMarkerStyle(20+i);
    hWidthRatio[i]->Draw("PE1same");
  }
  label->Draw("same");
  legWidthRatio->Draw("same");
  canNSig->SetLogx();
  canNSig->SaveAs("AsymmetryWidthRatios.eps","eps");


  TCanvas* canRMS = util::HistOps::createRatioTopCanvas();
  TPad *bottomPadRMS = util::HistOps::createRatioBottomPad();
  TH1 *topFrameRMS = util::HistOps::createRatioTopHist(hRMSMC);
  TH1 *bottomFrameRMS = util::HistOps::createRatioBottomFrame(hRMSMC,"p^{ave}_{T}","GeV",0.91,1.29);
  canRMS->cd();
  topFrameRMS->GetYaxis()->SetRangeUser(0.,0.3);
  topFrameRMS->Draw("HISTE");
  hRMSData->Draw("PE1same");
  label->Draw("same");
  legRMS->Draw("same");
  canRMS->SetLogx();
  bottomPadRMS->Draw();
  bottomPadRMS->cd();
  bottomFrameRMS->GetXaxis()->SetMoreLogLabels();
  bottomFrameRMS->Draw();
  hRMSRatio->Draw("PE1same");
  bottomPadRMS->SetLogx();
  canRMS->SaveAs("AsymmetryRMS.eps","eps");
}


