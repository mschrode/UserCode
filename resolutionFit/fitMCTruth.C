// $Id: $

//!  Fit mean response and resolution from
//!  Kalibri::ControlPlotsJetSmearing


#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


unsigned int gN_FILES = 0;


TH1* createTH1FromTH2(const TH2* h2, const TString &name, const TString &xTitle, const TString &xUnit, const TString &yTitle) {
  TH1* h = new TH1D(name,"",h2->GetNbinsX(),h2->GetXaxis()->GetXbins()->GetArray());
  util::HistOps::setAxisTitles(h,xTitle,xUnit,yTitle);

  return h;
}



void fitProfile(const TString &fileName, double nSigCore, TH1* &hMean, TH1* &hMeanGauss, TH1* &hRMS, TH1* &hSigmaGauss, const TString &outNamePrefix) {
  std::cout << "Fitting profiles (" << gN_FILES << ")" << std::endl;

  // Get 2D histogram resp vs pt from file
  std::cout << "  Getting 2D histogram 'response vs ptGen' from file '" << fileName << "'" << std::endl;
  TH2* hRespVsPt = util::FileOps::readTH2(fileName,"MeanResp_hRespVsPtGen","hRespVsPt_File"+util::toTString(gN_FILES));


  // Creating products
  std::cout << "  Creating result histograms" << std::endl;
  hMean = createTH1FromTH2(hRespVsPt,"hMean"+util::toTString(gN_FILES),"p^{gen}_{T}","GeV","Response Mean");
  hMean->SetMarkerStyle(22);
  hMeanGauss = createTH1FromTH2(hRespVsPt,"hMeanGauss"+util::toTString(gN_FILES),"p^{gen}_{T}","GeV","Response Gauss Mean");
  hMeanGauss->SetMarkerStyle(20);
  hRMS = createTH1FromTH2(hRespVsPt,"hRMS"+util::toTString(gN_FILES),"p^{gen}_{T}","GeV","Response Standard Deviation");
  hRMS->SetMarkerStyle(22);
  hSigmaGauss = createTH1FromTH2(hRespVsPt,"hSigmaGauss"+util::toTString(gN_FILES),"p^{gen}_{T}","GeV","Response Gaussian Width");
  hSigmaGauss->SetMarkerStyle(20);


  // Fill response distributions per pt bin
  std::cout << "  Filling response distributions per ptGen bin" << std::endl;
  util::HistVec hResp;
  util::HistOps::fillSlices(hRespVsPt,hResp,"hResp_File"+util::toTString(gN_FILES));
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt) {
    util::HistOps::setAxisTitles(*hIt,"Response","","jets");
    (*hIt)->SetMarkerStyle(20);
  }


  // Fitting response distributions
  std::cout << "  Fitting response distributions per ptGen bin" << std::endl;
  std::vector<TF1*> fGauss;
  int bin = 1;
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt, ++bin) {
    double width = 0.;
    double widthErr = 0.;
    double rms = 0.;
    double rmsErr = 0.;
    TF1* fit = 0;
    if( util::HistOps::fitCoreWidth(*hIt,nSigCore,fit,width,widthErr,rms,rmsErr) ) {
      hMean->SetBinContent(bin,(*hIt)->GetMean());
      hMean->SetBinError(bin,(*hIt)->GetMeanError());
      hMeanGauss->SetBinContent(bin,fit->GetParameter(1));
      hMeanGauss->SetBinError(bin,fit->GetParError(1));
      hRMS->SetBinContent(bin,rms);
      hRMS->SetBinError(bin,rmsErr);
      hSigmaGauss->SetBinContent(bin,width);
      hSigmaGauss->SetBinError(bin,widthErr);
    } else {
      fit = new TF1("fGauss"+util::toTString(gN_FILES)+"_Bin"+util::toTString(bin),"gaus");
    }
    fit->SetLineWidth(1);
    fit->SetLineColor(kRed);    
    fGauss.push_back(fit);
  }


  // Plotting resonse distributions and fits
  std::cout << "  Plotting resonse distributions and fits" << std::endl;
  bin = 1;
  for(util::HistItConst hIt = hResp.begin(); hIt != hResp.end(); ++hIt, ++bin) {
    double ptMin = hRespVsPt->GetXaxis()->GetBinLowEdge(bin);
    double ptMax = hRespVsPt->GetXaxis()->GetBinUpEdge(bin);
    TPaveText* label = util::LabelFactory::createPaveText(1);
    label->AddText(util::toTString(ptMin)+" < p^{gen}_{T} < "+util::toTString(ptMax)+" GeV");
    TCanvas* can = new TCanvas("can","Response PtBin "+util::toTString(bin),500,500);

    can->cd();
    util::HistOps::setYRange(*hIt,2);
    (*hIt)->GetXaxis()->SetRangeUser(0.55,1.45);
    (*hIt)->Draw("PE1");
    fGauss.at(bin-1)->Draw("same");
    label->Draw("same");
    can->SaveAs(outNamePrefix+"_Linear_PtBin"+util::toTString(bin-1)+".eps","eps");

    can->cd();
    util::HistOps::setYRange(*hIt,2,true);
    (*hIt)->GetXaxis()->SetRangeUser(0.,2.);
    (*hIt)->Draw("PE1");
    fGauss.at(bin-1)->Draw("same");
    label->Draw("same");
    can->SetLogy();
    can->SaveAs(outNamePrefix+"_Log_PtBin"+util::toTString(bin-1)+".eps","eps");

    delete label;
    delete can;
  }


  // Creating tex slides
  std::cout << "  Creating LaTeX slides" << std::endl;

  std::ofstream oFile(outNamePrefix+".tex");

  unsigned int nPtBins = hResp.size();
  
  oFile << "\n\n\n% ----- MC Truth Response ---------------------------" << std::endl;
  unsigned int nSlides = nPtBins/3;
  if( nPtBins%3 > 0 ) nSlides++;
  unsigned int ptBin = 0;
  for(unsigned int slide = 0; slide < nSlides; ++slide) {
    oFile << "\n% --------------------------------------------------\n";
    oFile << "\\begin{frame}\n";
    oFile << "  \\frametitle{MC Truth Response (" << slide+1 << "/" << nSlides << ")}\n";
    oFile << "  \\begin{columns}[T] \n";
    for(int col = 0; col < 3; ++col, ++ptBin) {
      oFile << "    \\begin{column}{0.3333\\textwidth} \n";
      oFile << "    \\centering\n";
      if( ptBin < nPtBins ) {
	oFile << "      \\includegraphics[width=\\textwidth]{" << outNamePrefix+"_Linear_PtBin"+util::toTString(ptBin)+"}\\\\ \n";
	oFile << "      \\includegraphics[width=\\textwidth]{" << outNamePrefix+"_Log_PtBin"+util::toTString(ptBin)+"}\\\\ \n";
      }
      oFile << "    \\end{column} \n";
    }
    oFile << "  \\end{columns} \n";
    oFile << "\\end{frame} \n";
  }
  oFile.close();


  gN_FILES++;
}



TF1* fitResolution(TH1* h, const TString &name, double min, double max) {
  TF1* fit = new TF1(name,"sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",min,max);
  fit->SetLineWidth(1);
  fit->SetLineColor(kRed);
  h->Fit(fit,"INR");

  return fit;
}



void fitMCTruth(const TString &fileName) {
  util::StyleSettings::presentationNoTitle();

  TString outNamePrefix = "MCTruthResponse_F"+util::toTString(gN_FILES);

  // Fit distributions
  TH1* hMean = 0;
  TH1* hMeanGauss = 0;
  TH1* hRMS = 0;
  TH1* hSigmaGauss = 0;
  fitProfile(fileName,2.5,hMean,hMeanGauss,hRMS,hSigmaGauss,outNamePrefix);

  // Fit resolution
  TF1* fit = fitResolution(hSigmaGauss,"Res_F"+util::toTString(gN_FILES),80.,2000.);
  
  // Plot response and resolution
  TLegend* legMean = util::LabelFactory::createLegendCol(2,0.7);
  legMean->AddEntry(hMean,"Arithmetic Mean","P");
  legMean->AddEntry(hMeanGauss,"Gaussian Mean","P");

  TCanvas* can = new TCanvas("can","",500,500);
  can->cd();
  TH1* hFrame = util::HistOps::createRatioFrame(hMeanGauss,"Mean Response",0.83,1.27);
  hFrame->Draw();
  hMean->Draw("PE1same");
  hMeanGauss->Draw("PE1same");
  legMean->Draw("same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_MeanResponse.eps","eps");

  TLegend* legRes = util::LabelFactory::createLegendCol(2,0.7);
  legRes->AddEntry(hRMS,"Standard Deviation","P");
  legRes->AddEntry(hSigmaGauss,"Gaussian Width","P");

  can->cd();
  hSigmaGauss->GetYaxis()->SetRangeUser(0.,0.33);
  hSigmaGauss->GetYaxis()->SetTitle("Resolution");
  hSigmaGauss->Draw("PE1");
  hRMS->Draw("PE1same");
  fit->Draw("same");
  legRes->Draw("same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_Resolution.eps","eps");
  
  delete legMean;
  delete legRes;
  delete can;
}
