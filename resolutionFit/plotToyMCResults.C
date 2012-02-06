// $Id: plotToyMCResults.C,v 1.1 2012/02/05 21:37:09 mschrode Exp $

#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"

#define UTILS_AS_HEADER_FILE
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


const double PT_CUT_MIN = 250.;
const double PT_CUT_MAX = 300.;
const TString LABEL_MC = "Toy Simulation";
const TString LABEL_TRUTH = "x^{true}";
const TString LABEL_PT_CUT = "x^{true}";
const TString LABEL_DELTA_PT = "|#Deltax|";
const TString LABEL_PT_BIN = util::toTString(PT_CUT_MIN)+" < "+LABEL_PT_CUT+" < "+util::toTString(PT_CUT_MAX);
const TString LABEL_ASYMMETRY = "Asymmetry";
const TString LABEL_RESPONSE = "Response";
const int MARKER_STYLE_DIS = 21;
const int COLOR_DIS = kBlack;
const int COLOR_PDF = kRed;
const int LINE_WIDTH = 2;


// ------------------------------------------------------------
void smallestInterval(const TH1* h, double frac, int &binMin, int &binMax, bool isSymmetric) {
  int startBin = h->GetMaximumBin();
  binMin = 1;
  binMax = h->GetNbinsX();
  double tot = h->Integral();
  if( frac < 1. && tot > 0. ) {
    int startBinTmpFirst = isSymmetric ? startBin : std::max(startBin-10,1);
    int startBinTmpLast = isSymmetric ? startBinTmpFirst : std::min(startBin+10,h->GetNbinsX());
    for(int startBinTmp = startBinTmpFirst; startBinTmp <= startBinTmpLast; ++startBinTmp) {
      int binMinTmp = startBinTmp;
      int binMaxTmp = startBinTmp;
      double currFrac = h->Integral(binMinTmp,binMaxTmp)/tot;
      while( currFrac < frac ) {
	if( binMinTmp > 1 ) --binMinTmp;
	if( binMaxTmp < h->GetNbinsX() ) ++binMaxTmp;
	currFrac = h->Integral(binMinTmp,binMaxTmp)/tot;
      }
      //std::cout << startBinTmp << ": " << binMinTmp << " (" << h->GetBinCenter(binMinTmp) << ") - " << binMaxTmp << " (" << h->GetBinCenter(binMaxTmp) << ")" << std::endl;
      if( binMaxTmp-binMinTmp < binMax-binMin ) {
	binMin = binMinTmp;
	binMax = binMaxTmp;
      }
    }
  }
}


// ------------------------------------------------------------
void chi2Test(const TH1* h, const TH1* f, bool isSymmetric, double frac, double &chi2, int &n, int &min, int &max) {
  chi2 = 0.;
  n = 0;
  min = 0;
  max = 0;
  smallestInterval(h,frac,min,max,isSymmetric);
  for(int i = min; i <= max; ++i) {
    if( h->GetBinContent(i) > 0. && h->GetBinError(i) > 0. ) {
      ++n;
      int binMin = f->GetXaxis()->FindBin(h->GetXaxis()->GetBinLowEdge(i));
      int binMax = f->GetXaxis()->FindBin(h->GetXaxis()->GetBinUpEdge(i));
      double t = f->Integral(binMin,binMax)/(binMax-binMin+1);
      chi2 += pow( (h->GetBinContent(i) - t)/h->GetBinError(i), 2. );
    }
  }
}


// ------------------------------------------------------------
void printChi2Test(const TString &var, const TH1* h, const TH1* f, bool isSymmetric = false) {
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout << std::endl << " "+var+":" << std::endl;
  for(int i = 0; i < 4; ++i) {
    double frac = 1.-i*0.1;
    double chi2 = 0.;
    int n = 0;
    int min = 1;
    int max = 1;
    chi2Test(h,f,isSymmetric,frac,chi2,n,min,max);
    std::cout << std::setprecision(2) << "  " << frac << " & ";
    std::cout << h->GetXaxis()->GetBinLowEdge(min) << " & " <<  h->GetXaxis()->GetBinUpEdge(max) << " & ";
    std::cout << std::setprecision(2) << chi2 << " & ";
    std::cout << std::setprecision(0) << n << " & ";
    std::cout << std::setprecision(2) << chi2/n << " \\\\ " << std::endl;
  }
  std::cout << std::endl;
}


// ------------------------------------------------------------
void plotDisVsPdf(TH1* hDis, TH1* hPdf, const TString &disVar, const TString &model, const TString &outName) {
  TPaveText* label = util::LabelFactory::createPaveText(2);
  label->AddText(LABEL_MC+", "+model);
  label->AddText(LABEL_PT_BIN);

  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.6,label->GetSize());
  leg->AddEntry(hDis,"Selected events","P");
  leg->AddEntry(hPdf,"Assumed pdf","L");

  hDis->SetTitle("");
  util::HistOps::setAxisTitles(hDis,disVar,"","events",true);
  hDis->SetMarkerStyle(MARKER_STYLE_DIS);
  hDis->SetMarkerColor(COLOR_DIS);
  hDis->SetLineWidth(LINE_WIDTH);
  hDis->SetLineColor(COLOR_DIS);
  hPdf->SetLineWidth(LINE_WIDTH);
  hPdf->SetLineColor(COLOR_PDF);

  TCanvas* can = new TCanvas("can"+outName,outName,500,500);
  can->cd();
  util::HistOps::setYRange(hDis,label->GetSize()+leg->GetNRows());
  hDis->Draw("PE1");
  hPdf->Draw("Lsame");
  label->Draw("same");
  leg->Draw("same");
  gPad->RedrawAxis();
  can->SaveAs(outName+".eps","eps");

  util::HistOps::setYRange(hDis,label->GetSize()+leg->GetNRows(),3E-6);
  if( outName.Contains("Asym") ) util::HistOps::setYRange(hDis,label->GetSize()+leg->GetNRows(),3E-4);
  else if( outName.Contains("Resp") ) util::HistOps::setYRange(hDis,label->GetSize()+leg->GetNRows(),3E-3);
  hDis->Draw("PE1");
  hPdf->Draw("Lsame");
  label->Draw("same");
  leg->Draw("same");
  can->SetLogy();
  gPad->RedrawAxis();
  can->SaveAs(outName+"Log.eps","eps");
}


// ------------------------------------------------------------
void plotToyMCResults(const TString &fileNameResults, const TString &model, const TString &outNamePrefix) {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::setStyleNoteNoTitle();

  // Fitted resolution
  TH1* hParameters = util::FileOps::readTH1(fileNameResults,"hAbsoluteParameters");
  const double fitRes = hParameters->GetBinContent(1);
  const double fitErr = hParameters->GetBinError(1);
  delete hParameters;

  // Plot spectrum
  TH1* hSpecDis = util::FileOps::readTH1(fileNameResults,"hPtGen");
  hSpecDis->SetBinContent(hSpecDis->FindBin(300.),0.);
  hSpecDis->SetBinError(hSpecDis->FindBin(300.),0.);
  const double meanPtGen = hSpecDis->GetMean();
  const double fitResRel = fitRes/meanPtGen;

  util::HistOps::setXRange(hSpecDis,10);
  TH1* hSpecPdf = util::FileOps::readTH1(fileNameResults,"hTruthPDF_0");
  util::HistOps::normHist(hSpecPdf,"width");
  plotDisVsPdf(hSpecDis,hSpecPdf,LABEL_TRUTH,model,outNamePrefix+"_Spectrum");

  // Plot DeltaX
  TH1* hDeltaXDis = util::FileOps::readTH1(fileNameResults,"hDeltaPtJet12_0");
  util::HistOps::setXRange(hDeltaXDis,10);
  TH1* hDeltaXPdf = new TH1D("hDeltaXPdf","",500,0.,hDeltaXDis->GetXaxis()->GetBinUpEdge(hDeltaXDis->GetNbinsX()));
  for(int i = 1; i <= hDeltaXPdf->GetNbinsX(); ++i) {
    hDeltaXPdf->SetBinContent(i,2.*TMath::Gaus(hDeltaXPdf->GetBinCenter(i),0.,fitRes/sqrt(2.),true));
  }
  plotDisVsPdf(hDeltaXDis,hDeltaXPdf,LABEL_DELTA_PT,model,outNamePrefix+"_DeltaX");

  // Plot asymmetry
  TH1* hAsymDis = util::FileOps::readTH1(fileNameResults,"hPtAsym_0");
  hAsymDis->GetXaxis()->SetRangeUser(-0.39,0.39);
  TH1* hAsymPdf = new TH1D("hAsymPdf","",10*hAsymDis->GetNbinsX(),-1.,1.);
  for(int i = 1; i <= hAsymPdf->GetNbinsX(); ++i) {
    hAsymPdf->SetBinContent(i,TMath::Gaus(hAsymPdf->GetBinCenter(i),0.,fitResRel/sqrt(2.),true));
  }
  plotDisVsPdf(hAsymDis,hAsymPdf,LABEL_ASYMMETRY,model,outNamePrefix+"_Asymmetry");

  // Plot response
  TH1* hRespDis = util::FileOps::readTH1(fileNameResults,"hRespMeas_0");
  util::HistOps::setXRange(hRespDis,10);
  TH1* hRespPdf = new TH1D("hRespPdf","",10*hRespDis->GetNbinsX(),0.,2.);
  for(int i = 1; i <= hRespPdf->GetNbinsX(); ++i) {
    hRespPdf->SetBinContent(i,TMath::Gaus(hRespPdf->GetBinCenter(i),1.,fitResRel,true));
  }
  plotDisVsPdf(hRespDis,hRespPdf,LABEL_RESPONSE,model,outNamePrefix+"_Response");

  // Plot likelihood
  TH1* hParScan = util::FileOps::readTH1(fileNameResults,"hParScan_Par0");
  TF1* parabola = util::FileOps::readTF1(fileNameResults,"ParabolaParScan_Par0");
  
  TPaveText* label = util::LabelFactory::createPaveText(1);
  label->AddText(LABEL_MC+", "+model);

  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.6,label->GetSize());
  leg->AddEntry(hParScan,"Likelihood L(#sigma)","L");
  leg->AddEntry(parabola,"Gaussian likelihood","L");

  hParScan->SetTitle("");
  util::HistOps::setAxisTitles(hParScan,"#sigma","","#DeltaF = -2#upoint[lnL(#sigma) - lnL(#hat{#sigma})]");
  hParScan->GetXaxis()->SetNdivisions(505);
  hParScan->SetLineWidth(LINE_WIDTH);
  hParScan->SetLineColor(COLOR_DIS);
  parabola->SetLineWidth(LINE_WIDTH);
  parabola->SetLineColor(kBlue);
  parabola->SetLineStyle(2);

  TCanvas* can = new TCanvas("canDeltaLikelihood","DeltaLikelihood",500,500);
  can->cd();
  util::HistOps::setYRange(hParScan,label->GetSize()+leg->GetNRows());
  hParScan->Draw("L");
  parabola->Draw("Lsame");
  label->Draw("same");
  leg->Draw("same");
  gPad->RedrawAxis();
  can->SaveAs(outNamePrefix+"_DeltaLikelihood.eps","eps");

  

  // Goodness-of-fit
  printChi2Test("Spectrum",hSpecDis,hSpecPdf);
  printChi2Test("DeltaX",hDeltaXDis,hDeltaXPdf);
  printChi2Test("Asymmetry",hAsymDis,hAsymPdf,true);
  printChi2Test("Response",hRespDis,hRespPdf,true);


  std::cout << std::endl << " Mean ptGen = " << meanPtGen << std::endl;
  std::cout << " Expected resolution = " << sqrt( 4.*4. + 1.2*1.2*meanPtGen + 0.05*0.05*meanPtGen*meanPtGen ) << std::endl;
  
}
