// $Id: plotToyMCResults.C,v 1.3 2012/03/22 19:24:52 mschrode Exp $

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
const TString LABEL_TRUTH = "p^{true}_{T} (GeV)";
const TString LABEL_PT_CUT = "p^{true}_{T}";
const TString LABEL_DELTA_PT = "|#Deltap_{T}| (GeV)";
const TString LABEL_PT_BIN = util::toTString(PT_CUT_MIN)+" < "+LABEL_PT_CUT+" < "+util::toTString(PT_CUT_MAX)+" GeV";
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
  leg->AddEntry(parabola,"Parabolic likelihood","L");

  hParScan->SetTitle("");
  util::HistOps::setAxisTitles(hParScan,"#sigma (GeV)","","#DeltaF = -2#upoint[lnL(#sigma) - lnL(#hat{#sigma})]");
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


// ------------------------------------------------------------
void plotPullDistribution(const TString &fileNamePrefix, unsigned int nFiles, const TString &model, const TString &outNamePrefix) {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::setStyleNoteNoTitle();

  // set of input files
  std::vector<TString> fileNames;
  for(unsigned int n = 0; n < nFiles; ++n) {
    TString tmp = fileNamePrefix+util::toTString(n)+".root";
    TFile f(tmp,"READ");
    if( f.IsOpen() ) {
      f.Close();
      fileNames.push_back(tmp);
    }
  }

  // determine expected resolution
  double expReso = 0.;
  if( model.Contains("1") ) {
    expReso = 20.;
  } else if( model.Contains("2") ) {
    TH1* h = new TH1D("hSpec","",1000,200.,500.);
    for(unsigned int n = 0; n < fileNames.size(); ++n) {
      // PtGenSpectrum
      TH1* hPtGen = util::FileOps::readTH1(fileNames.at(n),"hPtGen");
      h->Fill(hPtGen->GetMean());
      delete hPtGen;
    }
    double p = h->GetMean();
    delete h;
    expReso = p*sqrt( 4.*4./p/p + 1.2*1.2/p + 0.05*0.05 );
  }

  // fill pull and resolution
  TH1* hPull = new TH1D("hPull",";(#hat{#sigma} - #LT#sigma#GT) / #delta#hat{#sigma};Number of Fits",31,-4.5,4.5);
  hPull->SetMarkerStyle(21);
  hPull->SetMarkerSize(1.4);
  TH1* hReso = new TH1D("hReso",";#hat{#sigma} (GeV);Number of Fits",31,expReso-0.9,expReso+0.9);
  hReso->SetMarkerStyle(21);
  hReso->SetMarkerSize(1.4);

  for(unsigned int n = 0; n < fileNames.size(); ++n) {
    // Fitted resolution
    TH1* hParameters = util::FileOps::readTH1(fileNames.at(n),"hAbsoluteParameters");
    const double fitRes = hParameters->GetBinContent(1);
    const double fitErr = hParameters->GetBinError(1);
    delete hParameters;
    hPull->Fill((fitRes-expReso)/fitErr);
    hReso->Fill(fitRes);
  } // End of loop over files

  // Gaussian fit to distributions
  hReso->Fit("gaus","0ILL");
  TF1* fReso = hReso->GetFunction("gaus");
  fReso->SetName("fReso");
  fReso->SetLineColor(kBlue);
  fReso->SetLineWidth(2);

  hPull->Fit("gaus","0ILL");
  TF1* fPull = hPull->GetFunction("gaus");
  fPull->SetName("fPull");
  fPull->SetLineColor(kBlue);
  fPull->SetLineWidth(2);

  // position of expected resolution
  TLine* lExpRes = new TLine(expReso,0.,expReso,1.05*hReso->GetBinContent(hReso->GetMaximumBin()));
  lExpRes->SetLineWidth(2);
  lExpRes->SetLineStyle(2);
  lExpRes->SetLineColor(hReso->GetLineColor());

  // Labels
  TPaveText* label = util::LabelFactory::createPaveText(1);
  label->AddText(LABEL_MC+", "+model);

  TLegend* legReso = util::LabelFactory::createLegendWithOffset(5,label->GetSize());
  legReso->AddEntry(lExpRes,"Expected #LT#sigma#GT = "+util::toTString(expReso,util::firstSigDigit(fReso->GetParError(1)))+" GeV","L");
  legReso->AddEntry(hReso,"Fitted #hat{#sigma}","P");
  legReso->AddEntry(fReso,"Gaussian fit to distribution","L");
  util::LabelFactory::addExtraLegLine(legReso," mean = "+util::toTString(fReso->GetParameter(1),util::firstSigDigit(fReso->GetParError(1)))+" #pm "+util::toTString(fReso->GetParError(1),util::firstSigDigit(fReso->GetParError(1)))+" GeV");
  util::LabelFactory::addExtraLegLine(legReso," width = "+util::toTString(fReso->GetParameter(2),util::firstSigDigit(fReso->GetParError(2)))+" #pm "+util::toTString(fReso->GetParError(2),util::firstSigDigit(fReso->GetParError(2)))+" GeV");

  TLegend* legPull = util::LabelFactory::createLegendColWithOffset(4,-0.75,label->GetSize());
  legPull->AddEntry(hReso,"Fitted pull (#hat{#sigma} - #LT#sigma#GT) / #delta#hat{#sigma}","P");
  legPull->AddEntry(fReso,"Gaussian fit to distribution","L");
  util::LabelFactory::addExtraLegLine(legPull," mean = "+util::toTString(fPull->GetParameter(1),util::firstSigDigit(fPull->GetParError(1)))+" #pm "+util::toTString(fPull->GetParError(1),util::firstSigDigit(fPull->GetParError(1)))+" GeV");
  util::LabelFactory::addExtraLegLine(legPull," width = "+util::toTString(fPull->GetParameter(2),util::firstSigDigit(fPull->GetParError(2)))+" #pm "+util::toTString(fPull->GetParError(2),util::firstSigDigit(fPull->GetParError(2)))+" GeV");

  util::HistOps::setYRange(hReso,label->GetSize()+legReso->GetNRows());
  util::HistOps::setYRange(hPull,label->GetSize()+legPull->GetNRows());

  TCanvas* canPull = new TCanvas("canPull","Pull",500,500);
  canPull->cd();
  hPull->Draw("PE1");
  fPull->Draw("same");
  label->Draw("same");
  legPull->Draw("same");
  gPad->RedrawAxis();
  canPull->SaveAs(outNamePrefix+"_Pull.eps","eps");

  TCanvas* canReso = new TCanvas("canReso","Reso",500,500);
			     canReso->cd();
  hReso->Draw("PE1");
  lExpRes->Draw("same");
  fReso->Draw("same");
  label->Draw("same");
  legReso->Draw("same");
  gPad->RedrawAxis();
  canReso->SaveAs(outNamePrefix+"_Sigma.eps","eps");
}
