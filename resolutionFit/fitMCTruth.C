// $Id: fitMCTruth.C,v 1.13 2012/01/24 10:02:14 mschrode Exp $

//!  Fit mean response and resolution from MC truth
//!  histograms in Kalibri skims from
//!  sampleTools::writeDijetSkims.C

#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLine.h"
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


// Global style
const bool archivePlots_ = false;
const TString jetLabel_ = "Anti-k_{T} (R=0.5) PF Jets";
const double labelResponseMin_ = 0.61;
const double labelResponseMax_ = 1.39;
const double labelResponseLogMin_ = 0.09;
const double labelResponseLogMax_ = 1.91;


// ------------------------------------------------------------
void fitProfile(const TH2* h2, std::vector<TH1*> &hProf, double nSigCore, std::vector<TF1*> &fits, TH1* &hMean, TH1* &hSigma, TH1* &hChi2Ndof) {
  //std::cout << "Fitting profiles of '" << h2->GetName() << "'" << std::endl;

  // Histograms of mean values and standard deviations
  TString name = h2->GetName();
  TString xTitle = h2->GetXaxis()->GetTitle();
  TString yTitle = h2->GetYaxis()->GetTitle();
  hMean = new TH1D(name+"_YMean",";"+xTitle+";Gauss Mean",
		   h2->GetNbinsX(),h2->GetXaxis()->GetXbins()->GetArray());
  hMean->SetMarkerStyle(20);

  hSigma = static_cast<TH1*>(hMean->Clone(name+"_YSigma"));
  hSigma->SetYTitle("Gauss Width");

  hChi2Ndof = static_cast<TH1D*>(hMean->Clone(name+"_YChi2Ndof"));
  hChi2Ndof->SetYTitle("#chi^{2} / ndof");
  hChi2Ndof->SetMarkerColor(kBlack);

  // Fill y distributions per x bin
  util::HistOps::fillSlices(h2,hProf,name+"_YDistribution");
  for(util::HistItConst hIt = hProf.begin(); hIt != hProf.end(); ++hIt) {
    util::HistOps::setAxisTitles(*hIt,yTitle,"","jets");
    (*hIt)->UseCurrentStyle();
    (*hIt)->SetMarkerStyle(20);
    //    (*hIt)->Rebin(3);
  }

  // Fitting y distributions
  int bin = 1;
  for(util::HistItConst hIt = hProf.begin(); hIt != hProf.end(); ++hIt, ++bin) {
    double width = 0.;
    double widthErr = 0.;
    TF1* fit = 0;
    if( util::HistOps::fitCoreWidth(*hIt,nSigCore,fit,width,widthErr) ) {
      hMean->SetBinContent(bin,fit->GetParameter(1));
      hMean->SetBinError(bin,fit->GetParError(1));
      hSigma->SetBinContent(bin,width);
      hSigma->SetBinError(bin,widthErr);
      if( fit->GetNDF() > 0 ) {
	hChi2Ndof->SetBinContent(bin,fit->GetChisquare() / fit->GetNDF());
	hChi2Ndof->SetBinError(bin,0.);
      }
    } else {
      fit = new TF1(name+"_YGaussFit_XBin"+util::toTString(bin),"gaus");
      if( hSigma->Integral() ) fit->SetParameter(0,hSigma->Integral("width")/sqrt(2.*M_PI)/20.);
      fit->SetParameter(0,1./sqrt(2.*M_PI)/20.);
      fit->SetParameter(1,1.);
      fit->SetParameter(2,20.);
    }
    fit->SetLineWidth(1);
    fit->SetName(name+"_YGaussFit_XBin"+util::toTString(bin));
    fits.push_back(fit);
  }
}


// ------------------------------------------------------------
TF1* fitResolution(TH1* h, const TString &name, double min, double max) {
  TF1* fit = new TF1(name,"sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",min,max);
  fit->SetParameter(0,1.5);
  fit->SetParameter(1,0.5);
  fit->FixParameter(2,0.);
  fit->SetParameter(3,0.3);

  //    TF1* fit = new TF1(name,"sqrt( sq([0]/x) + sq([1])/x + sq([2]) )",min,max);
  //    fit->SetParameter(0,1.);
  //    fit->SetParameter(1,0.1);
  //    fit->SetParameter(2,0.01);

  fit->SetLineWidth(1);
  fit->SetLineColor(kRed);
  h->Fit(fit,"BINR");		// Option "B" to consider fixed parameters


  return fit;
}


// ------------------------------------------------------------
void fitResolutionVsPtGen(const TString &fileName,
			  const TString &histNamePrefix, unsigned int nEtaBins,
			  double nSigCore, double minPt,
			  std::vector< std::vector<TH1*> > &hProf,
			  std::vector< std::vector<TF1*> > &fProf,
			  std::vector<TH1*> &hChi2Ndof,
			  std::vector<TH1*> &hMean,
			  std::vector<TH1*> &hSigma,
			  std::vector<TF1*> &fSigma
			  ) {

  std::cout << "\nFitting resolution for '" << histNamePrefix << "'" << std::endl;
  std::cout << "  Eta bin: " << std::flush;
  for(unsigned int i = 0; i < nEtaBins; ++i) {
    std::cout << i << "   " << std::flush;
    // Get 2D histogram
    std::vector<TString> histName;
    // In case several histograms are to be added
    if( histNamePrefix.Contains(":") ) {
      int idxFirst = 0;
      int idxLast = histNamePrefix.Index(":");
      while( idxLast >= 0 ) {
	histName.push_back(histNamePrefix(idxFirst,idxLast-idxFirst));
	std::cout << ">>> " << idxFirst << " - " << idxLast << ": " << histName.back() << std::endl;
	idxFirst = idxLast+1;
	idxLast = histNamePrefix.Index(":",idxFirst);
	std::cout << "    " << idxFirst << " - " << idxLast << std::endl;
      }
      histName.push_back(histNamePrefix(idxFirst,histNamePrefix.Length()-idxFirst));
    } else {
      histName.push_back(histNamePrefix);
    }
    TH2* h2 = 0;
    for(std::vector<TString>::iterator it = histName.begin();
	it != histName.end(); ++it) {
      if( it->Contains("Eta?") ) it->ReplaceAll("Eta?","Eta"+util::toTString(i));
      else *it += "_Eta"+util::toTString(i);
      std::cout << "Getting '" << *it << "' from file" << std::endl;
      if( h2 == 0 ) {
	h2 = util::FileOps::readTH2(fileName,*it);
      } else {
	h2->Add(util::FileOps::readTH2(fileName,*it));
      }
    }

    std::vector<TH1*> hProfTmp;
    std::vector<TF1*> fProfTmp;
    TH1* hMeanTmp = 0;
    TH1* hSigmaTmp = 0;
    TH1* hChi2NdofTmp = 0;
    fitProfile(h2,hProfTmp,nSigCore,fProfTmp,hMeanTmp,hSigmaTmp,hChi2NdofTmp);
    for(int bin = 1; bin <= hMeanTmp->FindBin(minPt); ++bin) {
      hMeanTmp->SetBinContent(bin,0);
      hMeanTmp->SetBinError(bin,0);
      hSigmaTmp->SetBinContent(bin,0);
      hSigmaTmp->SetBinError(bin,0);
      hChi2NdofTmp->SetBinContent(bin,0);
      hChi2NdofTmp->SetBinError(bin,0);
    }
    int startOfPrecisionRegion = 1;
    int endOfPrecisionRegion = 10000;
    for(int bin = hMeanTmp->GetNbinsX()/2; bin >= 1; --bin) {
      if( hSigmaTmp->GetBinContent(bin) &&
	  hSigmaTmp->GetBinError(bin) / hSigmaTmp->GetBinContent(bin) > 0.15 ) {
	startOfPrecisionRegion = bin;
	break;
      }
    }
    for(int bin = hMeanTmp->GetNbinsX()/2; bin <= hMeanTmp->GetNbinsX(); ++bin) {
      if( hSigmaTmp->GetBinContent(bin) &&
	  hSigmaTmp->GetBinError(bin) / hSigmaTmp->GetBinContent(bin) > 0.15 ) {
	endOfPrecisionRegion = bin;
	break;
      }
    }
    for(int bin = 1; bin <= hMeanTmp->GetNbinsX(); ++bin) {
      if( bin <= startOfPrecisionRegion || bin >= endOfPrecisionRegion ) {
	hMeanTmp->SetBinContent(bin,0);
	hMeanTmp->SetBinError(bin,0);
	hSigmaTmp->SetBinContent(bin,0);
	hSigmaTmp->SetBinError(bin,0);
	hChi2NdofTmp->SetBinContent(bin,0);
	hChi2NdofTmp->SetBinError(bin,0);
      }
    }
    hProf.push_back(hProfTmp);
    fProf.push_back(fProfTmp);
    hMean.push_back(hMeanTmp);
    hSigma.push_back(hSigmaTmp);
    hChi2Ndof.push_back(hChi2NdofTmp);

    fSigma.push_back(fitResolution(hSigmaTmp,histNamePrefix+"_ResolutionFit",minPt,std::min(hSigmaTmp->GetBinCenter(endOfPrecisionRegion),hSigmaTmp->GetXaxis()->GetBinUpEdge(hSigmaTmp->GetNbinsX()))));

    delete h2;
  }
  std::cout << std::endl;
}


// ------------------------------------------------------------
void fitResolutionVsPtGen(const TString &fileName,
			  const TString &histNamePrefix, unsigned int nEtaBins,
			  double nSigCore, double minPt,
			  std::vector<TH1*> &hSigma,
			  std::vector<TF1*> &fSigma
			  ) {

  std::vector< std::vector<TH1*> > hProf;
  std::vector< std::vector<TF1*> > fProf;
  std::vector<TH1*> hChi2Ndof;
  std::vector<TH1*> hMean;
  fitResolutionVsPtGen(fileName,histNamePrefix,nEtaBins,nSigCore,minPt,hProf,fProf,hChi2Ndof,hMean,hSigma,fSigma);
  for(unsigned int i = 0; i < hProf.size(); ++i) {
    for(unsigned int j = 0; j < hProf.at(i).size(); ++j) {
      delete hProf.at(i).at(j);
      delete fProf.at(i).at(j);
    }
  }
  for(unsigned int i = 0; i < hChi2Ndof.size(); ++i) {
    delete hChi2Ndof.at(i);
    delete hMean.at(i);
  }
}


// Plot MC-truth resolution for one selection
// Different eta bins are compared
// ------------------------------------------------------------
void plotMCTruth(const TString &fileName, const TString &histNamePrefix, const sampleTools::BinningAdmin &binningAdmin, double nSigCore, double minPt, const TString &outNamePrefix, bool plotGaussFit = true) {
  
  // Fit response distributions in eta and ptGen bins,
  // obtain mean and width, and fit width vs ptGen
  std::vector< std::vector<TH1*> > hProf;
  std::vector< std::vector<TF1*> > fProf;
  std::vector<TH1*> hChi2Ndof;
  std::vector<TH1*> hMean;
  std::vector<TH1*> hSigma;
  std::vector<TF1*> fSigma;
  fitResolutionVsPtGen(fileName,histNamePrefix,binningAdmin.nEtaBins(),nSigCore,minPt,
		       hProf,fProf,hChi2Ndof,hMean,hSigma,fSigma);

  // Set style
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    hSigma.at(etaBin)->SetMarkerStyle(20+etaBin);
    hSigma.at(etaBin)->SetMarkerColor(util::StyleSettings::color(etaBin));
    hSigma.at(etaBin)->SetLineColor(hSigma.at(etaBin)->GetMarkerColor());
    fSigma.at(etaBin)->SetLineColor(hSigma.at(etaBin)->GetMarkerColor());
    fSigma.at(etaBin)->SetLineWidth(2);

    hMean.at(etaBin)->SetMarkerStyle(hSigma.at(etaBin)->GetMarkerStyle());
    hMean.at(etaBin)->SetMarkerColor(hSigma.at(etaBin)->GetMarkerColor());
    hMean.at(etaBin)->SetLineColor(hSigma.at(etaBin)->GetLineColor());

    hChi2Ndof.at(etaBin)->SetMarkerStyle(hSigma.at(etaBin)->GetMarkerStyle());
    hChi2Ndof.at(etaBin)->SetMarkerColor(hSigma.at(etaBin)->GetMarkerColor());
    hChi2Ndof.at(etaBin)->SetLineColor(hSigma.at(etaBin)->GetLineColor());

    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      hProf.at(etaBin).at(ptBin)->SetMarkerStyle(20);
      hProf.at(etaBin).at(ptBin)->SetMarkerColor(kBlack);
      hProf.at(etaBin).at(ptBin)->SetLineColor(hProf.at(etaBin).at(ptBin)->GetMarkerColor());
      hProf.at(etaBin).at(ptBin)->SetLineWidth(2);
      fProf.at(etaBin).at(ptBin)->SetLineColor(kRed);
      fProf.at(etaBin).at(ptBin)->SetLineWidth(2);
    }
  }
      
  // Create ratio plots of histogram and fit function
  std::vector< std::vector<TH1*> > hProfOverFit;
  std::vector<TH1*> hSigmaOverFit;
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    hSigmaOverFit.push_back(util::HistOps::createRatioPlot(hSigma.at(etaBin),fSigma.at(etaBin)));
    std::vector<TH1*> h;
    for(unsigned int ptBin = 0; ptBin < hProf.at(etaBin).size(); ++ptBin) {
      fProf.at(etaBin).at(ptBin)->SetRange(0.,2.);
      h.push_back(util::HistOps::createRatioPlot(hProf.at(etaBin).at(ptBin),fProf.at(etaBin).at(ptBin)));
    }
    hProfOverFit.push_back(h);
  }

  // Norm profile plots and fits
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < hProf.at(etaBin).size(); ++ptBin) {
      if( hProf.at(etaBin).at(ptBin)->Integral() ) {
 	double norm = 1./hProf.at(etaBin).at(ptBin)->Integral("width");
 	fProf.at(etaBin).at(ptBin)->SetParameter(0,norm*fProf.at(etaBin).at(ptBin)->GetParameter(0));
 	hProf.at(etaBin).at(ptBin)->Scale(norm);
      }
    }
  }


  // Labels for plots with all eta bins in one
  TPaveText* label = util::LabelFactory::createPaveText(2,-0.58);
  label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
  label->AddText(jetLabel_);

  TLegend* leg = util::LabelFactory::createLegendCol(binningAdmin.nEtaBins(),0.4);
  TLegend* legLines = util::LabelFactory::createLegendCol(binningAdmin.nEtaBins(),0.4);
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    leg->AddEntry(hSigma.at(etaBin),util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)),"P");
    legLines->AddEntry(fSigma.at(etaBin),util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)),"L");
  }

  // Canvases
  TCanvas* canPlain = new TCanvas("can","",500,500);
  TCanvas* canRatio = util::HistOps::createRatioTopCanvas();
  TPad *bPad = util::HistOps::createRatioBottomPad();

  // Resolution plots
  canPlain->cd();
  TH1* hFrame = new TH1D("hFrame",";p^{gen}_{T} (GeV);Resolution",10000,
			 std::max(minPt-0.99,hSigma.front()->GetXaxis()->GetBinLowEdge(1)),
			 1.2*hSigma.front()->GetXaxis()->GetBinUpEdge(hSigma.front()->GetNbinsX()));
  hFrame->GetYaxis()->SetRangeUser(3E-3,0.36);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->GetXaxis()->SetNoExponent();
  hFrame->Draw();
  for(int etaBin = static_cast<int>(binningAdmin.nEtaBins()-1); etaBin >=0; --etaBin) {
    fSigma.at(etaBin)->Draw("same");
    hSigma.at(etaBin)->Draw("PE1same");
  }
  label->Draw("same");
  leg->Draw("same");
  canPlain->SetLogx();
  canPlain->SaveAs(outNamePrefix+"_Resolution.eps","eps");

  canPlain->cd();
  hFrame->Draw();
  for(int etaBin = static_cast<int>(binningAdmin.nEtaBins()-1); etaBin >=0; --etaBin) {
    fSigma.at(etaBin)->SetLineStyle(1+etaBin);
    fSigma.at(etaBin)->Draw("same");
  }
  label->Draw("same");
  legLines->Draw("same");
  canPlain->SetLogx();
  canPlain->SaveAs(outNamePrefix+"_ResolutionFit.eps","eps");

  canRatio->cd();
  TH1 *tFrame = util::HistOps::createRatioTopHist(hFrame);
  tFrame->GetYaxis()->SetRangeUser(3E-3,0.39);
  tFrame->Draw();
  for(int etaBin = static_cast<int>(binningAdmin.nEtaBins()-1); etaBin >=0; --etaBin) {
    fSigma.at(etaBin)->Draw("same");
    hSigma.at(etaBin)->Draw("PE1same");
  }
  label->Draw("same");
  leg->Draw("same");
  bPad->Draw();
  bPad->cd();
  TH1 *bFrame = util::HistOps::createRatioBottomFrame(hFrame,0.81,1.19);
  bFrame->SetLineWidth(2);
  bFrame->Draw();
  for(int etaBin = static_cast<int>(binningAdmin.nEtaBins()-1); etaBin >=0; --etaBin) {
    hSigmaOverFit.at(etaBin)->Draw("PE1same");
  }
  canRatio->SetLogx();
  bPad->SetLogx();
  canRatio->SaveAs(outNamePrefix+"_ResolutionFitRatio.eps","eps");

  // Mean response plots
  canPlain->cd();
  hFrame->GetYaxis()->SetTitle("Mean Response");
  hFrame->GetYaxis()->SetRangeUser(0.71,1.59);
  for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
    hFrame->SetBinContent(bin,1.);
    hFrame->SetBinError(bin,0.);
  }
  hFrame->SetLineStyle(2);
  hFrame->SetLineWidth(2);
  hFrame->Draw();
  for(int etaBin = static_cast<int>(binningAdmin.nEtaBins()-1); etaBin >=0; --etaBin) {
    hMean.at(etaBin)->Draw("PE1same");
  }
  label->Draw("same");
  leg->Draw("same");
  canPlain->SetLogx();
  canPlain->SaveAs(outNamePrefix+"_MeanResponse.eps","eps");

  // Chi2/nodf resolution plots
  canPlain->cd();
  hFrame->GetYaxis()->SetTitle("#chi^{2} / ndof");
  hFrame->GetYaxis()->SetRangeUser(0.,8.9);
  hFrame->Draw();
  for(int etaBin = static_cast<int>(binningAdmin.nEtaBins()-1); etaBin >=0; --etaBin) {
    hChi2Ndof.at(etaBin)->Draw("PE1same");
  }
  label->Draw("same");
  leg->Draw("same");
  canPlain->SetLogx();
  canPlain->SaveAs(outNamePrefix+"_ResolutionFitChi2.eps","eps");

  // Resolution plots per etaBin
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    delete label;
    label = util::LabelFactory::createPaveText(3,-0.58);
    label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
    label->AddText(jetLabel_);
    label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)));
    
    canRatio->cd();
    tFrame->Draw();
    hSigma.at(etaBin)->SetMarkerStyle(20);
    hSigma.at(etaBin)->SetMarkerColor(kBlack);
    hSigma.at(etaBin)->SetLineColor(kBlack);
    fSigma.at(etaBin)->SetLineColor(kRed);
    fSigma.at(etaBin)->SetLineStyle(1);
    fSigma.at(etaBin)->Draw("same");
    hSigma.at(etaBin)->Draw("PE1same");
    label->Draw("same");
    bPad = util::HistOps::createRatioBottomPad();
    bPad->Draw();
    bPad->cd();
    bFrame->SetLineWidth(2);
    bFrame->SetLineColor(fSigma.at(etaBin)->GetLineColor());
    bFrame->Draw();
    hSigmaOverFit.at(etaBin)->SetMarkerStyle(20);
    hSigmaOverFit.at(etaBin)->SetMarkerColor(kBlack);
    hSigmaOverFit.at(etaBin)->SetLineColor(kBlack);
    hSigmaOverFit.at(etaBin)->Draw("PE1same");
    canRatio->SetLogx();
    bPad->SetLogx();
    canRatio->SaveAs(outNamePrefix+"_ResolutionFitRatio_EtaBin"+util::toTString(etaBin)+".eps","eps");
  }

  // Response distribution plots
  canPlain->SetLogx(0);
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < hProf.at(etaBin).size(); ++ptBin) {
      delete label;
      label = util::LabelFactory::createPaveText(3,-0.8);
      label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
      label->AddText(jetLabel_);
      label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin))+", "+util::toTString(hSigma.at(etaBin)->GetXaxis()->GetBinLowEdge(ptBin+1),0)+" < p^{gen}_{T} < "+util::toTString(hSigma.at(etaBin)->GetXaxis()->GetBinUpEdge(ptBin+1),0)+" GeV");

      TH1* hRespDist = hProf.at(etaBin).at(ptBin);
      canPlain->cd();
      util::HistOps::setAxisTitles(hRespDist,"Response","","jets",true);
      hRespDist->GetXaxis()->SetRangeUser(labelResponseMin_,labelResponseMax_);
      util::HistOps::setYRange(hRespDist,label->GetSize());
      hRespDist->SetLineWidth(2);
      hRespDist->Draw("HIST");
      label->Draw("same");
      canPlain->SetLogy(0);
      canPlain->SaveAs(outNamePrefix+"_Response_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");

      canPlain->cd();
      hRespDist->GetXaxis()->SetRangeUser(labelResponseLogMin_,labelResponseLogMax_);
      util::HistOps::setYRange(hRespDist,label->GetSize(),3E-5);
      hRespDist->Draw("HIST");
      label->Draw("same");
      canPlain->SetLogy(1);
      canPlain->SaveAs(outNamePrefix+"_ResponseLog_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");
    }
  }
  
  if( plotGaussFit ) {
    canRatio->SetLogx(0);
    for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
      for(unsigned int ptBin = 0; ptBin < hProf.at(etaBin).size(); ++ptBin) {
	delete tFrame;
	delete bFrame;
	delete label;
	delete leg;
      
 	TF1* fitRestricted = static_cast<TF1*>(fProf.at(etaBin).at(ptBin)->Clone("fitRestricted"));
 	double xMin = fitRestricted->GetParameter(1)-nSigCore*fitRestricted->GetParameter(2);
 	double xMax = fitRestricted->GetParameter(1)+nSigCore*fitRestricted->GetParameter(2);
 	fitRestricted->SetRange(xMin,xMax);

	label = util::LabelFactory::createPaveText(3,-0.8);
	label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
	label->AddText(jetLabel_);
	label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin))+", "+util::toTString(hSigma.at(etaBin)->GetXaxis()->GetBinLowEdge(ptBin+1),0)+" < p^{gen}_{T} < "+util::toTString(hSigma.at(etaBin)->GetXaxis()->GetBinUpEdge(ptBin+1),0)+" GeV");
	leg = util::LabelFactory::createLegendColWithOffset(1,-0.43,3);
	leg->AddEntry(fitRestricted,"#chi^{2} / ndof = "+util::toTString(hChi2Ndof.at(etaBin)->GetBinContent(ptBin+1),2),"L");
	
	canRatio->cd();
	canRatio->SetLogy(1);
	tFrame = util::HistOps::createRatioTopHist(hProf.at(etaBin).at(ptBin));

 	util::HistOps::setAxisTitles(tFrame,"","","jets",true);
 	util::HistOps::setYRange(tFrame,3,3E-4); 
 	tFrame->GetXaxis()->SetRangeUser(labelResponseLogMin_,labelResponseLogMax_);     
 	tFrame->Draw("HIST");
 	fProf.at(etaBin).at(ptBin)->SetLineStyle(2);      
 	fProf.at(etaBin).at(ptBin)->Draw("same");
 	fitRestricted->Draw("same");
 	label->Draw("same");
 	bPad  = util::HistOps::createRatioBottomPad();
 	bPad->Draw();
 	bPad->cd();
	bFrame = util::HistOps::createRatioBottomFrame(tFrame,"Response","",0.51,1.99);
 	bFrame->GetXaxis()->SetRangeUser(labelResponseLogMin_,labelResponseLogMax_);     
 	bFrame->Draw();
 	hProfOverFit.at(etaBin).at(ptBin)->Draw("PE1same");
 	label->Draw("same");
 	canRatio->SaveAs(outNamePrefix+"_ResponseFitLog_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");

	label = util::LabelFactory::createPaveText(3,-0.8);
	label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
	label->AddText(jetLabel_);
	label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin))+", "+util::toTString(hSigma.at(etaBin)->GetXaxis()->GetBinLowEdge(ptBin+1),0)+" < p^{gen}_{T} < "+util::toTString(hSigma.at(etaBin)->GetXaxis()->GetBinUpEdge(ptBin+1),0)+" GeV");
	canRatio->cd();
	canRatio->SetLogy(0);
	util::HistOps::setYRange(tFrame,4);
	tFrame->GetXaxis()->SetRangeUser(labelResponseMin_,labelResponseMax_);
	tFrame->Draw("HIST");
	fProf.at(etaBin).at(ptBin)->Draw("same");
	fitRestricted->Draw("same");
	label->Draw("same");
	leg->Draw("same");
	bPad  = util::HistOps::createRatioBottomPad();
	bPad->Draw();
	bPad->cd();
	bFrame->GetYaxis()->SetRangeUser(0.65,1.35);
	bFrame->GetXaxis()->SetRangeUser(labelResponseMin_,labelResponseMax_);
	bFrame->Draw();
	hProfOverFit.at(etaBin).at(ptBin)->Draw("PE1same");
	canRatio->SaveAs(outNamePrefix+"_ResponseFit_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");

	delete fitRestricted;
      }
    }
  }

  // Print fit parameters
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    std::cout << "$" << binningAdmin.etaMin(etaBin) << " - " << binningAdmin.etaMax(etaBin) << std::flush;
    for(int i = 0; i < fSigma.at(etaBin)->GetNpar(); ++i) {
      std::cout << "$ & $" << util::toTString(fSigma.at(etaBin)->GetParameter(i),3) << " \\pm " << util::toTString(fSigma.at(etaBin)->GetParError(i),3) << std::flush;
    }
    std::cout << "$ \\\\" << std::endl;
  }

  for(unsigned int i = 0; i < hProf.size(); ++i) {
    for(unsigned int j = 0; j < hProf.at(i).size(); ++j) {
      delete hProf.at(i).at(j);
      delete fProf.at(i).at(j);
      delete hProfOverFit.at(i).at(j);
    }
  }
  for(unsigned int i = 0; i < hChi2Ndof.size(); ++i) {
    delete hChi2Ndof.at(i);
    delete hMean.at(i);
    delete hSigma.at(i);
    delete fSigma.at(i);
    delete hSigmaOverFit.at(i);
  }
  delete label;
  delete leg;
  delete legLines;
  delete hFrame;
  delete tFrame;
  delete bFrame;
  delete canPlain;
  delete bPad;
  delete canRatio;

  if( archivePlots_ ) {
    gROOT->ProcessLine(".! tar -zcf "+outNamePrefix+".tar.gz "+outNamePrefix+"*.eps");
    gROOT->ProcessLine(".! rm "+outNamePrefix+"*.eps");
    std::cout << "Plots stored in "+outNamePrefix+".tar.gz" << std::endl;
  }
}


// Compare MC-truth resolution for different selections
// One plot per eta bin
// ------------------------------------------------------------
void plotMCTruth(const std::vector<TString> &fileName, const std::vector<TString> &histNamePrefix, const std::vector<TString> &legEntries, const sampleTools::BinningAdmin &binningAdmin, double nSigCore, double minPt, const TString &outNamePrefix, bool plotProfiles = false, bool normProfilesToMax = false, int colorOffset = 0) {

  assert( legEntries.size() == histNamePrefix.size() );
  if( fileName.size() > 1 ) assert( legEntries.size() == fileName.size() );

  std::vector< std::vector< std::vector<TH1*> > > hProf;
  std::vector< std::vector<TH1*> > hMean;
  std::vector< std::vector<TH1*> > hSigma;
  std::vector< std::vector<TF1*> > fSigma;
  for(unsigned int i = 0; i < histNamePrefix.size(); ++i) {
    std::vector<TH1*> tmpHMean;
    std::vector<TH1*> tmpHSigma;
    std::vector<TF1*> tmpFSigma;
    std::vector< std::vector<TH1*> > tmpHProf;
    std::vector< std::vector<TF1*> > fProf;
    std::vector<TH1*> hChi2Ndof;
    fitResolutionVsPtGen((fileName.size()>1?fileName.at(i):fileName.front()),histNamePrefix.at(i),binningAdmin.nEtaBins(),nSigCore,minPt,tmpHProf,fProf,hChi2Ndof,tmpHMean,tmpHSigma,tmpFSigma);
    hProf.push_back(tmpHProf);
    hMean.push_back(tmpHMean);
    hSigma.push_back(tmpHSigma);
    fSigma.push_back(tmpFSigma);
    for(unsigned int k = 0; k < fProf.size(); ++k) {
      for(unsigned int j = 0; j < fProf.at(k).size(); ++j) {
	delete fProf.at(k).at(j);
	tmpHProf.at(k).at(j)->Rebin(3);
      }
    }
    for(unsigned int k = 0; k < hChi2Ndof.size(); ++k) {
      delete hChi2Ndof.at(k);
    }
  }

  TCanvas* can = new TCanvas("can","",500,500);
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    TPaveText* label = util::LabelFactory::createPaveText(3,-0.58);
    label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
    label->AddText(jetLabel_);
    label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)));

    TLegend* leg = util::LabelFactory::createLegendCol(histNamePrefix.size(),0.4);
    for(unsigned int i = 0; i < histNamePrefix.size(); ++i) {
      leg->AddEntry(hSigma.at(i).at(etaBin),legEntries.at(i),"P");
      hSigma.at(i).at(etaBin)->SetMarkerStyle(20+i);
      hSigma.at(i).at(etaBin)->SetMarkerColor(util::StyleSettings::color(i+colorOffset));
      hSigma.at(i).at(etaBin)->SetLineColor(hSigma.at(i).at(etaBin)->GetMarkerColor());
      hSigma.at(i).at(etaBin)->SetLineWidth(2);
      fSigma.at(i).at(etaBin)->SetLineColor(hSigma.at(i).at(etaBin)->GetMarkerColor());
      fSigma.at(i).at(etaBin)->SetLineWidth(hSigma.at(i).at(etaBin)->GetLineWidth());
      hMean.at(i).at(etaBin)->SetLineColor(hSigma.at(i).at(etaBin)->GetMarkerColor());
      hMean.at(i).at(etaBin)->SetLineWidth(hSigma.at(i).at(etaBin)->GetLineWidth());
      hMean.at(i).at(etaBin)->SetMarkerColor(hSigma.at(i).at(etaBin)->GetMarkerColor());
      hMean.at(i).at(etaBin)->SetMarkerStyle(hSigma.at(i).at(etaBin)->GetMarkerStyle());
    }

    TH1* hFrame = new TH1D("hFrame",";p^{gen}_{T} (GeV);Resolution",1000,
			   std::max(minPt-0.99,hSigma.front().front()->GetXaxis()->GetBinLowEdge(1)),1.2*hSigma.front().front()->GetXaxis()->GetBinUpEdge(hSigma.front().front()->GetNbinsX()));
    hFrame->GetYaxis()->SetRangeUser(3E-3,0.39);
    hFrame->GetXaxis()->SetMoreLogLabels();
    hFrame->GetXaxis()->SetNoExponent();

    can->cd();
    hFrame->Draw();
    for(int i = static_cast<int>(histNamePrefix.size()-1); i >= 0; --i) {
      fSigma.at(i).at(etaBin)->Draw("same");
      hSigma.at(i).at(etaBin)->Draw("PE1same");
    }
    label->Draw("same");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(outNamePrefix+"_EtaBin"+util::toTString(etaBin)+".eps","eps");

    hFrame->GetYaxis()->SetTitle("Mean Response");
    hFrame->GetYaxis()->SetRangeUser(0.81,1.39);
    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
      hFrame->SetBinContent(bin,1.);
    }
    hFrame->SetLineStyle(2);
    hFrame->SetLineWidth(2);
    can->cd();
    hFrame->Draw();
    for(int i = static_cast<int>(histNamePrefix.size()-1); i >= 0; --i) {
      hMean.at(i).at(etaBin)->Draw("PE1same");
    }
    label->Draw("same");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(outNamePrefix+"_MeanResponse_EtaBin"+util::toTString(etaBin)+".eps","eps");
    delete hFrame;

    if( plotProfiles ) {
      // Response distribution plots
      can->SetLogx(0);
      for(unsigned int ptBin = 0; ptBin < hProf.front().at(etaBin).size(); ++ptBin) {
 	delete label;
	delete leg;

 	label = util::LabelFactory::createPaveText(3,-0.8);
 	label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
 	label->AddText(jetLabel_);
 	label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin))+", "+util::toTString(hSigma.front().at(etaBin)->GetXaxis()->GetBinLowEdge(ptBin+1),0)+" < p^{gen}_{T} < "+util::toTString(hSigma.front().at(etaBin)->GetXaxis()->GetBinUpEdge(ptBin+1),0)+" GeV");

	leg = util::LabelFactory::createLegendColWithOffset(histNamePrefix.size(),-0.4,label->GetSize());
	for(unsigned int i = 0; i < legEntries.size(); ++i) {
	  leg->AddEntry(hProf.at(i).at(etaBin).at(ptBin),legEntries.at(i),"L");
	}

	double totalProbNorm = 0;
 	for(unsigned int i = 0; i < hProf.size(); ++i) {
 	  totalProbNorm += hProf.at(i).at(etaBin).at(ptBin)->Integral("width");
	}
	totalProbNorm = 1./totalProbNorm;
 	for(unsigned int i = 0; i < hProf.size(); ++i) {
	  if( normProfilesToMax ) {
	    hProf.at(i).at(etaBin).at(ptBin)->Scale(totalProbNorm);
	  } else {
	    util::HistOps::normHist(hProf.at(i).at(etaBin).at(ptBin),"width");
	  }
 	  hProf.at(i).at(etaBin).at(ptBin)->SetLineColor(util::StyleSettings::color(i+colorOffset));
 	  hProf.at(i).at(etaBin).at(ptBin)->SetLineStyle(1+i);
 	  hProf.at(i).at(etaBin).at(ptBin)->SetLineWidth(2);
 	}

 	can->cd();
 	TH1* hRef = hProf.front().at(etaBin).at(ptBin);
 	util::HistOps::setAxisTitles(hRef,"Response","","jets",true);
 	hRef->GetXaxis()->SetRangeUser(labelResponseMin_,labelResponseMax_);
 	util::HistOps::normHist(hRef,"width");
 	util::HistOps::setYRange(hRef,label->GetSize());
 	hRef->Draw("HIST");
 	for(unsigned int i = 1; i < hProf.size(); ++i) {
 	  hProf.at(i).at(etaBin).at(ptBin)->Draw("HISTsame");
 	}	
 	label->Draw("same");
 	leg->Draw("same");
 	can->SetLogy(0);
 	can->SaveAs(outNamePrefix+"_Response_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");

 	can->cd();
 	hRef->GetXaxis()->SetRangeUser(labelResponseLogMin_,labelResponseLogMax_);
 	util::HistOps::setYRange(hRef,label->GetSize(),3E-5);
 	hRef->Draw("HIST");
 	for(unsigned int i = 1; i < hProf.size(); ++i) {
 	  hProf.at(i).at(etaBin).at(ptBin)->Draw("HISTsame");
 	}	
 	label->Draw("same");
 	leg->Draw("same");
 	can->SetLogy(1);
	gPad->RedrawAxis();
 	can->SaveAs(outNamePrefix+"_ResponseLog_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");
      }
    }
    
    delete label;
    delete leg;
  }

  for(unsigned int i = 0; i < hSigma.size(); ++i) {
    for(unsigned int j = 0; j < hSigma.at(i).size(); ++j) {
      delete hMean.at(i).at(j);
      delete hSigma.at(i).at(j);
      delete fSigma.at(i).at(j);
    }
  }
  for(unsigned int i = 0; i < hProf.size(); ++i) {
    for(unsigned int j = 0; j < hProf.at(i).size(); ++j) {
      for(unsigned int k = 0; k < hProf.at(i).at(j).size(); ++k) {
	delete hProf.at(i).at(j).at(k);
      }
    }
  }
  delete can;  

  if( archivePlots_ ) {
    gROOT->ProcessLine(".! tar -zcf "+outNamePrefix+".tar.gz "+outNamePrefix+"*.eps");
    gROOT->ProcessLine(".! rm "+outNamePrefix+"*.eps");
    std::cout << "Plots stored in "+outNamePrefix+".tar.gz" << std::endl;
  }
}

// ------------------------------------------------------------
void plotMCTruth(const TString &fileName, const std::vector<TString> &histNamePrefix, const std::vector<TString> &legEntries, const sampleTools::BinningAdmin &binningAdmin, double nSigCore, double minPt, const TString &outNamePrefix, bool plotProfiles = false, bool normProfilesToMax = false, int colorOffset = 0) {
  std::vector<TString> name;
  name.push_back(fileName);
  plotMCTruth(name,histNamePrefix,legEntries,binningAdmin,nSigCore,minPt,outNamePrefix,plotProfiles,normProfilesToMax,colorOffset);
}




// Plot MC-truth resolution for different
// Gaussian fit ranges (nSigCore)
// Different eta bins are compared
// ------------------------------------------------------------
void plotGaussianFitRange(const TString &fileName, const TString &histNamePrefix, const sampleTools::BinningAdmin &binningAdmin, double minPt, const TString &outNamePrefix) {

  // Tested nSigCore
  std::vector<double> nSigCore;
  nSigCore.push_back(1.5);
  nSigCore.push_back(2.0);
  nSigCore.push_back(2.5);
  
  // Fit response distributions in eta and ptGen bins,
  // obtain mean and width, and fit width vs ptGen
  std::vector< std::vector<TH1*> > hSigma;
  std::vector< std::vector<TF1*> > fSigma;
  for(unsigned int i = 0; i < nSigCore.size(); ++i) {
    std::vector<TH1*> hSigmaTmp;
    std::vector<TF1*> fSigmaTmp;
    fitResolutionVsPtGen(fileName,histNamePrefix,binningAdmin.nEtaBins(),nSigCore.at(i),minPt,hSigmaTmp,fSigmaTmp);
    hSigma.push_back(hSigmaTmp);
    fSigma.push_back(fSigmaTmp);

    for(unsigned int etaBin = 0; etaBin < hSigmaTmp.size(); ++etaBin) {
      hSigmaTmp.at(etaBin)->SetMarkerStyle(20+i);
      hSigmaTmp.at(etaBin)->SetMarkerColor(util::StyleSettings::color(i));
      hSigmaTmp.at(etaBin)->SetLineColor(hSigmaTmp.at(etaBin)->GetMarkerColor());
      fSigmaTmp.at(etaBin)->SetLineColor(hSigmaTmp.at(etaBin)->GetMarkerColor());
      fSigmaTmp.at(etaBin)->SetLineWidth(2);
    }
  }

  TCanvas* can = new TCanvas("can","",500,500);

  TLegend* leg = util::LabelFactory::createLegendCol(nSigCore.size(),0.4);
  for(unsigned int i = 0; i < nSigCore.size(); ++i) {
    leg->AddEntry(hSigma.at(i).front(),"#mu #pm "+util::toTString(nSigCore.at(i))+"#sigma","P");
  }
  
  TH1* hFrame = new TH1D("hFrame",";p^{gen}_{T} (GeV);Resolution",1000,
			 std::max(minPt-0.99,hSigma.front().front()->GetXaxis()->GetBinLowEdge(1)),1.2*hSigma.front().front()->GetXaxis()->GetBinUpEdge(hSigma.front().front()->GetNbinsX()));
  hFrame->GetYaxis()->SetRangeUser(3E-3,0.39);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->GetXaxis()->SetNoExponent();
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    TPaveText* label = util::LabelFactory::createPaveText(3,-0.58);
    label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
    label->AddText(jetLabel_);
    label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)));

    can->cd();
    hFrame->Draw();
    for(unsigned int i = 0; i < nSigCore.size(); ++i) {
      fSigma.at(i).at(etaBin)->Draw("same");
      hSigma.at(i).at(etaBin)->Draw("PE1same");
    }
    label->Draw("same");
    leg->Draw("same");
    can->SetLogx();
    can->SaveAs(outNamePrefix+"_NSigCore_EtaBin"+util::toTString(etaBin)+".eps","eps");
    
    delete label;
  }

  for(unsigned int i = 0; i < hSigma.size(); ++i) {
    for(unsigned int j = 0; j < hSigma.at(i).size(); ++j) {
      delete hSigma.at(i).at(j);
      delete fSigma.at(i).at(j);
    }
  }

  delete leg;
  delete hFrame;
  delete can;
}



// ------------------------------------------------------------
void plotDeltaR(const TString &fileName, const TString &histNamePrefix, const sampleTools::BinningAdmin &binningAdmin, double minPt, double maxPt, const TString &outNamePrefix) {

  
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    // Get DeltaR distributions for all pt from file
    std::vector<TH1*> hProf;
    TH2* h2 = util::FileOps::readTH2(fileName,histNamePrefix+"_Eta"+util::toTString(etaBin));
    util::HistOps::fillSlices(h2,hProf,histNamePrefix+"_Profile");

    // The plotted DeltaR distribution
    TH1* hDeltaR = static_cast<TH1*>(hProf.front()->Clone(outNamePrefix+"DeltaR_EtaBin"+util::toTString(etaBin)));
    hDeltaR->Reset();
    hDeltaR->Clear();

    // Select ptGen range
    unsigned int firstProf = h2->GetXaxis()->FindBin(minPt)-1;
    unsigned int lastProf = h2->GetXaxis()->FindBin(maxPt)-1;
    for(unsigned int i = firstProf; i < lastProf; ++i) {
      hDeltaR->Add(hProf.at(i));
    }
    util::HistOps::normHist(hDeltaR,"width");
    for(std::vector<TH1*>::iterator it = hProf.begin(); it != hProf.end(); ++it) {
      delete *it;
    }

    // Set style
    TPaveText* label = util::LabelFactory::createPaveText(3,-0.8);
    label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
    label->AddText(jetLabel_);
    label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin))+", "+util::toTString(minPt)+" < p^{gen}_{T} < "+util::toTString(maxPt)+" GeV");

    util::HistOps::setAxisTitles(hDeltaR,"#DeltaR(jet^{gen},jet^{det})","","jets",true);
    hDeltaR->UseCurrentStyle();
    util::HistOps::setYRange(hDeltaR,label->GetSize(),3E-4);
    hDeltaR->GetXaxis()->SetRangeUser(0.,0.49);
    hDeltaR->GetXaxis()->SetNdivisions(505);

    TH1* hDeltaRSel = static_cast<TH1*>(hDeltaR->Clone(outNamePrefix+"DeltaRSelected_EtaBin"+util::toTString(etaBin)));
    for(int bin = hDeltaRSel->FindBin(0.1); bin <= hDeltaRSel->GetNbinsX(); ++bin) {
      hDeltaRSel->SetBinContent(bin,0);
      hDeltaRSel->SetBinError(bin,0);
    }
    hDeltaRSel->SetFillColor(33);

    TLegend* leg = util::LabelFactory::createLegendColWithOffset(1,-0.4,3);
    leg->AddEntry(hDeltaRSel,"= "+util::toTString(100.*hDeltaRSel->Integral("width"),0)+"%","F");
    
    TCanvas* can = new TCanvas("can","DeltaR EtaBin"+util::toTString(etaBin),500,500);
    can->cd();
    hDeltaRSel->Draw("HIST");
    hDeltaR->Draw("HISTsame");
    label->Draw("same");
    leg->Draw("same");
    can->SetLogy();
    can->SaveAs(outNamePrefix+"_DeltaR_EtaBin"+util::toTString(etaBin)+".eps","eps");
    
    delete hDeltaR;
    delete hDeltaRSel;
    delete label;
    delete leg;
    delete can;
  }
}


// Compare one and two Gaussian fits
// ------------------------------------------------------------
void plotTwoGaussFit(const TString &fileName, const TString &histNamePrefix, const sampleTools::BinningAdmin &binningAdmin, double nSigCore, double minPt, const TString &outNamePrefix) {
  

  std::vector<TH1*> hChi2Ndof;
  std::vector<TH1*> hChi2Ndof2Gauss;

  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    // Obtain response distributions in ptGen bins,
    // perform standard, one-Gaussian fit
    TString histName = histNamePrefix;
    if( histName.Contains("Eta?") ) histName.ReplaceAll("Eta?","Eta"+util::toTString(etaBin));
    else histName += "_Eta"+util::toTString(etaBin);
    TH2* h2 = util::FileOps::readTH2(fileName,histName);

    std::vector<TH1*> hProf;
    std::vector<TF1*> fProf;
    std::vector<TF1*> fProf2Gauss;
    TH1* hChi2NdofTmp;
    TH1* hMean;
    TH1* hSigma;
    fitProfile(h2,hProf,nSigCore,fProf,hMean,hSigma,hChi2NdofTmp);

    hChi2NdofTmp->GetYaxis()->SetRangeUser(3E-4,8.9);
    hChi2NdofTmp->SetMarkerStyle(20+etaBin);
    hChi2NdofTmp->SetMarkerColor(util::StyleSettings::color(etaBin));
    hChi2NdofTmp->SetLineColor(hChi2NdofTmp->GetMarkerColor());
    util::HistOps::setAxisTitles(hChi2NdofTmp,"p^{gen}_{T}","GeV","#chi^{2} / ndof");
    TH1* hChi2Ndof2GaussTmp = static_cast<TH1*>(hChi2NdofTmp->Clone("hChi2Ndof2Gauss_Eta"+util::toTString(etaBin)));

    // Fit sum of two Gaussian within range of one-Gaus fit
    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      TString name = fProf.at(ptBin)->GetName();
      double min = 1.;
      double max = 1.;
      fProf.at(ptBin)->GetRange(min,max);
      TF1* f2Gauss = new TF1(name+"2Gauss","[0]*TMath::Exp(-0.5*sq(([1]-x[0])/[2])) + [0]*[3]*TMath::Exp(-0.5*sq(([1]-x[0])/[4]))",min,max);
      f2Gauss->SetParameter(0,fProf.at(ptBin)->GetParameter(0));
      f2Gauss->SetParameter(1,1.);
      f2Gauss->SetParameter(2,fProf.at(ptBin)->GetParameter(2));
      f2Gauss->SetParameter(3,0.);
      f2Gauss->SetParameter(4,f2Gauss->GetParameter(2));
      int fitRes = hProf.at(ptBin)->Fit(f2Gauss,"QBINR");
      fProf2Gauss.push_back(f2Gauss);

      double chi2ndof = f2Gauss->GetChisquare() / f2Gauss->GetNDF();
      if( fitRes == 0 && chi2ndof < 1000. ) {
	hChi2Ndof2GaussTmp->SetBinContent(ptBin+1,chi2ndof);
	hChi2Ndof2GaussTmp->SetBinError(ptBin+1,0.);
      }
    }
    int startOfPrecisionRegion = 1;
    int endOfPrecisionRegion = 10000;
    for(int bin = hSigma->GetNbinsX()/2; bin >= 1; --bin) {
      if( hSigma->GetBinContent(bin) &&
	  hSigma->GetBinError(bin) / hSigma->GetBinContent(bin) > 0.15 ) {
	startOfPrecisionRegion = bin;
	break;
      }
    }
    startOfPrecisionRegion = max(startOfPrecisionRegion,hSigma->FindBin(minPt));
    for(int bin = hSigma->GetNbinsX()/2; bin <= hSigma->GetNbinsX(); ++bin) {
      if( hSigma->GetBinContent(bin) &&
	  hSigma->GetBinError(bin) / hSigma->GetBinContent(bin) > 0.15 ) {
	endOfPrecisionRegion = bin;
	break;
      }
    }
    for(int bin = 1; bin <= hSigma->GetNbinsX(); ++bin) {
      if( bin <= startOfPrecisionRegion || bin >= endOfPrecisionRegion ) {
	hChi2NdofTmp->SetBinContent(bin,0);
	hChi2NdofTmp->SetBinError(bin,0);
	hChi2Ndof2GaussTmp->SetBinContent(bin,0);
	hChi2Ndof2GaussTmp->SetBinError(bin,0);
      }
    }
    hChi2Ndof.push_back(hChi2NdofTmp);
    hChi2Ndof2Gauss.push_back(hChi2Ndof2GaussTmp);


    // Create ratio plots of histogram and fit function
    std::vector<TH1*> hProfOverFit;
    std::vector<TH1*> hProfOver2Gauss;
    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      hProfOverFit.push_back(util::HistOps::createRatioPlot(hProf.at(ptBin),fProf.at(ptBin)));
      hProfOver2Gauss.push_back(util::HistOps::createRatioPlot(hProf.at(ptBin),fProf2Gauss.at(ptBin)));
    }

    // Norm profile plots and fits
    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      if( hProf.at(ptBin)->Integral() ) {
 	double norm = 1./hProf.at(ptBin)->Integral("width");
 	fProf.at(ptBin)->SetParameter(0,norm*fProf.at(ptBin)->GetParameter(0));
 	fProf2Gauss.at(ptBin)->SetParameter(0,norm*fProf2Gauss.at(ptBin)->GetParameter(0));
 	hProf.at(ptBin)->Scale(norm);
      }
    }

    // Set style
    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      hProf.at(ptBin)->SetMarkerStyle(20);
      hProf.at(ptBin)->SetMarkerColor(kBlack);
      hProf.at(ptBin)->SetLineColor(hProf.at(ptBin)->GetMarkerColor());
      hProf.at(ptBin)->SetLineWidth(2);

      fProf.at(ptBin)->SetLineColor(kRed);
      fProf.at(ptBin)->SetLineWidth(2);
      hProfOverFit.at(ptBin)->SetMarkerStyle(0);
      hProfOverFit.at(ptBin)->SetLineColor(fProf.at(ptBin)->GetLineColor());
      hProfOverFit.at(ptBin)->SetLineWidth(fProf.at(ptBin)->GetLineWidth());

      fProf2Gauss.at(ptBin)->SetLineColor(kBlue);
      fProf2Gauss.at(ptBin)->SetLineWidth(fProf.at(ptBin)->GetLineWidth());
      hProfOver2Gauss.at(ptBin)->SetMarkerStyle(0);
      hProfOver2Gauss.at(ptBin)->SetLineColor(fProf2Gauss.at(ptBin)->GetLineColor());
      hProfOver2Gauss.at(ptBin)->SetLineWidth(fProf2Gauss.at(ptBin)->GetLineWidth());
    }

    // Canvases
    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      // Indicate fit range
      fProf.at(ptBin)->SetRange(0.,2.);
      fProf2Gauss.at(ptBin)->SetRange(0.,2.);

      TF1* fitRestricted = static_cast<TF1*>(fProf.at(ptBin)->Clone("fitRestricted"));
      double xMin = fitRestricted->GetParameter(1)-nSigCore*fitRestricted->GetParameter(2);
      double xMax = fitRestricted->GetParameter(1)+nSigCore*fitRestricted->GetParameter(2);
      fitRestricted->SetRange(xMin,xMax);
      TF1* fitRestricted2Gauss = static_cast<TF1*>(fProf2Gauss.at(ptBin)->Clone("fitRestricted2Gauss"));
      fitRestricted2Gauss->SetRange(xMin,xMax);

      TH1* ratioRestricted = static_cast<TH1*>(hProfOverFit.at(ptBin)->Clone("ratioRestricted"));
      TH1* ratioRestricted2Gauss = static_cast<TH1*>(hProfOver2Gauss.at(ptBin)->Clone("ratioRestricted2Gauss"));
      for(int bin = 1; bin <= ratioRestricted->GetNbinsX(); ++bin) {
	if( ratioRestricted->GetXaxis()->GetBinCenter(bin) < xMin ||
	    ratioRestricted->GetXaxis()->GetBinCenter(bin) > xMax ) {
	  ratioRestricted->SetBinContent(bin,0.);
	  ratioRestricted->SetBinError(bin,0.);
	  ratioRestricted2Gauss->SetBinContent(bin,0.);
	  ratioRestricted2Gauss->SetBinError(bin,0.);
	}
      }

      fProf.at(ptBin)->SetLineStyle(2);
      fProf2Gauss.at(ptBin)->SetLineStyle(2);
      hProfOverFit.at(ptBin)->SetLineStyle(2);
      hProfOver2Gauss.at(ptBin)->SetLineStyle(2);

      // Labels
      TPaveText* label = util::LabelFactory::createPaveText(3,-0.8);
      label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
      label->AddText(jetLabel_);
      label->AddText(util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin))+", "+util::toTString(hChi2Ndof.at(etaBin)->GetXaxis()->GetBinLowEdge(ptBin+1),0)+" < p^{gen}_{T} < "+util::toTString(hChi2Ndof.at(etaBin)->GetXaxis()->GetBinUpEdge(ptBin+1),0)+" GeV");
      TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.9,3);
      leg->AddEntry(fitRestricted,"G_{A}(#mu,#sigma):  #chi^{2} / ndof = "+util::toTString(hChi2Ndof.at(etaBin)->GetBinContent(ptBin+1),2),"L");
      leg->AddEntry(fitRestricted2Gauss,"G_{A}(#mu,#sigma_{A}) + G_{B}(#mu,#sigma_{B}):  #chi^{2} / ndof = "+util::toTString(hChi2Ndof2Gauss.at(etaBin)->GetBinContent(ptBin+1),2),"L");

      // Plot
      TCanvas* canRatio = util::HistOps::createRatioTopCanvas();
      canRatio->cd();
      canRatio->SetLogy(1);
      TH1* tFrame = util::HistOps::createRatioTopHist(hProf.at(ptBin));
      util::HistOps::setAxisTitles(tFrame,"","","jets",true);
      util::HistOps::setYRange(tFrame,3,3E-4); 
      tFrame->GetXaxis()->SetRangeUser(labelResponseLogMin_,labelResponseLogMax_);     
      tFrame->Draw("HIST");
      fProf.at(ptBin)->SetLineStyle(2);      
      fProf.at(ptBin)->Draw("same");
      fitRestricted->Draw("same");
      fProf2Gauss.at(ptBin)->Draw("same");
      fitRestricted2Gauss->Draw("same");
      label->Draw("same");
      TPad* bPad  = util::HistOps::createRatioBottomPad();
      bPad->Draw();
      bPad->cd();
      TH1* bFrame = util::HistOps::createRatioBottomFrame(tFrame,"Response","",0.51,1.99);
      bFrame->GetXaxis()->SetRangeUser(labelResponseLogMin_,labelResponseLogMax_);     
      bFrame->Draw();
      hProfOverFit.at(ptBin)->Draw("HISTsameE");
      hProfOver2Gauss.at(ptBin)->Draw("HISTsameE");
      ratioRestricted->Draw("HISTsameE");
      ratioRestricted2Gauss->Draw("HISTsameE");
      label->Draw("same");
      gPad->RedrawAxis();
      canRatio->SaveAs(outNamePrefix+"_OneVsTwoGauss_ResponseFitLog_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");

      canRatio->cd();
      canRatio->SetLogy(0);
      util::HistOps::setYRange(tFrame,5);
      tFrame->GetXaxis()->SetRangeUser(labelResponseMin_,labelResponseMax_);
      tFrame->Draw("HIST");
      fProf.at(ptBin)->SetLineStyle(2);      
      fProf.at(ptBin)->Draw("same");
      fitRestricted->Draw("same");
      fitRestricted2Gauss->Draw("same");
      fProf2Gauss.at(ptBin)->Draw("same");
      label->Draw("same");
      leg->Draw("same");
      bPad  = util::HistOps::createRatioBottomPad();
      bPad->Draw();
      bPad->cd();
      delete bFrame;
      bFrame = util::HistOps::createRatioBottomFrame(tFrame,"Response","",0.61,1.39);
      bFrame->GetXaxis()->SetRangeUser(labelResponseMin_,labelResponseMax_);     
      bFrame->Draw();
      hProfOverFit.at(ptBin)->Draw("HISTsameE");
      hProfOver2Gauss.at(ptBin)->Draw("HISTsameE");
      ratioRestricted->Draw("HISTsameE");
      ratioRestricted2Gauss->Draw("HISTsameE");
      label->Draw("same");
      gPad->RedrawAxis();
      canRatio->SaveAs(outNamePrefix+"_OneVsTwoGauss_ResponseFit_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+".eps","eps");

      delete fitRestricted;
      delete fitRestricted2Gauss;
      delete ratioRestricted;
      delete ratioRestricted2Gauss;
      delete label;
      delete leg;
      delete tFrame;
      delete bFrame;
      delete canRatio;
    } // End of loop over pt bins

    delete hMean;
    delete hSigma;
    for(unsigned int ptBin = 0; ptBin < hProf.size(); ++ptBin) {
      delete hProf.at(ptBin);
      delete fProf.at(ptBin);
      delete fProf2Gauss.at(ptBin);
      delete hProfOverFit.at(ptBin);
      delete hProfOver2Gauss.at(ptBin);
    }
  } // End of loop over eta bins

  // Plot chi2/ndof fit probability
  TPaveText* label = util::LabelFactory::createPaveText(2,-0.58);
  label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
  label->AddText(jetLabel_);
  TLegend* leg = util::LabelFactory::createLegendCol(binningAdmin.nEtaBins(),0.4);
  for(unsigned int etaBin = 0; etaBin < binningAdmin.nEtaBins(); ++etaBin) {
    leg->AddEntry(hChi2Ndof.at(etaBin),util::LabelFactory::labelEtaGen(binningAdmin.etaMin(etaBin),binningAdmin.etaMax(etaBin)),"P");
  }

  TCanvas* can = new TCanvas("can","#chi^{2}/ndof",500,500);
  can->cd();
  hChi2Ndof.back()->Draw("P");
  for(int i = hChi2Ndof.size()-2; i >= 0; --i) {
    hChi2Ndof.at(i)->Draw("Psame");
  }
  label->Draw("same");
  leg->Draw("same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_OneVsTwoGauss_Chi2NdofOneGauss.eps","eps");

  can->cd();
  hChi2Ndof2Gauss.back()->Draw("P");
  for(int i = hChi2Ndof2Gauss.size()-2; i >= 0; --i) {
    hChi2Ndof2Gauss.at(i)->Draw("Psame");
  }
  label->Draw("same");
  leg->Draw("same");
  can->SetLogx();
  can->SaveAs(outNamePrefix+"_OneVsTwoGauss_Chi2NdofTwoGauss.eps","eps");  

  for(unsigned int i = 0; i < hChi2Ndof.size(); ++i) {
    delete hChi2Ndof.at(i);
    delete hChi2Ndof2Gauss.at(i);
  }
  delete leg;
  delete label;
  delete can;
}


// ------------------------------------------------------------
void fitMCTruth() {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::setStyleNoteNoTitle();
  
  sampleTools::BinningAdmin admin("config/Analysis2011/Binning/BinningAdmin2011_v2.cfg");
  
  TString file = "Kalibri_MCTruthResponse_Summer11_AK5PF.root";
  TString name = "MCTruthSummer11";

  //plotMCTruth(file,"hRespVsPtGen",admin,2.,10.,name);
  //plotDeltaR(file,"hDeltaRVsPtGen",admin,100.,200.,name);
  //  plotGaussianFitRange(file,"hRespVsPtGen",admin,10.,name);
  //plotTwoGaussFit(file,"hRespVsPtGenUncorrected",admin,2.,10.,name);
  //plotMCTruth(file,"hRespVsPtGenUncorrected",admin,2.,10.,name+"_Uncorr");


   std::vector<TString> histNamePrefix;
   std::vector<TString> legEntries;

//    histNamePrefix.push_back("hRespVsPtGenUncorrected");
//    histNamePrefix.push_back("hRespVsPtGen");
//    legEntries.push_back("Uncorrected");
//    legEntries.push_back("C_{off} #upoint C_{rel} #upoint C_{abs}");
//    plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_UncorrVsCorr",true,false,1);

//    histNamePrefix.clear();
//    legEntries.clear();
//    histNamePrefix.push_back("hRespVsPtGenUncorrected");
//    histNamePrefix.push_back("hRespVsPtGenL1Corrected");
//    histNamePrefix.push_back("hRespVsPtGen");
//    legEntries.push_back("Uncorrected");
//    legEntries.push_back("C_{off}");
//    legEntries.push_back("C_{off} #upoint C_{rel} #upoint C_{abs}");
//    plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_UncorrVsL1VsCorr");

//    histNamePrefix.clear();
//    legEntries.clear();
//    histNamePrefix.push_back("hRespVsPtGenUncorrected");
//    histNamePrefix.push_back("hRespVsPtGenUncorrectedLowPU");
//    histNamePrefix.push_back("hRespVsPtGenL1Corrected");
//    legEntries.push_back("Uncorrected  ");
//    legEntries.push_back("N(PU) < 2");
//    legEntries.push_back("C_{off}");
//    plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_UncorrVsLowPUVsL1Corr",true);

//    histNamePrefix.clear();
//    legEntries.clear();
//    histNamePrefix.push_back("hRespVsPtGenVsPdgId_Eta?_PdgId1:hRespVsPtGenVsPdgId_Eta?_PdgId2:hRespVsPtGenVsPdgId_Eta?_PdgId3");
//    histNamePrefix.push_back("hRespVsPtGenVsPdgId_Eta?_PdgId4");
//    histNamePrefix.push_back("hRespVsPtGenVsPdgId_Eta?_PdgId5");
//    histNamePrefix.push_back("hRespVsPtGenVsPdgId_Eta?_PdgId21");
//    legEntries.push_back("u,d,s");
//    legEntries.push_back("c");
//    legEntries.push_back("b");
//    legEntries.push_back("g");
//    plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_Flavour",true);

   histNamePrefix.clear();
   legEntries.clear();
   histNamePrefix.push_back("hRespVsPtGen");
   histNamePrefix.push_back("hRespVsPtGen");
   legEntries.push_back("Particle-flow");
   legEntries.push_back("Calorimeter");
   std::vector<TString> fileNames;
   fileNames.push_back(file);
   fileNames.push_back("Kalibri_MCTruthResponse_Summer11_AK5Calo.root");
   plotMCTruth(fileNames,histNamePrefix,legEntries,admin,2.,10.,name+"_PFVsCalo");

//      histNamePrefix.clear();
//      legEntries.clear();
//       histNamePrefix.push_back("hRespVsPtGenDeltaRLess05");
//       histNamePrefix.push_back("hRespVsPtGenDeltaRLess10");
//       histNamePrefix.push_back("hRespVsPtGenDeltaRLess15");
//       histNamePrefix.push_back("hRespVsPtGenDeltaRLess20");
//       histNamePrefix.push_back("hRespVsPtGenDeltaRLess25");
//       //  histNamePrefix.push_back("hRespVsPtGenDeltaRLess30");
//       legEntries.push_back("#DeltaR_{max} = 0.05");
//       legEntries.push_back("#DeltaR_{max} = 0.10");
//       legEntries.push_back("#DeltaR_{max} = 0.15");
//       legEntries.push_back("#DeltaR_{max} = 0.20");
//       legEntries.push_back("#DeltaR_{max} = 0.25");
//       //  legEntries.push_back("#DeltaR_{max} = 0.30");
//       plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_ResolutionVsDeltaR");

  //     histNamePrefix.clear();
  //     legEntries.clear();
  //     histNamePrefix.push_back("hRespVsPtGen");
  //     histNamePrefix.push_back("hRespVsPtGenNoJetID");
  //     legEntries.push_back("JetID applied");
  //     legEntries.push_back("JetID not applied");
  //     plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_ResolutionJetID");

//        histNamePrefix.clear();
//        legEntries.clear();
//        histNamePrefix.push_back("hRespVsPtGen");
//        histNamePrefix.push_back("hRespVsPtGenPUMC");
//        legEntries.push_back("With PU reweighting");
//        legEntries.push_back("No PU reweighting");
//        plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_ResolutionPU");

//        histNamePrefix.clear();
//        legEntries.clear();
//        histNamePrefix.push_back("hRespVsPtGen");
//        //   histNamePrefix.push_back("hRespVsPtGenVsPtSoft_Eta?_PtSoft4");
//        histNamePrefix.push_back("hRespVsPtGenVsPtSoft_Eta?_PtSoft6");
//        legEntries.push_back("QCD multijets");
//        legEntries.push_back("QCD dijets");
//        plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_ResolutionDijetSelection",true);
 
//      histNamePrefix.clear();
//      legEntries.clear();
//      histNamePrefix.push_back("hRespVsPtGenPULess05");
//      histNamePrefix.push_back("hRespVsPtGenPULess10");
//      histNamePrefix.push_back("hRespVsPtGenPULess15");
//      histNamePrefix.push_back("hRespVsPtGenPULess99");
//      legEntries.push_back("N(PU) #leq 4");
//      legEntries.push_back("4 < N(PU) #leq 9");
//      legEntries.push_back("9 < N(PU) #leq 14");
//      legEntries.push_back("N(PU) > 14");
//      plotMCTruth(file,histNamePrefix,legEntries,admin,2.,10.,name+"_ResolutionNPU");
}


// ------------------------------------------------------------
void plotAsymmetryComponents(const TString &file) {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::setStyleNoteNoTitle();

  TH1* hSigAsym = util::FileOps::readTH1(file,"Components:hSigAsym");
  TH1* hSigResp = util::FileOps::readTH1(file,"Components:hSigResp");
  TH1* hSigSoft = util::FileOps::readTH1(file,"Components:hSigSoft");
  TH1* hSigTotal = util::FileOps::readTH1(file,"Components:hSigTotal");

  hSigAsym->SetMarkerStyle(21);
  hSigAsym->SetMarkerColor(kBlack);
  hSigAsym->SetLineColor(kBlack);
  hSigAsym->SetLineWidth(2);

  hSigTotal->SetLineColor(kRed);
  hSigTotal->SetLineWidth(2);
  hSigTotal->GetYaxis()->SetRangeUser(3E-4,0.59);
  hSigTotal->GetXaxis()->SetRange(1,hSigTotal->GetNbinsX()-1);
  hSigTotal->GetXaxis()->SetMoreLogLabels();
  hSigTotal->GetXaxis()->SetNoExponent();
  util::HistOps::setAxisTitles(hSigTotal,"p^{gen}_{T}","GeV","Standard deviation");

  hSigResp->SetLineColor(kBlue);
  hSigResp->SetLineWidth(2);
  hSigResp->SetLineStyle(2);

  hSigSoft->SetLineColor(28);
  hSigSoft->SetLineWidth(2);
  hSigSoft->SetLineStyle(3);

  TPaveText* label = util::LabelFactory::createPaveText(3,-0.65);
  label->AddText("CMS Simulation,  #sqrt{s} = 7 TeV");
  label->AddText(jetLabel_);
  label->AddText("0.0 < |#eta^{gen}| < 0.5, |#Delta#phi^{gen}| > 2.7");
  TLegend* leg = util::LabelFactory::createLegendCol(4,0.35);
  leg->AddEntry(hSigAsym,"#sqrt{2} #upoint #sigma_{A}","P");
  leg->AddEntry(hSigResp,"#sigma/p_{T}","L");
  leg->AddEntry(hSigSoft,"#sigma_{soft}","L");
  leg->AddEntry(hSigTotal,"#sigma/p_{T} #oplus #sigma_{soft}","L");

  
  TCanvas* can = new TCanvas("can","Asymmetry Components",500,500);
  can->cd();
  hSigTotal->Draw("HIST");
  hSigSoft->Draw("HISTEsame");
  hSigResp->Draw("HISTEsame");
  hSigAsym->Draw("PE1same");
  label->Draw("same");
  leg->Draw("same");
  can->SetLogx();
  gPad->RedrawAxis();
  can->SaveAs("MCTruthSummer11_AsymmetryContributions.eps","eps");
}

