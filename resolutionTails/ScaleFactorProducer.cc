// $Id: ScaleFactorProducer.cc,v 1.3 2013/05/10 14:25:50 mschrode Exp $

#ifndef RESOLUTION_TAILS_SCALE_FACTOR_PRODUCER
#define RESOLUTION_TAILS_SCALE_FACTOR_PRODUCER

#include <vector>
#include <cstdlib>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "EtaPtBin.cc"
#include "FitParameters.h"
#include "LittleHelper.h"
#include "Output.h"
#include "Style.h"
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"



namespace resolutionTails {
  class ScaleFactorProducer {
  public:
    ScaleFactorProducer(const TString &fileNameData, const TString &fileNameMC, const TString &fileNameTailWindow, const resolutionTails::FitParameters* fitPars, const std::vector<double> &coreScaleFactors, const sampleTools::BinningAdmin* binAdm, Output* out, const Style* style);
    ~ScaleFactorProducer();

    void makeScaleFactors(bool useExtrapolatedValue, bool mcTruthRef);
    void makeControlPlots() const;

    void writeWindowBorders() const;
    void writeMCClosure() const;
    void writeExtrapolation() const;


  private:
    typedef std::vector<EtaPtBin*>::const_iterator EtaPtBinConstIt;

    const sampleTools::BinningAdmin* binAdm_;
    const resolutionTails::FitParameters* fitPars_;
    const Style* style_;

    const TString fileNameData_;
    const TString fileNameMC_;
    const TString fileNameTailWindow_;

    Output* out_;

    std::vector<EtaPtBin*> etaPtBins_;

    bool completePicName(int padIdx, TString &picName) const;
    void makeControlPlotsBins() const;
    void makeControlPlotsExtrapolation() const;
    void makeControlPlotsScaleFactors() const;
    void makeControlPlotsTailStart() const;
    void makeControlPlotsFAsym() const;
    void writeLaTeXSlides() const;
  };


  // ------------------------------------------------------------------------------------
  ScaleFactorProducer::ScaleFactorProducer(const TString &fileNameData, const TString &fileNameMC, const TString &fileNameTailWindow, const resolutionTails::FitParameters* fitPars, const std::vector<double> &coreScaleFactors, const sampleTools::BinningAdmin* binAdm, Output* out, const Style* style)
    : binAdm_(binAdm), fitPars_(fitPars), style_(style), fileNameData_(fileNameData), fileNameMC_(fileNameMC), fileNameTailWindow_(fileNameTailWindow), out_(out) {

    if( coreScaleFactors.size() < binAdm_->nEtaBins() ) {
      std::cerr << "ERROR: core scale factors not defined for every eta bin" << std::endl;
      exit(-1);
    }

    // Initialize bins
    std::cout << "Creating bins" << std::endl;
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      std::cout << "  Eta " << etaBin << ": Using core scale factor " << coreScaleFactors.at(etaBin) << std::endl;
      
      for(unsigned int ptBin = 0; ptBin < binAdm_->nPtBins(etaBin); ++ptBin) {
	if( Output::DEBUG ) std::cout << "Setting up bin (" << etaBin << ", " << ptBin << ")" << std::endl;
	
	// Create eta-pt bin
	EtaPtBin* bin = new EtaPtBin(etaBin,ptBin,fitPars_->nSigCore(),coreScaleFactors.at(etaBin),binAdm_,out_,style_);

	// Define window
	if( fileNameTailWindow_ == "" ) {
	  bin->findWindow(fileNameMC_,0,fitPars_->nSigTailStart(),fitPars_->nSigTailEnd());
	} else {
	  bin->setWindow(fileNameTailWindow_,fitPars_->nSigTailStart(),fitPars_->nSigTailEnd());
	}

	// Add pt3 bins
	for(unsigned int ptSoftBin = 0; ptSoftBin < binAdm_->nPtSoftBins(); ++ptSoftBin) {
	  bin->addPt3Bin(ptSoftBin,binAdm_->ptSoftMax(ptSoftBin),fileNameData_,fileNameMC_);
	}
	
	// Add mc truth response for toy asymmetry
	bin->addMCTruthForToyAsym(fileNameMC_);
	
	etaPtBins_.push_back(bin);
	if( Output::DEBUG ) std::cout << "Done setting up bin" << std::endl;
      }
    }
  }


  // ------------------------------------------------------------------------------------
  ScaleFactorProducer::~ScaleFactorProducer() {
    for(std::vector<EtaPtBin*>::iterator it = etaPtBins_.begin();
	it != etaPtBins_.end(); ++it) {
      delete *it;
    }
  }



  // Perform extrapolation
  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeScaleFactors(bool useExtrapolatedValue, bool mcTruthRef) {
    std::cout << "Performing extrapolation vs pt3rel" << std::endl;
    for(EtaPtBinConstIt it = etaPtBins_.begin();	
	it != etaPtBins_.end(); ++it ) {
      if( Output::DEBUG ) std::cout << "Extrapolation" << std::endl;
      
      // Extrapolation
      (*it)->extrapolate(fitPars_->minPt3Data(),fitPars_->fixDataShape(),useExtrapolatedValue,mcTruthRef);
    }
  }


  // Make all control plots
  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeControlPlots() const {
    std::cout << "Creating control plots" << std::endl;
    makeControlPlotsBins();
    makeControlPlotsExtrapolation();
    makeControlPlotsScaleFactors();
    makeControlPlotsTailStart();
    makeControlPlotsFAsym();
    writeLaTeXSlides();
  }



  // Make control plots
  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeControlPlotsBins() const {
    std::cout << "  Plots per pt3Bin" << std::endl;
    for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins_.begin();
	it != etaPtBins_.end(); ++it ) {
      if( Output::DEBUG ) std::cout << "ScaleFactorProducer::makeControlPlotsBins(): eta " << (*it)->etaBin() << ", pt " << (*it)->ptBin() << std::endl;
      
      (*it)->plotAsymmetryDistributions(fitPars_->nSigTailStart(),fitPars_->nSigTailEnd());
      (*it)->plotMCTruth();
      (*it)->plotSpectra();
      (*it)->plotTailStart();
    }
  }



  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeControlPlotsExtrapolation() const {
    std::cout << "  Extrapolation plots" << std::endl;
    for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins_.begin();
	it != etaPtBins_.end(); ++it ) {
      if( Output::DEBUG ) std::cout << "ScaleFactorProducer::makeControlPlotsExtrapolation(): eta " << (*it)->etaBin() << ", pt " << (*it)->ptBin() << std::endl;
      (*it)->plotExtrapolation();
    }
  }



  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeControlPlotsScaleFactors() const {
    std::cout << "  Scale factor summary plots" << std::endl;
    // Loop over eta bins
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      if( Output::DEBUG ) std::cout << "ScaleFactorProducer::makeControlPlotsScaleFactors(): eta " << etaBin << std::endl;

      //// Create histograms
      // Averave ptAve per bin in data per ptAve bin in this eta bin
      TH1* hPtAveMeanData = new TH1D("hPtAveMeanData_Eta"+util::toTString(etaBin),
				     style_->title()+";p^{ave}_{T} Bin Idx;<p^{ave}_{T}> (GeV)",
				     binAdm_->nPtBins(etaBin),-0.5,binAdm_->nPtBins(etaBin)-0.5);
      hPtAveMeanData->SetMarkerStyle(20);

      // Difference of extrapolated tail in data and MC as a function of
      // average ptAve in this eta bin
      TH1* hDelta = new TH1D("hDelta_Eta"+util::toTString(etaBin),
			     style_->title()+";p^{ave}_{T} (GeV);#Delta = "+style_->labelFAsymData()+"(0) - "+style_->labelFAsymMC()+"(0)",
			     binAdm_->nPtBins(etaBin),&(binAdm_->ptBinEdges(etaBin).front()));
      hDelta->SetMarkerStyle(20);

      // Scale factors as a function of average ptAve in this eta bin
      TH1* hScale = new TH1D("hScale_Eta"+util::toTString(etaBin),
			     style_->title()+";p^{ave}_{T} (GeV);"+style_->labelScaleFactor(),
			     binAdm_->nPtBins(etaBin),&(binAdm_->ptBinEdges(etaBin).front()));
      hScale->SetMarkerStyle(20);
      hScale->GetXaxis()->SetMoreLogLabels();
      hScale->GetXaxis()->SetNoExponent();


      //// Loop over pt bins in this eta bin
      for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins_.begin();
	  it != etaPtBins_.end(); ++it) {
	if( (*it)->etaBin() == etaBin ) {
	  int bin = 1+(*it)->ptBin();
	  hPtAveMeanData->SetBinContent(bin,(*it)->ptAveMeanData(0));
	  hPtAveMeanData->SetBinError(bin,(*it)->ptAveMeanDataErr(0));
	  hDelta->SetBinContent(bin,(*it)->deltaExtra());
	  hDelta->SetBinError(bin,(*it)->deltaExtraErr());
	  hScale->SetBinContent(bin,(*it)->scalingFactor());
	  hScale->SetBinError(bin,(*it)->scalingFactorErr());
	}
      } // End of loop over pt bins

      TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
    
      TH1* hDeltaFrame = static_cast<TH1D*>(hDelta->Clone("hDeltaFrame_Eta"+util::toTString(etaBin)));
      hDeltaFrame->SetTitle(style_->title());
      hDeltaFrame->SetLineColor(kBlack);
      hDeltaFrame->SetLineStyle(2);
      hDeltaFrame->SetMarkerStyle(1);
      for(int i = 1; i <= hDeltaFrame->GetNbinsX(); ++i) {
	hDeltaFrame->SetBinContent(i,0.);
	hDeltaFrame->SetBinError(i,0.);
      }
      hDeltaFrame->GetYaxis()->SetRangeUser(-0.0045,0.013);
      hDeltaFrame->GetXaxis()->SetMoreLogLabels();
      can->cd();
      hDeltaFrame->Draw("HIST");
      hDelta->Draw("PE1same");
      can->SetLogx();
      out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_Delta");
      
      TH1* hScaleFrame = static_cast<TH1D*>(hScale->Clone("hScaleFrame_Eta"+util::toTString(etaBin)));
      hScaleFrame->SetLineColor(kBlack);
      hScaleFrame->SetLineStyle(2);
      hScaleFrame->SetMarkerStyle(1);
      for(int i = 1; i <= hScaleFrame->GetNbinsX(); ++i) {
	hScaleFrame->SetBinContent(i,1.);
	hScaleFrame->SetBinError(i,0.);
      }
      hScaleFrame->GetYaxis()->SetRangeUser(0.01,2.99);
      hScaleFrame->GetXaxis()->SetMoreLogLabels();
      can->cd();
      hScaleFrame->Draw("HIST");
      hScale->Draw("PE1same");
      can->SetLogx();
      TPaveText* label = util::LabelFactory::createPaveText(1);
      label->AddText(style_->labelLumi()+",  "+style_->labelWindow(fitPars_->nSigTailStart(),fitPars_->nSigTailEnd())+",  "+util::LabelFactory::etaCut(binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin)));
      label->Draw("same");
      out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_ScaleFactors");
      
      out_->toFiles(hDelta);    
      out_->toFiles(hDeltaFrame);   
      out_->toFiles(hScale);    
      out_->toFiles(hScaleFrame);   
      out_->toFiles(hPtAveMeanData);
      
      delete hDelta;
      delete hScale;
      delete hDeltaFrame;
      delete hScaleFrame;
      delete hPtAveMeanData;
      delete label;
      delete can;
    } // End of loop over eta bins
  }



  // Plot the difference between the nominal tail start, defined in
  // multiples of sigmaAsym, and the effective tail start due to the
  // binning of the asymmetry histograms
  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeControlPlotsTailStart() const {
    std::cout << "  Tail start plots" << std::endl;
    //// Loop over eta bins
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      if( Output::DEBUG ) std::cout << "ScaleFactorProducer::makeControlPlotsTailStart(): eta " << etaBin << std::endl;

      //// Create histograms
      // Nominal start of tail region in asymmetry (same A for all pt3 bins)
      // for each pt bin in this eta bin
      TH1* hAsymTailStart = new TH1D("hAsymTailStart_EtaBin"+util::toTString(etaBin),style_->title()+";"+util::LabelFactory::ptAve()+" (GeV);Asymmetry at Tail Start",binAdm_->nPtBins(etaBin),&(binAdm_->ptBinEdges(etaBin).front()));
      hAsymTailStart->SetLineWidth(style_->lineWidth());
      hAsymTailStart->GetXaxis()->SetMoreLogLabels();
      hAsymTailStart->GetXaxis()->SetNoExponent();

      // Effective start of tail region in asymmetry (due to binning effects)
      // for each pt bin in this eta bin
      TH1* hAsymTailStartEff = static_cast<TH1*>(hAsymTailStart->Clone("hAsymTailStartEff_EtaBin"+util::toTString(etaBin)));
      hAsymTailStartEff->SetLineStyle(2);
      hAsymTailStartEff->SetLineColor(kBlue);


      //// Loop over pt bins in this eta bin
      for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins_.begin();
	  it != etaPtBins_.end(); ++it) {
	if( (*it)->etaBin() == etaBin ) {
	  int bin = 1+(*it)->ptBin();
	  hAsymTailStart->SetBinContent(bin,(*it)->tailWindowMin());
	  hAsymTailStartEff->SetBinContent(bin,(*it)->tailWindowMinEff());
	}
      }	// End of loop over pt bins

      // Labels
      TPaveText* label = util::LabelFactory::createPaveText(2);
      label->AddText(style_->labelMCSmear());
      label->AddText(util::LabelFactory::etaCut(binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin))+",  "+style_->labelWindow(fitPars_->nSigTailStart(),fitPars_->nSigTailEnd()));
      
      TLegend* leg = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
      leg->AddEntry(hAsymTailStart,style_->labelAsymTailStart(),"L");
      leg->AddEntry(hAsymTailStartEff,style_->labelAsymTailStartEff(),"L");

      // Plot
      TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
      can->cd();
      hAsymTailStart->GetYaxis()->SetRangeUser(0.071,0.169);
      hAsymTailStart->Draw("HIST");
      hAsymTailStartEff->Draw("HISTsame");
      hAsymTailStart->Draw("HISTsame");
      leg->Draw("same");
      label->Draw("same");
      gPad->RedrawAxis();
      can->SetLogx();
      out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_AsymTailStart");
      out_->toFiles(hAsymTailStart);
      out_->toFiles(hAsymTailStartEff);

      delete hAsymTailStart;
      delete hAsymTailStartEff;
      delete label;
      delete leg;
      delete can;
    } // End of loop over eta bins
  }



  // Plot fractional number of tail events
  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::makeControlPlotsFAsym() const {
    std::cout << "  Tail size plots" << std::endl;

    // Find min and max of fasym per pt3Bin
    std::vector<double> minFAsym(etaPtBins_.front()->nPt3Bins(),1000.);
    std::vector<double> maxFAsym(etaPtBins_.front()->nPt3Bins(),0.);
    for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins_.begin(); it != etaPtBins_.end(); ++it) {
      for(unsigned int i = 0; i < (*it)->nPt3Bins(); ++i) {
	if( (*it)->fTailMCSmearedGauss(i) < minFAsym.at(i) ) 
	  minFAsym.at(i) = (*it)->fTailMCSmearedGauss(i);
	if( (*it)->fTailData(i) > maxFAsym.at(i) )
	  maxFAsym.at(i) = (*it)->fTailData(i);
      }
    }
    
    // Loop over eta bins
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      if( Output::DEBUG ) std::cout << "ScaleFactorProducer::makeControlPlotsFAsym(): eta " << etaBin << std::endl;

      // Loop over pt3 bins
      for(unsigned int pt3Bin = 0; pt3Bin < binAdm_->nPtSoftBins(); ++pt3Bin) {
	  
	TH1* hFAsymMCSmeared = new TH1D("hFAsymMCSmeared_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin),style_->title()+";p^{ave}_{T} (GeV);"+style_->labelFAsym()+"  (%)",binAdm_->nPtBins(etaBin),&(binAdm_->ptBinEdges(etaBin).front()));
	hFAsymMCSmeared->GetXaxis()->SetMoreLogLabels();
	hFAsymMCSmeared->GetYaxis()->SetMoreLogLabels();
	hFAsymMCSmeared->GetXaxis()->SetNoExponent();
	hFAsymMCSmeared->GetYaxis()->SetNoExponent();
	hFAsymMCSmeared->GetYaxis()->SetNdivisions(505);
	hFAsymMCSmeared->SetLineWidth(style_->lineWidth());
	hFAsymMCSmeared->SetLineColor(style_->colorFilledAsymSmear());
	hFAsymMCSmeared->SetLineStyle(1);

	TH1* hFAsymData = static_cast<TH1*>(hFAsymMCSmeared->Clone("hFAsymData_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)));
	hFAsymData->SetMarkerStyle(20);
	hFAsymData->SetMarkerSize(style_->markerSize());
	hFAsymData->SetLineColor(kBlack);
      
	TH1* hFAsymMCSmearedGauss = static_cast<TH1*>(hFAsymMCSmeared->Clone("hFAsymMCSmearedGauss_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)));
	hFAsymMCSmearedGauss->SetLineColor(style_->colorGauss());
	hFAsymMCSmearedGauss->SetLineStyle(2);

	// Averave ptAve per bin in data per ptAve bin in this eta bin
	TH1* hPtAveMeanDataTmp = new TH1D("hPtAveMeanDataTmp","",binAdm_->nPtBins(etaBin),-0.5,binAdm_->nPtBins(etaBin)-0.5);
	
	// Loop over pt bins and fill histograms
	for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins_.begin();
	    it != etaPtBins_.end(); ++it) {
	  if( (*it)->etaBin() == etaBin ) {
	    int bin = 1+(*it)->ptBin();
	    hFAsymData->SetBinContent(bin,100.*(*it)->fTailData(pt3Bin));
	    hFAsymData->SetBinError(bin,100.*(*it)->fTailDataErr(pt3Bin));
	    hFAsymMCSmeared->SetBinContent(bin,100.*(*it)->fTailMCSmeared(pt3Bin));
	    hFAsymMCSmeared->SetBinError(bin,100.*(*it)->fTailMCSmearedErr(pt3Bin));
	    hFAsymMCSmearedGauss->SetBinContent(bin,100.*(*it)->fTailMCSmearedGauss(pt3Bin));
	    hPtAveMeanDataTmp->SetBinContent(bin,(*it)->ptAveMeanData(0));
	    hPtAveMeanDataTmp->SetBinError(bin,(*it)->ptAveMeanDataErr(0));
	  }
	} // End of loop over pt bins
	  
	  // Write basic histograms to ROOT file for further processing
	out_->toFiles(hFAsymData);
	out_->toFiles(hFAsymMCSmeared);
	out_->toFiles(hFAsymMCSmearedGauss);

	// Convert fAsymData to graph and set x value to mean
	// of ptAve spectrum in each bin
	TGraphAsymmErrors* gFAsymData = util::HistOps::createTGraph(hFAsymData);
	for(int n = 0; n < gFAsymData->GetN(); ++n) {
	  gFAsymData->GetX()[n] = hPtAveMeanDataTmp->GetBinContent(1+n);
	  gFAsymData->SetPointEXhigh(n,hPtAveMeanDataTmp->GetBinError(1+n));
	  gFAsymData->SetPointEXlow(n,hPtAveMeanDataTmp->GetBinError(1+n));
	}
	delete hPtAveMeanDataTmp;
      
	// Convert fAsymMC to graph to be able to plot error band
	TGraphAsymmErrors* gFAsymMCSmeared = util::HistOps::getUncertaintyBand(hFAsymMCSmeared,style_->colorFilledAsymSmear());
	gFAsymMCSmeared->SetLineColor(gFAsymMCSmeared->GetFillColor());
	gFAsymMCSmeared->SetLineStyle(1);
	gFAsymMCSmeared->SetLineWidth(style_->lineWidth());
	gFAsymMCSmeared->SetFillStyle(style_->hatchStyle());

	// Ratios to Gauss
	TGraphAsymmErrors* gFAsymDataRelToGauss = util::HistOps::createRatioGraph(gFAsymData,hFAsymMCSmearedGauss);
	TGraphAsymmErrors* gFAsymMCSmearedRelToGauss = util::HistOps::createRatioGraph(gFAsymMCSmeared,hFAsymMCSmearedGauss);
	gFAsymMCSmearedRelToGauss->SetFillStyle(gFAsymMCSmeared->GetFillStyle());
	gFAsymMCSmearedRelToGauss->SetLineColor(gFAsymMCSmeared->GetLineColor());
	gFAsymMCSmearedRelToGauss->SetLineWidth(style_->lineWidth());
	TH1* hFAsymMCSmearedRelToGauss = util::HistOps::createRatioPlot(hFAsymMCSmeared,hFAsymMCSmearedGauss);
	  
	// Labels
	TPaveText* label = util::LabelFactory::createPaveText(2);
	label->AddText(util::LabelFactory::etaCut(binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin))+",  "+style_->labelWindow(fitPars_->nSigTailStart(),fitPars_->nSigTailEnd()));
	label->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+util::LabelFactory::pt3RelCut(binAdm_->ptSoftMax(pt3Bin)));
	  
	TLegend* leg = util::LabelFactory::createLegendWithOffset(3,label->GetSize());
	leg->AddEntry(gFAsymData,style_->labelData(),"P");
	leg->AddEntry(gFAsymMCSmeared,style_->labelMCSmear(),"LF");
	leg->AddEntry(hFAsymMCSmearedGauss,style_->labelMCGauss(),"L");
	  
	// Plot
	TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
	can->cd();
	double min = 100.*minFAsym.at(pt3Bin);
	double max = 100.*maxFAsym.at(pt3Bin);
	double delta = max - min;
	hFAsymMCSmeared->GetYaxis()->SetRangeUser(0.9*min,max+1.2*delta);
	hFAsymMCSmeared->Draw("HIST");
	hFAsymMCSmearedGauss->Draw("HISTsame");
	gFAsymMCSmeared->Draw("E2same");
	hFAsymMCSmeared->Draw("HISTsame");
	gFAsymData->Draw("PE1same");
	leg->Draw("same");
	label->Draw("same");
	gPad->RedrawAxis();
	can->SetLogx();
	out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsym");
	  
	// Plot with ratio relative to Gauss
	delete can;
	can = util::HistOps::createRatioTopCanvas();
	TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
	//TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hFAsymData);
	TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hFAsymMCSmeared);
	TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hFAsymMCSmeared,util::LabelFactory::ptAve(),"GeV","#frac{Data, Sim.}{Gauss}",0.81,2.99);
	bRatioBottomFrame->GetYaxis()->SetRangeUser(0.95,1.1*(*std::max_element(gFAsymDataRelToGauss->GetY(),gFAsymDataRelToGauss->GetY()+gFAsymDataRelToGauss->GetN())));
	bRatioBottomFrame->SetLineStyle(hFAsymMCSmearedGauss->GetLineStyle());
	bRatioBottomFrame->SetLineColor(hFAsymMCSmearedGauss->GetLineColor());
	bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
	bRatioBottomFrame->GetXaxis()->SetNoExponent();
	bRatioTopFrame->GetYaxis()->SetRangeUser(0.9*min,max+1.2*delta);
	can->cd();
	bRatioTopFrame->Draw("HIST");
	hFAsymMCSmearedGauss->Draw("HISTsame");
	gFAsymMCSmeared->Draw("E2same");
	hFAsymMCSmeared->Draw("HISTsame");
	gFAsymData->Draw("PE1same");
	leg->Draw("same");
	label->Draw("same");
	gPad->RedrawAxis();
	can->SetLogx();
	bRatioBottomPad->Draw();
	bRatioBottomPad->cd();
	bRatioBottomFrame->Draw("HIST");
	gFAsymMCSmearedRelToGauss->Draw("E2same");
	hFAsymMCSmearedRelToGauss->Draw("HISTsame");
	gFAsymDataRelToGauss->Draw("PE1same");
	bRatioBottomPad->SetLogx();
	bRatioBottomPad->RedrawAxis();
	out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsymRelToGaussBottom");
	  
	delete hFAsymData;
	delete hFAsymMCSmeared;
	delete bRatioBottomPad;
	delete bRatioTopFrame;
	delete bRatioBottomFrame;
	delete label;
	delete leg;
	delete can;
      } // End of loop over pt3 bins
    }
  }



  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::writeWindowBorders() const {
    char tmp[50];
    
    TString text = "% ***  Window Borders ("+out_->namePrefix()+")  ***\n\n";
    text += "\\begin{tabular}{cr@{ -- }rcc}\n\\toprule\n";
    text += "\\multicolumn{3}{c}{Interval} & \\multicolumn{2}{c}{\\tailborder{";
    text += util::toTString(fitPars_->nSigTailStart(),1,true);
    text += "}} \\\\\n";
    text += "$|\\eta|$ & \\multicolumn{2}{c}{$\\ptave \\,(\\gevnospace)$} & \\asymtaileff & $\\asymtaileff/\\sigma_{\\asym}$ \\\\\n\\midrule\n";
    for(EtaPtBinConstIt it = etaPtBins_.begin(); it != etaPtBins_.end(); ++it) {
      EtaPtBin* bin = *it;
      text += bin->toTexTableCells();
      
      double min = bin->tailWindowMinEff();
      double sigSmear = bin->sigmaSmeared(0);
      sprintf(tmp," & $%.3lf$ & $%.3lf$ \\\\\n",min,min/sigSmear);
      text += tmp;
      if( bin->ptBin() == binAdm_->nPtBins(bin->etaBin())-1 ) text += "\\midrule\n";
    }
    text += "\\bottomrule\n\\end{tabular}\n\n\n\n";

    std::cout << "Writing window-border info to " << out_->toFile(text,"Info.tex") << std::endl;
  }


  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::writeMCClosure() const {
    // Print only in case of closure test...
    if( etaPtBins_.front()->hasToyMC() ) {
      char tmp[50];
      TString text = "% ***  MC Closure ("+out_->namePrefix()+")  ***\n\n";
      text += "\\begin{tabular}[ht]{ccccc}\n\\hline\n";
      text += "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\fasymmc(0)$ & \\fasymtoy & $\\fasymmc(0))/\\fasymtoy$ \\\\ \n\\hline \n";
      for(EtaPtBinConstIt it = etaPtBins_.begin(); it != etaPtBins_.end(); ++it) {
	EtaPtBin* bin = *it;
	text += bin->toTexTableCells();

	sprintf(tmp," & $%.3lf \\pm %.3lf$",bin->extraMC(),bin->extraMCErr());
	text += tmp;
	sprintf(tmp," & $%.3lf \\pm %.3lf$",bin->toyMC(),bin->toyMCErr());
	text += tmp;
	sprintf(tmp," & $%.3lf$ \\\\ \n",bin->extraMC()/bin->toyMC());

	if( bin->ptBin() == binAdm_->nPtBins(bin->etaBin())-1 ) text += "\\hline\n";
      }
      text += "  \\end{tabular}\n\n\n\n";

      std::cout << "Writing MC-closure info to " << out_->toFile(text,"Info.tex") << std::endl;
    }
  }


  // ------------------------------------------------------------------------------------
  void ScaleFactorProducer::writeExtrapolation() const {
    char tmp[50];
    TString text = "% ***  Extrapolation ("+out_->namePrefix()+")  ***\n";
    text += "\\begin{tabular}[ht]{ccccc}\n\\toprule\n";
    text += "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\mean{\\ptave} \\,(\\gevnospace)$ & $\\fasymdata(0)$ & $\\fasymmc(0)$ \\\\ \n\\midrule\n";
    for(EtaPtBinConstIt it = etaPtBins_.begin(); it != etaPtBins_.end(); ++it) {
      EtaPtBin* bin = *it;
      text += bin->toTexTableCells();

      sprintf(tmp," & $%.1lf \\pm %.1lf$",bin->ptAveMeanData(0),bin->ptAveMeanDataErr(0));
      text += tmp;
      sprintf(tmp," & $%.3lf \\pm %.3lf$",bin->extraData(),bin->extraDataErr());
      text += tmp;
      sprintf(tmp," & $%.3lf \\pm %.3lf$ \\\\ \n",bin->extraMC(),bin->extraMCErr());
      text += tmp;

      if( bin->ptBin() == binAdm_->nPtBins(bin->etaBin())-1 ) text += "\\midrule\n";
    }
    text += "  \\end{tabular}\n\n\n\n";

    std::cout << "Writing extrapolation-closure info to " << out_->toFile(text,"Info.tex") << std::endl;
  }


  // Create beamer-LaTeX source code for slides listing
  // created control plots.
  // -----------------------------------------------------------------------------------
  void ScaleFactorProducer::writeLaTeXSlides() const {
    char tmp[200];
    TString text = "";
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      for(unsigned int ptBin = 0; ptBin < binAdm_->nPtBins(etaBin); ++ptBin) {
	text += "\n\n\n% ==== Eta ";
	text += etaBin;
	text += ", Pt ";
	text += ptBin;
	text += " ================================================\n";
	
	int nPages = 1 + (1+2*binAdm_->nPtSoftBins())%12;
	for(int page = 0; page < nPages; ++page) {
	  text += "\n% -----------------------------------------------------------------\n";
	  sprintf(tmp,"\\begin{frame}\\frametitle{Asymmetry Tails $%.1lf < |\\eta| < %.1lf$}\n",binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin));
	  text += "  \\begin{columns}[T]\n";
	  for(int colIdx = 0; colIdx < 4; ++colIdx) {
	    text += "    \\begin{column}{0.25\\textwidth}\n";
	    text += "    \\centering\n";
	    for(int rowIdx = 0; rowIdx < 3; ++rowIdx) {
	      int padIdx = 12*page + 4*rowIdx + colIdx;
	      TString picName = out_->namePrefix()+"_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
	      if( completePicName(padIdx,picName) ) {
		text += "      \\includegraphics[width=\\textwidth]{figures/";
		text += picName;
		text += "}\\\\\n";
	      } else {
		continue;
	      }
	    }
	    text += "    \\end{column}\n";
	  }
	  text += "  \\end{columns}\n";
	  text += "\\end{frame}\n";
	}
      }
    }

    std::cout << "Writing LaTeX slides of control plots to file " << out_->toFile(text,".tex") << std::endl;
  }



  // ------------------------------------------------------------------------------------
  bool ScaleFactorProducer::completePicName(int padIdx, TString &picName) const {
    bool picExists = true;
    if( padIdx == 0 )
      picName += "_Extrapolation.pdf";
    else if( padIdx-1 < static_cast<int>(binAdm_->nPtSoftBins()) )
      picName += "_Pt3Bin"+util::toTString(padIdx-1)+"_PtSmearAsymTail.pdf";
    else if( padIdx-1-binAdm_->nPtSoftBins() < binAdm_->nPtSoftBins() )
      picName += "_Pt3Bin"+util::toTString(padIdx-1-binAdm_->nPtSoftBins())+"_PtAveSpectra.pdf";
    else
      picExists = false;
    
    return picExists;
  }
}
#endif
