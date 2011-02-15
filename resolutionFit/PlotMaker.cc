// $Id: $

#include "PlotMaker.h"

#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"

#include "FitResult.h"

#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  PlotMaker::PlotMaker(const Parameters *par, const EtaBins &etaBins)
    : par_(par), etaBins_(etaBins),
      xMinPt_(40.), xMaxPt_(1100.),
      yMinExtraRes_(1E-3), yMaxExtraRes_(0.38) {

    // Output manager
    out_ = OutputManager::createOutputManager(par_->outMode(),par_->outFilePrefix());

    // Style
    labelMk_ = new LabelMaker(par_);
  }


  // -------------------------------------------------------------------------------------
  PlotMaker::~PlotMaker() {
    delete out_;
    delete labelMk_;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotExtrapolation() const {
     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotExtrapolation(): Entering" << std::endl;
     std::cout << "Plotting extrapolation" << std::endl;
     
     // Loop over eta and pt bins
     for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
       out_->newPage("Extrapolation");

       for(PtBinIt ptBinIt = (*etaBinIt)->ptBinsBegin(); ptBinIt != (*etaBinIt)->ptBinsEnd(); ++ptBinIt) {
	 const PtBin* ptBin = *ptBinIt;
	 if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
       
	 // Loop over Samples
	 for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
	   const Sample* sample = sIt->second;
	   if( par_->verbosity() > 1 ) std::cout << "    " << sample->label() << std::endl;

	   // Loop over FitResultTypes
	   for(FitResultTypeIt rIt = (*etaBinIt)->fitResultTypesBegin(); 
	       rIt != (*etaBinIt)->fitResultTypesEnd(); ++rIt) {

	     if( par_->verbosity() > 1 ) std::cout << "      " << FitResult::toString(*rIt) << std::endl;
	 
	     // Graph of fitted resolutions
	     std::vector<double> val;
	     std::vector<double> uncert;
	     sample->values(*rIt,val,uncert);
	     std::vector<double> ptSoftx;
	     sample->ptSoft(ptSoftx);
	     std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
	     
	     TGraphAsymmErrors* gVals = 
	       new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
				     &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
				     &(uncert.front()),&(uncert.front()));
	     sample->setStyle(gVals);
	     
	     // Extrapolation function
	     TF1* fit = sample->extrapolationFunction(*rIt,"ExtrapolationForPlot");
	     fit->SetLineColor(kRed);

	     // Create frame
	     TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame","",
				    1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	     util::HistOps::setAxisTitles(hFrame,"p_{||,3} / p^{ref}_{T}","","#sigma / p^{ref}_{T}");
	     hFrame->GetYaxis()->SetRangeUser(0.8*val.front(),1.1*val.back());
	     
	     // Label
	     //TPaveText* label = labelMk_->binLabel(bin);
	     //TLegend* leg = util::LabelFactory::createLegendCol(1,0.5);
	     //leg->AddEntry(hPtGen,sample->label(),"P");
	   
	     // Linear scale
	     out_->nextMultiPad(sample->label()+": Extrapolation "+ptBin->toTString());
	     hFrame->Draw();
	     gVals->Draw("PE1same");
	     fit->Draw("same");
	     //label->Draw("same");
	     //leg->Draw("same");
	     out_->saveCurrentPad(histFileName("Extrapolation",ptBin,sample,*rIt));
	 
	     delete gVals;
	     delete fit;
	     delete hFrame;
	     //delete label;
	     //delete leg;
	   } // End of loop over FitResultTypes
	 } // End of loop over Samples
       } // End of loop over pt bins
     } // End of loop over eta bins
     
     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotExtrapolation(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  void PlotMaker::plotPtGenSpectra() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtGenSpectra(): Entering" << std::endl;
    std::cout << "Plotting ptGen spectra" << std::endl;

    // Loop over eta and pt bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      out_->newPage("PtGenSpectrum");

      for(PtBinIt ptBinIt = (*etaBinIt)->ptBinsBegin(); ptBinIt != (*etaBinIt)->ptBinsEnd(); ++ptBinIt) {
 	const PtBin* ptBin = *ptBinIt;
 	if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
	
 	// Loop over MCSamples in this bin
 	for(MCSampleIt sIt = ptBin->mcSamplesBegin(); sIt != ptBin->mcSamplesEnd(); ++sIt) {
 	  const MCSample* sample = sIt->second;
 	  if( par_->verbosity() > 1 ) std::cout << "    " << sample->label() << std::endl;
	  
 	  for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < sample->nPtSoftBins(); ++ptSoftBinIdx) {
 	    if( par_->verbosity() > 1 ) std::cout << "      PtSoftBin " << ptSoftBinIdx << std::endl;
	    
 	    // Get ptGen spectrum from sample and tweak style
 	    TH1* hPtGen = sample->histPtGen(ptSoftBinIdx);
 	    util::HistOps::setAxisTitles(hPtGen,"p^{gen}_{T}","GeV","jets",true);
	    
 	    // Label
 	    //TPaveText* label = labelMk_->binLabel(bin,ptSoftBinIdx);
 	    //TLegend* leg = util::LabelFactory::createLegendCol(1,0.5);
 	    //leg->AddEntry(hPtGen,sample->label(),"P");
	    
 	    // Linear scale
 	    util::HistOps::setYRange(hPtGen,3);
 	    out_->nextMultiPad(sample->label()+": PtGen Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
 	    hPtGen->Draw("PE1");
 	    //label->Draw("same");
 	    //leg->Draw("same");
 	    out_->saveCurrentPad(histFileName("PtGen",ptBin,sample,ptSoftBinIdx));
	    
 	    delete hPtGen;
 	    //delete label;
 	    //delete leg;
 	  } // End of loop over ptSoft bins	  
 	} // End of loop over MCSamples
       } // End of loop over pt bins
     } // End of loop over eta bins
    
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtGenSpectra(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  void PlotMaker::plotResolution() const {
     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotResolution(): Entering" << std::endl;
     std::cout << "Plotting resolution" << std::endl;
     out_->newPage("Resolution");

     // Loop over eta bins
     for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
       const EtaBin* etaBin = *etaBinIt;
       const unsigned int etaBinIdx = etaBin->etaBin();
       if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;

       // Loop over FitResultTypes
       for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); rIt != etaBin->fitResultTypesEnd(); ++rIt) {
	 FitResult::Type fitResType = *rIt;
	 if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(fitResType) << std::endl;

 	 // Loop over SampleLabels
 	 for(SampleTypeIt sTIt = etaBin->sampleTypesBegin(); sTIt != etaBin->sampleTypesEnd(); ++sTIt) {
	   SampleLabel sampleLabel = sTIt->first;
	   Sample::Type sampleType = sTIt->second;
 	   if( par_->verbosity() > 1 ) std::cout << "      " << sampleLabel << std::endl;

	   // Create graphs of extrapolated resolution
	   // and corrected for PLI
	   std::vector<double> pt;
	   std::vector<double> ptErr;
	   std::vector<double> res;
	   std::vector<double> resStatErr;
	   std::vector<double> resPLICorr;
	   std::vector<double> resPLICorrStatErr;
	   for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	     unsigned int ptBinIdx = (*ptBinIt)->ptBin();
	     pt.push_back(etaBin->meanPt(sampleLabel,fitResType,ptBinIdx));
	     ptErr.push_back(0.);
	     res.push_back(etaBin->extrapolatedValue(sampleLabel,fitResType,ptBinIdx));
	     resStatErr.push_back(etaBin->extrapolatedStatUncert(sampleLabel,fitResType,ptBinIdx));
	     resPLICorr.push_back(etaBin->correctedResolution(sampleLabel,fitResType,ptBinIdx));
	     resPLICorrStatErr.push_back(etaBin->correctedResolutionStatUncert(sampleLabel,fitResType,ptBinIdx));
	   }
	 
	   TGraphAsymmErrors* gRes = 
	     new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),&(ptErr.front()),&(ptErr.front()),
				   &(resStatErr.front()),&(resStatErr.front()));
	   gRes->SetMarkerStyle(markerStyleExtrapolatedResolution(sampleType));

	   TGraphAsymmErrors* gResCorr = 
	     new TGraphAsymmErrors(pt.size(),&(pt.front()),&(resPLICorr.front()),
				   &(ptErr.front()),&(ptErr.front()),
				   &(resPLICorrStatErr.front()),&(resPLICorrStatErr.front()));
	   gResCorr->SetMarkerStyle(markerStyleCorrectedResolution(sampleType));
	     
	   // MC truth resolution function
	   TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBinIdx));
	   mcTruth->SetLineColor(kRed);

	   // PLI function
	   TF1* pli = etaBin->pliFunc("PLI_Eta"+util::toTString(etaBinIdx));
	   pli->SetLineColor(kGreen);
	   pli->SetLineStyle(2);

	   // Create frame
	   TH1* hFrame = new TH1D("PlotMaker::plotResolution::hFrame","",1000,xMinPt_,xMaxPt_);
	   util::HistOps::setAxisTitles(hFrame,"p^{ref}_{T}","GeV","#sigma / p^{ref}_{T}");
	   hFrame->GetYaxis()->SetRangeUser(yMinExtraRes_,yMaxExtraRes_);
	   hFrame->GetXaxis()->SetMoreLogLabels();
	   hFrame->GetXaxis()->SetNoExponent();
	     
	   // Label
	   //TPaveText* label = labelMk_->binLabel(bin);
	   //TLegend* leg = util::LabelFactory::createLegendCol(1,0.5);
	   //leg->AddEntry(hPtGen,sample->label(),"P");
	   
	   out_->nextPad(sampleLabel+": Resolution "+etaBin->toString());
	   hFrame->Draw();
	   pli->Draw("same");
	   mcTruth->Draw("same");
	   gRes->Draw("PE1same");
	   gResCorr->Draw("PE1same");

	   //label->Draw("same");
	   //leg->Draw("same");
	   out_->logx();
	   out_->saveCurrentPad(histFileName("Resolution",etaBin,sampleLabel,fitResType));
	   
	   delete gRes;
	   delete mcTruth;
	   delete hFrame;
	   //delete label;
	   //delete leg;
	 } // End of loop over SampleLabels
       } // End of loop over FitResultTypes
     } // End of loop over eta bins
     
     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotResolution(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const EtaBin* etaBin, SampleLabel sampleLabel, FitResult::Type type) const {
    return par_->outFilePrefix()+"_"+id+"_"+(sampleLabel.ReplaceAll(" ","_"))+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(etaBin->etaBin())+".eps";
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, FitResult::Type type) const {
    return par_->outFilePrefix()+"_"+id+"_"+((sample->label()).ReplaceAll(" ","_"))+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+".eps";
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, unsigned int ptSoftBinIdx) const {
    return par_->outFilePrefix()+"_"+id+"_"+((sample->label()).ReplaceAll(" ","_"))+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx)+".eps";
  }



  // -------------------------------------------------------------------------------------
  int PlotMaker::markerStyleExtrapolatedResolution(Sample::Type type) const {
    return 27;
  }



  // -------------------------------------------------------------------------------------
  int PlotMaker::markerStyleCorrectedResolution(Sample::Type type) const {
    int style = 1;
    if( Sample::validType(type) ) {
      if( type == Sample::Data ) {
	style = 20;
      } else if( type == Sample::MC ) {
	style = 25;
      }
    }
    
    return style;
  }



  // -------------------------------------------------------------------------------------
  PlotMaker::LabelMaker::LabelMaker(const Parameters* par)
    : par_(par) {}



//   // -------------------------------------------------------------------------------------
//   TPaveText* PlotMaker::LabelMaker::binLabel(unsigned int etaBin, unsigned int ptBin) const {
//     TPaveText* label = util::LabelFactory::createPaveText(2,-0.6);
//     label->AddText(jetAlgo()+", "+etaRange(etaBin));
//     label->AddText(ptRange(etaBin,ptBin));

//     return label;
//   }



//   // -------------------------------------------------------------------------------------
//   TPaveText* PlotMaker::LabelMaker::binLabel(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBinIdx) const {
//     TPaveText* label = util::LabelFactory::createPaveText(3,-0.6);
//     label->AddText(jetAlgo()+", "+etaRange(etaBin));
//     label->AddText(ptSoftRange(ptSoftBinIdx));
//     label->AddText(ptRange(etaBin,ptBin));

//     return label;
//   }



//   // -------------------------------------------------------------------------------------
//   TString PlotMaker::LabelMaker::etaRange(unsigned int etaBin) const {
//     return util::toTString(par_->etaMin(etaBin))+" < |#eta| < "+util::toTString(par_->etaMax(etaBin));
//   }



//   // -------------------------------------------------------------------------------------
//   TString PlotMaker::LabelMaker::ptRange(unsigned int etaBin, unsigned int ptBin) const {
//     return util::toTString(par_->ptMin(etaBin,ptBin))+" < "+pt()+" < "+util::toTString(par_->ptMax(etaBin,ptBin))+" GeV";
//   }



//   // -------------------------------------------------------------------------------------
//   TString PlotMaker::LabelMaker::ptSoftRange(unsigned int ptSoftBinIdx) const {
//     TString label;
//     if( par_->extrapolationType() == Extrapolation::Threshold ) {
//       label = ptSoft()+" < "+util::toTString(par_->ptSoftMax(ptSoftBinIdx))+"#upoint#bar{p^{ave}_{T}}";
//     }
//     return label;
//   }



//   // -------------------------------------------------------------------------------------
//   TString PlotMaker::LabelMaker::jetAlgo() const {
//     TString algo = "default jets";
//     if( par_->jetAlgo() == AK5PF ) {
//       algo = "AK5 PF-Jets";
//     }
//     return algo;
//   }



//   // -------------------------------------------------------------------------------------
//   TString PlotMaker::LabelMaker::pt() const { 
//     return "p^{ave}_{T}";
//   }



//   // -------------------------------------------------------------------------------------
//   TString PlotMaker::LabelMaker::ptSoft() const { 
//     return "p_{||,3}";
//   }
}
