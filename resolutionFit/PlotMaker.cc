// $Id: PlotMaker.cc,v 1.9 2011/03/01 16:52:41 mschrode Exp $

#include "PlotMaker.h"

#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TPad.h"

#include "FitResult.h"
#include "SystematicUncertainty.h"

#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  PlotMaker::PlotMaker(const Parameters *par, const EtaBins &etaBins)
    : par_(par), etaBins_(etaBins),
      xMinPt_(40.), xMaxPt_(1100.),
      yMinExtraRes_(1E-3), yMaxExtraRes_(0.38), yMinResRatio_(0.74), yMaxResRatio_(1.26) {

    // Output manager
    out_ = OutputManager::createOutputManager(par_->outMode(),par_->outFilePrefix());

    // Style
    labelMk_ = new LabelMaker(par_);

    title_ = "";
  }


  // -------------------------------------------------------------------------------------
  PlotMaker::~PlotMaker() {
    delete out_;
    delete labelMk_;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotAsymmetry() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetry(): Entering" << std::endl;
    std::cout << "Plotting asymmetry" << std::endl;

    const double asymMax = 0.49;

    // +++++ Asymmetry plots per bin, with different FitResult types ++++++++++++++++++++++++++++++

    // Loop over SampleLabels
    for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      
      SampleLabel sampleLabel = sTIt->first;
      if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
      
      // Loop over ptSoft bins
      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	if( par_->verbosity() > 1 ) std::cout << "      PtSoftBin " << ptSoftBinIdx << std::endl;
	
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  
	  out_->newPage("Asymmetry");
	  
	  // Loop over pt bins
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "    " << ptBin->toTString() << std::endl;
	    
	    // Loop over Samples and select current one (by sampleLabel)
	    for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
	      const Sample* sample = sIt->second;
	      if( sample->label() != sampleLabel ) continue;
	      
	      // Asymmetry distribution
	      TH1* hPtAsym = sample->histPtAsym(ptSoftBinIdx);
	      util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","events");
	      hPtAsym->GetXaxis()->SetRangeUser(-asymMax,asymMax);
	      setStyle(sample->type(),hPtAsym);
	      
	      // Asymmetry fits
	      std::vector<TF1*> fits;
	      std::vector<TString> labels;

	      // Loop over FitResultTypes
	      for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin();
		  rIt != etaBin->fitResultTypesEnd(); ++rIt) {

		double sigma = sample->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
		fits.push_back(new TF1("AsymmetryFit_"+FitResult::toString(*rIt),"gaus",-1.,1.));
		fits.back()->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
		fits.back()->SetParameter(1,0.);
		fits.back()->SetParameter(2,sigma);

		labels.push_back("Fit: "+FitResult::toString(*rIt));
	      } // End of loop over FitResultTypes

	      for(unsigned int i = 0; i < fits.size(); ++i) {
		fits.at(i)->SetLineWidth(1);
		fits.at(i)->SetLineStyle(1+i);
		fits.at(i)->SetLineColor(1+i);
	      }
	      
	      // Labels
	      TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(1+labels.size(),labelMk_->start(),label->GetSize());
	      leg->AddEntry(hPtAsym,sample->label(),"P");
	      for(unsigned int i = 0; i < fits.size(); ++i) {
		leg->AddEntry(fits.at(i),labels.at(i),"L");
	      }
	    
	      util::HistOps::setYRange(hPtAsym,label->GetSize()+leg->GetNRows()+1);
	      out_->nextMultiPad(sample->label()+": PtAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hPtAsym->Draw("PE1");
	      for(unsigned int i = 0; i < fits.size(); ++i) {
		fits.at(i)->Draw("same");
	      }
	      label->Draw("same");
	      leg->Draw("same");
	      out_->saveCurrentPad(histFileName("PtAsym",ptBin,sample,ptSoftBinIdx));

	      delete hPtAsym;
	      for(unsigned int i = 0; i < fits.size(); ++i) {
		delete fits.at(i);
	      }
	      delete label;
	      delete leg;
	    } // End of loop over Samples
	  } // End of loop over pt bins
	} // End of loop over eta bins
      } // End of loop over ptSoft bins
    } // End of loop over SampleLabels



    // +++++ Asymmetry plots per bin, different Samples ++++++++++++++++++++++++++++++
    
    if( etaBins_.front()->nComparedSamples() > 0 ) {
      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); 
	  rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
	if( par_->verbosity() > 1 ) std::cout << "  " << FitResult::toString(*rIt) << std::endl;

	// Loop over ptSoft bins
	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	  
	  // Loop over eta and bins
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;
	    out_->newPage("Asymmetry");
	    
	    // Loop over pt bins
	    for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	      const PtBin* ptBin = *ptBinIt;
	      if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;

	      // Loop over to-be-compare Samples
	      for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
		  sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
 		SampleLabel sLabel1 = (*sCIt)->label1();
 		SampleLabel sLabel2 = (*sCIt)->label2();
		const Sample* sample1 = ptBin->findSample(sLabel1);
		const Sample* sample2 = ptBin->findSample(sLabel2);

		// Asymmetry distributions and fits
		TH1* hPtAsym1 = sample1->histPtAsym(ptSoftBinIdx);
		util::HistOps::setAxisTitles(hPtAsym1,"Asymmetry","","events");
		hPtAsym1->GetXaxis()->SetRangeUser(-asymMax,asymMax);
		setStyle(sample1->type(),hPtAsym1);

		double sigma = sample1->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
		TF1* fit1 = new TF1("AsymmetryFit_Sample1","gaus",-1.,1.);
		fit1->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
		fit1->SetParameter(1,0.);
		fit1->SetParameter(2,sigma);
		fit1->SetLineWidth(1);
		fit1->SetLineStyle(1);
		fit1->SetLineColor(hPtAsym1->GetLineColor());

		TH1* hPtAsym2 = sample2->histPtAsym(ptSoftBinIdx);
		util::HistOps::setAxisTitles(hPtAsym2,"Asymmetry","","events");
		hPtAsym2->GetXaxis()->SetRangeUser(-asymMax,asymMax);
		setStyle(sample2->type(),hPtAsym2);

		sigma = sample2->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
		TF1* fit2 = new TF1("AsymmetryFit_Sample1","gaus",-1.,1.);
		fit2->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
		fit2->SetParameter(1,0.);
		fit2->SetParameter(2,sigma);
		fit2->SetLineWidth(1);
		fit2->SetLineStyle(2);
		fit2->SetLineColor(hPtAsym2->GetLineColor());
			      
		// Labels
		TPaveText* label = labelMk_->ptSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
		TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
		leg->AddEntry(hPtAsym1,sLabel1,"P");
		leg->AddEntry(hPtAsym2,sLabel2,"P");
	    
		util::HistOps::setYRange(hPtAsym1,label->GetSize()+leg->GetNRows()+1);
		out_->nextMultiPad(sLabel1+" vs "+sLabel2+": PtAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		hPtAsym1->Draw("PE1");
		hPtAsym2->Draw("PE1same");
		fit1->Draw("same");
		fit2->Draw("same");
		label->Draw("same");
		leg->Draw("same");
		out_->saveCurrentPad(histFileName("PtAsym",ptBin,sLabel1,sLabel2,*rIt,ptSoftBinIdx));

		delete hPtAsym1;
		delete hPtAsym2;
		delete label;
		delete leg;
	      } // End of loop over to-be-compared samples
	    } // End of loop over pt bins
	  } // End of loop over eta bins
	} // End of loop over ptSoft bins
      } // End of loop over FitResult types
    } // End if samples to be compared

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetry(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotExtrapolation() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotExtrapolation(): Entering" << std::endl;
    std::cout << "Plotting extrapolation" << std::endl;
     

    // +++++ Extrapolation plots per Sample ++++++++++++++++++++++++++++++++++++++

    // Loop over FitResultTypes
    for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); 
	rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {

      // Loop over SampleLabels
      for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	  sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
	SampleLabel sampleLabel = sTIt->first;
	if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
	 
	// Loop over eta and pt bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  out_->newPage("Extrapolation");
	 
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "    " << ptBin->toTString() << std::endl;
	     
	    // Loop over Samples and select current one (by sampleLabel)
	    for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
	      const Sample* sample = sIt->second;
	      if( sample->label() != sampleLabel ) continue;
	       
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
	      setStyle(sample->type(),gVals);
	       
	      // Extrapolation function
	      TF1* fit = sample->extrapolationFunction(*rIt,"ExtrapolationForPlotFitRange");
	      fit->SetLineColor(kRed);
	      TF1* extra = static_cast<TF1*>(fit->Clone("ExtrapolationForPlotPlotRange"));
	      extra->SetRange(0.,1.4*ptSoftx.back());
	      extra->SetLineStyle(2);
	       
	      // Create frame
	      TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame","",
				     1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	      util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"","#sigma / p^{ref}_{T}");
	      hFrame->GetYaxis()->SetRangeUser(0.8*val.front(),1.3*val.back());
	       
	      // Labels
	      TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	      leg->AddEntry(gVals,labelMk_->label(*rIt),"P");
	      leg->AddEntry(fit,"Extrapolation","L");
	       
	      // Linear scale
	      out_->nextMultiPad(sample->label()+": Extrapolation "+ptBin->toTString());
	      hFrame->Draw();
	      gVals->Draw("PE1same");
	      extra->Draw("same");
	      fit->Draw("same");
	      label->Draw("same");
	      leg->Draw("same");
	      out_->saveCurrentPad(histFileName("Extrapolation",ptBin,sample,*rIt));
	       
	      delete gVals;
	      delete extra;
	      delete fit;
	      delete hFrame;
	      delete label;
	      delete leg;
	    } // End of loop over Samples
	  } // End of loop over pt bins
	} // End of loop over eta bins
      } // End of loop over SampleLabels
    } // End of loop over FitResultTypes
     



    // +++++ Extrapolation plots, comparison of different Sample +++++++++++++++++++++++++++++++

    if( etaBins_.front()->nComparedSamples() > 0 ) {
      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); 
	  rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
	if( par_->verbosity() > 1 ) std::cout << "  " << FitResult::toString(*rIt) << std::endl;
	 
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  out_->newPage("Extrapolation");
	  
	  // Loop over pt bins
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "    " << ptBin->toTString() << std::endl;

	    // Graphs and functions
	    TGraphAsymmErrors* gVals1 = 0;
	    TGraphAsymmErrors* gVals2 = 0;
	    TF1* fit1 = 0;
	    TF1* fit2 = 0;
	    TF1* extra1 = 0;
	    TF1* extra2 = 0;

	    // Labels
	    TPaveText* label = labelMk_->ptBin(etaBin->etaBin(),ptBin->ptBin());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());

	    // Loop over to-be-compare Samples
	    for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
		sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	      SampleLabel sLabel1 = (*sCIt)->label1();
	      SampleLabel sLabel2 = (*sCIt)->label2();
	      const Sample* sample1 = ptBin->findSample(sLabel1);
	      const Sample* sample2 = ptBin->findSample(sLabel2);
	     
	      // Graphs of fitted resolutions
	      std::vector<double> val;
	      std::vector<double> uncert;
	      sample1->values(*rIt,val,uncert);
	      std::vector<double> ptSoftx;
	      sample1->ptSoft(ptSoftx);
	      std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
	      gVals1 = new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
					     &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
					     &(uncert.front()),&(uncert.front()));
	      setStyle(sample1->type(),gVals1);

	      val.clear();
	      uncert.clear();
	      sample2->values(*rIt,val,uncert);
	      ptSoftx.clear();
	      sample2->ptSoft(ptSoftx);
	      ptSoftxUncert.clear();
	      ptSoftxUncert = std::vector<double>(ptSoftx.size(),0.);
	      gVals2 = new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
					     &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
					     &(uncert.front()),&(uncert.front()));
	      setStyle(sample2->type(),gVals2);

	      // Extrapolation functions
	      fit1 = sample1->extrapolationFunction(*rIt,"Extrapolation_"+((sLabel1).ReplaceAll(" ","_"))+"_FitRange");
	      fit1->SetLineColor(gVals1->GetMarkerColor());
	      extra1 = static_cast<TF1*>(fit1->Clone("Extrapolation_"+((sLabel1).ReplaceAll(" ","_"))+"_PlotRange"));
	      extra1->SetRange(0.,1.4*ptSoftx.back());
	      extra1->SetLineStyle(2);

	      fit2 = sample2->extrapolationFunction(*rIt,"Extrapolation_"+((sLabel2).ReplaceAll(" ","_"))+"_FitRange");
	      fit2->SetLineColor(gVals2->GetMarkerColor());
	      extra2 = static_cast<TF1*>(fit2->Clone("Extrapolation_"+((sLabel2).ReplaceAll(" ","_"))+"_PlotRange"));
	      extra2->SetRange(0.,1.4*ptSoftx.back());
	      extra2->SetLineStyle(2);

	      leg->AddEntry(gVals1,sLabel1,"P");
	      leg->AddEntry(gVals2,sLabel2,"P");
	    } // End of loop over Samples
	   
	      // Create frame
	    TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame","",
				   1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	    util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"","#sigma / p^{ref}_{T}");
	    hFrame->GetYaxis()->SetRangeUser(0.8*std::min(gVals1->GetY()[0],gVals2->GetY()[1]),
					     1.3*std::max(gVals1->GetY()[gVals1->GetN()-1],gVals2->GetY()[gVals2->GetN()-1]));
	     
	    // Linear scale
	    out_->nextMultiPad("Sample comparison: Extrapolation "+ptBin->toTString());
	    hFrame->Draw();
	    gVals1->Draw("PE1same");
	    extra1->Draw("same");
	    fit1->Draw("same");
	    gVals2->Draw("PE1same");
	    extra2->Draw("same");
	    fit2->Draw("same");
	    label->Draw("same");
	    leg->Draw("same");
	    out_->saveCurrentPad(histFileName("Extrapolation",ptBin,*rIt));
	 
	    delete gVals1;
	    delete extra1;
	    delete fit1;
	    delete gVals2;
	    delete extra2;
	    delete fit2;
	    delete hFrame;
	    delete label;
	    delete leg;
	  } // End of loop over FitResultTypes
	} // End of loop over pt bins
      } // End of loop over eta bins
    } // End if nSampleTypes() > 1
     



      // +++++ Extrapolation plots, comparison of different FitResult types +++++++++++++++++++++++++++++++

      // Loop over SampleLabels
    if( etaBins_.front()->nFitResultTypes() > 1 ) {
      for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	  sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
	SampleLabel sampleLabel = sTIt->first;
	if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
	 
	// Loop over eta and pt bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  out_->newPage("Extrapolation");
	 
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "    " << ptBin->toTString() << std::endl;
	     
	    // Loop over Samples and select current one (by sampleLabel)
	    for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
	      const Sample* sample = sIt->second;
	      if( sample->label() != sampleLabel ) continue;
   
	      std::vector<TGraphAsymmErrors*> gVals;
	      double valMin = 1000.;
	      double valMax = 0.;
	      std::vector<TF1*> fits;
	      std::vector<TF1*> extra;
	      std::vector<TString> labels;

	      // Loop over FitResultTypes
	      unsigned int frt = 0;
	      for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); 
		  rIt != etaBin->fitResultTypesEnd(); ++rIt,++frt) {

		labels.push_back(labelMk_->label(*rIt));
	     
		// Graph of fitted resolutions
		std::vector<double> val;
		std::vector<double> uncert;
		sample->values(*rIt,val,uncert);
		std::vector<double> ptSoftx;
		sample->ptSoft(ptSoftx);
		std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
	     
		gVals.push_back(new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
						      &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
						      &(uncert.front()),&(uncert.front())));
		gVals.back()->SetMarkerStyle(20+frt);
		gVals.back()->SetLineColor(1+frt);
		gVals.back()->SetMarkerColor(gVals.back()->GetLineColor());
	     
		// Extrapolation function
		fits.push_back(sample->extrapolationFunction(*rIt,"Extrapolation_"+((sample->label()).ReplaceAll(" ","_"))+"_FitRange"));
		fits.back()->SetLineColor(gVals.back()->GetLineColor());
		extra.push_back(static_cast<TF1*>(fits.back()->Clone("Extrapolation_"+((sample->label()).ReplaceAll(" ","_"))+"_PlotRange")));
		extra.back()->SetRange(0.,1.4*ptSoftx.back());
		extra.back()->SetLineStyle(2);
	       
		if( val.front() < valMin ) valMin = val.front();
		if( val.back() > valMax ) valMax = val.back();
	      } // End of loop over FitResultTypes
	     
	   
		// Create frame
	      TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame","",
				     1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	      util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"","#sigma / p^{ref}_{T}");
	      hFrame->GetYaxis()->SetRangeUser(0.8*valMin,1.3*valMax);
	     
	      // Labels
	      TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(labels.size(),labelMk_->start(),label->GetSize());
	      for(size_t i = 0; i < labels.size(); ++i) {
		leg->AddEntry(gVals.at(i),labels.at(i),"P");
	      }
	   
	      // Linear scale
	      out_->nextMultiPad("FitResult comparison: Extrapolation "+ptBin->toTString());
	      hFrame->Draw();
	      for(size_t i = 0; i < gVals.size(); ++i) {
		gVals.at(i)->Draw("PE1same");
		extra.at(i)->Draw("same");
		fits.at(i)->Draw("same");
	      }
	      label->Draw("same");
	      leg->Draw("same");
	      out_->saveCurrentPad(histFileName("Extrapolation",ptBin,sample));
	     
	      for(size_t i = 0; i < gVals.size(); ++i) {
		delete gVals.at(i);
		delete extra.at(i);
		delete fits.at(i);
	      }
	      delete hFrame;
	      delete label;
	      delete leg;
	    } // End of loop over Samples
	  } // End of loop over pt bins
	} // End of loop over eta bins
      } // End loop over SampleTypes
    } // End if nFitResultTypes() > 1
     
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotExtrapolation(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotParticleLevelImbalance() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetry(): Entering" << std::endl;
    if( etaBins_.front()->hasMCTruthSample() ) {

      std::cout << "Plotting particle level imabalance control plots" << std::endl;

      const double asymMax = 0.29;

      // +++++ PtGenAsymmetry per bin ++++++++++++++++++++++++++++++++++++++++++++++++

      // Loop over ptSoft bins
      std::cout << "  PtGen asymmetry" << std::endl;
      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  
	  out_->newPage("PtGenAsymmetry");
	  
	  // Loop over pt bins
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	    
	    const Sample* sample = ptBin->mcTruthSample();
	      
	    // PtGenAsymmetry distribution
	    TH1* hPtGenAsym = sample->histPtGenAsym(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hPtGenAsym,"p^{gen}_{T} Asymmetry","","events");
	    hPtGenAsym->SetTitle(title_);
	    hPtGenAsym->GetXaxis()->SetRangeUser(-asymMax,asymMax);
	    hPtGenAsym->SetMarkerStyle(24);
	    hPtGenAsym->SetMarkerColor(kBlack);
	    hPtGenAsym->SetLineColor(hPtGenAsym->GetMarkerColor());
	    hPtGenAsym->SetLineWidth(1);
	      
	    // Labels
	    TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
	    // Draw
	    util::HistOps::setYRange(hPtGenAsym,label->GetSize()+1);
	    out_->nextMultiPad(sample->label()+": PLIPtGenAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGenAsym->Draw("PE1");
	    label->Draw("same");
	    out_->saveCurrentPad(histFileName("PLIPtGenAsym",ptBin,sample,ptSoftBinIdx));
	      
	    delete hPtGenAsym;
	    delete label;
	  } // End of loop over pt bins
	} // End of loop over eta bins
      } // End of loop over ptSoft bins


      // +++++ PtGen spectra per bin ++++++++++++++++++++++++++++++++++++++++++++++++

      // Loop over ptSoft bins
      std::cout << "  PtGen spectra" << std::endl;
      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  
	  out_->newPage("PtGenSpectra");
	  
	  // Loop over pt bins
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	    
	    const Sample* sample = ptBin->mcTruthSample();
	      
	    // PtGen distribution
	    TH1* hPtGen = sample->histPtGen(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hPtGen,"p^{gen}_{T}","GeV","events");
	    hPtGen->SetTitle(title_);
	    hPtGen->SetMarkerStyle(24);
	    hPtGen->SetMarkerColor(kBlack);
	    hPtGen->SetLineColor(hPtGen->GetMarkerColor());
	    hPtGen->SetLineWidth(1);
	      
	    // Labels
	    TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
	    // Draw
	    util::HistOps::setYRange(hPtGen,label->GetSize()+1);
	    out_->nextMultiPad(sample->label()+": PLIPtGenSpectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGen->Draw("PE1");
	    label->Draw("same");
	    out_->saveCurrentPad(histFileName("PLIPtGenSpectrum",ptBin,sample,ptSoftBinIdx));
	      
	    delete hPtGen;
	    delete label;
	  } // End of loop over pt bins
	} // End of loop over eta bins
      } // End of loop over ptSoft bins



      // +++++ PtGenAsym Extrapolation +++++++++++++++++++++++++++++++++++++++++++++
      
      // Loop over eta and pt bins
      std::cout << "  Extrapolation" << std::endl;
      for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	const EtaBin* etaBin = *etaBinIt;
	out_->newPage("PtGenExtrapolation");
      
	for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	const PtBin* ptBin = *ptBinIt;
	if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;

	const Sample* sample = ptBin->mcTruthSample();
	
	// Graph of fitted resolutions
	std::vector<double> val;
	std::vector<double> uncert;
	sample->values(FitResult::PtGenAsym,val,uncert);
	std::vector<double> ptSoftx;
	sample->ptSoft(ptSoftx);
	std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
	       
	TGraphAsymmErrors* gVals = 
	  new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
				&(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
				&(uncert.front()),&(uncert.front()));
	gVals->SetMarkerStyle(24);
	       
	// Extrapolation function
	TF1* fit = sample->extrapolationFunction(FitResult::PtGenAsym,"ExtrapolationForPlotFitRange");
	fit->SetLineColor(kRed);
	TF1* extra = static_cast<TF1*>(fit->Clone("ExtrapolationForPlotPlotRange"));
	extra->SetRange(0.,1.4*ptSoftx.back());
	extra->SetLineStyle(2);
	       
	// Create frame
	TH1* hFrame = new TH1D("PlotMaker::plotParticleLevelImbalance::hFrame","",
			       1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	util::HistOps::setAxisTitles(hFrame,"p^{gen}_{3,T} / #LTp^{gen}_{T}#GT","","#sqrt{2}#upoint#sigma(p^{gen} Asymmetry) / #LTp^{gen}_{T}#GT");
	hFrame->GetYaxis()->SetRangeUser(0.,0.19);
	       
	// Labels
	TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	leg->AddEntry(gVals,"#sqrt{2}#upoint#sigma(p^{gen} Asymmetry)","P");
	leg->AddEntry(fit,"Extrapolation","L");
	       
	// Linear scale
	out_->nextMultiPad(sample->label()+": PLI PtGen Asymmetry Extrapolation "+ptBin->toTString());
	hFrame->Draw();
	gVals->Draw("PE1same");
	extra->Draw("same");
	fit->Draw("same");
	label->Draw("same");
	leg->Draw("same");
	out_->saveCurrentPad(histFileName("PLIPtGenAsymExtrapolation",ptBin,sample,FitResult::PtGenAsym));
	
	delete gVals;
	delete extra;
	delete fit;
	delete hFrame;
	delete label;
	delete leg;
	} // End of loop over pt bins
      } // End of loop over eta bins


      // +++++ Fitted PLI vs ptGenAve +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      // Loop over eta and pt bins
      std::cout << "  PLI fit" << std::endl;
      for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
 	const EtaBin* etaBin = *etaBinIt;
 	out_->newPage("PLIFit");
	  
	// Create graph of extrapolated ptGen asymmetry width
	std::vector<double> pt;
	std::vector<double> ptErr;
	std::vector<double> res;
	std::vector<double> resStatErr;
	for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	  const Sample* sample = (*ptBinIt)->mcTruthSample();
	  pt.push_back(sample->meanPt(FitResult::PtGenAsym));
	  ptErr.push_back(0.);
	  res.push_back(sample->extrapolatedValue(FitResult::PtGenAsym));
	  resStatErr.push_back(sample->extrapolatedStatUncert(FitResult::PtGenAsym));
	}
	TGraphAsymmErrors* gPLI = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),
							&(ptErr.front()),&(ptErr.front()),
							&(resStatErr.front()),&(resStatErr.front()));
	gPLI->SetMarkerStyle(24);
	
	// Fitted PLI function
	TF1* pli = etaBin->pliFunc("PLI_Eta"+util::toTString(etaBin->etaBin()));
	pli->SetRange(xMinPt_,xMaxPt_);
	pli->SetLineColor(kGreen);
	pli->SetLineStyle(2);
	
	// Ratio of measurement and fit
	TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gPLI,pli);
	
	// Create frames for main and ratio plots
	TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,0.109,"#sigma / #LTp^{gen}_{T}#GT");
	hFrameMain->GetXaxis()->SetMoreLogLabels();
	hFrameMain->GetXaxis()->SetNoExponent();
	
	TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"#LTp^{gen}_{T}#GT","GeV",yMinResRatio_,yMaxResRatio_);
	hFrameRatio->GetXaxis()->SetMoreLogLabels();
	hFrameRatio->GetXaxis()->SetNoExponent();
	
	// Labels
	TPaveText* label = labelMk_->etaBin(etaBin->mcTruthSampleLabel(),etaBin->etaBin());
	TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	leg->AddEntry(gPLI,"Particle Level Imbalance","P");
	leg->AddEntry(pli,"Fit","L");
	
	out_->nextMainPad(etaBin->mcTruthSampleLabel()+": PLI fit "+etaBin->toString());
	hFrameMain->Draw();
	gPLI->Draw("PE1same");
	pli->Draw("same");
	label->Draw("same");
	leg->Draw("same");
	out_->nextRatioPad();
	hFrameRatio->Draw();
	gRatio->Draw("PE1same");
	out_->logx();
	out_->saveCurrentPad(histFileName("PLIFit",etaBin,etaBin->mcTruthSampleLabel(),FitResult::PtGenAsym));
	
	delete gPLI;
	delete gRatio;
	delete pli;
	delete hFrameMain;
	delete hFrameRatio;
	delete label;
	delete leg;
      } // End of loop over eta bins
      
    } // End if hasMCTruthSample()

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetry(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotPtGenSpectra() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtGenSpectra(): Entering" << std::endl;
    std::cout << "Plotting ptGen spectra" << std::endl;

    // Loop over eta and pt bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
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
	    
	    // Get ptGen spectrum and fit from sample and tweak style
	    TH1* hPtGen = sample->histPtGen(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hPtGen,"p^{gen}_{T}","GeV","jets",true);
	    setStyle(sample->type(),hPtGen);

	    TH1* hPdf = sample->histPdfPtTrue(ptSoftBinIdx);
	    hPdf->SetLineColor(kRed);
	    
	    // Labels
	    TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	    leg->AddEntry(hPtGen,"Generated Spectrum","P");
	    leg->AddEntry(hPdf,"Assumed Spectrum #tilde{f}","L");
	    
	    // Linear scale
	    util::HistOps::setYRange(hPtGen,label->GetSize()+leg->GetNRows()+1);
	    out_->nextMultiPad(sample->label()+": PtGen Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGen->Draw("PE1");
	    hPdf->Draw("Lsame");
	    label->Draw("same");
	    leg->Draw("same");
	    out_->saveCurrentPad(histFileName("PtGen",ptBin,sample,ptSoftBinIdx));
	    
	    delete hPtGen;
	    delete label;
	    delete leg;
	  } // End of loop over ptSoft bins	  
	} // End of loop over MCSamples
      } // End of loop over pt bins
    } // End of loop over eta bins
    
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtGenSpectra(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotPtSpectra() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtSpectra(): Entering" << std::endl;
    std::cout << "Plotting pt spectra" << std::endl;

    // +++++ Pt spectra per bin ++++++++++++++++++++++++++++++++++++++++++++++

    // Loop over SampleLabels
    for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      
      SampleLabel sampleLabel = sTIt->first;
      if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
      
      // Loop over distributions
      for(unsigned int distIdx = 0; distIdx < 4; ++distIdx) {

	// Loop over ptSoft bins
	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	  
	  // Loop over eta bins
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;
	    
	    if( distIdx == 0 ) out_->newPage("PtAveSpectrum");
	    else if( distIdx == 1 ) out_->newPage("Pt1Spectrum");
	    else if( distIdx == 1 ) out_->newPage("Pt2Spectrum");
	    else if( distIdx == 1 ) out_->newPage("Pt3Spectrum");
	    
	    // Loop over pt bins
	    for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	      const PtBin* ptBin = *ptBinIt;
	      if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	    
	      // Loop over Samples and select current one (by sampleLabel)
	      for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
		const Sample* sample = sIt->second;
		if( sample->label() != sampleLabel ) continue;
	      
		// Spectrum
		TH1* h = 0;
		if( distIdx == 0 ) {
		  h = sample->histPtAve(ptSoftBinIdx);
		  util::HistOps::setAxisTitles(h,"p^{ave}_{T}","GeV","events");
		} else if( distIdx == 1 ) {
		  h = sample->histPt1(ptSoftBinIdx);
		  util::HistOps::setAxisTitles(h,"p_{T,1}","GeV","events");
		} else if( distIdx == 2 ) {
		  h = sample->histPt2(ptSoftBinIdx);
		  util::HistOps::setAxisTitles(h,"p_{T,2}","GeV","events");
		} else if( distIdx == 3 ) {
		  h = sample->histPt3(ptSoftBinIdx);
		  util::HistOps::setAxisTitles(h,"p_{T,3}","GeV","events");
		}
		//h->GetXaxis()->SetRangeUser(-asymMax,asymMax);
		setStyle(sample->type(),h);
	      
		// Labels
		TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
		util::HistOps::setYRange(h,label->GetSize()+1);
		out_->nextMultiPad(sample->label()+": Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		h->Draw("PE1");
		label->Draw("same");
		if( distIdx == 0 ) out_->saveCurrentPad(histFileName("PtAve",ptBin,sample,ptSoftBinIdx));
		else if( distIdx == 1 ) out_->saveCurrentPad(histFileName("PtJet1",ptBin,sample,ptSoftBinIdx));
		else if( distIdx == 2 ) out_->saveCurrentPad(histFileName("PtJet2",ptBin,sample,ptSoftBinIdx));
		else if( distIdx == 3 ) out_->saveCurrentPad(histFileName("PtJet3",ptBin,sample,ptSoftBinIdx));

		delete h;
		delete label;
	      } // End of loop over Samples
	    } // End of loop over pt bins
	  } // End of loop over eta bins
	} // End of loop over ptSoft bins
      } // End of loop over SampleLabels
    } // End of loop over spectra

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtSpectra(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotResolution() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotResolution(): Entering" << std::endl;
    std::cout << "Plotting resolution" << std::endl;
    out_->newPage("Resolution");

    // +++++ Resolution per Sample ++++++++++++++++++++++++++++++++++++++
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

	  TGraphAsymmErrors* gRes = etaBin->extrapolatedResolution(sampleLabel,fitResType);
	  setStyle(sampleType,gRes);
	  gRes->SetMarkerStyle(27);

	  TGraphAsymmErrors* gResCorr = etaBin->correctedResolution(sampleLabel,fitResType);
	  setStyle(sampleType,gResCorr);

	  // MC truth resolution function
	  TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBinIdx));
	  mcTruth->SetLineColor(kRed);

	  // Ratio of measurement and mc truth
	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gResCorr,mcTruth);

	  // PLI function
	  TF1* pli = etaBin->pliFunc("PLI_Eta"+util::toTString(etaBinIdx));
	  pli->SetRange(xMinPt_,xMaxPt_);
	  pli->SetLineColor(kGreen);
	  pli->SetLineStyle(2);

	  // Create frames for main and ratio plots
	  TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	  hFrameMain->SetTitle(title_);
	  hFrameMain->GetXaxis()->SetMoreLogLabels();
	  hFrameMain->GetXaxis()->SetNoExponent();

	  TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	  hFrameRatio->GetXaxis()->SetMoreLogLabels();
	  hFrameRatio->GetXaxis()->SetNoExponent();
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,labelMk_->start(),label->GetSize());
	  leg->AddEntry(gRes,"Resolution ("+FitResult::toString(fitResType)+")","P");
	  leg->AddEntry(pli,"Particle Level Imbalance (PLI)","L");
	  leg->AddEntry(gResCorr,"Resolution (PLI Corrected)","P");
	  leg->AddEntry(mcTruth,"Resolution (MC Truth)","L");

	  out_->nextMainPad(sampleLabel+": Resolution "+etaBin->toString());
	  hFrameMain->Draw();
	  pli->Draw("same");
	  mcTruth->Draw("same");
	  gRes->Draw("PE1same");
	  gResCorr->Draw("PE1same");
	  label->Draw("same");
	  leg->Draw("same");
	  out_->nextRatioPad();
	  hFrameRatio->Draw();
	  gRatio->Draw("PE1same");
	  out_->logx();
	  out_->saveCurrentPad(histFileName("Resolution",etaBin,sampleLabel,fitResType));

	  delete gRes;
	  delete gRatio;
	  delete mcTruth;
	  delete pli;
	  delete hFrameMain;
	  delete hFrameRatio;
	  delete label;
	  delete leg;
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins



    // +++++ Sample comparison ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over eta bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      const unsigned int etaBinIdx = etaBin->etaBin();
      if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;

      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); rIt != etaBin->fitResultTypesEnd(); ++rIt) {
	FitResult::Type fitResType = *rIt;
	if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(fitResType) << std::endl;

	// Loop over to-be-compare Samples
	for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	    sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	  SampleLabel sLabel1 = (*sCIt)->label1();
	  SampleLabel sLabel2 = (*sCIt)->label2();

	  // Create graphs of extrapolated resolution corrected for PLI
	  TGraphAsymmErrors* gRes1 = etaBin->correctedResolution(sLabel1,fitResType);
	  setStyle(etaBin->sampleType(sLabel1),gRes1);

	  TGraphAsymmErrors* gRes2 = etaBin->correctedResolution(sLabel2,fitResType);
	  setStyle(etaBin->sampleType(sLabel2),gRes2);

	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gRes1,gRes2);

	  // Create frames for main and ratio plots
	  TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	  hFrameMain->SetTitle(title_);
	  hFrameMain->GetXaxis()->SetMoreLogLabels();
	  hFrameMain->GetXaxis()->SetNoExponent();

	  TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	  hFrameRatio->GetXaxis()->SetMoreLogLabels();
	  hFrameRatio->GetXaxis()->SetNoExponent();
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	  leg->AddEntry(gRes1,sLabel1,"P");
	  leg->AddEntry(gRes2,sLabel2,"P");

	  out_->nextMainPad(sLabel1+"vs"+sLabel2+": Resolution "+etaBin->toString());
	  hFrameMain->Draw();
	  gRes1->Draw("PE1same");
	  gRes2->Draw("PE1same");
	  label->Draw("same");
	  leg->Draw("same");
	  out_->nextRatioPad();
	  hFrameRatio->Draw();
	  gRatio->Draw("PE1same");
	  out_->logx();
	  out_->saveCurrentPad(histFileName("Resolution",etaBin,sLabel1,sLabel2,fitResType));

	  delete gRes1;
	  delete gRes2;
	  delete gRatio;
	  delete hFrameMain;
	  delete hFrameRatio;
	  delete label;
	  delete leg;
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotResolution(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  void PlotMaker::plotScaledMCTruth() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotScaledMCTruth(): Entering" << std::endl;
    std::cout << "Plotting scaled MC truth" << std::endl;
    out_->newPage("MCTruthScaling");


    // +++++ Sample comparison (ratios) +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over eta bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      const unsigned int etaBinIdx = etaBin->etaBin();
      if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;

      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); rIt != etaBin->fitResultTypesEnd(); ++rIt) {
	FitResult::Type fitResType = *rIt;
	if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(fitResType) << std::endl;

	// Loop over to-be-compare Samples
	for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	    sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	  SampleLabel sLabel1 = (*sCIt)->label1();
	  SampleLabel sLabel2 = (*sCIt)->label2();

	  if( etaBin->hasKValue(sLabel1,sLabel2,fitResType) ) {

	    // Data-MC ratio (k-Value) with statistical und systematic uncertainty
	    TGraphAsymmErrors* gRatio = etaBin->ratioGraph(sLabel1,sLabel2,fitResType);
	    setStyle(etaBin->sampleType(sLabel1),gRatio);
	    TF1* kValueLine = etaBin->kValueLine(sLabel1,sLabel2,fitResType,"kValueLine",xMinPt_,xMaxPt_);
	    TGraphAsymmErrors* kStatBand = etaBin->kValueStatBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    TGraphAsymmErrors* kSystBand = etaBin->kValueSystBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    
	    // Create frame
	    TH1* hFrame = util::HistOps::createRatioFrame(xMinPt_,xMaxPt_,0.71,1.89,"p^{ref}_{T} (GeV)","#sigma("+sLabel1+") / #sigma("+sLabel2+")");
	    hFrame->SetTitle(title_);
	    hFrame->GetXaxis()->SetMoreLogLabels();
	    hFrame->GetXaxis()->SetNoExponent();
	     
	    // Labels
	    TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,-labelMk_->start(),label->GetSize());
	    leg->AddEntry(gRatio,"Measurement","P");
	    leg->AddEntry(kValueLine,"Fit","L");
	    leg->AddEntry(kStatBand,"Statistical Uncertainty","F");
	    leg->AddEntry(kSystBand,"Systematic Uncertainty","F");

	    out_->nextPad(sLabel1+"Over"+sLabel2+": Resolution "+etaBin->toString());
	    hFrame->Draw();
	    kSystBand->Draw("E2same");
	    kStatBand->Draw("E2same");
	    kValueLine->Draw("same");
	    gRatio->Draw("PE1same");
	    hFrame->Draw("same");
	    label->Draw("same");
	    leg->Draw("same");
	    out_->logx();
	    gPad->RedrawAxis();
	    out_->saveCurrentPad(histFileName("ResolutionRatio",etaBin,sLabel1,sLabel2,fitResType));

	    delete gRatio;
	    delete kValueLine;
	    delete kStatBand;
	    delete hFrame;
	    delete label;
	    delete leg;
	  } // End if hasKValue
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins



    // +++++ Sample comparison (scaled MC truth) +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over eta bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      const unsigned int etaBinIdx = etaBin->etaBin();
      if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;

      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); rIt != etaBin->fitResultTypesEnd(); ++rIt) {
	FitResult::Type fitResType = *rIt;
	if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(fitResType) << std::endl;

	// Loop over to-be-compare Samples
	for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	    sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	  SampleLabel sLabel1 = (*sCIt)->label1();
	  SampleLabel sLabel2 = (*sCIt)->label2();

	  if( etaBin->hasKValue(sLabel1,sLabel2,fitResType) ) {

	    // MC truth for comparison
	    TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthResoFunc");
	    mcTruth->SetLineColor(kRed);
	    mcTruth->SetLineStyle(2);

	    // Scaled MC truth and uncertainty bands
	    TF1* scaledMCT = etaBin->scaledMCTruthResoFunc("scaledMCTruthResoFunc");
	    scaledMCT->SetLineColor(kRed);
	    TGraphAsymmErrors* scaledMCTBand = etaBin->scaledMCTruthUncertaintyBand();
	    TGraphAsymmErrors* scaledMCTRatioBand = etaBin->scaledMCTruthRatioBand();

	    // Bias corrected measurement and ratio to scaled MC truth
	    TGraphAsymmErrors* biasCorrRes = etaBin->biasCorrectedResolution(sLabel1,sLabel2,fitResType);
	    setStyle(etaBin->sampleType(sLabel1),biasCorrRes);
	    TGraphAsymmErrors* biasCorrResRatio = util::HistOps::createRatioGraph(biasCorrRes,scaledMCT);

	    // Labels
	    TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,labelMk_->start(),label->GetSize());
	    leg->AddEntry(biasCorrRes,"Bias Corrected "+sLabel1,"P");
	    leg->AddEntry(scaledMCT,"Scaled MC Truth","L");
	    leg->AddEntry(scaledMCTBand,"Systematic Uncertainty","F");
	    leg->AddEntry(mcTruth,"MC Truth","L");

	    // Create frames for main and ratio plots
	    TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	    hFrameMain->SetTitle(title_);
	    hFrameMain->GetXaxis()->SetMoreLogLabels();
	    hFrameMain->GetXaxis()->SetNoExponent();

	    TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	    hFrameRatio->GetXaxis()->SetMoreLogLabels();
	    hFrameRatio->GetXaxis()->SetNoExponent();
	     
	    out_->nextMainPad("Bias-corrected "+sLabel1+" vs scaled MC truth "+etaBin->toString());
	    hFrameMain->Draw();
	    scaledMCTBand->Draw("E3same");
	    mcTruth->Draw("same");
	    scaledMCT->Draw("same");
	    biasCorrRes->Draw("PE1same");
	    label->Draw("same");
	    leg->Draw("same");
	    gPad->RedrawAxis();
	    out_->nextRatioPad();
	    hFrameRatio->Draw();
	    scaledMCTRatioBand->Draw("E3same");
	    hFrameRatio->Draw("same");
	    biasCorrResRatio->Draw("PE1same");
	    out_->logx();
	    gPad->RedrawAxis();
	    out_->saveCurrentPad(histFileName("ScaledMCTruthResolution",etaBin,sLabel1,fitResType));

	    delete biasCorrRes;
	    delete biasCorrResRatio;
	    delete mcTruth;
	    delete scaledMCT;
	    delete scaledMCTBand;
	    delete scaledMCTRatioBand;
	    delete label;
	    delete leg;
	    delete hFrameMain;
	    delete hFrameRatio;
	  } // End if hasKValue
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins
     
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotScaledMCTruth(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotSlopes() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotSlopes(): Entering" << std::endl;
    std::cout << "Plotting extrapolation slopes" << std::endl;
     
    // +++++ Slopes per Sample ++++++++++++++++++++++++++++++++++++++
    out_->newPage("Slope");

    unsigned int workingPointPtSoftIdx = 3;

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

	  // Graph of kSoft slopes
	  TGraphAsymmErrors* gSlope = etaBin->kSoftSlope(sampleLabel,fitResType);
	  setStyle(sampleType,gSlope);

	  // kSoft fitting function
	  TF1* fit = etaBin->kSoftFit(sampleLabel,*rIt,"PlotMaker::plotSlopes::KSoftFit");
	  fit->SetLineColor(kRed);

	  // Frame
	  TH1* hFrame = new TH1D("hFrame",title_+";p^{ref}_{T} (GeV);#sigma(p^{rel}_{T} = 0) / #sigma(p^{rel}_{T} = "+util::toTString(par_->ptSoftMax(workingPointPtSoftIdx))+")",1000,xMinPt_,xMaxPt_);
	  hFrame->GetYaxis()->SetRangeUser(0.61,1.29);
	  hFrame->GetXaxis()->SetMoreLogLabels();
	  hFrame->GetXaxis()->SetNoExponent();
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin());

	  out_->nextMultiPad(sampleLabel+": Slope "+etaBin->toString());
	  hFrame->Draw();
	  gSlope->Draw("PE1same");
	  fit->Draw("same");
	  label->Draw("same");
	  out_->logx();
	  out_->saveCurrentPad(histFileName("Slope",etaBin,sampleLabel,fitResType));

	  delete gSlope;
	  delete fit;
	  delete hFrame;
	  delete label;
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins


    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotSlopes(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotSystematicUncertainties() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotSystematicUncertainties(): Entering" << std::endl;
    std::cout << "Plotting systematic uncertainties" << std::endl;
    out_->newPage("Systematic Uncertainties");

    // +++++ Relative systematic uncertainties ++++++++++++++++++++++++++++++++++++++
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

	  // If uncertainties are assigned to this sample and 
	  // fit result type, get them
	  const SystematicUncertainty* uncert = 0;
	  if( etaBin->findSystematicUncertainty(sampleLabel,fitResType,uncert) ) {

	    // Uncertainty bands and labels for components
	    // and total uncertainty
	    std::vector<TGraphAsymmErrors*> bands;
	    TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(uncert->nComponents()+1,-labelMk_->start(),label->GetSize());
	    // Add total uncertainty
	    bands.push_back(uncert->relUncertSteps());
	    for(SystUncertIt it = uncert->componentsBegin(); it != uncert->componentsEnd(); ++it) {
	      // Add components
	      bands.push_back((*it)->relUncertSteps());
	      leg->AddEntry(bands.back(),(*it)->label(),"F");
	    }
	    leg->AddEntry(bands.front(),uncert->label(),"F");

	    // Create frame
	    TH1* hFrame = new TH1D("hFrame","",1000,xMinPt_,xMaxPt_);
	    hFrame->GetXaxis()->SetMoreLogLabels();
	    hFrame->GetXaxis()->SetNoExponent();
	    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
	      hFrame->SetBinContent(bin,0.);
	    }
	    hFrame->SetLineStyle(2);
	    hFrame->GetYaxis()->SetRangeUser(-0.39,0.89);
	    util::HistOps::setAxisTitles(hFrame,"p^{ref}_{T}","GeV","Relative Uncertainty");
	     
	    out_->nextPad(sampleLabel+": Systematic uncertainties "+etaBin->toString());
	    hFrame->Draw();
	    for(std::vector<TGraphAsymmErrors*>::iterator it = bands.begin();
		it != bands.end(); ++it) {
	      (*it)->Draw("E2same");
	    }
	    hFrame->Draw("same");
	    label->Draw("same");
	    leg->Draw("same");
	    out_->logx();
	    out_->saveCurrentPad(histFileName("RelativeSystematicUncertainty",etaBin,sampleLabel,fitResType));

	    for(std::vector<TGraphAsymmErrors*>::iterator it = bands.begin();
		it != bands.end(); ++it) {
	      delete *it;
	    }
	    delete hFrame;
	    delete label;
	    delete leg;
	  } // End if has uncertainties
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins


    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotSystematicUncertainties(): Leaving" << std::endl;    
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const EtaBin* etaBin, SampleLabel label1, SampleLabel label2, FitResult::Type type) const {
    return par_->outFilePrefix()+"_"+id+"_"+(label1.ReplaceAll(" ","_"))+"_vs_"+(label2.ReplaceAll(" ","_"))+FitResult::toString(type)+"_EtaBin"+util::toTString(etaBin->etaBin())+".eps";
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
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, FitResult::Type type) const {
    return par_->outFilePrefix()+"_"+id+"_AllSamples_"+FitResult::toString(type)+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+".eps";
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample) const {
    return par_->outFilePrefix()+"_"+id+"_AllTypes_"+((sample->label()).ReplaceAll(" ","_"))+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+".eps";
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, unsigned int ptSoftBinIdx) const {
    return par_->outFilePrefix()+"_"+id+"_"+((sample->label()).ReplaceAll(" ","_"))+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx)+".eps";
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, SampleLabel label1, SampleLabel label2, FitResult::Type type, unsigned int ptSoftBinIdx) const {
    return par_->outFilePrefix()+"_"+id+"_"+(label1.ReplaceAll(" ","_"))+"_vs_"+(label2.ReplaceAll(" ","_"))+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx)+".eps";
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::setStyle(Sample::Type type, TH1* &h) const {
    h->SetMarkerStyle(markerStyle(type));
    h->SetMarkerColor(color(type));
    h->SetLineColor(h->GetMarkerColor());
    h->SetTitle(title_);
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::setStyle(Sample::Type type, TGraphAsymmErrors* &g) const {
    g->SetMarkerStyle(markerStyle(type));
    g->SetMarkerColor(color(type));
    g->SetLineColor(g->GetMarkerColor());
    g->SetTitle(title_);
  }


  // -------------------------------------------------------------------------------------
  int PlotMaker::markerStyle(Sample::Type type) const {
    int style = 24;
    if( type == Sample::Data ) style = 20;
    else if( type == Sample::MC ) style = 25;

    return style;
  }


  // -------------------------------------------------------------------------------------
  int PlotMaker::color(Sample::Type type) const {
    int color = 1;
    if( type == Sample::Data ) color = 1;
    else if( type == Sample::MC ) color = 4;

    return color;
  }



  // -------------------------------------------------------------------------------------
  PlotMaker::LabelMaker::LabelMaker(const Parameters* par)
    : par_(par) {}


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptSoftBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = etaBin(label,etaBinIdx,1);
    lab->AddText(ptRange(etaBinIdx,ptBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptSoftBin(unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = etaBin(etaBinIdx,1);
    lab->AddText(ptRange(etaBinIdx,ptBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptBin(unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    TPaveText* lab = etaBin(etaBinIdx,1);
    lab->AddText(ptRange(etaBinIdx,ptBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    TPaveText* lab = etaBin(label,etaBinIdx,1);
    lab->AddText(ptRange(etaBinIdx,ptBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaBin(SampleLabel label, unsigned int etaBinIdx, unsigned int nExtraEntries) const {
    TPaveText* lab = util::LabelFactory::createPaveText(2+nExtraEntries);
    lab->AddText(label+",  #sqrt{s} = 7 TeV, L = "+util::toTString(33)+" pb^{-1}");
    lab->AddText(jets()+",  "+etaRange(etaBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaBin(unsigned int etaBinIdx, unsigned int nExtraEntries) const {
    TPaveText* lab = util::LabelFactory::createPaveText(2+nExtraEntries);
    lab->AddText("#sqrt{s} = 7 TeV, L = "+util::toTString(33)+" pb^{-1}");
    lab->AddText(jets()+",  "+etaRange(etaBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::etaRange(unsigned int etaBinIdx) const {
    return util::toTString(par_->etaMin(etaBinIdx))+" < |#eta| < "+util::toTString(par_->etaMax(etaBinIdx));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptRange(unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    return util::toTString(par_->ptMin(etaBinIdx,ptBinIdx))+" < "+pt()+" < "+util::toTString(par_->ptMax(etaBinIdx,ptBinIdx))+" GeV";
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptSoftRange(unsigned int ptSoftBinIdx) const {
    return ptSoft()+" < "+util::toTString(par_->ptSoftMax(ptSoftBinIdx));
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::jets() const {
    TString lab = "";
    if( par_->jetAlgo() == JetProperties::AK5 ) {
      lab += "Anti-k_{T} (R=0.5)";
    }
    lab += " ";
    if( par_->jetType() == JetProperties::Calo ) lab += "Calo Jets";
    else if( par_->jetType() == JetProperties::PF ) lab += "PF Jets";
    else if( par_->jetType() == JetProperties::JPT ) lab += "JPT Jets";

    return lab;
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::pt() const { 
    return "p^{ave}_{T}";
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptSoft() const { 
    return "p^{rel}_{T,3}";
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::label(FitResult::Type type) const { 
    TString lab = "";
    if( FitResult::validType(type) ) {
      if( type == FitResult::FullMaxLikeRel ) lab = "Max Like (Relative)";
      else if( type == FitResult::FullMaxLikeAbs ) lab = "Max Like (Absolute)";
      else if( type == FitResult::SimpleMaxLike ) lab = "Simplified Max Like";
      else if( type == FitResult::PtAsym ) lab = "Asymmetry";
      else if( type == FitResult::PtGenAsym ) lab = "PLI";
    }
    
    return lab;
  }
}
