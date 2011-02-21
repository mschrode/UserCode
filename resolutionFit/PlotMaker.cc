// $Id: PlotMaker.cc,v 1.3 2011/02/18 18:42:22 mschrode Exp $

#include "PlotMaker.h"

#include <cmath>
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
      yMinExtraRes_(1E-3), yMaxExtraRes_(0.38), yMinResRatio_(0.74), yMaxResRatio_(1.26) {

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
	      
	      // Asymmetry fits
	      std::vector<TF1*> fits;
	      std::vector<TString> labels;

	      // Loop over FitResultTypes
	      for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); 
		  rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {

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
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(3+labels.size(),labelMk_->start(),label->GetSize());
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
	       sample->setStyle(gVals);
	       
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

     if( etaBins_.front()->nSampleTypes() > 1 ) {
       // Loop over eta and pt bins
       for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	 const EtaBin* etaBin = *etaBinIt;
	 out_->newPage("Extrapolation");
	 
	 for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	   const PtBin* ptBin = *ptBinIt;
	   if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
	   
	   // Loop over FitResultTypes
	   for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); 
	       rIt != etaBin->fitResultTypesEnd(); ++rIt) {
	     
	     if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(*rIt) << std::endl;
	     
	     std::vector<TGraphAsymmErrors*> gVals;
	     double valMin = 1000.;
	     double valMax = 0.;
	     std::vector<TF1*> fits;
	     std::vector<TF1*> extra;
	     std::vector<TString> labels;
	     
	     // Loop over Samples
	     for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
	       const Sample* sample = sIt->second;

	       labels.push_back(sample->label());
	     
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
	       sample->setStyle(gVals.back());
	     
	       // Extrapolation function
	       fits.push_back(sample->extrapolationFunction(*rIt,"Extrapolation_"+((sample->label()).ReplaceAll(" ","_"))+"_FitRange"));
	       fits.back()->SetLineColor(gVals.back()->GetMarkerColor());
	       extra.push_back(static_cast<TF1*>(fits.back()->Clone("Extrapolation_"+((sample->label()).ReplaceAll(" ","_"))+"_PlotRange")));
	       extra.back()->SetRange(0.,1.4*ptSoftx.back());
	       extra.back()->SetLineStyle(2);

	       if( val.front() < valMin ) valMin = val.front();
	       if( val.back() > valMax ) valMax = val.back();
	     } // End of loop over Samples
	   
	     // Create frame
	     TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame","",
				    1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	     util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"","#sigma / p^{ref}_{T}");
	     hFrame->GetYaxis()->SetRangeUser(0.8*valMin,1.3*valMax);
	     
	     // Labels
	     TPaveText* label = labelMk_->ptBin(etaBin->etaBin(),ptBin->ptBin());
	     TLegend* leg = util::LabelFactory::createLegendColWithOffset(labels.size(),labelMk_->start(),label->GetSize());
	     for(size_t i = 0; i < labels.size(); ++i) {
	       leg->AddEntry(gVals.at(i),labels.at(i),"P");
	     }
	   
	     // Linear scale
	     out_->nextMultiPad("Sample comparison: Extrapolation "+ptBin->toTString());
	     hFrame->Draw();
	     for(size_t i = 0; i < gVals.size(); ++i) {
	       gVals.at(i)->Draw("PE1same");
	       extra.at(i)->Draw("same");
	       fits.at(i)->Draw("same");
	     }
	     label->Draw("same");
	     leg->Draw("same");
	     out_->saveCurrentPad(histFileName("Extrapolation",ptBin,*rIt));
	 
	     for(size_t i = 0; i < gVals.size(); ++i) {
	       delete gVals.at(i);
	       delete extra.at(i);
	       delete fits.at(i);
	     }
	     delete hFrame;
	     delete label;
	     delete leg;
	   } // End of loop over FitResultTypes
	 } // End of loop over pt bins
       } // End of loop over eta bins
     } // End if nSampleTypes() > 1
     



     // +++++ Extrapolation plots, comparison of different FitResult types +++++++++++++++++++++++++++++++

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
     } // End if nSampleTypes() > 1
     
     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotExtrapolation(): Leaving" << std::endl;
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

	   // Ratio of measurement and mc truth
	   TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gResCorr,mcTruth);

	   // PLI function
	   TF1* pli = etaBin->pliFunc("PLI_Eta"+util::toTString(etaBinIdx));
	   pli->SetLineColor(kGreen);
	   pli->SetLineStyle(2);

	   // Create frames for main and ratio plots
	   TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	   hFrameMain->GetXaxis()->SetMoreLogLabels();
	   hFrameMain->GetXaxis()->SetNoExponent();

	   TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	   hFrameRatio->GetXaxis()->SetMoreLogLabels();
	   hFrameRatio->GetXaxis()->SetNoExponent();
	     
	   // Labels
	   TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin());
	   TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,labelMk_->start(),label->GetSize());
	   leg->AddEntry(gRes,"Resolution ("+FitResult::toString(fitResType)+")","P");
	   leg->AddEntry(pli,"Particle Level Imbalance","L");
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
	   std::vector<double> pt1;
	   std::vector<double> ptErr1;
	   std::vector<double> resCorr1;
	   std::vector<double> resCorrStatErr1;
	   std::vector<double> pt2;
	   std::vector<double> ptErr2;
	   std::vector<double> resCorr2;
	   std::vector<double> resCorrStatErr2;
	   for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	     unsigned int ptBinIdx = (*ptBinIt)->ptBin();
	     pt1.push_back(etaBin->meanPt(sLabel1,fitResType,ptBinIdx));
	     ptErr1.push_back(0.);
	     resCorr1.push_back(etaBin->correctedResolution(sLabel1,fitResType,ptBinIdx));
	     resCorrStatErr1.push_back(etaBin->correctedResolutionStatUncert(sLabel1,fitResType,ptBinIdx));

	     pt2.push_back(etaBin->meanPt(sLabel2,fitResType,ptBinIdx));
	     ptErr2.push_back(0.);
	     resCorr2.push_back(etaBin->correctedResolution(sLabel2,fitResType,ptBinIdx));
	     resCorrStatErr2.push_back(etaBin->correctedResolutionStatUncert(sLabel2,fitResType,ptBinIdx));
	   }
	   TGraphAsymmErrors* gRes1 = new TGraphAsymmErrors(pt1.size(),&(pt1.front()),&(resCorr1.front()),
					 &(ptErr1.front()),&(ptErr1.front()),
					 &(resCorrStatErr1.front()),&(resCorrStatErr1.front()));
	   gRes1->SetMarkerStyle(markerStyleCorrectedResolution(etaBin->sampleType(sLabel1)));

	   TGraphAsymmErrors* gRes2 = new TGraphAsymmErrors(pt2.size(),&(pt2.front()),&(resCorr2.front()),
					 &(ptErr1.front()),&(ptErr2.front()),
					 &(resCorrStatErr2.front()),&(resCorrStatErr2.front()));
	   gRes2->SetMarkerStyle(markerStyleCorrectedResolution(etaBin->sampleType(sLabel2)));

	   TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gRes1,gRes2);

	   // Create frames for main and ratio plots
	   TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	   hFrameMain->GetXaxis()->SetMoreLogLabels();
	   hFrameMain->GetXaxis()->SetNoExponent();

	   TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	     
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
    return "p^{rel}_{||,3}";
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
