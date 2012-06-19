// $Id: PlotMaker.cc,v 1.38 2012/06/16 15:02:29 mschrode Exp $

#include "PlotMaker.h"

#include <cmath>
#include <iostream>

#include "TArrow.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TPad.h"

#include "FitResult.h"
#include "ResolutionFunction.h"
#include "SystematicUncertainty.h"

#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"
#include "../util/MultiCanvas.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  PlotMaker::PlotMaker(const Parameters *par, const EtaBins &etaBins)
    : par_(par), etaBins_(etaBins),
      xMinPt_(std::max(1.,0.9*par->ptMin(0))),
      xMaxPt_(1.1*par->ptMax(0)),
      yMinExtraRes_(1E-3),
      yMaxExtraRes_(util::StyleSettings::style() == util::StyleSettings::Presentation ? 0.47 : 0.28),
      yMinResRatio_(0.74),
      yMaxResRatio_(1.26) {

    // Output manager
    out_ = OutputManager::createOutputManager(par_->outMode(),par_->outFilePrefix());

    // Style
    labelMk_ = new LabelMaker(par_);

    lineWidth_ = 1;
    markerSize_ = 1.;
    if( util::StyleSettings::style() == util::StyleSettings::Presentation ) {
      if( par_->outMode() == OutputManager::EPSSingleFiles ) {
	lineWidth_ = 3;
	markerSize_ = 1.5;
      } else {
	lineWidth_ = 2;
	markerSize_ = 0.9;
      }
    } else if( util::StyleSettings::style() == util::StyleSettings::Note ) {
	lineWidth_ = 2;
	markerSize_ = 1.4;
    }

    if( util::StyleSettings::style() == util::StyleSettings::PAS )
      title_ = util::StyleSettings::title(par_->lumi(),true);
    else
      title_ = "";
  }


  // -------------------------------------------------------------------------------------
  PlotMaker::~PlotMaker() {
    delete out_;
    delete labelMk_;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::makeAllPlots() const {
//     plotControlDistributions();
//     plotPtGenSpectra();
    //    plotAsymmetry();
    //    plotExtrapolation();
    //    plotParticleLevelImbalance();
    plotResolution();
    plotScaledMCTruth();
    plotSystematicUncertainties();

//    plotAsymmetryComparison();

//    plotPtSpectra();
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::plotAsymmetry() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetry(): Entering" << std::endl;
    std::cout << "Plotting asymmetry" << std::endl;

    const double asymMax = 0.34;

    // +++++ Asymmetry plots per bin, with different FitResult types ++++++++++++++++++++++++++++++

//     // Loop over SampleLabels
//     for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
// 	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      
//       SampleLabel sampleLabel = sTIt->first;
//       if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
      
//       // Loop over ptSoft bins
//       for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
// 	if( par_->verbosity() > 1 ) std::cout << "      PtSoftBin " << ptSoftBinIdx << std::endl;
	
// 	// Loop over eta bins
// 	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
// 	  const EtaBin* etaBin = *etaBinIt;
	  
// 	  out_->newPage("Asymmetry");
	  
// 	  // Loop over pt bins
// 	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
// 	    const PtBin* ptBin = *ptBinIt;
// 	    if( par_->verbosity() > 1 ) std::cout << "    " << ptBin->toTString() << std::endl;
	    
// 	    // Loop over Samples and select current one (by sampleLabel)
// 	    for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
// 	      const Sample* sample = sIt->second;
// 	      if( sample->label() != sampleLabel ) continue;
	      
// 	      // Asymmetry distribution
// 	      TH1* hPtAsym = sample->histPtAsym(ptSoftBinIdx);
// 	      util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","events",false);
// 	      hPtAsym->GetXaxis()->SetRangeUser(-asymMax,asymMax);
// 	      setStyle(sample,hPtAsym);
	      
// 	      // Asymmetry fits
// 	      std::vector<TF1*> fits;
// 	      std::vector<TString> labels;

// 	      // Loop over FitResultTypes
// 	      for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin();
// 		  rIt != etaBin->fitResultTypesEnd(); ++rIt) {

// 		double sigma = sample->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
// 		fits.push_back(new TF1("AsymmetryFit_"+FitResult::toString(*rIt),"gaus",-1.,1.));
// 		fits.back()->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
// 		fits.back()->SetParameter(1,0.);
// 		fits.back()->SetParameter(2,sigma);

// 		labels.push_back("Fit: "+FitResult::toString(*rIt));
// 	      } // End of loop over FitResultTypes

// 	      for(unsigned int i = 0; i < fits.size(); ++i) {
// 		fits.at(i)->SetLineWidth(1);
// 		fits.at(i)->SetLineStyle(1+i);
// 		fits.at(i)->SetLineColor(1+i);
// 	      }
	      
// 	      // Labels
// 	      TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
// 	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(1,-labelMk_->start(),label->GetSize());
// 	      leg->AddEntry(hPtAsym,sample->label()+" (N = "+util::toTString(hPtAsym->Integral())+")","P");
// 	      // 	      for(unsigned int i = 0; i < fits.size(); ++i) {
// 	      // 		leg->AddEntry(fits.at(i),labels.at(i),"L");
// 	      // 	      }
	    
// 	      //util::HistOps::setYRange(hPtAsym,label->GetSize()+leg->GetNRows()+1);
// 	      util::HistOps::setYRange(hPtAsym,label->GetSize()+1);
// 	      out_->nextMultiPad(sample->label()+": PtAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
// 	      hPtAsym->Draw("PE1");
// 	      for(unsigned int i = 0; i < fits.size(); ++i) {
// 		//fits.at(i)->Draw("same");
// 	      }
// 	      label->Draw("same");
// 	      leg->Draw("same");
// 	      out_->saveCurrentPad(histFileName("PtAsym",ptBin,sample,ptSoftBinIdx));

// 	      delete hPtAsym;
// 	      for(unsigned int i = 0; i < fits.size(); ++i) {
// 		delete fits.at(i);
// 	      }
// 	      delete label;
// 	      delete leg;
// 	    } // End of loop over Samples
// 	  } // End of loop over pt bins
// 	} // End of loop over eta bins
//       } // End of loop over ptSoft bins
//     } // End of loop over SampleLabels



    // +++++ Asymmetry plots per bin for different pt3 ++++++++++++++++++++++++++++++

    // Loop over SampleLabels
    for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {

      SampleLabel sampleLabel = sTIt->first;

      // Plot spectra only for MC sample that is
      // compared to the data
      bool isToBeCompared = false;
      for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	  sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	if( (*sCIt)->contains(sampleLabel) ) {
	  isToBeCompared = true;
	  break;
	}
      }

      if( isToBeCompared ) {
	if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
      
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

	      // Asymmetry distribution for different ptSoft
	      std::vector<TH1*> hPtAsyms;
	      std::vector<TGraphAsymmErrors*> gUncerts;
	      std::vector<TString> labels;
	      unsigned int step = (par_->nPtSoftBins() > 5 ? 4 : 1);
	      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ptSoftBinIdx += step) {
		TH1* h = sample->histPtAsym(ptSoftBinIdx);
		util::HistOps::normHist(h,"width");
		hPtAsyms.push_back(h);
		gUncerts.push_back(util::HistOps::getUncertaintyBand(h));
		labels.push_back(labelMk_->ptSoftRange(ptSoftBinIdx));
	      }

	      // Labels and style
	      TPaveText* label = labelMk_->etaPtAveBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(labels.size(),-0.6,label->GetSize());

	      for(unsigned int i = 0; i < hPtAsyms.size(); ++i) {
		util::HistOps::setAxisTitles(hPtAsyms.at(i),"Asymmetry","","events",true);
		hPtAsyms.at(i)->GetXaxis()->SetRangeUser(-0.28,0.28);
		setStyle(sample,hPtAsyms.at(i));
		hPtAsyms.at(i)->SetLineColor(util::StyleSettings::color(i));
		hPtAsyms.at(i)->SetLineStyle(1+i);
		hPtAsyms.at(i)->SetMarkerStyle(0);
		util::HistOps::setYRange(hPtAsyms.at(i),label->GetSize()+leg->GetNRows()+2);
		leg->AddEntry(hPtAsyms.at(i),labels.at(i),"L");
		gUncerts.at(i)->SetFillColor(util::StyleSettings::colorLight(hPtAsyms.at(i)->GetLineColor()));
	      }

	      out_->nextMultiPad(sample->label()+": PtAsym "+ptBin->toTString());
	      hPtAsyms.front()->Draw("HIST");
	      for(unsigned int i = 0; i < hPtAsyms.size(); ++i) {
		gUncerts.at(i)->Draw("E3same");
		hPtAsyms.at(i)->Draw("HISTsame");
	      }
	      label->Draw("same");
	      leg->Draw("same");
	      gPad->RedrawAxis();
	      out_->saveCurrentPad(histFileName("PtAsym",ptBin,sample));
	    
	      for(unsigned int i = 0; i < hPtAsyms.size(); ++i) {
		delete hPtAsyms.at(i);
		delete gUncerts.at(i);
	      }
	      delete label;
	      delete leg;
	    } // End isToBeCompared
	  } // End of loop over Samples
	} // End of loop over pt bins
      } // End of loop over eta bins
    } // End of loop over SampleLabels



    // +++++ Asymmetry plots per bin, different Samples ++++++++++++++++++++++++++++++
    
 //    if( etaBins_.front()->nComparedSamples() > 0 ) {
//       // Loop over FitResultTypes
//       for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); 
// 	  rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
// 	if( par_->verbosity() > 1 ) std::cout << "  " << FitResult::toString(*rIt) << std::endl;

// 	// Loop over ptSoft bins
// 	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
// 	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	  
// 	  // Loop over eta and bins
// 	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
// 	    const EtaBin* etaBin = *etaBinIt;
// 	    out_->newPage("Asymmetry");
	    
// 	    // Loop over pt bins
// 	    for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
// 	      const PtBin* ptBin = *ptBinIt;
// 	      if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;

// 	      // Loop over to-be-compare Samples
// 	      for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
// 		  sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
//  		SampleLabel sLabel1 = (*sCIt)->label1();
//  		SampleLabel sLabel2 = (*sCIt)->label2();
// 		const Sample* sample1 = ptBin->findSample(sLabel1);
// 		const Sample* sample2 = ptBin->findSample(sLabel2);

// 		// Asymmetry distributions and fits
// 		TH1* hPtAsym1 = sample1->histPtAsym(ptSoftBinIdx);
// 		if( hPtAsym1->Integral() ) hPtAsym1->Scale(1./hPtAsym1->Integral("width"));
// 		util::HistOps::setAxisTitles(hPtAsym1,"Asymmetry","","events",true);
// 		hPtAsym1->GetXaxis()->SetRangeUser(-asymMax,asymMax);
// 		setStyle(sample1,hPtAsym1);

// 		double sigma = sample1->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
// 		TF1* fit1 = new TF1("AsymmetryFit_Sample1","gaus",-1.,1.);
// 		fit1->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
// 		fit1->SetParameter(1,0.);
// 		fit1->SetParameter(2,sigma);
// 		fit1->SetLineWidth(1);
// 		fit1->SetLineStyle(1);
// 		fit1->SetLineColor(hPtAsym1->GetLineColor());

// 		TH1* hPtAsym2 = sample2->histPtAsym(ptSoftBinIdx);
// 		if( hPtAsym2->Integral() ) hPtAsym2->Scale(1./hPtAsym2->Integral("width"));
// 		util::HistOps::setAxisTitles(hPtAsym2,"Asymmetry","","events",true);
// 		hPtAsym2->GetXaxis()->SetRangeUser(-asymMax,asymMax);
// 		hPtAsym2->SetLineColor(kBlack);
// 		hPtAsym2->SetLineWidth(lineWidth_);
// 		hPtAsym2->SetMarkerStyle(0);
// 		hPtAsym2->SetFillColor(38);


// 		sigma = sample2->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
// 		TF1* fit2 = new TF1("AsymmetryFit_Sample2","gaus",-1.,1.);
// 		fit2->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
// 		fit2->SetParameter(1,0.);
// 		fit2->SetParameter(2,sigma);
// 		fit2->SetLineWidth(1);
// 		fit2->SetLineStyle(2);
// 		fit2->SetLineColor(hPtAsym2->GetLineColor());
			      
// 		// Labels
// 		TPaveText* label = labelMk_->ptSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
// 		TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.65*labelMk_->start(),label->GetSize());
// 		leg->AddEntry(hPtAsym1,labelMk_->sample(sLabel1),"P");
// 		leg->AddEntry(hPtAsym2,labelMk_->sample(sLabel2),"F");
	    
// 		util::HistOps::setYRange(hPtAsym1,label->GetSize()+leg->GetNRows()+1);
// 		out_->nextMultiPad(sLabel1+" vs "+sLabel2+": PtAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
// 		hPtAsym1->Draw("PE1");
// 		hPtAsym2->Draw("HISTsame");
// 		hPtAsym1->Draw("PE1same");
// // 		fit1->Draw("same");
// // 		fit2->Draw("same");
// 		label->Draw("same");
// 		leg->Draw("same");
// 		gPad->RedrawAxis();
// 		out_->saveCurrentPad(histFileName("PtAsym",ptBin,sLabel1,sLabel2,*rIt,ptSoftBinIdx));

// 		delete hPtAsym1;
// 		delete hPtAsym2;
// 		delete fit1;
// 		delete fit2;
// 		delete label;
// 		delete leg;
// 	      } // End of loop over to-be-compared samples
// 	    } // End of loop over pt bins
// 	  } // End of loop over eta bins
// 	} // End of loop over ptSoft bins
//       } // End of loop over FitResult types
//     } // End if samples to be compared

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetry(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  void PlotMaker::plotAsymmetryComparison() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetryComparison(): Entering" << std::endl;
    std::cout << "Plotting comparison of asymmetry" << std::endl;

    const double asymMax = 0.44;

    // +++++ Asymmetry plots for different samples ++++++++++++++++++++++++++++++
    
    // Loop over eta bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      const unsigned int etaBinIdx = etaBin->etaBin();
      if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;
      
      // Loop over to-be-compared Samples
      for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	  sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	SampleLabel sLabel1 = (*sCIt)->label1();
	SampleLabel sLabel2 = (*sCIt)->label2();

	out_->newPage(("Asymmetry: "+sLabel1+" vs "+sLabel2));
	if( par_->verbosity() > 1 ) std::cout << "    Asymmetry: "+sLabel1+" vs "+sLabel2 << std::endl;

	// Loop over ptSoft bins
	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
      
	  // Loop over pt bins
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	    
	    const Sample* sample1 = ptBin->findSample(sLabel1);
	    const Sample* sample2 = ptBin->findSample(sLabel2);

	    // Asymmetry distributions and fits
	    TH1* hPtAsym1 = sample1->histPtAsym(ptSoftBinIdx);
	    if( hPtAsym1->Integral() ) hPtAsym1->Scale(1./hPtAsym1->Integral("width"));
	    util::HistOps::setAxisTitles(hPtAsym1,"Asymmetry","","events",true);
	    hPtAsym1->GetXaxis()->SetRangeUser(-asymMax,asymMax);
	    setStyle(sample1,hPtAsym1);

	    double sigma = sample1->asymmetryWidth(ptSoftBinIdx);
	    TF1* fit1 = new TF1("AsymmetryWidth_Sample1","gaus",-1.,1.);
	    fit1->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
	    fit1->SetParameter(1,0.);
	    fit1->SetParameter(2,sigma);
	    fit1->SetLineWidth(1);
	    fit1->SetLineStyle(1);
	    fit1->SetLineColor(hPtAsym1->GetLineColor());

	    TH1* hPtAsym2 = sample2->histPtAsym(ptSoftBinIdx);
	    if( hPtAsym2->Integral() ) hPtAsym2->Scale(1./hPtAsym2->Integral("width"));
	    util::HistOps::setAxisTitles(hPtAsym2,"Asymmetry","","events",true);
	    hPtAsym2->GetXaxis()->SetRangeUser(-asymMax,asymMax);
	    setStyle(sample2,hPtAsym2);

	    sigma = sample2->asymmetryWidth(ptSoftBinIdx);
	    TF1* fit2 = new TF1("AsymmetryWidth_Sample2","gaus",-1.,1.);
	    fit2->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
	    fit2->SetParameter(1,0.);
	    fit2->SetParameter(2,sigma);
	    fit2->SetLineWidth(1);
	    fit2->SetLineStyle(2);
	    fit2->SetLineColor(hPtAsym2->GetLineColor());
			      
	    // Labels
	    TPaveText* label = labelMk_->etaPtAvePtSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
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
	    out_->saveCurrentPad(histFileName("PtAsymWidth",ptBin,sLabel1,sLabel2,ptSoftBinIdx));
	    
	    delete hPtAsym1;
	    delete hPtAsym2;
	    delete fit1;
	    delete fit2;
	    delete label;
	    delete leg;
	  } // End of loop over eta bins
	} // End of loop over to-be-compared samples
      } // End of loop over ptSoft bins
    } // End of loop over pt bins


    // +++++ Asymmetry width for different samples ++++++++++++++++++++++++++++++
    
    out_->newPage("Asymmetry width vs pt");
    // Loop over eta bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      const unsigned int etaBinIdx = etaBin->etaBin();
      if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;
      
      // Loop over to-be-compared Samples
      for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	  sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	SampleLabel sLabel1 = (*sCIt)->label1();
	SampleLabel sLabel2 = (*sCIt)->label2();

	out_->newPage(("Asymmetry: "+sLabel1+" vs "+sLabel2));
	if( par_->verbosity() > 1 ) std::cout << "    Asymmetry: "+sLabel1+" vs "+sLabel2 << std::endl;

	// Loop over ptSoft bins
	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;

	  // Fitted standard deviation vs pt and quadratic difference
	  TGraphAsymmErrors* gWidthAsym1 = etaBin->asymmetryWidths(sLabel1,ptSoftBinIdx);
	  setStyle(sLabel1,gWidthAsym1);
	  TGraphAsymmErrors* gWidthAsym2 = etaBin->asymmetryWidths(sLabel2,ptSoftBinIdx);
	  setStyle(sLabel2,gWidthAsym2);
	  TGraphAsymmErrors* gDiff2 = util::HistOps::createQuadraticDifferenceGraph(gWidthAsym2,gWidthAsym1);
	  gDiff2->SetMarkerStyle(27);
	  gDiff2->SetMarkerColor(kRed);
	  gDiff2->SetLineColor(kRed);
	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gWidthAsym2,gWidthAsym1);
// 	  gRatio->SetMarkerStyle(27);
// 	  gRatio->SetMarkerColor(kBlack);
// 	  gRatio->SetLineColor(kBlack);
	  	  
	  // Fit to standard deviations vs pt
	  ResolutionFunction* rf = ResolutionFunction::fitTGraph(gWidthAsym1,ResolutionFunction::ModifiedNSC);
	  TF1* fWidthAsym1 = rf->func("Asym1_Fit");
	  delete rf;
	  fWidthAsym1->SetLineWidth(lineWidth_);
	  fWidthAsym1->SetLineColor(gWidthAsym1->GetLineColor());
	  fWidthAsym1->SetLineStyle(1);
	  rf = ResolutionFunction::fitTGraph(gWidthAsym2,ResolutionFunction::ModifiedNSC);
	  TF1* fWidthAsym2 = rf->func("Asym2_Fit");
	  delete rf;
	  fWidthAsym2->SetLineWidth(lineWidth_);
	  fWidthAsym2->SetLineColor(gWidthAsym2->GetLineColor());
	  fWidthAsym2->SetLineStyle(1);
	  rf = ResolutionFunction::fitTGraph(gDiff2,ResolutionFunction::ModifiedNSC);
	  TF1* fDiff2 = rf->func("AsymDiff2_Fit");
	  fDiff2->SetLineWidth(lineWidth_);
	  fDiff2->SetLineColor(gDiff2->GetLineColor());
	  delete rf;

	  TPaveText* label = labelMk_->etaPtSoftBin(etaBin->etaBin(),ptSoftBinIdx);
	  if( Sample::type(sLabel1) == Sample::type(sLabel2) ) {
	    delete label;
	    label = labelMk_->etaPtSoftBin(util::LabelFactory::data(labelMk_->lumi()),etaBin->etaBin(),ptSoftBinIdx);
	  }
	  TLegend* leg = util::LabelFactory::createLegendCol(2,0.3);
 	  leg->AddEntry(gWidthAsym1,sLabel1,"P");
 	  leg->AddEntry(gWidthAsym2,sLabel2,"P");
	  // 	  leg->AddEntry(gDiff2,"#sqrt{B^{2} - A^{2}}","P");
	  // 	  leg->AddEntry(gRatio,"B / A","P");
// 	  leg->AddEntry(gWidthAsym1,sLabel1,"P");
// 	  leg->AddEntry(gWidthAsym2,sLabel2,"P");

	  TH1* hFrame = new TH1D("hFrame",title_,1000,xMinPt_,xMaxPt_);
	  hFrame->GetXaxis()->SetMoreLogLabels();
	  hFrame->GetXaxis()->SetNoExponent();
//  	  for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
//  	    hFrame->SetBinContent(i,0.);
//  	  }
//  	  hFrame->SetLineStyle(2);
	  hFrame->SetLineWidth(lineWidth_);
	  hFrame->GetYaxis()->SetRangeUser(-0.09,0.34);
	  util::HistOps::setAxisTitles(hFrame,"p^{ave}_{T}","GeV","#sigma_{A}");
	  
// 	  out_->nextPad(sLabel1+" vs "+sLabel2+": PtAsym Width, "+etaBin->toString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
// 	  hFrame->Draw();
// 	  //fDiff2->Draw("same");
// 	  //gDiff2->Draw("PE1same");
// 	  //fWidthAsym2->Draw("same");
// 	  gWidthAsym2->Draw("PE1same");
// 	  //fWidthAsym1->Draw("same");
// 	  gWidthAsym1->Draw("PE1same");
// 	  label->Draw("same");
// 	  leg->Draw("same");
// 	  out_->logx();
// 	  out_->saveCurrentPad(histFileName("PtAsymWidth",sLabel1,sLabel2,etaBin,ptSoftBinIdx));

	  out_->nextMainPad(sLabel1+" vs "+sLabel2+": PtAsym Width, "+etaBin->toString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	  TH1* hFrameMain = out_->mainFrame(hFrame);
	  hFrameMain->GetYaxis()->SetRangeUser(0.005,0.27);
	  hFrameMain->SetTitle(title_);	  
	  TH1* hFrameRatio = out_->ratioFrame(hFrame,util::LabelFactory::ptAve(),"GeV","#frac{"+sLabel2+"}{"+sLabel1+"}",0.83,1.17);
	  hFrameRatio->SetLineWidth(lineWidth_);
  	  hFrameMain->GetXaxis()->SetRangeUser(27.,xMaxPt_);
  	  hFrameRatio->GetXaxis()->SetRangeUser(27.,xMaxPt_);
	  hFrameMain->Draw();
	  gWidthAsym2->Draw("PE1same");
	  gWidthAsym1->Draw("PE1same");
	  label->Draw("same");
	  leg->Draw("same");
	  gPad->RedrawAxis();
	  gPad->SetLogx();
	  out_->nextRatioPad();
	  hFrameRatio->Draw("][");
	  gRatio->Draw("PE1same");
	  gPad->SetLogx();
	  out_->saveCurrentPad(histFileName("PtAsymWidth",sLabel1,sLabel2,etaBin,ptSoftBinIdx));

	  delete gWidthAsym1;
	  delete fWidthAsym1;
	  delete gWidthAsym2;
	  delete fWidthAsym2;
	  delete gDiff2;
	  delete fDiff2;
	  delete gRatio;
	  delete hFrame;
	  delete label;
	  delete leg;
	  delete hFrameMain;
	  delete hFrameRatio;
	} // End of loop over ptSoft bins
      } // End of loop over to-be-compared samples
    } // End of loop over eta bins

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetryComparison(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  void PlotMaker::plotAsymmetryTails() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetryTails(): Entering" << std::endl;
    std::cout << "Plotting asymmetry tails" << std::endl;

    // Loop over eta bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      const unsigned int etaBinIdx = etaBin->etaBin();
      if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;
      bool printOutDone = false;

      // Loop over ptSoft bins
      //      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
      for(unsigned int ptSoftBinIdx = 4; ptSoftBinIdx < 5; ++ptSoftBinIdx) {
	if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;

	// Loop over pt bins
	for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	  const PtBin* ptBin = *ptBinIt;
	  if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;

	  // Loop over to-be-compare Samples
	  for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	      sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	    SampleLabel sLabelData = (*sCIt)->label1();
	    SampleLabel sLabelMC = (*sCIt)->label2();
	    const Sample* sData = ptBin->findSample(sLabelData);
	    const Sample* sMC = ptBin->findSample(sLabelMC);

	    // Scale of MC to data
	    double w = sMC->relativeWeightTo(sLabelData,ptSoftBinIdx);

	    // Asymmetry distributions
	    TH1* hAsymData = sData->histPtAsym(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hAsymData,"Asymmetry","","events");
	    
	    TH1* hAsymMC = sMC->histPtAsym(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hAsymMC,"Asymmetry","","events");

	    // Smear histogram by data / MC ratio
	    double kValue = 0.;
	    if( etaBin->hasKValue(sLabelData,sLabelMC,*(etaBin->fitResultTypesBegin())) ) {
	      kValue = etaBin->kValue(sLabelData,sLabelMC,*(etaBin->fitResultTypesBegin()));
	      kValue = std::max(0.,kValue-1.);
	    }	    
	    TH1* hAsymMCSmeared = 0;
	    double width = 0.;
	    double widthErr = 1000.;
	    if( util::HistOps::fitCoreWidth(hAsymMC,2.,width,widthErr) ) {
	      util::HistOps::smearHistogram(hAsymMC,hAsymMCSmeared,width,kValue);
	      if( !printOutDone ){
		std::cout << "  Smearing asymmetry in '" << sLabelMC << "' in " << etaBin->toString() << " with k = " << kValue << std::endl;
		printOutDone = true;
	      }
	    } else {
	      util::HistOps::smearHistogram(hAsymMC,hAsymMCSmeared,width,0.);
	    }
	    
	    hAsymMC->Scale(w);
	    hAsymMCSmeared->Scale(w);

	    setStyle(sLabelData,hAsymData);
	    setStyle(sLabelMC,hAsymMC);	      
	    hAsymMC->SetFillColor(38);
	    hAsymMC->SetLineColor(kBlack);
	    hAsymMC->SetLineWidth(1);
	    setStyle(sLabelMC,hAsymMCSmeared);	      
	    hAsymMCSmeared->SetFillColor(29);
	    hAsymMCSmeared->SetLineColor(kBlack);
	    hAsymMCSmeared->SetLineWidth(1);

	    TH1* hAsymRatio = util::HistOps::createRatioPlot(hAsymData,hAsymMC);
	    TH1* hAsymSmearedRatio = util::HistOps::createRatioPlot(hAsymData,hAsymMCSmeared);


	    // Cumulative distributions
	    TH1* hIntAsymData = util::HistOps::getCumulativeDistributionXToInf(hAsymData);
	    TH1* hIntAsymMC = util::HistOps::getCumulativeDistributionXToInf(hAsymMC);
	    TH1* hIntAsymMCSmeared = util::HistOps::getCumulativeDistributionXToInf(hAsymMCSmeared);
	    TH1* hIntAsymRatio = util::HistOps::createRatioPlot(hIntAsymData,hIntAsymMC);
	    TH1* hIntAsymSmearedRatio = util::HistOps::createRatioPlot(hIntAsymData,hIntAsymMCSmeared);



	    
	    // Labels
	    TPaveText* label = labelMk_->ptSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    TLegend* legNom = util::LabelFactory::createLegendColWithOffset(label->GetSize(),-labelMk_->start(),label->GetSize());
	    legNom->AddEntry(hAsymData,labelMk_->sample(sLabelData),"P");
	    legNom->AddEntry(hAsymMC,labelMk_->sample(sLabelMC)+" (Scaled)","F");
	    TLegend* legSmeared = util::LabelFactory::createLegendColWithOffset(label->GetSize(),-labelMk_->start(),label->GetSize()); 	   
	    legSmeared->AddEntry(hAsymData,labelMk_->sample(sLabelData),"P");
	    legSmeared->AddEntry(hAsymMCSmeared,labelMk_->sample(sLabelMC)+" Smeared (Scaled)","F");


	    if( par_->outMode() == OutputManager::EPSSingleFiles ) {

	      // Asymmetry, linear
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hAsymMC,label->GetSize()+legNom->GetNRows()+1);
	      TH1* hFrameMain = out_->mainFrame(hAsymMC);
	      hFrameMain->SetTitle(title_);	  
	      TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.6);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.6);

	      hFrameMain->Draw();
	      hAsymMC->Draw("HISTsame");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hAsymRatio->Draw("PE1same");
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsym",ptBin,sMC,ptSoftBinIdx));

	      // Asymmetry, log
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hAsymMC,label->GetSize()+legNom->GetNRows(),3E-1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hAsymMC);
	      hFrameMain->SetTitle(title_);	  
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.6);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.6);

	      hFrameMain->Draw();
	      hAsymMC->Draw("HISTsame");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hAsymRatio->Draw("PE1same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymLog",ptBin,sMC,ptSoftBinIdx));

	      // Cumulative distributions
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hIntAsymMC,label->GetSize()+legNom->GetNRows(),3E-1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hIntAsymMC);
	      hFrameMain->SetTitle(title_);	  
	      delete hFrameRatio;
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.6);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.6);

	      hFrameMain->Draw();
	      hIntAsymMC->Draw("HISTsame");
	      hIntAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hIntAsymRatio->Draw("PE1same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymInt",ptBin,sMC,ptSoftBinIdx));



	      // Smeared asymmetry, linear
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hAsymMCSmeared,label->GetSize()+legNom->GetNRows()+1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hAsymMCSmeared);
	      hFrameMain->SetTitle(title_);	  
	      delete hFrameRatio;
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.6);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.6);

	      hFrameMain->Draw();
	      hAsymMCSmeared->Draw("HISTsame");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hAsymSmearedRatio->Draw("PE1same");
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmeared",ptBin,sMC,ptSoftBinIdx));


	      // Smeared asymmetry, log
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hAsymMCSmeared,label->GetSize()+legNom->GetNRows(),3E-1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hAsymMCSmeared);
	      hFrameMain->SetTitle(title_);	  
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.6);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.6);

	      hFrameMain->Draw();
	      hAsymMCSmeared->Draw("HISTsame");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hAsymSmearedRatio->Draw("PE1same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmearedLog",ptBin,sMC,ptSoftBinIdx));


	      // Cumulative distributions
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hIntAsymMCSmeared,label->GetSize()+legNom->GetNRows()+1,3E-1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hIntAsymMCSmeared);
	      hFrameMain->SetTitle(title_);	  
	      delete hFrameRatio;
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.6);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.6);

	      hFrameMain->Draw();
	      hIntAsymMCSmeared->Draw("HISTsame");
	      hIntAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hIntAsymSmearedRatio->Draw("PE1same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmearedInt",ptBin,sMC,ptSoftBinIdx));

	      delete hFrameMain;
	      delete hFrameRatio;

	    } else {

	      // Asymmetry, linear
	      hAsymMC->GetXaxis()->SetRangeUser(0.,0.4);
	      util::HistOps::setYRange(hAsymMC,label->GetSize()+legNom->GetNRows()+1);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hAsymMC->Draw("HIST");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsym",ptBin,sMC,ptSoftBinIdx));
	    
	      // Asymmetry, log
	      hAsymMC->GetXaxis()->SetRangeUser(0.,1.);
	      util::HistOps::setYRange(hAsymMC,label->GetSize()+legNom->GetNRows(),1E-3);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hAsymMC->Draw("HIST");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymLog",ptBin,sMC,ptSoftBinIdx));

	      // Cumulative asymmetry
	      hIntAsymMC->GetXaxis()->SetRangeUser(0.,0.6);
	      util::HistOps::setYRange(hIntAsymMC,label->GetSize()+legNom->GetNRows(),1E-3);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hIntAsymMC->Draw("HIST");
	      hIntAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymInt",ptBin,sMC,ptSoftBinIdx));
	    
	      // Smeared asymmetry, linear
	      hAsymMCSmeared->GetXaxis()->SetRangeUser(0.,0.4);
	      util::HistOps::setYRange(hAsymMCSmeared,label->GetSize()+legSmeared->GetNRows()+1);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hAsymMCSmeared->Draw("HIST");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmeared",ptBin,sMC,ptSoftBinIdx));
	    
	      // Smeared asymmetry, log
	      hAsymMCSmeared->GetXaxis()->SetRangeUser(0.,1.);
	      util::HistOps::setYRange(hAsymMCSmeared,label->GetSize()+legNom->GetNRows(),1E-3);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hAsymMCSmeared->Draw("HIST");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmearedLog",ptBin,sMC,ptSoftBinIdx));

	      // Cumulative smeared asymmetry
	      hIntAsymMCSmeared->GetXaxis()->SetRangeUser(0.,0.6);
	      util::HistOps::setYRange(hIntAsymMCSmeared,label->GetSize()+legNom->GetNRows(),1E-3);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hIntAsymMCSmeared->Draw("HIST");
	      hIntAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmearedInt",ptBin,sMC,ptSoftBinIdx));

	    }
	    
	    
	    // Clean up
	    delete hAsymData;
	    delete hAsymMC;
	    delete hAsymMCSmeared;
	    delete hIntAsymData;
	    delete hIntAsymMC;
	    delete hIntAsymMCSmeared;
	    delete hIntAsymRatio;
	    delete hIntAsymSmearedRatio;
	    delete label;
	    delete legNom;
	    delete legSmeared;
	  } // End of loop over to-be-compared samples
	} // End of loop over pt bins

	out_->newPage("Data vs MC");      

	// Loop over pt bins
	for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	  const PtBin* ptBin = *ptBinIt;
	  if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	  
	  // Loop over to-be-compare Samples
	  for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	      sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	    SampleLabel sLabelData = (*sCIt)->label1();
	    SampleLabel sLabelMC = (*sCIt)->label2();
	    const Sample* sData = ptBin->findSample(sLabelData);
	    const Sample* sMC = ptBin->findSample(sLabelMC);
	    
	    // Scale of MC to data
	    double w = sMC->relativeWeightTo(sLabelData,ptSoftBinIdx);
	    
	    // PtAve spectra
	    TH1* hPtAveData = sData->histPtAve(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hPtAveData,"p^{ave}_{T}","GeV","events");

	    TH1* hPtAveMC = sMC->histPtAve(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hPtAveMC,"p^{ave}_{T}","GeV","events");
	    hPtAveMC->Scale(w);

	    setStyle(sLabelData,hPtAveData);
	    setStyle(sLabelMC,hPtAveMC);
	    hPtAveMC->SetFillColor(38);
	    hPtAveMC->SetLineColor(kBlack);
	    hPtAveMC->SetLineWidth(1);

	    TH1* hPtAveRatio = util::HistOps::createRatioPlot(hPtAveData,hPtAveMC);
	    
	    // Labels
	    TPaveText* label = labelMk_->ptSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    TLegend* legNom = util::LabelFactory::createLegendColWithOffset(label->GetSize(),-labelMk_->start(),label->GetSize());
	    legNom->AddEntry(hPtAveData,labelMk_->sample(sLabelData),"P");
	    legNom->AddEntry(hPtAveMC,labelMk_->sample(sLabelMC)+" (Scaled)","F");

	    if( par_->outMode() == OutputManager::EPSSingleFiles ) {
	      // PtAve spectrum
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hPtAveMC,label->GetSize()+legNom->GetNRows()+2);
	      TH1* hFrameMain = out_->mainFrame(hPtAveMC);
	      hFrameMain->SetTitle(title_);	    
	      TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ave}_{T}","GeV",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      int xMin = 1;
	      int xMax = hPtAveData->GetNbinsX();
	      util::HistOps::findXRange(hPtAveData,xMin,xMax);
	      xMin = std::max(1,xMin-5);
	      xMax = std::min(hPtAveData->GetNbinsX(),xMax+5);
	      hFrameMain->GetXaxis()->SetRangeUser(hPtAveData->GetXaxis()->GetBinLowEdge(xMin),hPtAveData->GetXaxis()->GetBinUpEdge(xMax));
	      hFrameRatio->GetXaxis()->SetRangeUser(hPtAveData->GetXaxis()->GetBinLowEdge(xMin),hPtAveData->GetXaxis()->GetBinUpEdge(xMax));
	      hFrameMain->Draw();
	      hPtAveMC->Draw("HISTsame");
	      hPtAveData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      gPad->RedrawAxis();
	      out_->nextRatioPad();
	      hFrameRatio->Draw("][");
	      hPtAveRatio->Draw("PE1same");
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAve",ptBin,sMC,ptSoftBinIdx));
	      delete hFrameMain;
	      delete hFrameRatio;

	    } else {

	      // PtAve spectra
	      util::HistOps::setYRange(hPtAveMC,label->GetSize()+legNom->GetNRows()+1);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hPtAveMC->Draw("HIST");
	      hPtAveData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAve",ptBin,sMC,ptSoftBinIdx));
	    }
	    
	    
	    // Clean up
	    delete hPtAveData;
	    delete hPtAveMC;
	    delete hPtAveRatio;
	    delete label;
	    delete legNom;
	  } // End of loop over to-be-compared samples
	} // End of loop over pt bins

	out_->newPage("Data vs MC");      

      } // End of loop over ptSoft bins
    } // End of loop over eta bins

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotAsymmetryTails: Leaving" << std::endl;
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

	// Plot spectra only for sample that is
	// compared to the data
	bool isToBeCompared = false;
	for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	    sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	  if( (*sCIt)->contains(sampleLabel) ) {
	    isToBeCompared = true;
	    break;
	  }
	}

	if( isToBeCompared ) {

	  if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
	  
  	// Loop over eta and pt bins
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;
	    out_->newPage("Extrapolation");
	    
	    std::vector<double> errExtrapolation;
	    std::vector<double> errFirstPoint;
	    std::vector<double> errLastPoint;	    

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
		sample->valuesInExtrapolation(*rIt,val,uncert);
		std::vector<double> ptSoftx;
		sample->ptSoft(ptSoftx);
		std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
		
		TGraphAsymmErrors* gVals = 
		  new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
					&(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
					&(uncert.front()),&(uncert.front()));
		setStyle(sample,gVals);
		
		// Extrapolation function
		TF1* fit = sample->extrapolationFunction(*rIt,"ExtrapolationForPlotFitRange");
		fit->SetLineColor(kRed);
		fit->SetLineWidth(lineWidth_);
		TF1* extra = static_cast<TF1*>(fit->Clone("ExtrapolationForPlotPlotRange"));
		extra->SetRange(0.,1.4*ptSoftx.back());
		extra->SetLineStyle(2);
		
		// Create frame
		TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame",title_,
				       1000,0.,1.3*par_->ptSoftMax(par_->nPtSoftBins()-1));
		util::HistOps::setAxisTitles(hFrame,"Threshold "+labelMk_->ptSoftMax(),"",sample->labelQuantityInExtrapolation(*rIt));
		hFrame->GetYaxis()->SetRangeUser(0.7*val.front(),1.7*val.back());
		
		
		// Labels
		TPaveText* label = labelMk_->etaPtAveBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
		TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-labelMk_->start(),label->GetSize());
		//leg->AddEntry(gVals,labelMk_->type(*rIt),"P");
		leg->AddEntry(gVals,"Fitted Resolution","P");
		leg->AddEntry(fit,"Extrapolation","L");
		
		// Linear scale
		out_->nextMultiPad(sample->label()+": Extrapolation "+ptBin->toTString());
		hFrame->Draw();
		label->Draw("same");
		leg->Draw("same");
		gVals->Draw("PE1same");
		extra->Draw("same");
		fit->Draw("same");
		//out_->saveCurrentPad(histFileName("Extrapolation",ptBin,sample,*rIt));
		
		// Uncertainties
		errExtrapolation.push_back(fit->GetParError(0));
		errFirstPoint.push_back(sample->uncertInExtrapolation(*rIt,sample->firstPointInExtrapolation(*rIt)));
		errLastPoint.push_back(sample->uncertInExtrapolation(*rIt,sample->nPtSoftBins()-1));

		// Illustration for systematic uncertainties
		TLine* hor = new TLine(0.,extra->Eval(0.),fit->GetXmin(),extra->Eval(0.));
		hor->SetLineWidth(4);
		hor->SetLineColor(kRed+1);
		hor->SetLineStyle(2);
		TArrow* ver = new TArrow(fit->GetXmin(),extra->Eval(0.),fit->GetXmin(),extra->Eval(fit->GetXmin()),0.04,"<>");
		ver->SetLineWidth(hor->GetLineWidth());
		ver->SetLineColor(hor->GetLineColor());
		
		TPaveText *txt = new TPaveText(fit->GetXmin()+0.01,extra->Eval(0.),fit->GetXmin()+0.08,extra->Eval(fit->GetXmin()));
		txt->SetBorderSize(0);
		txt->SetFillColor(0);
		txt->SetTextFont(42);
		txt->SetTextAlign(12);
		txt->SetTextSize(0.08);
		txt->SetTextColor(hor->GetLineColor());
		txt->AddText("50%");
		
		hFrame->GetXaxis()->SetNdivisions(505);
		hFrame->GetXaxis()->SetRangeUser(0.,0.19);
		hFrame->GetYaxis()->SetRangeUser(0.8*val.front(),1.1*val.back());
		fit->SetLineColor(gVals->GetLineColor());
		extra->SetLineColor(gVals->GetLineColor());
		
		out_->nextMultiPad("Illustration: Extrapolation uncertainty "+ptBin->toTString());
		hFrame->Draw();
		label->Draw("same");
		gVals->Draw("PE1same");
		extra->Draw("same");
		fit->Draw("same");
		hor->Draw("same");
		ver->Draw();
		txt->Draw("same");
		//		out_->saveCurrentPad(histFileName("ExtrapolationUncertaintyIllustration",ptBin,sample,*rIt));
		
		
		delete gVals;
		delete extra;
		delete fit;
		delete hFrame;
		delete label;
		delete leg;
		delete hor;
		delete ver;
		delete txt;
	      } // End of loop over Samples
	    } // End of loop over pt bins

	    TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame",title_+";p^{ave}_{T} (GeV);#delta#sigma_{ex} / #delta#sigma_{stat}",par_->nPtBins(etaBin->etaBin()),&(par_->ptBinEdges(etaBin->etaBin()).front()));
	    hFrame->GetXaxis()->SetNoExponent();
	    hFrame->GetXaxis()->SetMoreLogLabels();
	    hFrame->SetLineWidth(lineWidth_);
	    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
	      hFrame->SetBinContent(bin,1.);
	    }

	    TH1* hErrExtrapolation = static_cast<TH1*>(hFrame->Clone("PlotMaker::plotExtrapolation::hErrExtrapolation"));
	    TH1* hErrFirstPoint = static_cast<TH1*>(hErrExtrapolation->Clone("PlotMaker::plotExtrapolation::hErrFirstPoint"));
	    hErrFirstPoint->SetLineColor(kBlue);
	    hErrFirstPoint->SetLineStyle(2);
	    hErrFirstPoint->SetLineWidth(3);
	    TH1* hErrLastPoint = static_cast<TH1*>(hErrExtrapolation->Clone("PlotMaker::plotExtrapolation::hErrLastPoint"));
	    hErrLastPoint->SetLineColor(kRed);
	    hErrLastPoint->SetLineStyle(3);
	    hErrLastPoint->SetLineWidth(3);
	    for(unsigned int i = 0; i < errExtrapolation.size(); ++i) {
	      hErrExtrapolation->SetBinContent(1+i,errExtrapolation.at(i));
	      hErrFirstPoint->SetBinContent(1+i,errFirstPoint.at(i));
	      hErrLastPoint->SetBinContent(1+i,errLastPoint.at(i));
	    }
	    hErrFirstPoint->Divide(hErrExtrapolation,hErrFirstPoint);
	    hErrLastPoint->Divide(hErrExtrapolation,hErrLastPoint);

	    TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin(),0.8);
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,0.8,label->GetSize());
	    leg->AddEntry(hErrLastPoint,"#delta#sigma_{ex} / #delta#sigma_{stat}(Last Point)","L");
	    leg->AddEntry(hErrFirstPoint,"#delta#sigma_{ex} / #delta#sigma_{stat}(First Point)","L");

	    out_->nextPad("Extrapolation uncertainties comparison "+etaBin->toString());
	    hFrame->GetYaxis()->SetRangeUser(0.,5.2);
	    hFrame->Draw();
	    hErrLastPoint->Draw("HIST][same");
	    hErrFirstPoint->Draw("HIS][Tsame");
	    label->Draw("same");
	    leg->Draw("same");
	    out_->logx();
	    gPad->RedrawAxis();
	    //out_->saveCurrentPad(histFileName("ExtrapolationUncertaintiesComparison",etaBin,sampleLabel,*rIt));

	    delete hFrame;
	    delete hErrExtrapolation;
	    delete hErrFirstPoint;
	    delete hErrLastPoint;
	    delete label;
	    delete leg;
	  } // End of loop over eta bins
	} // End if isComparedSample
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
 	    // Loop over to-be-compare Samples
 	    TString yAxisLabel = "";
 	    for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
 		sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
 	      SampleLabel sLabel1 = (*sCIt)->label1();
 	      SampleLabel sLabel2 = (*sCIt)->label2();
 	      const Sample* sample1 = ptBin->findSample(sLabel1);
 	      const Sample* sample2 = ptBin->findSample(sLabel2);
	     
 	      // Graphs of fitted resolutions
 	      std::vector<double> val;
 	      std::vector<double> uncert;
 	      sample1->valuesInExtrapolation(*rIt,val,uncert);
 	      std::vector<double> ptSoftx;
 	      sample1->ptSoft(ptSoftx);
 	      std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
 	      TGraphAsymmErrors* gVals1 = new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),&(ptSoftxUncert.front()),&(ptSoftxUncert.front()),&(uncert.front()),&(uncert.front()));
 	      setStyle(sample1,gVals1);

 	      val.clear();
 	      uncert.clear();
 	      sample2->valuesInExtrapolation(*rIt,val,uncert);
 	      ptSoftx.clear();
 	      sample2->ptSoft(ptSoftx);
 	      ptSoftxUncert.clear();
 	      ptSoftxUncert = std::vector<double>(ptSoftx.size(),0.);
 	      TGraphAsymmErrors* gVals2 = new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),&(ptSoftxUncert.front()),&(ptSoftxUncert.front()),&(uncert.front()),&(uncert.front()));
 	      setStyle(sample2,gVals2);

 	      // Deal with failed fits
 	      while( gVals1->GetY()[0] == gVals1->GetY()[1] ) {
 		gVals1->RemovePoint(0);
 		gVals2->RemovePoint(0);
 	      }
 	      while( gVals2->GetY()[0] == gVals2->GetY()[1] ) {
 		gVals1->RemovePoint(0);
 		gVals2->RemovePoint(0);
 	      }
 	      while( gVals1->GetEYhigh()[0] == 0. ) {
 		gVals1->RemovePoint(0);
 	      }
 	      while( gVals2->GetEYhigh()[0] == 0. ) {
 		gVals2->RemovePoint(0);
 	      }

 	      // Extrapolation functions
 	      TF1* fit1 = sample1->extrapolationFunction(*rIt,"Extrapolation_"+sLabel1+"_FitRange");
 	      fit1->SetLineColor(gVals1->GetMarkerColor());
 	      fit1->SetLineWidth(lineWidth_);
 	      TF1* extra1 = static_cast<TF1*>(fit1->Clone("Extrapolation_"+sLabel1+"_PlotRange"));
 	      extra1->SetRange(0.,1.4*ptSoftx.back());
 	      extra1->SetLineStyle(2);

 	      TF1* fit2 = sample2->extrapolationFunction(*rIt,"Extrapolation_"+sLabel2+"_FitRange");
 	      fit2->SetLineColor(gVals2->GetMarkerColor());
 	      fit2->SetLineWidth(lineWidth_);
 	      TF1* extra2 = static_cast<TF1*>(fit2->Clone("Extrapolation_"+sLabel2+"_PlotRange"));
 	      extra2->SetRange(0.,1.4*ptSoftx.back());
 	      extra2->SetLineStyle(2);

 	      TPaveText* label = labelMk_->etaPtAveBin(etaBin->etaBin(),ptBin->ptBin());
 	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,-1.,label->GetSize());
 	      leg->AddEntry(gVals1,labelMk_->sample(sLabel1),"PL");
 	      leg->AddEntry(gVals1,"#sigma("+labelMk_->ptSoftMax()+"#rightarrow 0) = "+util::toTString(fit1->GetParameter(0),1)+" #pm "+util::toTString(fit1->GetParError(0),1)+" GeV","");
 	      leg->AddEntry(gVals2,labelMk_->sample(sLabel2),"PL");
 	      leg->AddEntry(gVals2,"#sigma("+labelMk_->ptSoftMax()+"#rightarrow 0) = "+util::toTString(fit2->GetParameter(0),1)+" #pm "+util::toTString(fit2->GetParError(0),1)+" GeV","");

 	      TString yAxisLabel = sample1->labelQuantityInExtrapolation(*rIt);
	   
 	      // Create frame
 	      TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame",title_,
 				     1000,0.,1.3*par_->ptSoftMax(par_->nPtSoftBins()-1));
 	      util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"",yAxisLabel);
 	      hFrame->GetXaxis()->SetTitle("Threshold "+labelMk_->ptSoftMax());
	      double yMin1 = 0.;
	      double yMax1 = 0.;
	      util::HistOps::findYRangeInclErrors(gVals1,yMin1,yMax1);
	      double yMin2 = 0.;
	      double yMax2 = 0.;
	      util::HistOps::findYRangeInclErrors(gVals2,yMin2,yMax2);
	      double yMin = std::min(std::min(std::min(yMin1,yMin2),extra1->Eval(0.)),extra2->Eval(0.));
	      double yMax = std::max(yMax1,yMax2);
	      if( yMax > 100. ) yMax = 25.; // One corrupted fit point
	      double delta = yMax - yMin;
	      hFrame->GetYaxis()->SetRangeUser(std::max(0.,yMin-0.3*delta),yMax+1.5*delta);
 	      out_->nextMultiPad("Sample comparison: Extrapolation "+ptBin->toTString());
 	      hFrame->Draw();
 	      label->Draw("same");
 	      leg->Draw("same");
 	      gVals1->Draw("PE1same");
 	      extra1->Draw("same");
 	      fit1->Draw("same");
 	      gVals2->Draw("PE1same");
 	      extra2->Draw("same");
 	      fit2->Draw("same");
 	      out_->saveCurrentPad(histFileName("Extrapolation",ptBin,sLabel1,sLabel2,*rIt));
	      
 	      delete gVals1;
 	      delete extra1;
 	      delete fit1;
 	      delete gVals2;
 	      delete extra2;
 	      delete fit2;
 	      delete hFrame;
 	      delete label;
 	      delete leg;
 	    } // End of loop over Samples
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

		labels.push_back(labelMk_->type(*rIt));
	     
		// Graph of fitted resolutions
		std::vector<double> val;
		std::vector<double> uncert;
		sample->valuesInExtrapolation(*rIt,val,uncert);
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
	      TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame",title_,
				     1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	      util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"","Resolution");
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
	    util::HistOps::normHist(hPtGenAsym,"width");
 	    util::HistOps::setAxisTitles(hPtGenAsym,"p^{gen}_{T} Asymmetry","","events",true);
 	    hPtGenAsym->SetTitle(title_);
	    hPtGenAsym->GetXaxis()->SetNdivisions(505);
 	    hPtGenAsym->GetXaxis()->SetRangeUser(-0.24,0.24);
 	    hPtGenAsym->SetMarkerStyle(24);
 	    hPtGenAsym->SetMarkerColor(kBlack);
 	    hPtGenAsym->SetLineColor(hPtGenAsym->GetMarkerColor());
 	    hPtGenAsym->SetLineWidth(lineWidth_);
	      
 	    // Labels
 	    TPaveText* label = labelMk_->etaPtAvePtSoftGenBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
 	    // Draw
 	    util::HistOps::setYRange(hPtGenAsym,label->GetSize()-1);
 	    out_->nextMultiPad(sample->label()+": PLIPtGenAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
 	    hPtGenAsym->Draw("HIST");
 	    label->Draw("same");
 	    out_->saveCurrentPad(histFileName("PLI_PtGenAsym",ptBin,sample,ptSoftBinIdx));
	      
 	    delete hPtGenAsym;
 	    delete label;
 	  } // End of loop over pt bins
 	} // End of loop over eta bins
       } // End of loop over ptSoft bins


//       // +++++ PtGen spectra per bin ++++++++++++++++++++++++++++++++++++++++++++++++

//       // Loop over ptSoft bins
//       std::cout << "  PtGen spectra" << std::endl;
//       for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
// 	if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	
// 	// Loop over eta bins
// 	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
// 	  const EtaBin* etaBin = *etaBinIt;
	  
// 	  out_->newPage("PtGenSpectra");
	  
// 	  // Loop over pt bins
// 	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
// 	    const PtBin* ptBin = *ptBinIt;
// 	    if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	    
// 	    const Sample* sample = ptBin->mcTruthSample();
	      
// 	    // PtGen distribution
// 	    TH1* hPtGen = sample->histPtGen(ptSoftBinIdx);
// 	    util::HistOps::setAxisTitles(hPtGen,"p^{gen}_{T}","GeV","events");
// 	    hPtGen->SetTitle(title_);
// 	    hPtGen->SetMarkerStyle(24);
// 	    hPtGen->SetMarkerColor(kBlack);
// 	    hPtGen->SetLineColor(hPtGen->GetMarkerColor());
// 	    hPtGen->SetLineWidth(1);
	      
// 	    // Labels
// 	    TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
// 	    // Draw
// 	    util::HistOps::setYRange(hPtGen,label->GetSize()+1);
// 	    out_->nextMultiPad(sample->label()+": PLIPtGenSpectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
// 	    hPtGen->Draw("PE1");
// 	    label->Draw("same");
// 	    out_->saveCurrentPad(histFileName("PLI_PtGenSpectrum",ptBin,sample,ptSoftBinIdx));
	      
// 	    delete hPtGen;
// 	    delete label;
// 	  } // End of loop over pt bins
// 	} // End of loop over eta bins
//       } // End of loop over ptSoft bins



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
	  sample->valuesInExtrapolation(FitResult::PtGenAsym,val,uncert);
	  std::vector<double> ptSoftx;
	  sample->ptSoft(ptSoftx);
	  std::vector<double> ptSoftxUncert(ptSoftx.size(),0.);
	       
	  TGraphAsymmErrors* gVals = 
	    new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
				  &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
				  &(uncert.front()),&(uncert.front()));
	  gVals->SetMarkerStyle(24);
	  gVals->SetMarkerColor(kBlue);
	  gVals->SetLineColor(kBlue);
	       
	  // Extrapolation function
	  TF1* fit = sample->extrapolationFunction(FitResult::PtGenAsym,"ExtrapolationForPlotFitRange");
	  fit->SetLineColor(kBlue);
	  TF1* extra = static_cast<TF1*>(fit->Clone("ExtrapolationForPlotPlotRange"));
	  extra->SetRange(0.,1.4*ptSoftx.back());
	  extra->SetLineStyle(2);
	       
	  // Create frame
	  TH1* hFrame = new TH1D("PlotMaker::plotParticleLevelImbalance::hFrame",title_,
				 1000,0.,1.3*par_->ptSoftMax(par_->nPtSoftBins()-1));
	  util::HistOps::setAxisTitles(hFrame,"Threshold "+labelMk_->ptSoftGenMax(),"","#sqrt{2}#upoint#sigma_{A^{gen}}");
	  hFrame->GetYaxis()->SetRangeUser(0.,0.19);
	       
	  // Labels
	  TPaveText* label = labelMk_->etaPtAveBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-labelMk_->start(),label->GetSize());
	  leg->AddEntry(fit,"Fit  #sqrt{2}#upoint#sigma_{A^{gen}}("+labelMk_->ptSoftGenMax()+"#rightarrow 0)","PL");
	  leg->AddEntry(fit," = "+util::toTString(fit->GetParameter(0),5)+" #pm "+util::toTString(fit->GetParError(0),5),"");

	  // Linear scale
	  out_->nextMultiPad(sample->label()+": PLI PtGen Asymmetry Extrapolation "+ptBin->toTString());
	  hFrame->Draw();
	  gVals->Draw("PE1same");
	  extra->Draw("same");
	  fit->Draw("same");
	  label->Draw("same");
	  leg->Draw("same");
	  out_->saveCurrentPad(histFileName("PLI_PtGenAsymExtrapolation",ptBin,sample,FitResult::PtGenAsym));
	
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

      std::vector<TGraphAsymmErrors*> plis;
      std::vector<TF1*> pliFits;
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
	pli->SetRange(std::max(xMinPt_,0.75*gPLI->GetX()[0]),std::min(xMaxPt_,gPLI->GetX()[gPLI->GetN()-1]*1.5));
	pli->SetLineColor(kGreen);
	pli->SetLineStyle(2);
	
	// Ratio of measurement and fit
	TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gPLI,pli);
	
	// Create frames for main and ratio plots
	TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,0.109,"#sigma_{PLI}");
	hFrameMain->GetXaxis()->SetMoreLogLabels();
	hFrameMain->GetXaxis()->SetNoExponent();
	
	TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{gen}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	hFrameRatio->GetXaxis()->SetMoreLogLabels();
	hFrameRatio->GetXaxis()->SetNoExponent();
	
	// Labels
	TPaveText* label = labelMk_->etaBin(etaBin->mcTruthSampleLabel(),etaBin->etaBin());
	TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	leg->AddEntry(gPLI,"Particle-level imbalance","P");
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
	out_->saveCurrentPad(histFileName("PLI_Fit",etaBin,etaBin->mcTruthSampleLabel(),FitResult::PtGenAsym));

	plis.push_back(gPLI);
	pliFits.push_back(pli);
	
	delete gRatio;
	delete hFrameMain;
	delete hFrameRatio;
	delete label;
	delete leg;
      } // End of loop over eta bins

      // PLI for all eta bins in one plot
      out_->newPage("PLIFit");
      // Labels
      TPaveText* label = labelMk_->sampleBin(etaBins_.front()->mcTruthSampleLabel());
      TLegend* leg = util::LabelFactory::createLegendColWithOffset(par_->nEtaBins(),0.7*labelMk_->start(),label->GetSize());
      for(unsigned int i = 0; i < plis.size(); ++i) {
	leg->AddEntry(plis.at(i),labelMk_->etaRange(i),"PL");
      }

      TH1* hFrame = new TH1D("hFrame","",1000,xMinPt_,xMaxPt_);
      hFrame->GetYaxis()->SetRangeUser(1E-4,0.139);
      util::HistOps::setAxisTitles(hFrame,"p^{gen}_{T}","GeV","#sigma_{PLI}");
      hFrame->GetXaxis()->SetMoreLogLabels();
      hFrame->GetXaxis()->SetNoExponent();
      for(unsigned int i = 0; i < plis.size(); ++i) {
	plis.at(i)->SetMarkerStyle(24+i);
	plis.at(i)->SetMarkerColor(util::StyleSettings::color(i));
	plis.at(i)->SetLineColor(plis.at(i)->GetMarkerColor());
	plis.at(i)->SetLineWidth(2);
	pliFits.at(i)->SetLineColor(plis.at(i)->GetLineColor());
	pliFits.at(i)->SetLineWidth(2);
	pliFits.at(i)->SetLineStyle(1+i);
      }

      out_->nextPad(etaBins_.front()->mcTruthSampleLabel()+": PLI fits");
      hFrame->Draw();
      for(int i = static_cast<int>(plis.size())-1; i >= 0; --i) {
	pliFits.at(i)->Draw("same");
	plis.at(i)->Draw("PE1same");
      }
      label->Draw("same");
      leg->Draw("same");
      out_->logx();
      out_->saveCurrentPad(histFileName("PLI_Fits",etaBins_.front()->mcTruthSampleLabel()));

      delete leg;
      delete label;
      delete hFrame;
      for(unsigned int i = 0; i < plis.size(); ++i) {
	delete plis.at(i);
	delete pliFits.at(i);
      }      
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
	  
	  // Plot spectra only for MC sample that is
	  // compared to the data
	  bool isToBeCompared = false;
	  for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	      sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	    if( (*sCIt)->contains(sample->label()) ) {
	      isToBeCompared = true;
	      break;
	    }
	  }

	  for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < sample->nPtSoftBins(); ++ptSoftBinIdx) {
	    if( par_->verbosity() > 1 ) std::cout << "      PtSoftBin " << ptSoftBinIdx << std::endl;
	    
	    // Get ptGen spectrum and fit from sample and tweak style
	    TH1* hPtGen = sample->histPtGen(ptSoftBinIdx);
	    util::HistOps::setAxisTitles(hPtGen,"p^{gen}_{T}","GeV","jets",true);
	    setStyle(sample,hPtGen);

	    TH1* hPdf = sample->histPdfPtTrue(ptSoftBinIdx);
	    hPdf->SetLineColor(kRed);
	    hPdf->SetLineWidth(lineWidth_);
	    
	    // Labels
	    TPaveText* label = labelMk_->etaPtAvePtSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,0.8*labelMk_->start(),label->GetSize());
	    leg->AddEntry(hPtGen,labelMk_->sample(sample->label()),"P");
	    leg->AddEntry(hPdf,"Estimate","L");
	    
	    // Linear scale
	    util::HistOps::setYRange(hPtGen,label->GetSize()+leg->GetNRows());
	    out_->nextMultiPad(sample->label()+": PtGen Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGen->Draw("PE1");
	    hPdf->Draw("Lsame");
	    label->Draw("same");
	    leg->Draw("same");
	    gPad->RedrawAxis();
	    out_->saveCurrentPad(histFileName("PtGen",ptBin,sample,ptSoftBinIdx));

	    // Illustration of migration effects
	    TLine* ll = new TLine(ptBin->min(),0.,ptBin->min(),hPtGen->GetBinContent(hPtGen->GetMaximumBin()));
	    ll->SetLineWidth(4);
	    ll->SetLineStyle(3);
	    ll->SetLineColor(kOrange+1);
	    TLine* lr = new TLine(ptBin->max(),0.,ptBin->max(),hPtGen->GetBinContent(hPtGen->GetMaximumBin()));
	    lr->SetLineWidth(4);
	    lr->SetLineStyle(3);
	    lr->SetLineColor(kOrange+1);
	    TLegend* legl = util::LabelFactory::createLegendWithOffset(1,label->GetSize());
	    legl->AddEntry(ll,"p^{ave}_{T} Interval Boundaries","L");
	    util::HistOps::setYRange(hPtGen,label->GetSize()+legl->GetNRows());
	    out_->nextMultiPad(sample->label()+": PtGen Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGen->Draw("PE1");
	    ll->Draw("same");
	    lr->Draw("same");
	    label->Draw("same");
	    legl->Draw("same");
	    gPad->RedrawAxis();
	    out_->saveCurrentPad(histFileName("PtGen-IllustrationOfMigration",ptBin,sample,ptSoftBinIdx));
	    
	    // Log scale
	    util::HistOps::setYRange(hPtGen,label->GetSize()+leg->GetNRows(),3E-8);
	    out_->nextMultiPad(sample->label()+": PtGen Spectrum (log) "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGen->Draw("PE1");
	    hPdf->Draw("Lsame");
	    label->Draw("same");
	    leg->Draw("same");
	    gPad->RedrawAxis();
	    out_->logy();
	    out_->saveCurrentPad(histFileName("PtGenLog",ptBin,sample,ptSoftBinIdx));
	    
	    delete ll;
	    delete lr;
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

//     // Loop over SampleLabels
//     for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
// 	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      
//       SampleLabel sampleLabel = sTIt->first;
//       if( par_->verbosity() > 1 ) std::cout << "  " << sampleLabel << std::endl;
      
//       // Loop over distributions
//       for(unsigned int distIdx = 0; distIdx < 4; ++distIdx) {

// 	// Loop over ptSoft bins
// 	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
// 	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	  
// 	  // Loop over eta bins
// 	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
// 	    const EtaBin* etaBin = *etaBinIt;
	    
// 	    if( distIdx == 0 ) out_->newPage("PtAveSpectrum");
// 	    else if( distIdx == 1 ) out_->newPage("Pt1Spectrum");
// 	    else if( distIdx == 2 ) out_->newPage("Pt2Spectrum");
// 	    else if( distIdx == 3 ) out_->newPage("Pt3Spectrum");
	    
// 	    // Loop over pt bins
// 	    for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
// 	      const PtBin* ptBin = *ptBinIt;
// 	      if( par_->verbosity() > 1 ) std::cout << "      " << ptBin->toTString() << std::endl;
	    
// 	      // Loop over Samples and select current one (by sampleLabel)
// 	      for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
// 		const Sample* sample = sIt->second;
// 		if( sample->label() != sampleLabel ) continue;
	      
// 		// Spectrum
// 		TH1* h = 0;
// 		int minBin = 1;
// 		int maxBin = 1;
// 		if( distIdx == 0 ) {
// 		  h = sample->histPtAve(ptSoftBinIdx);
// 		  minBin = std::max(h->FindBin(0.8*ptBin->min()),1);
// 		  maxBin = std::min(h->FindBin(1.2*ptBin->max()),h->GetNbinsX());
// 		  util::HistOps::setAxisTitles(h,"p^{ave}_{T}","GeV","events");
// 		} else if( distIdx == 1 ) {
// 		  h = sample->histPt1(ptSoftBinIdx);
// 		  minBin = std::max(h->FindBin(0.5*ptBin->min()),1);
// 		  maxBin = std::min(h->FindBin(1.5*ptBin->max()),h->GetNbinsX());
// 		  util::HistOps::setAxisTitles(h,"p_{T,1}","GeV","events");
// 		} else if( distIdx == 2 ) {
// 		  h = sample->histPt2(ptSoftBinIdx);
// 		  minBin = std::max(h->FindBin(0.5*ptBin->min()),1);
// 		  maxBin = std::min(h->FindBin(1.5*ptBin->max()),h->GetNbinsX());
// 		  util::HistOps::setAxisTitles(h,"p_{T,2}","GeV","events");
// 		} else if( distIdx == 3 ) {
// 		  h = sample->histPt3(ptSoftBinIdx);
// 		  minBin = 1;
// 		  maxBin = std::min(h->FindBin(0.3*ptBin->min()),h->GetNbinsX());
// 		  util::HistOps::setAxisTitles(h,"p_{T,3}","GeV","events");
// 		}
// 		h->GetXaxis()->SetRange(minBin,maxBin);
// 		h->Rebin(10);
// 		setStyle(sample,h);
	      
// 		// Labels
// 		TPaveText* label = labelMk_->etaPtAvePtSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
// 		util::HistOps::setYRange(h,label->GetSize()+1);
// 		out_->nextMultiPad(sample->label()+": Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
// 		h->Draw("HIST");//h->Draw("PE1");
// 		label->Draw("same");
// 		if( distIdx == 0 ) out_->saveCurrentPad(histFileName("PtAve",ptBin,sample,ptSoftBinIdx));
// 		else if( distIdx == 1 ) out_->saveCurrentPad(histFileName("PtJet1",ptBin,sample,ptSoftBinIdx));
// 		else if( distIdx == 2 ) out_->saveCurrentPad(histFileName("PtJet2",ptBin,sample,ptSoftBinIdx));
// 		else if( distIdx == 3 ) out_->saveCurrentPad(histFileName("PtJet3",ptBin,sample,ptSoftBinIdx));

// 		delete h;
// 		delete label;
// 	      } // End of loop over Samples
// 	    } // End of loop over pt bins
// 	  } // End of loop over eta bins
// 	} // End of loop over ptSoft bins
//       } // End of loop over SampleLabels
//     } // End of loop over spectra


    if( etaBins_.front()->nComparedSamples() > 0 ) {
      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); 
	  rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
	if( par_->verbosity() > 1 ) std::cout << "  " << FitResult::toString(*rIt) << std::endl;
	
	// Loop over ptSoft bins
	for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	  if( par_->verbosity() > 1 ) std::cout << "    PtSoftBin " << ptSoftBinIdx << std::endl;
	  
	  // Loop over eta and bins
	  std::vector<TH1*> hPtSpecInclS1;
	  std::vector<TH1*> hPtSpecInclS2;
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;
	    out_->newPage("Pt Spectra");
	    
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

		// Spectra
		std::vector<TH1*> hPtSpecS1;
		std::vector<TH1*> hPtSpecS2;		
	 	hPtSpecS1.push_back(sample1->histPtAve(ptSoftBinIdx));
		hPtSpecS1.push_back(sample1->histPt1(ptSoftBinIdx));
		hPtSpecS1.push_back(sample1->histPt2(ptSoftBinIdx));
		hPtSpecS1.push_back(sample1->histPt3(ptSoftBinIdx));
		hPtSpecS2.push_back(sample2->histPtAve(ptSoftBinIdx));
		hPtSpecS2.push_back(sample2->histPt1(ptSoftBinIdx));
		hPtSpecS2.push_back(sample2->histPt2(ptSoftBinIdx));
		hPtSpecS2.push_back(sample2->histPt3(ptSoftBinIdx));

		// Norm sample2 to sample1 from ptAve spectrum
		double norm = 1.;
		if( hPtSpecS2.front()->Integral() > 0. ) {
		  norm = hPtSpecS1.front()->Integral()/hPtSpecS2.front()->Integral();
		}
		for(unsigned int i = 0; i < hPtSpecS2.size(); ++i) {
		  hPtSpecS2.at(i)->Scale(norm);
		}

		// Set style
		for(unsigned int i = 0; i < hPtSpecS1.size(); ++i) {
		  hPtSpecS1.at(i)->SetMarkerStyle(20);
		  hPtSpecS1.at(i)->SetMarkerColor(kBlack);
		  hPtSpecS1.at(i)->SetLineColor(kBlack);
		  hPtSpecS1.at(i)->SetLineWidth(2);

		  hPtSpecS2.at(i)->SetMarkerStyle(1);
		  hPtSpecS2.at(i)->SetLineColor(kBlack);
		  hPtSpecS2.at(i)->SetLineWidth(2);
		  hPtSpecS2.at(i)->SetFillColor(38);
		}
		
		// Add to inclusive spectra
		for(unsigned int i = 0; i < hPtSpecS1.size(); ++i) {
		  if( ptBinIt == etaBin->ptBinsBegin() ) {
		    // Clone hists to be able to have independent binning etc
		    hPtSpecInclS1.push_back(static_cast<TH1*>(hPtSpecS1.at(i)->Clone("hPtSpec1Incl"+util::toTString(i))));
		    hPtSpecInclS2.push_back(static_cast<TH1*>(hPtSpecS2.at(i)->Clone("hPtSpec2Incl"+util::toTString(i))));
		  } else {
		    hPtSpecInclS1.at(i)->Add(hPtSpecS1.at(i));
		    hPtSpecInclS2.at(i)->Add(hPtSpecS2.at(i));
		  }
		}
		
		// tweak x range
		for(unsigned int i = 0; i < hPtSpecS2.size(); ++i) {
		  int rebin = 2;
		  if( i == 1 || i == 2 ) rebin = 5;
		  hPtSpecS1.at(i)->Rebin(rebin);
		  hPtSpecS2.at(i)->Rebin(rebin);
		  
		  TH1* h = hPtSpecS2.at(i);
		  int minBin = std::max(h->FindBin(0.95*ptBin->min()),1);
		  int maxBin = std::min(h->FindBin(1.05*ptBin->max()),h->GetNbinsX());
		  if( i == 1 || i == 2 ) {
		    minBin = std::max(h->FindBin(0.5*ptBin->min()),1);
		    maxBin = std::min(h->FindBin(1.5*ptBin->max()),h->GetNbinsX());
		  } else if( i == 3 ) {
		    minBin = 1;
		    maxBin = std::min(h->FindBin(0.3*ptBin->min()),h->GetNbinsX());
		  }
		  hPtSpecS1.at(i)->GetXaxis()->SetRange(minBin,maxBin);
		  hPtSpecS2.at(i)->GetXaxis()->SetRange(minBin,maxBin);

		  if( i == 0 ) {
		    util::HistOps::setAxisTitles(hPtSpecS1.at(i),"p^{ave}_{T}","GeV","events");		 
		    util::HistOps::setAxisTitles(hPtSpecS2.at(i),"p^{ave}_{T}","GeV","events");		 
		  } else {
		    util::HistOps::setAxisTitles(hPtSpecS1.at(i),"p_{T,"+util::toTString(i)+"}","GeV","events");		 
		    util::HistOps::setAxisTitles(hPtSpecS2.at(i),"p_{T,"+util::toTString(i)+"}","GeV","events");		 
		  }
		}
		  
		// Labels
		TPaveText* label = labelMk_->etaPtAvePtSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);	    
		TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.65*labelMk_->start(),label->GetSize());
		leg->AddEntry(hPtSpecS1.front(),labelMk_->sample(sLabel1),"P");
		leg->AddEntry(hPtSpecS2.front(),labelMk_->sample(sLabel2),"F");


		// Plot
		for(unsigned int i = 0; i < hPtSpecS2.size(); ++i) {
		  if( i > 0 ) 
		    util::HistOps::setYRange(hPtSpecS2.at(i),label->GetSize()+leg->GetNRows(),3E-1);
		  else
		    util::HistOps::setYRange(hPtSpecS2.at(i),label->GetSize()+leg->GetNRows());

		  out_->nextMultiPad(sample1->label()+" vs "+sample2->label()+": Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		  hPtSpecS2.at(i)->Draw("HIST");
		  hPtSpecS1.at(i)->Draw("PE1same");
		  label->Draw("same");
		  leg->Draw("same");
		  if( i == 0 ) {
		    out_->saveCurrentPad(histFileName("PtAve",ptBin,sLabel1,sLabel2,*rIt,ptSoftBinIdx));
		  } else {
		    out_->logy();
		    if( i == 1 ) out_->saveCurrentPad(histFileName("PtJet1",ptBin,sLabel1,sLabel2,*rIt,ptSoftBinIdx));
		    else if( i == 2 ) out_->saveCurrentPad(histFileName("PtJet2",ptBin,sLabel1,sLabel2,*rIt,ptSoftBinIdx));
		    else if( i == 3 ) out_->saveCurrentPad(histFileName("PtJet3",ptBin,sLabel1,sLabel2,*rIt,ptSoftBinIdx));
		  }
		}

		for(unsigned int i = 0; i < hPtSpecS2.size(); ++i) {
		  delete hPtSpecS1.at(i);
		  delete hPtSpecS2.at(i);
		}
		delete label;
		delete leg;
	      } // End of loop over to-be-compared samples
	    } // End of loop over pt bins

	    // Plot inclusive spectra
	    out_->newPage("Pt Spectra");

	    // Currently only implemented if two samples are to be compared
	    if( etaBin->nComparedSamples() == 1 ) {
	      ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
	      SampleLabel sLabel1 = (*sCIt)->label1();
	      SampleLabel sLabel2 = (*sCIt)->label2();

	      // tweak x range
	      for(unsigned int i = 0; i < hPtSpecInclS2.size(); ++i) {
		int minBin = 1;
		int maxBin = 1;
		if( i == 3 ) {
		  hPtSpecInclS1.at(i)->Rebin(8);
		  hPtSpecInclS2.at(i)->Rebin(8);
		  maxBin = std::min(hPtSpecInclS2.at(i)->FindBin(245),hPtSpecInclS2.at(i)->GetNbinsX());
		} else {
		  hPtSpecInclS1.at(i)->Rebin(25);
		  hPtSpecInclS2.at(i)->Rebin(25);
		  maxBin = std::min(hPtSpecInclS2.at(i)->FindBin(1690),hPtSpecInclS2.at(i)->GetNbinsX());
		}
		hPtSpecInclS1.at(i)->GetXaxis()->SetNdivisions(505);
		hPtSpecInclS2.at(i)->GetXaxis()->SetNdivisions(505);
		hPtSpecInclS1.at(i)->GetXaxis()->SetRange(minBin,maxBin);
		hPtSpecInclS2.at(i)->GetXaxis()->SetRange(minBin,maxBin);
		if( i == 0 ) {
		  util::HistOps::setAxisTitles(hPtSpecInclS1.at(i),"p^{ave}_{T}","GeV","events");		 
		  util::HistOps::setAxisTitles(hPtSpecInclS2.at(i),"p^{ave}_{T}","GeV","events");		 
		} else {
		  util::HistOps::setAxisTitles(hPtSpecInclS1.at(i),"p_{T,"+util::toTString(i)+"}","GeV","events");		 
		  util::HistOps::setAxisTitles(hPtSpecInclS2.at(i),"p_{T,"+util::toTString(i)+"}","GeV","events");		 
		}
	      }

	      // Labels
	      TPaveText* label = labelMk_->etaPtSoftBin(etaBin->etaBin(),ptSoftBinIdx);
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-0.65*labelMk_->start(),label->GetSize());
	      leg->AddEntry(hPtSpecInclS1.front(),labelMk_->sample(sLabel1),"P");
	      leg->AddEntry(hPtSpecInclS2.front(),labelMk_->sample(sLabel2),"F");
	      
	      // Plot
	      for(unsigned int i = 0; i < hPtSpecInclS2.size(); ++i) {
		util::HistOps::setYRange(hPtSpecInclS2.at(i),label->GetSize()+leg->GetNRows(),3E-1);
		  
		out_->nextMultiPad(sLabel1+" vs "+sLabel2+": Inclusive Spectrum "+etaBin->toString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		hPtSpecInclS2.at(i)->Draw("HIST");
		hPtSpecInclS1.at(i)->Draw("PE1same");
		label->Draw("same");
		leg->Draw("same");
		out_->logy();
		if( i == 0 )
		  out_->saveCurrentPad(histFileName("PtAve",sLabel1,sLabel2,etaBin,ptSoftBinIdx));
		else if( i == 1 ) out_->saveCurrentPad(histFileName("PtJet1",sLabel1,sLabel2,etaBin,ptSoftBinIdx));
		else if( i == 2 ) out_->saveCurrentPad(histFileName("PtJet2",sLabel1,sLabel2,etaBin,ptSoftBinIdx));
		else if( i == 3 ) out_->saveCurrentPad(histFileName("PtJet3",sLabel1,sLabel2,etaBin,ptSoftBinIdx));
	      }

	      for(unsigned int i = 0; i < hPtSpecInclS2.size(); ++i) {
		delete hPtSpecInclS1.at(i);
		delete hPtSpecInclS2.at(i);
	      }
	      delete label;
	      delete leg;
	    }
	  } // End of loop over eta bins
	} // End of loop over ptSoft bins
      } // End of loop over FitResult types
    } // End if samples to be compared


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

	  if( par_->verbosity() > 1 ) std::cout << "      " << sampleLabel << std::endl;

	  TGraphAsymmErrors* gRes = etaBin->extrapolatedResolution(sampleLabel,fitResType);
	  setStyle(sampleLabel,gRes);
	  gRes->SetMarkerStyle(27);

	  TGraphAsymmErrors* gResCorr = etaBin->correctedResolution(sampleLabel,fitResType);
	  setStyle(sampleLabel,gResCorr);
	  //gResCorr->SetMarkerStyle(21);

	  // MC truth resolution function
	  TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBinIdx));
	  mcTruth->SetLineColor(kRed);
	  mcTruth->SetLineWidth(lineWidth_);

	  // Ratio of measurement and mc truth
	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gResCorr,mcTruth);
	  gRatio->SetLineWidth(lineWidth_);

	  // PLI function
	  TF1* pli = etaBin->pliFunc("PLI_Eta"+util::toTString(etaBinIdx));
	  pli->SetRange(xMinPt_,xMaxPt_);
	  pli->SetLineColor(kGreen);
	  pli->SetLineStyle(2);
	  pli->SetLineWidth(lineWidth_);

	  // Create frames for main and ratio plots
	  TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{true}_{T}");
	  hFrameMain->SetTitle(title_);
	  hFrameMain->GetXaxis()->SetMoreLogLabels();
	  hFrameMain->GetXaxis()->SetNoExponent();

	  double yMinRatio = yMinResRatio_;
	  double yMaxRatio = yMaxResRatio_;
	  if( sTIt->second == Sample::Data && etaBin->etaBin() == 4 ) {
	    yMinRatio = 0.91;
	    yMaxRatio = 1.68;
	  }
	  TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{true}_{T}","GeV","#frac{Meas.}{Truth}",yMinRatio,yMaxRatio);
	  hFrameRatio->GetXaxis()->SetMoreLogLabels();
	  hFrameRatio->GetXaxis()->SetNoExponent();
	  hFrameRatio->SetLineWidth(lineWidth_);
	  hFrameRatio->SetLineColor(mcTruth->GetLineColor());
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,0.8*labelMk_->start(),label->GetSize());
	  // 	  leg->AddEntry(gRes,"Resolution (Extrapolated)","P");
	  // 	  leg->AddEntry(pli,"Particle Level Imbalance (PLI)","L");
	  // 	  leg->AddEntry(gResCorr,"Resolution (PLI Corrected)","P");
	  // 	  leg->AddEntry(mcTruth,"Resolution (MC Truth)","L");

	  leg->AddEntry(gResCorr,"Measurement","P");
	  leg->AddEntry(mcTruth,"MC Truth","L");

	  out_->nextMainPad(sampleLabel+": Resolution "+etaBin->toString());
	  hFrameMain->Draw();
	  //pli->Draw("same");
	  //gRes->Draw("PE1same");
	  mcTruth->Draw("same");
	  gResCorr->Draw("PE1same");
	  label->Draw("same");
	  leg->Draw("same");
	  out_->nextRatioPad();
	  hFrameRatio->Draw("][");
	  //fitToRatio->Draw("same");
	  gRatio->Draw("PE1same");
	  out_->logx();
	  out_->saveCurrentPad(histFileName("Resolution",etaBin,sampleLabel,fitResType));

// 	  if( sTIt->second == Sample::Data ) outDirData->WriteTObject(gResCorr);
// 	  else if( sTIt->second == Sample::MC ) outDirMC->WriteTObject(gResCorr);

	  delete gRes;
	  delete gResCorr;
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
//     fileData->Close();
//     fileMC->Close();
//     delete fileData;
//     delete fileMC;


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
	  setStyle(sLabel1,gRes1);
	  TGraphAsymmErrors* gRes2 = etaBin->correctedResolution(sLabel2,fitResType);
	  setStyle(sLabel2,gRes2);
	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gRes1,gRes2);
	  gRatio->SetLineWidth(lineWidth_);

	  TGraphAsymmErrors* gSystErr1 = etaBin->correctedResolutionSystUncert(sLabel1,sLabel2,fitResType);
	  gSystErr1->SetFillColor(kYellow);
	  gSystErr1->SetMarkerStyle(1);
	  gSystErr1->SetMarkerColor(gSystErr1->GetFillColor());
	  gSystErr1->SetLineColor(gSystErr1->GetFillColor());
	  gSystErr1->SetLineWidth(1);
	  gSystErr1->SetFillStyle(1001);

	  TGraphAsymmErrors* gRatioSystErr = static_cast<TGraphAsymmErrors*>(gRatio->Clone());
	  for(int i = 0; i < gRatioSystErr->GetN(); ++i) {
	    double res2 = gRes2->GetY()[i];
	    gRatioSystErr->GetEYlow()[i] = gSystErr1->GetEYlow()[i]/res2;
	    gRatioSystErr->GetEYhigh()[i] = gSystErr1->GetEYhigh()[i]/res2;
	  }
	  gRatioSystErr->SetFillColor(gSystErr1->GetFillColor());
	  gRatioSystErr->SetMarkerStyle(1);
	  gRatioSystErr->SetMarkerColor(gRatioSystErr->GetFillColor());
	  gRatioSystErr->SetLineColor(gRatioSystErr->GetFillColor());
	  gRatioSystErr->SetLineWidth(1);
	  gRatioSystErr->SetFillStyle(gSystErr1->GetFillStyle());

	  // Hack to remove first point in case of failed fits
	  if( gRes1->GetY()[0] < 0.02 ) {
	    gRes1->RemovePoint(0);
	    gRatio->RemovePoint(0);
	    gSystErr1->RemovePoint(0);
	    gRatioSystErr->RemovePoint(0);
	  }

	  // If data / MC comparison for scaling factor is shown,
	  // plot also MC truth
	  // MC truth resolution function
	  bool showMCTruth = false;//etaBin->hasKValue(sLabel1,sLabel2,fitResType);
	  TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBinIdx));
	  mcTruth->SetLineColor(kRed);
	  mcTruth->SetLineWidth(lineWidth_);

	  // Create frames for main and ratio plots
	  TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{true}_{T}");
	  hFrameMain->SetTitle(title_);
	  hFrameMain->GetXaxis()->SetMoreLogLabels();
	  hFrameMain->GetXaxis()->SetNoExponent();

	  TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{true}_{T}","GeV","#frac{Data}{Sim.}",yMinResRatio_,yMaxResRatio_);
	  hFrameRatio->GetXaxis()->SetMoreLogLabels();
	  hFrameRatio->GetXaxis()->SetNoExponent();
	  hFrameRatio->SetLineWidth(lineWidth_);
	  hFrameRatio->SetLineColor(gRes2->GetLineColor());
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	  TLegend* leg = 0;
	  if( showMCTruth ) 
	    leg = util::LabelFactory::createLegendCol(5,0.75*(labelMk_->start()));
	  else 
	    leg = util::LabelFactory::createLegendCol(4,0.75*(labelMk_->start()));
	  leg->AddEntry(gRes1,labelMk_->sample(sLabel1),"P");
	  leg->AddEntry(gRes2,labelMk_->sample(sLabel2),"P");
	  leg->AddEntry(gRes1,"#delta#sigma_{ex}","L");
	  leg->AddEntry(gSystErr1,"Syst. Uncert.","F");
	  if( showMCTruth ) 
	    leg->AddEntry(mcTruth,"MC truth","L");

	  out_->nextMainPad(sLabel1+"vs"+sLabel2+": Resolution "+etaBin->toString());
	  hFrameMain->Draw();
	  if( showMCTruth) mcTruth->Draw("same");
	  gSystErr1->Draw("PE3same");
	  gRes1->Draw("PE1same");
	  gRes2->Draw("PE1same");
	  label->Draw("same");
	  leg->Draw("same");
	  out_->nextRatioPad();
	  hFrameRatio->Draw("][");
	  gRatioSystErr->Draw("PE3same");
	  hFrameRatio->Draw("][same");
	  gRatio->Draw("PE1same");
	  out_->logx();
	  gPad->RedrawAxis();
	  out_->saveCurrentPad(histFileName("Resolution",etaBin,sLabel1,sLabel2,fitResType));

	  std::cout << "\nEta: " << etaBin->etaBin() << std::endl;
	  for(int i = 0; i < gRes1->GetN(); ++i) {
	    std::cout << i << ": " << gRes1->GetX()[i] << ": " << 100.*gRes1->GetY()[i] << " +/- " << 100.*gRes1->GetEYhigh()[i] << " +" << 100.*gSystErr1->GetEYhigh()[i] << " -" << 100.*gSystErr1->GetEYlow()[i] << std::endl;
	    std::cout << i << ": " << gRes2->GetX()[i] << ": " << 100.*gRes2->GetY()[i] << " +/- " << 100.*gRes2->GetEYhigh()[i] << std::endl;
	  }

	  delete gRes1;
	  delete gRes2;
	  delete gSystErr1;
	  delete gRatio;
	  delete gRatioSystErr;
	  delete mcTruth;
	  delete hFrameMain;
	  delete hFrameRatio;
	  delete label;
	  delete leg;
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins


    // +++++ Multi Canvas ++++++++++++++++++++++++++++++++++++++++++++++++

    if( par_->verbosity() > 1 ) std::cout << "  MultiCanvas plots" << std::endl;
    TH1* hFrame = new TH1D("hFrame",";p^{true}_{T} (GeV);#sigma / p^{true}_{T}  ",10000,xMinPt_,xMaxPt_);
    hFrame->GetYaxis()->SetRangeUser(yMinExtraRes_,yMaxExtraRes_);

    // Loop over SampleLabels
    if( par_->verbosity() > 1 ) std::cout << "    Resolution per sample" << std::endl;
    for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
 	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      SampleLabel sampleLabel = sTIt->first;
      if( par_->verbosity() > 1 ) std::cout << "      " << sampleLabel << std::endl;	

      util::MultiCanvas* mc = new util::MultiCanvas("Resolution",3,2,5,true);
      TLegend* leg = util::LabelFactory::createLegendWithOffset(2,0.25+util::LabelFactory::lineHeightMultiCan(),util::LabelFactory::lineHeightMultiCan());
      mc->adjustLegend(leg);

      // Loop over eta bins
      for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	const EtaBin* etaBin = *etaBinIt;
	if( par_->verbosity() > 1 ) std::cout << "        Eta " << etaBin->etaBin() << "   " << std::flush;	

	// Only first FitResultType for the time being
	FitResult::Type fitResType = *(etaBin->fitResultTypesBegin());
	
	TGraphAsymmErrors* gResCorr = etaBin->correctedResolution(sampleLabel,fitResType);
	setStyle(sampleLabel,gResCorr);
	mc->markForDeletion(gResCorr);
	
	// MC truth resolution function
	TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBin->etaBin()));
	mcTruth->SetLineColor(kRed);
	mcTruth->SetLineWidth(lineWidth_);
	mc->markForDeletion(mcTruth);

	// Ratio of measurement and mc truth
	TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gResCorr,mcTruth);
	gRatio->SetLineWidth(lineWidth_);
	mc->markForDeletion(gRatio);
	
	// Create frames for main and ratio plots
	TH1* mFrame = mc->mainFrame(hFrame,etaBin->etaBin());
	TH1* rFrame = mc->ratioFrame(hFrame,"#frac{Meas.}{Truth}",yMinResRatio_,yMaxResRatio_,etaBin->etaBin());
	rFrame->SetLineWidth(lineWidth_);
	rFrame->SetLineColor(mcTruth->GetLineColor());

	// Labels
	TPaveText* label = labelMk_->etaBin(etaBin->etaBin(),0.45);
	mc->markForDeletion(label);
	mc->adjustPaveText(label);
	mc->canvas()->cd();
	TPad* padt = mc->mainPad(etaBin->etaBin());
	padt->Draw();
	padt->cd();
	mFrame->Draw();
	mcTruth->Draw("same");
	gResCorr->Draw("PE1same");
	label->Draw("same");
	TPad* padr = mc->ratioPad(etaBin->etaBin());
	padr->Draw();
	padr->cd();
	rFrame->Draw();
	gRatio->Draw("PE1same");
	gPad->RedrawAxis();
 	if( etaBin->etaBin() == 0 ) {
 	  leg->AddEntry(gResCorr,"Measurement","P");
 	  leg->AddEntry(mcTruth,"MC Truth","L");
 	}
	if( par_->verbosity() > 1 ) std::cout << "ok" << std::endl;
      } // End of loop over eta bins

      TPaveText* label = util::LabelFactory::createPaveTextWithOffset(1,1.,0.25,util::LabelFactory::lineHeightMultiCan());
      label->AddText(labelMk_->sample(sampleLabel));
      mc->markForDeletion(label);
      mc->adjustPaveText(label);
      
      mc->canvas()->cd();
      TPad* padt = mc->mainPad(5);
      padt->Draw();
      padt->cd();
      label->Draw();
      leg->Draw("same");      
      mc->moreLogLabelsX();
      mc->noExponentX();
      mc->setLogx();
      out_->saveCanvas(mc->canvas(),histFileName("Resolution",sampleLabel));
      delete mc;
      delete leg;
    } // End of loop over SampleTypes


    if( par_->verbosity() > 1 ) std::cout << "  Sample comparison plots" << std::endl;
    // Only first FitResultType for the time being
    FitResult::Type fitResType = *(etaBins_.front()->fitResultTypesBegin());

    // Loop over to-be-compare Samples
    for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
      SampleLabel sLabel1 = (*sCIt)->label1();
      SampleLabel sLabel2 = (*sCIt)->label2();
      if( par_->verbosity() > 1 ) std::cout << "    " << sLabel1 << " vs " << sLabel2 << ":" <<  std::endl;

      util::MultiCanvas* mc = new util::MultiCanvas("Resolution",3,2,5,true);
      TLegend* leg = util::LabelFactory::createLegendWithOffset(4,0.2,util::LabelFactory::lineHeightMultiCan());
      mc->adjustLegend(leg);
      
      // Loop over eta bins
      for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	const EtaBin* etaBin = *etaBinIt;
	if( par_->verbosity() > 1 ) std::cout << "      Eta " << etaBin->etaBin() << std::endl;

	// Create graphs of extrapolated resolution corrected for PLI
	TGraphAsymmErrors* gRes1 = etaBin->correctedResolution(sLabel1,fitResType);
	setStyle(sLabel1,gRes1);
	mc->markForDeletion(gRes1);
	TGraphAsymmErrors* gRes2 = etaBin->correctedResolution(sLabel2,fitResType);
	setStyle(sLabel2,gRes2);
	mc->markForDeletion(gRes2);
	TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gRes1,gRes2);
	gRatio->SetLineWidth(lineWidth_);
	mc->markForDeletion(gRatio);

	TGraphAsymmErrors* gSystErr1 = etaBin->correctedResolutionSystUncert(sLabel1,sLabel2,fitResType);
	gSystErr1->SetFillColor(kYellow);
	gSystErr1->SetMarkerStyle(1);
	gSystErr1->SetMarkerColor(gSystErr1->GetFillColor());
	gSystErr1->SetLineColor(gSystErr1->GetFillColor());
	gSystErr1->SetLineWidth(1);
	gSystErr1->SetFillStyle(1001);
	mc->markForDeletion(gSystErr1);
	
	TGraphAsymmErrors* gRatioSystErr = static_cast<TGraphAsymmErrors*>(gRatio->Clone());
	for(int i = 0; i < gRatioSystErr->GetN(); ++i) {
	  double res2 = gRes2->GetY()[i];
	  gRatioSystErr->GetEYlow()[i] = gSystErr1->GetEYlow()[i]/res2;
	  gRatioSystErr->GetEYhigh()[i] = gSystErr1->GetEYhigh()[i]/res2;
	}
	gRatioSystErr->SetFillColor(gSystErr1->GetFillColor());
	gRatioSystErr->SetMarkerStyle(1);
	gRatioSystErr->SetMarkerColor(gRatioSystErr->GetFillColor());
	gRatioSystErr->SetLineColor(gRatioSystErr->GetFillColor());
	gRatioSystErr->SetLineWidth(1);
	gRatioSystErr->SetFillStyle(gSystErr1->GetFillStyle());
	mc->markForDeletion(gRatioSystErr);

	// Hack to remove first point in case of failed fits
	if( gRes1->GetY()[0] < 0.02 ) {
	  gRes1->RemovePoint(0);
	  gRatio->RemovePoint(0);
	  gSystErr1->RemovePoint(0);
	  gRatioSystErr->RemovePoint(0);
	}

	// Create frames for main and ratio plots
	TH1* mFrame = mc->mainFrame(hFrame,etaBin->etaBin());
	TH1* rFrame = mc->ratioFrame(hFrame,"#frac{Data}{Sim.}",0.82,1.38,etaBin->etaBin());
	rFrame->SetLineWidth(lineWidth_);
	rFrame->SetLineColor(gRes2->GetLineColor());
	
	// Labels
	TPaveText* label = labelMk_->etaBin(etaBin->etaBin(),0.45);
	mc->markForDeletion(label);
	mc->adjustPaveText(label);

	mc->canvas()->cd();
	TPad* padt = mc->mainPad(etaBin->etaBin());
	padt->Draw();
	padt->cd();
	mFrame->Draw();
	gSystErr1->Draw("PE3same");
	gRes1->Draw("PE1same");
	gRes2->Draw("PE1same");
	label->Draw("same");
	TPad* padr = mc->ratioPad(etaBin->etaBin());
	padr->Draw();
	padr->cd();
	rFrame->Draw("][");
	gRatioSystErr->Draw("PE3same");
	rFrame->Draw("][same");
	gRatio->Draw("PE1same");
	gPad->RedrawAxis();

	if( etaBin->etaBin() == 0 ) {
	  leg->AddEntry(gRes1,labelMk_->sample(sLabel1),"P");
	  leg->AddEntry(gRes2,labelMk_->sample(sLabel2),"P");
	  leg->AddEntry(gRes1,"Extrapolation Uncertainty #delta#sigma_{ex}","L");
	  leg->AddEntry(gSystErr1,"Systematic Uncertainty","F");
	}
      } // End of loop over eta bins
      
      if( par_->verbosity() > 1 ) std::cout << "    Storing MultiCanvas" << std::endl;
      mc->canvas()->cd();
      TPad* padt = mc->mainPad(5);
      padt->Draw();
      padt->cd();
      leg->Draw();      
      mc->moreLogLabelsX();
      mc->noExponentX();
      mc->setLogx();
      out_->saveCanvas(mc->canvas(),histFileName("Resolution",sLabel1,sLabel2,fitResType));
      delete mc;
      delete leg;
    } // End of loop over to-be-compared samples

    delete hFrame;
    
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotResolution(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  void PlotMaker::plotScaledMCTruth() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotScaledMCTruth(): Entering" << std::endl;
    std::cout << "Plotting scaled MC truth" << std::endl;
    
    unsigned int nComparedSamples = etaBins_.front()->nComparedSamples();
    if( nComparedSamples > 0 ) {
      out_->newPage("MCTruthScaling");

      // +++++ Sample comparison (ratios) +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      // Loop over to-be-compare Samples
      for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin(); sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	SampleLabel sLabel1 = (*sCIt)->label1();
	SampleLabel sLabel2 = (*sCIt)->label2();
	
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  const unsigned int etaBinIdx = etaBin->etaBin();
	  if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;

	  // Loop over FitResultTypes
	  for(FitResultTypeIt rIt = etaBin->fitResultTypesBegin(); rIt != etaBin->fitResultTypesEnd(); ++rIt) {
	    FitResult::Type fitResType = *rIt;
	    if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(fitResType) << std::endl;
	    
	    if( etaBin->hasKValue(sLabel1,sLabel2,fitResType) ) {

	      // Data-MC ratio (k-Value) with statistical und systematic uncertainty
	      TGraphAsymmErrors* gRatio = etaBin->ratioGraph(sLabel1,sLabel2,fitResType);
	      setStyle(sLabel1,gRatio);
	      TF1* kValueLine = etaBin->kValueLine(sLabel1,sLabel2,fitResType,"kValueLine",xMinPt_,xMaxPt_);
	      kValueLine->SetLineWidth(lineWidth_);
	      TGraphAsymmErrors* kStatBand = etaBin->kValueStatBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	      kStatBand->SetFillStyle(3445);
	      kStatBand->SetFillColor(kValueLine->GetLineColor());
	      kStatBand->SetLineColor(kValueLine->GetLineColor());
	      kStatBand->SetLineWidth(lineWidth_);
	      TGraphAsymmErrors* kSystBand = etaBin->kValueSystBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    
	      // Create frame
	      TH1* hFrame = util::HistOps::createRatioFrame(xMinPt_,xMaxPt_,0.71,2.09,"p^{true}_{T} (GeV)","#sigma(Data) / #sigma(Simulation)");
	      hFrame->GetYaxis()->SetRangeUser(0.71,2.1);
	      if( etaBinIdx >= 3 ) hFrame->GetYaxis()->SetRangeUser(0.51,3.6);
	      hFrame->SetTitle(title_);
	      hFrame->GetXaxis()->SetMoreLogLabels();
	      hFrame->GetXaxis()->SetNoExponent();
	      hFrame->SetLineWidth(lineWidth_);
	     
	      // Labels
	      TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,0.95,label->GetSize());
	      if( util::StyleSettings::style() == util::StyleSettings::PAS )
		leg->AddEntry(gRatio,"Measurement","P");
	      else
		leg->AddEntry(gRatio,"Measurement ("+labelMk_->lumi()+")","P");
	      leg->AddEntry(kValueLine,"Fitted #rho_{res}","L");
	      //	      leg->AddEntry(kStatBand,"Statistical Uncertainty","F");
	      leg->AddEntry(kStatBand,"Extrapolation Uncertainty #delta#sigma_{ex}","F");
	      leg->AddEntry(kSystBand,"Systematic Uncertainty","F");

	      out_->nextPad(sLabel1+"Over"+sLabel2+": Resolution "+etaBin->toString());
	      hFrame->Draw("][");
	      kSystBand->Draw("E2same");
	      kStatBand->Draw("E2same");
	      kValueLine->Draw("same");
	      gRatio->Draw("PE1same");
	      hFrame->Draw("][same");
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
	  } // End of loop over FitResultTypes
	} // End of loop over eta bins
      } // End of loop over SampleLabels
      
      
      if( par_->verbosity() > 1 ) std::cout << "  Sample comparison plots" << std::endl;
      // Only first FitResultType for the time being
      FitResult::Type fitResType = *(etaBins_.front()->fitResultTypesBegin());
      
      // Loop over to-be-compare Samples
      for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	  sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	SampleLabel sLabel1 = (*sCIt)->label1();
	SampleLabel sLabel2 = (*sCIt)->label2();
	if( par_->verbosity() > 1 ) std::cout << "    " << sLabel1 << " vs " << sLabel2 << ":" <<  std::endl;
      
	util::MultiCanvas* mc = new util::MultiCanvas("ResolutionRatio",3,2,5,false);
	TLegend* leg = util::LabelFactory::createLegendWithOffset(4,0.2,util::LabelFactory::lineHeightMultiCan());
	mc->adjustLegend(leg);
      
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  if( par_->verbosity() > 1 ) std::cout << "      Eta " << etaBin->etaBin() << std::endl;
	  if( etaBin->hasKValue(sLabel1,sLabel2,fitResType) ) {
	  
	    // Data-MC ratio (k-Value) with statistical und systematic uncertainty
	    TGraphAsymmErrors* gRatio = etaBin->ratioGraph(sLabel1,sLabel2,fitResType);
	    setStyle(sLabel1,gRatio);
	    mc->markForDeletion(gRatio);
	    TF1* kValueLine = etaBin->kValueLine(sLabel1,sLabel2,fitResType,"kValueLine",xMinPt_,xMaxPt_);
	    kValueLine->SetLineWidth(lineWidth_);
	    mc->markForDeletion(kValueLine);
	    TGraphAsymmErrors* kStatBand = etaBin->kValueStatBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    kStatBand->SetFillStyle(3445);
	    kStatBand->SetFillColor(kValueLine->GetLineColor());
	    kStatBand->SetLineColor(kValueLine->GetLineColor());
	    kStatBand->SetLineWidth(lineWidth_);
	    mc->markForDeletion(kStatBand);
	    TGraphAsymmErrors* kSystBand = etaBin->kValueSystBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    mc->markForDeletion(kSystBand);
	  
	    // Labels
	    TPaveText* label = labelMk_->etaBin(etaBin->etaBin(),0.45);
	    mc->adjustPaveText(label);
	    if( etaBin->etaBin() == 0 ) {
	      leg->AddEntry(gRatio,"Measurement ("+labelMk_->lumi()+")","P");
	      leg->AddEntry(kValueLine,"Fitted #rho_{res}","L");
	      leg->AddEntry(kStatBand,"Extrapolation Uncertainty #delta#sigma_{ex}","F");
	      leg->AddEntry(kSystBand,"Systematic Uncertainty","F");
	    }

	    TH1* hFrame = util::HistOps::createRatioFrame(xMinPt_,xMaxPt_,0.75,1.95,"p^{true}_{T} (GeV)","#sigma(Data) / #sigma(Simulation)  ");
	    hFrame->GetYaxis()->SetRangeUser(0.71,2.1);
	    if( etaBin->etaBin() > 3 ) hFrame->GetYaxis()->SetRangeUser(0.5,2.5);
	    hFrame->SetTitle(title_);
	    hFrame->SetLineWidth(lineWidth_);
	    mc->markForDeletion(hFrame);
	    TH1* mFrame = mc->mainFrame(hFrame,etaBin->etaBin());

	    mc->canvas()->cd();
	    TPad* padt = mc->mainPad(etaBin->etaBin());
	    padt->Draw();
	    padt->cd();
	    mFrame->Draw("][");
	    kSystBand->Draw("E2same");
	    kStatBand->Draw("E2same");
	    kValueLine->Draw("same");
	    gRatio->Draw("PE1same");
	    mFrame->Draw("][same");
	    label->Draw("same");
	    gPad->RedrawAxis();
	  } // End hasKValue
	} // End of loop over eta bins
      
	if( par_->verbosity() > 1 ) std::cout << "    Storing MultiCanvas" << std::endl;
	mc->canvas()->cd();
	TPad* padt = mc->mainPad(5);
	padt->Draw();
	padt->cd();
	leg->Draw();      
	mc->moreLogLabelsX();
	mc->noExponentX();
	mc->setLogx();
	out_->saveCanvas(mc->canvas(),histFileName("ResolutionRatio",sLabel1,sLabel2,fitResType));
	delete mc;
	delete leg;
      }


    
      // +++++ kValues vs eta +++++++++++++++++++++++++++++++++++++++++++++++++
      // Loop over FitResultTypes
      for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
	FitResult::Type fitResType = *rIt;

	// Loop over compared samples
	for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin(); sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	  SampleLabel sLabel1 = (*sCIt)->label1();
	  SampleLabel sLabel2 = (*sCIt)->label2();

	  if( etaBins_.front()->hasKValue(sLabel1,sLabel2,fitResType) ) {

	    std::vector<double> etaBinEdges;
	    etaBinEdges.push_back(par_->etaMin(0));
	    for(int i = 0; i < par_->nEtaBins(); ++i) {
	      etaBinEdges.push_back(par_->etaMax(i));
	    }
	    TH1* hKValueVsEta = new TH1D("hKValueVsEta",";|#eta|;#rho_{res}",etaBinEdges.size()-1,&(etaBinEdges.front()));
	    hKValueVsEta->GetXaxis()->SetNdivisions(505);
	    setStyle(sLabel1,hKValueVsEta);
	    hKValueVsEta->SetTitle(title_);
	    hKValueVsEta->SetLineWidth(lineWidth_);
	    TH1* hKValueDnVsEta = static_cast<TH1*>(hKValueVsEta->Clone("hKValueDnVsEta"));
	    TH1* hKValueUpVsEta = static_cast<TH1*>(hKValueVsEta->Clone("hKValueUpVsEta"));

	    // Loop over eta bins
	    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	      const EtaBin* etaBin = *etaBinIt;
	      const unsigned int etaBinIdx = etaBin->etaBin();
	      
	      hKValueVsEta->SetBinContent(etaBin->etaBin()+1,etaBin->kValue(sLabel1,sLabel2,fitResType));
	      hKValueVsEta->SetBinError(etaBin->etaBin()+1,etaBin->kStat(sLabel1,sLabel2,fitResType));
	      hKValueDnVsEta->SetBinContent(etaBin->etaBin()+1,etaBin->kValue(sLabel1,sLabel2,fitResType)-etaBin->kSystDown(sLabel1,sLabel2,fitResType));
	      hKValueUpVsEta->SetBinContent(etaBin->etaBin()+1,etaBin->kValue(sLabel1,sLabel2,fitResType)+etaBin->kSystUp(sLabel1,sLabel2,fitResType));
	    } // End of loop over eta bins

	    TGraphAsymmErrors* gKSyst = util::HistOps::getUncertaintyBand(hKValueVsEta,hKValueUpVsEta,hKValueDnVsEta,kYellow,1001);
	    gKSyst->SetLineColor(gKSyst->GetFillColor());

	    TH1* hFrame = util::HistOps::createRatioFrame(hKValueVsEta,"#rho_{res}",0.77,2.03);
	    hFrame->SetTitle(title_);

	    TPaveText* label = util::LabelFactory::createPaveText(1);
	    label->AddText(labelMk_->lumi());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,1,label->GetSize());
	    leg->AddEntry(hKValueVsEta,"Extrapolation Uncertainty #delta#sigma_{ex}","L");
	    leg->AddEntry(gKSyst,"Systematic Uncertainty","F");

	    out_->nextPad(sLabel1+"Over"+sLabel2+": Resolution");
	    hFrame->Draw();
	    gKSyst->Draw("E2same");
	    hKValueVsEta->Draw("PE1same");
	    hFrame->Draw("same");
	    label->Draw("same");
	    leg->Draw("same");
	    gPad->RedrawAxis();
	    out_->saveCurrentPad(histFileName("Rho",sLabel1,sLabel2,fitResType));
    
	    delete hFrame;
	    delete hKValueVsEta;
	    delete hKValueUpVsEta;
	    delete hKValueDnVsEta;
	    delete gKSyst;
	    delete leg;
	    delete label;
	  } // End if has kValue
	} // End of loop over compared samples
      } // End of loop over fit-result type



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
	      mcTruth->SetLineWidth(lineWidth_);

	      // Scaled MC truth and uncertainty bands
	      TF1* scaledMCT = etaBin->scaledMCTruthResoFunc("scaledMCTruthResoFunc");
	      scaledMCT->SetLineColor(kRed);
	      scaledMCT->SetLineWidth(lineWidth_);
	      TGraphAsymmErrors* scaledMCTBand = etaBin->scaledMCTruthUncertaintyBand();
	      TGraphAsymmErrors* scaledMCTRatioBand = etaBin->scaledMCTruthRatioBand();

	      // Bias corrected measurement and ratio to scaled MC truth
	      TGraphAsymmErrors* biasCorrRes = etaBin->biasCorrectedResolution(sLabel1,sLabel2,fitResType);
	      setStyle(sLabel1,biasCorrRes);
	      TGraphAsymmErrors* biasCorrResRatio = util::HistOps::createRatioGraph(biasCorrRes,scaledMCT);
	      setStyle(sLabel1,biasCorrResRatio);

	      // Labels
	      TPaveText* label = util::LabelFactory::createPaveText(1);
	      label->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+labelMk_->lumi());
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,labelMk_->start(),1);
	      if( util::StyleSettings::style() == util::StyleSettings::PAS )
		leg->AddEntry(biasCorrRes,"Bias corrected "+sLabel1,"P");
	      else
		leg->AddEntry(biasCorrRes,"Data (Bias Corrected)","P");
	      leg->AddEntry(mcTruth,"MC Truth","L");
	      leg->AddEntry(scaledMCT,"Scaled MC Truth","L");
	      leg->AddEntry(scaledMCTBand,"Systematic Uncertainty","F");

	      // Create frames for main and ratio plots
	      TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,0.32,"#sigma / p^{true}_{T}");
	      hFrameMain->SetTitle(title_);
	      hFrameMain->GetXaxis()->SetMoreLogLabels();
	      hFrameMain->GetXaxis()->SetNoExponent();

	      TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{true}_{T}","GeV","#frac{Data}{Scaled MC}",yMinResRatio_,yMaxResRatio_);
	      hFrameRatio->GetXaxis()->SetMoreLogLabels();
	      hFrameRatio->GetXaxis()->SetNoExponent();
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameRatio->SetLineColor(scaledMCT->GetLineColor());
	      hFrameRatio->SetLineStyle(scaledMCT->GetLineStyle());
	     
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
	      hFrameRatio->Draw("][");
	      scaledMCTRatioBand->Draw("E3same");
	      hFrameRatio->Draw("][same");
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
    } // End if compared samples
     
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
	    int nLegEntries = uncert->nComponents();
	    TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(std::min(nLegEntries,3),-0.5,label->GetSize());
	    TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(std::max(nLegEntries-3,0),0.5,label->GetSize());
	    
	    for(SystUncertIt it = uncert->componentsBegin(); it != uncert->componentsEnd(); ++it) {
	      TGraphAsymmErrors* b = (*it)->relUncertSteps();
		// Scale to percent
		for(int i = 0; i < b->GetN(); ++i) {
		  b->GetEYlow()[i] = 100.*b->GetEYlow()[i];
		  b->GetEYhigh()[i] = 100.*b->GetEYhigh()[i];
		}

	      
	      // 	      // Add components
	      // 	      if( (*it)->label() != "Extrapolation" && (*it)->label() != "PLI" ) {
// 		b->SetLineColor(b->GetFillColor());
// 		b->SetLineWidth(4);
// 		b->SetFillStyle(0);
// 		if( leg1->GetNRows() < 3 ) leg1->AddEntry(b,(*it)->label(),"L");
// 		else leg2->AddEntry(b,(*it)->label(),"L");
// 	      } else {
// 		if( leg1->GetNRows() < 3 ) leg1->AddEntry(b,(*it)->label(),"F");
// 		else leg2->AddEntry(b,(*it)->label(),"F");
// 	      }
// 	      bands.push_back(b);

	      // Add components in quadrature
	      if( bands.size() > 0 ) {
		for(int i = 0; i < b->GetN(); ++i) {
		  double errPrev = bands.back()->GetEYlow()[i];
		  double err = b->GetEYlow()[i];
		  b->GetEYlow()[i] = sqrt( err*err + errPrev*errPrev );
		  errPrev = bands.back()->GetEYhigh()[i];
		  err = b->GetEYhigh()[i];
		  b->GetEYhigh()[i] = sqrt( err*err + errPrev*errPrev );
		}
	      }
	      bands.push_back(b);
	      if( leg1->GetNRows() < 3 ) leg1->AddEntry(b,(*it)->label(),"F");
	      else leg2->AddEntry(b,(*it)->label(),"F");
	    }


	    // hack for the time being to catch failed fit
	    if( bands.back()->GetEYhigh()[0] > 60. ) {
	      for(std::vector<TGraphAsymmErrors*>::iterator it = bands.begin();
		  it != bands.end(); ++it) {
		(*it)->RemovePoint(0);
		(*it)->RemovePoint(1);
	      }
	    }


	    // Create frame
	    TH1* hFrame = new TH1D("hFrame",title_,1000,xMinPt_,xMaxPt_);
	    hFrame->GetXaxis()->SetMoreLogLabels();
	    hFrame->GetXaxis()->SetNoExponent();
	    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
	      hFrame->SetBinContent(bin,0.);
	    }
	    hFrame->SetLineWidth(lineWidth_);
	    hFrame->SetLineStyle(2);
	    if( etaBinIdx == 4 ) {
	      hFrame->GetYaxis()->SetRangeUser(-58,109);
	    } else {
	      hFrame->GetYaxis()->SetRangeUser(-28,54);
	    }
	    util::HistOps::setAxisTitles(hFrame,"p^{ave}_{T}","GeV","Relative Uncertainty on #sigma / p^{true}_{T} (%)");
	     
	    out_->nextPad(sampleLabel+": Systematic uncertainties "+etaBin->toString());
	    hFrame->Draw("][");
	    for(std::vector<TGraphAsymmErrors*>::reverse_iterator it = bands.rbegin();
		it != bands.rend(); ++it) {
	      (*it)->Draw("E3same");
	    }
	    hFrame->Draw("][same");
	    label->Draw("same");
	    leg1->Draw("same");
	    if( leg2->GetNRows() > 0 ) leg2->Draw("same");
	    out_->logx();
	    out_->saveCurrentPad(histFileName("RelativeSystematicUncertainty",etaBin,sampleLabel,fitResType));

	    for(std::vector<TGraphAsymmErrors*>::iterator it = bands.begin();
		it != bands.end(); ++it) {
	      delete *it;
	    }
	    delete hFrame;
	    delete label;
	    delete leg1;
	    delete leg2;
	  } // End if has uncertainties
	} // End of loop over SampleLabels
      } // End of loop over FitResultTypes
    } // End of loop over eta bins


    // +++++ MultiCanvas ++++++++++++++++++++++++++++++++++++

    // Loop over FitResultTypes
    for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin(); rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
      FitResult::Type fitResType = *rIt;
      if( par_->verbosity() > 1 ) std::cout << "    " << FitResult::toString(fitResType) << std::endl;
      
      // Loop over SampleLabels
      for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin(); sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
	SampleLabel sampleLabel = sTIt->first;
	Sample::Type sampleType = sTIt->second;
	if( par_->verbosity() > 1 ) std::cout << "      " << sampleLabel << std::endl;

	const SystematicUncertainty* unct = 0;
	if( etaBins_.front()->findSystematicUncertainty(sampleLabel,fitResType,unct) ) {
	  
	  util::MultiCanvas* mc = new util::MultiCanvas("RelativeSystematicUncertainty",3,2,5,false);
	  TLegend* leg = util::LabelFactory::createLegendWithOffset(unct->nComponents(),0.1+util::LabelFactory::lineHeightMultiCan(),util::LabelFactory::lineHeightMultiCan());
	  mc->adjustLegend(leg);
	  mc->markForDeletion(leg);
	  
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;
	    const unsigned int etaBinIdx = etaBin->etaBin();
	    if( par_->verbosity() > 1 ) std::cout << "  Eta " << etaBinIdx << std::endl;
	    
	    // If uncertainties are assigned to this sample and 
	    // fit result type, get them
	    const SystematicUncertainty* uncert = 0;
	    if( etaBin->findSystematicUncertainty(sampleLabel,fitResType,uncert) ) {
	      
	      // Uncertainty bands and labels for components
	      // and total uncertainty
	      std::vector<TGraphAsymmErrors*> bands;
	      for(SystUncertIt it = uncert->componentsBegin(); it != uncert->componentsEnd(); ++it) {
		TGraphAsymmErrors* b = (*it)->relUncertSteps();
		// Scale to percent
		for(int i = 0; i < b->GetN(); ++i) {
		  b->GetEYlow()[i] = 100.*b->GetEYlow()[i];
		  b->GetEYhigh()[i] = 100.*b->GetEYhigh()[i];
		}
		
		// Add components in quadrature
		if( bands.size() > 0 ) {
		  for(int i = 0; i < b->GetN(); ++i) {
		    double errPrev = bands.back()->GetEYlow()[i];
		    double err = b->GetEYlow()[i];
		    b->GetEYlow()[i] = sqrt( err*err + errPrev*errPrev );
		    errPrev = bands.back()->GetEYhigh()[i];
		    err = b->GetEYhigh()[i];
		    b->GetEYhigh()[i] = sqrt( err*err + errPrev*errPrev );
		  }
		}
		bands.push_back(b);
		mc->markForDeletion(b);
		if( etaBin->etaBin() == 0 ) {
		  leg->AddEntry(b,(*it)->label(),"F");
		}
	      }

	      
	      // hack for the time being to catch failed fit
	      if( bands.back()->GetEYhigh()[0] > 60. ) {
		for(std::vector<TGraphAsymmErrors*>::iterator it = bands.begin();
		    it != bands.end(); ++it) {
		  (*it)->RemovePoint(0);
		  (*it)->RemovePoint(1);
		}
	      }
	      
	      TPaveText* label = labelMk_->etaBin(etaBin->etaBin(),0.45);
	      mc->adjustPaveText(label);
	      mc->markForDeletion(label);
	      
	      // Create frame
	      TH1* hFrame = util::HistOps::createRatioFrame(xMinPt_,xMaxPt_,-38,38,"p^{ave}_{T} (GeV)","Rel. Uncertainty on #rho_{res} (%)");
	      for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
		hFrame->SetBinContent(bin,0.);
	      }
	      hFrame->SetLineWidth(lineWidth_);
	      if( etaBinIdx == 4 ) {
		hFrame->GetYaxis()->SetRangeUser(-68,68);
	      }
	      util::HistOps::setAxisTitles(hFrame,"p^{ave}_{T}","GeV","Rel. Uncertainty on #sigma / p^{true}_{T} (%)");
	      mc->markForDeletion(hFrame);
	      TH1* mFrame = mc->mainFrame(hFrame,etaBin->etaBin());
	      
	      mc->canvas()->cd();
	      TPad* padt = mc->mainPad(etaBin->etaBin());
	      padt->Draw();
	      padt->cd();
	      mFrame->Draw("][");
	      for(std::vector<TGraphAsymmErrors*>::reverse_iterator it = bands.rbegin();
		  it != bands.rend(); ++it) {
		(*it)->Draw("E3same");
	      }
	      mFrame->Draw("][same");
	      label->Draw("same");
	      gPad->RedrawAxis();
	    } // End has systematic uncertainty
	  } // End of loop over eta bins
	  
	  if( par_->verbosity() > 1 ) std::cout << "    Storing MultiCanvas" << std::endl;
	  mc->canvas()->cd();
	  TPad* padt = mc->mainPad(5);
	  padt->Draw();
	  padt->cd();
	  leg->Draw();      
	  mc->moreLogLabelsX();
	  mc->noExponentX();
	  mc->setLogx();
	  out_->saveCanvas(mc->canvas(),histFileName("RelativeSystematicUncertainty",sampleLabel));
	  delete mc;
	}
      }	// End of loop over sample labels
    } // End of loop over fit result types

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotSystematicUncertainties(): Leaving" << std::endl;    
  }


  // MC event information:
  // - Event weights;
  // - Number of added PU interactions
  // -------------------------------------------------------------------------------------
  void PlotMaker::plotMCEventInfo() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotMCEventInfo(): Entering" << std::endl;
    std::cout << "Plotting MC event information" << std::endl;

    // +++++ MC Weights ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Loop over eta and pt bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;

      out_->newPage("MCEventWeights");

      for(PtBinIt ptBinIt = (*etaBinIt)->ptBinsBegin(); ptBinIt != (*etaBinIt)->ptBinsEnd(); ++ptBinIt) {
	const PtBin* ptBin = *ptBinIt;
	if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
	
	// Loop over MCSamples in this bin
	for(MCSampleIt sIt = ptBin->mcSamplesBegin(); sIt != ptBin->mcSamplesEnd(); ++sIt) {
	  const MCSample* sample = sIt->second;
	  if( par_->verbosity() > 1 ) std::cout << "    " << sample->label() << std::endl;
	  
	  // Labels
	  TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(1+label->GetSize(),labelMk_->start()-0.2,label->GetSize());	    

	  // Weights in different ptSoft bins
	  std::vector<TH1*> hWeights;
	  unsigned int step = (par_->nPtSoftBins() > 3 ? 2 : 1);
	  unsigned int i = 0;
	  for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ptSoftBinIdx += step, ++i) {
	    TH1* h = sample->histMCWeight(ptSoftBinIdx);
	    setStyle(sample,h);
	    h->SetLineColor(util::StyleSettings::color(i));
	    util::HistOps::normHist(h,"width");
	    hWeights.push_back(h);
	    leg->AddEntry(h,labelMk_->ptSoftRange(ptSoftBinIdx),"L");
	  }

	  // Style
	  double min = hWeights.back()->GetMinimum(1E-15);
	  for(std::vector<TH1*>::iterator it = hWeights.begin(); it != hWeights.end(); ++it) {
	    util::HistOps::setAxisTitles(*it,"Event Weights","","events",true);
	    (*it)->GetXaxis()->SetRangeUser(0.,1.3);
	    util::HistOps::setYRange(*it,label->GetSize()+leg->GetNRows(),min);
	  }

	  // Plots
	  out_->nextMultiPad(sample->label()+": MC Event Weights "+ptBin->toTString());
	  hWeights.front()->Draw("HIST");
	  for(std::vector<TH1*>::iterator it = hWeights.begin()+1; it != hWeights.end(); ++it) {
	    (*it)->Draw("HISTsame");
	  }
	  label->Draw("same");
	  leg->Draw("same");
	  out_->logy();
	  gPad->RedrawAxis();
	  out_->saveCurrentPad(histFileName("MCWeights",ptBin,sample));
	  
	  for(std::vector<TH1*>::iterator it = hWeights.begin(); it != hWeights.end(); ++it) {
	    delete *it;
	  }
	  delete leg;
	  delete label;
	} // End of loop over Samples
      } // End of loop over pt bins
    } // End of loop over eta bins



    // +++++ MC PU Distribution ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Loop over eta and pt bins
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;

      out_->newPage("MCPUDistributions");

      for(PtBinIt ptBinIt = (*etaBinIt)->ptBinsBegin(); ptBinIt != (*etaBinIt)->ptBinsEnd(); ++ptBinIt) {
	const PtBin* ptBin = *ptBinIt;
	if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
	
	// Loop over MCSamples in this bin
	for(MCSampleIt sIt = ptBin->mcSamplesBegin(); sIt != ptBin->mcSamplesEnd(); ++sIt) {
	  const MCSample* sample = sIt->second;
	  if( par_->verbosity() > 1 ) std::cout << "    " << sample->label() << std::endl;
	  
	  // Labels
	  TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(1+label->GetSize(),labelMk_->start()-0.2,label->GetSize());	    

	  // NPU in different ptSoft bins
	  std::vector<TH1*> hNPU;
	  unsigned int step = (par_->nPtSoftBins() > 3 ? 2 : 1);
	  unsigned int i = 0;
	  for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ptSoftBinIdx += step, ++i) {
	    TH1* h = sample->histMCNumPU(ptSoftBinIdx);
	    setStyle(sample,h);
	    h->SetLineColor(util::StyleSettings::color(i));
	    hNPU.push_back(h);
	    leg->AddEntry(h,labelMk_->ptSoftRange(ptSoftBinIdx),"L");
	  }

	  // Style
	  for(std::vector<TH1*>::iterator it = hNPU.begin(); it != hNPU.end(); ++it) {
	    util::HistOps::normHist(*it,"width");
	    util::HistOps::setAxisTitles(*it,"N(PU) Interactions","","events",true);
	    (*it)->GetXaxis()->SetRangeUser(0.,25.);
	    util::HistOps::setYRange(*it,label->GetSize()+leg->GetNRows()+1);
	  }

	  // Plots
	  out_->nextMultiPad(sample->label()+": MC Num PU "+ptBin->toTString());
	  hNPU.front()->Draw("HIST");
	  for(std::vector<TH1*>::iterator it = hNPU.begin()+1; it != hNPU.end(); ++it) {
	    (*it)->Draw("HISTsame");
	  }
	  label->Draw("same");
	  leg->Draw("same");
	  gPad->RedrawAxis();
	  out_->saveCurrentPad(histFileName("MCNPU",ptBin,sample));
	  
	  for(std::vector<TH1*>::iterator it = hNPU.begin(); it != hNPU.end(); ++it) {
	    delete *it;
	  }
	  delete leg;
	  delete label;
	} // End of loop over Samples
      } // End of loop over pt bins
    } // End of loop over eta bins


    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotMCEventInfo(): Leaving" << std::endl;
  }


  //! \brief Plot miscellaneous control distributions
  //!
  //! Produces the following control distributions
  //! - Number of reconstructed vertices:
  //!    comparison of all samples marked for direct comparison
  //!    (\sa CommanderCool::compareSamples)
  // -------------------------------------------------------------------------------------
  void PlotMaker::plotControlDistributions() const {
    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotControlDistributions(): Entering" << std::endl;

    // +++++ Number of reconstructed vertices per (ptSoft,eta,pt) bin, different Samples ++
    
    if( etaBins_.front()->nComparedSamples() > 0 ) {
      // Loop over ptSoft bins
      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
	if( par_->verbosity() > 1 ) std::cout << "  PtSoftBin " << ptSoftBinIdx << std::endl;
	  
	// Loop over eta bins
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  out_->newPage("Number of Vertices");
	    
	  // Loop over pt bins
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    const PtBin* ptBin = *ptBinIt;
	    if( par_->verbosity() > 1 ) std::cout << "    " << ptBin->toTString() << std::endl;

	    // Loop over to-be-compare Samples
	    for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
		sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
	      SampleLabel sLabel1 = (*sCIt)->label1();
	      SampleLabel sLabel2 = (*sCIt)->label2();
	      const Sample* sample1 = ptBin->findSample(sLabel1);
	      const Sample* sample2 = ptBin->findSample(sLabel2);

	      // NVtx distributions and fits
	      TH1* hNVtx1 = sample1->histNumVtx(ptSoftBinIdx);
	      hNVtx1->Sumw2();
	      if( hNVtx1->Integral() ) hNVtx1->Scale(1./hNVtx1->Integral("width"));
	      util::HistOps::setAxisTitles(hNVtx1,"N_{Vtx}","","events",true);
	      setStyle(sample1,hNVtx1);

	      TH1* hNVtx2 = sample2->histNumVtx(ptSoftBinIdx);
	      hNVtx2->Sumw2();
	      if( hNVtx2->Integral() ) hNVtx2->Scale(1./hNVtx2->Integral("width"));
	      util::HistOps::setAxisTitles(hNVtx2,"N_{Vtx}","","events",true);
	      setStyle(sample2,hNVtx2);

	      // Labels
	      TPaveText* label = labelMk_->etaPtAvePtSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	      TLegend* leg = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
	      leg->AddEntry(hNVtx1,"#LT"+labelMk_->sample(sLabel1)+"#GT = "+util::toTString(hNVtx1->GetMean(),3)+" #pm "+util::toTString(hNVtx1->GetMeanError(),3),"P");
	      leg->AddEntry(hNVtx2,"#LT"+labelMk_->sample(sLabel2)+"#GT = "+util::toTString(hNVtx2->GetMean(),3)+" #pm "+util::toTString(hNVtx2->GetMeanError(),3),"L");
	    
	      util::HistOps::setYRange(hNVtx2,label->GetSize()+leg->GetNRows());
	      hNVtx2->GetXaxis()->SetRange(1,19);
	      out_->nextMultiPad(sLabel1+" vs "+sLabel2+": NVtx "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hNVtx2->Draw("HIST");
	      hNVtx1->Draw("PE1same");
	      label->Draw("same");
	      leg->Draw("same");
	      out_->saveCurrentPad(histFileName("NumVtx",ptBin,sLabel1,sLabel2,ptSoftBinIdx));

	      delete hNVtx1;
	      delete hNVtx2;
	      delete label;
	      delete leg;
	    } // End of loop over to-be-compared samples
	  } // End of loop over pt bins
	} // End of loop over eta bins
      } // End of loop over ptSoft bins
    } // End if samples to be compared


//     // +++++ Number of reconstructed vertices per (ptSoft,eta) bin, different Samples +++++++++++++
    
//     if( etaBins_.front()->nComparedSamples() > 0 ) {
//       // Loop over ptSoft bins
//       for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
// 	if( par_->verbosity() > 1 ) std::cout << "  PtSoftBin " << ptSoftBinIdx << std::endl;
	  
// 	// Loop over eta bins
// 	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
// 	  const EtaBin* etaBin = *etaBinIt;
// 	  if( par_->verbosity() > 1 ) std::cout << "    " << etaBin->toString() << std::endl;
// 	  out_->newPage("Number of Vertices");
	    
// 	  // Loop over to-be-compare Samples
// 	  for(ComparedSamplesIt sCIt = etaBin->comparedSamplesBegin();
// 	      sCIt != etaBin->comparedSamplesEnd(); ++sCIt) {
// 	    SampleLabel sLabel1 = (*sCIt)->label1();
// 	    SampleLabel sLabel2 = (*sCIt)->label2();

// 	    TH1* hNVtx1 = 0;
// 	    TH1* hNVtx2 = 0;

// 	    // Loop over pt bins
// 	    for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
// 	      const PtBin* ptBin = *ptBinIt;
// 	      const Sample* sample1 = ptBin->findSample(sLabel1);
// 	      const Sample* sample2 = ptBin->findSample(sLabel2);
	      
// 	      // NVtx distributions and fits
// 	      // Sample1
// 	      if( hNVtx1 == 0 ) {
// 		hNVtx1 = sample1->histNumVtx(ptSoftBinIdx);
// 		hNVtx1->Sumw2();
// 		util::HistOps::setAxisTitles(hNVtx1,"N(reconstructed vertices)","","events",true);
// 		setStyle(sample1,hNVtx1);
// 	      } else {
// 		TH1* h = sample1->histNumVtx(ptSoftBinIdx);
// 		hNVtx1->Add(h);
// 		delete h;
// 	      }

// 	      // Sample2
// 	      // Scale of MC to data
// 	      double w = 1.;
// 	      if( sample1->type() == Sample::Data && sample2->type() == Sample::MC ) {
// 		w = sample2->relativeWeightTo(sLabel1,ptSoftBinIdx);
// 	      }
	      
// 	      if( hNVtx2 == 0 ) {
// 		hNVtx2 = sample2->histNumVtx(ptSoftBinIdx);
// 		hNVtx2->Sumw2();
// 		hNVtx2->Scale(w);
// 		util::HistOps::setAxisTitles(hNVtx2,"N(reconstructed vertices)","","events",true);
// 		setStyle(sample2,hNVtx2);
// 	      } else {
// 		TH1* h = sample2->histNumVtx(ptSoftBinIdx);
// 		h->Scale(w);
// 		hNVtx2->Add(h);
// 		delete h;
// 	      }
// 	    } // End of loop over pt bins

// 	    if( hNVtx1 && hNVtx2 ) {
// 	      if( hNVtx1->Integral() ) hNVtx1->Scale(1./hNVtx1->Integral("width"));
// 	      if( hNVtx2->Integral() ) hNVtx2->Scale(1./hNVtx2->Integral("width")); 

// 	      // Labels
// 	      TPaveText* label = labelMk_->ptSoftBin(etaBin->etaBin(),ptSoftBinIdx);
// 	      TLegend* leg = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
// 	      leg->AddEntry(hNVtx1,sLabel1+": #LTN#GT = "+util::toTString(hNVtx1->GetMean(),3)+" #pm "+util::toTString(hNVtx1->GetMeanError(),3),"P");
// 	      leg->AddEntry(hNVtx2,sLabel2+": #LTN#GT = "+util::toTString(hNVtx2->GetMean(),3)+" #pm "+util::toTString(hNVtx2->GetMeanError(),3),"L");
	    
// 	      util::HistOps::setYRange(hNVtx2,label->GetSize()+leg->GetNRows());
// 	      hNVtx2->GetXaxis()->SetRange(1,25);
// 	      out_->nextMultiPad(sLabel1+" vs "+sLabel2+": NVtx PtSoftBin "+util::toTString(ptSoftBinIdx));
// 	      hNVtx2->Draw("HIST");
// 	      hNVtx1->Draw("PE1same");
// 	      label->Draw("same");
// 	      leg->Draw("same");
// 	      out_->saveCurrentPad(histFileName("NumVtx",sLabel1,sLabel2,etaBin,ptSoftBinIdx));

// 	      delete hNVtx1;
// 	      delete hNVtx2;
// 	      delete label;
// 	      delete leg;
// 	    }
// 	  } // End of loop over to-be-compared samples
// 	} // End of loop over eta bins
//       } // End of loop over ptSoft bins
//     } // End if samples to be compared


    // +++++ Delta Pt +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if( etaBins_.front()->nComparedSamples() > 0 ) {
      // Loop over SampleLabels
      for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	  sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
	SampleLabel sampleLabel = sTIt->first;
	if( par_->verbosity() > 1 ) std::cout << sampleLabel << std::endl;

	// Plot spectra only for samples that are compared
	bool isToBeCompared = false;
	for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	    sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	  if( (*sCIt)->contains(sampleLabel) ) {
	    isToBeCompared = true;
	    break;
	  }
	}

	if( isToBeCompared ) {
	  // Loop over eta and pt bins
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;
	    out_->newPage("Delta Pt");
	  
	    for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	      const PtBin* ptBin = *ptBinIt;
	      if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
	    
	      // Loop over Samples and select current one (by sampleLabel)
	      for(SampleIt sIt = ptBin->samplesBegin(); sIt != ptBin->samplesEnd(); ++sIt) {
		const Sample* sample = sIt->second;
		if( sample->label() != sampleLabel ) continue;
	      
		for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < sample->nPtSoftBins(); ++ptSoftBinIdx) {
		  if( par_->verbosity() > 1 ) std::cout << "      PtSoftBin " << ptSoftBinIdx << std::endl;
		
		  // delta pt of first and second jet
		  TH1* hDeltaPt = sample->histDeltaPt(ptSoftBinIdx);
		  hDeltaPt->GetXaxis()->SetNdivisions(505);
		  hDeltaPt->GetXaxis()->SetRangeUser(0.,0.7*hDeltaPt->GetXaxis()->GetBinUpEdge(hDeltaPt->GetNbinsX()));
		  util::HistOps::normHist(hDeltaPt,"width");
		  util::HistOps::setAxisTitles(hDeltaPt,"#Deltap_{T} = #frac{1}{2}#upoint(p_{T,1} - p_{T,2})","GeV","events",true);
		  setStyle(sample,hDeltaPt);
		  double yMin = 7E-6;
		
		  // Illustrate start of tails and motivate
		  // DeltaPtMax cut
		
		  // Iterative Gaussian fit
		  TF1* fit = new TF1("fit","gaus",0.,hDeltaPt->GetXaxis()->GetBinUpEdge(hDeltaPt->GetNbinsX()));
		  hDeltaPt->Fit(fit,"Q0IR");
		  fit->SetRange(0.,2.*fit->GetParameter(2));
		  hDeltaPt->Fit(fit,"Q0IR");
		  fit->SetRange(0.,2.*fit->GetParameter(2));
		  fit->SetLineWidth(2);
		  fit->SetLineColor(hDeltaPt->GetLineColor());
		  TF1* fitEx = static_cast<TF1*>(fit->Clone("fitEx"));
		  fitEx->SetLineStyle(2);
		  fitEx->SetRange(0.,fitEx->GetX(yMin,10.,200.));
		
		  // n sig lines
		  std::vector<TLine*> lines;
		  for(unsigned int i = 0; i < 4; ++i) {
		    double x = (2+i)*fit->GetParameter(2);
		    lines.push_back(new TLine(x,yMin,x,0.1));
		    lines.back()->SetLineWidth(4);
		    lines.back()->SetLineStyle(3);
		    lines.back()->SetLineColor(kOrange+1);
		  }	      
		
		  // Labels
		  TPaveText* label = labelMk_->etaPtAvePtSoftBin(etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
		
		  TLegend* leg = util::LabelFactory::createLegendColWithOffset(3,0.7*labelMk_->start(),label->GetSize());
		  leg->AddEntry(hDeltaPt,labelMk_->sample(sample->label()),"P");
		  leg->AddEntry(fit,"Gaussian Fit","L");
		  leg->AddEntry(lines.front(),"n#upoint#sigma',  n = 2,3,...","L");
		
		  util::HistOps::setYRange(hDeltaPt,label->GetSize()+leg->GetNRows()-1,yMin);
		  out_->nextMultiPad(sample->label()+": Delta Pt "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		  hDeltaPt->Draw("PE1");
		  fitEx->Draw("same");
		  fit->Draw("same");
		  label->Draw("same");
		  for(std::vector<TLine*>::iterator it = lines.begin(); it != lines.end(); ++it) {
		    (*it)->Draw("same");
		  }
		  leg->Draw("same");
		  gPad->RedrawAxis();
		  out_->logy();
		  out_->saveCurrentPad(histFileName("DeltaPt-CutIllustration",ptBin,sample,ptSoftBinIdx));

		
		  // Superimpose fit result
		  // (first FitResult type only)
		  FitResult::Type fitResType = *(etaBin->fitResultTypesBegin());
// 		  double fittedSigma = sample->meanPt(fitResType)*sample->fittedValue(fitResType,ptSoftBinIdx)/sqrt(2.);
// 		  double uncertSigma = sample->meanPt(fitResType)*sample->fittedUncert(fitResType,ptSoftBinIdx)/sqrt(2.);
		  double fittedSigma = sample->fittedValue(fitResType,ptSoftBinIdx)/sqrt(2.);
		  double uncertSigma = sample->fittedUncert(fitResType,ptSoftBinIdx)/sqrt(2.);
		  TF1* fitResult = static_cast<TF1*>(fit->Clone("fitResult"));
		  // Norm to data distribution: the data distribution is normalised
		  // to its integral, the pdf from 0 to 2sigma.
		  fitResult->SetParameter(0,fitEx->Integral(0.,2.*fittedSigma)/(sqrt(M_PI/2.)*fittedSigma*erf(sqrt(2.))));
		  fitResult->SetParameter(1,0.);
		  fitResult->SetParameter(2,fittedSigma);
		  fitResult->SetRange(0.,2.*fittedSigma);
		  fitResult->SetLineWidth(2);
		  fitResult->SetLineColor(kRed);
		  TF1* fitResultEx = static_cast<TF1*>(fitResult->Clone("fitResultEx"));
		  fitResultEx->SetLineStyle(2);
		  fitResultEx->SetRange(0.,fitResultEx->GetX(yMin,10.,200.));
		
		  TLegend* legResult = util::LabelFactory::createLegendColWithOffset(3,labelMk_->start(),label->GetSize());
		  legResult->AddEntry(hDeltaPt,labelMk_->sample(sample->label()),"P");
		  legResult->AddEntry(fitResult,"Estimate","L");
		  legResult->AddEntry(fitResult,"#hat{#sigma}' = "+util::toTString(fittedSigma,2)+" #pm "+util::toTString(uncertSigma,2)+" GeV","");
		
		  util::HistOps::setYRange(hDeltaPt,label->GetSize());
		  out_->nextMultiPad(sample->label()+": Delta Pt "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		  hDeltaPt->Draw("PE1");
		  fitResultEx->Draw("same");
		  fitResult->Draw("same");
		  label->Draw("same");
		  legResult->Draw("same");
		  gPad->RedrawAxis();
		  out_->saveCurrentPad(histFileName("DeltaPt",ptBin,sample,ptSoftBinIdx));
		  util::HistOps::setYRange(hDeltaPt,label->GetSize(),yMin);
		  out_->nextMultiPad(sample->label()+": Delta Pt "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		  hDeltaPt->Draw("PE1");
		  fitResultEx->Draw("same");
		  fitResult->Draw("same");
		  label->Draw("same");
		  legResult->Draw("same");
		  gPad->RedrawAxis();
		  out_->logy();
		  out_->saveCurrentPad(histFileName("DeltaPtLog",ptBin,sample,ptSoftBinIdx));
		
		
		  delete hDeltaPt;
		  delete fit;
		  delete fitEx;
		  delete fitResult;
		  delete fitResultEx;
		  delete leg;
		  delete legResult;
		  delete label;
		  for(std::vector<TLine*>::iterator it = lines.begin();
		      it != lines.end(); ++it) {
		    delete *it;
		  }		
		} // End isToBeCompared
	      } // End of loop over pt soft bins
	    } // End of loop over SampleLabels
	  } // End of loop over pt bins
	} // End of loop over eta bins
      }
    }

    if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotControlDistributions(): Leaving" << std::endl;
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, SampleLabel label) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label);
  }

  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, SampleLabel label1, SampleLabel label2, const EtaBin* etaBin, unsigned int ptSoftBin) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_EtaBin"+util::toTString(etaBin->etaBin())+"_PtSoftBin"+util::toTString(ptSoftBin));
  }

  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, SampleLabel label1, SampleLabel label2, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+histFileNameFitResultTypeLabel(type));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const EtaBin* etaBin, SampleLabel label1, SampleLabel label2, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_"+histFileNameFitResultTypeLabel(type)+"EtaBin"+util::toTString(etaBin->etaBin()));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const EtaBin* etaBin, SampleLabel sampleLabel, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sampleLabel+"_"+histFileNameFitResultTypeLabel(type)+"EtaBin"+util::toTString(etaBin->etaBin()));
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sample->label()+"_"+histFileNameFitResultTypeLabel(type)+"EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin()));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, SampleLabel label1, SampleLabel label2, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_"+histFileNameFitResultTypeLabel(type)+"EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin()));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sample->label()+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin()));
  }
  

  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, unsigned int ptSoftBinIdx) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sample->label()+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, SampleLabel label1, SampleLabel label2, unsigned int ptSoftBinIdx) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, SampleLabel label1, SampleLabel label2, FitResult::Type type, unsigned int ptSoftBinIdx) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_"+histFileNameFitResultTypeLabel(type)+"EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileNameFitResultTypeLabel(FitResult::Type type) const {
    return etaBins_.front()->nFitResultTypes() > 1 ? FitResult::toString(type)+"_" : "";
  }

  // -------------------------------------------------------------------------------------
  TString PlotMaker::cleanFileName(TString str) const {
    str.ReplaceAll(" ","");
    str.ReplaceAll("(","");
    str.ReplaceAll(")","");
    str.ReplaceAll("#","");
    str.ReplaceAll("<","");
    str.ReplaceAll(">","");
    str.ReplaceAll("{","");
    str.ReplaceAll("}","");
    str.ReplaceAll("__","_");

    return str;
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::setStyle(const SampleLabel &label, TH1* &h) const {
    h->SetMarkerStyle(markerStyle(label));
    h->SetMarkerColor(color(label));
    h->SetLineColor(h->GetMarkerColor());
    h->SetTitle(title_);
    h->SetLineWidth(lineWidth_);
    h->SetMarkerSize(markerSize_);
  }


  // -------------------------------------------------------------------------------------
  void PlotMaker::setStyle(const SampleLabel &label, TGraphAsymmErrors* &g) const {
    g->SetMarkerStyle(markerStyle(label));
    g->SetMarkerColor(color(label));
    g->SetLineColor(g->GetMarkerColor());
    g->SetTitle(title_);
    g->SetLineWidth(lineWidth_);
    g->SetMarkerSize(markerSize_);
  }



  // -------------------------------------------------------------------------------------
  PlotMaker::LabelMaker::LabelMaker(const Parameters* par)
    : par_(par) {}


  // -------------------------------------------------------------------------------------
  double PlotMaker::LabelMaker::start() const {
    double result = 0.75;
    if( util::StyleSettings::style() == util::StyleSettings::Presentation ) result = 0.9;

    return result;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtAvePtSoftBin(unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(2);
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+util::LabelFactory::ptAveCut(par_->ptMin(etaBinIdx,ptBinIdx),par_->ptMax(etaBinIdx,ptBinIdx)));
    lab->AddText(util::LabelFactory::deltaPhiCut(par_->deltaPhi())+",  "+util::LabelFactory::pt3RelCut(par_->ptSoftMax(ptSoftBinIdx)));

    return lab;
  }

  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtAvePtSoftBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(3);
    lab->AddText(sample(label));
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+util::LabelFactory::ptAveCut(par_->ptMin(etaBinIdx,ptBinIdx),par_->ptMax(etaBinIdx,ptBinIdx)));
    lab->AddText(util::LabelFactory::deltaPhiCut(par_->deltaPhi())+",  "+util::LabelFactory::pt3RelCut(par_->ptSoftMax(ptSoftBinIdx)));

    return lab;
  }

  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtAvePtSoftGenBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(3);
    lab->AddText(sample(label));
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+util::LabelFactory::ptAveCut(par_->ptMin(etaBinIdx,ptBinIdx),par_->ptMax(etaBinIdx,ptBinIdx)));
    lab->AddText(util::LabelFactory::deltaPhiCut(par_->deltaPhi())+",  "+util::LabelFactory::pt3RelGenCut(par_->ptSoftMax(ptSoftBinIdx)));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtAveBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(3);
    lab->AddText(sample(label));
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+util::LabelFactory::ptAveCut(par_->ptMin(etaBinIdx,ptBinIdx),par_->ptMax(etaBinIdx,ptBinIdx)));
    lab->AddText(util::LabelFactory::deltaPhiCut(par_->deltaPhi()));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtAveBin(unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(2);
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+util::LabelFactory::ptAveCut(par_->ptMin(etaBinIdx,ptBinIdx),par_->ptMax(etaBinIdx,ptBinIdx)));
    lab->AddText(util::LabelFactory::deltaPhiCut(par_->deltaPhi()));

    return lab;
  }
 
  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtSoftBin(unsigned int etaBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(1);
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx))+",  "+util::LabelFactory::deltaPhiCut(par_->deltaPhi())+",  "+util::LabelFactory::pt3RelCut(par_->ptSoftMax(ptSoftBinIdx)));

    return lab;
  }

  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaPtSoftBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(3);
    lab->AddText(label);
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx)));
    lab->AddText(util::LabelFactory::deltaPhiCut(par_->deltaPhi())+",  "+util::LabelFactory::pt3RelCut(par_->ptSoftMax(ptSoftBinIdx)));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaBin(SampleLabel label, unsigned int etaBinIdx, double pos) const {
    TPaveText* lab = util::LabelFactory::createPaveText(1,pos);
    lab->AddText(sample(label)+",  "+util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx)));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaBin(unsigned int etaBinIdx, double pos) const {
    TPaveText* lab = util::LabelFactory::createPaveText(1,pos);
    lab->AddText(util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx)));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::sampleBin(SampleLabel label) const {
    TPaveText* lab = util::LabelFactory::createPaveText(1);
    lab->AddText(sample(label)+",  "+util::LabelFactory::deltaPhiCut(par_->deltaPhi()));

    return lab;
  }




  // deprecated style functions
  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptSoftBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx, unsigned int nExtraEntries) const {
    TPaveText* lab = etaBin(label,etaBinIdx,1+nExtraEntries);
    lab->AddText(ptRange(etaBinIdx,ptBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptSoftBin(unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(3);
    lab->AddText("#sqrt{s} = 7 TeV,  "+jets());
    lab->AddText(etaRange(etaBinIdx)+",  "+ptRange(etaBinIdx,ptBinIdx));
    lab->AddText(ptSoftRange(ptSoftBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptSoftBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = 0;
    if( util::StyleSettings::style() != util::StyleSettings::PAS ) {
      if( Sample::type(label) == Sample::Data ) {
	lab = util::LabelFactory::createPaveText(3);
	lab->AddText("#sqrt{s} = 7 TeV,  L = "+lumi());
	lab->AddText(jets());
	lab->AddText(etaRange(etaBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));
      } else {
	lab = util::LabelFactory::createPaveText(2);
	lab->AddText("#sqrt{s} = 7 TeV,  "+jets());
	lab->AddText(etaRange(etaBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));
      }
    } else {
      lab = util::LabelFactory::createPaveText(2);
      lab->AddText("#sqrt{s} = 7 TeV,  "+jets());
      lab->AddText(etaRange(etaBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));
    }

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptSoftBin(unsigned int etaBinIdx, unsigned int ptSoftBinIdx) const {
    TPaveText* lab = util::LabelFactory::createPaveText(2);
    lab->AddText("#sqrt{s} = 7 TeV,  "+jets());
    lab->AddText(etaRange(etaBinIdx)+",  "+ptSoftRange(ptSoftBinIdx));

    return lab;
  }



  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptBin(unsigned int etaBinIdx, unsigned int ptBinIdx, unsigned int nExtraEntries) const {
    TPaveText* lab = util::LabelFactory::createPaveText(2+nExtraEntries);
    lab->AddText("#sqrt{s} = 7 TeV,  "+jets());
    lab->AddText(etaRange(etaBinIdx)+", "+ptRange(etaBinIdx,ptBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::ptBin(SampleLabel label, unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    TPaveText* lab = etaBin(label,etaBinIdx,(unsigned int)1);
    lab->AddText(ptRange(etaBinIdx,ptBinIdx));

    return lab;
  }


   // -------------------------------------------------------------------------------------
   TPaveText* PlotMaker::LabelMaker::etaBin(SampleLabel label, unsigned int etaBinIdx, unsigned int nExtraEntries) const {
     TPaveText* lab = 0;
     if( util::StyleSettings::style() != util::StyleSettings::PAS ) {
       lab = util::LabelFactory::createPaveText(2+nExtraEntries);
       if( Sample::type(label) == Sample::Data ) {
 	lab->AddText(label+",  #sqrt{s} = 7 TeV, L = "+lumi());
       } else {
 	lab->AddText(label+",  #sqrt{s} = 7 TeV");
       }
     } else {
       lab = util::LabelFactory::createPaveText(1+nExtraEntries);
     }
     lab->AddText(jets()+",  "+etaRange(etaBinIdx));

     return lab;
   }


   // -------------------------------------------------------------------------------------
   TPaveText* PlotMaker::LabelMaker::inclusiveLabel(SampleLabel label, unsigned int nExtraEntries) const {
     TPaveText* lab = 0;
     if( util::StyleSettings::style() != util::StyleSettings::PAS ) {
       lab = util::LabelFactory::createPaveText(2+nExtraEntries);
       if( Sample::type(label) == Sample::Data ) {
 	lab->AddText(label+",  #sqrt{s} = 7 TeV, L = "+lumi());
       } else {
 	lab->AddText(label+",  #sqrt{s} = 7 TeV");
       }
     } else {
       lab = util::LabelFactory::createPaveText(1+nExtraEntries);
     }
     lab->AddText(jets());

     return lab;
   }


   // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaBin(unsigned int etaBinIdx, unsigned int nExtraEntries) const {
    TPaveText* lab = util::LabelFactory::createPaveText(1+nExtraEntries);
    lab->AddText("#sqrt{s} = 7 TeV,  "+jets());
    lab->AddText(etaRange(etaBinIdx));
    
    return lab;
  }


   // -------------------------------------------------------------------------------------
   TString PlotMaker::LabelMaker::etaRange(unsigned int etaBinIdx, const TString &superscript) const {
     return util::LabelFactory::etaCut(par_->etaMin(etaBinIdx),par_->etaMax(etaBinIdx),superscript);
   }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptRange(unsigned int etaBinIdx, unsigned int ptBinIdx) const {
    return util::LabelFactory::ptAveCut(par_->ptMin(etaBinIdx,ptBinIdx),par_->ptMax(etaBinIdx,ptBinIdx));
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptSoftRange(unsigned int ptSoftBinIdx) const {
    return util::LabelFactory::pt3RelCut(par_->ptSoftMax(ptSoftBinIdx));
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::jets() const {
    TString id = "";
    if( par_->jetAlgo() == JetProperties::AK5 ) id += "AK5";
    
    if( par_->jetType() == JetProperties::Calo ) id += "Calo";
    else if( par_->jetType() == JetProperties::PF ) id += "PF";
    else if( par_->jetType() == JetProperties::JPT ) id += "JPT";
    

    return util::LabelFactory::labelJet(id);
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::pt() const { 
    return util::LabelFactory::ptAve();
  }

  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptSoft() const { 
    return util::LabelFactory::pt3Rel();
  }
  TString PlotMaker::LabelMaker::ptSoftMax() const { 
    return util::LabelFactory::pt3RelMax();
  }

  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::ptSoftGen() const { 
    return util::LabelFactory::pt3RelGen();
  }
  TString PlotMaker::LabelMaker::ptSoftGenMax() const { 
    return util::LabelFactory::pt3RelGenMax();
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::lumi() const {
    return util::StyleSettings::luminosity(par_->lumi());
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::type(FitResult::Type type) const { 
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


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::sample(const SampleLabel &label) const {
    TString lab = label;
    if( util::StyleSettings::style() != util::StyleSettings::PAS &&
	Sample::type(label) == Sample::Data )
      lab += " ("+util::StyleSettings::luminosity(par_->lumi())+")";

    return lab;    
  }
}
