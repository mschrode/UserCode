// $Id: PlotMaker.cc,v 1.14 2011/06/08 06:12:15 mschrode Exp $

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
      yMinExtraRes_(1E-3),
      yMaxExtraRes_(util::StyleSettings::style() == util::StyleSettings::Presentation ? 0.47 : 0.38),
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
    }

    //gStyle->SetEndErrorSize(6);

    if( false ) {
      title_ = "CMS preliminary";
//       // For paper-style
//       gStyle->SetTitleAlign(13);
//       gStyle->SetTitleFontSize(0.1);
//       gStyle->SetTitleX(0.7);
//       gStyle->SetTitleH(0.038);

      // For presentation-style
      gStyle->SetTitleAlign(13);
      gStyle->SetTitleFontSize(0.15);
      gStyle->SetTitleX(0.67);
      gStyle->SetTitleH(0.041);
    } else {
      title_ = "";
    }
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
	      util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","events",true);
	      hPtAsym->GetXaxis()->SetRangeUser(-asymMax,asymMax);
	      setStyle(sample,hPtAsym);
	      
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
	    
	      //util::HistOps::setYRange(hPtAsym,label->GetSize()+leg->GetNRows()+1);
	      util::HistOps::setYRange(hPtAsym,label->GetSize()+1);
	      out_->nextMultiPad(sample->label()+": PtAsym "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hPtAsym->Draw("PE1");
	      for(unsigned int i = 0; i < fits.size(); ++i) {
		//fits.at(i)->Draw("same");
	      }
	      label->Draw("same");
	      //leg->Draw("same");
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



    // +++++ Asymmetry plots per bin for different pt3 ++++++++++++++++++++++++++++++

    // Loop over SampleLabels
    for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      
      SampleLabel sampleLabel = sTIt->first;
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
	    std::vector<TString> labels;
	    unsigned int step = (par_->nPtSoftBins() > 3 ? 2 : 1);
	    for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ptSoftBinIdx += step) {
	      hPtAsyms.push_back(sample->histPtAsym(ptSoftBinIdx));
	      labels.push_back(labelMk_->ptSoftRange(ptSoftBinIdx));
	    }

	    // Labels and style
	    TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(labels.size(),-0.35,label->GetSize());

	    for(unsigned int i = 0; i < hPtAsyms.size(); ++i) {
	      util::HistOps::setAxisTitles(hPtAsyms.at(i),"Asymmetry","","events",true);
	      hPtAsyms.at(i)->GetXaxis()->SetRangeUser(-asymMax,asymMax);
	      setStyle(sample,hPtAsyms.at(i));
	      hPtAsyms.at(i)->SetLineColor(util::StyleSettings::color(i));
	      util::HistOps::setYRange(hPtAsyms.at(i),label->GetSize()+1);
	      leg->AddEntry(hPtAsyms.at(i),labels.at(i),"L");
	    }

	    out_->nextMultiPad(sample->label()+": PtAsym "+ptBin->toTString());
	    hPtAsyms.front()->Draw("HIST");
	    for(unsigned int i = 1; i < hPtAsyms.size(); ++i) {
	      hPtAsyms.at(i)->Draw("HISTsame");
	    }
	    label->Draw("same");
	    leg->Draw("same");
	    gPad->RedrawAxis();
	    out_->saveCurrentPad(histFileName("PtAsym",ptBin,sample));
	    
	    for(unsigned int i = 0; i < hPtAsyms.size(); ++i) {
	      delete hPtAsyms.at(i);
	    }
	    delete label;
	    delete leg;
	  } // End of loop over Samples
	} // End of loop over pt bins
      } // End of loop over eta bins
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
		if( hPtAsym1->Integral() ) hPtAsym1->Scale(1./hPtAsym1->Integral("width"));
		util::HistOps::setAxisTitles(hPtAsym1,"Asymmetry","","events",true);
		hPtAsym1->GetXaxis()->SetRangeUser(-asymMax,asymMax);
		setStyle(sample1,hPtAsym1);

		double sigma = sample1->fittedValue(*rIt,ptSoftBinIdx)/sqrt(2.);
		TF1* fit1 = new TF1("AsymmetryFit_Sample1","gaus",-1.,1.);
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
      for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ++ptSoftBinIdx) {
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
	    legNom->AddEntry(hAsymData,labelMk_->label(sLabelData),"P");
	    legNom->AddEntry(hAsymMC,labelMk_->label(sLabelMC)+" (Scaled)","F");
	    TLegend* legSmeared = util::LabelFactory::createLegendColWithOffset(label->GetSize(),-labelMk_->start(),label->GetSize()); 	   
	    legSmeared->AddEntry(hAsymData,labelMk_->label(sLabelData),"P");
	    legSmeared->AddEntry(hAsymMCSmeared,labelMk_->label(sLabelMC)+" Smeared (Scaled)","F");


	    if( par_->outMode() == OutputManager::EPSSingleFiles ) {

	      // Asymmetry, linear
	      out_->newPage("Data vs MC");      
	      out_->nextMainPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));

	      util::HistOps::setYRange(hAsymMC,label->GetSize()+legNom->GetNRows()+1);
	      TH1* hFrameMain = out_->mainFrame(hAsymMC);
	      hFrameMain->SetTitle(title_);	  
	      TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.4);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.4);

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

	      util::HistOps::setYRange(hAsymMC,50,3E-1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hAsymMC);
	      hFrameMain->SetTitle(title_);	  
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,1.);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,1.);

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

	      util::HistOps::setYRange(hIntAsymMC,50,3E-1);
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
	      hFrameMain->GetXaxis()->SetRangeUser(0.,0.4);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,0.4);

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

	      util::HistOps::setYRange(hAsymMCSmeared,80,3E-1);
	      delete hFrameMain;
	      hFrameMain = out_->mainFrame(hAsymMCSmeared);
	      hFrameMain->SetTitle(title_);	  
	      hFrameRatio = out_->ratioFrame(hFrameMain,"Asymmetry","",0.53,1.47);
	      hFrameRatio->SetLineWidth(lineWidth_);
	      hFrameMain->GetXaxis()->SetRangeUser(0.,1.);
	      hFrameRatio->GetXaxis()->SetRangeUser(0.,1.);

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

	      util::HistOps::setYRange(hIntAsymMCSmeared,50,3E-1);
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
	      util::HistOps::setYRange(hAsymMC,50,1E-3);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hAsymMC->Draw("HIST");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legNom->Draw("same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymLog",ptBin,sMC,ptSoftBinIdx));

	      // Cumulative asymmetry
	      hIntAsymMC->GetXaxis()->SetRangeUser(0.,0.6);
	      util::HistOps::setYRange(hIntAsymMC,50,1E-3);
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
	      util::HistOps::setYRange(hAsymMCSmeared,50,1E-3);
	      out_->nextMultiPad("DataVsMC "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	      hAsymMCSmeared->Draw("HIST");
	      hAsymData->Draw("PE1same");
	      label->Draw("same");
	      legSmeared->Draw("same");
	      out_->logy();
	      out_->saveCurrentPad(histFileName("DataVsMC_PtAsymSmearedLog",ptBin,sMC,ptSoftBinIdx));

	      // Cumulative smeared asymmetry
	      hIntAsymMCSmeared->GetXaxis()->SetRangeUser(0.,0.6);
	      util::HistOps::setYRange(hIntAsymMCSmeared,50,1E-3);
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
	    legNom->AddEntry(hPtAveData,labelMk_->label(sLabelData),"P");
	    legNom->AddEntry(hPtAveMC,labelMk_->label(sLabelMC)+" (Scaled)","F");

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
				     1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	      util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"",
					   sample->labelQuantityInExtrapolation(*rIt));
	      hFrame->GetYaxis()->SetRangeUser(0.8*val.front(),1.5*val.back());
	       
	      // Labels
	      TPaveText* label = labelMk_->ptBin(sample->label(),etaBin->etaBin(),ptBin->ptBin());
	      TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-labelMk_->start(),label->GetSize());
	      //leg->AddEntry(gVals,labelMk_->label(*rIt),"P");
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
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,-labelMk_->start(),label->GetSize());

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
	      gVals1 = new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
					     &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
					     &(uncert.front()),&(uncert.front()));
	      setStyle(sample1,gVals1);

	      val.clear();
	      uncert.clear();
	      sample2->valuesInExtrapolation(*rIt,val,uncert);
	      ptSoftx.clear();
	      sample2->ptSoft(ptSoftx);
	      ptSoftxUncert.clear();
	      ptSoftxUncert = std::vector<double>(ptSoftx.size(),0.);
	      gVals2 = new TGraphAsymmErrors(ptSoftx.size(),&(ptSoftx.front()),&(val.front()),
					     &(ptSoftxUncert.front()),&(ptSoftxUncert.front()),
					     &(uncert.front()),&(uncert.front()));
	      setStyle(sample2,gVals2);

	      // Extrapolation functions
	      fit1 = sample1->extrapolationFunction(*rIt,"Extrapolation_"+sLabel1+"_FitRange");
	      fit1->SetLineColor(gVals1->GetMarkerColor());
	      fit1->SetLineWidth(lineWidth_);
	      extra1 = static_cast<TF1*>(fit1->Clone("Extrapolation_"+sLabel1+"_PlotRange"));
	      extra1->SetRange(0.,1.4*ptSoftx.back());
	      extra1->SetLineStyle(2);

	      fit2 = sample2->extrapolationFunction(*rIt,"Extrapolation_"+sLabel2+"_FitRange");
	      fit2->SetLineColor(gVals2->GetMarkerColor());
	      fit2->SetLineWidth(lineWidth_);
	      extra2 = static_cast<TF1*>(fit2->Clone("Extrapolation_"+sLabel2+"_PlotRange"));
	      extra2->SetRange(0.,1.4*ptSoftx.back());
	      extra2->SetLineStyle(2);

	      leg->AddEntry(gVals1,labelMk_->label(sLabel1),"PL");
	      leg->AddEntry(gVals2,labelMk_->label(sLabel2),"PL");

	      yAxisLabel = sample1->labelQuantityInExtrapolation(*rIt);
	    } // End of loop over Samples
	   
	      // Create frame
	    TH1* hFrame = new TH1D("PlotMaker::plotExtrapolation::hFrame",title_,
				   1000,0.,1.4*par_->ptSoftMax(par_->nPtSoftBins()-1));
	    util::HistOps::setAxisTitles(hFrame,labelMk_->ptSoft(),"",yAxisLabel);
	    hFrame->GetYaxis()->SetRangeUser(0.8*std::min(gVals1->GetY()[0],gVals2->GetY()[1]),
					     1.5*std::max(gVals1->GetY()[gVals1->GetN()-1],gVals2->GetY()[gVals2->GetN()-1]));
	     
	    // Linear scale
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
	sample->valuesInExtrapolation(FitResult::PtGenAsym,val,uncert);
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
	TH1* hFrame = new TH1D("PlotMaker::plotParticleLevelImbalance::hFrame",title_,
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
	    setStyle(sample,hPtGen);

	    TH1* hPdf = sample->histPdfPtTrue(ptSoftBinIdx);
	    hPdf->SetLineColor(kRed);
	    hPdf->SetLineWidth(lineWidth_);
	    
	    // Labels
	    TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,labelMk_->start(),label->GetSize());
	    leg->AddEntry(hPtGen,"Generated Spectrum","P");
	    leg->AddEntry(hPdf,"Assumed Spectrum f","L");
	    
	    // Linear scale
	    util::HistOps::setYRange(hPtGen,label->GetSize()+leg->GetNRows()+1);
	    out_->nextMultiPad(sample->label()+": PtGen Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
	    hPtGen->Draw("PE1");
	    hPdf->Draw("Lsame");
	    label->Draw("same");
	    leg->Draw("same");
	    gPad->RedrawAxis();
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
		  
		  if( ptSoftBinIdx == 2 && sample->label() == "Data" ) {
		    std::cout << "**** " << ptBin->toTString() << ": " << h->GetMean() << " \\pm " << h->GetMeanError() << std::endl;
		  }
		}
		//h->GetXaxis()->SetRangeUser(-asymMax,asymMax);
		setStyle(sample,h);
	      
		// Labels
		TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
		util::HistOps::setYRange(h,label->GetSize()+1);
		out_->nextMultiPad(sample->label()+": Spectrum "+ptBin->toTString()+", PtSoftBin "+util::toTString(ptSoftBinIdx));
		h->Draw("HIST");//h->Draw("PE1");
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
	  setStyle(sampleLabel,gRes);
	  gRes->SetMarkerStyle(27);

	  TGraphAsymmErrors* gResCorr = etaBin->correctedResolution(sampleLabel,fitResType);
	  setStyle(sampleLabel,gResCorr);

	  // MC truth resolution function
	  TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBinIdx));
	  mcTruth->SetLineColor(kRed);
	  mcTruth->SetLineWidth(lineWidth_);

	  // Ratio of measurement and mc truth
	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gResCorr,mcTruth);

	  // PLI function
	  TF1* pli = etaBin->pliFunc("PLI_Eta"+util::toTString(etaBinIdx));
	  pli->SetRange(xMinPt_,xMaxPt_);
	  pli->SetLineColor(kGreen);
	  pli->SetLineStyle(2);
	  pli->SetLineWidth(lineWidth_);

	  // Create frames for main and ratio plots
	  TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	  hFrameMain->SetTitle(title_);
	  hFrameMain->GetXaxis()->SetMoreLogLabels();
	  hFrameMain->GetXaxis()->SetNoExponent();

	  TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	  hFrameRatio->GetXaxis()->SetMoreLogLabels();
	  hFrameRatio->GetXaxis()->SetNoExponent();
	  hFrameRatio->SetLineWidth(2);
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(sampleLabel,etaBin->etaBin());
	  TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,labelMk_->start(),label->GetSize());
	  leg->AddEntry(gRes,"Resolution (Extrapolated)","P");
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
	  hFrameRatio->Draw("][");
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
	  setStyle(sLabel1,gRes1);

	  TGraphAsymmErrors* gRes2 = etaBin->correctedResolution(sLabel2,fitResType);
	  setStyle(sLabel2,gRes2);

	  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gRes1,gRes2);

	  // If data / MC comparison for scaling factor is shown,
	  // plot also MC truth
	  // MC truth resolution function
	  bool showMCTruth = etaBin->hasKValue(sLabel1,sLabel2,fitResType);
	  TF1* mcTruth = etaBin->mcTruthResoFunc("MCTruthReso_Eta"+util::toTString(etaBinIdx));
	  mcTruth->SetLineColor(kRed);
	  mcTruth->SetLineWidth(lineWidth_);

	  // Create frames for main and ratio plots
	  TH1* hFrameMain = out_->mainFrame(xMinPt_,xMaxPt_,yMinExtraRes_,yMaxExtraRes_,"#sigma / p^{ref}_{T}");
	  hFrameMain->SetTitle(title_);
	  hFrameMain->GetXaxis()->SetMoreLogLabels();
	  hFrameMain->GetXaxis()->SetNoExponent();

	  TH1* hFrameRatio = out_->ratioFrame(hFrameMain,"p^{ref}_{T}","GeV",yMinResRatio_,yMaxResRatio_);
	  hFrameRatio->GetXaxis()->SetMoreLogLabels();
	  hFrameRatio->GetXaxis()->SetNoExponent();
	  hFrameRatio->SetLineWidth(lineWidth_);
	     
	  // Labels
	  TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	  TLegend* leg = 0;
	  if( showMCTruth ) 
	    leg = util::LabelFactory::createLegendColWithOffset(3,-1.*(labelMk_->start()),label->GetSize());
	  else 
	    leg = util::LabelFactory::createLegendColWithOffset(2,-1.*(labelMk_->start()),label->GetSize());
	  leg->AddEntry(gRes1,labelMk_->label(sLabel1),"P");
	  leg->AddEntry(gRes2,labelMk_->label(sLabel2),"P");
	  if( showMCTruth ) 
	    leg->AddEntry(mcTruth,"MC Truth","L");

	  out_->nextMainPad(sLabel1+"vs"+sLabel2+": Resolution "+etaBin->toString());
	  hFrameMain->Draw();
	  mcTruth->Draw("same");
	  gRes1->Draw("PE1same");
	  gRes2->Draw("PE1same");
	  label->Draw("same");
	  leg->Draw("same");
	  out_->nextRatioPad();
	  hFrameRatio->Draw("][");
	  gRatio->Draw("PE1same");
	  out_->logx();
	  out_->saveCurrentPad(histFileName("Resolution",etaBin,sLabel1,sLabel2,fitResType));

	  delete gRes1;
	  delete gRes2;
	  delete gRatio;
	  delete mcTruth;
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
	    setStyle(sLabel1,gRatio);
	    TF1* kValueLine = etaBin->kValueLine(sLabel1,sLabel2,fitResType,"kValueLine",xMinPt_,xMaxPt_);
	    kValueLine->SetLineWidth(lineWidth_);
	    TGraphAsymmErrors* kStatBand = etaBin->kValueStatBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    kStatBand->SetLineWidth(lineWidth_);
	    TGraphAsymmErrors* kSystBand = etaBin->kValueSystBand(sLabel1,sLabel2,fitResType,xMinPt_,xMaxPt_);
	    
	    // Create frame
	    TH1* hFrame = util::HistOps::createRatioFrame(xMinPt_,xMaxPt_,0.71,1.89,"p^{ref}_{T} (GeV)","#sigma("+sLabel1+") / #sigma("+sLabel2+")");
	    hFrame->SetTitle(title_);
	    hFrame->GetXaxis()->SetMoreLogLabels();
	    hFrame->GetXaxis()->SetNoExponent();
	    hFrame->SetLineWidth(lineWidth_);
	     
	    // Labels
	    TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	    //TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,-std::min(0.8,labelMk_->start()),label->GetSize());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,-1.,label->GetSize());
	    leg->AddEntry(gRatio,"Measurement","P");
	    leg->AddEntry(kValueLine,"Fit","L");
	    leg->AddEntry(kStatBand,"Stat. Uncertainty","F");
	    leg->AddEntry(kSystBand,"Syst. Uncertainty (Incomplete)","F");

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
	    TPaveText* label = labelMk_->etaBin(etaBin->etaBin());
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,std::min(0.85,labelMk_->start()),label->GetSize());
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
	    hFrameRatio->SetLineWidth(lineWidth_);
	     
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
	    TLegend* leg = util::LabelFactory::createLegendColWithOffset(uncert->nComponents()+1,-0.88*labelMk_->start(),label->GetSize());
	    // Add total uncertainty
	    bands.push_back(uncert->relUncertSteps());
	    for(SystUncertIt it = uncert->componentsBegin(); it != uncert->componentsEnd(); ++it) {
	      // Add components
	      bands.push_back((*it)->relUncertSteps());
	      leg->AddEntry(bands.back(),(*it)->label(),"F");
	    }
	    leg->AddEntry(bands.front(),uncert->label(),"F");

	    // Create frame
	    TH1* hFrame = new TH1D("hFrame",title_,1000,xMinPt_,xMaxPt_);
	    hFrame->GetXaxis()->SetMoreLogLabels();
	    hFrame->GetXaxis()->SetNoExponent();
	    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
	      hFrame->SetBinContent(bin,0.);
	    }
	    hFrame->SetLineWidth(lineWidth_);
	    hFrame->SetLineStyle(2);
	    hFrame->GetYaxis()->SetRangeUser(-0.32,0.96);
	    util::HistOps::setAxisTitles(hFrame,"p^{ref}_{T}","GeV","Relative Uncertainty");
	     
	    out_->nextPad(sampleLabel+": Systematic uncertainties "+etaBin->toString());
	    hFrame->Draw("][");
	    for(std::vector<TGraphAsymmErrors*>::iterator it = bands.begin();
		it != bands.end(); ++it) {
	      (*it)->Draw("E2same");
	    }
	    hFrame->Draw("][same");
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


  // MC event information:
  // - Event weights;
  // - Number of added PU interactions
  // -------------------------------------------------------------------------------------
  void PlotMaker::plotMCEventInfo() const {
//     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotMCEventInfo(): Entering" << std::endl;
//     std::cout << "Plotting MC event information" << std::endl;

//     // Loop over eta and pt bins
//     for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
//       const EtaBin* etaBin = *etaBinIt;

//       out_->newPage("MCEventWeights");

//       for(PtBinIt ptBinIt = (*etaBinIt)->ptBinsBegin(); ptBinIt != (*etaBinIt)->ptBinsEnd(); ++ptBinIt) {
// 	const PtBin* ptBin = *ptBinIt;
// 	if( par_->verbosity() > 1 ) std::cout << "  " << ptBin->toTString() << std::endl;
	
// 	// Loop over MCSamples in this bin
// 	for(MCSampleIt sIt = ptBin->mcSamplesBegin(); sIt != ptBin->mcSamplesEnd(); ++sIt) {
// 	  const MCSample* sample = sIt->second;
// 	  if( par_->verbosity() > 1 ) std::cout << "    " << sample->label() << std::endl;
	  
// 	  // Weight distribution for different ptSoft
// 	  std::vector<TH1*> hWeights;
// 	  std::vector<TString> labels;
// 	  unsigned int step = (par_->nPtSoftBins() > 3 ? 2 : 1);
// 	  for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < par_->nPtSoftBins(); ptSoftBinIdx += step) {
// 	    hWeights.push_back(sample->histPtAsym(ptSoftBinIdx));
// 	    labels.push_back(labelMk_->ptSoftRange(ptSoftBinIdx));
// 	  }


// 	  for(unsigned int ptSoftBinIdx = 0; ptSoftBinIdx < sample->nPtSoftBins(); ++ptSoftBinIdx) {
// 	    if( par_->verbosity() > 1 ) std::cout << "      PtSoftBin " << ptSoftBinIdx << std::endl;
	    
// 	    // Get ptGen spectrum and fit from sample and tweak style
// 	    TH1* hPtGen = sample->histPtGen(ptSoftBinIdx);



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
// 	    else if( distIdx == 1 ) out_->newPage("Pt2Spectrum");
// 	    else if( distIdx == 1 ) out_->newPage("Pt3Spectrum");
	    
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
// 		if( distIdx == 0 ) {
// 		  h = sample->histPtAve(ptSoftBinIdx);
// 		  util::HistOps::setAxisTitles(h,"p^{ave}_{T}","GeV","events");
// 		} else if( distIdx == 1 ) {
// 		  h = sample->histPt1(ptSoftBinIdx);
// 		  util::HistOps::setAxisTitles(h,"p_{T,1}","GeV","events");
// 		} else if( distIdx == 2 ) {
// 		  h = sample->histPt2(ptSoftBinIdx);
// 		  util::HistOps::setAxisTitles(h,"p_{T,2}","GeV","events");
// 		} else if( distIdx == 3 ) {
// 		  h = sample->histPt3(ptSoftBinIdx);
// 		  util::HistOps::setAxisTitles(h,"p_{T,3}","GeV","events");
		  
// 		  if( ptSoftBinIdx == 2 && sample->label() == "Data" ) {
// 		    std::cout << "**** " << ptBin->toTString() << ": " << h->GetMean() << " \\pm " << h->GetMeanError() << std::endl;
// 		  }
// 		}
// 		//h->GetXaxis()->SetRangeUser(-asymMax,asymMax);
// 		setStyle(sample,h);
	      
// 		// Labels
// 		TPaveText* label = labelMk_->ptSoftBin(sample->label(),etaBin->etaBin(),ptBin->ptBin(),ptSoftBinIdx);
	    
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

//     if( par_->verbosity() > 1 ) std::cout << "PlotMaker::plotPtSpectra(): Leaving" << std::endl;
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const EtaBin* etaBin, SampleLabel label1, SampleLabel label2, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(etaBin->etaBin())+".eps");
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const EtaBin* etaBin, SampleLabel sampleLabel, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sampleLabel+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(etaBin->etaBin())+".eps");
  }



  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sample->label()+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+".eps");
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, FitResult::Type type) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_AllSamples_"+FitResult::toString(type)+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+".eps");
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sample->label()+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+".eps");
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, const Sample* sample, unsigned int ptSoftBinIdx) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+sample->label()+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx)+".eps");
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::histFileName(const TString &id, const PtBin* ptBin, SampleLabel label1, SampleLabel label2, FitResult::Type type, unsigned int ptSoftBinIdx) const {
    return cleanFileName(par_->outFilePrefix()+"_"+id+"_"+label1+"_vs_"+label2+"_"+FitResult::toString(type)+"_EtaBin"+util::toTString(ptBin->etaBin())+"_PtBin"+util::toTString(ptBin->ptBin())+"_PtSoftBin"+util::toTString(ptSoftBinIdx)+".eps");
  }


  // -------------------------------------------------------------------------------------
  TString PlotMaker::cleanFileName(TString str) const {
    str.ReplaceAll(" ","");
    str.ReplaceAll("(","");
    str.ReplaceAll(")","");
    str.ReplaceAll("#","");

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
    if( Sample::type(label) == Sample::Data ) {
      lab->AddText(label+",  #sqrt{s} = 7 TeV, L = "+util::toTString(par_->lumi())+" pb^{-1}");
    } else {
      lab->AddText(label+",  #sqrt{s} = 7 TeV");
    }
    lab->AddText(jets()+",  "+etaRange(etaBinIdx));

    return lab;
  }


  // -------------------------------------------------------------------------------------
  TPaveText* PlotMaker::LabelMaker::etaBin(unsigned int etaBinIdx, unsigned int nExtraEntries) const {
    TPaveText* lab = util::LabelFactory::createPaveText(1+nExtraEntries);
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
    return "p_{T,3} / p^{ave}_{T} < "+util::toTString(par_->ptSoftMax(ptSoftBinIdx));
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
    return "p^{rel}_{T,3} = p_{T,3} / p^{ave}_{T}";
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


  // -------------------------------------------------------------------------------------
  TString PlotMaker::LabelMaker::label(const SampleLabel &label) const {
    TString lab = label;
    if( Sample::type(label) == Sample::Data )
      lab += " (L = "+util::toTString(par_->lumi())+" pb^{-1})";

    return lab;    
  }
}
