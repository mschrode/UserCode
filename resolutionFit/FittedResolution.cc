#include "FittedResolution.h"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h" 
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TString.h"

#include "KalibriFileParser.h"
#include "ResponseFunction.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"
#include "../util/utils.h"


namespace resolutionFit {
  // -------------------------------------------------------------------------------------
  FittedResolution::FittedResolution(const std::vector<PtBin*> &ptBins, const Parameters *par) 
    : par_(par), ptBins_(ptBins) {

    lineWidth_ = util::StyleSettings::lineWidth();
    lineHeight_ = util::LabelFactory::lineHeight();

    if( par_->hasMCTruthBins() ) std::cout << "Adding pseudo gamma+jet measurement from MC truth" << std::endl;

    ptMin_ = 0.8*meanPt(0);
    if( par_->hasMCTruthBins() ) ptMin_ = 0.8*par_->mcTruthPtMin(0);
    ptMax_ = 1.1*meanPt(nPtBins()-1);

    // MC truth resolution
//       trueRes_ = new TF1("FittedResolution:trueRes",
//   		       "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
//   		       ptMin_,ptMax_);
//       for(int i = 0; i < 3; i++) {
//         trueRes_->SetParameter(i,par_->trueGaussResPar(i));
//       }


     ////////////// HACK //////////////
    trueRes_ = new TF1("FittedResolution:trueRes","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",ptMin_,ptMax_);
    for(int i = 0; i < 4; i++) {
      trueRes_->SetParameter(i,par_->trueGaussResPar(i));
    }


    // PF Eta 0 - 1.1
//     trueRes_->SetParameter(0,-1.18591);
//     trueRes_->SetParameter(1,0.40573);
//     trueRes_->FixParameter(2,0.);
//     trueRes_->SetParameter(3,0.352069);

// PF Eta 1.1 - 1.7
//     trueRes_->SetParameter(0,-1.49592);
//     trueRes_->SetParameter(1,0.627368);
//     trueRes_->SetParameter(2,0);
//     trueRes_->SetParameter(3,0.218365);

// // PF eta 1.7 - 2.3
// trueRes_->SetParameter(0,-1.51996);
// trueRes_->SetParameter(1,0.766658);
// trueRes_->SetParameter(2,0);
// trueRes_->SetParameter(3,0.0228755);

// // PF Eta 2.3 - 5.0
// trueRes_->SetParameter(0,-0.336561);
// trueRes_->SetParameter(1,0.572859);
// trueRes_->SetParameter(2,0);
// trueRes_->SetParameter(3,0.144683);




//      // Calo 30 - 2000 GeV
//      trueRes_->SetParameter(0,3.8663);
//      trueRes_->SetParameter(1,0.728714);
//      trueRes_->FixParameter(2,0.);
//      trueRes_->SetParameter(3,0.224013);
    
    ///////////////////////

    trueRes_->SetLineWidth(lineWidth_);
    trueRes_->SetLineColor(4);
    trueRes_->SetLineStyle(1);

    // Correction from ptGen asymmetry
    ptGenAsym_ = new TF1("FittedResolution:ptGenAsym",
		       "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		       ptMin_,ptMax_);
    if( par_->hasCorrPtGenAsym() ) {
      std::cout << "Correction of ptGen asymmetry " << std::flush;
      for(int i = 0; i < 3; ++i) {
	ptGenAsym_->SetParameter(i,par_->trueGaussResPar(i));
      }
      if( par_->fitPtGenAsym() ) {
	std::cout << "from fit of distributions" << std::endl;
	TGraphAsymmErrors *g = getTGraphOfResolution("PtGenAsym","Statistic");
	g->Fit(ptGenAsym_,"0QR");
	delete g;
      } else {
	std::cout << "with specified function" << std::endl;
	for(int i = 0; i < 3; i++) {
	  ptGenAsym_->SetParameter(i,par_->corrPtGenAsymPar(i));
	}
      }
    }
    ptGenAsym_->SetLineWidth(lineWidth_);
    ptGenAsym_->SetLineColor(6);
    ptGenAsym_->SetLineStyle(2);

    // Fit extrapolated resolution using
    // statistic uncertainty
    TGraphAsymmErrors *gMaxLike = getTGraphOfResolution("MaxLike","Statistic",true);
    fittedRes_ = new TF1("FittedResolution::fittedRes_","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
			 ptMin_,ptMax_);
    TGraphAsymmErrors *gPtAsym = getTGraphOfResolution("PtAsym","Statistic",true);
    fittedResAsym_ = new TF1("FittedResolution::fittedResAsym_","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
			     ptMin_,ptMax_);
    for(int i = 0; i < 3; ++i) {
      fittedRes_->SetParameter(i,par_->trueGaussResPar(i));
      fittedResAsym_->SetParameter(i,par_->trueGaussResPar(i));
    }
    fittedRes_->SetLineColor(2);
    fittedRes_->SetLineWidth(lineWidth_);
    fittedResAsym_->SetLineColor(2);
    fittedResAsym_->SetLineWidth(lineWidth_);
    if( par_->fitExtrapolatedSigma() ) {
      gMaxLike->Fit(fittedRes_,"0QR");
      gPtAsym->Fit(fittedResAsym_,"0QR");
    }
    delete gMaxLike;
    delete gPtAsym;

    if( par_->hasMCClosure() ) fitMCClosure();
  }

  // -------------------------------------------------------------------------------------
  FittedResolution::~FittedResolution() {
    if( par_->verbosity() > 1 ) std::cout << "FittedResolution::~FittedResolution(): Entering\n";
    delete trueRes_;
    delete ptGenAsym_;
    delete fittedRes_;
    delete fittedResAsym_;
    for(int i = 0; i < nMCClosureResFits(); ++i) {
      delete mcClosureGReso_[i];
      delete mcClosureGScale_[i];
      for(int p = 0; p < nPtBins(); ++p) {
	delete mcClosureFits_[i][p];
      }
    }
    if( par_->verbosity() > 1 ) std::cout << "FittedResolution::~FittedResolution(): Leaving\n";
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotPtAsymmetry() const {
    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotPtAsymmetry(): Entering method" << std::endl;
    }

    // Loop over pt3 cuts
    for(int c = 0; c < par_->nPt3Cuts(); ++c) {
      std::vector<double> ptMeanMaxLike(nPtBins());
      std::vector<double> ptMeanMaxLikeErr(nPtBins());
      std::vector<double> ptMeanPtAsym(nPtBins());
      std::vector<double> ptMeanPtAsymErr(nPtBins());
      std::vector<double> resMaxLike(nPtBins());
      std::vector<double> resMaxLikeErr(nPtBins());
      std::vector<double> resPtAsym(nPtBins());
      std::vector<double> resPtAsymErr(nPtBins());
      std::vector<double> resRatio(nPtBins());
      std::vector<double> resRatioErr(nPtBins());

      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	  it != ptBins_.end(); it++) {
	int ptBin = (it-ptBins_.begin());

	// Store fitted widths
	ptMeanMaxLike[ptBin] = (*it)->meanPt();
	ptMeanMaxLikeErr[ptBin] = (*it)->meanPtUncert();
	ptMeanPtAsym[ptBin] = (*it)->meanPtAsym();
	ptMeanPtAsymErr[ptBin] = (*it)->meanPtAsymUncert();
	resMaxLike[ptBin] = (*it)->fittedValue(0,c)/sqrt(2.);
	resMaxLikeErr[ptBin] = (*it)->fittedValueUncert(0,c)/sqrt(2.);
	resPtAsym[ptBin] = (*it)->fittedAsym(c)/sqrt(2.);
	resPtAsymErr[ptBin] = (*it)->fittedAsymUncert(c)/sqrt(2.);
	resRatio[ptBin] = resMaxLike[ptBin]/resPtAsym[ptBin];
	resRatioErr[ptBin] = sqrt(pow(resMaxLikeErr[ptBin]/resPtAsym[ptBin],2) + pow(resPtAsymErr[ptBin]*resMaxLike[ptBin]/resPtAsym[ptBin]/resPtAsym[ptBin],2));

	double xMin = -0.4;
	double xMax = 0.4;

	// Pt asymmetry distribution
	TH1 *hPtAsym = (*it)->getHistPtAsym(c,"PtAsymmetryBins:hPtAsym");
	util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","events",true);
	hPtAsym->GetXaxis()->SetRangeUser(xMin,xMax);
	util::HistOps::setYRange(hPtAsym,5);
	hPtAsym->SetMarkerStyle(20);

	// Fitted pt asymmetry
	TF1 *fPtAsym = new TF1("PtAsymmetryBins:fPtAsym","gaus",hPtAsym->GetMean()-2.*hPtAsym->GetRMS(),hPtAsym->GetMean()+2.*hPtAsym->GetRMS());
	fPtAsym->SetParameter(1,0.);
	fPtAsym->SetParameter(2,resPtAsym[ptBin]);
	fPtAsym->SetParameter(0,1./sqrt(2.*M_PI)/fPtAsym->GetParameter(2));	
	fPtAsym->SetLineWidth(lineWidth_);
	fPtAsym->SetLineStyle(2);
	fPtAsym->SetLineColor(4);
	
	// Pt asymmetry from MaxLike
	TF1 *fMaxLike = new TF1("PtAsymmetryBins:fMaxLike","gaus",-1.,1);
	fMaxLike->SetParameter(1,0.);
	fMaxLike->SetParameter(2,resMaxLike[ptBin]);
	fMaxLike->SetParameter(0,1./sqrt(2.*M_PI)/fMaxLike->GetParameter(2));	
	fMaxLike->SetLineWidth(lineWidth_);
	fMaxLike->SetLineColor(2);

	// Labels
	TPaveText *label = 0;
	TLegend *leg = 0;
	if( util::StyleSettings::style() == util::StyleSettings::Presentation ) {
	  label = util::LabelFactory::createPaveText(1);
	  label->AddText(par_->labelPtBin(ptBin)+",  "+par_->labelPt3Cut(c));
	  leg = util::LabelFactory::createLegendWithOffset(3,1);
	  leg->SetTextSize(0.04);
	} else {
	  label = util::LabelFactory::createPaveText(4,-0.4);	  
	  label->AddText(par_->labelLumi()+",  "+par_->labelEtaBin());	
	  label->AddText(par_->labelPtSoftCut());
	  label->AddText(par_->labelPt3Cut(c));
	  label->AddText(par_->labelPtBin(ptBin));
	  leg = util::LabelFactory::createLegendCol(3,0.5);
	}
	if( par_->isData() ) leg->AddEntry(hPtAsym,"Data","P");
	else leg->AddEntry(hPtAsym,"Measurement","P");
	leg->AddEntry(fPtAsym,"Binned Fit","L");
	leg->AddEntry(fMaxLike,"Likelihood Fit","L");

	// Draw
	TCanvas *can = new TCanvas("PtAsymmetry","PtAsymmetry ("+util::toTString(ptBin)+" "+util::toTString(c)+")",500,500);
	can->cd();
	hPtAsym->Draw("PE1");
	fPtAsym->Draw("same");
	fMaxLike->Draw("same");
	label->Draw("same");
	leg->Draw("same");
	can->SaveAs(par_->outNamePrefix()+"PtAsymmetry_PtBin"+util::toTString(ptBin)+"_Pt3Cut"+util::toTString(c)+".eps","eps");

	can->cd();
	hPtAsym->GetXaxis()->SetRangeUser(-1.,1.);
	util::HistOps::setYRange(hPtAsym,3,3E-3);
	hPtAsym->Draw("PE1");
	fPtAsym->Draw("same");
	fMaxLike->Draw("same");
	label->Draw("same");
	leg->Draw("same");
	can->SetLogy();
	can->SaveAs(par_->outNamePrefix()+"PtAsymmetryLog_PtBin"+util::toTString(ptBin)+"_Pt3Cut"+util::toTString(c)+".eps","eps");


	// Scan for startpoint of tails
	hPtAsym->Fit("gaus","I0QR","",hPtAsym->GetMean()-2.*hPtAsym->GetRMS(),hPtAsym->GetMean()+2.*hPtAsym->GetRMS());
	TF1 *fit = hPtAsym->GetFunction("gaus");
	fit->SetRange(-1.,1.);
	fit->SetLineWidth(2);
	fit->SetLineColor(2);
	double mean = fit->GetParameter(1);
	double sigma = fit->GetParameter(2);
	double min = 0.;
	double max = 0.;
	util::HistOps::findYRange(hPtAsym,min,max);
	min = 3E-3;
	TLine *lFitRangeLeft = new TLine(hPtAsym->GetMean()-2.*hPtAsym->GetRMS(),min,
					 hPtAsym->GetMean()-2.*hPtAsym->GetRMS(),max);
	lFitRangeLeft->SetLineWidth(1);
	lFitRangeLeft->SetLineColor(4);
	TLine *lFitRangeRight = new TLine(hPtAsym->GetMean()+2.*hPtAsym->GetRMS(),min,
					  hPtAsym->GetMean()+2.*hPtAsym->GetRMS(),max);
	lFitRangeRight->SetLineWidth(1);
	lFitRangeRight->SetLineColor(4);
	std::vector<TLine*> lines;
	for(int k = -6; k < 6; ++k) {
	  if( k == 0 ) continue;
	  double x = mean - (2.-k)*sigma;
	  if( k > 0 ) x = mean + (2.+k)*sigma;
	  if( std::abs(x) > 1. ) continue;
	  TLine *line = new TLine(x,min,x,max);
	  line->SetLineStyle(2);
	  line->SetLineWidth(1);
	  line->SetLineColor(4);
	  lines.push_back(line);
	}
	hPtAsym->Draw("PE1");
	fit->Draw("same");
	lFitRangeLeft->Draw("same");
	lFitRangeRight->Draw("same");
	for(size_t k = 0; k < lines.size(); ++k) {
	  lines[k]->Draw("same");
	}
	leg->Draw("same");
	can->SetLogy();
	can->SaveAs(par_->outNamePrefix()+"PtAsymmetryTails_PtBin"+util::toTString(ptBin)+"_Pt3Cut"+util::toTString(c)+".eps","eps");

	// Clean up
	delete hPtAsym;
	delete fMaxLike;
	delete fPtAsym;
	delete lFitRangeLeft;
	delete lFitRangeRight;
	for(size_t k = 0; k < lines.size(); ++k) {
	  delete lines[k];
	}
	delete label;
	delete leg;
	delete can;
      } // End of loop over pt bins


      // Compare fitted widths
      TGraphErrors *gResMaxLike = new TGraphErrors(ptMeanMaxLike.size(),&(ptMeanMaxLike.front()),&(resMaxLike.front()),&(ptMeanMaxLikeErr.front()),&(resMaxLikeErr.front()));
      gResMaxLike->SetMarkerStyle(20);
      gResMaxLike->SetMarkerColor(2);
      gResMaxLike->SetLineColor(2);
      TGraphErrors *gResPtAsym = new TGraphErrors(ptMeanPtAsym.size(),&(ptMeanPtAsym.front()),&(resPtAsym.front()),&(ptMeanPtAsymErr.front()),&(resPtAsymErr.front()));
      gResPtAsym->SetMarkerStyle(25);
      gResPtAsym->SetMarkerColor(4);
      gResPtAsym->SetLineColor(4);
      TGraphErrors *gRatio = new TGraphErrors(ptMeanMaxLike.size(),&(ptMeanMaxLike.front()),&(resRatio.front()),&(ptMeanMaxLikeErr.front()),&(resRatioErr.front()));
      gRatio->SetMarkerStyle(20);
      gRatio->SetMarkerColor(2);
      gRatio->SetLineColor(gRatio->GetMarkerColor());

      TPaveText *label = 0;
      TLegend *leg = 0;
      if( util::StyleSettings::style() == util::StyleSettings::Presentation ) {
	label = util::LabelFactory::createPaveText(1);
	label->AddText(par_->labelPtSoftCut()+",  "+par_->labelPt3Cut(c));
	leg = util::LabelFactory::createLegendWithOffset(2,1);
	leg->SetTextSize(0.045);
      } else {
	label = util::LabelFactory::createPaveText(3,-0.4);
	label->AddText(par_->labelLumi()+",  "+par_->labelEtaBin());
	label->AddText(par_->labelPt3Cut(c));
	label->AddText(par_->labelPtSoftCut());
	leg = util::LabelFactory::createLegendCol(3,0.5);
      }
      leg->AddEntry(gResPtAsym,"Binned Fit","P");
      leg->AddEntry(gResMaxLike,"Likelihood Fit","P");

      TH1 *hFrame = new TH1D("FittedResolutions",";"+par_->labelPtRef()+" (GeV);#sigma(Asymmetry) / "+par_->labelPtRef(),1000,0.9*ptMin(),1.1*ptMax());
      double yMin = *(std::min_element(gResMaxLike->GetY(),gResMaxLike->GetY()+gResMaxLike->GetN()));
      double yMax = *(std::max_element(gResMaxLike->GetY(),gResMaxLike->GetY()+gResMaxLike->GetN()));
      double tmp =  *(std::min_element(gResPtAsym->GetY(),gResPtAsym->GetY()+gResPtAsym->GetN()));
      if( tmp < yMin ) yMin = tmp;
      tmp =  *(std::max_element(gResPtAsym->GetY(),gResPtAsym->GetY()+gResPtAsym->GetN()));
      if( tmp > yMax ) yMax = tmp;
      hFrame->GetYaxis()->SetRangeUser(0.,2.5*yMax);
      
      TCanvas *can = new TCanvas("PtAsymmetry","PtAsymmetry Width ("+util::toTString(c)+")",500,500);
      can->cd();
      hFrame->Draw();
      gResMaxLike->Draw("PE1same");
      gResPtAsym->Draw("PE1same");      
      label->Draw("same");
      leg->Draw("same");
      can->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"PtAsymmetryWidth_Pt3Cut"+util::toTString(c)+".eps","eps");

      TH1 *hRatioFrame = util::HistOps::createRatioFrame(hFrame,"#sigma(MaxLike) / #sigma(PtAsym)",0.6,1.5);
      hRatioFrame->Draw();
      gRatio->Draw("PE1same");
      label->Draw("same");
      can->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"PtAsymmetryWidthRatio_Pt3Cut"+util::toTString(c)+".eps","eps");
      
      // Ratio at bottom
      TCanvas *bRatioTopCan = util::HistOps::createRatioTopCanvas();
      TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
      TH1 *bRatioTopFrame = util::HistOps::createRatioTopFrame(hFrame);
      TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hFrame,par_->labelPtRef(),"GeV",0.81,1.19);
    
      bRatioTopCan->cd();
      bRatioTopFrame->GetYaxis()->SetRangeUser(1E-3,2.5*yMax);
      bRatioTopFrame->Draw();
      gResPtAsym->Draw("PE1same");
      gResMaxLike->Draw("PE1same");
      label->Draw("same");
      leg->Draw("same");
      gPad->SetLogx();
      bRatioBottomPad->Draw();
      bRatioBottomPad->cd();
      bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
      bRatioBottomFrame->Draw();
      gRatio->Draw("PE1same");
      gPad->SetLogx();
      bRatioTopCan->SaveAs(par_->outNamePrefix()+"PtAsymmetryWidthBottomRatio_Pt3Cut"+util::toTString(c)+".eps","eps");

      delete hFrame;
      delete hRatioFrame;
      delete gResMaxLike;
      delete gResPtAsym;
      delete gRatio;
      delete bRatioTopFrame;
      delete bRatioBottomFrame;
      delete bRatioBottomPad;
      delete bRatioTopCan;
      delete label;
      delete leg;
      delete can;
    } // End of loop over pt3 cuts

    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotPtAsymmetry(): Leaving method" << std::endl;
    }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotExtrapolation() const {
    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotExtrapolation(): Entering method" << std::endl;
    }

    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());

      // Loop over fitted parameters
      for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {

	// Create a frame
	TH1 *h = (*it)->getFrameOfVariation(parIdx,"FrameExtrapolation");
	//	if( !par_->isData() ) {
// 	  for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
// 	    h->SetBinContent(bin,par_->trueGaussSigma((*it)->meanPt()));
// 	  }
//	}
	//	h->GetYaxis()->SetRangeUser(0.8*par_->trueGaussSigma((*it)->meanPt()),2.*par_->trueGaussSigma((*it)->meanPt()));
	h->SetLineStyle(2);

	TH1 *hAsym = (*it)->getFrameOfVariationAsym("FrameExtrapolationAsym");
// 	if( !par_->isData() ) {
// 	  for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
// 	    hAsym->SetBinContent(bin,par_->trueGaussSigma((*it)->meanPt()));
// 	  }
// 	}
// 	hAsym->SetLineStyle(2);

	// Get graph and extrapolation
	TGraphAsymmErrors *g = (*it)->getTGraphOfVariation(parIdx);
	TF1 *f = (*it)->getTF1OfVariation(parIdx,"FitExtrapolation");
	f->SetLineWidth(lineWidth_);

	TGraphAsymmErrors *gAsym = (*it)->getTGraphOfVariationAsym();
	gAsym->SetMarkerStyle(24);
	TF1 *fAsym = (*it)->getTF1OfVariationAsym("FitExtrapolation");
	fAsym->SetLineWidth(lineWidth_);
	fAsym->SetLineStyle(2);
	
	// Draw label
	TPaveText *txt = 0;
	TLegend *leg = 0;
	TLegend *legAsym = 0;
	if( util::StyleSettings::style() == util::StyleSettings::Presentation ) {
	  txt = util::LabelFactory::createPaveText(1);
	  txt->AddText(par_->labelPtBin(ptBin));
	  leg = util::LabelFactory::createLegendWithOffset(2,1);
	  legAsym = util::LabelFactory::createLegendWithOffset(2,1);
	} else {
	  txt = util::LabelFactory::createPaveText(3,-0.4);
	  txt->AddText(par_->labelLumi()+", "+par_->labelEtaBin());
	  txt->AddText(par_->labelPtSoftCut());
	  txt->AddText(par_->labelPtBin(ptBin));
	  leg = util::LabelFactory::createLegendCol(2,0.45);
	  legAsym = util::LabelFactory::createLegendCol(2,0.45);
	}
	leg->AddEntry(g,"Likelihood Fit","P");
	leg->AddEntry(f,"Extrapolation","L");
	legAsym->AddEntry(g,"Binned Fit","P");
	legAsym->AddEntry(f,"Extrapolation","L");

	TLegend *legComp = util::LabelFactory::createLegendCol(4,0.47);
	legComp->AddEntry(g,"Maximum Likelihood Fit","P");
	legComp->AddEntry(f,"Fit to MaxLike","L");
	legComp->AddEntry(gAsym,"Binned Fit","P");
	legComp->AddEntry(fAsym,"Fit to PtAsym","L");

	// Draw
	TCanvas *can = new TCanvas("PlotExtrapolation","Extrapolation ("+util::toTString(ptBin)+")",500,500);
	can->cd();
	h->SetXTitle((par_->fitMode()==FitModeMaxLikeFull) ? "p_{||,3} / <p^{true}_{T}>" : "p_{||,3} / <p^{ave}_{T}>");
	h->SetYTitle((par_->fitMode()==FitModeMaxLikeFull) ? "#sigma / <p^{true}_{T}>" : "#sigma / <p^{ave}_{T}>");
	h->Draw();
	f->Draw("same");
	g->Draw("PE1same");
	txt->Draw("same");
	leg->Draw("same");
	can->SaveAs(par_->outNamePrefix()+"ExtrapolatedPar"+util::toTString(parIdx)+"_PtBin"+util::toTString(ptBin)+".eps","eps");
	if( parIdx == 0 ) {
	  can->cd();
	  hAsym->SetYTitle("#sigma / <p^{ave}_{T}>");
	  hAsym->Draw();
	  gAsym->Draw("PE1same");
	  fAsym->Draw("same");
	  txt->Draw("same");
	  legAsym->Draw("same");
	  can->SaveAs(par_->outNamePrefix()+"ExtrapolatedPar"+util::toTString(parIdx)+"Asym_PtBin"+util::toTString(ptBin)+".eps","eps");

	  can->cd();
	  h->SetYTitle("#sigma / <p^{ref}_{T}>");
	  h->Draw();
	  g->Draw("PE1same");
	  gAsym->Draw("PE1same");
	  f->Draw("same");
	  fAsym->Draw("same");
	  txt->Draw("same");
	  legComp->Draw("same");
	  can->SaveAs(par_->outNamePrefix()+"ExtrapolatedPar"+util::toTString(parIdx)+"CompAsym_PtBin"+util::toTString(ptBin)+".eps","eps");
	}

	// Clean up
	delete h;
	delete hAsym;
	delete f;
	delete fAsym;
	delete g;
	delete gAsym;
	delete txt;
	delete leg;
	delete legAsym;
	delete legComp;
	delete can;
      } // End of loop over fitted parameters



      // Draw pt asymmetry for different pt3 cuts
      TCanvas *can = new TCanvas("PlotVariedAsym","Varied Pt Asymmetry",500,500);
      can->cd();
      std::vector<TH1*> hAsym;
      //TLegend *leg = util::LabelFactory::createLegendCol(3,0.5);
      TLegend *leg = util::LabelFactory::createLegendWithOffset(3,1);
      leg->SetTextSize(0.04);
      for(int i = 0; i < 3; ++i) {
	if( par_->nPt3Cuts() < 3 ) continue;
	int pt3Idx = 0;
	if( i == 1 ) pt3Idx = par_->nPt3Cuts()/2;
	else if( i == 2 ) pt3Idx = par_->nPt3Cuts()-1;

	// Draw pt asymmetry distribution	
	if( par_->verbosity() == 3 ) {
	  std::cout << "FittedResolution::plotExtrapolation(): Drawing pt asymmetry for pt3 bin " << pt3Idx << std::endl;
	}
	TString name = "PlotExtrapolation:VariedPtAsymmetry";
	name += pt3Idx;
	TH1 *h = (*it)->getHistPtAsym(pt3Idx,name);
	h->Rebin();
	util::HistOps::setAxisTitles(h,"Asymmetry","","events",true);
	h->SetMarkerStyle(1);
	h->SetMarkerColor(util::StyleSettings::color(i));
	h->SetLineColor(h->GetMarkerColor());
	h->GetXaxis()->SetRangeUser(-0.3,0.3);
	if( i == 0 ) {
	  util::HistOps::setYRange(h,4);
	  h->Draw("Hist");
	} else {
	  h->Draw("HistSame");
	}
	leg->AddEntry(h,par_->labelPt3Cut(pt3Idx),"L");
	hAsym.push_back(h);
      }
      
      // Draw label
      //TPaveText *txt = util::LabelFactory::createPaveText(2,-0.5);
      TPaveText *txt = util::LabelFactory::createPaveText(1);
      //      txt->AddText(par_->labelLumi()+", "+par_->labelEtaBin());
      txt->AddText(par_->labelPtBin(ptBin));
      txt->Draw("same");
      leg->Draw("same");
      
      // Write canvas to file
      TString name = par_->outNamePrefix();
      name += "VariedPtAsymmetry_PtBin";
      name += ptBin;
      name += ".eps";
      can->SaveAs(name,"eps");

      // Clean up
      for(size_t i = 0; i < hAsym.size(); ++i) {
	delete hAsym[i];
      }
      delete txt;
      delete leg;
      delete can;
    } // End of loop over pt bins

    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotExtrapolation(): Leaving method" << std::endl;
    }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotResolution() const {
    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotResolution(): Entering method" << std::endl;
    }


    // ----- Plot relative resolution sigma / pt -----
    // Create graphs with statistical uncertainties
    TGraphAsymmErrors *gMaxLikeStat = getTGraphOfResolution("MaxLike","Statistic");
    TGraphAsymmErrors *gPtAsymStat = getTGraphOfResolution("PtAsym","Statistic");
    TGraphAsymmErrors *gPtGenAsym = getTGraphOfResolution("PtGenAsym","Statistic");
    TGraphAsymmErrors *gMaxLikeStatCorr = getTGraphOfResolution("MaxLike","Statistic",true);
    TGraphAsymmErrors *gPtAsymStatCorr = getTGraphOfResolution("PtAsym","Statistic",true);

    // Create a frame
    TH1 *h = new TH1D("FrameExtrapolatedResolution","",1000,ptMin_,ptMax_);
    h->GetXaxis()->SetMoreLogLabels();
    double min = 0.5*(*std::min_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
    //double max = 1.6*(*std::max_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
    double max = 2.2*(*std::max_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
    h->GetYaxis()->SetRangeUser(0.,max);

    // Create labels
    int nLegEntries = 2;
    if( par_->hasCorrPtGenAsym() ) nLegEntries += 2;
    if( par_->fitPtGenAsym() ) nLegEntries++;
    if( par_->fitExtrapolatedSigma() ) nLegEntries++;

    TPaveText *txt = 0;
    TLegend *legMaxLike = 0;
    TLegend *legPtAsym = 0;
    if( util::StyleSettings::style() == util::StyleSettings::Presentation ) {
      txt = util::LabelFactory::createPaveText(1,0.8);
      legMaxLike = util::LabelFactory::createLegendColWithOffset(nLegEntries,0.8,1);
      legPtAsym = util::LabelFactory::createLegendColWithOffset(nLegEntries,0.8,1);
    } else {
      txt = util::LabelFactory::createPaveText(1,0.6);
      legMaxLike = util::LabelFactory::createLegendColWithOffset(nLegEntries,0.6,1);
      legPtAsym = util::LabelFactory::createLegendColWithOffset(nLegEntries,0.6,1);
    }
    txt->AddText(par_->labelLumi()+", "+par_->labelEtaBin());

    legMaxLike->AddEntry(gMaxLikeStat,"Likelihood Fit","P");
    if( par_->hasCorrPtGenAsym() ) {
      if( par_->fitPtGenAsym() ) {
	legMaxLike->AddEntry(gPtGenAsym,par_->labelPtGen()+" Asymmetry","P");
	legMaxLike->AddEntry(ptGenAsym_,par_->labelPtGen()+" Asymmetry (Fit)","L");
      } else {
	legMaxLike->AddEntry(ptGenAsym_,par_->labelPtGen()+" Asymmetry","L");
      }
      legMaxLike->AddEntry(gMaxLikeStatCorr,"Likelihood Fit (Corrected)","P");
    }
    if( par_->fitExtrapolatedSigma() ) legMaxLike->AddEntry(fittedRes_,"Interpolated Resolution","L");
    legMaxLike->AddEntry(trueRes_,"MC Truth Resolution","L");

    legPtAsym->AddEntry(gPtAsymStat,"Asymmetry","P");
    if( par_->hasCorrPtGenAsym() ) {
      if( par_->fitPtGenAsym() ) {
	legPtAsym->AddEntry(gPtGenAsym,par_->labelPtGen()+" Asymmetry","P");
	legPtAsym->AddEntry(ptGenAsym_,par_->labelPtGen()+" Asymmetry (Fit)","L");
      } else {
	legPtAsym->AddEntry(ptGenAsym_,par_->labelPtGen()+" Asymmetry","L");
      }
      legPtAsym->AddEntry(gPtAsymStatCorr,"Asymmetry (Corrected)","P");
    }
    if( par_->fitExtrapolatedSigma() ) legPtAsym->AddEntry(fittedRes_,"Interpolated Resolution","L");
    legPtAsym->AddEntry(trueRes_,"MC Truth Resolution","L");

    TLegend *legComp = util::LabelFactory::createLegendCol(3,0.6);
    if( par_->hasCorrPtGenAsym() ) {
      legComp->AddEntry(gMaxLikeStatCorr,"Likelihood Fit (Corrected)","P");
      legComp->AddEntry(gPtAsymStatCorr,"Asymmetry (Corrected)","P");
    } else {
      legComp->AddEntry(gMaxLikeStat,"Likelihood Fit","P");
      legComp->AddEntry(gPtAsymStat,"Asymmetry","P");
    }
    legComp->AddEntry(trueRes_,"MC Truth Resolution","L");
    
    // Create Canvas
    TCanvas *can = new TCanvas("CanExtrapolatedResolution","Extrapolated Resolution",500,500);
    can->cd();

    // Draw MaxLike results
    h->SetXTitle("p^{ref}_{T} (GeV)");
    h->SetYTitle("#sigma / p^{ref}_{T}");
    h->Draw();
    trueRes_->Draw("same");
    gMaxLikeStat->Draw("PE1same");
    if( par_->hasCorrPtGenAsym() ) {
      gMaxLikeStatCorr->Draw("PE1same");
      ptGenAsym_->Draw("same");
      if( par_->fitPtGenAsym() ) gPtGenAsym->Draw("PE1same");
    }
    if( par_->fitExtrapolatedSigma() ) fittedRes_->Draw("same");
    txt->Draw("same");
    legMaxLike->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraReso.eps","eps");

    // Draw PtAsym results
    h->Draw();
    trueRes_->Draw("same");
    gPtAsymStat->Draw("PE1same");
    if( par_->hasCorrPtGenAsym() ) {
      gPtAsymStatCorr->Draw("PE1same");
      ptGenAsym_->Draw("same");
      if( par_->fitPtGenAsym() ) gPtGenAsym->Draw("PE1same");
    }
    if( par_->fitExtrapolatedSigma() ) fittedResAsym_->Draw("same");
    txt->Draw("same");
    legPtAsym->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoPtAsym.eps","eps");

    // Draw comparison
    h->Draw();
    trueRes_->Draw("same");
    if( par_->hasCorrPtGenAsym() ) {
      gMaxLikeStatCorr->Draw("PE1same");
      gPtAsymStatCorr->Draw("PE1same");
    } else {
      gMaxLikeStat->Draw("PE1same");
      gPtAsymStat->Draw("PE1same");
    }
    txt->Draw("same");
    legComp->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoComp.eps","eps");

    if( par_->hasMCClosure() ) {
      TH1 *hMCClosureReso = new TH1D("hMCClosureReso","",nPtBins(),&(par_->ptBinEdges()->front()));
      for(int i = 0; i < nPtBins(); ++i) {
	hMCClosureReso->SetBinContent(1+i,mcClosureReso(0,i));
	hMCClosureReso->SetBinError(1+i,mcClosureResoErr(0,i));
      }
      hMCClosureReso->SetLineColor(2);
      hMCClosureReso->SetLineWidth(1);
      hMCClosureReso->GetYaxis()->SetRangeUser(min,max);
      hMCClosureReso->SetXTitle(h->GetXaxis()->GetTitle());
      hMCClosureReso->SetYTitle(h->GetYaxis()->GetTitle());
      hMCClosureReso->Draw("HISTE");
      trueRes_->Draw("same");
      gMaxLikeStat->Draw("PE1same");
      gPtAsymStat->Draw("PE1same");
      txt->Draw("same");
      TLegend *legCompMCClosure = util::LabelFactory::createLegendWithOffset(nLegEntries+1,txt->GetSize());
      legCompMCClosure->AddEntry(gMaxLikeStat,"MaxLike: Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
      legCompMCClosure->AddEntry(gPtAsymStat,"PtAsym: Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
      legCompMCClosure->AddEntry(trueRes_,"MC truth resolution","L");
      legCompMCClosure->AddEntry(hMCClosureReso,"MC truth distributions","L");
      legCompMCClosure->Draw("same");
      gPad->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"ExtraResoMCClosure.eps","eps");
      delete hMCClosureReso;
      delete legCompMCClosure;
    }


    // ----- Plot relative deviation (sigma(fit)-sigma(true) ) / sigma(true)  -----

    // Create ratio graphs with statistical uncertainties
    TGraphAsymmErrors *gMaxLikeRatioStat = util::HistOps::createRatioGraph(gMaxLikeStat,trueRes_);
    TGraphAsymmErrors *gMaxLikeRatioStatCorr = util::HistOps::createRatioGraph(gMaxLikeStatCorr,trueRes_);
    TGraphAsymmErrors* gMaxLikeRatioMCClosure = getTGraphOfResolution("MaxLike","Statistic");
    if( par_->hasMCClosure() ) {
      for(int i = 0; i < gMaxLikeRatioStat->GetN(); ++i) {
	double x = gMaxLikeStat->GetX()[i];
	double y = gMaxLikeStat->GetY()[i];
	double yTrue = trueRes_->Eval(x);
	double exh = gMaxLikeRatioStat->GetEXhigh()[i];
	double eyh = gMaxLikeRatioStat->GetEYhigh()[i];
	y = y / mcClosureReso(0,i);
	eyh = util::ratioError(y,eyh,mcClosureReso(0,i),mcClosureResoErr(0,i));
	gMaxLikeRatioMCClosure->SetPoint(i,x,y);
	gMaxLikeRatioMCClosure->SetPointError(i,exh,exh,eyh,eyh);
      }	
    }

    TF1 *lineFitRatio = new TF1("LineFitRatioExtraplolatedResolution","pol0",ptMin_,ptMax_);
    lineFitRatio->SetLineWidth(lineWidth_);
    lineFitRatio->SetLineStyle(2);
    lineFitRatio->SetLineColor(2);
    if( par_->fitRatio() ) {
      gMaxLikeRatioStat->Fit(lineFitRatio,"0QR");
      std::cout << "Fitted ratio " << lineFitRatio->GetParameter(0) << std::flush;
      std::cout << " +/- " << lineFitRatio->GetParError(0) << std::endl;
    }
    TF1 *lineStartRes = static_cast<TF1*>(lineFitRatio->Clone("LineStartResExtrapolatedResolution"));
    lineStartRes->SetLineColor(8);
    lineStartRes->SetParameter(0,par_->relStartOffset());

    TGraphAsymmErrors *gPtAsymRatioStat = util::HistOps::createRatioGraph(gPtAsymStat,trueRes_);
    TGraphAsymmErrors *gPtAsymRatioStatCorr = util::HistOps::createRatioGraph(gPtAsymStatCorr,trueRes_);

    TH1 *hRatioFittedRes = util::HistOps::createRatio(fittedRes_,trueRes_,ptMin_,ptMax_,"","");
    hRatioFittedRes->SetLineColor(2);
    TH1 *hRatioFittedResAsym = util::HistOps::createRatio(fittedResAsym_,trueRes_,ptMin_,ptMax_,"","");
    hRatioFittedResAsym->SetLineColor(2);


    // Resolution plots with ratio at bottom
    TCanvas *bRatioTopCan = util::HistOps::createRatioTopCanvas();
    TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
    TH1 *bRatioTopFrame = util::HistOps::createRatioTopFrame(h);
    TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(h,"p^{ref}_{T}","GeV",0.91,1.29);
    
    // MaxLike
    bRatioTopCan->cd();
    bRatioTopFrame->GetYaxis()->SetRangeUser(1E-3,max);
    bRatioTopFrame->Draw();
    trueRes_->Draw("same");
    gMaxLikeStat->Draw("PE1same");
    if( par_->hasCorrPtGenAsym() ) {
      gMaxLikeStatCorr->Draw("PE1same");
      ptGenAsym_->Draw("same");
      if( par_->fitPtGenAsym() ) gPtGenAsym->Draw("PE1same");
    }
    if( par_->fitExtrapolatedSigma() ) fittedRes_->Draw("same");
    txt->Draw("same");
    legMaxLike->Draw("same");
    gPad->SetLogx();
    bRatioBottomPad->Draw();
    bRatioBottomPad->cd();
    bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
    bRatioBottomFrame->Draw();
    if( par_->hasCorrPtGenAsym() ) gMaxLikeRatioStatCorr->Draw("PE1same");
    else gMaxLikeRatioStat->Draw("PE1same");
    gPad->SetLogx();
    bRatioTopCan->SaveAs(par_->outNamePrefix()+"ExtraResoBottomRatio.eps","eps");

    // PtAsym
    TCanvas *bRatioTopCanAsym = util::HistOps::createRatioTopCanvas();
    TPad *bRatioBottomPadAsym = util::HistOps::createRatioBottomPad();
    bRatioTopCanAsym->cd();
    bRatioTopFrame->Draw();
    trueRes_->Draw("same");
    gPtAsymStat->Draw("PE1same");
    if( par_->hasCorrPtGenAsym() ) {
      gPtAsymStatCorr->Draw("PE1same");
      ptGenAsym_->Draw("same");
      if( par_->fitPtGenAsym() ) gPtGenAsym->Draw("PE1same");
    }
    if( par_->fitExtrapolatedSigma() ) fittedResAsym_->Draw("same");
    txt->Draw("same");
    legPtAsym->Draw("same");
    gPad->SetLogx();
    bRatioBottomPadAsym->Draw();
    bRatioBottomPadAsym->cd();
    bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
    bRatioBottomFrame->Draw();
    if( par_->hasCorrPtGenAsym() ) gPtAsymRatioStatCorr->Draw("PE1same");
    else gPtAsymRatioStat->Draw("PE1same");
    gPad->SetLogx();
    bRatioTopCanAsym->SaveAs(par_->outNamePrefix()+"ExtraResoPtAsymBottomRatio.eps","eps");

    // Adjust frame
    nLegEntries = 1;
    if( par_->fitRatio() ) nLegEntries++;
    if( par_->hasStartOffset() ) nLegEntries++;
    if( par_->fitExtrapolatedSigma() ) nLegEntries++;
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      h->SetBinContent(bin,1.);
    }
    h->SetLineStyle(trueRes_->GetLineStyle());
    h->SetLineColor(trueRes_->GetLineColor());
    h->GetYaxis()->SetTitle("#sigma_{Fit} / #sigma_{MC}");
    h->GetYaxis()->SetRangeUser(0.75,1.45+0.8*nLegEntries*lineHeight_);
    //h->GetYaxis()->SetRangeUser(0.85,1.85);


    // Recreate labels
    delete legMaxLike;
    delete legPtAsym;
    delete legComp;

    legPtAsym = util::LabelFactory::createLegendCol(nLegEntries,0.6);
    if( par_->hasCorrPtGenAsym() ) legPtAsym->AddEntry(gPtAsymRatioStatCorr,"p_{T} Asymmetry (Corrected)","P");
    else legPtAsym->AddEntry(gPtAsymRatioStat,"p_{T} Asymmetry","P");
    if( par_->fitExtrapolatedSigma() ) legPtAsym->AddEntry(hRatioFittedResAsym,"Interpolated Resolution","L");

    legMaxLike = util::LabelFactory::createLegendCol(nLegEntries,0.6);
    if( par_->hasCorrPtGenAsym() ) legMaxLike->AddEntry(gMaxLikeRatioStatCorr,"Likelihood (Corrected)","P");
    else legMaxLike->AddEntry(gMaxLikeRatioStat,"Likelihood Fit","P");
    if( par_->fitRatio() ) legMaxLike->AddEntry(lineFitRatio,"Mean fitted #bar{#sigma}","L");
    if( par_->hasStartOffset() ) legMaxLike->AddEntry(lineStartRes,"Resolution in spectrum","L");
    if( par_->fitExtrapolatedSigma() ) legMaxLike->AddEntry(hRatioFittedRes,"Interpolated Resolution","L");

    legComp = util::LabelFactory::createLegendCol(2,0.6);
    if( par_->hasCorrPtGenAsym() ) {
      legComp->AddEntry(gMaxLikeRatioStatCorr,"Likelihood (Corrected)","P");
      legComp->AddEntry(gPtAsymRatioStatCorr,"p_{T} Asymmetry (Corrected)","P");
    } else {
      legComp->AddEntry(gMaxLikeRatioStat,"Likelihood Fit","P");
      legComp->AddEntry(gPtAsymRatioStat,"p_{T} Asymmetry","P");
    }

    // Draw MaxLike results
    h->SetXTitle(par_->labelPtRef("MaxLike")+" (GeV)");
    h->SetYTitle("#sigma_{MaxLike} / #sigma_{MCTruth}");
    h->Draw();
    if( par_->fitRatio() ) lineFitRatio->Draw("same");
    if( par_->hasStartOffset() ) lineStartRes->Draw("same");
    if( par_->hasCorrPtGenAsym() ) gMaxLikeRatioStatCorr->Draw("PE1same");
    else gMaxLikeRatioStat->Draw("PE1same");
    if( par_->fitExtrapolatedSigma() ) hRatioFittedRes->Draw("Lsame");
    
    txt->Draw("same");
    legMaxLike->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoRatio.eps","eps");

    // Draw PtAsym results
    h->SetXTitle(par_->labelPtRef("PtAsym")+" (GeV)");
    h->SetYTitle("#sigma_{PtAsym} / #sigma_{MCTruth}");
    h->Draw();
    if( par_->hasCorrPtGenAsym() ) gPtAsymRatioStatCorr->Draw("PE1same");
    else gPtAsymRatioStat->Draw("PE1same");
    if( par_->fitExtrapolatedSigma() ) hRatioFittedResAsym->Draw("Lsame");
    txt->Draw("same");
    legPtAsym->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoPtAsymRatio.eps","eps");

    // Draw comparison
    h->SetXTitle(par_->labelPtRef()+" (GeV)");
    h->SetYTitle("#sigma_{Fit} / #sigma_{MCTruth}");
    h->Draw();
    if( par_->hasCorrPtGenAsym() ) {
      gMaxLikeRatioStatCorr->Draw("PE1same");
      gPtAsymRatioStatCorr->Draw("PE1same");
    } else {
      gMaxLikeRatioStat->Draw("PE1same");
      gPtAsymRatioStat->Draw("PE1same");
    }
    txt->Draw("same");
    legComp->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoCompRatio.eps","eps");

    if( par_->hasMCClosure() ) {
      h->GetYaxis()->SetTitle("#sigma_{fit} / #sigma_{MCClosure}");
      h->Draw();
      gMaxLikeRatioMCClosure->Draw("PE1same");
      txt->Draw("same");
      legComp->Draw("same");
      gPad->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"ExtraResoMCClosureRatio.eps","eps");
    }


    // Clean up
    delete h;
    delete gMaxLikeStat;
    delete gPtAsymStat;
    delete gMaxLikeStatCorr;
    delete gPtAsymStatCorr;
    delete gPtGenAsym;
    delete gMaxLikeRatioStat;
    delete gPtAsymRatioStat;
    delete gMaxLikeRatioMCClosure;
    delete lineFitRatio;
    delete lineStartRes;
    delete legMaxLike;
    delete legPtAsym;
    delete legComp;
    delete txt;
    delete can;
    delete bRatioTopFrame;
    delete bRatioBottomFrame;
    delete bRatioBottomPad;
    delete bRatioTopCan;
    delete bRatioBottomPadAsym;
    delete bRatioTopCanAsym;


    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotResolution(): Leaving method" << std::endl;
    }
  }
  
  

  
  // -------------------------------------------------------------------------------------
  void FittedResolution::plotSpectra() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());
      
      // PtGen spectrum
      TH1 *hPtGen = (*it)->getHistPtGen("PlotPtGen");
      util::HistOps::setAxisTitles(hPtGen,"p^{"+par_->labelPtGen()+"}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPtGen,5);
      hPtGen->SetMarkerStyle(20);

      TH1 *hPtGenJet1 = (*it)->getHistPtGenJet1("PlotPtGenJet1");
      util::HistOps::setAxisTitles(hPtGenJet1,par_->labelPtGen(),"GeV","jets",true);
      util::HistOps::setYRange(hPtGenJet1,5);
      hPtGenJet1->SetMarkerStyle(20);

      // PtAve spectrum
      TH1 *hPtAve = (*it)->getHistPtAve("PlotPtAve");
      util::HistOps::setAxisTitles(hPtAve,"p^{ave}_{T}","GeV","events");
      util::HistOps::setYRange(hPtAve,5);
      hPtAve->SetMarkerStyle(20);

      // Pdf
      TH1 *hPdfPtTrue = (*it)->getHistPdfPtTrue("PlotPdfPtTrue");
      hPdfPtTrue->SetLineColor(2);
      hPdfPtTrue->SetLineWidth(lineWidth_);
      util::HistOps::setAxisTitles(hPdfPtTrue,"p^{"+par_->labelTruth()+"}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPdfPtTrue,5);

      // Labels
      TPaveText *txt = 0;
      TLegend *legPtGen = 0;
      TLegend *legPtGenJet1 = 0;
      if( util::StyleSettings::style() == util::StyleSettings::Presentation ) {
	txt = util::LabelFactory::createPaveText(1);
	txt->AddText(par_->labelPtBin(ptBin)+",  "+par_->labelPt3Cut());
	legPtGen = util::LabelFactory::createLegendCol(2,0.45);
	legPtGenJet1 = util::LabelFactory::createLegendColWithOffset(2,-0.6,1);
      } else {
	txt = util::LabelFactory::createPaveText(4,-0.4);
	txt->AddText(par_->labelEtaBin());	
	txt->AddText(par_->labelPtSoftCut());
	txt->AddText(par_->labelPt3Cut());
	txt->AddText(par_->labelPtBin(ptBin));
	legPtGen = util::LabelFactory::createLegendCol(2,0.5);
	legPtGenJet1 = util::LabelFactory::createLegendCol(2,0.45);
      }
      legPtGen->AddEntry(hPtGen,"MC truth: p^{gen}_{T,1+2}","P");
      legPtGen->AddEntry(hPdfPtTrue,"Spectrum  f(p^{"+par_->labelTruth()+"}_{T})","L");
      legPtGenJet1->AddEntry(hPtGenJet1,"MC truth: p^{gen}_{T,1}","P");
      legPtGenJet1->AddEntry(hPdfPtTrue,"Spectrum  f(p^{"+par_->labelTruth()+"}_{T})","L");

      // Plot
      TCanvas *can = new TCanvas("PlotSpectrum","Spectrum ("+util::toTString(ptBin)+")",500,500);
      can->cd();

      if( !par_->isData() ) {
	hPtGen->Draw("PE1");
	if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
	txt->Draw("same");
	if( par_->fitMode() == FitModeMaxLikeFull ) legPtGen->Draw("same");
	can->SaveAs(par_->outNamePrefix()+"Spectrum_PtBin"+util::toTString(ptBin)+".eps","eps");

	hPtGenJet1->Draw("PE1");
	if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
	txt->Draw("same");
	if( par_->fitMode() == FitModeMaxLikeFull ) legPtGenJet1->Draw("same");
	can->SaveAs(par_->outNamePrefix()+"SpectrumJet1_PtBin"+util::toTString(ptBin)+".eps","eps");

	util::HistOps::setYRange(hPtGen,5,true);
	hPtGen->Draw("PE1");
	if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
	txt->Draw("same");
	if( par_->fitMode() == FitModeMaxLikeFull ) legPtGen->Draw("same");
	can->SetLogy();
	can->SaveAs(par_->outNamePrefix()+"SpectrumLog_PtBin"+util::toTString(ptBin)+".eps","eps");

	util::HistOps::setYRange(hPtGenJet1,5,true);
	hPtGenJet1->Draw("PE1");
	if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
	txt->Draw("same");
	if( par_->fitMode() == FitModeMaxLikeFull ) legPtGenJet1->Draw("same");
	can->SetLogy();
	can->SaveAs(par_->outNamePrefix()+"SpectrumJet1Log_PtBin"+util::toTString(ptBin)+".eps","eps");
      }

      can->SetLogy(0);
      hPtAve->Draw("PE1");
      txt->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"PtAveSpectrum_PtBin"+util::toTString(ptBin)+".eps","eps");

      util::HistOps::setYRange(hPtAve,5,true);
      hPtAve->Draw("PE1");
      txt->Draw("same");
      can->SetLogy();
      can->SaveAs(par_->outNamePrefix()+"PtAveSpectrumLog_PtBin"+util::toTString(ptBin)+".eps","eps");

      // Clean up
      delete hPtGen;
      delete hPtGenJet1;
      delete hPtAve;
      delete hPdfPtTrue;
      delete txt;
      delete legPtGen;
      delete legPtGenJet1;
      delete can;
    }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotResolutionDistributions() const {
    if( !par_->isData() ) {
      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  	  it != ptBins_.end(); it++) {
  	int ptBin = (it-ptBins_.begin());

	std::vector<TH1*> hPtGenRes;
	TCanvas *canGenRes = new TCanvas("canGenRes","PtGen Resolution (12)",500,500);
	TLegend *leg = util::LabelFactory::createLegendCol(3,par_->pt3Bins() ? 0.47 : 0.35);
	for(int i = 0; i < 3; ++i) {
	  if( par_->nPt3Cuts() < 3 ) continue;
	  int pt3Idx = 0;
	  if( i == 1 ) pt3Idx = par_->nPt3Cuts()/2;
	  else if( i == 2 ) pt3Idx = par_->nPt3Cuts()-1;

	  canGenRes->cd();
	  TH1 *h = (*it)->getHistMCRes(pt3Idx,"ResolutionDistributions:hResGen12"+util::toTString(pt3Idx));
	  util::HistOps::setAxisTitles(h,"p^{"+par_->labelMeas()+"}_{T,12} / p^{"+par_->labelTruth()+"}_{T,12}","GeV","jets",true);
	  h->SetLineColor(util::StyleSettings::color(ptBin));
	  h->GetXaxis()->SetRangeUser(-0.4,0.4);
	  if( i == 0 ) {
	    util::HistOps::setYRange(h,3);
	    h->Draw("HistE");
	  } else {
	    h->Draw("HistEsame");
	  }
	  leg->AddEntry(h,par_->labelPt3Cut(pt3Idx),"L");
	  hPtGenRes.push_back(h);
	}
	canGenRes->cd();
	leg->Draw("same");
	
	for(size_t i = 0; i < hPtGenRes.size(); ++i) {
	  delete hPtGenRes[i];
	}
	delete leg;
	delete canGenRes;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  void FittedResolution::plotControlDistributions() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());

      TCanvas *can = new TCanvas("canControlDistributions","Control Distributions",500,500);      
      can->cd();

      // Label
      TPaveText *label = util::LabelFactory::createPaveText(2);
      label->AddText(par_->labelLumi()+", "+par_->labelEtaBin()+", "+par_->labelDeltaPhiCut()+", "+par_->labelPt3Cut());
      label->AddText(par_->labelPtBin(ptBin));

      // Jet pt spectra
      TH1 *hPtJet1 = (*it)->getHistPtJet1("ControlDistributions:hPtJet1");
      util::HistOps::setYRange(hPtJet1,2,true);
      hPtJet1->SetMarkerStyle(20);
      hPtJet1->Draw("PE1");
      label->Draw("same");
      can->SetLogy();
      can->SaveAs(par_->outNamePrefix()+"PtJet1_PtBin"+util::toTString(ptBin)+".eps","eps");

      TH1 *hPtJet2 = (*it)->getHistPtJet2("ControlDistributions:hPtJet2");
      util::HistOps::setYRange(hPtJet2,2,true);
      hPtJet2->SetMarkerStyle(20);
      hPtJet2->Draw("PE1");
      label->Draw("same");
      can->SetLogy();
      can->SaveAs(par_->outNamePrefix()+"PtJet2_PtBin"+util::toTString(ptBin)+".eps","eps");

      TH1 *hPtJet3 = (*it)->getHistPtJet3("ControlDistributions:hPtJet3");
      util::HistOps::setYRange(hPtJet3,2,true);
      hPtJet3->SetMarkerStyle(20);
      hPtJet3->Draw("PE1");
      label->Draw("same");
      can->SetLogy();
      can->SaveAs(par_->outNamePrefix()+"PtJet3_PtBin"+util::toTString(ptBin)+".eps","eps");

      TH1 *hPtJet4 = (*it)->getHistPtJet4("ControlDistributions:hPtJet4");
      util::HistOps::setYRange(hPtJet4,2,true);
      hPtJet4->SetMarkerStyle(20);
      hPtJet4->Draw("PE1");
      label->Draw("same");
      can->SetLogy();
      can->SaveAs(par_->outNamePrefix()+"PtJet4_PtBin"+util::toTString(ptBin)+".eps","eps");

      TH1 *hEta = (*it)->getHistEta("ControlDistributions:hEta");
      hEta->SetMarkerStyle(20);
      util::HistOps::setYRange(hEta,2);
      hEta->Draw("PE1");
      label->Draw("same");
      can->SetLogy(0);
      can->SaveAs(par_->outNamePrefix()+"Eta_PtBin"+util::toTString(ptBin)+".eps","eps");

      TH1 *hDeltaPhi12 = (*it)->getHistDeltaPhi12("ControlDistributions:hDeltaPhi12");
      util::HistOps::setYRange(hDeltaPhi12,2);
      hDeltaPhi12->SetMarkerStyle(20);
      hDeltaPhi12->Draw("PE1");
      label->Draw("same");
      can->SetLogy(0);
      can->SaveAs(par_->outNamePrefix()+"DeltaPhi12_PtBin"+util::toTString(ptBin)+".eps","eps");

      // DeltaPt
      TH1 *hDeltaPtJet12 = (*it)->getHistDeltaPtJet12("ControlDistributions:hDeltaPtJet12");
      util::HistOps::setAxisTitles(hDeltaPtJet12,"#frac{1}{2} |p_{T,1} - p_{T,2}|","GeV","events",true);
      util::HistOps::setYRange(hDeltaPtJet12,3);
      //hDeltaPtJet12->GetYaxis()->SetRangeUser(3E-5,10.);
      hDeltaPtJet12->SetMarkerStyle(20);

      TF1* pdf = new TF1("pdf","gaus",0.,hDeltaPtJet12->GetXaxis()->GetBinUpEdge(hDeltaPtJet12->GetNbinsX()));
      pdf->SetLineColor(kRed);
      TH1 *htmp = (*it)->getHistPtGen("htmp");
      double sig = (*it)->fittedValue(0,0)*htmp->GetMean();
      delete htmp;
      //      double sig = (*it)->fittedValue(0,0)*(*it)->meanPt();
      double norm = hDeltaPtJet12->Integral("width");
      pdf->SetParameter(0,2.*norm/sqrt(2.*M_PI)/sig);
      pdf->SetParameter(1,0.);
      pdf->SetParameter(2,sig);
      pdf->SetLineWidth(1);
      std::cout << "NORM " << norm << ": " << pdf->Integral(0.,1000.) << std::endl;


      hDeltaPtJet12->Fit("gaus","I0QR","",0.,2.*hDeltaPtJet12->GetRMS());
      TF1 *fit = hDeltaPtJet12->GetFunction("gaus");
      fit->SetRange(0.,hDeltaPtJet12->GetXaxis()->GetBinUpEdge(hDeltaPtJet12->GetNbinsX()));
      fit->SetLineWidth(lineWidth_);
      fit->SetLineColor(4);
      double min = 3E-5;
      double max = 1E-1;
      std::vector<TLine*> lines;
      for(int k = 0; k < 4; ++k) {
	double x = (2.+k)*std::abs(fit->GetParameter(2));
	if( std::abs(x) > hDeltaPtJet12->GetXaxis()->GetBinUpEdge(hDeltaPtJet12->GetNbinsX()) ) continue;
	TLine *line = new TLine(x,min,x,max);
	line->SetLineStyle(2);
	line->SetLineWidth(lineWidth_);
	line->SetLineColor(4);
	lines.push_back(line);
      }
      hDeltaPtJet12->Draw("PE1");
      //fit->Draw("same");
      for(size_t k = 0; k < lines.size(); ++k) {
	//lines[k]->Draw("same");
      }
      pdf->Draw("same");
      label->Draw("same");
      //can->SetLogy(1);
      can->SaveAs(par_->outNamePrefix()+"DeltaPtJet12_PtBin"+util::toTString(ptBin)+".eps","eps");


      // Clean up
      delete hPtJet1;
      delete hPtJet2;
      delete hPtJet3;
      delete hPtJet4;      
      delete hEta;
      delete hDeltaPhi12;
      delete hDeltaPtJet12;
      for(size_t k = 0; k < lines.size(); ++k) {
	delete lines[k];
      }
      delete pdf;
      delete label;
      delete can;
    } // End of loop over ptbins
  }


  // -------------------------------------------------------------------------------------
  void FittedResolution::plotAdditionalJetActivity() const {
    
    // Components of asymmetry
    TH1 *hSumComps = new TH1D("AddJetAct:hSumComps",";"+par_->labelPtRef("PtAsym")+" (GeV);Standard Deviation",
			    nPtBins(),&(par_->ptBinEdges()->front()));
    hSumComps->GetXaxis()->SetMoreLogLabels();

    TH1 *hPtAsym = static_cast<TH1D*>(hSumComps->Clone("AddJetAct:hCompPtAsym"));
    hPtAsym->SetMarkerStyle(20);

    TH1 *hCompMCRes = static_cast<TH1D*>(hPtAsym->Clone("AddJetAct:hCompMCRes"));    
    hCompMCRes->SetMarkerStyle(21);
    util::HistOps::setColor(hCompMCRes,2);

    TH1 *hCompPJ3 = static_cast<TH1D*>(hPtAsym->Clone("AddJetAct:hCompPJ3"));
    hCompPJ3->SetMarkerStyle(22);
    util::HistOps::setColor(hCompPJ3,3);

    TH1 *hCompPSJ = static_cast<TH1D*>(hPtAsym->Clone("AddJetAct:hCompPSJ"));
    hCompPSJ->SetMarkerStyle(24);
    util::HistOps::setColor(hCompPSJ,4);

    TCanvas *can = new TCanvas("canAdditionalJetActivity","Additional Jet Activity",500,500);      
    can->cd();

    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());
      double sum2 = 0.;
      double sum2Err = 0.;

      // Label
      TPaveText *label = util::LabelFactory::createPaveText(2);
      label->AddText(par_->labelLumi()+", "+par_->labelEtaBin()+", "+par_->labelDeltaPhiCut()+", "+par_->labelPt3Cut());
      label->AddText(par_->labelPtBin(ptBin));

      TH1 *h = (*it)->getHistPJet3("AddJetAct:hPJ3");
      util::HistOps::setYRange(h,2);
      h->SetMarkerStyle(20);
      h->Draw("PE1");
      label->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"PJ3_PtBin"+util::toTString(ptBin)+".eps","eps");
      delete h;

      h = (*it)->getHistPJet3Rel("AddJetAct:hPJ3Rel");
      util::HistOps::setYRange(h,2);
      h->SetMarkerStyle(20);
      h->SetNdivisions(505);
      h->GetXaxis()->SetTitle("p_{||,3} / <p^{ave}_{T}>");
      h->Draw("PE1");
      label->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"PJ3Rel_PtBin"+util::toTString(ptBin)+".eps","eps");
      hCompPJ3->SetBinContent(ptBin+1,h->GetRMS());
      hCompPJ3->SetBinError(ptBin+1,h->GetRMSError());
      sum2 += h->GetRMS()*h->GetRMS();
      sum2Err += h->GetRMS()*h->GetRMS()*h->GetRMSError()*h->GetRMSError();
      delete h;

      h = (*it)->getHistPSJ("AddJetAct:hPSJ");
      util::HistOps::setYRange(h,2);
      h->SetMarkerStyle(20);
      h->GetXaxis()->SetTitle("p_{||,>3} (GeV)");
      h->Draw("PE1");
      label->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"PSJ_PtBin"+util::toTString(ptBin)+".eps","eps");
      delete h;

      h = (*it)->getHistPSJRel("AddJetAct:hPSJRel");
      util::HistOps::setYRange(h,2);
      h->SetMarkerStyle(20);
      h->GetXaxis()->SetTitle("p_{||,>3} / <p^{ave}_{T}>");
      h->SetNdivisions(505);
      h->Draw("PE1");
      label->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"PSJRel_PtBin"+util::toTString(ptBin)+".eps","eps");
      hCompPSJ->SetBinContent(ptBin+1,h->GetRMS());
      hCompPSJ->SetBinError(ptBin+1,h->GetRMSError());
      sum2 += h->GetRMS()*h->GetRMS();
      sum2Err += h->GetRMS()*h->GetRMS()*h->GetRMSError()*h->GetRMSError();
      delete h;

      h = (*it)->getHistPtAsym("AddJetAct:hPtAsym");
      double val = sqrt(2.)*h->GetRMS();
      val = sqrt( val*val - ptGenAsym_->Eval((*it)->meanPtAve())*ptGenAsym_->Eval((*it)->meanPtAve()) );
      hPtAsym->SetBinContent(ptBin+1,val);
      hPtAsym->SetBinError(ptBin+1,sqrt(2.)*h->GetRMSError());
      delete h;

      if( par_->hasMCClosure() ) {
  	h = (*it)->getHistMCRes("AddJetAct:MCRes");
	hCompMCRes->SetBinContent(ptBin+1,h->GetRMS());
	hCompMCRes->SetBinError(ptBin+1,h->GetRMSError());
	sum2 += h->GetRMS()*h->GetRMS();
	sum2Err += h->GetRMS()*h->GetRMS()*h->GetRMSError()*h->GetRMSError();
	delete h;

	hSumComps->SetBinContent(ptBin+1,sqrt(sum2));
	hSumComps->SetBinError(ptBin+1,sqrt(sum2Err/sum2));
      }

      // Clean up
      delete label;
    } // End of loop over pt bins

    TLegend *leg = util::LabelFactory::createLegendCol(5,0.5);
    leg->AddEntry(hPtAsym,"p_{T} asymmetry (#upoint#sqrt{2})","P");
    leg->AddEntry(hCompMCRes,"1: MC truth response","P");
    leg->AddEntry(hCompPJ3,"2: p_{||,3}","P");
    leg->AddEntry(hCompPSJ,"3: p_{||,>3}","P");
    leg->AddEntry(hSumComps,"1 #oplus 2 #oplus 3","L");

    if( par_->hasMCClosure() ) {
      can->cd();
      hSumComps->GetYaxis()->SetRangeUser(0.,0.2);
      hSumComps->Draw("HISTE");
      hPtAsym->Draw("PE1same");
      hCompMCRes->Draw("PE1same");
      hCompPJ3->Draw("PE1same");
      hCompPSJ->Draw("PE1same");
      leg->Draw("same");
      can->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"PtAsymmetryComponents.eps","eps");

      TH1 *hRatio = util::HistOps::createRatioPlot(hSumComps,hPtAsym,"StdDev(Comp/Asym))",0.7,1.3);
      TH1 *hRatioFrame = util::HistOps::createRatioFrame(hRatio,"StdDev(Comp/Asym)",0.7,1.3);
      hRatio->GetXaxis()->SetMoreLogLabels();
      hRatioFrame->Draw();
      hRatio->Draw("PE1same");
      can->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"PtAsymmetryComponentsRatio.eps","eps");
      delete hRatio;
      delete hRatioFrame;
    }

    delete leg;
    delete hPtAsym;
    delete hSumComps;
    delete hCompMCRes;
    delete hCompPJ3;
    delete hCompPSJ;
    delete can;
  }


  // -------------------------------------------------------------------------------------
  void FittedResolution::plotMCClosure() const {
    if( par_->hasMCClosure() ) {

      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  	  it != ptBins_.end(); it++) {
  	int ptBin = (it-ptBins_.begin());

  	double rMin = 0.4;
  	double rMax = 1.6;
  	TH1 *hMCRes = (*it)->getHistMCRes("PlotMCClosureMCRes");
  	util::HistOps::setAxisTitles(hMCRes,par_->xAxisTitleResponse(),"","jets",true);
  	hMCRes->SetMarkerStyle(20);
  	hMCRes->GetXaxis()->SetRangeUser(rMin,rMax);

  	// Fill histogram of extrapolated response
  	TH1 *hFitRes = new TH1D("PlotMCClosureFitRes","",5000,rMin,rMax);
  	hFitRes->SetLineWidth(lineWidth_);
  	hFitRes->SetLineColor(2);
  	hFitRes->SetXTitle(hMCRes->GetXaxis()->GetTitle());
  	hFitRes->SetYTitle(hMCRes->GetYaxis()->GetTitle());
  	hFitRes->GetXaxis()->SetRangeUser(rMin,rMax);
  	TH1 *hFitResCore = 0;
  	if( par_->respFuncType() == ResponseFunction::CrystalBall ) {
  	  hFitResCore = static_cast<TH1D*>(hFitRes->Clone("PlotMCClosureFitResCore"));
  	  hFitResCore->SetLineStyle(2);
  	}
  	TH1 *hFitResBins = static_cast<TH1D*>(hMCRes->Clone("PlotMCClosureFitResBins"));
  	hFitResBins->SetLineWidth(lineWidth_);
  	hFitResBins->SetLineColor(2);
  	hFitResBins->GetXaxis()->SetRangeUser(rMin,rMax);
  	std::vector<double> pars;
  	pars.push_back(1.);
  	for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
  	  pars.push_back((*it)->extrapolatedValue(parIdx));
  	}
  	for(int rBin = 1; rBin <= hFitRes->GetNbinsX(); ++rBin) {
  	  double res = hFitRes->GetBinCenter(rBin);
  	  hFitRes->SetBinContent(rBin,(*(par_->respFunc()))(res,pars));
  	  if( hFitResCore ) hFitResCore->SetBinContent(rBin,par_->respFunc()->pdfGauss(res,pars));
  	}
  	for(int rBin = 1; rBin <= hFitResBins->GetNbinsX(); ++rBin) {
  	  int minBin = hFitRes->FindBin(hFitResBins->GetXaxis()->GetBinLowEdge(rBin));
  	  int maxBin = hFitRes->FindBin(hFitResBins->GetXaxis()->GetBinUpEdge(rBin));
  	  hFitResBins->SetBinContent(rBin,hFitRes->Integral(minBin,maxBin,"width"));
  	  hFitResBins->SetBinError(rBin,0.);
  	}
  	hFitResBins->Scale(1./hFitResBins->GetBinWidth(1));

	
  	// Labels
  	TPaveText *txt = util::LabelFactory::createPaveText(2,-0.3);
	txt->AddText(par_->labelLumi());
	txt->AddText(par_->labelEtaBin());
  	TPaveText *txtLarge = util::LabelFactory::createPaveText(3);
	txtLarge->AddText("CMS QCD Simulation,  #sqrt{s} = 7 TeV");
	txtLarge->AddText("Anti-k_{T} 0.5 Calorimeter Dijets");
	txtLarge->AddText(par_->labelEtaBin()+", "+par_->labelPtBin(ptBin));
  	TLegend *leg = util::LabelFactory::createLegendCol(2,0.7);
  	leg->AddEntry(hMCRes,"MC truth, "+util::toTString((*it)->ptMin())+" < p^{"+par_->labelTruth()+"}_{T} < "+util::toTString((*it)->ptMax())+" GeV","P");
  	leg->AddEntry(hFitRes,"Fit result, "+par_->labelPtBin(ptBin),"L");
  	TLegend *legFits = util::LabelFactory::createLegendColWithOffset(nMCClosureResFits(),0.7,txt->GetSize());
	for(int i = 0; i < nMCClosureResFits(); ++i) {
	  legFits->AddEntry(mcClosureFits_[i][ptBin],mcClosureResoLabels_[i],"L");
	}

	hMCRes->GetYaxis()->SetRangeUser(3E-3,95);
  	util::HistOps::setYRange(hFitRes,txt->GetSize()+3);
  	util::HistOps::setYRange(hFitResBins,txt->GetSize()+3);

  	// Create Canvas
  	TCanvas *can = new TCanvas("PlotMCClosure","MCClosure ("+util::toTString(ptBin)+")",500,500);
  	can->cd();

  	hMCRes->Draw("PE1");
  	txtLarge->Draw("same");
	can->SetLogy();
  	can->SaveAs(par_->outNamePrefix()+"MCClosureResponse_PtBin"+util::toTString(ptBin)+".eps","eps");

  	hFitRes->Draw("L");	
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
	can->SetLogy(0);
  	can->SaveAs(par_->outNamePrefix()+"MCClosure_PtBin"+util::toTString(ptBin)+".eps","eps");

  	hFitResBins->Draw("HIST");	
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
  	can->SaveAs(par_->outNamePrefix()+"MCClosureHist_PtBin"+util::toTString(ptBin)+".eps","eps");

  	hFitRes->Draw("L");	
  	for(int i = 0; i < nMCClosureResFits(); ++i) {
  	  mcClosureFits_[i][ptBin]->Draw("same");
  	}
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
	legFits->Draw("same");
  	can->SaveAs(par_->outNamePrefix()+"MCClosureFits_PtBin"+util::toTString(ptBin)+".eps","eps");

	util::HistOps::setYRange(hFitRes,2,3E-4);
	hFitRes->GetXaxis()->SetRangeUser(0.,2.);
  	hFitRes->Draw("L");	
  	if( hFitResCore ) hFitResCore->Draw("Lsame");
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
  	gPad->SetLogy();
  	can->SaveAs(par_->outNamePrefix()+"MCClosureLog_PtBin"+util::toTString(ptBin)+".eps","eps");

	// Scan for startpoint of tails
	hMCRes->Fit("gaus","I0QR","",hMCRes->GetMean()-2.*hMCRes->GetRMS(),hMCRes->GetMean()+2.*hMCRes->GetRMS());
	TF1 *fit = hMCRes->GetFunction("gaus");
	fit->SetRange(0.,2.);
	fit->SetLineWidth(2);
	fit->SetLineColor(2);
	double mean = fit->GetParameter(1);
	double sigma = fit->GetParameter(2);
	double min = 0.;
	double max = 0.;
	util::HistOps::findYRange(hMCRes,min,max);
	min = 3E-4;
	TLine *lFitRangeLeft = new TLine(hMCRes->GetMean()-2.*hMCRes->GetRMS(),min,
					 hMCRes->GetMean()-2.*hMCRes->GetRMS(),max);
	lFitRangeLeft->SetLineWidth(1);
	lFitRangeLeft->SetLineColor(4);
	TLine *lFitRangeRight = new TLine(hMCRes->GetMean()+2.*hMCRes->GetRMS(),min,
					  hMCRes->GetMean()+2.*hMCRes->GetRMS(),max);
	lFitRangeRight->SetLineWidth(1);
	lFitRangeRight->SetLineColor(4);
	std::vector<TLine*> lines;
	for(int k = -6; k < 6; ++k) {
	  if( k == 0 ) continue;
	  double x = mean - (2.-k)*sigma;
	  if( k > 0 ) x = mean + (2.+k)*sigma;
	  if( x < 0. || x > 2. ) continue;
	  TLine *line = new TLine(x,min,x,max);
	  line->SetLineStyle(2);
	  line->SetLineWidth(1);
	  line->SetLineColor(4);
	  lines.push_back(line);
	}
	util::HistOps::setYRange(hMCRes,2,min);
	hMCRes->GetXaxis()->SetRangeUser(0.,2.);
	hMCRes->Draw("PE1");
	fit->Draw("same");
	lFitRangeLeft->Draw("same");
	lFitRangeRight->Draw("same");
	for(size_t k = 0; k < lines.size(); ++k) {
	  lines[k]->Draw("same");
	}
  	txt->Draw("same");
  	leg->Draw("same");
  	gPad->SetLogy();
  	can->SaveAs(par_->outNamePrefix()+"MCClosureTails_PtBin"+util::toTString(ptBin)+".eps","eps");


  	gPad->SetLogy(0);
  	TH1 *hFitResBinsRatioFrame = util::HistOps::createRatioFrame(hMCRes,"Fit / MC truth",0.5,2.0);
  	hFitResBinsRatioFrame->GetXaxis()->SetRangeUser(rMin,rMax);
  	TH1 *hFitResBinsRatio = util::HistOps::createRatioPlot(hFitResBins,hMCRes);
  	util::HistOps::setColor(hFitResBinsRatio,1);
  	hFitResBinsRatioFrame->Draw();	
  	hFitResBinsRatio->Draw("PE1same");
  	leg->Draw("same");
  	can->SaveAs(par_->outNamePrefix()+"MCClosureRatio_PtBin"+util::toTString(ptBin)+".eps","eps");
	
	
  	// Clean up
  	delete hMCRes;
  	delete hFitRes;
  	delete hFitResBins;
  	delete hFitResBinsRatio;
  	delete hFitResBinsRatioFrame;
  	if( hFitResCore ) delete hFitResCore;
	delete lFitRangeLeft;
	delete lFitRangeRight;
	for(size_t k = 0; k < lines.size(); ++k) {
	  delete lines[k];
	}
  	delete txt;
	delete txtLarge;
  	delete leg;
	delete legFits;
  	delete can;
      } // End of loop over pt bins


      // Compare width to MC truth resolution and plot scale
      std::vector<TF1*> fitsMCTruthReso;
      //       fitsMCTruthReso.push_back(new TF1("fitsMCTruthReso:1","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",ptMin(),ptMax()));
      //       fitsMCTruthReso.back()->SetParameter(0,3.60458);
      //       fitsMCTruthReso.back()->SetParameter(1,1.12229);
      //       fitsMCTruthReso.back()->SetParameter(2,0.0363994);
      for(size_t i = 0; i < fitsMCTruthReso.size(); ++i) {
	fitsMCTruthReso[i]->SetLineColor(util::StyleSettings::color(3+i));
	fitsMCTruthReso[i]->SetLineWidth(lineWidth_);
	fitsMCTruthReso[i]->SetLineStyle(2+i);
      }      
      TH1 *hFrameReso = new TH1D("hFrameReso",";"+par_->labelPtGen()+" (GeV);#sigma / "+par_->labelPtGen(),1000,ptMin(),ptMax());
      hFrameReso->GetYaxis()->SetRangeUser(0.,0.3);
      TLegend *leg = util::LabelFactory::createLegendCol(nMCClosureResFits(),0.6);
      for(int i = 0; i < nMCClosureResFits(); ++i) {
	leg->AddEntry(mcClosureGReso_[i],mcClosureResoLabels_[i],"P");
      }
      TCanvas *can1 = new TCanvas("PlotMCClosureResolution","MCClosure Resolution",500,500);
      can1->cd();
      hFrameReso->Draw();
      trueRes_->Draw("same");
      for(size_t i = 0; i < fitsMCTruthReso.size(); ++i) {
	fitsMCTruthReso[i]->Draw("same");
      }
      for(int i = 0; i < nMCClosureResFits(); ++i) {
	mcClosureGReso_[i]->Draw("PE1same");
      }
      leg->Draw("same");
      can1->SetLogx();
      can1->SaveAs(par_->outNamePrefix()+"MCClosureReso.eps","eps");


      TH1 *hFrameScale = util::HistOps::createRatioFrame(ptMin(),ptMax(),0.9,1.2,par_->labelPtGen()+" (GeV)","< R >");
      TCanvas *can2 = new TCanvas("PlotMCClosureScale","MCClosure Scale",500,500);
      can2->cd();
      hFrameScale->Draw();
      for(int i = 0; i < nMCClosureResFits(); ++i) {
	mcClosureGScale_[i]->Draw("PE1same");
      }
      leg->Draw("same");
      can2->SetLogx();
      can2->SaveAs(par_->outNamePrefix()+"MCClosureScale.eps","eps");
      

      for(size_t i = 0; i < fitsMCTruthReso.size(); ++i) {
	delete fitsMCTruthReso[i];
      }
      delete hFrameReso;
      delete hFrameScale;
      delete leg;
      delete can1;
      delete can2;

    } else {
      std::cerr << "No MCClosure response distribution available." << std::endl;
    }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::print() const {
    std::cout << std::endl;
    std::cout << "Fitted Resolution (MaxLike)" << std::endl;
    for(int k = 0; k < fittedRes_->GetNpar(); ++k) {
      std::cout << " & $" << fittedRes_->GetParameter(k) << " \\pm " << fittedRes_->GetParError(k) << "$" << std::flush;
    }
    std::cout << " \\\\" << std::endl;
    std::cout << "Fitted Resolution (PtAsym)" << std::endl;
    for(int k = 0; k < fittedResAsym_->GetNpar(); ++k) {
      std::cout << " & $" << fittedResAsym_->GetParameter(k) << " \\pm " << fittedResAsym_->GetParError(k) << "$" << std::flush;
    }
    std::cout << " \\\\" << std::endl;
    
    if( par_->hasCorrPtGenAsym() ) {
      std::cout << std::endl;
      std::cout << "PtGen Asymmetry" << std::endl;
      for(int k = 0; k < ptGenAsym_->GetNpar(); ++k) {
	std::cout << " & $" << ptGenAsym_->GetParameter(k) << " \\pm " << ptGenAsym_->GetParError(k) << "$" << std::flush;
      }
      std::cout << " \\\\" << std::endl;
    }
      

    std::cout << "\n\nFITTED RESOLUTION\n";
    // Loop over ptbins
    for(size_t bin = 0; bin < nPtBins(); bin++) {
      std::cout << bin << ": " << meanPt(bin) << std::flush;
      std::cout << " (" << ptMin(bin) << " - " << ptMax(bin) << "): " << std::flush;
      for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
	if( parIdx == 0 ) {
	  double val = extrapolatedValue(bin,0);
	  if( par_->hasCorrPtGenAsym() ) val = sqrt( val*val - ptGenAsym_->Eval(meanPt(bin))*ptGenAsym_->Eval(meanPt(bin)) );
	  std::cout << val << " (" << uncertStat(bin,parIdx) << ")" << std::flush;
	  std::cout << ", " << trueRes_->Eval(meanPt(bin)) << std::flush;
	}
	if( parIdx < par_->nFittedPars()-1 ) std::cout << " |" << std::flush;
	else std::cout << std::endl;
      }
    } 


//     std::cout << "\n\nFITTED RESOLUTION\n";
//     // Loop over ptbins
//     for(size_t bin = 0; bin < nPtBins(); bin++) {
//       std::cout << bin << ": " << meanPt(bin) << std::flush;
//       std::cout << " (" << ptMin(bin) << " - " << ptMax(bin) << "): " << std::flush;
//       std::cout << ptBins_[bin]->fittedValue(0,par_->stdSelIdx()) << ", " << std::flush;
//       std::cout << ptBins_[bin]->fittedValue(0,par_->stdSelIdx())*meanPt(bin) << std::endl;
//     }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::printLaTeX() const {
    for(size_t bin = 0; bin < nPtBins(); bin++) {
      std::cout << "$ " << bin << " $ & $ ";
      std::cout << ptMin(bin) << " - " << ptMax(bin) << " $ & $ ";
      std::cout << meanPt(bin) << " $ & $";
      for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
	std::cout << extrapolatedValue(bin,parIdx) << " \\pm " << uncertStat(bin,parIdx) << " $ ";
	if( parIdx == 0 ) std::cout << " & $ " << trueRes_->Eval(meanPt(bin)) << " $ ";
      }
      std::cout << " \\\\ \n";
    } 
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::printPoints() const {
    TString x = "ptMC";
    TString xe = "ptMCErr";
    TString y = "resMC";
    TString ye = "resMCErr";
    if( par_->isData() ) {
      x = "ptData";
      xe = "ptDataErr";
      y = "resData";
      ye = "resDataErr";
    }
    std::cout << std::endl;
    for(size_t bin = 0; bin < nPtBins(); bin++) {
      std::cout << x << ".push_back(" << meanPt(bin) << ");" << std::endl;
      std::cout << xe << ".push_back(" << meanPtUncert(bin) << ");" << std::endl;
      double val = extrapolatedValue(bin,0);


// //       // +/- variation on linearity
//       double diff = std::abs(fittedValue(bin,0,0)-val);
//       //val = val - 0.5*diff;
//       //val = val + 0.5*diff;
				 

      double pli = 0.;
      if( par_->hasCorrPtGenAsym() ) pli = ptGenAsym_->Eval(meanPt(bin));
//       // +/- 0.25 variation on PLI
//       //pli = 0.75*pli;
//       //pli = 1.25*pli;



      val = sqrt( val*val - pli*pli );
      if( val == val ) {
	std::cout << y << ".push_back(" << val << ");" << std::endl;
	std::cout << ye << ".push_back(" << uncertStat(bin,0) << ");" << std::endl;
      } else {
	std::cout << y << ".push_back(" << 0 << ");" << std::endl;
	std::cout << ye << ".push_back(" << 100000. << ");" << std::endl;
      }
    } 

//     for(size_t bin = 0; bin < nPtBins(); bin++) {
//       std::cout << meanPt(bin) << "\t   " << meanPtUncert(bin) << "\t   " << std::flush;
//       //std::cout << extrapolatedValue(bin,0) << "\t   " << std::flush;
//       double val = extrapolatedValue(bin,0);
//       if( par_->hasCorrPtGenAsym() ) val = sqrt( val*val - ptGenAsym_->Eval(meanPt(bin))*ptGenAsym_->Eval(meanPt(bin)) );
//       std::cout << val << "\t   " << uncertStat(bin,0) << std::endl;
//     } 
  }

    

  //! Write resolution to root file (via histograms; how to 
  //! write a TGraph to file??)
  // -------------------------------------------------------------------------------------
  void FittedResolution::writeRootOutput() const {
    TGraphAsymmErrors *g = getTGraphOfResolution("MaxLike","Statistic",true);
    TGraphAsymmErrors *gUncorr = getTGraphOfResolution("MaxLike","Statistic");
    
    TH1* hPt = new TH1D("hPtMean","",g->GetN(),&(par_->ptBinEdges()->front()));
    TH1* hRes = new TH1D("hResolution","",g->GetN(),&(par_->ptBinEdges()->front()));
    TH1* hExt = new TH1D("hExtrapolationResult","",g->GetN(),&(par_->ptBinEdges()->front()));
    for(int i = 0; i < g->GetN(); ++i) {
      int bin = 1+i;
      hPt->SetBinContent(bin,g->GetX()[i]);
      hPt->SetBinError(bin,g->GetEXhigh()[i]);
      hRes->SetBinContent(bin,g->GetY()[i]);
      hRes->SetBinError(bin,g->GetEYhigh()[i]);
      hExt->SetBinContent(bin,gUncorr->GetY()[i]);
      hExt->SetBinError(bin,gUncorr->GetEYhigh()[i]);
    }


    TFile file(par_->outNamePrefix()+".root","RECREATE");
    file.WriteTObject(hPt);
    file.WriteTObject(hRes);
    file.WriteTObject(hExt);
    file.WriteTObject(trueRes_);
    file.WriteTObject(ptGenAsym_);
    file.Close();

    delete hPt;
    delete hRes;
    delete g;
  }


    
  // -------------------------------------------------------------------------------------
  void FittedResolution::createSlides() const {
      
    // ----- Create slides with MC closure plots ------------------------------------
      
    // Open file
    TString name = par_->outNamePrefix();
    name += "SlidesMCClosure.tex";
    std::ofstream oFile(name);
    // Parse tex-code
    oFile << "% ----- MC closure plots ---------------------------" << std::endl;
    int nSlides = nPtBins()/6;
    if( nPtBins()%6 > 0 ) nSlides++;
    for(int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{MC closure in various \\pt bins (" << slide+1 << "/" << nSlides << ")}\n";
      oFile << "  \\vskip-0.7cm\n";
      oFile << "  \\begin{columns}[t] \n";
      for(int col = 0; col < 3; ++col) {
	oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	oFile << "      \\begin{center} \n";
	for(int row = 0; row < 2; ++row) {
	  int ptBin = 6*slide + 3*row + col;
	  if( ptBin < nPtBins() ) {
	    oFile << "        \\includegraphics[width=\\textwidth]{figures/ResFit_Spring10QCDFlat_CB_Eta0_MCClosure_PtBin" << ptBin << "}\\\\ \n";
	  }
	}
	oFile << "      \\end{center} \n";
	oFile << "    \\end{column} \n";
      }
      oFile << "  \\end{columns} \n";
      oFile << "\\end{frame} \n";
    }
    oFile.close();



    // ----- Create slides of all fit results ----------------------------------------

    if( par_->respFuncType() == ResponseFunction::CrystalBall ) {
      name = par_->outNamePrefix();
      name += "SlidesAllResults.tex";
      std::ofstream oFile(name);

      oFile << "\n\n\n% ----- Fit result plots ---------------------------" << std::endl;
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	  it != ptBins_.end(); it++) {
	int ptBin = (it-ptBins_.begin());
	oFile << "\n% --------------------------------------------------\n";
	oFile << "\\begin{frame}\n";
	oFile << "  \\frametitle{Fit results and MC closure $" << (*it)->ptMinStr() << " < \\pt < " << (*it)->ptMaxStr() << "\\gev$}\n";
	oFile << "  \\vskip-0.5cm\n";
	oFile << "  \\begin{columns}[t]\n";
	for(int col = 0; col < 3; ++col) {
	  oFile << "    \\begin{column}{0.3333\\textwidth}\n";
	  oFile << "      \\begin{center}\n";
	  oFile << "        \\includegraphics[width=\\textwidth]{figures/ResFit_Spring10QCDFlat_CB_Eta0_ExtrapolatedPar" << col << "_PtBin" << ptBin << "}\\\\ \n";
	  if( col == 0 ) {
	    oFile << "        \\includegraphics[width=\\textwidth]{figures/ResFit_Spring10QCDFlat_CB_Eta0_Spectrum_PtBin" << ptBin << "} \n";
	  } else if( col == 1 ) {
	    oFile << "        \\includegraphics[width=\\textwidth]{figures/ResFit_Spring10QCDFlat_CB_Eta0_MCClosure_PtBin" << ptBin << "} \n";
	  }
	  oFile << "      \\end{center} \n";
	  oFile << "    \\end{column} \n";
	}
	oFile << "  \\end{columns} \n";
	oFile << "\\end{frame} \n";
      }
      oFile.close();
    }
    else if( par_->respFuncType() == ResponseFunction::Gauss ) {
      // Create slides of all fit results
      name = par_->outNamePrefix();
      name += "SlidesAllResults.tex";
      std::ofstream oFile(name);

      // Create slides with result plots
      oFile << "% ----- All fit results ---------------------------" << std::endl;
      int nSlides = nPtBins()/3;
      if( nPtBins()%3 > 0 ) nSlides++;
      for(int slide = 0; slide < nSlides; ++slide) {
	oFile << "\n% --------------------------------------------------\n";
	oFile << "\\begin{frame}\n";
	oFile << "  \\frametitle{Fit results in various \\pt bins (" << slide+1 << "/" << nSlides << ")}\n";
	oFile << "  \\vskip-1cm\n";
	oFile << "  \\begin{columns}[T] \n";
	for(int col = 0; col < 3; ++col) {
	  oFile << "    \\begin{column}{0.3333\\textwidth} \n";
	  oFile << "      \\begin{center} \n";
	  int ptBin = 3*slide + col;
	  if( ptBin < nPtBins() ) {
	    oFile << "        \\includegraphics[width=\\textwidth]{figures/" << par_->outNamePrefix() << "ExtrapolatedPar0_PtBin" << ptBin << "}\\\\ \n";
	    oFile << "        \\includegraphics[width=\\textwidth]{figures/" << par_->outNamePrefix() << "Spectrum_PtBin" << ptBin << "} \n";
	  }
	  oFile << "      \\end{center} \n";
	  oFile << "    \\end{column} \n";
	}
	oFile << "  \\end{columns} \n";
	oFile << "\\end{frame} \n";
      }
      oFile.close();
    }




    // ----- Create slides with asymmetry distributions ---------------------------------

    // Create slides of asymmetry distributions
    name = par_->outNamePrefix();
    name += "SlidesAsymDistributions.tex";
    oFile.open(name);

    oFile << "\n\n\n% ----- Asymmetry distributions --------------------" << std::endl;
    for(int d = 0; d < par_->nFittedPars()+2; d++) {
      if( d == 0 ) oFile << "\n% ------- Linear scale -----------------------------" << std::endl;
      else if( d == 1 ) oFile << "\n% ------- Log scale --------------------------------" << std::endl;
      else oFile << "\n% ------- Variation par " << d-2 << "--------------------------" << std::endl;
      nSlides = nPtBins()/3;
      if( nPtBins()%3 > 0 ) nSlides++;
      for(int slide = 0; slide < nSlides; ++slide) {
	oFile << "\n% --------------------------------------------------\n";
	oFile << "\\begin{frame}\n";
	oFile << "  \\frametitle{Measured and predicted asymmetry ";
	if( d == 1 ) oFile << "-- log ";
	else if( d > 1 ) oFile << " -- varied $" << par_->parLabelTex(d-2) << "$ ";
	oFile << "(" << slide+1 << "/" << nSlides << ")}\n";
	oFile << "  \\begin{columns}[T]\n";
	for(int col = 0; col < 3; ++col) {
	  oFile << "    \\begin{column}{0.3333\\textwidth}\n";
	  oFile << "    \\centering\n";
	  int ptBin = 3*slide + col;
	  if( ptBin < nPtBins() ) {
	    // Asymmetry distribution
	    oFile << "      \\includegraphics[width=\\textwidth]{figures/" << par_->outNamePrefix() << "PtAsymmetry";
	    if( d == 1 || (d > 1 && par_->respFuncType() == ResponseFunction::CrystalBall) ) oFile << "Log";
	    if( d > 1 ) oFile << "_Var" << d-2;
	    oFile << "_PtBin" << ptBin << "}\\\\ \n";
	    // Ratio
	    oFile << "      \\includegraphics[width=\\textwidth]{figures/" << par_->outNamePrefix() << "PtAsymmetryRatio";
	    if( d > 1 ) oFile << "_Var" << d-2;
	    oFile << "_PtBin" << ptBin << "}\\\\ \n";
	  }
	  oFile << "    \\end{column} \n";
	}
	oFile << "  \\end{columns} \n";
	oFile << "\\end{frame} \n";
      }
    }
    oFile.close();
  }

    

  //! Fit MC closure response distributions with different
  //! constraints and create graphs of scale and resolution
  // -------------------------------------------------------------------------------------
  void FittedResolution::fitMCClosure() {
    if( par_->verbosity() > 1 ) std::cout << "FittedResolution::fitMCClosure() Entering\n";
    // Different fits of MC truth resolution
    mcClosureResoLabels_.push_back("<R> = free parameter");
    mcClosureResoLabels_.push_back("<R> = 1");

    mcClosureFits_ = std::vector< std::vector<TF1*> >(nMCClosureResFits());
    
    std::vector<double> ptMean;
    std::vector<double> ptMeanErr;
    std::vector< std::vector<double> > mcClosureReso(nMCClosureResFits());
    std::vector< std::vector<double> > mcClosureResoErr(nMCClosureResFits());     
    std::vector< std::vector<double> > mcClosureScale(nMCClosureResFits());
    std::vector< std::vector<double> > mcClosureScaleErr(nMCClosureResFits());     
    
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());
      ptMean.push_back((*it)->meanPt());
      ptMeanErr.push_back((*it)->meanPtUncert());

      // Fit MC truth resolution
      TH1 *hMCRes = (*it)->getHistMCRes("PlotMCClosureMCRes");
      double min = hMCRes->GetMean() - 3.*hMCRes->GetRMS();
      double max = hMCRes->GetMean() + 3.*hMCRes->GetRMS();
      hMCRes->Fit("gaus","I0Q","",min,max);
      TF1 *fit = hMCRes->GetFunction("gaus");
      min = fit->GetParameter(1) - 1.5*fit->GetParameter(2);
      max = fit->GetParameter(1) + 1.5*fit->GetParameter(2);
      for(int i = 0; i < nMCClosureResFits(); ++i) {
	TString name = "mcClosureFits"+util::toTString(i)+"_PtBin"+util::toTString(ptBin);
	mcClosureFits_[i].push_back(new TF1(name,"gaus",min,max));
	mcClosureFits_[i].back()->SetParameter(0,5.);
	if( i == 0 ) mcClosureFits_[i].back()->SetParameter(1,1.);
	else if( i == 1 ) mcClosureFits_[i].back()->FixParameter(1,1.);
	mcClosureFits_[i].back()->SetParameter(2,0.1);
	mcClosureFits_[i].back()->SetLineStyle(2);	  
	mcClosureFits_[i].back()->SetLineColor(util::StyleSettings::color(2+i));
	mcClosureFits_[i].back()->SetLineWidth(1);

	if( hMCRes->Fit(mcClosureFits_[i].back(),"I0QBR") == 0 ) {
	  mcClosureReso[i].push_back(std::abs(mcClosureFits_[i].back()->GetParameter(2))); 
	  mcClosureResoErr[i].push_back(mcClosureFits_[i].back()->GetParError(2));
	  mcClosureScale[i].push_back(mcClosureFits_[i].back()->GetParameter(1)); 
	  mcClosureScaleErr[i].push_back(mcClosureFits_[i].back()->GetParError(1));
	} else {
	  mcClosureReso[i].push_back(0.); 
	  mcClosureResoErr[i].push_back(0.);
	  mcClosureScale[i].push_back(0.); 
	  mcClosureScaleErr[i].push_back(0.);
	}
      }
      delete hMCRes;
    } // End of loop over ptbins

    // Create graphs
    mcClosureGReso_ = std::vector<TGraphAsymmErrors*>(nMCClosureResFits());
    mcClosureGScale_ = std::vector<TGraphAsymmErrors*>(nMCClosureResFits());
    for(int i = 0; i < nMCClosureResFits(); ++i) {
      mcClosureGReso_[i] = new TGraphAsymmErrors(ptMean.size(),&(ptMean.front()),&(mcClosureReso[i].front()),
						 &(ptMeanErr.front()),&(ptMeanErr.front()),
						 &(mcClosureResoErr[i].front()),&(mcClosureResoErr[i].front()));
      mcClosureGReso_[i]->SetMarkerStyle(20+i);
      mcClosureGReso_[i]->SetMarkerColor(util::StyleSettings::color(i));
      mcClosureGReso_[i]->SetLineColor(mcClosureGReso_[i]->GetMarkerColor());
      
      mcClosureGScale_[i] = new TGraphAsymmErrors(ptMean.size(),&(ptMean.front()),&(mcClosureScale[i].front()),
						  &(ptMeanErr.front()),&(ptMeanErr.front()),
						  &(mcClosureScaleErr[i].front()),&(mcClosureScaleErr[i].front()));
      mcClosureGScale_[i]->SetMarkerStyle(20+i);
      mcClosureGScale_[i]->SetMarkerColor(util::StyleSettings::color(i));
      mcClosureGScale_[i]->SetLineColor(mcClosureGScale_[i]->GetMarkerColor());
    }

    if( par_->verbosity() > 1 ) std::cout << "FittedResolution::fitMCClosure() Leaving\n";
  }




  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors *FittedResolution::getTGraphOfResolution(const TString &method, const TString &uncertainty, bool corrected) const {
    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution: Creating graph of resolution" << std::endl;
    }
    std::vector<double> x;
    std::vector<double> ex;
    std::vector<double> y;
    std::vector<double> eyDown;
    std::vector<double> eyUp;

    for(int i = 0; i < nPtBins(); ++i) {
      if( method == "MaxLike" ) {
	x.push_back(ptBins_[i]->meanPt());
	ex.push_back(ptBins_[i]->meanPtUncert());
	double val = ptBins_[i]->extrapolatedValue(0);
	if( corrected ) {
	  double pli = par_->scalePli()*ptGenAsym_->Eval(x[i]);
	  val = sqrt( val*val - pli*pli );
	}
	y.push_back(val);
	eyDown.push_back(ptBins_[i]->uncertStatDown(0));
	eyUp.push_back(ptBins_[i]->uncertStatUp(0));
      } else if( method == "PtAsym" ) {
	x.push_back(ptBins_[i]->meanPtAsym());
	ex.push_back(ptBins_[i]->meanPtAsymUncert());
	double val = ptBins_[i]->extrapolatedAsym();
	if( corrected ) val = sqrt( val*val - ptGenAsym_->Eval(x[i])*ptGenAsym_->Eval(x[i]) );
	y.push_back(val);
	eyDown.push_back(ptBins_[i]->uncertDownAsym());
	eyUp.push_back(ptBins_[i]->uncertUpAsym());
      } else if( method == "PtGenAsym" ) {
	x.push_back(ptBins_[i]->meanPtAsym());
	ex.push_back(ptBins_[i]->meanPtAsymUncert());
	y.push_back(ptBins_[i]->extrapolatedGenAsym());
	eyDown.push_back(ptBins_[i]->uncertDownGenAsym());
	eyUp.push_back(ptBins_[i]->uncertUpGenAsym());
      } else {
	std::cerr << "ERROR: FittedResolution::getTGraphOfResolution(): Unknown resolution method '" << method << "'" << std::endl;
	exit(-1);
      }
    }
    TGraphAsymmErrors *graph = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						     &(ex.front()),&(ex.front()),
						     &(eyDown.front()),&(eyUp.front()));

    if( method == "PtGenAsym" ) {
      graph->SetMarkerStyle(26);
    } else {
      if( corrected ) {
	if( util::StyleSettings::style() == util::StyleSettings::Presentation ) graph->SetMarkerStyle(21);
	else graph->SetMarkerStyle(25);
	if( par_->isData() ) {
	  graph->SetMarkerColor(2);
	  graph->SetLineColor(graph->GetMarkerColor());
	}
      } else {
	if( util::StyleSettings::style() == util::StyleSettings::Presentation ) graph->SetMarkerStyle(20);
	else graph->SetMarkerStyle(24);
	if( par_->isData() ) {
	  graph->SetMarkerColor(2);
	  graph->SetLineColor(graph->GetMarkerColor());
	}
      }
    }    
    return graph;
  }
  
  
  
  // -------------------------------------------------------------------------------------
  double FittedResolution::gaussian(double *x, double *par) {
    double u = (x[0]-1.)/par[0];
    return par[0] > 0. ? exp(-0.5*u*u)/sqrt(2.*M_PI)/par[0] : 0.;
  }
}
  
