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

    trueRes_ = new TF1("FittedResolution::trueRes",
		       "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		       ptMin_,ptMax_);
    for(int i = 0; i < 3; i++) {
      trueRes_->SetParameter(i,par_->trueGaussResPar(i));
    }
    trueRes_->SetLineWidth(lineWidth_);
    trueRes_->SetLineColor(4);
    trueRes_->SetLineStyle(1);

    // Fit extrapolated resolution using
    // statistic uncertainty
    TGraphAsymmErrors *gStat = getTGraphOfResolution("MaxLike","Statistic");
    fittedRes_ = new TF1("FittedResolution::fittedRes_",
			 "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
			 ptMin_,ptMax_);
    for(int i = 0; i < 3; ++i) {
      fittedRes_->SetParameter(i,par_->trueGaussResPar(i));
    }
    if( !par_->hasMCTruthBins() ) fittedRes_->FixParameter(0,par_->trueGaussResPar(0));
    fittedRes_->SetLineColor(2);
    fittedRes_->SetLineWidth(lineWidth_);
    if( par_->fitExtrapolatedSigma() ) gStat->Fit(fittedRes_,"0R");
    delete gStat;
  }

  // -------------------------------------------------------------------------------------
  FittedResolution::~FittedResolution() {
    delete trueRes_;
    delete fittedRes_;
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotPtAsymmetryBins() const {
    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotPtAsymmetryBins(): Entering method" << std::endl;
    }

    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());
      for(int c = 0; c < par_->nPt3Cuts(); ++c) {
	double xMin = -0.4;
	double xMax = 0.4;

	// Pt asymmetry distribution
	TH1 *hPtAsym = (*it)->getHistPtAsym(c,"PtAsymmetryBins:hPtAsym");
	util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","events",true);
	hPtAsym->GetXaxis()->SetRangeUser(xMin,xMax);
	util::HistOps::setYRange(hPtAsym,4);
	hPtAsym->SetMarkerStyle(20);

	// Fitted pt asymmetry
	TF1 *fPtAsym = new TF1("PtAsymmetryBins:fPtAsym","gaus",hPtAsym->GetMean()-2.*hPtAsym->GetRMS(),hPtAsym->GetMean()+2.*hPtAsym->GetRMS());
	fPtAsym->SetParameter(1,0.);
	fPtAsym->SetParameter(2,(*it)->fittedAsym(c)/sqrt(2.));
	fPtAsym->SetParameter(0,1./sqrt(2.*M_PI)/fPtAsym->GetParameter(2));	
	fPtAsym->SetLineWidth(lineWidth_);
	fPtAsym->SetLineColor(4);
	
	// Pt asymmetry from MaxLike
	TF1 *fMaxLike = new TF1("PtAsymmetryBins:fMaxLike","gaus",xMin,xMax);
	fMaxLike->SetParameter(1,0.);
	fMaxLike->SetParameter(2,(*it)->fittedValue(0,c)/sqrt(2.));
	fMaxLike->SetParameter(0,1./sqrt(2.*M_PI)/fMaxLike->GetParameter(2));	
	fMaxLike->SetLineWidth(lineWidth_);
	fMaxLike->SetLineColor(2);

	// Labels
	TPaveText *label = util::LabelFactory::createPaveText(3,-0.45);
	label->AddText(par_->labelLumi()+",  "+par_->labelEtaBin());
	label->AddText(par_->labelPtBin(ptBin));
	label->AddText(par_->labelPt3Cut(c));

	TLegend *leg = util::LabelFactory::createLegendCol(3,0.55);
	leg->AddEntry(hPtAsym,"Measurement","P");
	leg->AddEntry(fPtAsym,"Gaussian Fit","L");
	leg->AddEntry(fMaxLike,"MaxLike Prediction","L");

	// Draw
	TCanvas *can = new TCanvas("PtAsymmetry","PtAsymmetry ("+util::toTString(ptBin)+" "+util::toTString(c)+")",500,500);
	can->cd();
	hPtAsym->Draw("PE1");
	fPtAsym->Draw("same");
	fMaxLike->Draw("same");
	label->Draw("same");
	leg->Draw("same");
	can->SaveAs(par_->outNamePrefix()+"PtAsymmetry_PtBin"+util::toTString(ptBin)+"_Pt3Cut"+util::toTString(c)+".eps","eps");

	// Clean up
	delete hPtAsym;
	delete fMaxLike;
	delete fPtAsym;
	delete label;
	delete leg;
	delete can;
      }
    }

    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotPtAsymmetryBins(): Leaving method" << std::endl;
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

	// Create Canvas
	TCanvas *can = new TCanvas("PlotExtrapolation","Extrapolation ("+util::toTString(ptBin)+")",500,500);
	can->cd();

	// Draw a frame
	TH1 *h = (*it)->getFrameOfVariation(parIdx,"FrameExtrapolation");
	for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
	  h->SetBinContent(bin,par_->trueGaussSigma((*it)->meanPt()));
	}
	h->GetYaxis()->SetRangeUser(0.9*par_->trueGaussSigma((*it)->meanPt()),1.5*par_->trueGaussSigma((*it)->meanPt()));
	h->SetLineStyle(2);
	h->Draw();

	// Draw graph and extrapolation
	TGraphAsymmErrors *g = (*it)->getTGraphOfVariation(parIdx);
	g->Draw("PE1same");
	TF1 *f = (*it)->getTF1OfVariation(parIdx,"FitExtrapolation");
	f->SetLineWidth(lineWidth_);
	f->Draw("same");
	g->Draw("PE1same");

	TGraphAsymmErrors *gAsym = (*it)->getTGraphOfVariationAsym();
	gAsym->SetMarkerStyle(24);
	TF1 *fAsym = (*it)->getTF1OfVariationAsym("FitExtrapolation");
	fAsym->SetLineWidth(lineWidth_);
	fAsym->SetLineStyle(2);
	
	// Draw label
	TPaveText *txt = util::LabelFactory::createPaveText(2,-0.5);
	txt->AddText(par_->labelLumi()+", "+par_->labelEtaBin());
	txt->AddText(par_->labelPtBin(ptBin));
	txt->Draw("same");

	TLegend *leg = util::LabelFactory::createLegendCol(3,0.47);
	leg->AddEntry(g,"Maximum likelihood","P");
	leg->AddEntry(gAsym,"p_{T} asymmetry","P");
	leg->AddEntry(h,"MC truth","L");

	// Write canvas to file
	TString name = par_->outNamePrefix();
	name += "ExtrapolatedPar";
	name += parIdx;
	name += "_PtBin";
	name += ptBin;
	name += ".eps";
	can->SaveAs(name,"eps");

	if( parIdx == 0 ) {
	  can->cd();
	  h->Draw();
	  g->Draw("PE1same");
	  gAsym->Draw("PE1same");
	  f->Draw("same");
	  fAsym->Draw("same");
	  txt->Draw("same");
	  leg->Draw("same");
	  name = par_->outNamePrefix();
	  name += "ExtrapolatedPar";
	  name += parIdx;
	  name += "CompAsym_PtBin";
	  name += ptBin;
	  name += ".eps";
	  can->SaveAs(name,"eps");
	}

	// Clean up
	delete h;
	delete f;
	delete fAsym;
	delete g;
	delete gAsym;
	delete txt;
	delete leg;
	delete can;
      } // End of loop over fitted parameters



      // Draw pt asymmetry for different pt3 cuts
      TCanvas *can = new TCanvas("PlotVariedAsym","Varied Pt Asymmetry",500,500);
      can->cd();
      std::vector<TH1*> hAsym;
      TLegend *leg = util::LabelFactory::createLegendCol(3,par_->pt3Bins() ? 0.47 : 0.35);
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
	util::HistOps::setAxisTitles(h,"p_{T} asymmetry","","events");
	h->SetMarkerStyle(20+2*i);
	h->SetMarkerColor(util::StyleSettings::color(i));
	h->SetLineColor(h->GetMarkerColor());
	h->GetXaxis()->SetRangeUser(-0.4,0.4);
	if( i == 0 ) {
	  util::HistOps::setYRange(h,3);
	  h->Draw("PE1");
	} else {
	  h->Draw("PE1same");
	}
	leg->AddEntry(h,par_->labelPt3Cut(pt3Idx),"P");
	hAsym.push_back(h);
      }
      
      // Draw label
      TPaveText *txt = util::LabelFactory::createPaveText(2,-0.5);
      txt->AddText(par_->labelLumi()+", "+par_->labelEtaBin());
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


    // Create Canvas
    TCanvas *can = new TCanvas("CanExtrapolatedResolution","Extrapolated Resolution",500,500);
    can->cd();


    // ----- Plot relative resolution sigma / pt -----
    // Create graphs with statistical uncertainties
    TGraphAsymmErrors *gMaxLikeStat = getTGraphOfResolution("MaxLike","Statistic");
    TGraphAsymmErrors *gPtAsymStat = getTGraphOfResolution("PtAsym","Statistic");

    // Create a frame
    TH1 *h = new TH1D("FrameExtrapolatedResolution","",1000,ptMin_,ptMax_);
    h->SetXTitle("p^{ref}_{T} (GeV)");
    h->SetYTitle("#sigma / p^{ref}_{T}");
    h->SetNdivisions(510);
    double min = 0.7*(*std::min_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
    double max = 1.4*(*std::max_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
    h->GetYaxis()->SetRangeUser(min,max);

    // Create labels
    TPaveText *txt = 0;
    if( par_->extendedLegend() ) {
      txt = util::LabelFactory::createPaveText(2);
      txt->AddText("PYTHIA, #sqrt{s} = 7 TeV, "+par_->labelLumi());
      txt->AddText("Anti-k_{T} d = 0.5 jets, "+par_->labelEtaBin());
    } else {
      txt = util::LabelFactory::createPaveText(1);
      txt->AddText(par_->labelLumi()+", "+par_->labelEtaBin());
    }

    int nLegEntries = 2;

    TLegend *legMaxLike = util::LabelFactory::createLegendWithOffset(nLegEntries,txt->GetSize());
    legMaxLike->AddEntry(gMaxLikeStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    legMaxLike->AddEntry(trueRes_,"MC truth resolution","L");

    TLegend *legPtAsym = util::LabelFactory::createLegendWithOffset(nLegEntries,txt->GetSize());
    legPtAsym->AddEntry(gMaxLikeStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    legPtAsym->AddEntry(trueRes_,"MC truth resolution","L");

    TLegend *legComp = util::LabelFactory::createLegendWithOffset(nLegEntries,txt->GetSize());
    legComp->AddEntry(gMaxLikeStat,"MaxLike: Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    legComp->AddEntry(gPtAsymStat,"PtAsym: Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    legComp->AddEntry(trueRes_,"MC truth resolution","L");

    // Draw MaxLike results
    h->Draw();
    trueRes_->Draw("same");
    gMaxLikeStat->Draw("PE1same");
    txt->Draw("same");
    legMaxLike->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraReso.eps","eps");

    // Draw PtAsym results
    h->Draw();
    trueRes_->Draw("same");
    gPtAsymStat->Draw("PE1same");
    txt->Draw("same");
    legPtAsym->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoPtAsym.eps","eps");

    // Draw comparison
    h->Draw();
    trueRes_->Draw("same");
    gMaxLikeStat->Draw("PE1same");
    gPtAsymStat->Draw("PE1same");
    txt->Draw("same");
    legComp->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoComp.eps","eps");



    // ----- Plot relative deviation (sigma(fit)-sigma(true) ) / sigma(true)  -----

    // Fit MC truth response distributions
    std::vector< std::vector<double> > mcTruthReso(2);
    std::vector< std::vector<double> > mcTruthResoErr(mcTruthReso.size());
    if( par_->hasMCClosure() ) {
      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  	  it != ptBins_.end(); it++) {
  	TH1 *hMCRes = (*it)->getHistMCRes("PlotMCClosureMCRes");
	double min = hMCRes->GetMean() - 3.*hMCRes->GetRMS();
  	double max = hMCRes->GetMean() + 3.*hMCRes->GetRMS();
	hMCRes->Fit("gaus","I0Q","",min,max);
	TF1 *f = hMCRes->GetFunction("gaus");
	double reso = f->GetParameter(2);
	double resoErr = f->GetParError(2);
	min = f->GetParameter(1) - 2.*f->GetParameter(2);
	max = f->GetParameter(1) + 2.*f->GetParameter(2);
  	for(size_t i = 0; i < mcTruthReso.size(); ++i) {
	  TF1 *fit = 0;
  	  if( i == 0 ) {
  	    fit = new TF1("Fit","gaus",min,max);
  	    fit->SetParameter(0,5.);
  	    fit->SetParameter(1,1.);
  	    fit->SetParameter(2,0.1);
  	  } else if( i == 1 ) {
  	    fit = new TF1("Fit","gaus",min,max);
  	    fit->SetParameter(0,5.);
  	    fit->FixParameter(1,1.);
  	    fit->SetParameter(2,0.1);
  	  }
	  if( hMCRes->Fit(fit,"I0QBR") == 0 ) {
	    reso = fit->GetParameter(2);
	    resoErr = fit->GetParError(2);
	  }
	  mcTruthReso[i].push_back(reso);
	  mcTruthResoErr[i].push_back(resoErr);
	  delete fit;
  	}
	delete hMCRes;
      }
    }

    // Create ratio graphs with statistical uncertainties
    TGraphAsymmErrors *gMaxLikeRatioStat = getTGraphOfResolution("MaxLike","Statistic");
    std::vector<TGraphAsymmErrors*> gRatioMCClosure(mcTruthReso.size());
    for(int i = 0; i < gMaxLikeRatioStat->GetN(); ++i) {
      double x = gMaxLikeRatioStat->GetX()[i];
      double y = gMaxLikeRatioStat->GetY()[i];
      double yTrue = trueRes_->Eval(x);
      double exh = gMaxLikeRatioStat->GetEXhigh()[i];
      double exl = gMaxLikeRatioStat->GetEXlow()[i];
      double eyh = gMaxLikeRatioStat->GetEYhigh()[i];
      double eyl = gMaxLikeRatioStat->GetEYlow()[i];	
      gMaxLikeRatioStat->SetPoint(i,x,y/yTrue);
      gMaxLikeRatioStat->SetPointError(i,exl,exh,eyl/yTrue,eyh/yTrue);

      for(size_t t = 0; t < gRatioMCClosure.size(); ++t) {
	if( i == 0 ) gRatioMCClosure[t] = getTGraphOfResolution("MaxLike","Statistic");
	y = mcTruthReso[t][i] / yTrue;
	eyh = mcTruthResoErr[t][i] / yTrue;
	gRatioMCClosure[t]->SetPoint(i,x,y);
	gRatioMCClosure[t]->SetPointError(i,exl,exh,eyh,eyh);
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

    TGraphAsymmErrors *gPtAsymRatioStat = getTGraphOfResolution("PtAsym","Statistic");
    for(int i = 0; i < gPtAsymRatioStat->GetN(); ++i) {
      double x = gPtAsymRatioStat->GetX()[i];
      double y = gPtAsymRatioStat->GetY()[i];
      double yTrue = trueRes_->Eval(x);
      double exh = gPtAsymRatioStat->GetEXhigh()[i];
      double exl = gPtAsymRatioStat->GetEXlow()[i];
      double eyh = gPtAsymRatioStat->GetEYhigh()[i];
      double eyl = gPtAsymRatioStat->GetEYlow()[i];	
      gPtAsymRatioStat->SetPoint(i,x,y/yTrue);
      gPtAsymRatioStat->SetPointError(i,exl,exh,eyl/yTrue,eyh/yTrue);
    }


    // Adjust frame
    nLegEntries = 1;
    if( par_->fitRatio() ) nLegEntries++;
    if( par_->hasStartOffset() ) nLegEntries++;
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      h->SetBinContent(bin,1.);
    }
    h->SetLineStyle(trueRes_->GetLineStyle());
    h->SetLineColor(trueRes_->GetLineColor());
    h->GetYaxis()->SetTitle("#sigma_{fit} / #sigma_{MC}");
    h->GetYaxis()->SetRangeUser(0.65,1.45+0.8*nLegEntries*lineHeight_);


    // Recreate labels
    delete legMaxLike;
    delete legPtAsym;
    delete legComp;

    legMaxLike = util::LabelFactory::createLegendWithOffset(nLegEntries,txt->GetSize());
    legMaxLike->AddEntry(gMaxLikeRatioStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    if( par_->fitRatio() ) legMaxLike->AddEntry(lineFitRatio,"Mean fitted #bar{#sigma}","L");
    if( par_->hasStartOffset() ) legMaxLike->AddEntry(lineStartRes,"Resolution in spectrum","L");

    legPtAsym = util::LabelFactory::createLegendWithOffset(1,txt->GetSize());
    legPtAsym->AddEntry(gPtAsymRatioStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");

    legComp = util::LabelFactory::createLegendWithOffset(2,txt->GetSize());
    legComp->AddEntry(gMaxLikeRatioStat,"MaxLike: Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    legComp->AddEntry(gPtAsymRatioStat,"PtAsym: Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");


    // Draw MaxLike results
    h->Draw();
    if( par_->fitRatio() ) lineFitRatio->Draw("same");
    if( par_->hasStartOffset() ) lineStartRes->Draw("same");
    gMaxLikeRatioStat->Draw("PE1same");
    txt->Draw("same");
    legMaxLike->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoRatio.eps","eps");

    // Draw PtAsym results
    h->Draw();
    gPtAsymStat->Draw("PE1same");
    txt->Draw("same");
    legPtAsym->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoPtAsymRatio.eps","eps");

    // Draw comparison
    h->Draw();
    gMaxLikeRatioStat->Draw("PE1same");
    gPtAsymRatioStat->Draw("PE1same");
    txt->Draw("same");
    legComp->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoCompRatio.eps","eps");

    h->Draw();
    gMaxLikeRatioStat->Draw("PE1same");
    gPtAsymRatioStat->Draw("PE1same");
    for(size_t i = 0; i < gRatioMCClosure.size(); ++i) {
      gRatioMCClosure[i]->SetMarkerStyle(25+i);
      gRatioMCClosure[i]->SetMarkerColor(util::StyleSettings::color(1+i));
      gRatioMCClosure[i]->SetLineColor(util::StyleSettings::color(1+i));
      gRatioMCClosure[i]->Draw("PE1same");
    }
    txt->Draw("same");
    legComp->Draw("same");
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtraResoMCClosureRatio.eps","eps");


    // Clean up
    delete h;
    delete gMaxLikeStat;
    delete gPtAsymStat;
    delete gMaxLikeRatioStat;
    delete gPtAsymRatioStat;
    for(size_t i = 0; i < gRatioMCClosure.size(); ++i) {
      delete gRatioMCClosure[i];
    }
    delete lineFitRatio;
    delete lineStartRes;
    delete legMaxLike;
    delete legPtAsym;
    delete legComp;
    delete txt;
    delete can;

    if( par_->verbosity() >= 2 ) {
      std::cout << "FittedResolution::plotResolution(): Leaviing method" << std::endl;
    }
  }
  
  

  
  // -------------------------------------------------------------------------------------
  void FittedResolution::plotSpectra() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());
      
      // Create Canvas
      TCanvas *can = new TCanvas("PlotSpectrum","Spectrum ("+util::toTString(ptBin)+")",500,500);
      can->cd();
      
      // PtGen spectrum
      TH1 *hPtGen = (*it)->getHistPtGen("PlotPtGen");
      util::HistOps::setAxisTitles(hPtGen,"p^{"+par_->labelTruth()+"}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPtGen,5);
      hPtGen->SetMarkerStyle(24);

      // Pdf
      TH1 *hPdfPtTrue = (*it)->getHistPdfPtTrue("PlotPdfPtTrue");
      hPdfPtTrue->SetLineColor(2);
      hPdfPtTrue->SetLineWidth(lineWidth_);
      util::HistOps::setAxisTitles(hPdfPtTrue,"p^{"+par_->labelTruth()+"}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPdfPtTrue,5);

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(2);
      txt->AddText(par_->labelPtBin(ptBin));
      txt->AddText(par_->labelEtaBin()+",  "+par_->labelPt3Cut());
      TLegend *leg = util::LabelFactory::createLegendWithOffset(2,2.5*lineHeight_);
      leg->AddEntry(hPtGen,"MC truth: p^{"+par_->labelTruth()+"}_{T,1+2}","P");
      if( par_->fitMode() == FitModeMaxLikeFull ) leg->AddEntry(hPdfPtTrue,"Spectrum  #tilde{f}(p^{true}_{T})","L");

      // Plot
      hPtGen->Draw("PE1");
      if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
      txt->Draw("same");
      leg->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"Spectrum_PtBin"+util::toTString(ptBin)+".eps","eps");

      util::HistOps::setYRange(hPtGen,5,true);
      hPtGen->Draw("PE1");
      if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
      txt->Draw("same");
      leg->Draw("same");
      can->SetLogy();
      can->SaveAs(par_->outNamePrefix()+"Spectrum_PtBin"+util::toTString(ptBin)+"_Log.eps","eps");

      // Clean up
      delete hPtGen;
      delete hPdfPtTrue;
      delete txt;
      delete leg;
      delete can;
    }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotMCClosure() const {
    if( par_->hasMCClosure() ) {

      // Different fits of MC truth resolution
      std::vector<TString> mcTruthResoLabels;
      mcTruthResoLabels.push_back("#sigma,  <R> = free parameter");
      mcTruthResoLabels.push_back("#sigma,  <R> = 1");
      
      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  	  it != ptBins_.end(); it++) {
  	int ptBin = (it-ptBins_.begin());

  	TH1 *hMCRes = (*it)->getHistMCRes("PlotMCClosureMCRes");
  	util::HistOps::setAxisTitles(hMCRes,par_->xAxisTitleResponse(),"","jets",true);
  	hMCRes->SetMarkerStyle(20);
	
  	// Fit MC truth resolution
  	std::vector<TF1*> gaussFits(mcTruthResoLabels.size());
  	double min = hMCRes->GetMean() - 3.*hMCRes->GetRMS();
  	double max = hMCRes->GetMean() + 3.*hMCRes->GetRMS();
	hMCRes->Fit("gaus","I0Q","",min,max);
	TF1 *f = hMCRes->GetFunction("gaus");
	min = f->GetParameter(1) - 1.5*f->GetParameter(2);
	max = f->GetParameter(1) + 1.5*f->GetParameter(2);
  	for(size_t i = 0; i < gaussFits.size(); ++i) {
  	  TString name = "gaussFits"+util::toTString(i);
  	  if( i == 0 ) {
  	    gaussFits[i] = new TF1(name,"gaus",min,max);
  	    gaussFits[i]->SetParameter(0,5.);
  	    gaussFits[i]->SetParameter(1,1.);
  	    gaussFits[i]->SetParameter(2,0.1);
  	  } else if( i == 1 ) {
  	    gaussFits[i] = new TF1(name,"gaus",min,max);
  	    gaussFits[i]->SetParameter(0,5.);
  	    gaussFits[i]->FixParameter(1,1.);
  	    gaussFits[i]->SetParameter(2,0.1);
  	  }
  	  hMCRes->Fit(gaussFits[i],"I0QBR");
  	  gaussFits[i]->SetLineStyle(2);	  
  	  gaussFits[i]->SetLineColor(util::StyleSettings::color(2+i));
  	  gaussFits[i]->SetLineWidth(1);
  	}


  	// Fill histogram of extrapolated response
  	double rMin = 0.4;
  	double rMax = 1.6;
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
  	TPaveText *txt = util::LabelFactory::createPaveText(2,-0.45);
	txt->AddText(par_->labelLumi());
	txt->AddText(par_->labelEtaBin());
  	TLegend *leg = util::LabelFactory::createLegendCol(2,0.55);
  	leg->AddEntry(hMCRes,"MC truth, "+par_->labelPtBin(ptBin),"P");
  	leg->AddEntry(hFitRes,"Fit result, "+par_->labelPtBin(ptBin),"L");
  	TLegend *legFits = util::LabelFactory::createLegendColWithOffset(gaussFits.size(),0.55,txt->GetSize());
	for(size_t i = 0; i < gaussFits.size(); ++i) {
	  legFits->AddEntry(gaussFits[i],mcTruthResoLabels[i],"L");
	}

  	util::HistOps::setYRange(hFitRes,txt->GetSize()+3);
  	util::HistOps::setYRange(hFitResBins,txt->GetSize()+3);

  	// Create Canvas
  	TCanvas *can = new TCanvas("PlotMCClosure","MCClosure ("+util::toTString(ptBin)+")",500,500);
  	can->cd();

  	hFitRes->Draw("L");	
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
  	can->SaveAs(par_->outNamePrefix()+"MCClosure_PtBin"+util::toTString(ptBin)+".eps","eps");

  	hFitResBins->Draw("HIST");	
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
  	can->SaveAs(par_->outNamePrefix()+"MCClosureHist_PtBin"+util::toTString(ptBin)+".eps","eps");

  	hFitRes->Draw("L");	
  	for(size_t i = 0; i < gaussFits.size(); ++i) {
  	  gaussFits[i]->Draw("same");
  	}
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
	legFits->Draw("same");
  	can->SaveAs(par_->outNamePrefix()+"MCClosureFits_PtBin"+util::toTString(ptBin)+".eps","eps");
	
	hFitRes->GetXaxis()->SetRangeUser(0.,2.);
  	hFitRes->Draw("L");	
  	if( hFitResCore ) hFitResCore->Draw("Lsame");
  	hMCRes->Draw("PE1same");
  	txt->Draw("same");
  	leg->Draw("same");
  	gPad->SetLogy();
  	can->SaveAs(par_->outNamePrefix()+"MCClosureLog_PtBin"+util::toTString(ptBin)+".eps","eps");

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
  	for(size_t i = 0; i < gaussFits.size(); ++i) {
  	  delete gaussFits[i];
  	}
  	if( hFitResCore ) delete hFitResCore;
  	delete txt;
  	delete leg;
	delete legFits;
  	delete can;
      }
    } else {
      std::cerr << "No MCClosure response distribution available." << std::endl;
    }
  }





  //   void FittedResolution::plotResolutionBins() const {
  //     // Loop over ptbins
  //     for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  // 	it != ptBins_.end(); it++) {
  //       int bin = (it-ptBins_.begin());

  //       // Create Canvas
  //       TCanvas *can = new TCanvas("PlotResolution","PlotResolution",500,500);
  //       can->cd();

  //       // Draw MC truth resolution
  //       TH1 *hResGen = (*it)->getHistResGen("PlotResolution");
  //       hResGen->SetXTitle(par_->xAxisTitleResponse());
  //       hResGen->SetYTitle(par_->yAxisTitleResponse());
  //       hResGen->SetMarkerStyle(20);
  //       hResGen->GetXaxis()->SetRangeUser(0.4,1.6);
  //       hResGen->GetYaxis()->SetRangeUser(0.,(1.4+3.*lineHeight_)*hResGen->GetMaximum());
  //       hResGen->Draw("PE1");

  //       // Draw pdf
  //       TH1 *hPdfRes = (*it)->getHistPdfRes("PlotPdfResolution");
  //       hPdfRes->SetLineColor(2);
  //       hPdfRes->SetLineWidth(lineWidth_);
  //       hPdfRes->Draw("Lsame");
  //       hResGen->Draw("same");

  //       // Labels
  //       TPaveText *txt = util::LabelFactory::createPaveText(1);
  //       TString name = (*it)->ptMinStr();
  //       name += " < p^{";
  //       name += par_->labelMeas();
  //       name += "}_{T} < ";
  //       name += (*it)->ptMaxStr();
  //       name += " GeV";
  //       txt->AddText(name);
  //       txt->Draw("same");

  //       TLegend *leg = util::LabelFactory::createLegendWithOffset(2,1);
  //       leg->AddEntry(hResGen,"MC truth","P");
  //       char entry[100];
  //       sprintf(entry,"Fitted #sigma/p^{true}_{T}, p^{true}_{T} = %.1f GeV",(*it)->meanPt());
  //       leg->AddEntry(hPdfRes,entry,"L");
  //       leg->Draw("same");

  //       hResGen->Draw("PE1same");

  //       // Write Canvas to fiel
  //       name = par_->outNamePrefix();
  //       name += "Resolution_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");

  //       // Clean up
  //       delete hResGen;
  //       delete hPdfRes;
  //       delete txt;
  //       delete leg;
  //       delete can;
  //     }
  //   }


  //   void FittedResolution::plotPtAsymmetryBins() const {
  //     std::cout << "Plotting pt asymmetry distributions" << std::endl;

  //     TRandom3 *rand = new TRandom3(0);
  //     bool hasGaussCore = (par_->respFuncType() == ResponseFunction::CrystalBall) ||
  //       (par_->respFuncType() == ResponseFunction::TruncCrystalBall);

  //     // Loop over ptbins
  //     for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  // 	it != ptBins_.end(); it++) {
  //       int bin = (it-ptBins_.begin());

  //       // Measured asymmetry distribution
  //       // for default cut on pt3rel
  //       TH1 *hPtAsym = (*it)->getHistPtAsym("plotPtAsymmetryBins:PtAsymmetry");
  //       double xMin = -0.6;
  //       double xMax = 0.6;
  //       double yMin = 0.;
  //       double yMinLog = 3E-4;
  //       double yMax = (1.6+4.*lineHeight_)*hPtAsym->GetMaximum();
  //       double yMaxLog = 800.*hPtAsym->GetMaximum();
  //       util::HistOps::setAxisTitles(hPtAsym,"p_{T} asymmetry","","events",true);
  //       hPtAsym->SetMarkerStyle(20);
  //       hPtAsym->GetYaxis()->SetRangeUser(yMin,yMax);

  //       // Simulate asymmetry from fitted response
  //       int nDists = 1+2*par_->nFittedPars();
  //       std::vector<TH1*> hAsymPred(nDists);
  //       std::vector<TH1*> hAsymPredRatio(nDists);
  //       for(int d = 0; d < nDists; ++d) {
  // 	TString name = "PtAsymPredBin";
  // 	name += d;
  // 	hAsymPred[d] = static_cast<TH1D*>(hPtAsym->Clone(name));
  // 	hAsymPred[d]->SetMarkerStyle(0);
  // 	if( d == 0 ) hAsymPred[d]->SetLineColor(2);
  // 	else hAsymPred[d]->SetLineColor(util::StyleSettings::color(d%2+1));
  // 	hAsymPred[d]->SetLineWidth(lineWidth_);
  // 	hAsymPred[d]->GetXaxis()->SetRangeUser(xMin,xMax);

  // 	name = "hPtAsymPredRatio";
  // 	name += d;
  // 	hAsymPredRatio[d] = static_cast<TH1D*>(hAsymPred[d]->Clone(name));
  // 	hAsymPredRatio[d]->SetMarkerStyle(20);
  // 	if( d == 0 ) {
  // 	  hAsymPredRatio[d]->SetMarkerColor(2);
  // 	  hAsymPredRatio[d]->SetLineColor(2);
  // 	} else {
  // 	  hAsymPredRatio[d]->SetMarkerColor(hAsymPred[d]->GetLineColor());
  // 	  hAsymPredRatio[d]->SetLineColor(hAsymPred[d]->GetLineColor());
  // 	}
  // 	hAsymPredRatio[d]->Reset();
  // 	hAsymPredRatio[d]->SetYTitle("Prediction / Measurement");
  //       }
  //       TH1 *hAsymPredGauss = static_cast<TH1D*>(hAsymPred[0]->Clone("PtAsymPredGauss"));
  //       hAsymPredGauss->SetLineColor(4);
  //       hAsymPredGauss->SetLineStyle(2);
  //       TH1 *hAsymPredGaussRatio = static_cast<TH1D*>(hAsymPredRatio[0]->Clone("hPtAsymPredGaussRatio"));
  //       hAsymPredGaussRatio->SetMarkerColor(hAsymPredGauss->GetLineColor());
  //       hAsymPredGaussRatio->SetLineColor(hAsymPredGauss->GetLineColor());
  //       hAsymPredGaussRatio->SetMarkerStyle(24);

  //       // Fitted parameters
  //       std::vector<double> pars;
  //       pars.push_back(1.); // Correct scale assumed
  //       for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
  // 	pars.push_back((*it)->fittedValue(parIdx,par_->stdSelIdx()));
  //       }
  //       std::vector<double> parsGauss(2,0.);
  //       parsGauss[1] = pars[1]/sqrt(2.);

  //       // Variation factors
  //       std::vector<double> varFac(nDists,1.);
  //       for(int d = 1; d < nDists; ++d) {
  // 	if( d == 1 )      varFac[d] = 1./(1.-0.2);
  // 	else if( d == 2 ) varFac[d] = (1.-0.2);
  //  	else if( d == 3 ) varFac[d] = 1./(1.-0.5);
  //  	else if( d == 4 ) varFac[d] = (1.-0.5);
  // 	else if( d == 5 ) varFac[d] = 1./(1.-0.5);
  // 	else if( d == 6 ) varFac[d] = (1.-0.5);
  //       }

  //       // Fill asymmetry prediction
  //       // Loop over bins in asymmetry distributions
  //       for(int aBin = 1; aBin <= hAsymPred[0]->GetNbinsX(); ++aBin) {
  // 	double aMin = hAsymPred[0]->GetXaxis()->GetBinLowEdge(aBin);
  // 	double aMax = hAsymPred[0]->GetXaxis()->GetBinUpEdge(aBin);
  // 	int nIntBins = 50;
  // 	double deltaA = (aMax-aMin)/nIntBins;
  // 	std::vector<double> pdfA(nDists,0.);
  // 	double pdfAGauss = 0.;
  // 	// Integrate pdf over aBin
  // 	for(int intBin = 0; intBin < nIntBins; ++intBin) {
  // 	  double a = aMin + (intBin+0.5)*deltaA;

  // 	  // Asymmetry from chosen response
  // 	  pdfA[0] += par_->respFunc()->pdfAsymmetry(a,pars);
  // 	  // Asymmetry from varied response
  // 	  for(int d = 1; d < nDists; ++d) {
  // 	    int parIdx = (d-1)/2 + 1;
  // 	    std::vector<double> parsTmp = pars;
  // 	    parsTmp[parIdx] = varFac[d]*pars[parIdx];
  // 	    pdfA[d] += par_->respFunc()->pdfAsymmetry(a,parsTmp);
  // 	  }
  // 	  // Asymmetry assuming Gaussian response
  // 	  pdfAGauss += par_->respFunc()->pdfGauss(a,parsGauss);
  // 	}
  // 	for(int d = 0; d < nDists; ++d) {
  // 	  pdfA[d] *= deltaA/(aMax-aMin);
  // 	}
  // 	pdfAGauss *= deltaA/(aMax-aMin);

  // 	// Fill plots
  // 	for(int d = 0; d < nDists; ++d) {
  // 	  hAsymPred[d]->SetBinContent(aBin,pdfA[d]);
  // 	  hAsymPredRatio[d]->SetBinContent(aBin,pdfA[d]);
  // 	  hAsymPred[d]->SetBinError(aBin,0.);
  // 	  hAsymPredRatio[d]->SetBinError(aBin,0.);
  // 	}
  // 	hAsymPredGauss->SetBinContent(aBin,pdfAGauss);
  // 	hAsymPredGaussRatio->SetBinContent(aBin,pdfAGauss);
  // 	hAsymPredGauss->SetBinError(aBin,0.);
  // 	hAsymPredGaussRatio->SetBinError(aBin,0.);
  //       } // End of loop over asymmetry bins
  //       for(int d = 0; d < nDists; ++d) {
  // 	hAsymPredRatio[d]->Divide(hPtAsym);
  //       }
  //       hAsymPredGauss->Divide(hPtAsym);

  //       for(int d = 0; d < nDists; ++d) {
  // 	if( hAsymPred[d]->Integral("width") < 0.999 ) {
  // 	  std::cerr << "  WARNING: Int(" << bin << "," << d << ") = " << hAsymPred[d]->Integral("width") << std::endl;
  // 	}
  //       }



  //       // ----- Smearing ------------------------------------------------

  //       // Histograms
  //       TH1 *hAsymSmear = static_cast<TH1D*>(hAsymPred[0]->Clone("PtAsymSmear"));
  //       hAsymSmear->Reset();
  //       hAsymSmear->SetLineColor(4);
  //       hAsymSmear->SetLineWidth(lineWidth_);
  //       TH1 *hAsymSmearRatio = static_cast<TH1D*>(hAsymPred[0]->Clone("hPtAsymSmearRatio"));
  //       hAsymSmearRatio->Reset();
  //       hAsymSmearRatio->SetMarkerStyle(20);
  //       hAsymSmearRatio->SetMarkerColor(hAsymSmear->GetLineColor());
  //       hAsymSmearRatio->SetLineColor(hAsymSmear->GetLineColor());
  //       hAsymSmearRatio->SetYTitle("Prediction / Measurement");
      
  //       // Underlying truth spectrum
  //       TH1 *h = 0;
  //       TFile file("~/Kalibri/input/Spring10_TruthSpectrum_Eta0.root","READ");
  //       file.GetObject("hPtGen",h);
  //       if( !h ) {
  // 	std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
  // 	exit(1);
  //       }
  //       h->SetDirectory(0);
  //       file.Close();

  //       int minBin = h->FindBin(0.6*(*it)->ptMin());
  //       if( minBin < 1 ) minBin = 1;
  //       int maxBin = h->FindBin(1.6*(*it)->ptMax());
  //       if( maxBin > h->GetNbinsX() ) maxBin = h->GetNbinsX();
  //       TH1 *hGenSpectrum = new TH1D("plotAsymmetryBins:hGenSpectrum","",
  // 				   1+maxBin-minBin,h->GetXaxis()->GetBinLowEdge(minBin),
  // 				   h->GetXaxis()->GetBinUpEdge(maxBin));
  //       hGenSpectrum->SetLineColor(4);
  //       for(int xBin = 1; xBin <= hGenSpectrum->GetNbinsX(); ++xBin) {
  // 	hGenSpectrum->SetBinContent(xBin,h->GetBinContent(minBin+xBin-1));
  // 	hGenSpectrum->SetBinError(xBin,h->GetBinError(minBin+xBin-1));
  //       }
  //       delete h;

  //       double t = (*it)->meanPt();
  //       for(int n = 0; n < 500000; ++n) {
  // 	//double t = hGenSpectrum->GetRandom();
  // 	double m1 = t*par_->respFunc()->random(pars);
  // 	double m2 = t*par_->respFunc()->random(pars);
  // 	//	if( m1 > (*it)->ptMin() && m1 < (*it)->ptMax() ) {
  // 	  if( m1 + m2 > 0 ) {
  // 	    double x1 = m1;
  // 	    double x2 = m2;
  // 	    if( rand->Uniform() > 0.5 ) {
  // 	      x1 = m2;
  // 	      x2 = m1;
  // 	    }
  // 	    hAsymSmear->Fill((x1-x2)/(x1+x2));
  // 	    hAsymSmearRatio->Fill((x1-x2)/(x1+x2));
  // 	  }
  // // 	}
  // // 	if( m2 > (*it)->ptMin() && m2 < (*it)->ptMax() ) {
  // // 	  if( m1 + m2 > 0 ) {
  // // 	    double x1 = m1;
  // // 	    double x2 = m2;
  // // 	    if( rand->Uniform() > 0.5 ) {
  // // 	      x1 = m2;
  // // 	      x2 = m1;
  // // 	    }
  // // 	    hAsymSmear->Fill((x1-x2)/(x1+x2));
  // // 	    hAsymSmearRatio->Fill((x1-x2)/(x1+x2));
  // // 	  }
  // // 	}
  //       }
  //       if( hAsymSmear->Integral() ) hAsymSmear->Scale(1./hAsymSmear->Integral("width"));
  //       if( hAsymSmearRatio->Integral() ) hAsymSmearRatio->Scale(1./hAsymSmearRatio->Integral("width"));      
  //       hAsymSmearRatio->Divide(hPtAsym);



  //       // ----- Gaussian fit --------------------------------------------------------

  //       // Histograms
  //       TH1 *hAsymGaussFit = static_cast<TH1D*>(hAsymPred[0]->Clone("PtAsymGaussFit"));
  //       hAsymGaussFit->Reset();
  //       hAsymGaussFit->SetLineColor(4);
  //       hAsymGaussFit->SetLineStyle(2);
  //       hAsymGaussFit->SetLineWidth(lineWidth_);
  //       TH1 *hAsymGaussFitRatio = static_cast<TH1D*>(hAsymGaussFit->Clone("hPtAsymGaussFitRatio"));
  //       hAsymGaussFitRatio->SetMarkerStyle(20);
  //       hAsymGaussFitRatio->SetMarkerColor(hAsymGaussFit->GetLineColor());
  //       hAsymGaussFitRatio->SetYTitle("Prediction / Measurement");
      
  //       hPtAsym->Fit("gaus","I0Q");
  //       TF1 *fit = hPtAsym->GetFunction("gaus");
  //       for(int aBin = 1; aBin <= hAsymGaussFit->GetNbinsX(); ++aBin) {
  // 	double val = fit->Integral(hAsymGaussFit->GetXaxis()->GetBinLowEdge(aBin),
  // 				   hAsymGaussFit->GetXaxis()->GetBinUpEdge(aBin));
  // 	hAsymGaussFit->SetBinContent(aBin,val);
  // 	hAsymGaussFitRatio->SetBinContent(aBin,val);
  //       }
  //       hAsymGaussFit->Scale(1./hAsymGaussFit->GetBinWidth(1));
  //       hAsymGaussFitRatio->Scale(1./hAsymGaussFitRatio->GetBinWidth(1));
  //       hAsymGaussFitRatio->Divide(hPtAsym);



  //       // ----- Plotting -------------------------------------------------------------

  //       // Labels

  //       // For all
  //       TPaveText *txt = util::LabelFactory::createPaveText(2,-0.47);
  //       txt->AddText(par_->labelEtaBin()+",  "+par_->labelPt3Cut());
  //       txt->AddText(par_->labelPtBin(bin,0));

  //       // For prediction with fitted parameters only
  //       TLegend *legComp = util::LabelFactory::createLegendCol(4,0.5);
  //       legComp->AddEntry(hPtAsym,"Measurement","P");
  //       util::LabelFactory::addExtraLegLine(legComp,"#sqrt{V} = "+util::toTString(hPtAsym->GetRMS(),4));
  //       legComp->AddEntry(hAsymPred[0],"Prediction ("+par_->respFunc()->typeLabel()+")","L");
  //       util::LabelFactory::addExtraLegLine(legComp,"#sqrt{V} = "+util::toTString(hAsymPred[0]->GetRMS(),4));


  //       // Ratio with fitted parameters only
  //       TLegend *legRatio = util::LabelFactory::createLegendCol(1,0.5);
  //       legRatio->AddEntry(hAsymPredRatio[0],par_->respFunc()->typeLabel(),"P");

  //       // For prediction and ratio with varied parameters
  //       int nVars = par_->nFittedPars();
  //       TLegend *legMeas = util::LabelFactory::createLegendColWithOffset(1,-0.5,2);
  //       legMeas->AddEntry(hPtAsym,"Measurement","P");
  //       std::vector<TLegend*> legsPred(nVars);
  //       std::vector<TLegend*> legsRatio(nVars);
  //       for(int v = 0; v < nVars; ++v) {
  // 	double varUp = 100.*(1. - 1./varFac[1+2*v]);
  // 	double varDown = 100.*(1. - varFac[2+2*v]);
  // 	legsPred[v] = util::LabelFactory::createLegendCol(4,0.5);
  // 	util::LabelFactory::addExtraLegLine(legsPred[v],"Prediction ("+par_->respFunc()->typeLabel()+")");
  // 	legsPred[v]->AddEntry(hAsymPred[0],"Fitted parameters","L");
  // 	legsPred[v]->AddEntry(hAsymPred[1+2*v],par_->parLabel(v)+": + "+util::toTString(varUp,0)+"%","L");
  // 	legsPred[v]->AddEntry(hAsymPred[1+2*v],par_->parLabel(v)+": - "+util::toTString(varUp,0)+"%","L");

  // 	legsRatio[v] = util::LabelFactory::createLegendCol(4,0.5);
  // 	util::LabelFactory::addExtraLegLine(legsRatio[v],par_->respFunc()->typeLabel());
  // 	legsRatio[v]->AddEntry(hAsymPredRatio[0],"Fitted parameters","P");
  // 	legsRatio[v]->AddEntry(hAsymPredRatio[1+2*v],par_->parLabel(v)+": + "+util::toTString(varUp,0)+"%","L");
  // 	legsRatio[v]->AddEntry(hAsymPredRatio[1+2*v],par_->parLabel(v)+": - "+util::toTString(varUp,0)+"%","L");
  //       }

  //       // For prediction from smearing
  //       TLegend *legSmear = util::LabelFactory::createLegendCol(2,0.5);
  //       legSmear->AddEntry(hPtAsym,"Measurement","P");
  //       legSmear->AddEntry(hAsymSmear,"Smearing ("+par_->respFunc()->typeLabel()+")","L");

  //       // For Gaussian fit
  //       TLegend *legGaussFit = util::LabelFactory::createLegendCol(2,0.5);
  //       legGaussFit->AddEntry(hPtAsym,"Measurement","P");
  //       legGaussFit->AddEntry(hAsymGaussFit,"Gaussian fit","L");

  //       // For simple Gaussian prediction
  //       TLegend *legGauss = util::LabelFactory::createLegendCol(2,0.5);
  //       legGauss->AddEntry(hPtAsym,"Measurement","P");
  //       legGauss->AddEntry(hAsymPredGauss,"Gaussian prediction","L");


  //       // Create Canvas and write to file
  //       TCanvas *can = new TCanvas("PlotPtAsymmetry","PlotPtAsymmetry",500,500);
  //       can->cd();

  //       hGenSpectrum->GetYaxis()->SetRangeUser(3E-10,8.);
  //       hGenSpectrum->Draw("PE1");
  //       can->SetLogy(1);
  //       TString name = par_->outNamePrefix();
  //       name += "GenSpectrum1_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       //can->SaveAs(name,"eps");
  //       can->SetLogy(0);

  //       // Linear y axis scale
  //       hAsymPred[0]->GetYaxis()->SetRangeUser(yMin,yMax);
  //       hAsymPredGauss->GetYaxis()->SetRangeUser(yMin,yMax);
  //       hAsymPred[0]->Draw("H");
  //       hPtAsym->Draw("PE1same");
  //       txt->Draw("same");
  //       legComp->Draw("same");
  //       name = par_->outNamePrefix();
  //       name += "PtAsymmetry_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");

  //       hAsymPredGauss->GetYaxis()->SetRangeUser(yMin,yMax);
  //       hAsymPredGauss->Draw("H");
  //       hPtAsym->Draw("PE1same");
  //       txt->Draw("same");
  //       legGauss->Draw("same");
  //       name = par_->outNamePrefix();
  //       name += "PtAsymmetryGaussSimple_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");

  //       hAsymGaussFit->GetYaxis()->SetRangeUser(yMin,yMax);
  //       hAsymGaussFit->Draw("H");
  //       hPtAsym->Draw("PE1same");
  //       txt->Draw("same");
  //       legGaussFit->Draw("same");
  //       name = par_->outNamePrefix();
  //       name += "PtAsymmetryGaussFit_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");

  //       hAsymSmear->GetYaxis()->SetRangeUser(yMin,yMax);
  //       hAsymSmear->Draw("HIST E");
  //       hPtAsym->Draw("PE1same");
  //       txt->Draw("same");
  //       legSmear->Draw("same");
  //       name = par_->outNamePrefix();
  //       name += "PtAsymmetrySmear_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");
  //       for(int v = 0; v < nVars; ++v) {
  // 	double tmp1 = hAsymPred[2*v+1]->GetBinContent(hAsymPred[2*v+1]->GetMaximumBin());
  // 	double tmp2 = hAsymPred[2*v+2]->GetBinContent(hAsymPred[2*v+2]->GetMaximumBin());
  // 	if( tmp2 > tmp1 ) tmp1 = tmp2;
  // 	hAsymPred[0]->GetYaxis()->SetRangeUser(yMin,yMax*tmp1/hAsymPred[0]->GetBinContent(hAsymPred[0]->GetMaximumBin()));
  // 	hAsymPred[0]->SetLineColor(1);
  // 	hAsymPred[0]->Draw("H");
  // 	hAsymPred[2*v+1]->Draw("Hsame");
  // 	hAsymPred[2*v+2]->Draw("Hsame");
  // 	hAsymPred[0]->Draw("Hsame");
  // 	hPtAsym->Draw("PE1same");
  // 	txt->Draw("same");
  // 	legMeas->Draw("same");
  // 	legsPred[v]->Draw("same");
  // 	name = par_->outNamePrefix();
  // 	name += "PtAsymmetry_Var";
  // 	name += v;
  // 	name += "_PtBin";
  // 	name += bin;
  // 	name += ".eps";
  // 	can->SaveAs(name,"eps");
  //       }

  //       // Log y axis scale
  //       can->Clear();
  //       hPtAsym->GetYaxis()->SetRangeUser(yMinLog,yMaxLog);
  //       hAsymPred[0]->SetLineColor(2);
  //       for(int d = 0; d < nDists; ++d) {
  // 	hAsymPred[d]->GetYaxis()->SetRangeUser(yMinLog,yMaxLog);
  //       }
  //       hAsymPred[0]->Draw("H");
  //       hPtAsym->Draw("PE1same");
  //       txt->Draw("same");
  //       legComp->Draw("same");
  //       can->SetLogy();
  //       name = par_->outNamePrefix();
  //       name += "PtAsymmetryLog_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");
  //       for(int v = 0; v < nVars; ++v) {
  // 	hAsymPred[0]->SetLineColor(1);
  // 	hAsymPred[0]->Draw("H");
  // 	hAsymPred[2*v+1]->Draw("Hsame");
  // 	hAsymPred[2*v+2]->Draw("Hsame");
  // 	hAsymPred[0]->Draw("Hsame");
  // 	hPtAsym->Draw("PE1same");
  // 	txt->Draw("same");
  // 	legMeas->Draw("same");
  // 	legsPred[v]->Draw("same");
  // 	name = par_->outNamePrefix();
  // 	name += "PtAsymmetryLog_Var";
  // 	name += v;
  // 	name += "_PtBin";
  // 	name += bin;
  // 	name += ".eps";
  // 	can->SaveAs(name,"eps");
  //       }

  //       TH1 *hFrame = util::HistOps::createRatioFrame(hAsymPredRatio[0],"Prediction / Measurement",0.,2.5);
  //       hFrame->GetXaxis()->SetRangeUser(xMin,xMax);
  //       hFrame->Draw("H");
  //       hAsymPredRatio[0]->Draw("PE1same");
  //       txt->Draw("same");
  //       legRatio->Draw("same");
  //       can->SetLogy(0);
  //       name = par_->outNamePrefix();
  //       name += "PtAsymmetryRatio_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");
  //       for(int v = 0; v < nVars; ++v) {
  // 	hFrame->Draw("H");
  // 	hAsymPredRatio[2*v+1]->Draw("PE1same");
  // 	hAsymPredRatio[2*v+2]->Draw("PE1same");
  // 	hAsymPredRatio[0]->SetMarkerColor(1);
  // 	hAsymPredRatio[0]->SetLineColor(1);
  // 	hAsymPredRatio[0]->Draw("PE1same");
  // 	txt->Draw("same");
  // 	legsRatio[v]->Draw("same");
  // 	name = par_->outNamePrefix();
  // 	name += "PtAsymmetryRatio_Var";
  // 	name += v;
  // 	name += "_PtBin";
  // 	name += bin;
  // 	name += ".eps";
  // 	can->SaveAs(name,"eps");
  //       }

  //       // Clean up
  //       delete hPtAsym;
  //       for(int d = 0; d < nDists; ++d) {
  // 	delete hAsymPred[d];
  // 	delete hAsymPredRatio[d];
  //       }
  //       delete hAsymPredGauss;
  //       delete hAsymPredGaussRatio;
  //       delete hAsymSmear;
  //       delete hAsymSmearRatio;
  //       delete hGenSpectrum;
  //       delete hAsymGaussFit;
  //       delete hAsymGaussFitRatio;
  //       delete hFrame;
  //       delete txt;
  //       delete legComp;
  //       delete legRatio;
  //       delete legMeas;
  //       delete legSmear;
  //       delete legGaussFit;
  //       delete legGauss;
  //       for(int v = 0; v < nVars; ++v) {
  // 	delete legsPred[v];
  // 	delete legsRatio[v];
  //       }
  //       util::LabelFactory::deleteDummies();
  //       delete can;
  //     }
  //     delete rand;
  //   }


  //   void FittedResolution::plotPtAsymmetryAndResponseWidth() const {

  //     // Loop over ptbins
  //     for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  // 	it != ptBins_.end(); it++) {
  //       int ptBin = (it-ptBins_.begin());

  //       // Graph of resolution from pt asymmetry
  //       // for different cuts on pt3
  //       std::vector<double> x(par_->nPt3Cuts());
  //       std::vector<double> y(par_->nPt3Cuts());
  //       std::vector<double> xu(par_->nPt3Cuts(),0.);
  //       std::vector<double> yu(par_->nPt3Cuts());
  //       for(int c = 0; c < par_->nPt3Cuts(); ++c) {
  // 	x[c] = (*it)->pt3Max(c);
  // 	TH1 *h = (*it)->getHistPtAsym(c,"h");
  // 	h->Fit("gaus","I0Q");
  // 	TF1 *f = h->GetFunction("gaus");
  // 	y[c] = sqrt(2.)*f->GetParameter(2);
  // 	yu[c] = sqrt(2.)*f->GetParError(2);		  
  // 	delete h;
  //       }
  //       TGraphAsymmErrors *gResAsym = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
  // 							  &(xu.front()),&(xu.front()),
  // 							  &(yu.front()),&(yu.front()));
  //       gResAsym->SetMarkerStyle(24);


  //       // Graph of fitted resolution
  //       TGraphAsymmErrors *gResFit = (*it)->getTGraphOfVariation(0);

  // //       TF1 *f = (*it)->getTF1OfVariation(parIdx,"FitExtrapolation");
  // //       f->SetLineWidth(lineWidth_);


  //       // Frame and MC truth resolution
  //       TH1 *hFrame = (*it)->getFrameOfVariation(0,"FrameComparisonAsymFit");
  //       hFrame->GetYaxis()->SetRangeUser(hFrame->GetMinimum(),(1.+2*lineHeight_)*hFrame->GetMaximum());
  //       for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
  // 	hFrame->SetBinContent(bin,trueRes_->Eval((*it)->meanPt()));
  // 	hFrame->SetBinError(bin,0.);
  //       }
  //       hFrame->SetLineStyle(2);


  //       // Label
  //       TPaveText *txt = util::LabelFactory::createPaveText(2,-0.5);
  //       txt->AddText(par_->labelEtaBin());
  //       txt->AddText(par_->labelPtBin(ptBin,0));
      
  //       TLegend *leg = util::LabelFactory::createLegendCol(3,0.5);
  //       leg->AddEntry(gResFit,"#sigma(Response)","P");
  //       leg->AddEntry(gResAsym,"#sqrt{2}#upoint#sigma(Asymmetry)","P");
  //       leg->AddEntry(hFrame,"MC truth","L");

  //       TCanvas *can = new TCanvas("PlotPtAsymmetryAndResponse","PlotPtAsymmetryAndResponse",500,500);
  //       can->cd();
  //       hFrame->Draw("HIST");
  //       gResFit->Draw("PE1same");
  //       gResAsym->Draw("PE1same");
  //       txt->Draw("same");
  //       leg->Draw("same");
  //       TString name = par_->outNamePrefix();
  //       name += "PtAsymAndResp_PtBin";
  //       name += ptBin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");

  //       delete gResAsym;
  //       delete gResFit;
  //       delete hFrame;
  //       delete txt;
  //       delete leg;
  //       delete can;
  //     }
  //   }


  //   void FittedResolution::plotPtGenAsymmetry() const {
  //     std::vector<double> meanPt;
  //     std::vector<double> meanPtErr;
  //     std::vector<double> gaussAsym;
  //     std::vector<double> gaussAsymErr;
  //     std::vector<double> stdevAsym;
  //     std::vector<double> stdevAsymErr;

  //     // Loop over ptbins
  //     for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  // 	it != ptBins_.end(); it++) {
  //       int bin = (it-ptBins_.begin());

  //       // ptGen asymmetry distribution
  //       // for default cut on pt3rel
  //       TH1 *hPtGenAsym = (*it)->getHistPtGenAsym("plotPtGenAsymmetry:PtGenAsymmetry");
  //       hPtGenAsym->SetMarkerStyle(20);
  //       double xMin = -0.3;
  //       double xMax = 0.3;
  //       double yMin = 0.;
  //       double yMax = (1.4+2.*lineHeight_)*hPtGenAsym->GetMaximum();
  //       hPtGenAsym->GetXaxis()->SetRangeUser(xMin,xMax);
  //       hPtGenAsym->GetYaxis()->SetRangeUser(yMin,yMax);
  //       hPtGenAsym->GetYaxis()->SetTitle("1 / N  dN / dA");
  //       hPtGenAsym->GetXaxis()->SetTitle("p^{gen}_{T} asymmetry");

  //       // Gaussian fit vs standard deviation
  //       hPtGenAsym->Fit("gaus","0QI");
  //       TF1 *fitPtGenAsym = hPtGenAsym->GetFunction("gaus");
  //       fitPtGenAsym->SetLineWidth(lineWidth_);
  //       fitPtGenAsym->SetLineStyle(2);
  //       gaussAsym.push_back(fitPtGenAsym->GetParameter(2));
  //       gaussAsymErr.push_back(fitPtGenAsym->GetParError(2));
  //       stdevAsym.push_back(hPtGenAsym->GetRMS());
  //       stdevAsymErr.push_back(hPtGenAsym->GetRMSError());
  //       meanPt.push_back((*it)->meanPt());
  //       meanPtErr.push_back(0.);

  //       // Labels
  //       TPaveText *txt = util::LabelFactory::createPaveText(2);
  //       txt->AddText(par_->labelEtaBin()+",  "+par_->labelPt3Cut());
  //       txt->AddText(par_->labelPtBin(bin,0));

  //       TCanvas *can = new TCanvas("PlotPtGenAsymmetry","PlotPtGenAsymmetry",500,500);
  //       can->cd();
  //       hPtGenAsym->Draw("PE1");
  //       fitPtGenAsym->Draw("same");
  //       txt->Draw("same");
  //       TString name = par_->outNamePrefix();
  //       name += "PtGenAsymmetry_PtBin";
  //       name += bin;
  //       name += ".eps";
  //       can->SaveAs(name,"eps");

  //       delete hPtGenAsym;
  //       delete txt;
  //       delete can;
  //     }

  //     // Gaussian width vs standard deviation
  //     TGraphErrors *gGauss = new TGraphErrors(meanPt.size(),&(meanPt.front()),&(gaussAsym.front()),
  // 					    &(meanPtErr.front()),&(gaussAsymErr.front()));
  //     gGauss->SetMarkerStyle(24);
  //     TGraphErrors *gStdev = new TGraphErrors(meanPt.size(),&(meanPt.front()),&(stdevAsym.front()),
  // 					    &(meanPtErr.front()),&(stdevAsymErr.front()));
  //     gStdev->SetMarkerStyle(20);

  //     TH1 *hFrame = new TH1D("plotPtGenAsymmetry:hFrame",";p_{T} (GeV);p^{gen}_{T} asymmetry",
  // 			   100,ptMin(),ptMax());
  //     hFrame->GetYaxis()->SetRangeUser(0.,0.08);

  //     TCanvas *can = new TCanvas("PlotPtGenAsymmetry","PlotPtGenAsymmetry",500,500);
  //     can->cd();
  //     TPaveText *txt = util::LabelFactory::createPaveText(1);
  //     txt->AddText(par_->labelEtaBin()+",  "+par_->labelPt3Cut());
  //     TLegend *leg = util::LabelFactory::createLegendWithOffset(2,1);
  //     leg->AddEntry(gGauss,"#sigma Gauss fit","P");
  //     leg->AddEntry(gStdev,"Standard deviation","P");
  //     hFrame->Draw();
  //     gGauss->Draw("PE1same");
  //     gStdev->Draw("PE1same");
  //     txt->Draw("same");
  //     leg->Draw("same");
  //     can->SetLogx();
  //     TString name = par_->outNamePrefix();
  //     name += "PtGenAsymmetry.eps";
  //     can->SaveAs(name,"eps");
   
  //     delete txt;
  //     delete hFrame;
  //     delete gGauss;
  //     delete gStdev;
  //     delete can;
  //   }



  //     // Draw underlying truth spectra for different pt3 cuts
  //     if( par_->hasTruthSpectra() ) {
  //       std::vector<TH1*> hSpectra;
  //       TPaveText *txt = util::LabelFactory::createPaveText(1,0.47);
  //       txt->AddText(par_->labelEtaBin());
  //       TLegend *leg = util::LabelFactory::createLegendColWithOffset(par_->nPt3Cuts()-1,0.47,1);
  //       for(int cutIdx = 0; cutIdx < par_->nPt3Cuts(); ++cutIdx) {
  // 	KalibriFileParser *parser = new KalibriFileParser(par_->fileNameTruthSpectrum(cutIdx),par_->verbosity(),false);
  // 	TString name = "hTruthSpectrum";
  // 	name += cutIdx;
  // 	TH1 *h = parser->hist("hPtGen",name);
  // 	util::HistOps::normHist(h,"width");
  // 	util::HistOps::setAxisTitles(h,"p^{"+par_->labelTruth()+"}_{T}","GeV","jets",true);
  // 	h->SetLineColor(util::StyleSettings::color(cutIdx));      
  // 	h->SetLineWidth(1);
  // 	hSpectra.push_back(h);
  // 	leg->AddEntry(h,par_->labelPt3Cut(cutIdx),"L");
  // 	delete parser;
  //       }
  //       hSpectra[0]->GetYaxis()->SetRangeUser(3E-20,8);
  //       TCanvas *can = new TCanvas("canTruthSpectra","Truth Spectra",500,500);
  //       can->cd();
  //       hSpectra[0]->Draw("HISTC");
  //       for(size_t i = 1; i < hSpectra.size(); ++i) {
  // 	hSpectra[i]->Draw("HISTCsame");
  //       }
  //       txt->Draw("same");
  //       leg->Draw("same");
  //       can->SetLogy();
  //       can->SaveAs(par_->outNamePrefix()+"TruthSpectra.eps","eps");
  //       for(size_t i = 0; i < hSpectra.size(); ++i) {
  // 	delete hSpectra[i];
  //       }
  //       delete leg;
  //       delete txt;
  //       delete can;
  //     }    
  //   }


  //   void FittedResolution::plotSystematicUncertainties() const {
  // //     // Create graphs of systematic uncertainties
  // //     std::vector<double> pt(nPtBins());
  // //     for(int bin = 0; bin < nPtBins(); bin++) {
  // //       pt[bin] = meanPt(bin);
  // //     }
  // //     std::vector<TGraph*> gSystUp(nUncertSyst());
  // //     std::vector<TGraph*> gSystDown(nUncertSyst());
  // //     TLegend *leg = util::LabelFactory::createLegend(gSystUp.size());
  // //     // Loop over systematic uncertainties
  // //     for(int n = 0; n < nUncertSyst(); n++) {
  // //       std::vector<double> systUp(nPtBins());
  // //       std::vector<double> systDown(nPtBins());
  // //       for(int bin = 0; bin < nPtBins(); bin++) {
  // // 	const Uncertainty *uncertSyst = ptBins_[bin]->uncertSyst();
  // // 	systUp[bin] = uncertSyst->up(n)/relSigma(bin);
  // // 	systDown[bin] = uncertSyst->down(n)/relSigma(bin);
  // //       }
  // //       gSystUp[n] = new TGraph(pt.size(),&(pt.front()),&(systUp.front()));
  // //       gSystDown[n] = new TGraph(pt.size(),&(pt.front()),&(systDown.front()));
  // //       gSystUp[n]->SetLineColor(util::StyleSettings::color(n));
  // //       gSystDown[n]->SetLineColor(util::StyleSettings::color(n));
  // //       gSystUp[n]->SetLineWidth(lineWidth_);
  // //       gSystDown[n]->SetLineWidth(lineWidth_);

  // //       leg->AddEntry(gSystUp[n],ptBins_[0]->uncertSyst()->label(n),"L");
  // //     }

  // //     // Create Canvas
  // //     TCanvas *can = new TCanvas("CanSystematicUncertainties","Systematic uncertainties",500,500);
  // //     can->cd();

  // //     // Create frame histogram
  // //     double xMin = 0.8*meanPt(0);
  // //     double xMax = 1.1*meanPt(nPtBins()-1);
  // //     TH1 *h = new TH1D("FrameSystematicUncertainties",";p_{T} (GeV);#Delta#sigma / #sigma",
  // // 		       1000,xMin,xMax);
  // //     for(int bin = 1; bin < h->GetNbinsX(); bin++) {
  // //       h->SetBinContent(bin,0.);
  // //     }
  // //     h->GetYaxis()->SetRangeUser(-0.1,0.2);
  // //     h->SetLineStyle(2);
  // //     h->Draw();

  // //     // Plot systematic uncertainties
  // //     for(int n = 0; n < gSystUp.size(); n++) {
  // //       gSystUp[n]->Draw("Lsame");
  // //       //gSystDown[n]->Draw("Lsame");
  // //     }

  // //     // Add legend
  // //     leg->Draw("same");

  // //     // Write canvas to file
  // //     can->SaveAs(par_->outNamePrefix()+"SystematicUncertainties.eps","eps");

  // //     // Clean up
  // //     delete h;
  // //     for(int n = 0; n < gSystUp.size(); n++) {
  // //       delete gSystUp[n];
  // //       delete gSystDown[n];
  // //     }
  // //     delete leg;
  // //     delete can;
  //   }




  //   void FittedResolution::plotCrystalBallTest() const {
  //     if( par_->respFuncType() == ResponseFunction::CrystalBall ) {
  //       if( par_->hasMCClosure() ) {
  // 	// Loop over ptbins
  // 	for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  // 	    it != ptBins_.end(); it++) {
  // 	  int bin = (it-ptBins_.begin());
	
  // 	  double rMin = 0.;
  // 	  double rMax = 2.;
  // 	  double yMin = 6E-4;
  // 	  double yMax = 4E5;

  // 	  // Create Canvas
  // 	  TCanvas *can = new TCanvas("PlotCrystalBallTest","PlotCrystalBallTest",500,500);
  // 	  can->cd();

  // 	  // Draw MC truth resolution
  // 	  TH1 *hMCRes = (*it)->getHistMCRes("PlotCrystalBallTestMCRes");
  // 	  hMCRes->SetMarkerStyle(20);
  // 	  hMCRes->SetXTitle("Response R = p^{reco}_{T} / p^{particle}_{T}");
  // 	  hMCRes->SetYTitle("1 / N  dN / dR");
  // 	  hMCRes->GetXaxis()->SetRangeUser(rMin,rMax);
  // 	  hMCRes->GetYaxis()->SetRangeUser(yMin,yMax);
  // 	  hMCRes->Draw("PE1");
	
  // 	  // Fill histogram of extrapolated response
  // 	  TH1 *hExRes = new TH1D("PlotCrystalBallTestExRes","",500,rMin,rMax);
  // 	  hExRes->SetLineWidth(lineWidth_);
  // 	  hExRes->SetLineStyle(2);
  // 	  hExRes->SetLineColor(1);
  // 	  std::vector<double> pars;
  // 	  pars.push_back(1.);
  // 	  for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
  // 	    pars.push_back((*it)->extrapolatedValue(parIdx));
  // 	  }
  // 	  for(int rBin = 1; rBin <= hExRes->GetNbinsX(); ++rBin) {
  // 	    double res = hExRes->GetBinCenter(rBin);
  // 	    hExRes->SetBinContent(rBin,(*(par_->respFunc()))(res,pars));
  // 	  }
  // 	  hExRes->Draw("Lsame");	

  // 	  std::vector<TH1*> hCutVarRes;
  // 	  for(int c = 0; c < (*it)->nPt3Cuts(); ++c) {
  // 	    pars[2] = (*it)->fittedValue(1,c);
  // 	    pars[3] = (*it)->fittedValue(2,c);

  // 	    TH1 *h = static_cast<TH1D*>(hExRes->Clone("PlotCrystalBallTestH"));
  // 	    h->Reset();
  // 	    h->SetLineStyle(1);
  // 	    h->SetLineColor(util::StyleSettings::color(c));
  // 	    for(int rBin = 1; rBin <= hExRes->GetNbinsX(); ++rBin) {
  // 	      double res = h->GetBinCenter(rBin);
  // 	      h->SetBinContent(rBin,(*(par_->respFunc()))(res,pars));
  // 	    }
  // 	    h->Draw("Lsame");
  // 	    hCutVarRes.push_back(h);	  
  // 	  }
	
  // 	  gPad->RedrawAxis(); // Leads to slightly thicker font!!
		
  // 	  // Labels
  // 	  TPaveText *txt = util::LabelFactory::createPaveText(1,1.,0.04);
  // 	  txt->AddText((*it)->ptMinStr()+" < p_{T} < "+(*it)->ptMaxStr()+" GeV,  "+par_->labelEtaBin());
  // 	  txt->Draw("same");

  // 	  TLegend *leg = util::LabelFactory::createLegendWithOffset((*it)->nPt3Cuts()+2,0.05,0.038);
  // 	  leg->AddEntry(hMCRes,"MC truth","P");
  // 	  leg->AddEntry(hExRes,"Full extrapolation","L");
  // 	  for(int c = 0; c < (*it)->nPt3Cuts(); ++c) {
  // 	    leg->AddEntry(hCutVarRes[c],par_->labelPt3Cut(c),"L");
  // 	  }
  // 	  leg->Draw("same");
	
  // 	  gPad->SetLogy();
	
  // 	  // Write Canvas to fiel
  // 	  TString name = par_->outNamePrefix();
  // 	  name += "CrystalBallTest_PtBin";
  // 	  name += bin;
  // 	  name += ".eps";
  // 	  can->SaveAs(name,"eps");
	
  // 	  // Clean up
  // 	  delete hMCRes;
  // 	  delete hExRes;
  // 	  for(std::vector<TH1*>::iterator it = hCutVarRes.begin();
  // 	      it != hCutVarRes.end(); ++it) {
  // 	    delete *it;
  // 	  }
  // 	  delete txt;
  // 	  delete leg;
  // 	  delete can;
  // 	}
  //       } else {
  // 	std::cerr << "No MCClosure response distribution available." << std::endl;
  //       }
  //     }
  //   }



  // -------------------------------------------------------------------------------------
  void FittedResolution::print() const {
    std::cout << std::endl;
    for(int k = 0; k < fittedRes_->GetNpar(); ++k) {
      std::cout << " & $" << fittedRes_->GetParameter(k) << " \\pm " << fittedRes_->GetParError(k) << "$" << std::flush;
    }
    std::cout << " \\\\" << std::endl;
      

    std::cout << "\n\nFITTED RESOLUTION\n";
    // Loop over ptbins
    for(size_t bin = 0; bin < nPtBins(); bin++) {
      std::cout << bin << ": " << meanPt(bin) << std::flush;
      std::cout << " (" << ptMin(bin) << " - " << ptMax(bin) << "): " << std::flush;
      for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
	if( parIdx == 0 ) std::cout << extrapolatedValue(bin,parIdx)*meanPt(bin) << "," << std::flush;
	std::cout << " " << extrapolatedValue(bin,parIdx) << " (" << uncertStat(bin,parIdx) << ")" << std::flush;
	if( parIdx == 0 ) std::cout << ", " << trueRes_->Eval(meanPt(bin)) << std::flush;
	if( parIdx < par_->nFittedPars()-1 ) std::cout << " |" << std::flush;
	else std::cout << std::endl;
      }
    } 
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

    

  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors *FittedResolution::getTGraphOfResolution(const TString &method, const TString &uncertainty) const {
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
	y.push_back(ptBins_[i]->extrapolatedValue(0));
	eyDown.push_back(ptBins_[i]->uncertStatDown(0));
	eyUp.push_back(ptBins_[i]->uncertStatUp(0));
      } else if( method == "PtAsym" ) {
	x.push_back(ptBins_[i]->meanPtAsym());
	ex.push_back(ptBins_[i]->meanPtAsymUncert());
	y.push_back(ptBins_[i]->extrapolatedAsym());
	eyDown.push_back(ptBins_[i]->uncertDownAsym());
	eyUp.push_back(ptBins_[i]->uncertUpAsym());
      } else {
	std::cerr << "ERROR: FittedResolution::getTGraphOfResolution(): Unknown resolution method '" << method << "'" << std::endl;
	exit(-1);
      }
    }
    TGraphAsymmErrors *graph = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						     &(ex.front()),&(ex.front()),
						     &(eyDown.front()),&(eyUp.front()));
    graph->SetMarkerStyle(20);
    if( method == "PtAsym" ) graph->SetMarkerStyle(24);

    return graph;
  }



  // -------------------------------------------------------------------------------------
  double FittedResolution::gaussian(double *x, double *par) {
    double u = (x[0]-1.)/par[0];
    return par[0] > 0. ? exp(-0.5*u*u)/sqrt(2.*M_PI)/par[0] : 0.;
  }
}
  
