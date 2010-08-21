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

    if( par_->hasMCClosure() ) fitMCClosure();
  }

  // -------------------------------------------------------------------------------------
  FittedResolution::~FittedResolution() {
    delete trueRes_;
    delete fittedRes_;
    for(int i = 0; i < nMCClosureResFits(); ++i) {
      delete mcClosureGReso_[i];
      delete mcClosureGScale_[i];
      for(int p = 0; p < nPtBins(); ++p) {
	delete mcClosureFits_[i][p];
      }
    }
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
	util::HistOps::setYRange(hPtAsym,4);
	hPtAsym->SetMarkerStyle(20);

	// Fitted pt asymmetry
	TF1 *fPtAsym = new TF1("PtAsymmetryBins:fPtAsym","gaus",hPtAsym->GetMean()-2.*hPtAsym->GetRMS(),hPtAsym->GetMean()+2.*hPtAsym->GetRMS());
	fPtAsym->SetParameter(1,0.);
	fPtAsym->SetParameter(2,resPtAsym[ptBin]);
	fPtAsym->SetParameter(0,1./sqrt(2.*M_PI)/fPtAsym->GetParameter(2));	
	fPtAsym->SetLineWidth(lineWidth_);
	fPtAsym->SetLineColor(4);
	
	// Pt asymmetry from MaxLike
	TF1 *fMaxLike = new TF1("PtAsymmetryBins:fMaxLike","gaus",xMin,xMax);
	fMaxLike->SetParameter(1,0.);
	fMaxLike->SetParameter(2,resMaxLike[ptBin]);
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
      } // End of loop over pt bins


      // Compare fitted widths
      TGraphErrors *gResMaxLike = new TGraphErrors(ptMeanMaxLike.size(),&(ptMeanMaxLike.front()),&(resMaxLike.front()),&(ptMeanMaxLikeErr.front()),&(resMaxLikeErr.front()));
      gResMaxLike->SetMarkerStyle(20);
      gResMaxLike->SetMarkerColor(2);
      gResMaxLike->SetLineColor(2);
      TGraphErrors *gResPtAsym = new TGraphErrors(ptMeanPtAsym.size(),&(ptMeanPtAsym.front()),&(resPtAsym.front()),&(ptMeanPtAsymErr.front()),&(resPtAsymErr.front()));
      gResPtAsym->SetMarkerStyle(24);
      gResPtAsym->SetMarkerColor(4);
      gResPtAsym->SetLineColor(4);
      TGraphErrors *gRatio = new TGraphErrors(ptMeanMaxLike.size(),&(ptMeanMaxLike.front()),&(resRatio.front()),&(ptMeanMaxLikeErr.front()),&(resRatioErr.front()));
      gRatio->SetMarkerStyle(20);

      TPaveText *label = util::LabelFactory::createPaveText(2,-0.5);
      label->AddText(par_->labelLumi()+",  "+par_->labelEtaBin());
      label->AddText(par_->labelPt3Cut(c));

      TLegend *leg = util::LabelFactory::createLegendCol(2,0.5);
      leg->AddEntry(gResMaxLike,"MaxLike","P");
      leg->AddEntry(gResPtAsym,"PtAsym","P");

      TH1 *hFrame = new TH1D("FittedResolutions",";p^{ref}_{T} (GeV);#sigma / p^{ref}_{T}",1000,0.9*ptMin(),1.1*ptMax());
      double yMin = *(std::min_element(gResMaxLike->GetY(),gResMaxLike->GetY()+gResMaxLike->GetN()));
      double yMax = *(std::max_element(gResMaxLike->GetY(),gResMaxLike->GetY()+gResMaxLike->GetN()));
      double tmp =  *(std::min_element(gResPtAsym->GetY(),gResPtAsym->GetY()+gResPtAsym->GetN()));
      if( tmp < yMin ) yMin = tmp;
      tmp =  *(std::max_element(gResPtAsym->GetY(),gResPtAsym->GetY()+gResPtAsym->GetN()));
      if( tmp > yMax ) yMax = tmp;
      hFrame->GetYaxis()->SetRangeUser(0.,2.5*yMax);
      
      TCanvas *can = new TCanvas("PtAsymmetry","Fitted Resolutions ("+util::toTString(c)+")",500,500);
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
      
      
      delete hFrame;
      delete hRatioFrame;
      delete gResMaxLike;
      delete gResPtAsym;
      delete gRatio;
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
    double min = 0.5*(*std::min_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
    double max = 1.6*(*std::max_element(gMaxLikeStat->GetY(),gMaxLikeStat->GetY()+gMaxLikeStat->GetN()));
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
    TGraphAsymmErrors *gMaxLikeRatioStat = getTGraphOfResolution("MaxLike","Statistic");
    TGraphAsymmErrors* gMaxLikeRatioMCClosure = getTGraphOfResolution("MaxLike","Statistic");
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

      if( par_->hasMCClosure() ) {
	y = y / mcClosureReso(0,i);
	eyh = util::ratioError(y,eyh,mcClosureReso(0,i),mcClosureResoErr(0,i));
	gMaxLikeRatioMCClosure->SetPoint(i,x,y);
	gMaxLikeRatioMCClosure->SetPointError(i,exl,exh,eyh,eyh);
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
    TGraphAsymmErrors *gPtAsymRatioMCClosure = getTGraphOfResolution("PtAsym","Statistic");
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

      if( par_->hasMCClosure() ) {
	y = y / mcClosureReso(0,i);
	eyh = util::ratioError(y,eyh,mcClosureReso(0,i),mcClosureResoErr(0,i));
	gPtAsymRatioMCClosure->SetPoint(i,x,y);
	gPtAsymRatioMCClosure->SetPointError(i,exl,exh,eyh,eyh);
      }	
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
    gPtAsymRatioStat->Draw("PE1same");
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

    if( par_->hasMCClosure() ) {
      h->GetYaxis()->SetTitle("#sigma_{fit} / #sigma_{MCClosure}");
      h->Draw();
      gMaxLikeRatioMCClosure->Draw("PE1same");
      gPtAsymRatioMCClosure->Draw("PE1same");
      txt->Draw("same");
      legComp->Draw("same");
      gPad->SetLogx();
      can->SaveAs(par_->outNamePrefix()+"ExtraResoMCClosureRatio.eps","eps");
    }


    // Clean up
    delete h;
    delete gMaxLikeStat;
    delete gPtAsymStat;
    delete gMaxLikeRatioStat;
    delete gPtAsymRatioStat;
    delete gMaxLikeRatioMCClosure;
    delete gPtAsymRatioMCClosure;
    delete lineFitRatio;
    delete lineStartRes;
    delete legMaxLike;
    delete legPtAsym;
    delete legComp;
    delete txt;
    delete can;

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
      util::HistOps::setAxisTitles(hPtGen,"p^{"+par_->labelTruth()+"}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPtGen,5);
      hPtGen->SetMarkerStyle(20);

      TH1 *hPtGenJet1 = (*it)->getHistPtGenJet1("PlotPtGenJet1");
      util::HistOps::setAxisTitles(hPtGenJet1,"p^{"+par_->labelTruth()+"}_{T}","GeV","jets",true);
      util::HistOps::setYRange(hPtGenJet1,5);
      hPtGenJet1->SetMarkerStyle(20);

      // PtAve spectrum
      TH1 *hPtAve = (*it)->getHistPtAve("PlotPtAve");
      util::HistOps::setAxisTitles(hPtAve,"p^{ave}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPtAve,5);
      hPtAve->SetMarkerStyle(20);

      // Pdf
      TH1 *hPdfPtTrue = (*it)->getHistPdfPtTrue("PlotPdfPtTrue");
      hPdfPtTrue->SetLineColor(2);
      hPdfPtTrue->SetLineWidth(lineWidth_);
      util::HistOps::setAxisTitles(hPdfPtTrue,"p^{"+par_->labelTruth()+"}_{T}","GeV","events",true);
      util::HistOps::setYRange(hPdfPtTrue,5);

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(2,-0.55);
      txt->AddText(par_->labelPtBin(ptBin));
      txt->AddText(par_->labelEtaBin()+",  "+par_->labelPt3Cut());
      TLegend *legPtGen = util::LabelFactory::createLegendCol(2,0.45);
      legPtGen->AddEntry(hPtGen,"MC truth: p^{"+par_->labelTruth()+"}_{T,1+2}","P");
      legPtGen->AddEntry(hPdfPtTrue,"Spectrum  #tilde{f}(p^{true}_{T})","L");
      TLegend *legPtGenJet1 = util::LabelFactory::createLegendCol(2,0.45);
      legPtGenJet1->AddEntry(hPtGenJet1,"MC truth: p^{"+par_->labelTruth()+"}_{T,1}","P");
      legPtGenJet1->AddEntry(hPdfPtTrue,"Spectrum  #tilde{f}(p^{true}_{T})","L");
      TLegend *legPtAve = util::LabelFactory::createLegendCol(2,0.45);
      legPtAve->AddEntry(hPtAve,"p^{ave}_{T}","P");
      legPtAve->AddEntry(hPdfPtTrue,"Spectrum  #tilde{f}(p^{true}_{T})","L");

      // Plot
      TCanvas *can = new TCanvas("PlotSpectrum","Spectrum ("+util::toTString(ptBin)+")",500,500);
      can->cd();

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

      hPtAve->Draw("PE1");
      if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
      txt->Draw("same");
      if( par_->fitMode() == FitModeMaxLikeFull ) legPtAve->Draw("same");
      can->SaveAs(par_->outNamePrefix()+"PtAveSpectrum_PtBin"+util::toTString(ptBin)+".eps","eps");

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
      util::HistOps::setYRange(hPtAve,5,true);
      hPtAve->Draw("PE1");
      if( par_->fitMode() == FitModeMaxLikeFull ) hPdfPtTrue->Draw("Lsame");
      txt->Draw("same");
      if( par_->fitMode() == FitModeMaxLikeFull ) legPtAve->Draw("same");
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
      delete legPtAve;
      delete can;
    }
  }



  // -------------------------------------------------------------------------------------
  void FittedResolution::plotMCClosure() const {
    if( par_->hasMCClosure() ) {

      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
  	  it != ptBins_.end(); it++) {
  	int ptBin = (it-ptBins_.begin());
  	TH1 *hMCRes = (*it)->getHistMCRes("PlotMCClosureMCRes");
  	util::HistOps::setAxisTitles(hMCRes,par_->xAxisTitleResponse(),"","jets",true);
  	hMCRes->SetMarkerStyle(20);

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
  	TLegend *legFits = util::LabelFactory::createLegendColWithOffset(nMCClosureResFits(),0.55,txt->GetSize());
	for(int i = 0; i < nMCClosureResFits(); ++i) {
	  legFits->AddEntry(mcClosureFits_[i][ptBin],mcClosureResoLabels_[i],"L");
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
  	for(int i = 0; i < nMCClosureResFits(); ++i) {
  	  mcClosureFits_[i][ptBin]->Draw("same");
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
  	if( hFitResCore ) delete hFitResCore;
  	delete txt;
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

    

  //! Fit MC closure response distributions with different
  //! constraints and create graphs of scale and resolution
  // -------------------------------------------------------------------------------------
  void FittedResolution::fitMCClosure() {
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
    }

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
  
