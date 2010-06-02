#include "FittedResolution.h"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h" 
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"

#include "ResponseFunction.h"

#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


namespace resolutionFit {
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
    trueRes_->SetLineStyle(2);

    // Fit extrapolated resolution using
    // statistic uncertainty
    TGraphAsymmErrors *gStat = getTGraphOfResolution("Statistic+MCTruth");
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


  FittedResolution::~FittedResolution() {
    delete trueRes_;
    delete fittedRes_;
  }


  void FittedResolution::plotExtrapolation() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int ptBin = (it-ptBins_.begin());

      // Loop over fitted parameters
      for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {

	// Create Canvas
	TCanvas *can = new TCanvas("PlotExtrapolation","PlotExtrapolation",500,500);
	can->cd();

	// Draw a frame
	TH1 *h = (*it)->getFrameOfVariation(parIdx,"FrameExtrapolation");
	h->GetYaxis()->SetRangeUser(h->GetMinimum(),(1.+2*lineHeight_)*h->GetMaximum());
	h->Draw();

	// Draw graph and extrapolation
	TGraphAsymmErrors *g = (*it)->getTGraphOfVariation(parIdx);
	g->Draw("PE1same");
	TF1 *f = (*it)->getTF1OfVariation(parIdx,"FitExtrapolation");
	f->SetLineWidth(lineWidth_);
	f->Draw("same");
	g->Draw("PE1same");

	// Draw label
	TPaveText *txt = util::LabelFactory::createPaveText(2,1.);
	txt->AddText(par_->labelEtaBin());
	txt->AddText(par_->labelPtBin(ptBin,0));
	txt->Draw("same");
      
	// Write canvas to file
	TString name = par_->outNamePrefix();
	name += "ExtrapolatedPar";
	name += parIdx;
	name += "_PtBin";
	name += ptBin;
	name += ".eps";
	can->SaveAs(name,"eps");

	// Clean up
	delete h;
	delete f;
	delete g;
	delete txt;
	delete can;
      }
    }
  }


  void FittedResolution::plotResolutionBins() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int bin = (it-ptBins_.begin());

      // Create Canvas
      TCanvas *can = new TCanvas("PlotResolution","PlotResolution",500,500);
      can->cd();

      // Draw MC truth resolution
      TH1 *hResGen = (*it)->getHistResGen("PlotResolution");
      hResGen->SetMarkerStyle(20);
      hResGen->GetXaxis()->SetRangeUser(0.4,1.6);
      hResGen->GetYaxis()->SetRangeUser(0.,(1.4+3.*lineHeight_)*hResGen->GetMaximum());
      hResGen->Draw("PE1");

      // Draw pdf
      TH1 *hPdfRes = (*it)->getHistPdfRes("PlotPdfResolution");
      hPdfRes->SetLineColor(2);
      hPdfRes->SetLineWidth(lineWidth_);
      hPdfRes->Draw("Lsame");
      hResGen->Draw("same");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1);
      TString name = (*it)->ptMinStr();
      name += " < p^{recoJet}_{T} < ";
      name += (*it)->ptMaxStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,1,1);
      leg->AddEntry(hResGen,"MC truth","P");
      char entry[100];
      sprintf(entry,"Fitted #sigma/p^{true}_{T}, p^{true}_{T} = %.1f GeV",(*it)->meanPt());
      leg->AddEntry(hPdfRes,entry,"L");
      leg->Draw("same");

      hResGen->Draw("PE1same");

      // Write Canvas to fiel
      name = par_->outNamePrefix();
      name += "Resolution_PtBin";
      name += bin;
      name += ".eps";
      can->SaveAs(name,"eps");

      // Clean up
      delete hResGen;
      delete hPdfRes;
      delete txt;
      delete leg;
      delete can;
    }
  }


  void FittedResolution::plotPtAsymmetryBins() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int bin = (it-ptBins_.begin());

      // Create Canvas
      TCanvas *can = new TCanvas("PlotPtAsymmetry","PlotPtAsymmetry",500,500);
      can->cd();

      // Draw MC truth resolution
      TH1 *hPtAsym = (*it)->getHistPtAsym("PlotPtAsymmetry");
      hPtAsym->SetMarkerStyle(20);
      hPtAsym->GetXaxis()->SetRangeUser(0.4,1.6);
      hPtAsym->GetYaxis()->SetRangeUser(0.,(1.4+3.*lineHeight_)*hPtAsym->GetMaximum());
      hPtAsym->GetYaxis()->SetTitle("1 / N  dN / dA");
      hPtAsym->Draw("PE1");

      // Draw pdf
      TH1 *hPdfPtAsym = (*it)->getHistPdfPtAsym("PlotPdfPtAsymmetry");
      hPdfPtAsym->SetLineColor(2);
      hPdfPtAsym->SetLineWidth(lineWidth_);
      hPdfPtAsym->Draw("Lsame");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1);
      TString name = (*it)->ptMinStr();
      name += " < p^{recoJet}_{T} < ";
      name += (*it)->ptMaxStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,1,1);
      leg->AddEntry(hPtAsym,"MC truth","P");
      char entry[100];
      sprintf(entry,"Fitted  #sqrt{2}#sigma/p^{true}_{T}, p^{true}_{T} = %.1f GeV",(*it)->meanPt());
      leg->AddEntry(hPdfPtAsym,entry,"L");
      leg->Draw("same");

      hPtAsym->Draw("PE1same");

      // Write Canvas to fiel
      name = par_->outNamePrefix();
      name += "PtAsymmetry_PtBin";
      name += bin;
      name += ".eps";
      can->SaveAs(name,"eps");

      // Clean up
      delete hPtAsym;
      delete hPdfPtAsym;
      delete txt;
      delete leg;
      delete can;
    }
  }


  void FittedResolution::plotResolution() const {
    // Create Canvas
    TCanvas *can = new TCanvas("CanExtrapolatedResolution","Extrapolated Resolution",500,500);
    can->cd();


    // ----- Plot relative resolution sigma / pt -----
    // Create graph with statistical uncertainties
    TGraphAsymmErrors *gStat = getTGraphOfResolution("Statistic");
    TGraphAsymmErrors *gPseudo = 0;
    if( par_->hasMCTruthBins() ) gPseudo = getTGraphOfResolution("MCTruth");

    // Draw a frame
    int nLegEntries = 2;
    if( par_->fitExtrapolatedSigma() ) nLegEntries++;
    if( gPseudo ) nLegEntries++;
    TH1 *h = new TH1D("FrameExtrapolatedResolution",";p_{T} (GeV);#sigma / p_{T}",
		      1000,ptMin_,ptMax_);
    h->SetNdivisions(510);
    double min = 0.7*(*std::min_element(gStat->GetY(),gStat->GetY()+gStat->GetN()));
    double max = 1.3*(*std::max_element(gStat->GetY(),gStat->GetY()+gStat->GetN()));
    if( gPseudo ) max = 1.8*(*std::max_element(gPseudo->GetY(),gPseudo->GetY()+gPseudo->GetN()));
    h->GetYaxis()->SetRangeUser(min,max);
    h->Draw();

    // Draw systematic uncertainty band
    //TGraphAsymmErrors *gSyst = getTGraphOfResolution("Systematic");
    //gSyst->Draw("E3same");

    // Draw true and fitted resolution functions
    trueRes_->Draw("same");
    if( par_->fitExtrapolatedSigma() ) fittedRes_->Draw("same");

    // Draw fitted mean resolution
    gStat->Draw("PE1same");
    if( gPseudo ) gPseudo->Draw("PE1same");

    // Add a legend
    TLegend *leg = util::LabelFactory::createLegend(nLegEntries);
    leg->AddEntry(gStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    //leg->AddEntry(gSyst,"Uncertainty from spectrum","F");
    if( gPseudo ) leg->AddEntry(gPseudo,"Pseudo #gamma+jet measurement","P");
    leg->AddEntry(trueRes_,"MC truth resolution","L");
    if( par_->fitExtrapolatedSigma() ) leg->AddEntry(fittedRes_,"Fit #sigma(p_{T}) to #bar{#sigma}","L");
    leg->Draw("same");

    // Write canvas to file
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtrapolatedResolution.eps","eps");


    // ----- Plot relative deviation (sigma(fit)-sigma(true) ) / sigma(true)  -----

    // Fitted resolution function
    TH1 *hTrueRes = trueRes_->GetHistogram();
    TH1 *hFittedRatio = fittedRes_->GetHistogram();
    hFittedRatio->Divide(hTrueRes);

    // Draw graph with statistical uncertainties
    TGraphAsymmErrors *gStatRatio = getTGraphOfResolution("Statistic");
    for(int i = 0; i < gStatRatio->GetN(); ++i) {
      double x = gStatRatio->GetX()[i];
      double y = gStatRatio->GetY()[i];
      double yTrue = trueRes_->Eval(x);
      double exh = gStatRatio->GetEXhigh()[i];
      double exl = gStatRatio->GetEXlow()[i];
      double eyh = gStatRatio->GetEYhigh()[i];
      double eyl = gStatRatio->GetEYlow()[i];	
      gStatRatio->SetPoint(i,x,y/yTrue);
      gStatRatio->SetPointError(i,exl,exh,eyl/yTrue,eyh/yTrue);
    }

    TF1 *lineFitRatio = new TF1("LineFitRatioExtraplolatedResolution","pol0",ptMin_,ptMax_);
    lineFitRatio->SetLineWidth(lineWidth_);
    lineFitRatio->SetLineStyle(2);
    lineFitRatio->SetLineColor(2);
    if( par_->fitRatio() ) {
      gStatRatio->Fit(lineFitRatio,"0QR");
      std::cout << "Fitted ratio " << lineFitRatio->GetParameter(0) << std::flush;
      std::cout << " +/- " << lineFitRatio->GetParError(0) << std::endl;
    }

    TF1 *lineStartRes = static_cast<TF1*>(lineFitRatio->Clone("LineStartResExtrapolatedResolution"));
    lineStartRes->SetLineColor(8);
    lineStartRes->SetParameter(0,par_->relStartOffset());

    TGraphAsymmErrors *gPseudoRatio = 0;
    if( gPseudo ) {
      gPseudoRatio = getTGraphOfResolution("MCTruth");
      for(int i = 0; i < gPseudoRatio->GetN(); ++i) {
	double x = gPseudoRatio->GetX()[i];
	double y = gPseudoRatio->GetY()[i];
	double yTrue = trueRes_->Eval(x);
	double exh = gPseudoRatio->GetEXhigh()[i];
	double exl = gPseudoRatio->GetEXlow()[i];
	double eyh = gPseudoRatio->GetEYhigh()[i];
	double eyl = gPseudoRatio->GetEYlow()[i];	
	gPseudoRatio->SetPoint(i,x,y/yTrue);
	gPseudoRatio->SetPointError(i,exl,exh,eyl/yTrue,eyh/yTrue);
      }
    }

    // Frame
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      h->SetBinContent(bin,1.);
    }
    h->SetLineStyle(2);
    h->SetLineColor(4);
    h->GetYaxis()->SetTitle("#sigma_{fit} / #sigma_{MC}");
    h->GetYaxis()->SetRangeUser(0.65,1.45+0.8*nLegEntries*lineHeight_);

    nLegEntries = 2;
    if( par_->fitExtrapolatedSigma() ) nLegEntries++;
    if( gPseudo ) nLegEntries++;
    if( par_->fitRatio() ) nLegEntries++;
    if( par_->hasStartOffset() ) nLegEntries++;
    
    max = *(std::max_element(gStatRatio->GetY(),gStatRatio->GetY()+gStat->GetN()));
    if( max > 1.4 ) h->GetYaxis()->SetRangeUser(0.65,1.75+1.1*nLegEntries*lineHeight_);
    h->Draw();
    if( par_->fitExtrapolatedSigma() ) hFittedRatio->Draw("Lsame");
    if( par_->fitRatio() ) lineFitRatio->Draw("same");
    if( par_->hasStartOffset() ) lineStartRes->Draw("same");
    gStatRatio->Draw("PE1same");
    if( gPseudoRatio ) gPseudoRatio->Draw("PE1same");

    // Add a legend
    delete leg;
    leg = util::LabelFactory::createLegend(nLegEntries);
    leg->AddEntry(gStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    //leg->AddEntry(gSyst,"Uncertainty from spectrum","F");
    if( gPseudo ) leg->AddEntry(gPseudo,"Pseudo #gamma+jet measurement","P");
    leg->AddEntry(trueRes_,"MC truth resolution","L");
    if( par_->fitRatio() ) leg->AddEntry(lineFitRatio,"Mean fitted #bar{#sigma}","L");
    if( par_->hasStartOffset() ) leg->AddEntry(lineStartRes,"Resolution in spectrum","L");
    if( par_->fitExtrapolatedSigma() ) leg->AddEntry(fittedRes_,"Fit #sigma(p_{T}) to #bar{#sigma}","L");
    leg->Draw("same");

    // Write canvas to file
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtrapolatedResolutionRatio.eps","eps");


    // ----- Plot quadratic difference sigma(fit) - sigma(MC) -----

    TGraphAsymmErrors *gStatDiff = getTGraphOfResolution("Statistic");
    for(int i = 0; i < gStatDiff->GetN(); ++i) {
      double x = gStatDiff->GetX()[i];
      double y = gStatDiff->GetY()[i];
      double yTrue = trueRes_->Eval(x);
      y = sqrt( y*y - yTrue*yTrue );
      double exh = gStatDiff->GetEXhigh()[i];
      double exl = gStatDiff->GetEXlow()[i];
//       double eyh = gStatDiff->GetEYhigh()[i];
//       double eyl = gStatDiff->GetEYlow()[i];	
      gStatDiff->SetPoint(i,x,y);
      gStatDiff->SetPointError(i,exl,exh,0.,0.);
    }
    gStatDiff->SetMarkerStyle(21);

    h->Reset();
    h->GetYaxis()->SetRangeUser(0.,0.2);
    h->Draw();
    gStatDiff->Draw("PE1same");
    gStat->Draw("PE1same");
    trueRes_->Draw("same");

    // Write canvas to file
    gPad->SetLogx();
    can->SaveAs(par_->outNamePrefix()+"ExtrapolatedResolutionDifference.eps","eps");



    // Clean up
    delete h;
    //    delete gSyst;
    delete gStat;
    if( gPseudo ) delete gPseudo;
    if( gPseudoRatio ) delete gPseudoRatio;
    delete gStatRatio;
    delete gStatDiff;
    delete lineFitRatio;
    delete lineStartRes;
    delete leg;
    delete can;
  }


  void FittedResolution::plotSpectra() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int bin = (it-ptBins_.begin());
      bool logY = par_->respFuncType() == ResponseFunction::CrystalBall ? true : false;

      // Create Canvas
      TCanvas *can = new TCanvas("PlotSpectrum","PlotSpectrum",500,500);
      can->cd();

      // Draw MC truth spectra
      TH1 *hPtGen = (*it)->getHistPtGen("PlotPtGen");

      double yMin = logY ? 6E-6 : 0.;
      double yMax = logY ? 8E1 : (1.6+2*lineHeight_)*hPtGen->GetMaximum();

      hPtGen->SetXTitle("p^{gen}_{T}  (GeV)");
      hPtGen->SetYTitle("1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
      hPtGen->GetYaxis()->SetRangeUser(yMin,yMax);
      hPtGen->SetMarkerStyle(24);

      TH1 *hPtGenJet1 = (*it)->getHistPtGenJet1("PlotPtGenJet1");
      hPtGenJet1->SetXTitle("p^{particle}_{T}  (GeV)");
      hPtGenJet1->SetYTitle("1 / N  dN / dp^{particle}_{T}  1 / (GeV)");
      hPtGenJet1->GetYaxis()->SetRangeUser(yMin,yMax);
      hPtGenJet1->SetMarkerStyle(20);

      // Draw pdf
      TH1 *hPdfPtTrue = (*it)->getHistPdfPtTrue("PlotPdfPtTrue");
      hPdfPtTrue->SetLineColor(2);
      hPdfPtTrue->SetLineWidth(lineWidth_);
      hPdfPtTrue->SetXTitle("p^{gen}_{T}  (GeV)");
      hPdfPtTrue->SetYTitle("1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
      hPdfPtTrue->GetYaxis()->SetRangeUser(yMin,yMax);
      hPdfPtTrue->Draw("L");
      hPtGenJet1->Draw("PE1same");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(2);
      txt->AddText(par_->labelPtBin(bin,0));
      txt->AddText(par_->labelEtaBin()+",  p^{rel}_{T,3} < 0.1");
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,1,2);
      leg->AddEntry(hPtGenJet1,"MC truth: p^{particle}_{T}","P");
      //leg->AddEntry(hPtGen,"MC truth: p^{particleJet}_{T,1+2}","P");
      leg->AddEntry(hPdfPtTrue,"Spectrum  #tilde{f}(p_{T})","L");
      leg->Draw("same");

      if( logY ) gPad->SetLogy();

      // Write Canvas to fiel
      TString name = par_->outNamePrefix();
      name += "Spectrum_PtBin";
      name += bin;
      name += ".eps";
      can->SaveAs(name,"eps");

      // Clean up
      delete hPtGen;
      delete hPtGenJet1;
      delete hPdfPtTrue;
      delete txt;
      delete leg;
      delete can;
    }
  }


  void FittedResolution::plotSystematicUncertainties() const {
//     // Create graphs of systematic uncertainties
//     std::vector<double> pt(nPtBins());
//     for(int bin = 0; bin < nPtBins(); bin++) {
//       pt[bin] = meanPt(bin);
//     }
//     std::vector<TGraph*> gSystUp(nUncertSyst());
//     std::vector<TGraph*> gSystDown(nUncertSyst());
//     TLegend *leg = util::LabelFactory::createLegend(gSystUp.size());
//     // Loop over systematic uncertainties
//     for(int n = 0; n < nUncertSyst(); n++) {
//       std::vector<double> systUp(nPtBins());
//       std::vector<double> systDown(nPtBins());
//       for(int bin = 0; bin < nPtBins(); bin++) {
// 	const Uncertainty *uncertSyst = ptBins_[bin]->uncertSyst();
// 	systUp[bin] = uncertSyst->up(n)/relSigma(bin);
// 	systDown[bin] = uncertSyst->down(n)/relSigma(bin);
//       }
//       gSystUp[n] = new TGraph(pt.size(),&(pt.front()),&(systUp.front()));
//       gSystDown[n] = new TGraph(pt.size(),&(pt.front()),&(systDown.front()));
//       gSystUp[n]->SetLineColor(util::StyleSettings::color(n));
//       gSystDown[n]->SetLineColor(util::StyleSettings::color(n));
//       gSystUp[n]->SetLineWidth(lineWidth_);
//       gSystDown[n]->SetLineWidth(lineWidth_);

//       leg->AddEntry(gSystUp[n],ptBins_[0]->uncertSyst()->label(n),"L");
//     }

//     // Create Canvas
//     TCanvas *can = new TCanvas("CanSystematicUncertainties","Systematic uncertainties",500,500);
//     can->cd();

//     // Create frame histogram
//     double xMin = 0.8*meanPt(0);
//     double xMax = 1.1*meanPt(nPtBins()-1);
//     TH1 *h = new TH1D("FrameSystematicUncertainties",";p_{T} (GeV);#Delta#sigma / #sigma",
// 		       1000,xMin,xMax);
//     for(int bin = 1; bin < h->GetNbinsX(); bin++) {
//       h->SetBinContent(bin,0.);
//     }
//     h->GetYaxis()->SetRangeUser(-0.1,0.2);
//     h->SetLineStyle(2);
//     h->Draw();

//     // Plot systematic uncertainties
//     for(int n = 0; n < gSystUp.size(); n++) {
//       gSystUp[n]->Draw("Lsame");
//       //gSystDown[n]->Draw("Lsame");
//     }

//     // Add legend
//     leg->Draw("same");

//     // Write canvas to file
//     can->SaveAs(par_->outNamePrefix()+"SystematicUncertainties.eps","eps");

//     // Clean up
//     delete h;
//     for(int n = 0; n < gSystUp.size(); n++) {
//       delete gSystUp[n];
//       delete gSystDown[n];
//     }
//     delete leg;
//     delete can;
  }


  void FittedResolution::plotMCClosure() const {
    if( par_->hasMCClosure() ) {
      // Loop over ptbins
      for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	  it != ptBins_.end(); it++) {
	int bin = (it-ptBins_.begin());
	
	bool logY = par_->respFuncType() == ResponseFunction::CrystalBall ? true : false;
	double rMin = 0.;
	double rMax = 2.;
	double yMin = logY ? 6E-4 : 0.;
	double yMax = logY ? 4E3 : 12.;

	// Create Canvas
	TCanvas *can = new TCanvas("PlotMCClosure","PlotMCClosure",500,500);
	can->cd();

	// Draw MC truth resolution
	TH1 *hMCRes = (*it)->getHistMCRes("PlotMCClosureMCRes");
	hMCRes->SetMarkerStyle(20);
	hMCRes->SetXTitle("Response R = p^{reco}_{T} / p^{particle}_{T}");
	hMCRes->SetYTitle("1 / N  dN / dR");
	hMCRes->GetXaxis()->SetRangeUser(rMin,rMax);
	hMCRes->GetYaxis()->SetRangeUser(yMin,yMax);
	hMCRes->Draw("PE1");
	
	// Fill histogram of extrapolated response
	TH1 *hFitRes = new TH1D("PlotMCClosureFitRes","",500,rMin,rMax);
	hFitRes->SetLineWidth(lineWidth_);
	hFitRes->SetLineColor(2);
	std::vector<double> pars;
	for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
	  pars.push_back((*it)->extrapolatedValue(parIdx));
	}
	hFitRes->SetXTitle("Response R = p^{reco}_{T} / p^{particle}_{T}");
	hFitRes->SetYTitle("1 / N  dN / dR");
 	hFitRes->GetXaxis()->SetRangeUser(rMin,rMax);
 	hFitRes->GetYaxis()->SetRangeUser(yMin,yMax);

	TH1 *hFitResCore = 0;
	if( par_->respFuncType() == ResponseFunction::CrystalBall ) {
	  hFitResCore = static_cast<TH1D*>(hFitRes->Clone("PlotMCClosureFitResCore"));
	  hFitResCore->SetLineStyle(2);
	}
	for(int rBin = 1; rBin <= hFitRes->GetNbinsX(); ++rBin) {
	  double res = hFitRes->GetBinCenter(rBin);
	  hFitRes->SetBinContent(rBin,(*(par_->respFunc()))(res,pars));
	  if( hFitResCore ) hFitResCore->SetBinContent(rBin,par_->respFunc()->pdfGauss(res,1.,pars[0]));
	}
	if( hFitResCore ) hFitResCore->Draw("Lsame");
	hFitRes->Draw("Lsame");	
	
	gPad->RedrawAxis(); // Leads to slightly thicker font!!
	
	
	// Labels
	TPaveText *txt = util::LabelFactory::createPaveText(2);
	txt->AddText("PYTHIA, #sqrt{s} = 7 TeV, L = 50 pb^{-1}");
	txt->AddText("Anti-k_{T} d = 0.5 jets, "+par_->labelEtaBin());
	txt->Draw("same");

	TLegend *leg = util::LabelFactory::createLegend(2,1.,2);
	leg->AddEntry(hMCRes,"MC truth, "+par_->labelPtBin(bin,1),"P");
	leg->AddEntry(hFitRes,"Fit result, "+par_->labelPtBin(bin,0),"L");
	leg->Draw("same");
	
	if( logY ) gPad->SetLogy();
	
	// Write Canvas to fiel
	TString name = par_->outNamePrefix();
	name += "MCClosure_PtBin";
	name += bin;
	name += ".eps";
	can->SaveAs(name,"eps");
	
	// Clean up
	delete hMCRes;
	delete hFitRes;
	if( hFitResCore ) delete hFitResCore;
	delete txt;
	delete leg;
	delete can;
      }
    } else {
      std::cerr << "No MCClosure response distribution available." << std::endl;
    }
  }



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


  void FittedResolution::createSlides() const {
    // Open file
    TString name = par_->outNamePrefix();
    name += "SlidesMCClosure.tex";
    std::ofstream oFile(name);

    // Create slides with MC closure plots
    oFile << "% ----- MC closure plots ---------------------------" << std::endl;
    int nSlides = nPtBins()/6;
    if( nPtBins()%6 > 0 ) nSlides++;
    for(int slide = 0; slide < nSlides; ++slide) {
      oFile << "\n% --------------------------------------------------\n";
      oFile << "\\begin{frame}\n";
      oFile << "  \\frametitle{MC closure in various \\pt bins (" << slide << ")}\n";
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

    if( par_->respFuncType() == ResponseFunction::CrystalBall ) {
      // Create slides of all fit results
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
  }

    
  TGraphAsymmErrors *FittedResolution::getTGraphOfResolution(const TString &type) const {
    if( par_->verbosity() == 2 ) {
      std::cout << "FittedResolution: Creating graph of resolution" << std::endl;
    }
    std::vector<double> x;
    std::vector<double> ex;
    std::vector<double> y;
    std::vector<double> eyDown;
    std::vector<double> eyUp;

    if( type == "MCTruth" || type == "Statistic+MCTruth" ) {
      for(int i = 0; i < par_->nMCTruthPtBins(); ++i) {
	x.push_back(par_->mcTruthPtBinCenter(i));
	ex.push_back(0.);
	y.push_back( (trueRes_->Eval(x[i]))*(par_->mcTruthPseudoMeas(i)) );
	eyDown.push_back(y[i]*par_->mcTruthRelUncert());
	eyUp.push_back(y[i]*par_->mcTruthRelUncert());
      }
    }
    if( type != "MCTruth" ) {
      for(int i = 0; i < nPtBins(); ++i) {
	x.push_back(ptBins_[i]->meanPt());
	ex.push_back(ptBins_[i]->meanPtUncert());
	y.push_back(ptBins_[i]->extrapolatedValue(0));
	if( type == "Total" ) {
	  eyDown.push_back(ptBins_[i]->uncertDown(0));
	  eyUp.push_back(ptBins_[i]->uncertUp(0));
	} else if( type == "Statistic" || type == "Statistic+MCTruth" ) {
	  eyDown.push_back(ptBins_[i]->uncertStatDown(0));
	  eyUp.push_back(ptBins_[i]->uncertStatUp(0));
	} else if( type == "Systematic" ) {
	  eyDown.push_back(ptBins_[i]->uncertSystDown(0));
	  eyUp.push_back(ptBins_[i]->uncertSystUp(0));
	}
      }   
    }
    TGraphAsymmErrors *graph = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						     &(ex.front()),&(ex.front()),
						     &(eyDown.front()),&(eyUp.front()));
    graph->SetMarkerStyle(20);
    if( type == "Systematic" ) {
      graph->SetFillColor(43);
      graph->SetLineColor(43);
    } else if( type == "MCTruth" ) {
      graph->SetMarkerStyle(24);
      graph->SetMarkerColor(14);
      graph->SetLineColor(14);
    }
    return graph;
  }
}

