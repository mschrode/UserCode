#include "FittedResolution.h"

#include <algorithm>
#include <cmath>

#include "TCanvas.h"
#include "TError.h"
#include "TH1.h"
#include "TF1.h" 
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


namespace resolutionFit {
  FittedResolution::FittedResolution(const std::vector<PtBin*> &ptBins, const Parameters *par) 
    : par_(par), ptBins_(ptBins) {

    // Do not print ROOT message if eps file has been created
    gErrorIgnoreLevel = 1001;

    if( par_->hasMCTruthBins() ) std::cout << "Adding pseudo gamma+jet measurement from MC truth" << std::endl;

    xMin_ = 0.8*meanPt(0);
    if( par_->hasMCTruthBins() ) xMin_ = 0.8*par_->mcTruthPtMin(0);
    xMax_ = 1.1*meanPt(nPtBins()-1);

    trueRes_ = new TF1("FittedResolution::trueRes",
		       "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		       xMin_,xMax_);
    for(int i = 0; i < 3; i++) {
      trueRes_->SetParameter(i,par_->trueGaussResPar(i));
    }
    trueRes_->SetLineColor(4);
    trueRes_->SetLineStyle(2);

    // Fit extrapolated resolution using
    // statistic uncertainty
    TGraphAsymmErrors *gStat = getTGraphOfResolution("Statistic+MCTruth");
    fittedRes_ = new TF1("FittedResolution::fittedRes_",
			 "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
			 xMin_,xMax_);
    for(int i = 0; i < 3; ++i) {
      fittedRes_->SetParameter(i,par_->trueGaussResPar(i));
    }
    //    fittedRes_->FixParameter(0,par_->trueGaussResPar(0));
    fittedRes_->SetLineColor(2);
    if( par_->fitExtrapolatedSigma() ) gStat->Fit(fittedRes_,"0R");
    delete gStat;

    // Set style
    util::StyleSettings::presentationNoTitle();
  }


  FittedResolution::~FittedResolution() {
    delete trueRes_;
    delete fittedRes_;
  }


  void FittedResolution::plotExtrapolation() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int bin = (it-ptBins_.begin());

      // Create Canvas
      TString name = "PlotExtrapolationPtBin";
      name += bin;
      TCanvas *can = new TCanvas(name,name,500,500);
      can->cd();

      // Draw a frame
      name = "FrameExtrapolationPtBin";
      name += bin;
      TH1 *h = (*it)->getFrameOfVariation(name);
      h->Draw();

      // Draw graph and extrapolation
      TGraphAsymmErrors *g = (*it)->getTGraphOfVariation();
      g->Draw("PE1same");
      name = "FitExtrapolationPtBin";
      name += bin;
      TF1 *f = (*it)->getTF1OfVariation(name);
      f->Draw("same");
      g->Draw("PE1same");

      // Draw label
      TPaveText *txt = util::LabelFactory::createPaveText(1,1.,0.06);
      name = (*it)->ptMinStr();
      name += " < p^{recoJet}_{T} < ";
      name += (*it)->ptMaxStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");
      
      // Write canvas to file
      name = par_->outNamePrefix();
      name += "ExtrapolatedSigma_PtBin";
      name += bin;
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


  void FittedResolution::plotResolutionBins() const {
    // Loop over ptbins
    for(std::vector<PtBin*>::const_iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      int bin = (it-ptBins_.begin());

      // Create Canvas
      TString name = "PlotResolutionPtBin";
      name += bin;
      TCanvas *can = new TCanvas(name,name,500,500);
      can->cd();

      // Draw MC truth resolution
      name = "PlotResolutionBin";
      name += bin;
      TH1 *hResGen = (*it)->getHistResGen(name);
      hResGen->UseCurrentStyle();
      hResGen->SetMarkerStyle(20);
      hResGen->GetXaxis()->SetRangeUser(0.4,1.6);
      hResGen->GetYaxis()->SetRangeUser(0.,1.8*hResGen->GetMaximum());
      hResGen->Draw("PE1");

      // Draw pdf
      name = "PlotPdfResolutionBin";
      name += bin;
      TH1 *hPdfRes = (*it)->getHistPdfRes(name);
      hPdfRes->Draw("Lsame");
      hResGen->Draw("same");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1,0.9);
      name = (*it)->ptMinStr();
      name += " < p^{recoJet}_{T} < ";
      name += (*it)->ptMaxStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,1.,util::LabelFactory::lineHeight(),util::LabelFactory::lineHeight());
      leg->AddEntry(hResGen,"MC truth","P");
      char entry[100];
      sprintf(entry,"Fitted #sigma/p^{true}_{T}, p^{true}_{T} = %.1f GeV",(*it)->meanPt());
      leg->AddEntry(hPdfRes,entry,"L");
      leg->Draw("same");

      gPad->RedrawAxis();

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
      TString name = "PlotPtAsymmetryPtBin";
      name += bin;
      TCanvas *can = new TCanvas(name,name,500,500);
      can->cd();

      // Draw MC truth resolution
      name = "PlotPtAsymmetryBin";
      name += bin;
      TH1 *hPtAsym = (*it)->getHistPtAsym(name);
      hPtAsym->UseCurrentStyle();
      hPtAsym->SetMarkerStyle(20);
      hPtAsym->GetXaxis()->SetRangeUser(0.4,1.6);
      hPtAsym->GetYaxis()->SetRangeUser(0.,1.8*hPtAsym->GetMaximum());
      hPtAsym->GetYaxis()->SetTitle("1 / N  dN / dA");
      hPtAsym->Draw("PE1");

      // Draw pdf
      name = "PlotPdfPtAsymmetryBin";
      name += bin;
      TH1 *hPdfPtAsym = (*it)->getHistPdfPtAsym(name);
      hPdfPtAsym->Draw("Lsame");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1,0.9);
      name = (*it)->ptMinStr();
      name += " < p^{recoJet}_{T} < ";
      name += (*it)->ptMaxStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,1.,util::LabelFactory::lineHeight(),util::LabelFactory::lineHeight());
      leg->AddEntry(hPtAsym,"MC truth","P");
      char entry[100];
      sprintf(entry,"Fitted  #sqrt{2}#sigma/p^{true}_{T}, p^{true}_{T} = %.1f GeV",(*it)->meanPt());
      leg->AddEntry(hPdfPtAsym,entry,"L");
      leg->Draw("same");

      gPad->RedrawAxis();

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
    TH1 *h = new TH1D("FrameExtrapolatedResolution",";p_{T} (GeV);#sigma / p_{T}",
		      1000,xMin_,xMax_);
    double min = *std::min_element(gStat->GetY(),gStat->GetY()+gStat->GetN());
    double max = *std::max_element(gStat->GetY(),gStat->GetY()+gStat->GetN());
    if( gPseudo ) max = *std::max_element(gPseudo->GetY(),gPseudo->GetY()+gPseudo->GetN());
    h->GetYaxis()->SetRangeUser(0.7*min,1.2*max);
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
    int nLegEntries = 2;
    if( par_->fitExtrapolatedSigma() ) nLegEntries++;
    if( gPseudo ) nLegEntries++;
    TLegend *leg = util::LabelFactory::createLegend(nLegEntries);
    leg->AddEntry(gStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    //leg->AddEntry(gSyst,"Uncertainty from spectrum","F");
    if( gPseudo ) leg->AddEntry(gPseudo,"Pseudo #gamma+jet measurement","P");
    leg->AddEntry(trueRes_,"MC truth resolution","L");
    if( par_->fitExtrapolatedSigma() ) leg->AddEntry(fittedRes_,"Fit #sigma(p_{T}) to #bar{#sigma}","L");
    leg->Draw("same");

    // Write canvas to file
    gPad->SetLogx();
    gPad->RedrawAxis();
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

    TF1 *lineFitRatio = new TF1("LineFitRatioExtraplolatedResolution","pol0",xMin_,xMax_);
    lineFitRatio->SetLineWidth(2);
    lineFitRatio->SetLineStyle(2);
    lineFitRatio->SetLineColor(2);
    if( par_->fitRatio() ) gStatRatio->Fit(lineFitRatio,"0QR");

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
    h->GetYaxis()->SetRangeUser(0.65,1.75);
    
    max = *(std::max_element(gStatRatio->GetY(),gStatRatio->GetY()+gStat->GetN()));
    if( max > 1.45 ) h->GetYaxis()->SetRangeUser(0.65,1.25*max);
    h->Draw();
    if( par_->fitExtrapolatedSigma() ) hFittedRatio->Draw("Lsame");
    if( par_->fitRatio() ) lineFitRatio->Draw("same");
    if( par_->hasStartOffset() ) lineStartRes->Draw("same");
    gStatRatio->Draw("PE1same");
    if( gPseudoRatio ) gPseudoRatio->Draw("PE1same");

    // Add a legend
    delete leg;
    nLegEntries = 2;
    if( par_->fitExtrapolatedSigma() ) nLegEntries++;
    if( gPseudo ) nLegEntries++;
    if( par_->fitRatio() ) nLegEntries++;
    if( par_->hasStartOffset() ) nLegEntries++;
    leg = util::LabelFactory::createLegend(nLegEntries);
    leg->AddEntry(gStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    //leg->AddEntry(gSyst,"Uncertainty from spectrum","F");
    if( gPseudo ) leg->AddEntry(gPseudo,"Pseudo #gamma+jet measurement","P");
    leg->AddEntry(trueRes_,"MC truth resolution","L");
    if( par_->fitRatio() ) leg->AddEntry(lineFitRatio,"Mean deviation","L");
    if( par_->hasStartOffset() ) leg->AddEntry(lineStartRes,"Resolution in spectrum","L");
    if( par_->fitExtrapolatedSigma() ) leg->AddEntry(fittedRes_,"Fit #sigma(p_{T}) to #bar{#sigma}","L");
    leg->Draw("same");

    // Write canvas to file
    gPad->SetLogx();
    gPad->RedrawAxis();
    can->SaveAs(par_->outNamePrefix()+"ExtrapolatedResolutionRatio.eps","eps");


    // Clean up
    delete h;
    //    delete gSyst;
    delete gStat;
    if( gPseudo ) delete gPseudo;
    if( gPseudoRatio ) delete gPseudoRatio;
    delete gStatRatio;
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

      // Create Canvas
      TString name = "PlotSpectrumPtBin";
      name += bin;
      TCanvas *can = new TCanvas(name,name,500,500);
      can->cd();

      // Draw MC truth spectra
      name = "PlotPtGenPtBin";
      name += bin;
      TH1 *hPtGen = (*it)->getHistPtGen(name);
      hPtGen->UseCurrentStyle();
      hPtGen->SetLineWidth(1);
      hPtGen->SetXTitle("p^{gen}_{T}  (GeV)");
      hPtGen->SetYTitle("1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
      hPtGen->GetYaxis()->SetRangeUser(0.,2.4*hPtGen->GetMaximum());
      hPtGen->SetMarkerStyle(24);
      //hPtGen->Draw("PE1");

      name = "PlotPtGenJet1PtBin";
      name += bin;
      TH1 *hPtGenJet1 = (*it)->getHistPtGenJet1(name);
      hPtGenJet1->UseCurrentStyle();
      hPtGenJet1->SetLineWidth(1);
      hPtGenJet1->SetXTitle("p^{gen}_{T,1}  (GeV)");
      hPtGenJet1->SetYTitle("1 / N  dN / dp^{gen}_{T,1}  1 / (GeV)");
      hPtGenJet1->GetYaxis()->SetRangeUser(0.,2.*hPtGenJet1->GetMaximum());
      hPtGenJet1->SetMarkerStyle(20);
      hPtGenJet1->Draw("PE1");

      // Draw pdf
      name = "PlotPdfPtTrueBin";
      name += bin;
      TH1 *hPdfPtTrue = (*it)->getHistPdfPtTrue(name);
      hPdfPtTrue->Draw("Lsame");
      gPad->RedrawAxis();

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1,-0.7,0.06);
      name = (*it)->ptMinStr();
      name += " < p^{recoJet}_{T} < ";
      name += (*it)->ptMaxStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,-0.6,0.07,0.07);
      leg->AddEntry(hPtGenJet1,"MC truth: p^{gen}_{T,1}","P");
      //leg->AddEntry(hPtGen,"MC truth: p^{particleJet}_{T,1+2}","P");
      leg->AddEntry(hPdfPtTrue,"Spectrum  f(p_{T})","L");
      leg->Draw("same");

      // Write Canvas to fiel
      name = par_->outNamePrefix();
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
    // Create graphs of systematic uncertainties
    std::vector<double> pt(nPtBins());
    for(int bin = 0; bin < nPtBins(); bin++) {
      pt[bin] = meanPt(bin);
    }
    std::vector<TGraph*> gSystUp(nUncertSyst());
    std::vector<TGraph*> gSystDown(nUncertSyst());
    TLegend *leg = util::LabelFactory::createLegend(gSystUp.size());
    // Loop over systematic uncertainties
    for(int n = 0; n < nUncertSyst(); n++) {
      std::vector<double> systUp(nPtBins());
      std::vector<double> systDown(nPtBins());
      for(int bin = 0; bin < nPtBins(); bin++) {
	const Uncertainty *uncertSyst = ptBins_[bin]->uncertSyst();
	systUp[bin] = uncertSyst->up(n)/relSigma(bin);
	systDown[bin] = uncertSyst->down(n)/relSigma(bin);
      }
      gSystUp[n] = new TGraph(pt.size(),&(pt.front()),&(systUp.front()));
      gSystDown[n] = new TGraph(pt.size(),&(pt.front()),&(systDown.front()));
      gSystUp[n]->SetLineColor(util::StyleSettings::color(n));
      gSystDown[n]->SetLineColor(util::StyleSettings::color(n));
      gSystUp[n]->SetLineWidth(2);
      gSystDown[n]->SetLineWidth(2);

      leg->AddEntry(gSystUp[n],ptBins_[0]->uncertSyst()->label(n),"L");
    }

    // Create Canvas
    TCanvas *can = new TCanvas("CanSystematicUncertainties","Systematic uncertainties",500,500);
    can->cd();

    // Create frame histogram
    double xMin = 0.8*meanPt(0);
    double xMax = 1.1*meanPt(nPtBins()-1);
    TH1 *h = new TH1D("FrameSystematicUncertainties",";p_{T} (GeV);#Delta#sigma / #sigma",
		       1000,xMin,xMax);
    for(int bin = 1; bin < h->GetNbinsX(); bin++) {
      h->SetBinContent(bin,0.);
    }
    h->GetYaxis()->SetRangeUser(-0.1,0.2);
    h->SetLineStyle(2);
    h->Draw();

    // Plot systematic uncertainties
    for(int n = 0; n < gSystUp.size(); n++) {
      gSystUp[n]->Draw("Lsame");
      //gSystDown[n]->Draw("Lsame");
    }

    // Add legend
    leg->Draw("same");

    // Write canvas to file
    can->SaveAs(par_->outNamePrefix()+"SystematicUncertainties.eps","eps");

    // Clean up
    delete h;
    for(int n = 0; n < gSystUp.size(); n++) {
      delete gSystUp[n];
      delete gSystDown[n];
    }
    delete leg;
    delete can;
  }


  void FittedResolution::print() const {
    std::cout << "\n\nFITTED RESOLUTION\n";
    // Loop over ptbins
    for(size_t bin = 0; bin < nPtBins(); bin++) {
      std::cout << bin << ": " << meanPt(bin) << std::flush;
      std::cout << " (" << ptMin(bin) << " - " << ptMax(bin) << "): " << std::flush;
      std::cout << relSigma(bin)*meanPt(bin) << ", " << std::flush;
      std::cout << relSigma(bin) << " (" << uncertStat(bin) << std::flush;
      std::cout << ", +" << uncertSystUp(bin) << ", -" << uncertSystDown(bin) << std::flush;
      std::cout << "),  " << fittedRes_->Eval(meanPt(bin)) << std::flush;
      std::cout << ",  " << trueRes_->Eval(meanPt(bin)) << std::endl;
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
	y.push_back(ptBins_[i]->relSigma());
	if( type == "Total" ) {
	  eyDown.push_back(ptBins_[i]->uncertDown());
	  eyUp.push_back(ptBins_[i]->uncertUp());
	} else if( type == "Statistic" || type == "Statistic+MCTruth" ) {
	  eyDown.push_back(ptBins_[i]->uncertStatDown());
	  eyUp.push_back(ptBins_[i]->uncertStatUp());
	} else if( type == "Systematic" ) {
	  eyDown.push_back(ptBins_[i]->uncertSystDown());
	  eyUp.push_back(ptBins_[i]->uncertSystUp());
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

