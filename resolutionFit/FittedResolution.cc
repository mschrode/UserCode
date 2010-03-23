#include "FittedResolution.h"

#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h" 
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


namespace resolutionFit {
  FittedResolution::FittedResolution(const std::vector<PtBin*> &ptBins, const std::vector<double> &trueResPar)
    : verbose_(true) {
    // Set ptbins
    ptBins_ = ptBins;

    // Create true resolution function
    assert( trueResPar.size() >= 3 );

    double xMin = 0.8*meanPt(0);
    double xMax = 1.1*meanPt(nPtBins()-1);

    trueRes_ = new TF1("FittedResolution::trueRes",
		       "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		       xMin,xMax);
    for(int i = 0; i < 3; i++) {
      trueRes_->SetParameter(i,trueResPar[i]);
    }
    trueRes_->SetLineColor(4);
    trueRes_->SetLineStyle(2);

    // Fit extrapolated resolution using
    // total uncertainty
    TGraphAsymmErrors *gStat = getTGraphOfResolution("Statistic");
    fittedRes_ = new TF1("FittedResolution::fittedRes_",
			 "sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
			 xMin,xMax);
    fittedRes_->SetParameter(0,0.);
    fittedRes_->SetParameter(0,1.);
    fittedRes_->SetParameter(0,0.03);
    fittedRes_->SetLineColor(2);
    gStat->Fit(fittedRes_,"0R");
    delete gStat;

    // Set style
    util::StyleSettings::presentationNoTitle();
  }


  FittedResolution::~FittedResolution() {
    for(std::vector<PtBin*>::iterator it = ptBins_.begin();
	it != ptBins_.end(); it++) {
      delete *it;
    }
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
      TH1D *h = (*it)->getFrameOfVariation(name);
      h->SetXTitle("p^{3}_{T,rel} - Schnitt");
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
      name = (*it)->minPtStr();
      name += " < p_{T} < ";
      name += (*it)->maxPtStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");
      
      // Write canvas to file
      name = "ExtrapolatedSigma_PtBin";
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
      TH1F *hResGen = (*it)->getHistResGen(name);
      hResGen->UseCurrentStyle();
      hResGen->GetXaxis()->SetRangeUser(0.4,1.6);
      hResGen->GetYaxis()->SetRangeUser(0.,1.8*hResGen->GetMaximum());
      hResGen->Draw();

      // Draw pdf
      name = "PlotPdfResolutionBin";
      name += bin;
      TH1F *hPdfRes = (*it)->getHistPdfRes(name);
      hPdfRes->Draw("Lsame");
      hResGen->Draw("same");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1,0.9);
      name = (*it)->minPtStr();
      name += " < p_{T} < ";
      name += (*it)->maxPtStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(2,0.9,util::LabelFactory::lineHeight(),util::LabelFactory::lineHeight());
      leg->AddEntry(hResGen,"MC truth","L");
      char entry[100];
      sprintf(entry,"Fit #sigma/p_{T}, p_{T} = %.1f GeV",(*it)->meanPt());
      leg->AddEntry(hPdfRes,entry,"L");
      leg->Draw("same");

      // Write Canvas to fiel
      name = "Resolution_PtBin";
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


  void FittedResolution::plotResolution() const {
    // Create Canvas
    TCanvas *can = new TCanvas("CanExtrapolatedResolution","Extrapolated Resolution",500,500);
    can->cd();

    // Draw a frame
    double xMin = 0.8*meanPt(0);
    double xMax = 1.1*meanPt(nPtBins()-1);
    TH1D *h = new TH1D("FrameExtrapolatedResolution",";p_{T} (GeV);#sigma / p_{T}",
		       1000,xMin,xMax);
    h->GetYaxis()->SetRangeUser(0.7*relSigma(nPtBins()-1),0.185);
    h->Draw();

    // Draw systematic uncertainty band
    TGraphAsymmErrors *gSyst = getTGraphOfResolution("Systematic");
    gSyst->Draw("E3same");

    // Draw true and fitted resolution functions
    trueRes_->Draw("same");
    fittedRes_->Draw("same");

    // Draw graph with statistical uncertainties
    TGraphAsymmErrors *gStat = getTGraphOfResolution("Statistic");
    gStat->Draw("PE1same");

    // Add a legend
    TLegend *leg = util::LabelFactory::createLegend(5);
    leg->AddEntry(gStat,"Extrapolated #bar{#sigma} (p^{3}_{T,rel} #rightarrow 0)","P");
    leg->AddEntry(gStat,"Statistical uncertainty","L");
    leg->AddEntry(gSyst,"Uncertainty from spectrum","F");
    leg->AddEntry(trueRes_,"MC truth","L");
    leg->AddEntry(fittedRes_,"Fit #sigma(p_{T}) to #bar{#sigma}","L");
    leg->Draw("same");

    // Write canvas to file
    can->SaveAs("ExtrapolatedResolution.eps","eps");

    // Clean up
    delete h;
    delete gSyst;
    delete gStat;
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

      // Draw MC truth spectrum
      name = "PlotPtGenPtBin";
      name += bin;
      TH1F *hPtGen = (*it)->getHistPtGen(name);
      hPtGen->UseCurrentStyle();
      hPtGen->GetYaxis()->SetRangeUser(0.,2.4*hPtGen->GetMaximum());
      hPtGen->Draw();

      // Draw pdf
      name = "PlotPdfPtTrueBin";
      name += bin;
      TH1F *hPdfPtTrue = (*it)->getHistPdfPtTrue(name);
      hPdfPtTrue->Draw("Lsame");
      hPtGen->Draw("same");

      // Labels
      TPaveText *txt = util::LabelFactory::createPaveText(1,1.,0.06);
      name = (*it)->minPtStr();
      name += " < p_{T} < ";
      name += (*it)->maxPtStr();
      name += " GeV";
      txt->AddText(name);
      txt->Draw("same");

      TLegend *leg = util::LabelFactory::createLegend(3,1.,0.09,0.06);
      leg->AddEntry(hPtGen,"MC truth","L");
      leg->AddEntry(hPdfPtTrue,"Spectrum #tilde{f}(p_{T})","L");
      char entry[100];
      sprintf(entry,"#frac{#sigma_{MC}(p_{T})}{p_{T}} = #frac{%.3f}{#sqrt{p_{T} / GeV}} #oplus %.3f",
	      trueRes_->GetParameter(1),trueRes_->GetParameter(2));
      TH1F *hDummy = new TH1F("hDummy","",10,0,1);
      hDummy->SetLineColor(0);
      leg->AddEntry(hDummy,entry,"L");
      leg->Draw("same");

      // Write Canvas to fiel
      name = "Spectrum_PtBin";
      name += bin;
      name += ".eps";
      can->SaveAs(name,"eps");

      // Clean up
      delete hPtGen;
      delete hPdfPtTrue;
      delete hDummy;
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
    TH1D *h = new TH1D("FrameSystematicUncertainties",";p_{T} (GeV);#Delta#sigma / #sigma",
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
      gSystDown[n]->Draw("Lsame");
    }

    // Add legend
    leg->Draw("same");

    // Write canvas to file
    can->SaveAs("SystematicUncertainties.eps","eps");

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
      std::cout << " (" << minPt(bin) << " - " << maxPt(bin) << "): " << std::flush;
      std::cout << relSigma(bin) << " (" << uncertStat(bin) << std::flush;
      std::cout << ", +" << uncertSystUp(bin) << ", -" << uncertSystDown(bin) << std::flush;
      std::cout << "),  " << fittedRes_->Eval(meanPt(bin)) << std::flush;
      std::cout << ",  " << trueRes_->Eval(meanPt(bin)) << std::endl;
    } 
  }

    
  TGraphAsymmErrors *FittedResolution::getTGraphOfResolution(const TString &uncertType) const {
    if( verbose_ ) {
      std::cout << "FittedResolution: Creating graph of resolution" << std::endl;
    }
    std::vector<double> x(nPtBins());
    std::vector<double> ex(nPtBins(),0.);
    std::vector<double> y(nPtBins());
    std::vector<double> eyDown(nPtBins());
    std::vector<double> eyUp(nPtBins());
    for(int i = 0; i < nPtBins(); i++) {
      x[i] = ptBins_[i]->meanPt();
      y[i] = ptBins_[i]->relSigma();
      if( uncertType == "Total" ) {
	eyDown[i] = ptBins_[i]->uncertDown();
	eyUp[i] = ptBins_[i]->uncertUp();
      } else if( uncertType == "Statistic" ) {
	eyDown[i] = ptBins_[i]->uncertStatDown();
	eyUp[i] = ptBins_[i]->uncertStatUp();
      } else if( uncertType == "Systematic" ) {
	eyDown[i] = ptBins_[i]->uncertSystDown();
	eyUp[i] = ptBins_[i]->uncertSystUp();
      } else {
	eyDown[i] = 0.;
	eyUp[i] = 0.;
	std::cerr << "WARNING (FittedResolution): Unknown uncertainty type '" << uncertType << "'" << std::endl;
      }
    }   
    TGraphAsymmErrors *graph = new TGraphAsymmErrors(nPtBins(),&(x.front()),&(y.front()),
						     &(ex.front()),&(ex.front()),
						     &(eyDown.front()),&(eyUp.front()));
    graph->SetMarkerStyle(20);
    if( uncertType == "Systematic" ) {
      graph->SetFillColor(43);
      graph->SetLineColor(43);
    }
    return graph;
  }
}

