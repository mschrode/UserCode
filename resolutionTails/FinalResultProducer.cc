// $Id: FinalResultProducer.cc,v 1.2 2013/05/10 13:22:19 mschrode Exp $

#ifndef RESOLUTION_TAILS_FINAL_RESULT_PRODUCER
#define RESOLUTION_TAILS_FINAL_RESULT_PRODUCER

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "FitParameters.h"
#include "Output.h"
#include "Style.h"
#include "Uncertainty.h"

#include "../sampleTools/BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"


// Produce final plots
// Merge scale factors from nominal and systematic variations into final plots
// Compute total systematic uncertainties

namespace resolutionTails {
  class FinalResultProducer {
  public:
    FinalResultProducer(const TString &fileName, const resolutionTails::FitParameters* fitPars, const sampleTools::BinningAdmin* binAdm, const std::vector<resolutionTails::Uncertainty::SystematicVariation> &variations, const Style* style);
    ~FinalResultProducer();

    void makeScaleFactorPlots();
    void makeUncertaintyPlots();
    void print() const;


  private:
    const sampleTools::BinningAdmin* binAdm_;
    const FitParameters* fitPars_;
    const Style* style_;
    Output* out_;

    std::vector<TH1*> hFrames_;
    std::vector<TGraphAsymmErrors*> gScaleFactors_;
    std::vector<TGraphAsymmErrors*> gAbsUncertTot_;
    std::vector< std::vector<TGraphAsymmErrors*> > gStackedRelUncerts_;
    std::vector<TString> uncertLabels_;

    void init(const TString &inFilePrefix, const std::vector<resolutionTails::Uncertainty::SystematicVariation> &variations);
    TGraphAsymmErrors* correctBinCenter(const TH1* h1, const TH1* h2) const;
    TGraphAsymmErrors* relativeDifference(const TH1* hNom, int color, double weight, const TH1* hVarUp, const TH1* hVarDn = 0) const;
    TGraphAsymmErrors* sumErrorsInQuadrature(const std::vector<TGraphAsymmErrors*> &graphs) const;
    TGraphAsymmErrors* uncertaintyBand(const TH1* hNom, const TGraphAsymmErrors* gRelUncert, int color) const;
    std::vector<TGraphAsymmErrors*> uncertaintyStack(const std::vector<TGraphAsymmErrors*> &gUncerts, unsigned int etaBin) const;
  };



  // fileName: output from ScaleFactor producer for *nominal* variation
  // ------------------------------------------------------------------------------------
  FinalResultProducer::FinalResultProducer(const TString &fileName, const FitParameters* fitPars, const sampleTools::BinningAdmin* binAdm, const std::vector<resolutionTails::Uncertainty::SystematicVariation> &variations, const Style* style) 
    : binAdm_(binAdm), fitPars_(fitPars), style_(style) {
    TString inFilePrefix = util::baseName(fileName);
    out_ = new Output(inFilePrefix+"_Final",true,true,false);
    init(inFilePrefix,variations);
  }


  // ------------------------------------------------------------------------------------
  void FinalResultProducer::init(const TString &inFilePrefix, const std::vector<resolutionTails::Uncertainty::SystematicVariation> &variations) {
    if( Output::DEBUG ) std::cout << "FinalResultProducer::init() Entering" << std::endl;
    // Nominal results
    TString inFileNameNom = inFilePrefix+".root";

    // Loop over eta bins
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      if( Output::DEBUG ) std::cout << "  EtaBin " << etaBin << std::endl;
    
      //// 1. Frame for scale factor and uncertainty plots
      //      (has the correct binning)
      if( Output::DEBUG ) std::cout << "    Frame" << std::endl;
      TString histName = "hScaleFrame_Eta"+util::toTString(etaBin);
      TH1* hScaleFrame = util::FileOps::readTH1(inFileNameNom,histName,histName);
      hScaleFrame->GetXaxis()->SetTitle("p^{ave}_{T} (GeV)");
      hScaleFrame->SetLineWidth(style_->lineWidth());
      hScaleFrame->SetLineStyle(2);
      hFrames_.push_back(hScaleFrame);


      //// 2. Nominal scale factors
      if( Output::DEBUG ) std::cout << "    Scale factors" << std::endl;
      // Read mean ptAve per ptAve bin in this eta bin
      histName = "hPtAveMeanData_Eta"+util::toTString(etaBin);
      TH1* hPtAveMean = util::FileOps::readTH1(inFileNameNom,histName,histName);

      // Read nominal scaling factors per ptAve bin in this eta bin
      histName = "hScale_Eta"+util::toTString(etaBin);
      TH1* hScaleNom = util::FileOps::readTH1(inFilePrefix+".root",histName,histName+"_Nominal");

      // Graph of nominal scaling factors with ptAveMean as x
      gScaleFactors_.push_back(correctBinCenter(hScaleNom,hPtAveMean));


      //// 3. Systematic uncertainties
      if( Output::DEBUG ) std::cout << "    Systematic uncertainties" << std::endl;
      TH1* hVarCoreUp = 0;
      TH1* hVarCoreDn = 0;
      TH1* hVarClosure = 0;
      TH1* hVarPUUp = 0;
      TH1* hVarPUDn = 0;
      TH1* hVarExtra = 0;

      // Loop over systematic variations
      std::vector<Uncertainty::SystematicVariation>::const_iterator vit = variations.begin();
      for(;vit != variations.end(); ++vit) {
	if( *vit != Uncertainty::Nominal ) {
	  TH1* h = util::FileOps::readTH1(inFilePrefix+"_"+Uncertainty::id(*vit)+".root",histName,histName+"_"+Uncertainty::id(*vit));
	
	  if( *vit == Uncertainty::CoreUp ) hVarCoreUp = h;
	  else if( *vit == Uncertainty::CoreDn ) hVarCoreDn = h;
	  else if( *vit == Uncertainty::Closure ) hVarClosure = h;
	  else if( *vit == Uncertainty::PUUp ) hVarPUUp = h;
	  else if( *vit == Uncertainty::PUDn ) hVarPUDn = h;
	  else if( *vit == Uncertainty::Extrapolation ) hVarExtra = h;
	}
      }

      // Get relative uncertainties
      std::vector<TGraphAsymmErrors*> gRelUncerts;
      if( hVarPUUp && hVarPUDn ) {
	gRelUncerts.push_back(relativeDifference(hScaleNom,11,1.,hVarPUUp,hVarPUDn));
	if( uncertLabels_.size() < gRelUncerts.size() ) {
	  uncertLabels_.push_back("PU");
	}
	delete hVarPUUp;
	delete hVarPUDn;
      }
      if( hVarExtra ) {
	gRelUncerts.push_back(relativeDifference(hScaleNom,38,1.,hVarExtra));
	if( uncertLabels_.size() < gRelUncerts.size() ) {
	  uncertLabels_.push_back("Additional Jets");
	}
	delete hVarExtra;
      }
      if( hVarClosure ) {
	gRelUncerts.push_back(relativeDifference(hScaleNom,8,0.5,hVarClosure));
	if( uncertLabels_.size() < gRelUncerts.size() ) {
	  uncertLabels_.push_back("Bias");
	}
	delete hVarClosure;
      }
      if( hVarCoreUp && hVarCoreDn ) {
	gRelUncerts.push_back(relativeDifference(hScaleNom,46,1.,hVarCoreUp,hVarCoreDn));
	if( uncertLabels_.size() < gRelUncerts.size() ) {
	  uncertLabels_.push_back("#sigma_{A} Correction");
	}
	delete hVarCoreUp;
	delete hVarCoreDn;
      }

      // Compute total error band (sum errors in quadrature)
      TGraphAsymmErrors* gRelUncertTot = sumErrorsInQuadrature(gRelUncerts);
      gAbsUncertTot_.push_back(uncertaintyBand(hScaleNom,gRelUncertTot,5));
      delete gRelUncertTot;
      delete hScaleNom;
      delete hPtAveMean;


      // Compute stacked (quadratic addition) uncertainties for plot of
      // relative systematic uncertainties
      gStackedRelUncerts_.push_back(uncertaintyStack(gRelUncerts,etaBin));
      for(unsigned int i = 0; i < gRelUncerts.size(); ++i) {
	delete gRelUncerts.at(i);
      }
    } // End of loop over eta bins
    if( Output::DEBUG ) std::cout << "FinalResultProducer::init() Leaving" << std::endl;
  }


  // ------------------------------------------------------------------------------------
  FinalResultProducer::~FinalResultProducer() {
    for(std::vector<TH1*>::iterator it = hFrames_.begin();
	it != hFrames_.end(); ++it) {
      delete *it;
    }
    for(std::vector<TGraphAsymmErrors*>::iterator it = gScaleFactors_.begin();
	it != gScaleFactors_.end(); ++it) {
      delete *it;
    }
    for(std::vector<TGraphAsymmErrors*>::iterator it = gAbsUncertTot_.begin();
	it != gAbsUncertTot_.end(); ++it) {
      delete *it;
    }
    for(std::vector< std::vector<TGraphAsymmErrors*> >::iterator it = gStackedRelUncerts_.begin();
	it != gStackedRelUncerts_.end(); ++it) {
      for(std::vector<TGraphAsymmErrors*>::iterator it2 = it->begin();
	  it2 != it->end(); ++it2) {
	delete *it2;
      }
    }
    delete out_;
  }


  // ------------------------------------------------------------------------------------
  void FinalResultProducer::makeScaleFactorPlots() {
    if( Output::DEBUG ) std::cout << "FinalResultProducer::makeScaleFactorPlots() Entering" << std::endl;

    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      if( Output::DEBUG ) std::cout << "  EtaBin " << etaBin << std::endl;

      // Labels
      TPaveText* label = util::LabelFactory::createPaveText(1);
      label->AddText(style_->labelLumi()+",  "+style_->labelWindow(fitPars_->nSigTailStart(),fitPars_->nSigTailEnd())+",  "+util::LabelFactory::etaCut(binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin)));
      TLegend* leg = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
      leg->AddEntry(gScaleFactors_.at(etaBin),"Extrapolation Uncertainty #deltaf_{ex}","L");
      leg->AddEntry(gAbsUncertTot_.at(etaBin),"Systematic Uncertainty","F");

      // Prepare frame
      hFrames_.at(etaBin)->SetTitle(style_->title());
      hFrames_.at(etaBin)->GetYaxis()->SetRangeUser(0.01,3.3);
      hFrames_.at(etaBin)->GetYaxis()->SetTitle(style_->labelScaleFactor());
      hFrames_.at(etaBin)->GetXaxis()->SetMoreLogLabels();
      hFrames_.at(etaBin)->GetXaxis()->SetNoExponent();
      for(int bin = 1; bin <= hFrames_.at(etaBin)->GetNbinsX(); ++bin) {
	hFrames_.at(etaBin)->SetBinContent(bin,1.);
      }


      // Plot scale factors with total uncertainty
      TCanvas* can = util::HistOps::createTCanvas("can","can",500,500);
      can->cd();
      hFrames_.at(etaBin)->Draw("HIST");
      gAbsUncertTot_.at(etaBin)->Draw("E2same");
      hFrames_.at(etaBin)->Draw("HISTsame");
      gScaleFactors_.at(etaBin)->Draw("PE1same");
      label->Draw("same");
      leg->Draw("same");
      can->SetLogx();
      gPad->RedrawAxis();
      out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_ScaleFactors");
      delete label;
      delete leg;
      delete can;
    } // End of loop over eta bins
    if( Output::DEBUG ) std::cout << "FinalResultProducer::makeScaleFactorPlots() Leaving" << std::endl;
  }


  // ------------------------------------------------------------------------------------
  void FinalResultProducer::makeUncertaintyPlots() {
    if( Output::DEBUG ) std::cout << "FinalResultProducer::makeUncertaintyPlots() Entering" << std::endl;

    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      if( Output::DEBUG ) std::cout << "  EtaBin " << etaBin << std::endl;

      // Labels
      TPaveText* label = util::LabelFactory::createPaveText(1);
      label->AddText(style_->labelWindow(fitPars_->nSigTailStart(),fitPars_->nSigTailEnd())+",   "+util::LabelFactory::etaCut(binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin)));
      TLegend* legUncert1 = util::LabelFactory::createLegendColWithOffset(2,-0.5,label->GetSize());
      TLegend* legUncert2 = util::LabelFactory::createLegendColWithOffset(2,0.5,label->GetSize());
      for(int i = static_cast<int>(gStackedRelUncerts_.at(etaBin).size())-1; i >= 0; --i) {
	if( i > 1 ) legUncert1->AddEntry(gStackedRelUncerts_.at(etaBin).at(i),uncertLabels_.at(i),"F");
	else        legUncert2->AddEntry(gStackedRelUncerts_.at(etaBin).at(i),uncertLabels_.at(i),"F");
      }

      // Frame
      hFrames_.at(etaBin)->GetYaxis()->SetRangeUser(0.,65);
      hFrames_.at(etaBin)->GetXaxis()->SetMoreLogLabels();
      hFrames_.at(etaBin)->GetXaxis()->SetNoExponent();
      hFrames_.at(etaBin)->GetYaxis()->SetTitle("Relative Uncertainty on "+style_->labelScaleFactor()+"  (%)");
      for(int bin = 1; bin <= hFrames_.at(etaBin)->GetNbinsX(); ++bin) {
	hFrames_.at(etaBin)->SetBinContent(bin,-1.);
      }

      // Plot relative (stacked) uncertainties
      TCanvas* can = util::HistOps::createTCanvas("can","can",500,500);
      can->cd();
      hFrames_.at(etaBin)->Draw("HIST");
      for(std::vector<TGraphAsymmErrors*>::reverse_iterator it = gStackedRelUncerts_.at(etaBin).rbegin();
	  it != gStackedRelUncerts_.at(etaBin).rend(); ++it) {
	(*it)->Draw("E3same");
      }
      label->Draw("same");
      legUncert1->Draw("same");
      legUncert2->Draw("same");
      can->SetLogx();
      gPad->RedrawAxis();
      out_->toFiles(can,"EtaBin"+util::toTString(etaBin)+"_Uncertainties");
      delete label;
      delete legUncert1;
      delete legUncert2;
      delete can;
    }
    if( Output::DEBUG ) std::cout << "FinalResultProducer::makeUncertaintyPlots() Leaving" << std::endl;
  }



  // Print scale factors with total uncertainty (stat + syst)
  // to screen and in LaTeX table format to file
  // ------------------------------------------------------------------------------------
  void FinalResultProducer::print() const {
    if( Output::DEBUG ) std::cout << "FinalResultProducer::print() Entering" << std::endl;

    char tmp[100];
    TString text = "";

    // Table body
    for(unsigned int etaBin = 0; etaBin < binAdm_->nEtaBins(); ++etaBin) {
      for(unsigned int ptBin = 0; ptBin < binAdm_->nPtBins(etaBin); ++ptBin) {

	// Eta bin borders
	sprintf(tmp,"  $%.1lf$ -- $%.1lf$ & $",binAdm_->etaMin(etaBin),binAdm_->etaMax(etaBin));
	text += tmp;
	// PtAve bin borders
	double x = binAdm_->ptMin(etaBin,ptBin);
	if( x < 1000 ) text += " ";
	if( x < 100  ) text += " ";
	sprintf(tmp,"%.0lf$ -+- $",x);
	text += tmp;
	x = binAdm_->ptMax(etaBin,ptBin);
	if( x < 1000 ) text += " ";
	if( x < 100  ) text += " ";
	sprintf(tmp,"%.0lf$ & $",x);
	text += tmp;
	// PtAve mean
	x = gScaleFactors_.at(etaBin)->GetX()[ptBin];
	if( x < 1000 ) text += " ";
	if( x < 100  ) text += " ";
	sprintf(tmp,"%.1lf \\pm %.1lf$ &  $",x,gScaleFactors_.at(etaBin)->GetEXlow()[ptBin]);
	text += tmp;
	// Scale factors
	sprintf(tmp,"%.3lf",gScaleFactors_.at(etaBin)->GetY()[ptBin]);
	text += tmp;
	double estat = gScaleFactors_.at(etaBin)->GetEYlow()[ptBin];
	double esystd = gAbsUncertTot_.at(etaBin)->GetEYlow()[ptBin];
// 	double esystu = gAbsUncertTot_.at(etaBin)->GetEYhigh()[ptBin];
// 	double etotd = sqrt( estat*estat + esystd*esystd );
// 	double etotu = sqrt( estat*estat + esystu*esystu );
	sprintf(tmp," \\pm %.3lf \\pm %.3lf$ \\\\\n",estat,esystd);
	text += tmp;
      }	// End of loop over ptAve bins
      if( etaBin < binAdm_->nEtaBins()-1 ) text += "\\midrule\n";
    } // End of loop over eta bins

    // LaTeX output
    TString latex = "% ***  Final Scale Factors ("+out_->namePrefix()+")  ***\n\n";
    latex += "\\renewcommand{\\arraystretch}{1.5}\n";
    latex += "\\begin{tabular}{cr@{ -- }rcc}\n\\toprule\n";
    latex += "$|\\eta|$ & \\multicolumn{2}{c}{$\\ptave\\,(\\gevnospace)$} & $\\mean{\\ptave}\\,(\\gevnospace)$ & $";
    sprintf(tmp,"\\tailratio(\\tailborder{%.1lf})$ \\\\\n\\midrule\n",fitPars_->nSigTailStart());
    latex += tmp;
    latex += text;
    latex += "\\bottomrule\n\\end{tabular}\n\\renewcommand{\\arraystretch}{1.}\n";
    latex.ReplaceAll("-+-","&");
    std::cout << "Writing final scale factors to " << out_->toFile(latex,"ScaleFactors.tex",false) << std::endl;
    
    // Screen output
    text.ReplaceAll("$","");
    text.ReplaceAll("&","");
    text.ReplaceAll("\\pm","+/-");
    text.ReplaceAll("\\\\","");
    text.ReplaceAll("-+-","--");
    text.ReplaceAll("\\midrule"," ----------------------------------------------------------------------");

    std::cout << "\n\n   Eta Bin      PtAve Bin     Mean PtAve     Scale     Stat      Syst" << std::endl;
    std::cout << " ----------------------------------------------------------------------" << std::endl;
    std::cout << text;
    std::cout << " ----------------------------------------------------------------------\n" << std::endl;

    if( Output::DEBUG ) std::cout << "FinalResultProducer::print() Leaving" << std::endl;
  }



  // Returns graph where the coordinates of each point have:
  // - y value: bin content of h1;
  // - x value: bin content of h2; and
  // - the errors, respectively.
  // Requires h1 and h2 to have the same number of bins.
  // ------------------------------------------------------------------------------------
  TGraphAsymmErrors* FinalResultProducer::correctBinCenter(const TH1* h1, const TH1* h2) const {
    
    std::vector<double> x;
    std::vector<double> xe;
    std::vector<double> y;
    std::vector<double> ye;
    for(int bin = 1; bin <= h1->GetNbinsX(); ++bin) {
      x.push_back(h2->GetBinContent(bin));
      xe.push_back(h2->GetBinError(bin));
      y.push_back(h1->GetBinContent(bin));
      ye.push_back(h1->GetBinError(bin));
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						 &(xe.front()),&(xe.front()),&(ye.front()),&(ye.front()));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(style_->markerSize());
    g->SetLineWidth(style_->lineWidth());
    
    return g;
  }
  


  // Returns graph of the relative difference (relD) between the bin contents
  // of hNom and hVarUp, hVarDn. The y value of each point is 0, the y error is
  // relD, the x value corresponds to the bin center of hNom, the x error to
  // half the bin width.
  // If
  // - hVarDn == 0: relD = |( Y(hNom) - Y(hVarUp) ) / Y(hNom)|, hVarDn is ignored
  // - hVarDn != 0: relD = 0.5 ( |( Y(hNom) - Y(hVarUp) ) / Y(hNom)| + |( Y(hNom) - Y(hVarDn) ) / Y(hNom)|)
  // Then, relD is multiplied by 'weight'.
  // 'color' specifies the fill color of the graph.
  // ------------------------------------------------------------------------------------
  TGraphAsymmErrors* FinalResultProducer::relativeDifference(const TH1* hNom, int color, double weight, const TH1* hVarUp, const TH1* hVarDn) const {
    
    std::vector<double> x;
    std::vector<double> xe;
    std::vector<double> y;
    std::vector<double> yeu;
    std::vector<double> yed;
    for(int bin = 1; bin <= hNom->GetNbinsX(); ++bin) {
      x.push_back(hNom->GetBinCenter(bin));
      xe.push_back(0.5*hNom->GetBinWidth(bin));
      y.push_back(0.);
      double nom = hNom->GetBinContent(bin);
      double var = hVarUp->GetBinContent(bin);
      double err = std::abs((var-nom)/nom);
      if( hVarDn ) {
	var = hVarDn->GetBinContent(bin);
	err += std::abs((var-nom)/nom);
	err /= 2.;
      }
      err *= weight;
      yeu.push_back(err);
      yed.push_back(err);
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						 &(xe.front()),&(xe.front()),&(yed.front()),&(yeu.front()));
    g->SetFillStyle(1001);
    g->SetFillColor(color);
    g->SetLineColor(color);
    
    return g;
  }



  // Returns graph with y errors summed in quadrature
  // ------------------------------------------------------------------------------------
  TGraphAsymmErrors* FinalResultProducer::sumErrorsInQuadrature(const std::vector<TGraphAsymmErrors*> &graphs) const {

    TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(graphs.at(0)->Clone());
    for(int i = 0; i < g->GetN(); ++i) {
      double eu2 = 0.;
      double ed2 = 0.;
      for(unsigned int j = 0; j < graphs.size(); ++j) {
	eu2 += pow(graphs.at(j)->GetEYhigh()[i],2); 
	ed2 += pow(graphs.at(j)->GetEYlow()[i],2);
      }
      if( ed2 == 0. ) ed2 = eu2;
      g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],sqrt(ed2),sqrt(eu2));
    }
    g->SetFillStyle(3444);
    g->SetFillColor(kBlack);
    g->SetLineColor(kBlack);
    
    
    return g;
  }



  // Create an error band for a given histogram and a given graph of the
  // relative uncertainties on the bin content.
  //
  // Returns a graph with each point having
  //          y       : bin content of 'hNom'
  // error on y       : bin content of 'hNom' multiplied by bin content of 'gRelUncert'
  // x and error on x : as in 'gRelUncert'
  // 'color' specifies the fill color of the error band.
  // ------------------------------------------------------------------------------------
  TGraphAsymmErrors* FinalResultProducer::uncertaintyBand(const TH1* hNom, const TGraphAsymmErrors* gRelUncert, int color) const {
    TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(gRelUncert->Clone());
    for(int i = 0; i < g->GetN(); ++i) {
      double y = hNom->GetBinContent(1+i);
      g->SetPoint(i,g->GetX()[i],y);
      double eu = (g->GetEYhigh()[i])*y;
      double ed = (g->GetEYlow()[i])*y;
      g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],ed,eu);
    }
    g->SetFillStyle(1001);
    g->SetFillColor(color);
    g->SetLineColor(color);
    
    return g;
  }



  // ------------------------------------------------------------------------------------
  std::vector<TGraphAsymmErrors*> FinalResultProducer::uncertaintyStack(const std::vector<TGraphAsymmErrors*> &gUncerts, unsigned int etaBin) const {
    std::vector<TGraphAsymmErrors*> stack;
    for(unsigned int i = 0; i < gUncerts.size(); ++i) {
      std::vector<double> x;
      std::vector<double> xed;
      std::vector<double> xeu;
      std::vector<double> y;
      std::vector<double> yed;
      std::vector<double> yeu;
      for(int n = 0; n < gUncerts.at(i)->GetN(); ++n) {
	x.push_back(gUncerts.at(i)->GetX()[n]);
	xed.push_back(gUncerts.at(i)->GetEXlow()[n]);
	xeu.push_back(gUncerts.at(i)->GetEXhigh()[n]);
	y.push_back(gUncerts.at(i)->GetY()[n]);
	yed.push_back(0.);
	double err = 100.*gUncerts.at(i)->GetEYhigh()[n];
	if( i > 0 ) {
	  double prevErr = stack.back()->GetEYhigh()[n];
	  err = sqrt( err*err + prevErr*prevErr );
	}
	yeu.push_back(err);
      }
      x.insert(x.begin(),binAdm_->ptMin(etaBin,0));
      xed.insert(xed.begin(),0.);
      xeu.insert(xeu.begin(),0.);
      y.insert(y.begin(),y.front());
      yed.insert(yed.begin(),0.);
      yeu.insert(yeu.begin(),yeu.front());
	
      x.push_back(binAdm_->ptMax(etaBin,binAdm_->nPtBins(etaBin)-1));
      xed.push_back(0.);
      xeu.push_back(0.);
      y.push_back(y.back());
      yed.push_back(0.);
      yeu.push_back(yeu.back());
	
      TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),&(xed.front()),&(xeu.front()),&(yed.front()),&(yeu.front()));
      g->SetLineWidth(0);
      g->SetFillStyle(1001);
      g->SetFillColor(gUncerts.at(i)->GetFillColor());
      g->SetLineColor(g->GetFillColor());
      stack.push_back(g);
    }

    return stack;
  }
}
#endif
