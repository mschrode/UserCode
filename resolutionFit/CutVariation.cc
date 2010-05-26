// $Id: CutVariation.cc,v 1.9 2010/05/18 12:04:45 mschrode Exp $

#include "CutVariation.h"
#include "KalibriFileParser.h"

#include <algorithm>
#include <cmath>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"

namespace resolutionFit {
  CutVariation::CutVariation(const Parameters::PtBinParameters *par, int parIndex)
    : par_(par), parIdx_(parIndex) {
    // Read uncertainty from MC statistics
    mcStatUncert_ = 0.;
    if( par_->hasMCStatUncert() ) {
      KalibriFileParser *parser = new KalibriFileParser(par_->fileNameMCStatUncert(),par_->verbosity());
      mcStatUncert_ = parser->statUncert(parIdx());
      delete parser;
    }
    // Read values from file
    varPoints_ = std::vector<VariationPoint*>(nCutValues());
    for(int i = 0; i < nCutValues(); i++) {
      KalibriFileParser *parser = new KalibriFileParser(par_->fileNamePt3CutVariations(i),par_->verbosity());
      // Mean pt for all varied cuts (set only once)
      if( i == 0 ) {
	meanPt_ = parser->meanPt();
	if( isRelValue() ) mcStatUncert_ /= meanPt_;
      }
      // Create value at cut variation i
      double fittedValue = parser->value(parIdx());
      double statUncert = parser->statUncert(parIdx());
      if( isRelValue() ) {
	fittedValue /= meanPt_;
	statUncert /= meanPt_;
      }
      statUncert = sqrt( statUncert*statUncert + mcStatUncert_*mcStatUncert_ );
      Uncertainty *uncert = new Uncertainty("Statistical uncertainty",statUncert);
      delete parser;
      varPoints_[i] = new VariationPoint(fittedValue,uncert,cutValue(i));
    }
    // Create graph of varied values
    createTGraph();
    // Default value for extrapolation
    extrapolatedPoint_ = new VariationPoint();
    // Initialize fit function
    TString name = "resolutionFit::CutVariationFit_PtBin";
    name += par_->ptBinIdx();
    name += "_Par";
    name += parIdx();
    fit_ = new TF1(name,"pol1",cutValue(0),cutValue(nCutValues()-1));
    fit_->SetParameter(0,0.);
    fit_->SetParameter(1,0.);
    fit_->SetLineColor(2);
  }


  CutVariation::~CutVariation() {
    for(std::vector<VariationPoint*>::iterator it = varPoints_.begin();
	it != varPoints_.end(); it++) {
      delete *it;
    }
    delete fit_;
    delete graph_;
    delete extrapolatedPoint_;
  }

  TF1 *CutVariation::getTF1(const TString &name) const {
    TF1 *fit = new TF1(name,"pol1",0.,1.4*maxCutValue());
    for(int i = 0; i < 2; i++) {
      fit->SetParameter(i,fit_->GetParameter(i));
      fit->SetParError(i,fit_->GetParError(i));
    }
    fit->SetLineColor(2);
    return fit;
  }

  TGraphAsymmErrors *CutVariation::getTGraph() const {
    TGraphAsymmErrors *graph = static_cast<TGraphAsymmErrors*>(graph_->Clone());
    return graph;
  }

  TH1 *CutVariation::getFrame(const TString &name) const { 
    double xMin = 0.;
    double xMax = 1.4*maxCutValue();
    double yMin = fit_->GetParameter(0);
    double yMax = fit_->GetParameter(1)*xMax + fit_->GetParameter(0);
    if( fit_->GetParameter(1) < 0 ) {
      yMin = fit_->GetParameter(1)*xMax + fit_->GetParameter(0);
      yMax = fit_->GetParameter(0);
    }

    double gMin = *(std::min_element(graph_->GetY(),graph_->GetY()+graph_->GetN()));
    if( gMin < yMin ) yMin = gMin;
    double gMax = *(std::max_element(graph_->GetY(),graph_->GetY()+graph_->GetN()));
    if( gMax > yMax ) yMax = gMax;    

    double deltaY = yMax - yMin;
    yMin -= 0.2*deltaY;
    if( yMin < 0. ) yMin = 0.;
    yMax += 0.8*deltaY;

    TH1 *hFrame = new TH1D(name,";p^{rel}_{T,3} cut;",1000,xMin,xMax);
    hFrame->SetYTitle(par_->parAxisLabel(parIdx()));
    hFrame->GetYaxis()->SetRangeUser(yMin,yMax);
    hFrame->SetNdivisions(505);
    return hFrame;
  }

  void CutVariation::extrapolate() {
    if( par_->verbosity() == 2 ) {
      std::cout << "CutVariation: Fitting extrapolation" << std::endl;
    }
    if( nCutValues() >= 2 ) {
      // Fit varied values
      graph_->Fit(fit_,"0Q");
      // Replace extrapolated value with fit results
      delete extrapolatedPoint_;
      Uncertainty *uncert = new Uncertainty("Extrapolation",fit_->GetParError(0));
      extrapolatedPoint_ = new VariationPoint(fit_->GetParameter(0),uncert,0.);
    } else {
      std::cerr << "  WARNING: Less than 2 cut variations" << std::endl;
    }
  }


  void CutVariation::createTGraph() {
    if( par_->verbosity() == 2 ) {
      std::cout << "CutVariation: Creating graph of varied values" << std::endl;
    }
    std::vector<double> x(nCutValues());
    std::vector<double> ex(nCutValues(),0.);
    std::vector<double> y(nCutValues());
    std::vector<double> ey(nCutValues());
    for(int i = 0; i < nCutValues(); i++) {
      x[i] = cutValue(i);
      y[i] = fittedValue(i);
      ey[i] = uncert(i);
    }   
    graph_ = new TGraphAsymmErrors(nCutValues(),&(x.front()),&(y.front()),
				   &(ex.front()),&(ex.front()),
				   &(ey.front()),&(ey.front()));
    graph_->SetMarkerStyle(20);
  }


  CutVariation::VariationPoint::VariationPoint(double fitValue, const Uncertainty *uncert, double cutValue)
    : fitValue_(fitValue), uncert_(uncert), cutValue_(cutValue) {};

  CutVariation::VariationPoint::VariationPoint()
    : fitValue_(0.), uncert_(new Uncertainty()), cutValue_(0.) {};

  CutVariation::VariationPoint::~VariationPoint() {
    delete uncert_;
  }
}
