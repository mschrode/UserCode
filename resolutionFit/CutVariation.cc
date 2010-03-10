// $Id: $

#include "CutVariation.h"

#include "KalibriFileParser.h"

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"

namespace resolutionFit {
  int CutVariation::nObjs = 0;

  CutVariation::CutVariation(const std::vector<TString> &fileNames, const std::vector<double> &cutValues, int verbose)
    : verbose_(verbose) {
    // Read values from file
    varPoints_ = std::vector<VariationPoint*>(fileNames.size());
    for(int i = 0; i < nCutValues(); i++) {
      KalibriFileParser *parser = new KalibriFileParser(fileNames[i],verbose_);
      // Mean pt for all varied cuts (set only once)
      if( i == 0 ) meanPt_ = parser->meanPtGen();//parser->meanPdfPtTrue();
      // Create value at cut variation i
      double relSigma = parser->value()/meanPt_;
      Uncertainty *uncert = new Uncertainty("Statistical uncertainty",parser->statUncert()/meanPt_);
      varPoints_[i] = new VariationPoint(relSigma,uncert,cutValues[i]);
      delete parser;
    }
    // Create graph of varied values
    createTGraph();
    // Default value for extrapolation
    extrapolatedPoint_ = new VariationPoint();
    // Initialize fit function
    TString name = "resolutionFit::CutVariationFit";
    name += CutVariation::nObjs;
    fit_ = new TF1(name,"pol1",cutValue(0),cutValue(nCutValues()-1));
    fit_->SetParameter(0,0.);
    fit_->SetParameter(1,0.);
    fit_->SetLineColor(2);

    // Update number of created CutVariation objects
    CutVariation::nObjs++;
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

  TH1D *CutVariation::getFrame(const TString &name) const { 
    TH1D *hFrame = new TH1D(name,";p^{3}_{T,rel};#sigma / p_{T}",1000,0.,1.4*maxCutValue());
    hFrame->GetYaxis()->SetRangeUser(0.8*relSigma(0),1.2*relSigma(nCutValues()-1));
    return hFrame;
  }

  void CutVariation::extrapolate() {
    if( verbose_ == 1 ) {
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
    if( verbose_ == 2 ) {
      std::cout << "CutVariation: Creating graph of varied values" << std::endl;
    }
    std::vector<double> x(nCutValues());
    std::vector<double> ex(nCutValues(),0.);
    std::vector<double> y(nCutValues());
    std::vector<double> ey(nCutValues());
    for(int i = 0; i < nCutValues(); i++) {
      x[i] = cutValue(i);
      y[i] = relSigma(i);
      ey[i] = uncert(i);
    }   
    graph_ = new TGraphAsymmErrors(nCutValues(),&(x.front()),&(y.front()),
				   &(ex.front()),&(ex.front()),
				   &(ey.front()),&(ey.front()));
    graph_->SetMarkerStyle(20);
  }


  CutVariation::VariationPoint::VariationPoint(double val, const Uncertainty *uncert, double cutValue)
    : relSigma_(val), uncert_(uncert), cutValue_(cutValue) {};

  CutVariation::VariationPoint::VariationPoint()
    : relSigma_(0.), uncert_(new Uncertainty()), cutValue_(0.) {};

  CutVariation::VariationPoint::~VariationPoint() {
    delete uncert_;
  }
}
