// $Id: CutVariation.cc,v 1.24 2010/12/27 13:34:05 mschrode Exp $

#include "CutVariation.h"
#include "KalibriFileParser.h"

#include <algorithm>
#include <cmath>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"

#include "../util/utils.h"


namespace resolutionFit {
  CutVariation::CutVariation(const Parameters::PtBinParameters *par, int parIndex, const TString &quantity)
    : par_(par), parIdx_(parIndex) {

    mcStatUncert_ = 0.;

    if( quantity == "MaxLikeFit" ) {
      // Read uncertainty from MC statistics
      if( par_->hasMCStatUncert() ) {
	KalibriFileParser *parser = new KalibriFileParser(par_->fileNameMCStatUncert(),par_->ptBinIdx(),par_->verbosity());
	mcStatUncert_ = parser->statUncert(parIdx());
	delete parser;
      }
      // Read values from file
      varPoints_ = std::vector<VariationPoint*>(nPt3Cuts());
      for(int i = 0; i < nPt3Cuts(); i++) {
	KalibriFileParser *parser = new KalibriFileParser(par_->fileNamePt3CutVariations(i),par_->ptBinIdx(),par_->verbosity());
	// Mean pt for all varied cuts
	// Set only once from lowest ptsoft bin
	// to have smallest bias
	if( i == 0 ) {
	  if( par_->binPt() == BinPtGen ) { 
	    meanPt_ = parser->meanPtGen();
	    meanPtUncert_ = parser->meanPtGenUncert();
	  } else if( par_->binPt() == BinPtAve ) {
	    if( par_->fitMode() == FitModeMaxLikeFull ) {
	      meanPt_ = parser->meanPt();
	      meanPtUncert_ = parser->meanPtAveUncert(); // TODO: Use correct uncertainty
	    } else if( par_->fitMode() == FitModeMaxLikeSimple ) {
	      meanPt_ = parser->meanPtAve();
	      meanPtUncert_ = parser->meanPtAveUncert();
	    }
	  } else {
	    meanPt_ = 1.;
	    meanPtUncert_ = 1.;
	  }
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
	Uncertainty *uncert = 0;
// 	if( par_->isData() ) uncert = new Uncertainty("Statistical uncertainty",statUncert*sqrt(32.));
// 	else uncert = new Uncertainty("Statistical uncertainty",statUncert);
	uncert = new Uncertainty("Statistical uncertainty",statUncert);
	
	double pt3Val = pt3Max(i);
	if( pt3Bins() ) pt3Val = pt3Mean(i);

	TH1 *hMCRes = parser->hist("hRespMeas","resolutionFit::CutVariationMCRes_PtBin"+util::toTString(par_->ptBinIdx())+"_Var"+util::toTString(i));
	delete parser;

	varPoints_[i] = new VariationPoint(fittedValue,uncert,pt3Val,hMCRes);
      }
    } else if( quantity == "PtAsym" ) {
      // Read asymmetry distributions and mean pt from file
      varPoints_ = std::vector<VariationPoint*>(nPt3Cuts());
      for(int i = 0; i < nPt3Cuts(); i++) {
	KalibriFileParser *parser = new KalibriFileParser(par_->fileNamePt3CutVariations(i),par_->ptBinIdx(),par_->verbosity(),false);
	// Mean pt for all varied cuts (set only once)
	if( i == 0 ) {
	  if( par_->binPt() == BinPtGen ) { 
	    meanPt_ = parser->meanPtGen();
	    meanPtUncert_ = parser->meanPtGenUncert();
	  } else {
	    meanPt_ = parser->meanPtAve();
	    meanPtUncert_ = parser->meanPtAveUncert();
	  }
	}	      
	
	// Fit asymmetry distribution
	TString name = "resolutionFit::CutVariationPtAsym_PtBin";
	name += par_->ptBinIdx();
	name += "_Var";
	name += i;
	TH1 *hAsym = parser->hist("hPtAsym",name);
	double pt3Val = pt3Max(i);
	if( pt3Bins() ) pt3Val = pt3Mean(i);
	varPoints_[i] = new VariationPoint(hAsym,true,pt3Val);
	
	delete parser;
      }
    } else if( quantity == "PtGenAsym" ) {
      // Read ptGen asymmetry distributions and mean pt from file
      varPoints_ = std::vector<VariationPoint*>(nPt3Cuts());
      for(int i = 0; i < nPt3Cuts(); i++) {
	KalibriFileParser *parser = new KalibriFileParser(par_->fileNamePt3CutVariations(i),par_->ptBinIdx(),par_->verbosity(),false);
	// Mean pt for all varied cuts (set only once)
	if( i == 0 ) {
	  meanPt_ = parser->meanPtGen();
	  meanPtUncert_ = parser->meanPtGenUncert();
	}	      
	
	// Fit ptGen asymmetry distribution
	TString name = "resolutionFit::CutVariationPtGenAsym_PtBin";
	name += par_->ptBinIdx();
	name += "_Var";
	name += i;
	TH1 *hAsym = parser->hist("hPtGenAsym",name);
	double pt3Val = pt3Max(i);
	if( pt3Bins() ) pt3Val = pt3Mean(i);
	varPoints_[i] = new VariationPoint(hAsym,false,pt3Val);
	
	delete parser;
      }
    } else {
      std::cerr << "ERROR: Unknown extrapolation quantity '" << quantity << "'" << std::endl;
      exit(-1);
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

//     if( pt3Bins() ) fit_ = new TF1(name,"pol1",pt3Mean(0),pt3Mean(nPt3Cuts()-1));
//     else fit_ = new TF1(name,"pol1",pt3Max(0),pt3Max(nPt3Cuts()-1));

    fit_ = new TF1(name,"pol1",pt3Max(0),pt3Max(nPt3Cuts()-1));
    fit_->SetParameter(0,0.);
    fit_->SetParameter(1,0.);
    fit_->SetLineColor(2);

//     fit_ = new TF1(name,"sqrt( sq([0]) + sq([1]*x) )",pt3Max(0),pt3Max(nPt3Cuts()-1));
//     fit_->SetParameter(0,0.1);
//     fit_->SetParameter(1,0.01);
//     fit_->SetLineColor(2);
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
//     TF1 *fit = new TF1(name,"pol1",0.,1.4*pt3Max(nPt3Cuts()-1));
//     for(int i = 0; i < 2; i++) {
//       fit->SetParameter(i,fit_->GetParameter(i));
//       fit->SetParError(i,fit_->GetParError(i));
//     }

    TF1* fit = static_cast<TF1*>(fit_->Clone(name));
    fit->SetRange(0.,1.4*pt3Max(nPt3Cuts()-1));
    fit->SetLineColor(2);
    return fit;
  }

  TGraphAsymmErrors *CutVariation::getTGraph() const {
    TGraphAsymmErrors *graph = static_cast<TGraphAsymmErrors*>(graph_->Clone());
    return graph;
  }

  TH1 *CutVariation::getFrame(const TString &name) const { 
    double xMin = 0.;
    double xMax = 1.4*pt3Max(nPt3Cuts()-1);
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
    yMin -= 0.4*deltaY;
    if( yMin < 0. ) yMin = 0.;
    yMax += 1.2*deltaY;

    TH1 *hFrame = new TH1D(name,"",1000,xMin,xMax);
    hFrame->SetNdivisions(505);
    TString xTitle;
    if( par_->pt3Var() == Pt3Rel ) xTitle = "p^{rel}_{T,3}";
    else if( par_->pt3Var() == Pt3Abs ) xTitle = "p_{T,3}";
    else if( par_->pt3Var() == Pp3Rel ) xTitle = "p_{||,3} / <p^{ave}_{T}>";
    if( !pt3Bins() ) xTitle += " threshold";
    if( par_->pt3Var() == Pt3Abs ) xTitle += " (GeV)";
    hFrame->SetXTitle(xTitle);
    hFrame->SetYTitle(par_->parAxisLabel(parIdx()));
    hFrame->GetYaxis()->SetRangeUser(yMin,yMax);

    return hFrame;
  }

  void CutVariation::extrapolate() {
    if( par_->verbosity() == 2 ) {
      std::cout << "CutVariation: Fitting extrapolation" << std::endl;
    }
    if( nPt3Cuts() >= 3 ) {
      // Fit varied values
      graph_->Fit(fit_,"0Q");
      // Replace extrapolated value with fit results
      delete extrapolatedPoint_;
      Uncertainty *uncert = new Uncertainty("Extrapolation",fit_->GetParError(0));

      double val = fit_->GetParameter(0);
      if( par_->systExtrapolation() == "low" ) {
 	double diff = std::abs(fittedValue(0)-val);
 	val = val - 0.5*diff;
      } else if( par_->systExtrapolation() == "high" ) {
 	double diff = std::abs(fittedValue(0)-val);
 	val = val + 0.5*diff;
      }
      extrapolatedPoint_ = new VariationPoint(val,uncert,0.,0);

    } else {
      std::cerr << "  WARNING: Less than 3 cut variations" << std::endl;
    }
  }


  void CutVariation::createTGraph() {
    if( par_->verbosity() == 2 ) {
      std::cout << "CutVariation: Creating graph of varied values" << std::endl;
    }
    std::vector<double> x(nPt3Cuts());
    std::vector<double> ex(nPt3Cuts(),0.);
    std::vector<double> y(nPt3Cuts());
    std::vector<double> ey(nPt3Cuts());
    for(int i = 0; i < nPt3Cuts(); i++) {
      if( pt3Bins() ) x[i] = pt3Mean(i);
      else x[i] = pt3Max(i);
      y[i] = fittedValue(i);
      ey[i] = uncert(i);
    }   
    graph_ = new TGraphAsymmErrors(nPt3Cuts(),&(x.front()),&(y.front()),
				   &(ex.front()),&(ex.front()),
				   &(ey.front()),&(ey.front()));
    graph_->SetMarkerStyle(20);
  }


  CutVariation::VariationPoint::VariationPoint(double fitValue, Uncertainty *uncert, double cutValue, TH1 *hMCRes)
    : fitValue_(fitValue), uncert_(uncert), cutValue_(cutValue), hPtAsym_(0), hMCRes_(hMCRes) {};

  CutVariation::VariationPoint::VariationPoint(TH1 *hPtAsym, bool fitGauss, double cutValue)
    : cutValue_(cutValue), hPtAsym_(hPtAsym), hMCRes_(0) {
    TF1 *fit = new TF1("fit","gaus",hPtAsym_->GetMean()-2.*hPtAsym_->GetRMS(),hPtAsym_->GetMean()+2.*hPtAsym_->GetRMS());
    fit->SetParameter(0,5.*hPtAsym_->GetRMS()*hPtAsym_->GetEntries());
    fit->SetParameter(1,0.);
    fit->SetParameter(2,hPtAsym_->GetRMS());
    
    if( hPtAsym_->Fit(fit,"Q0IRB") != 0 ) {
      std::cerr << "WARNING: No convergence" << std::endl;
    }
    fitValue_ = sqrt(2.)*std::abs(fit->GetParameter(2));
    uncert_ = new Uncertainty("Statistical uncertainty",sqrt(2.)*fit->GetParError(2));
    delete fit;


//     double rms = hPtAsym_->GetRMS();
//     std::cout << "Before: " << sqrt(2)*rms << std::endl;
//     rms *= 2.7;
//     hPtAsym_->GetXaxis()->SetRangeUser(-rms,rms);
//     fitValue_ = sqrt(2.)*hPtAsym_->GetRMS();
//     std::cout << "After: " << fitValue_ << std::endl;
//     uncert_ = new Uncertainty("Statistical uncertainty",sqrt(2.)*hPtAsym_->GetRMSError());
  }

  CutVariation::VariationPoint::VariationPoint()
    : fitValue_(0.), uncert_(new Uncertainty()), cutValue_(0.), hPtAsym_(0), hMCRes_(0) {};

  CutVariation::VariationPoint::~VariationPoint() {
    delete uncert_;
    if( hPtAsym_ ) delete hPtAsym_;
    if( hMCRes_ ) delete hMCRes_;
  }


  TH1 *CutVariation::VariationPoint::histPtAsym(const TString &name) const {
    return hPtAsym_ ? static_cast<TH1D*>(hPtAsym_->Clone(name)) : 0;
  }


  TH1 *CutVariation::VariationPoint::histMCRes(const TString &name) const {
    return hMCRes_ ? static_cast<TH1D*>(hMCRes_->Clone(name)) : 0;
  }
}
