// $Id: Extrapolation.cc,v 1.16 2012/06/07 21:10:55 mschrode Exp $

#include "Extrapolation.h"

#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TGraphAsymmErrors.h"


namespace resolutionFit {

  unsigned int Extrapolation::NUM_EXTRAPOLATION_FUNCTIONS = 0;
  unsigned int Extrapolation::NUM_INSTANCES = 0;
  

  //! @param wpIdx Index of stable working-point. If negative, use simple method
  //!              where y-axis intercept is taken as resolution and error on this
  //!              as uncertainty.
  // -------------------------------------------------------------------------------------
  Extrapolation::Extrapolation(double minPt, double maxPt, double minPt3, int wpIdx)
    : minPt_(minPt), maxPt_(maxPt), minPt3_(minPt3),
      useWPExtrapolation_(wpIdx<0 ? false : true), wpIdx_(wpIdx_<0 ? 1000 : static_cast<int>(wpIdx))
 {
    
    if( NUM_INSTANCES == 0 ) {
      std::cout << "\n\n++++++ Extrapolation +++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "\nLinear extrapolation\nMinimum ptSoft = " << minPt3_ << " GeV" << std::endl;
      if( useWPExtrapolation_ ) std::cout << "Using 'Working-Point Extrapolation' with WP = " << wpIdx_ << std::endl;
      std::cout << std::endl;
    }
    
    ++NUM_INSTANCES;
  }


  // -------------------------------------------------------------------------------------
  bool Extrapolation::operator()(const std::vector<double> &ptSoft,
				 const std::vector<double> &values,
				 const std::vector<double> &uncerts,
				 TF1* &fit, unsigned int &firstPointInExtrapolation,
				 double &extra, double &statUncert, double &systUncert) const {
    bool result = true;
    
    unsigned int effWPIdx = wpIdx_;
    TGraphAsymmErrors* g = getGraph(ptSoft,values,uncerts,effWPIdx);
    
    // Setup linear fit
    TString name = "LinearExtrapolation";
    name += NUM_EXTRAPOLATION_FUNCTIONS;
    ++NUM_EXTRAPOLATION_FUNCTIONS;
    fit = new TF1(name,"pol1",g->GetX()[0],g->GetX()[g->GetN()-1]);
    fit->SetLineWidth(1);
    bool fitResult = !(g->Fit(fit,"0QR"));

    // HACK for Eta4 bin
    // Should be replace by proper chi2 criterion
    if( fit->GetParameter(0) > 34 && g->GetY()[2] > 45 ) {
      std::cout << "\n\n************* HACK IN EXTRAPOLATION *********************\n\n" << std::endl;
      for(int i = 0; i < 7; ++i) {
	std::cout << "  ptSoftRel = " << g->GetX()[0] << " ("  << g->GetX()[0]*minPt_ << " GeV)" << std::endl;
	g->RemovePoint(0);
	if( effWPIdx > 0 ) effWPIdx--;
      }
      fitResult = !(g->Fit(fit,"0QR"));
    }

    firstPointInExtrapolation = ptSoft.size() - g->GetN();
     
    // Work-around to keep the program running in case of 
    // fitting failure
    if( !fitResult ) {
      double mean = 0.;
      for(int i = 0; i < g->GetN(); ++i) {
	mean += g->GetY()[i];
      }
      mean /= g->GetN();
      fit->SetParameter(0,mean);
      fit->SetParError(0,100.);
      fit->SetParameter(1,0.);
      fit->SetParError(1,10.);
      fitResult = true;
      std::cerr << "\n\n+++++ WARNING: Fit failed in " << bin() << " +++++++++++++++++++++++++++++++++++++++\n";
      std::cerr << "  Cured by setting fit result to mean " << mean << " with large error " << 100 << std::endl << std::endl;
    }

    if( useWPExtrapolation_ ) {
      // Working point * slope
      //std::cout << "   WP = " << g->GetX()[effWPIdx] << std::endl;
      extra = g->GetY()[effWPIdx] - fit->GetParameter(1)*g->GetX()[effWPIdx];
      double statPoint = g->GetEYhigh()[effWPIdx];
      double statSlope = g->GetX()[effWPIdx]*fit->GetParError(1);
      statUncert = sqrt( statPoint*statPoint + statSlope*statSlope );
      systUncert = 0.5*std::abs(fit->Eval(g->GetX()[0])-fit->Eval(0.));
    } else {
      // sigma = y(0)
      // syst uncertainty wrt to extrapolation fit
      extra = fit->GetParameter(0);
      statUncert = fit->GetParError(0);
      systUncert = 0.5*std::abs(fit->Eval(g->GetX()[0])-fit->Eval(0));
    }

    return fitResult;
  }



  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* Extrapolation::getGraph(const std::vector<double> &ptSoft,
					     const std::vector<double> &values,
					     const std::vector<double> &uncerts,
					     unsigned int &effectiveWPIdx) const {
    std::vector<double> ptSoftUncert(ptSoft.size(),0.);
    TGraphAsymmErrors* g = new TGraphAsymmErrors(ptSoft.size(),&(ptSoft.front()),&(values.front()),
						 &(ptSoftUncert.front()),&(ptSoftUncert.front()),
						 &(uncerts.front()),&(uncerts.front()));
    
    // Remove points that correspond to bins with
    // ptAveMin below 'minPt3_'
    effectiveWPIdx = wpIdx_;
    double minPt3 = minPt3_;
    if( minPt3 > 6. && minPt_ < 70. ) minPt3 = 6.; // Otherwise no points left in first bin
    int nPointsToDelete = 0;
    for(int i = 0; i < g->GetN(); ++i) {
      if( g->GetX()[i]*minPt_ < minPt3 ) nPointsToDelete++;
    }
    if( nPointsToDelete ) {
      std::cout << "  Extrapolation in " << bin() << ": Removing " << nPointsToDelete << " points" << std::endl;
      for(int i = 0; i < nPointsToDelete; ++i) {
	std::cout << "    ptSoftRel = " << g->GetX()[0] << " ("  << g->GetX()[0]*minPt_ << " GeV)" << std::endl;
	g->RemovePoint(0);
      }
      effectiveWPIdx -= nPointsToDelete;
      std::cout << "    " << g->GetN() << " points remaining" << std::endl;
    }

    return g;
  }
  

  // -------------------------------------------------------------------------------------
  TString Extrapolation::bin() const {
    TString str = "pt bin ";
    str += minPt_;
    str += " - ";
    str += maxPt_;
    str += " GeV";

    return str;
  }
}
