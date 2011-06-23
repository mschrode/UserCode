// $Id: Extrapolation.cc,v 1.11 2011/06/08 16:58:02 mschrode Exp $

#include "Extrapolation.h"

#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TGraphAsymmErrors.h"


namespace resolutionFit {

  unsigned int Extrapolation::NUM_EXTRAPOLATION_FUNCTIONS = 0;
  
  
  // -------------------------------------------------------------------------------------
  bool Extrapolation::operator()(const std::vector<double> &ptSoft,
				 const std::vector<double> &values,
				 const std::vector<double> &uncerts,
				 TF1* &fit, double &systUncert) const {
    bool result = true;
    
    TGraphAsymmErrors* g = getGraph(ptSoft,values,uncerts);

    // Setup linear fit
    TString name = "LinearExtrapolation";
    name += NUM_EXTRAPOLATION_FUNCTIONS;
    ++NUM_EXTRAPOLATION_FUNCTIONS;
    fit = new TF1(name,"pol1",g->GetX()[0],g->GetX()[g->GetN()-1]);
    fit->SetLineWidth(1);
    bool fitResult = !(g->Fit(fit,"0QR"));
     
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

    systUncert = 0.5*std::abs(g->GetY()[0]-fit->GetParameter(0));


    return fitResult;
  }



  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* Extrapolation::getGraph(const std::vector<double> &ptSoft,
					     const std::vector<double> &values,
					     const std::vector<double> &uncerts) const {
    std::vector<double> ptSoftUncert(ptSoft.size(),0.);
    TGraphAsymmErrors* g = new TGraphAsymmErrors(ptSoft.size(),&(ptSoft.front()),&(values.front()),
						 &(ptSoftUncert.front()),&(ptSoftUncert.front()),
						 &(uncerts.front()),&(uncerts.front()));
    
    // Remove points that correspond to bins with an
    // average ptAve below 10 GeV
    double minPtAve = minPt_ > 70. ? 10. : 6.;
    int nPointsToDelete = 0;
    for(int i = 0; i < g->GetN(); ++i) {
      if( g->GetX()[i]*minPt_ < minPtAve ) nPointsToDelete++;
    }
    if( nPointsToDelete ) {
      std::cout << "Extrapolation in " << bin() << ": Removing " << nPointsToDelete << " points" << std::endl;
      for(int i = 0; i < nPointsToDelete; ++i) {
	std::cout << "  ptSoftRel = " << g->GetX()[0] << " ("  << g->GetX()[0]*minPt_ << " GeV)" << std::endl;
	g->RemovePoint(0);
      }
      std::cout << "  " << g->GetN() << " points remaining" << std::endl;
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
