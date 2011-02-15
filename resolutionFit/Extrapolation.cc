// $Id: $

#include "Extrapolation.h"

#include <iostream>

#include "TGraphAsymmErrors.h"
#include "TString.h"


namespace resolutionFit {

  unsigned int Extrapolation::NUM_EXTRAPOLATION_FUNCTIONS = 0;

  // -------------------------------------------------------------------------------------
  bool Extrapolation::operator()(double meanPt,
				 const std::vector<double> &ptSoft,
				 const std::vector<double> &values,
				 const std::vector<double> &uncerts,
				 TF1* &fit) const {
    bool result = true;

    std::vector<double> ptSoftUncert(ptSoft.size(),0.);
    TGraphAsymmErrors* g = new TGraphAsymmErrors(ptSoft.size(),&(ptSoft.front()),&(values.front()),
						 &(ptSoftUncert.front()),&(ptSoftUncert.front()),
						 &(uncerts.front()),&(uncerts.front()));

//     int nPointsToDelete = 0;
//     for(int i = 0; i < g->GetN(); ++i) {
//       if( g->GetX()[i]*meanPt < 6. ) nPointsToDelete++;
//     }
//     for(int i = 0; i < nPointsToDelete; ++i) {
//       g->RemovePoint(0);
//     }
//    if( nPointsToDelete ) std::cout << "Deleted " << nPointsToDelete << " points: " << g->GetN() << std::endl;

    TString name = "LinearExtrapolation";
    name += NUM_EXTRAPOLATION_FUNCTIONS;
    ++NUM_EXTRAPOLATION_FUNCTIONS;
    fit = new TF1(name,"pol1",0.,1.4*ptSoft.back());
    fit->SetLineWidth(1);

    return !( g->Fit(fit,"0QR") );
  }
}
