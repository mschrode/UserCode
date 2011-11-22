// $Id: Extrapolation.h,v 1.5 2011/11/21 17:18:04 mschrode Exp $

#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include <vector>

#include "TString.h"

class TF1;
class TGraphAsymmErrors;

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class Extrapolation {
  public:
    Extrapolation(double minPt, double maxPt, double minPt3);
    
    bool operator()(const std::vector<double> &values,
		    const std::vector<double> &uncerts,
		    const std::vector<double> &ptSoftconst,
		    TF1* &fit, double &systUncert) const;


  private:
    static unsigned int NUM_EXTRAPOLATION_FUNCTIONS;
    static unsigned int NUM_INSTANCES;

    const double minPt_;
    const double maxPt_;
    const double minPt3_;

    TGraphAsymmErrors* getGraph(const std::vector<double> &ptSoft,
				const std::vector<double> &values,
				const std::vector<double> &uncerts) const;
    TString bin() const;
  };
}
#endif
