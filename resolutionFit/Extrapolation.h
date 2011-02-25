// $Id: Extrapolation.h,v 1.2 2011/02/17 13:42:32 mschrode Exp $

#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include <vector>

#include "TF1.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class Extrapolation {
  public:
    Extrapolation(double meanPt);

    bool operator()(const std::vector<double> &values,
		    const std::vector<double> &uncerts,
		    const std::vector<double> &ptSoftconst,
		    TF1* &fit, double &systUncert) const;


  private:
    static unsigned int NUM_EXTRAPOLATION_FUNCTIONS;

    const double meanPt_;
  };
}
#endif
