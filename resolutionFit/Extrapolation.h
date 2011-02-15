// $Id: $

#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include <vector>

#include "TF1.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class Extrapolation {
  public:
    Extrapolation() {};

    virtual bool operator()(double meanPt,
			    const std::vector<double> &values,
			    const std::vector<double> &uncerts,
			    const std::vector<double> &ptSoftconst,
			    TF1* &fit) const;


  private:
    static unsigned int NUM_EXTRAPOLATION_FUNCTIONS;
  };
}
#endif
