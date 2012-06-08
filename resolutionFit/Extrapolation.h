// $Id: Extrapolation.h,v 1.7 2012/06/07 21:10:55 mschrode Exp $

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
    Extrapolation(double minPt, double maxPt, double minPt3, int wpIdx);
    
    bool operator()(const std::vector<double> &values,
		    const std::vector<double> &uncerts,
		    const std::vector<double> &ptSoftconst,
		    TF1* &fit, unsigned int &firstPointInExtrapolation,
		    double &extra, double &statUncert, double &systUncert) const;


  private:
    static unsigned int NUM_EXTRAPOLATION_FUNCTIONS;
    static unsigned int NUM_INSTANCES;

    const double minPt_;
    const double maxPt_;
    const double minPt3_;
    const bool useWPExtrapolation_;
    const unsigned int wpIdx_;

    TGraphAsymmErrors* getGraph(const std::vector<double> &ptSoft,
				const std::vector<double> &values,
				const std::vector<double> &uncerts,
				unsigned int &effectiveWPIdx) const;
    TString bin() const;
  };
}
#endif
