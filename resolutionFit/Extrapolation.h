// $Id: Extrapolation.h,v 1.8 2012/06/08 21:14:44 mschrode Exp $

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
    Extrapolation(double pt3Min, int wpIdx, const TString &sampleLabel, double etaMin, double etaMax, double ptMin, double ptMax);
    
    bool operator()(const std::vector<double> &values,
		    const std::vector<double> &uncerts,
		    const std::vector<double> &ptSoftconst,
		    TF1* &fit, unsigned int &firstPointInExtrapolation,
		    double &extra, double &statUncert, double &systUncert) const;


  private:
    static unsigned int NUM_EXTRAPOLATION_FUNCTIONS;
    static unsigned int NUM_INSTANCES;

    const TString sampleLabel_;
    const double etaMin_;
    const double etaMax_;
    const double ptMin_;
    const double ptMax_;
    const double pt3Min_;
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
