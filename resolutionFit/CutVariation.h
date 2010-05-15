#ifndef CUT_VARIATION_H
#define CUT_VARIATION_H

#include <cassert>
#include <vector>

#include "TString.h"

#include "Parameters.h"
#include "Uncertainty.h"

class TF1;
class TGraphAsymmErrors;
class TH1;

namespace resolutionFit {
  class CutVariation {
  public:
    CutVariation(const Parameters::PtBinParameters *par);
    ~CutVariation();

    int nCutValues() const { return par_->nPt3CutVariations(); }
    double minCutValue() const { return par_->pt3CutValue(0); }
    double maxCutValue() const { return par_->pt3CutValue(nCutValues()-1); }
    double cutValue(int i) const { return par_->pt3CutValue(i); }
    double relSigma(int i) const { return (*(varPoints_.at(i)))(); }
    double uncert(int i) const { return varPoints_.at(i)->uncertUp(); }
    double extrapolatedRelSigma() const { return (*extrapolatedPoint_)(); }
    double extrapolatedUncert() const { return extrapolatedPoint_->uncertUp(); }
    double meanPt() const { return meanPt_; }

    TF1 *getTF1(const TString &name) const;
    TGraphAsymmErrors *getTGraph() const;
    TH1 *getFrame(const TString &name) const;

    void extrapolate();

  private:
    class VariationPoint {
    public:
      VariationPoint();
      VariationPoint(double relSigma, const Uncertainty *uncert, double cutValue);
      ~VariationPoint();
      
      double operator()() const { return relSigma_; }
      double uncertDown() const { return uncert_->down(); }
      double uncertUp() const { return uncert_->up(); }
      double cutValue() const { return cutValue_; }
      
    private:
      const double relSigma_;
      const Uncertainty *uncert_;
      const double cutValue_;
    };

    const Parameters::PtBinParameters *par_;

    double meanPt_;    
    double mcStatUncert_;
    std::vector<VariationPoint*> varPoints_;
    VariationPoint* extrapolatedPoint_;
    TGraphAsymmErrors *graph_;
    TF1 *fit_;

    void createTGraph();
  };
}
#endif
