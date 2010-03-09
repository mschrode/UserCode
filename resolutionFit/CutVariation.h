#ifndef CUT_VARIATION_H
#define CUT_VARIATION_H

#include <cassert>
#include <vector>

#include "TString.h"

#include "Uncertainty.h"

class TF1;
class TGraphAsymmErrors;
class TH1D;

namespace resolutionFit {
  class CutVariation {
  public:
    CutVariation(const std::vector<TString> &fileNames, const std::vector<double> &cutValues, int verbose = 1);
    ~CutVariation();

    int nCutValues() const { return static_cast<int>(varPoints_.size()); }
    double minCutValue() const { return varPoints_.front()->cutValue(); }
    double maxCutValue() const { return varPoints_.back()->cutValue(); }
    double cutValue(int i) const { assert( i>=0 && i<nCutValues() ); return varPoints_[i]->cutValue(); }
    double relSigma(int i) const { assert( i>=0 && i<nCutValues() ); return (*varPoints_[i])(); }
    double uncert(int i) const { assert( i>=0 && i<nCutValues() ); return varPoints_[i]->uncertUp(); }
    double extrapolatedRelSigma() const { return (*extrapolatedPoint_)(); }
    double extrapolatedUncert() const { return extrapolatedPoint_->uncertUp(); }
    double meanPt() const { return meanPt_; }

    TF1 *getTF1(const TString &name) const;
    TGraphAsymmErrors *getTGraph() const;
    TH1D *getFrame(const TString &name) const;

    void extrapolate();

  private:
    static int nObjs;

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

    const int verbose_;

    double meanPt_;    
    std::vector<VariationPoint*> varPoints_;
    VariationPoint* extrapolatedPoint_;
    TGraphAsymmErrors *graph_;
    TF1 *fit_;

    void createTGraph();
  };
}
#endif
