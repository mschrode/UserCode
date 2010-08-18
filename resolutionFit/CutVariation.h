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
    CutVariation(const Parameters::PtBinParameters *par, int parIndex, bool maxLikeFit = true);
    ~CutVariation();

    int nPt3Cuts() const { return par_->nPt3Cuts(); }
    bool pt3Bins() const { return par_->pt3Bins(); }
    double pt3Min(int i) const { return par_->pt3Min(i); }
    double pt3Max(int i) const { return par_->pt3Max(i); }
    double pt3Mean(int i) const { return par_->pt3Mean(i); }

    int parIdx() const { return parIdx_; }
    bool isRelValue() const { return par_->isRelParValue(parIdx()); }
    double fittedValue(int i) const { return (*(varPoints_.at(i)))(); }
    double uncert(int i) const { return varPoints_.at(i)->uncertUp(); }
    double extrapolatedValue() const { return (*extrapolatedPoint_)(); }
    double extrapolatedUncert() const { return extrapolatedPoint_->uncertUp(); }
    double meanPt() const { return meanPt_; }
    double meanPtUncert() const { return meanPtUncert_; }

    TF1 *getTF1(const TString &name) const;
    TGraphAsymmErrors *getTGraph() const;
    TH1 *getFrame(const TString &name) const;
    TH1 *getPtAsymmetry(int i, const TString &name) const {
      return varPoints_.at(i)->histPtAsym(name);
    }

    void extrapolate();

  private:
    class VariationPoint {
    public:
      VariationPoint();
      VariationPoint(double fitValue, Uncertainty *uncert, double cutValue, TH1 *hPtAsym);
      VariationPoint(TH1 *hPtAsym, double cutValue);
      ~VariationPoint();
      
      double operator()() const { return fitValue_; }
      double uncertDown() const { return uncert_->down(); }
      double uncertUp() const { return uncert_->up(); }
      double cutValue() const { return cutValue_; }
      TH1 *histPtAsym(const TString &name) const;
      
    private:
      double fitValue_;
      Uncertainty *uncert_;
      const double cutValue_;
      TH1 *hPtAsym_;
    };

    const Parameters::PtBinParameters *par_;
    const int parIdx_;

    double meanPt_;    
    double meanPtUncert_;
    double mcStatUncert_;
    std::vector<VariationPoint*> varPoints_;
    VariationPoint* extrapolatedPoint_;
    TGraphAsymmErrors *graph_;
    TF1 *fit_;

    void createTGraph();
  };
}
#endif
