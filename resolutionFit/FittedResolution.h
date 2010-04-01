#ifndef FITTED_RESOLUTION_H
#define FITTED_RESOLUTION_H

#include <cassert>
#include <vector>

#include "TString.h"

#include "PtBin.h"

class TF1;
class TGraphAsymmErrors;

namespace resolutionFit {
  class FittedResolution {
  public:
    FittedResolution(const std::vector<PtBin*> &ptBins, const std::vector<double> &trueResPar, const TString &outNamePrefix = "");
    ~FittedResolution();

    int nPtBins() const { return static_cast<int>(ptBins_.size()); }
    double minPt() const { return ptBins_.front()->minPt(); }
    double maxPt() const { return ptBins_.back()->maxPt(); }
    double minPt(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->minPt(); }
    double maxPt(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->maxPt(); }
    double meanPt(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->meanPt(); }
    double relSigma(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->relSigma(); }
    double uncertStat(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->uncertStatDown(); }
    int nUncertSyst() const { return ptBins_[0]->nUncertSyst(); }
    double uncertSystDown(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->uncertSystDown(); }
    double uncertSystUp(int i) const { assert( i>=0 && i < nPtBins() ); return ptBins_[i]->uncertSystUp(); }


    void plotExtrapolation() const;
    void plotResolution() const;
    void plotResolutionBins() const;
    void plotSpectra() const;
    void plotSystematicUncertainties() const;
    void print() const;

  private:
    const bool verbose_;
    const TString outNamePrefix_;

    std::vector<PtBin*> ptBins_;
    TF1 *trueRes_;
    TF1 *fittedRes_;

    TGraphAsymmErrors *getTGraphOfResolution(const TString &uncertType) const;
  };
}
#endif
