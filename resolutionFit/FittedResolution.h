#ifndef FITTED_RESOLUTION_H
#define FITTED_RESOLUTION_H

#include <cassert>
#include <vector>

#include "TString.h"

#include "Parameters.h"
#include "PtBin.h"

class TF1;
class TGraphAsymmErrors;

namespace resolutionFit {
  class FittedResolution {
  public:
    FittedResolution(const std::vector<PtBin*> &ptBins, const Parameters *par);
    ~FittedResolution();

    int nPtBins() const { return static_cast<int>(ptBins_.size()); }
    double ptMin() const { return ptBins_.front()->ptMin(); }
    double ptMax() const { return ptBins_.back()->ptMax(); }
    double ptMin(int ptBin) const { return par_->ptMin(ptBin); }
    double ptMax(int ptBin) const { return par_->ptMax(ptBin); }
    double meanPt(int ptBin) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->meanPt();
    }
    double extrapolatedValue(int ptBin, int parIdx) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->extrapolatedValue(parIdx);
    }
    double uncertStat(int ptBin, int parIdx) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->uncertStatDown(parIdx);
    }
    int nUncertSyst(int parIdx) const { return ptBins_[0]->nUncertSyst(parIdx); }
    double uncertSystDown(int ptBin, int parIdx) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->uncertSystDown(parIdx);
    }
    double uncertSystUp(int ptBin, int parIdx) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->uncertSystUp(parIdx);
    }


    void plotExtrapolation() const;
    void plotPtAsymmetryBins() const;
    void plotResolution() const;
    void plotResolutionBins() const;
    void plotSpectra() const;
    void plotSystematicUncertainties() const;
    void print() const;

  private:
    const Parameters *par_;
    const std::vector<PtBin*> ptBins_;

    double ptMin_;
    double ptMax_;

    TF1 *trueRes_;
    TF1 *fittedRes_;

    TGraphAsymmErrors *getTGraphOfResolution(const TString &uncertType) const;
  };
}
#endif
