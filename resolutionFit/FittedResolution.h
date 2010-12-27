#ifndef FITTED_RESOLUTION_H
#define FITTED_RESOLUTION_H

#include <cassert>
#include <vector>

#include "TGraphAsymmErrors.h"
#include "TString.h"

#include "Parameters.h"
#include "PtBin.h"

class TF1;

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
    double meanPtUncert(int ptBin) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->meanPtUncert();
    }
    double meanPtAve(int ptBin) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->meanPtAve();
    }
    double meanPtAsym(int ptBin) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->meanPtAsym();
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
    double fittedValue(int ptBin, int parIdx, int cutIdx) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->fittedValue(parIdx,cutIdx);
    }

    double extrapolatedAsym(int ptBin) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->extrapolatedAsym();
    }
    double uncertStatAsym(int ptBin) const {
      assert( ptBin>=0 && ptBin<nPtBins() );
      return ptBins_[ptBin]->uncertDownAsym();
    }

    int nMCClosureResFits() const { return static_cast<int>(mcClosureResoLabels_.size()); }
    TString mcClosureLabel(int i) const { return mcClosureResoLabels_.at(i); }
    double mcClosureReso(int i, int ptBin) const { return mcClosureGReso_.at(i)->GetY()[ptBin]; }
    double mcClosureResoErr(int i, int ptBin) const { return mcClosureGReso_.at(i)->GetEYhigh()[ptBin]; }

    void plotExtrapolation() const;
    void plotPtAsymmetry() const;
    void plotResolution() const;
    void plotResolutionDistributions() const;
    void plotSpectra() const;
    void plotAdditionalJetActivity() const;
    void plotControlDistributions() const;
    void plotMCClosure() const;
    void print() const;
    void printLaTeX() const;
    void printPoints() const;
    void createSlides() const;
    void writeRootOutput() const;
	

  private:
    const Parameters *par_;
    const std::vector<PtBin*> ptBins_;

    double ptMin_;
    double ptMax_;

    TF1 *trueRes_;
    TF1 *fittedRes_;
    TF1 *fittedResAsym_;
    TF1 *ptGenAsym_;

    std::vector<TString> mcClosureResoLabels_;
    std::vector< std::vector<TF1*> > mcClosureFits_;
    std::vector<TGraphAsymmErrors*> mcClosureGReso_;
    std::vector<TGraphAsymmErrors*> mcClosureGScale_;

    int lineWidth_;
    double lineHeight_;

    static double gaussian(double *x, double *par);
    void fitMCClosure();
    TGraphAsymmErrors *getTGraphOfResolution(const TString &method, const TString &uncertainty, bool corrected = false) const;
  };
}
#endif
