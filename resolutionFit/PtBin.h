#ifndef PT_BIN_H
#define PT_BIN_H

#include <vector>

#include "TString.h"

#include "CutVariation.h"
#include "Parameters.h"
#include "Uncertainty.h"

class TF1;
class TGraphAsymmErrors;
class TH1;

namespace resolutionFit {
  class PtBin {
  public:
    PtBin(const Parameters::PtBinParameters *par);
    ~PtBin();

    int ptBinIdx() const { return par_->ptBinIdx(); }

    double meanPt() const { return cutVar_.at(0)->meanPt(); }
    double meanPtUncert() const { return cutVar_.at(0)->meanPtUncert(); }
    double meanPtAve() const { return cutVarAsym_->meanPt(); }
    double meanPtAveUncert() const { return cutVarAsym_->meanPtUncert(); }
    double meanPtAsym() const { return meanPtAve(); }
    double meanPtAsymUncert() const { return meanPtAveUncert(); }
    double ptMin() const { return par_->ptMin(); }
    double ptMax() const { return par_->ptMax(); }
    TString ptMinStr() const { char str[10]; sprintf(str,"%.0f",ptMin()); return str; }
    TString ptMaxStr() const { char str[10]; sprintf(str,"%.0f",ptMax()); return str; }

    int nPt3Cuts() const { return cutVar_.at(0)->nPt3Cuts(); }
    bool pt3Bins() const { return cutVar_.at(0)->pt3Bins(); }
    double pt3Min(int cutVarIdx) const { return cutVar_.at(0)->pt3Min(cutVarIdx); }
    double pt3Max(int cutVarIdx) const { return cutVar_.at(0)->pt3Max(cutVarIdx); }
    double pt3Mean(int cutVarIdx) const { return cutVar_.at(0)->pt3Mean(cutVarIdx); }

    double fittedValue(int parIdx, int cutVarIdx) const { return cutVar_.at(parIdx)->fittedValue(cutVarIdx); }
    double fittedValueUncert(int parIdx, int cutVarIdx) const { return cutVar_.at(parIdx)->uncert(cutVarIdx); }
    double extrapolatedValue(int parIdx) const { return extrapolatedVal_.at(parIdx); }
    double uncertDown(int parIdx) const { return uncert_.at(parIdx)->down(); }
    double uncertUp(int parIdx) const { return uncert_.at(parIdx)->up(); }
    double uncertStatDown(int parIdx) const { return uncert_.at(parIdx)->down(0); }
    double uncertStatUp(int parIdx) const { return uncert_.at(parIdx)->up(0); }
    const Uncertainty *uncertStat(int parIdx) const { return uncert_.at(parIdx)->uncert(0); }
    double uncertSystDown(int parIdx) const { return uncert_.at(parIdx)->down(1); }
    double uncertSystUp(int parIdx) const { return uncert_.at(parIdx)->up(1); }
    const Uncertainty *uncertSyst(int parIdx) const { return uncert_.at(parIdx)->uncert(1); }
    int nUncertSyst(int parIdx) const { return uncert_.at(parIdx)->uncert(1)->nUncerts(); }

    double fittedAsym(int cutVarIdx) const {
      return cutVarAsym_->fittedValue(cutVarIdx);
    }
    double fittedAsymUncert(int cutVarIdx) const {
      return cutVarAsym_->uncert(cutVarIdx);
    }
    double extrapolatedAsym() const { return extrapolatedAsym_; }
    double uncertDownAsym() const { return uncertAsym_->down(); }
    double uncertUpAsym() const { return uncertAsym_->up(); }

    TH1 *getHistPtGen(const TString &newName) const { return getHist("hPtGen",newName); }
    TH1 *getHistPtGenJet1(const TString &newName) const { return getHist("hPtGenJet1",newName); }
    TH1 *getHistPtAve(const TString &newName) const { return getHist("hPtAve",newName); }
    TH1 *getHistPtAveAbs(const TString &newName) const { return getHist("hPtAveAbs",newName); }
    TH1 *getHistPdfPtTrue(const TString &newName) const { return getHist("hPdfPtTrue",newName); }
    TH1 *getHistResGen(const TString &newName) const { return getHist("hResGen",newName); }
    TH1 *getHistPdfRes(const TString &newName) const { return getHist("hPdfRes",newName); }
    TH1 *getHistPtAsym(const TString &newName) const { return getHistPtAsym(par_->stdSelIdx(),newName); }
    TH1 *getHistPtAsym(int cutVarIdx, const TString &newName) const {
      return cutVarAsym_->getPtAsymmetry(cutVarIdx,newName);
    }
    TH1 *getHistPtGenAsym(const TString &newName) const { return getHist("hPtGenAsym",newName); }
    TH1 *getHistMCRes(const TString &newName) const { return getHist("hMCRes",newName); }

    TF1 *getTF1OfVariation(int parIdx, const TString &name) const { return cutVar_.at(parIdx)->getTF1(name); }
    TF1 *getTF1OfVariationAsym(const TString &name) const { return cutVarAsym_->getTF1(name); }
    TGraphAsymmErrors *getTGraphOfVariation(int parIdx) const { return cutVar_.at(parIdx)->getTGraph(); }
    TGraphAsymmErrors *getTGraphOfVariationAsym() const { return cutVarAsym_->getTGraph(); }
    TH1 *getFrameOfVariation(int parIdx, const TString &name) const { return cutVar_.at(parIdx)->getFrame(name); }


  private:
    const Parameters::PtBinParameters *par_;

    std::vector<double> extrapolatedVal_;
    std::vector<CutVariation*> cutVar_;
    std::vector<Uncertainty*> uncert_;

    double extrapolatedAsym_;
    CutVariation* cutVarAsym_;
    Uncertainty* uncertAsym_;


    TH1* hPtGen_;
    TH1* hPtGenJet1_;
    TH1* hPdfPtTrue_;
    TH1 *hPtAve_;
    TH1 *hPtAveAbs_;
    TH1* hResGen_;
    TH1* hPdfRes_;
    TH1* hPtGenAsym_;
    TH1 *hMCRes_;

    TH1 *getHist(const TString &name, const TString &newName) const;
  };
}
#endif
