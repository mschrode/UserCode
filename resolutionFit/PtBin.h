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
    double extrapolatedValue(int parIdx, bool corrected = false) const;
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
    double extrapolatedAsym(bool corrected = false) const;
    double uncertDownAsym() const { return uncertAsym_->down(); }
    double uncertUpAsym() const { return uncertAsym_->up(); }

    double fittedGenAsym(int cutVarIdx) const {
      return cutVarGenAsym_->fittedValue(cutVarIdx);
    }
    double fittedGenAsymUncert(int cutVarIdx) const {
      return cutVarGenAsym_->uncert(cutVarIdx);
    }
    double extrapolatedGenAsym() const { return extrapolatedGenAsym_; }
    double uncertDownGenAsym() const { return uncertGenAsym_->down(); }
    double uncertUpGenAsym() const { return uncertGenAsym_->up(); }

    TH1 *getHistPtJet1(const TString &newName) const { return getHist("hPtJet1",newName); }
    TH1 *getHistPtJet2(const TString &newName) const { return getHist("hPtJet2",newName); }
    TH1 *getHistPtJet3(const TString &newName) const { return getHist("hPtJet3",newName); }
    TH1 *getHistPtJet4(const TString &newName) const { return getHist("hPtJet4",newName); }
    TH1 *getHistPtGen(const TString &newName) const { return getHist("hPtGen",newName); }
    TH1 *getHistPtGenJet1(const TString &newName) const { return getHist("hPtGenJet1",newName); }
    TH1 *getHistPtAve(const TString &newName) const { return getHist("hPtAve",newName); }
    TH1 *getHistPtAveAbs(const TString &newName) const { return getHist("hPtAve",newName); }
    TH1 *getHistPdfPtTrue(const TString &newName) const { return getHist("hPdfPtTrue",newName); }
    TH1 *getHistMCRes(const TString &newName) const { return getHistMCRes(par_->stdSelIdx(),newName); }
    TH1 *getHistMCRes(int cutVarIdx, const TString &newName) const { 
      return cutVar_[0]->getMCRes(cutVarIdx,newName); 
    }
    TH1 *getHistPdfRes(const TString &newName) const { return getHist("hPdfRes",newName); }
    TH1 *getHistPtAsym(const TString &newName) const { return getHistPtAsym(par_->stdSelIdx(),newName); }
    TH1 *getHistPtAsym(int cutVarIdx, const TString &newName) const {
      return cutVarAsym_->getPtAsymmetry(cutVarIdx,newName);
    }
    TH1 *getHistPtGenAsym(const TString &newName) const { return getHist("hPtGenAsym",newName); }
    TH1 *getHistPJet3(const TString &newName) const { return getHist("hPJet",newName); }
    TH1 *getHistPJet3Rel(const TString &newName) const { return getHist("hPJet3Rel",newName); }
    TH1 *getHistPJet3GenRel(const TString &newName) const { return getHist("hPJet3GenRel",newName); }
    TH1 *getHistPSJ(const TString &newName) const { return getHist("hPSJ",newName); }
    TH1 *getHistPSJRel(const TString &newName) const { return getHist("hPSJRel",newName); }
    TH1 *getHistPSJGenRel(const TString &newName) const { return getHist("hPSJGenRel",newName); }
    TH1 *getHistEta(const TString &newName) const { return getHist("hEta",newName); }
    TH1 *getHistDeltaPhi12(const TString &newName) const { return getHist("hDeltaPhi12",newName); }


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

    double extrapolatedGenAsym_;
    CutVariation* cutVarGenAsym_;
    Uncertainty* uncertGenAsym_;

    TH1* hPtJet1_;
    TH1* hPtJet2_;
    TH1* hPtJet3_;
    TH1* hPtJet4_;
    TH1* hPtGen_;
    TH1* hPtGenJet1_;
    TH1* hPdfPtTrue_;
    TH1 *hPtAve_;
    TH1* hResGen_;
    TH1* hPdfRes_;
    TH1* hPtGenAsym_;
    TH1 *hPJet3_;
    TH1 *hPJet3Rel_;
    TH1 *hPJet3GenRel_;
    TH1 *hPSJ_;
    TH1 *hPSJRel_;
    TH1 *hPSJGenRel_;
    TH1 *hEta_;
    TH1 *hDeltaPhi12_;


    TH1 *getHist(const TString &name, const TString &newName) const;
  };
}
#endif
