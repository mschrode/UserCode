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
    int nCutValues() const { return cutVar_.at(0)->nCutValues(); }
    double cutValue(int cutVarIdx) const { return cutVar_.at(0)->cutValue(cutVarIdx); }
    double fittedValue(int parIdx, int cutVarIdx) const { return cutVar_.at(parIdx)->fittedValue(cutVarIdx); }
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
    double meanPt() const { return meanPt_; }
    double meanPtUncert() const { return meanPtUncert_; }
    double ptMin() const { return par_->ptMin(); }
    double ptMax() const { return par_->ptMax(); }
    TString ptMinStr() const { char str[10]; sprintf(str,"%.0f",ptMin()); return str; }
    TString ptMaxStr() const { char str[10]; sprintf(str,"%.0f",ptMax()); return str; }

    TH1 *getHistPtGen(const TString &newName) const { return getHist("hPtGen",newName); }
    TH1 *getHistPtGenJet1(const TString &newName) const { return getHist("hPtGenJet1",newName); }
    TH1 *getHistPdfPtTrue(const TString &newName) const { return getHist("hPdfPtTrue",newName); }
    TH1 *getHistResGen(const TString &newName) const { return getHist("hResGen",newName); }
    TH1 *getHistPdfRes(const TString &newName) const { return getHist("hPdfRes",newName); }
    TH1 *getHistPtAsym(const TString &newName) const { return getHist("hPtAsym",newName); }
    TH1 *getHistPdfPtAsym(const TString &newName) const { return getHist("hPdfPtAsym",newName); }
    TH1 *getHistMCRes(const TString &newName) const { return getHist("hMCRes",newName); }

    TF1 *getTF1OfVariation(int parIdx, const TString &name) const { return cutVar_.at(parIdx)->getTF1(name); }
    TGraphAsymmErrors *getTGraphOfVariation(int parIdx) const { return cutVar_.at(parIdx)->getTGraph(); }
    TH1 *getFrameOfVariation(int parIdx, const TString &name) const { return cutVar_.at(parIdx)->getFrame(name); }


  private:
    const Parameters::PtBinParameters *par_;

    double meanPt_;
    double meanPtUncert_;
    std::vector<double> extrapolatedVal_;
    std::vector<CutVariation*> cutVar_;
    std::vector<Uncertainty*> uncert_;

    TH1* hPtGen_;
    TH1* hPtGenJet1_;
    TH1* hPdfPtTrue_;
    TH1* hResGen_;
    TH1* hPdfRes_;
    TH1* hPtAsym_;
    TH1* hPdfPtAsym_;
    TH1 *hMCRes_;

    TH1 *getHist(const TString &name, const TString &newName) const;
  };
}
#endif
