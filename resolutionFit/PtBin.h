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

    int ptBinIdx() const { return par_->idx(); }
    double relSigma() const { return relSigma_; }
    double uncertDown() const { return uncert_->down(); }
    double uncertUp() const { return uncert_->up(); }
    double uncertStatDown() const { return uncert_->down(0); }
    double uncertStatUp() const { return uncert_->up(0); }
    const Uncertainty *uncertStat() const { return uncert_->uncert(0); }
    double uncertSystDown() const { return uncert_->down(1); }
    double uncertSystUp() const { return uncert_->up(1); }
    const Uncertainty *uncertSyst() const { return uncert_->uncert(1); }
    int nUncertSyst() const { return uncert_->uncert(1)->nUncerts(); }
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

    TF1 *getTF1OfVariation(const TString &name) const { return cutVar_->getTF1(name); }
    TGraphAsymmErrors *getTGraphOfVariation() const { return cutVar_->getTGraph(); }
    TH1 *getFrameOfVariation(const TString &name) const { return cutVar_->getFrame(name); }


  private:
    const Parameters::PtBinParameters *par_;

    double meanPt_;
    double meanPtUncert_;
    double relSigma_;
    CutVariation *cutVar_;
    Uncertainty *uncert_;

    TH1* hPtGen_;
    TH1* hPtGenJet1_;
    TH1* hPdfPtTrue_;
    TH1* hResGen_;
    TH1* hPdfRes_;
    TH1* hPtAsym_;
    TH1* hPdfPtAsym_;

    TH1 *getHist(const TString &name, const TString &newName) const;
  };
}
#endif
