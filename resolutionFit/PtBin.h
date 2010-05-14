#ifndef PT_BIN_H
#define PT_BIN_H

#include <vector>

#include "TString.h"

#include "Uncertainty.h"
#include "CutVariation.h"

class TF1;
class TGraphAsymmErrors;
class TH1;

namespace resolutionFit {
  class PtBin {
  public:
    PtBin(const TString &fileNameStdSel,
	  const std::vector<TString> &fileNamesCutVariation, const std::vector<double> &cutValues,
	  const TString &fileNameMCStatUncert,
	  const std::vector<TString> &fileNamesSystUncertUp,
	  const std::vector<TString> &fileNamesSystUncertDown,
	  const std::vector<TString> &labelsSystUncertainties, double minPt, double maxPt, int verbose = 1);
    ~PtBin();

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
    double minPt() const { return minPt_; }
    double maxPt() const { return maxPt_; }
    TString minPtStr() const { char str[10]; sprintf(str,"%.0f",minPt()); return str; }
    TString maxPtStr() const { char str[10]; sprintf(str,"%.0f",maxPt()); return str; }

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
    static int nPtBins;

    const int verbose_;

    double minPt_;
    double maxPt_;
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
