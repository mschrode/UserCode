#ifndef KALIBRI_FILE_PARSER_H
#define KALIBRI_FILE_PARSER_H

#include <cassert>
#include <vector>
#include <map>

#include "TString.h"

class TH1F;

namespace resolutionFit {
  class KalibriFileParser {
  public:
    KalibriFileParser(const TString &fileName, int verbose = 1);
    ~KalibriFileParser();

    int nValues() const { return static_cast<int>(values_.size()); }
    double value(int i = 0) const { assert( i>=0 && i < nValues() ); return values_[i]; }
    double statUncert(int i = 0) const { assert( i>=0 && i < nValues() ); return statUncert_[i]; }
    double meanPtGen() const { return meanPtGen_; }
    double meanPtDijet() const { return meanPtDijet_; }
    double meanPdfPtTrue() const { return meanPdfPtTrue_; }
    TH1F *hist(const TString &name, const TString &newName) const;

  private:
    typedef std::map<TString,TH1F*>::const_iterator HistIt;

    const int verbose_;

    std::vector<double> values_;
    std::vector<double> statUncert_;
    std::map<TString,TH1F*> hists_;

    double meanPtGen_;
    double meanPtDijet_;
    double meanPdfPtTrue_;

    int parse(const TString &fileName);
    void setMeanPt();
  };
}
#endif

