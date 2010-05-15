// $Id:  $

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <list>
#include <vector>

#include "TString.h"


namespace resolutionFit {

  //! \brief Store parameters like binning and file names centrally
  //!
  //! \author Matthias Schroeder
  //! \date 2010/05/15
  //! $Id:  $
  // --------------------------------------------
  class Parameters {
  public:
    
    class PtBinParameters {
    public:
      PtBinParameters(int idx, const Parameters *par)
	: idx_(idx), par_(par) {};
      ~PtBinParameters() {};

      int idx() const { return idx_; }
      double ptMin() const { return par_->ptMin(idx()); }
      double ptMax() const { return par_->ptMax(idx()); }
      double etaMin() const { return par_->etaMin(); }
      double etaMax() const { return par_->etaMax(); }
      TString fileNameStdSel() const { return par_->fileNameStdSel(idx()); }
      
      bool hasMCStatUncert() const { return par_->hasMCStatUncert(); }
      TString fileNameMCStatUncert() const { return par_->fileNameMCStatUncert(idx()); }
      
      int nPt3CutVariations() const { return par_->nPt3CutVariations(); }
      double pt3CutValue(int k) const { return par_->pt3CutValue(k); }
      TString fileNamePt3CutVariations(int k) const { return par_->fileNamePt3CutVariations(idx(),k); }

      int nSystUncerts() const { return par_->nSystUncerts(); }
      TString labelSystUncert(int k) const { return par_->labelSystUncert(k); }
      TString fileNameSystUncertDown(int k) const { return par_->fileNameSystUncertDown(idx(),k); }
      TString fileNameSystUncertUp(int k) const { return par_->fileNameSystUncertUp(idx(),k); }

      int verbosity() const { return par_->verbosity(); }
      
    private:
      const int idx_;
      const Parameters *par_;
    };

    
    Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, int startIdx, int endIdx, const TString &outNamePrefix, int verbosity);
    ~Parameters();

    const PtBinParameters *createPtBinParameters(int ptBinIdx) const;

    double etaMin() const { return etaMin_; }
    double etaMax() const { return etaMax_; }
    TString outNamePrefix() const { return outNamePrefix_; }
    int nPtBins() const { return static_cast<int>(ptBinEdges_.size())-1; }
    double ptMin(int ptBin) const { return ptBinEdges_.at(ptBin); }
    double ptMax(int ptBin) const { return ptBinEdges_.at(ptBin+1); }
    TString fileNameStdSel(int i) const { return namesStdSel_.at(i); }

    bool hasMCStatUncert() const { return namesMCStat_.size() > 0 ? true : false; }
    TString fileNameMCStatUncert(int i) const { return namesMCStat_.at(i); }

    int nPt3CutVariations() const { return static_cast<int>(pt3RelCutValues_.size()); }
    double pt3CutValue(int k) const { return pt3RelCutValues_.at(k); }
    TString fileNamePt3CutVariations(int i, int k) const { return namesCutVars_[k][i]; }

    int nSystUncerts() const { return static_cast<int>(labelSyst_.size()); }
    TString labelSystUncert(int k) const { return labelSyst_.at(k); }
    TString fileNameSystUncertDown(int i, int k) const { return namesSystDown_[k][i]; }
    TString fileNameSystUncertUp(int i, int k) const { return namesSystUp_[k][i]; }

    double trueGaussResPar(int i) const { return trueResPar_.at(i); }

    int verbosity() const { return verbosity_; }

    void addPt3Cut(double pt3RelCutValue, const TString &fileBaseName);
    void addFileBaseNameMCStat(const TString &name) { writeFileNames(namesMCStat_,name); }
    void addSystUncert(const TString &label, const TString &fileBaseNameDown, const TString &fileBaseNameUp);
    void setTrueGaussResPar(double a0, double a1, double a2) {
      trueResPar_.at(0) = a0;
      trueResPar_.at(1) = a1;
      trueResPar_.at(2) = a2;
    }


  private:
    const double etaMin_;
    const double etaMax_;
    const int startIdx_;
    const int endIdx_;
    const TString outNamePrefix_;
    const int verbosity_;

    std::vector<double> ptBinEdges_;
    std::vector<TString> namesStdSel_;
    std::vector<TString> namesMCStat_;

    std::vector<double> pt3RelCutValues_;
    std::vector< std::vector<TString> > namesCutVars_;  // [ cut variation ] [ pt bin ]

    std::vector<TString> labelSyst_;
    std::vector< std::vector<TString> > namesSystDown_;  // [ uncertainty ] [ pt bin ]
    std::vector< std::vector<TString> > namesSystUp_;  // [ uncertainty ] [ pt bin ]

    std::vector<double> trueResPar_;

    mutable std::list<PtBinParameters*> listOfPtBinParameters_;
    
    void writeFileNames(std::vector<TString> &names, const TString &baseName) const;
  };
}
#endif
