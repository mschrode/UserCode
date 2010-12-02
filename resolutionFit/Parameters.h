// $Id: Parameters.h,v 1.22 2010/11/11 12:57:04 mschrode Exp $

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <list>
#include <vector>

#include "TString.h"

#include "ResponseFunction.h"
#include "../util/utils.h"


class TRandom3;
namespace resolutionFit {

  enum FitMode { FitModeMaxLikeFull, FitModeMaxLikeSimple };
  enum BinPt { BinPtGen, BinPtAve, BinPtJet };
  enum Pt3Var { Pt3Rel, Pt3Abs, Pp3Rel };


  //! \brief Store parameters like binning and file names centrally
  //!
  //! \author Matthias Schroeder
  //! \date 2010/05/15
  //! $Id: Parameters.h,v 1.22 2010/11/11 12:57:04 mschrode Exp $
  // --------------------------------------------
  class Parameters {
  public:
    class PtBinParameters {
    public:
      PtBinParameters(int ptBinIdx, const Parameters *par)
	: ptBinIdx_(ptBinIdx), par_(par) {};
      ~PtBinParameters() {};

      FitMode fitMode() const { return par_->fitMode(); }
      ResponseFunction::Type respFuncType() const { return par_->respFuncType(); }
      int nFittedPars() const { return par_->nFittedPars(); }
      bool isRelParValue(int parIdx) const { return par_->isRelParValue(parIdx); }
      TString parLabel(int parIdx) const { return par_->parLabel(parIdx); }
      TString parAxisLabel(int parIdx) const { return par_->parAxisLabel(parIdx); }

      int ptBinIdx() const { return ptBinIdx_; }
      BinPt binPt() const { return par_->binPt(); }
      double ptMin() const { return par_->ptMin(ptBinIdx()); }
      double ptMax() const { return par_->ptMax(ptBinIdx()); }
      double etaMin() const { return par_->etaMin(); }
      double etaMax() const { return par_->etaMax(); }
      TString fileNameStdSel() const { return par_->fileNameStdSel(); }
      
      bool hasMCStatUncert() const { return par_->hasMCStatUncert(); }
      TString fileNameMCStatUncert() const { return par_->fileNameMCStatUncert(); }
      
      int nPt3Cuts() const { return par_->nPt3Cuts(); }
      bool pt3Bins() const { return par_->pt3Bins(); }
      Pt3Var pt3Var() const { return par_->pt3Var(); }
      double pt3Min(int k) const { return pt3Bins() ? par_->pt3Min(k) : 0.; }
      double pt3Max(int k) const { return par_->pt3Max(k); }
      double pt3Mean(int k) const { return pt3Bins() ? par_->pt3Mean(k) : 0.; }
      TString fileNamePt3CutVariations(int k) const { return par_->fileNamePt3CutVariations(k); }
      int stdSelIdx() const { return par_->stdSelIdx(); }

      bool hasMCClosure() const { return par_->hasMCClosure(); }
      TString fileNameMCClosure() const { return par_->fileNameMCClosure(); }

      bool fitPtGenAsym() const { return par_->fitPtGenAsym(); }
      bool hasCorrPtGenAsym() const { return par_->hasCorrPtGenAsym(); }

      int nSystUncerts() const { return par_->nSystUncerts(); }
      TString labelSystUncert(int k) const { return par_->labelSystUncert(k); }
      TString fileNameSystUncertDown(int k) const { return par_->fileNameSystUncertDown(k); }
      TString fileNameSystUncertUp(int k) const { return par_->fileNameSystUncertUp(k); }

      int verbosity() const { return par_->verbosity(); }
      
    private:
      const int ptBinIdx_;
      const Parameters *par_;
    };


    Parameters(double etaMin, double etaMax, double deltaPhi12, const TString &fileNameStdSel, const std::vector<double> &ptBinEdges, const TString &outNamePrefix, ResponseFunction::Type type, FitMode fitMode, BinPt binPt, int verbosity);
    ~Parameters();

    const PtBinParameters *createPtBinParameters(int ptBinIdx) const;


    bool isData() const { return isData_; }
    FitMode fitMode() const { return fitMode_; }
    ResponseFunction::Type respFuncType() const { return respFunc_->type(); }
    const ResponseFunction *respFunc() const { return respFunc_; }
    int nFittedPars() const { return respFunc_->nPars()-1; }
    bool isRelParValue(int parIdx) const;
    TString parLabel(int parIdx, bool maxLikeFit = true) const;
    TString parLabelTex(int parIdx, bool maxLikeFit = true) const;
    TString parAxisLabel(int parIdx) const;

    double etaMin() const { return etaMin_; }
    double etaMax() const { return etaMax_; }
    double deltaPhi12() const { return deltaPhi12_; }
    TString outNamePrefix() const { return outNamePrefix_; }
    BinPt binPt() const { return binPt_; }
    int nPtBins() const { return static_cast<int>(ptBinEdges_.size())-1; }
    double ptMin(int ptBin) const { return ptBinEdges_.at(ptBin); }
    double ptMax(int ptBin) const { return ptBinEdges_.at(ptBin+1); }
    const std::vector<double> *ptBinEdges() const { return &ptBinEdges_; }
    TString fileNameStdSel() const { return nameStdSel_; }

    bool hasMCStatUncert() const { return (nameMCStat_=="" ? false : true); }
    TString fileNameMCStatUncert() const { return nameMCStat_; }

    int nPt3Cuts() const { return static_cast<int>(pt3Max_.size()); }
    bool pt3Bins() const { return pt3Bins_; }
    Pt3Var pt3Var() const { return pt3Var_; }
    double pt3Min(int k) const { return pt3Min_.at(k); }
    double pt3Max(int k) const { return pt3Max_.at(k); }
    double pt3Mean(int k) const { return pt3Mean_.at(k); }
    TString fileNamePt3CutVariations(int k) const { return namesCutVars_[k]; }
    bool hasTruthSpectra() const { return hasTruthSpectra_; }
    TString fileNameTruthSpectrum(int varIdx) const { return namesTruthSpectra_[varIdx]; }
    int stdSelIdx() const { return stdSelIdx_; }

    int nSystUncerts() const { return static_cast<int>(labelSyst_.size()); }
    TString labelSystUncert(int k) const { return labelSyst_.at(k); }
    TString fileNameSystUncertDown(int k) const { return namesSystDown_[k]; }
    TString fileNameSystUncertUp(int k) const { return namesSystUp_[k]; }

    bool hasMCClosure() const { return (nameMCClosure_=="" ? false : false); }
    TString fileNameMCClosure() const { return nameMCClosure_; }

    TString labelMeas() const;
    TString labelTruth() const;
    TString labelEtaBin() const;
    TString labelJetAlgo() const;
    TString labelLumi() const;
    TString labelDeltaPhiCut() const { return "|#Delta#phi| < "+util::toTString(deltaPhi12()); }
    TString labelPtBin(int ptBin) const;
    TString labelPt3Cut() const { return labelPt3Cut(stdSelIdx()); }
    TString labelPt3Cut(int ptBin) const;
    TString labelPtSoftCut() const;
    TString labelPtGen() const;
    TString labelPtMeas() const;
    TString labelPtTrue() const;
    TString labelPtRef(const TString &plot = "") const;
    TString xAxisTitleResponse() const;
    TString yAxisTitleResponse() const;
    bool extendedLegend() const { return extendedLegend_; }
    double lumi() const { return lumi_; }

    double trueGaussResPar(int i) const { return trueResPar_.at(i); }
    double trueGaussSigma(double pt) const { 
      return sqrt( trueResPar_[0]*trueResPar_[0]/pt/pt + trueResPar_[1]*trueResPar_[1]/pt + trueResPar_[2]*trueResPar_[2] ); 
    }
    bool hasMCTruthBins() const { return mcTruthPtBinEdges_.size() > 0 ? true : false; }
    int nMCTruthPtBins() const { return static_cast<int>(mcTruthPtBinEdges_.size())-1; }
    double mcTruthPtMin(int i) const { return mcTruthPtBinEdges_.at(i); }
    double mcTruthPtMax(int i) const { return mcTruthPtBinEdges_.at(i+1); }
    double mcTruthPtBinCenter(int i) const { return 0.5*(mcTruthPtMin(i)+mcTruthPtMax(i)); }
    double mcTruthPseudoMeas(int i) const { return mcTruthPseudoMeas_.at(i); }
    double mcTruthRelUncert() const { return mcTruthRelUncert_; }

    bool fitPtGenAsym() const { return fitPtGenAsym_; }
    bool hasCorrPtGenAsym() const { return hasCorrPtGenAsym_; }
    double corrPtGenAsymPar(int i) const { return ptGenAsymPar_.at(i); }

    bool fitExtrapolatedSigma() const { return fitExtrapolatedSigma_; }
    bool fitRatio() const { return fitRatio_; }
    bool hasStartOffset() const { return startResOffset_ > 0 ? true : false; }
    double relStartOffset() const { return startResOffset_; }

    TString styleMode() const { return styleMode_; }
    int verbosity() const { return verbosity_; }


    void addPt3Threshold(Pt3Var pt3Variable, double pt3RelMax, const TString &fileName, const TString &spectrumName = "");
    void addPt3Bin(Pt3Var pt3Variable, double pt3RelMin, double pt3RelMax, double pt3RelMean, const TString &fileName, const TString &spectrumName = "");
    void addFileNameMCStat(const TString &name) { nameMCStat_ = name; }
    void addFileNameMCClosure(const TString &name) { nameMCClosure_ = name; }
    void addSystUncert(const TString &label, const TString &fileNameDown, const TString &fileNameUp);
    void addMCTruthBins(int nBins, double min, double max, double relUncert);
    void setTrueGaussResPar(double a0, double a1, double a2) {
      trueResPar_.at(0) = a0;
      trueResPar_.at(1) = a1;
      trueResPar_.at(2) = a2;
    }
    void addStartOffset(double delta) { startResOffset_ = delta; }
    void fitExtrapolatedSigma(bool fit) { fitExtrapolatedSigma_ = fit; }
    void fitRatio(bool fit) { fitRatio_ = fit; }
    void fitPtGenAsym(bool fit);
    void setParPtGenAsym(double a0, double a1, double a2);
    void setLumi(double lumi) { lumi_ = lumi; }
    void isData(bool data) { isData_ = data; }
    void extendedLegend(bool ext) { extendedLegend_ = ext; }


  private:
    const double etaMin_;
    const double etaMax_;
    const double deltaPhi12_;
    const FitMode fitMode_;
    const BinPt binPt_;
    const TString styleMode_;
    const int verbosity_;

    TString outNamePrefix_;
    bool isData_;
    double lumi_;

    std::vector<double> ptBinEdges_;
    const TString nameStdSel_;
    TString nameMCStat_;
    TString nameMCClosure_;

    bool pt3Bins_;
    int stdSelIdx_;
    Pt3Var pt3Var_;
    std::vector<double> pt3Min_;
    std::vector<double> pt3Max_;
    std::vector<double> pt3Mean_;
    std::vector<TString> namesCutVars_;  // [ cut variation ]
    bool hasTruthSpectra_;
    std::vector<TString> namesTruthSpectra_;

    std::vector<TString> labelSyst_;
    std::vector<TString> namesSystDown_;  // [ uncertainty ]
    std::vector<TString> namesSystUp_;  // [ uncertainty ]

    std::vector<double> trueResPar_;
    std::vector<double> mcTruthPtBinEdges_;
    std::vector<double> mcTruthPseudoMeas_;
    double mcTruthRelUncert_;

    bool fitPtGenAsym_;
    bool hasCorrPtGenAsym_;
    std::vector<double> ptGenAsymPar_;

    bool fitExtrapolatedSigma_;
    bool fitRatio_;
    double startResOffset_;

    bool extendedLegend_;

    ResponseFunction *respFunc_;
    TRandom3 *rand_;

    mutable std::list<PtBinParameters*> listOfPtBinParameters_;
    
    void init(const std::vector<double> &ptBinEdges, ResponseFunction::Type type);
    void print() const;
  };
}
#endif
