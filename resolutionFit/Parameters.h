// $Id: Parameters.h,v 1.7 2010/07/20 14:08:38 mschrode Exp $

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <list>
#include <vector>

#include "TString.h"

#include "ResponseFunction.h"

class TRandom3;
namespace resolutionFit {

  //! \brief Store parameters like binning and file names centrally
  //!
  //! \author Matthias Schroeder
  //! \date 2010/05/15
  //! $Id: Parameters.h,v 1.7 2010/07/20 14:08:38 mschrode Exp $
  // --------------------------------------------
  class Parameters {
  public:
    class PtBinParameters {
    public:
      PtBinParameters(int ptBinIdx, const Parameters *par)
	: ptBinIdx_(ptBinIdx), par_(par) {};
      ~PtBinParameters() {};

      ResponseFunction::Type respFuncType() const { return par_->respFuncType(); }
      int nFittedPars() const { return par_->nFittedPars(); }
      bool isRelParValue(int parIdx) const { return par_->isRelParValue(parIdx); }
      TString parLabel(int parIdx) const { return par_->parLabel(parIdx); }
      TString parAxisLabel(int parIdx) const { return par_->parAxisLabel(parIdx); }
    
      int ptBinIdx() const { return ptBinIdx_; }
      double ptMin() const { return par_->ptMin(ptBinIdx()); }
      double ptMax() const { return par_->ptMax(ptBinIdx()); }
      double etaMin() const { return par_->etaMin(); }
      double etaMax() const { return par_->etaMax(); }
      TString fileNameStdSel() const { return par_->fileNameStdSel(ptBinIdx()); }
      
      bool hasMCStatUncert() const { return par_->hasMCStatUncert(); }
      TString fileNameMCStatUncert() const { return par_->fileNameMCStatUncert(ptBinIdx()); }
      
      int nPt3CutVariations() const { return par_->nPt3CutVariations(); }
      double pt3CutValue(int k) const { return par_->pt3CutValue(k); }
      TString fileNamePt3CutVariations(int k) const { return par_->fileNamePt3CutVariations(ptBinIdx(),k); }

      bool hasMCClosure() const { return par_->hasMCClosure(); }
      TString fileNameMCClosure() const { return par_->fileNameMCClosure(ptBinIdx()); }

      int nSystUncerts() const { return par_->nSystUncerts(); }
      TString labelSystUncert(int k) const { return par_->labelSystUncert(k); }
      TString fileNameSystUncertDown(int k) const { return par_->fileNameSystUncertDown(ptBinIdx(),k); }
      TString fileNameSystUncertUp(int k) const { return par_->fileNameSystUncertUp(ptBinIdx(),k); }

      int verbosity() const { return par_->verbosity(); }
      
    private:
      const int ptBinIdx_;
      const Parameters *par_;
    };

    
    Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, int startIdx, int endIdx, const TString &outNamePrefix, ResponseFunction::Type type, int verbosity);
    Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, const std::vector<int> fileNameIdx, const TString &outNamePrefix, ResponseFunction::Type type, int verbosity);

    ~Parameters();

    const PtBinParameters *createPtBinParameters(int ptBinIdx) const;

    ResponseFunction::Type respFuncType() const { return respFunc_->type(); }
    const ResponseFunction *respFunc() const { return respFunc_; }
    int nFittedPars() const { return respFunc_->nPars()-1; }
    bool isRelParValue(int parIdx) const;
    TString parLabel(int parIdx) const;
    TString parLabelTex(int parIdx) const;
    TString parAxisLabel(int parIdx) const;

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

    bool hasMCClosure() const { namesMCClosure_.size() > 0 ? true : false; }
    TString fileNameMCClosure(int ptBin) const { return namesMCClosure_.at(ptBin); }

    TString labelMeas() const;
    TString labelTruth() const;
    TString labelEtaBin() const;
    TString labelJetAlgo() const;
    TString labelLumi() const;
    TString labelPtBin(int ptBin, int type) const;
    TString labelPt3Cut(int ptBin) const;
    TString xAxisTitleResponse() const;
    TString yAxisTitleResponse() const;
    bool extendedLegend() const { return extendedLegend_; }


    double trueGaussResPar(int i) const { return trueResPar_.at(i); }
    bool hasMCTruthBins() const { return mcTruthPtBinEdges_.size() > 0 ? true : false; }
    int nMCTruthPtBins() const { return static_cast<int>(mcTruthPtBinEdges_.size())-1; }
    double mcTruthPtMin(int i) const { return mcTruthPtBinEdges_.at(i); }
    double mcTruthPtMax(int i) const { return mcTruthPtBinEdges_.at(i+1); }
    double mcTruthPtBinCenter(int i) const { return 0.5*(mcTruthPtMin(i)+mcTruthPtMax(i)); }
    double mcTruthPseudoMeas(int i) const { return mcTruthPseudoMeas_.at(i); }
    double mcTruthRelUncert() const { return mcTruthRelUncert_; }

    bool fitExtrapolatedSigma() const { return fitExtrapolatedSigma_; }
    bool fitRatio() const { return fitRatio_; }
    bool hasStartOffset() const { return startResOffset_ > 0 ? true : false; }
    double relStartOffset() const { return startResOffset_; }

    TString styleMode() const { return styleMode_; }
    int verbosity() const { return verbosity_; }


    void addPt3Cut(double pt3RelCutValue, const TString &fileBaseName);
    void addFileBaseNameMCStat(const TString &name) { writeFileNames(namesMCStat_,name); }
    void addFileBaseNameMCClosure(const TString &name) { writeFileNames(namesMCClosure_,name); }
    void addSystUncert(const TString &label, const TString &fileBaseNameDown, const TString &fileBaseNameUp);
    void addMCTruthBins(int nBins, double min, double max, double relUncert);
    void setTrueGaussResPar(double a0, double a1, double a2) {
      trueResPar_.at(0) = a0;
      trueResPar_.at(1) = a1;
      trueResPar_.at(2) = a2;
    }
    void addStartOffset(double delta) { startResOffset_ = delta; }
    void fitExtrapolatedSigma(bool fit) { fitExtrapolatedSigma_ = fit; }
    void fitRatio(bool fit) { fitRatio_ = fit; }

    void extendedLegend(bool ext) { extendedLegend_ = ext; }


  private:
    const double etaMin_;
    const double etaMax_;
    const TString outNamePrefix_;
    const TString styleMode_;
    const int verbosity_;

    std::vector<int> fileNameIdx_;
    std::vector<double> ptBinEdges_;
    std::vector<TString> namesStdSel_;
    std::vector<TString> namesMCStat_;
    std::vector<TString> namesMCClosure_;

    std::vector<double> pt3RelCutValues_;
    std::vector< std::vector<TString> > namesCutVars_;  // [ cut variation ] [ pt bin ]

    std::vector<TString> labelSyst_;
    std::vector< std::vector<TString> > namesSystDown_;  // [ uncertainty ] [ pt bin ]
    std::vector< std::vector<TString> > namesSystUp_;  // [ uncertainty ] [ pt bin ]

    std::vector<double> trueResPar_;
    std::vector<double> mcTruthPtBinEdges_;
    std::vector<double> mcTruthPseudoMeas_;
    double mcTruthRelUncert_;

    bool fitExtrapolatedSigma_;
    bool fitRatio_;
    double startResOffset_;

    bool extendedLegend_;

    ResponseFunction *respFunc_;
    TRandom3 *rand_;

    mutable std::list<PtBinParameters*> listOfPtBinParameters_;
    
    void init(const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, ResponseFunction::Type type);
    void writeFileNames(std::vector<TString> &names, const TString &baseName) const;
  };
}
#endif
