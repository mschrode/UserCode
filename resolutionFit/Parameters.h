// $ Id: $

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "TString.h"

#include "BinningAdmin.h"
#include "../util/utils.h"
#include "OutputManager.h"

namespace resolutionFit {

  class JetProperties {
  public:
    enum Algo { AK5 };
    enum Type { Calo, JPT, PF };
  };


  // Store parameters
  //
  // Interface to eta,pt,ptSoft binning information,
  // file name conventions, jet types, labels
  //
  // Make Singleton
  // -------------------------------------------------------------------------------------
  class Parameters {
  public:
    Parameters(const TString &id, const TString &binAdmCfg, unsigned int verbosity = 0);
    ~Parameters();

    void setJetProperties(JetProperties::Algo algo, JetProperties::Type type);

    unsigned int verbosity() const { return verbosity_; }

    JetProperties::Algo jetAlgo() const { return jetAlgo_; }
    JetProperties::Type jetType() const { return jetType_; }

    unsigned int nEtaBins() const { return binAdm_->nEtaBins(); }
    unsigned int nPtBins(unsigned int etaBin) const { return binAdm_->nPtBins(etaBin); }
    unsigned int nPtSoftBins() const { return binAdm_->nPtSoftBins(); }
    double etaMin(unsigned int etaBin) const { return binAdm_->etaMin(etaBin); }
    double etaMax(unsigned int etaBin) const { return binAdm_->etaMax(etaBin); }
    double ptMin(unsigned int etaBin, unsigned int ptBin) const { return binAdm_->ptMin(etaBin,ptBin); }
    double ptMax(unsigned int etaBin, unsigned int ptBin) const { return binAdm_->ptMax(etaBin,ptBin); }
    double ptMin(unsigned int etaBin) const { return binAdm_->ptMin(etaBin); }
    double ptMax(unsigned int etaBin) const { return binAdm_->ptMax(etaBin); }
    const std::vector<double> ptBinEdges(unsigned int etaBin) const { return binAdm_->ptBinEdges(etaBin); }
    const std::vector<double>& ptSoft() const { return binAdm_->ptSoftMax(); }
    double ptSoftMin(unsigned int ptSoftBin) const { return binAdm_->ptSoftMin(ptSoftBin); }
    double ptSoftMax(unsigned int ptSoftBin) const { return binAdm_->ptSoftMax(ptSoftBin); }

    void printBinning() const { binAdm_->printBinning(); }

    TString histNameSuffix(unsigned int etaBin, unsigned int ptBin) const {
      return etaHistNameTag(etaBin)+ptHistNameTag(ptBin);
    }
    TString etaHistNameTag(unsigned int etaBin) const { return "_Eta"+util::toTString(etaBin); }
    TString ptHistNameTag(unsigned int ptBin) const { return "_Pt"+util::toTString(ptBin); }

    TString fileNameSuffix(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) const {
      return etaFileNameTag(etaBin)+ptSoftFileNameTag(ptSoftBin)+fileNameEnding();
    }
    TString etaFileNameTag(unsigned int etaBin) const { return "_Eta"+util::toTString(etaBin); }
    TString ptFileNameTag(unsigned int ptBin) const { return "_Pt"+util::toTString(ptBin); }
    TString ptSoftFileNameTag(unsigned int ptSoftBin) const { return "_PtSoft"+util::toTString(ptSoftBin+1); }
    TString fileNameEnding() const { return ".root"; }

    OutputManager::Mode outMode() const { return OutputManager::PSAllInOne; }
    TString outFilePrefix() const { return id_; }
      

  private:
    const unsigned int verbosity_;
    const TString id_;

    JetProperties::Algo jetAlgo_;
    JetProperties::Type jetType_;

    sampleTools::BinningAdmin* binAdm_;
  };
}
#endif
