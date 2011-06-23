// $Id: Parameters.h,v 1.28 2011/06/07 18:23:30 mschrode Exp $

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

    static bool isValidType(Type type);
    static TString toString(Type type);
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
    void setOutMode(OutputManager::Mode mode);
    void setNEtaBinsUser(unsigned int nBins);

    unsigned int verbosity() const { return verbosity_; }

    JetProperties::Algo jetAlgo() const { return jetAlgo_; }
    JetProperties::Type jetType() const { return jetType_; }
    double lumi() const { return 295.; }

    unsigned int nEtaBins() const { return nEtaBinsUser_>0 ? nEtaBinsUser_ : binAdm_->nEtaBins(); }
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

    OutputManager::Mode outMode() const { return outPutMode_; }
    TString outFilePrefix() const { return id_; }
      

  private:
    const unsigned int verbosity_;

    TString id_;
    unsigned int nEtaBinsUser_;

    JetProperties::Algo jetAlgo_;
    JetProperties::Type jetType_;
    OutputManager::Mode outPutMode_;

    sampleTools::BinningAdmin* binAdm_;
  };
}
#endif
