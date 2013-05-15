// $Id: Parameters.h,v 1.1 2013/05/08 13:23:32 mschrode Exp $

#ifndef SAMPLETOOLS_PARAMETERS_H
#define SAMPLETOOLS_PARAMETERS_H

#include "TString.h"

#include "WeightProducer.h"


namespace sampleTools {

  class Parameters {
  public:
    Parameters();

    // Input
    TString inFileListData() const { return inFileListData_; }
    TString inFileListMC() const { return inFileListMC_; }
    bool runOnData() const { return runOnData_; }
    sampleTools::WeightProducer::QCDSample qcdSampleType() const { return qcdSampleType_; }
    sampleTools::WeightProducer::PUScenario puScenario() const { return puScenario_; }
    TString puHistFile() const { return puHistFile_; }
    int nEvts() const { return nEvts_; }

    // Output
    TString histFilePath() const { return histFilePath_; }
    TString histFilePrefix() const { return histFilePrefix_; }
    TString histFileData() const { return histFileData_; }
    TString histFileMC() const { return histFileMC_; }

    // Binning (eta,pt,pt3rel)
    TString configAdm() const { return configAdm_; }
    TString configBin() const { return configBin_; }
    TString configBinMerged() const { return configBinMerged_; }


    // Set values
    void setInFileListData(const TString &val) { inFileListData_ = val; }
    void setInFileListMC(const TString &val) { inFileListMC_ = val; }
    void setQCDSampleType(sampleTools::WeightProducer::QCDSample val) { qcdSampleType_ = val; }
    void setPUScenario(sampleTools::WeightProducer::PUScenario val) { puScenario_ = val; }
    void setPUHistFile(const TString &val) { puHistFile_ = val; }
    void setHistFilePath(const TString &val) { histFilePath_ = val; }
    void setHistFilePrefix(const TString &val) { histFilePrefix_ = val; }
    void setHistFileData(const TString &val) { histFileData_ = val; }
    void setHistFileMC(const TString &val) { histFileMC_ = val; }
    void setConfigAdm(const TString &val) { configAdm_ = val; }
    void setConfigBin(const TString &val) { configBin_ = val; }
    void setConfigBinMerged(const TString &val) { configBinMerged_ = val; }


  private:
    TString inFileListData_;
    TString inFileListMC_;
    bool runOnData_;
    sampleTools::WeightProducer::QCDSample qcdSampleType_;
    sampleTools::WeightProducer::PUScenario puScenario_;
    TString puHistFile_;
    int nEvts_;

    TString histFilePath_;
    TString histFilePrefix_;
    TString histFileData_;
    TString histFileMC_;

    TString configAdm_;
    TString configBin_;
    TString configBinMerged_;
  };


  Parameters::Parameters() {
    inFileListData_ = "";
    inFileListMC_ = "";
    runOnData_ = "";
    qcdSampleType_ = sampleTools::WeightProducer::Flat;
    puScenario_ = sampleTools::WeightProducer::Summer12S10;
    puHistFile_ = "";
    nEvts_ = 0;
    histFilePath_ = "";
    histFilePrefix_ = "";
    histFileData_ = "";
    histFileMC_ = "";
    configAdm_ = "";
    configBin_ = "";
    configBinMerged_ = "";
  }
}
#endif
