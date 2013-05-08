// $Id: $

#ifndef SAMPLETOOLS_PARAMETERS_H
#define SAMPLETOOLS_PARAMETERS_H

#include "TString.h"

#include "WeightProducer.h"


namespace sampleTools {

class Parameters {
public:
  // Input
  static TString inFileListData;
  static TString inFileListMC;
  static bool runOnData;
  static sampleTools::WeightProducer::QCDSample qcdSampleType;
  static sampleTools::WeightProducer::PUScenario puScenario;
  static TString puHistFile;
  static int nEvts;

  // Output
  static TString histFilePath;
  static TString histFilePrefix;
  static TString histFileData;
  static TString histFileMC;

  // Binning (eta,pt,pt3rel)
  static TString configAdm;
  static TString configBin;
  static TString configBinMerged;

  // Histogram names
  static TString hNamePtAve;
  static TString hNamePtAsym;
  static TString hNamePtAsymAbs;
  static TString hNameResp;


  // Checks
  static bool isInit;
};


bool Parameters::isInit = false;
TString Parameters::inFileListData = "";
TString Parameters::inFileListMC = "";
bool Parameters::runOnData = "";
sampleTools::WeightProducer::QCDSample Parameters::qcdSampleType = sampleTools::WeightProducer::Flat;
sampleTools::WeightProducer::PUScenario Parameters::puScenario = sampleTools::WeightProducer::Summer12S10;
TString Parameters::puHistFile = "";
int Parameters::nEvts = 0;
TString Parameters::histFilePath = "";
TString Parameters::histFilePrefix = "";
TString Parameters::histFileData = "";
TString Parameters::histFileMC = "";
TString Parameters::configAdm = "";
TString Parameters::configBin = "";
TString Parameters::configBinMerged = "";

TString Parameters::hNamePtAve = "";
TString Parameters::hNamePtAsym = "";
TString Parameters::hNamePtAsymAbs = "";
TString Parameters::hNameResp = "";

}
#endif
