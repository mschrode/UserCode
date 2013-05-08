// $Id: $

#include "WeightProducer.h"
#include "Parameters.h"

#include "fillResTailInputHists.C"
#include "mergePtBins.C"


TString fileNameId(const TString &inFileList, bool isData) {
  TString fileNameId = "";

  if( inFileList.Contains("Calo") ) fileNameId += "_Calo";
  else if( inFileList.Contains("PF") ) fileNameId += "_PF";  
  else if( inFileList.Contains("JPT") ) fileNameId += "_JPT";  

  if( inFileList.Contains("L1Offset") ) fileNameId += "_L1Offset";
  else if( inFileList.Contains("L1FastJet") ) fileNameId += "_L1FastJet";

  if( isData ) fileNameId += "_Data";
  else if( inFileList.Contains("Fall10") ) fileNameId += "_MCFall10";  
  else if( inFileList.Contains("Spring11") ) fileNameId += "_MCSpring11";  
  else if( inFileList.Contains("Summer11") ) fileNameId += "_MCSummer11";  
  else if( inFileList.Contains("Fall11") ) fileNameId += "_MCFall11";  
  else if( inFileList.Contains("Summer12") ) fileNameId += "_MCSummer12";  
  else fileNameId += "_MC";  

  return fileNameId;
}



void produceResolutionTailInput() {

  // ++++ PARAMETER DEFINITIONS ++++++++++++++++++++++++++++++++
  
   // Define input
   sampleTools::Parameters::inFileListData = "input/Analysis2012/Kalibri_Run2012ABC-ReReco13Jul-PromptRecoC_L1FastJet_AK5PF";
   sampleTools::Parameters::inFileListMC = "input/Analysis2012/Kalibri_QCD_Pt15-3000-Flat_pythia_Z2_Summer12_53X_S10_L1FastJet_AK5PF";
   sampleTools::Parameters::runOnData = true;
   sampleTools::Parameters::qcdSampleType = sampleTools::WeightProducer::Flat;
   sampleTools::Parameters::puScenario = sampleTools::WeightProducer::Summer12S10;
   sampleTools::Parameters::puHistFile = "/afs/naf.desy.de/user/m/mschrode/CMSSW_5_3_5/src/RA2Classic/AdditionalInputFiles/DataPileupHistogram_RA2Summer12_190456-203002_ABC.root";
   sampleTools::Parameters::nEvts = 300000;

   // Set ouput names
   // Note: the values of sampleTools::Parameters::histFileData and sampleTools::Parameters::histFileMC
   // will be set automatically by fillResTailInputHists(). They are the input
   // for mergeBins().
   sampleTools::Parameters::histFilePath = ".";
   sampleTools::Parameters::histFilePrefix = "TEST";

   sampleTools::Parameters::histFileData = sampleTools::Parameters::histFilePath+"/"+sampleTools::Parameters::histFilePrefix+"_"+fileNameId(sampleTools::Parameters::inFileListData,true)+".root";
   sampleTools::Parameters::histFileMC = sampleTools::Parameters::histFilePath+"/"+sampleTools::Parameters::histFilePrefix+"_"+fileNameId(sampleTools::Parameters::inFileListMC,false)+".root";


   // Set binning
   sampleTools::Parameters::configAdm = "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/resolutionFit/config/Analysis2012/Binning/BinningAdmin2012_v1.cfg";
   sampleTools::Parameters::configBin = "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/resolutionFit/config/Analysis2012/Binning/Binning2012_v1_skims.cfg";
   sampleTools::Parameters::configBinMerged = "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/resolutionFit/config/Analysis2012/Binning/Binning2012_v1_mergedPtBins.cfg";

   // Histogram names
   sampleTools::Parameters::hNamePtAve = "hPtAve";
   sampleTools::Parameters::hNamePtAsym = "hPtAsym";
   sampleTools::Parameters::hNamePtAsymAbs = "hPtAsymAbs";
   sampleTools::Parameters::hNameResp = "hResp";



   // Declare parameters initialized
   sampleTools::Parameters::isInit = true;
  
  


  
   // ++++ PRODUCE HISTOGRAMS +++++++++++++++++++++++++++++++++++
   sampleTools::Parameters::runOnData = true;
   sampleTools::fillResTailInputHists();
   sampleTools::Parameters::runOnData = false;
   sampleTools::fillResTailInputHists();

   sampleTools::mergePtBins();
}

