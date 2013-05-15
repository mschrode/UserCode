// $Id: produceResolutionTailInput.C,v 1.1 2013/05/08 13:31:25 mschrode Exp $

#include <cassert>

#include "TString.h"

#include "WeightProducer.h"
#include "Parameters.h"
#include "fillResTailInputHists.C"
#include "mergePtBins.C"


using namespace sampleTools;

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
  
  Parameters par;
  
  // Define input
  par.setInFileListData("input/Analysis2012/Kalibri_Run2012ABC-ReReco13Jul-PromptRecoC_L1FastJet_AK5PF");
  par.setInFileListMC("input/Analysis2012/Kalibri_QCD_Pt15-3000-Flat_pythia_Z2_Summer12_53X_S10_L1FastJet_AK5PF");
  par.setQCDSampleType(WeightProducer::Flat);
  par.setPUScenario(WeightProducer::Summer12S10);
  par.setPUHistFile("/afs/naf.desy.de/user/m/mschrode/CMSSW_5_3_5/src/RA2Classic/AdditionalInputFiles/DataPileupHistogram_RA2Summer12_190456-203002_ABC.root");
  int nEvts = 300000;


  // Set ouput names
  par.setHistFilePath(".");
  par.setHistFilePrefix("TEST");

  par.setHistFileData(par.histFilePath()+"/"+par.histFilePrefix()+"_"+fileNameId(par.inFileListData(),true)+".root");
  par.setHistFileMC(par.histFilePath()+"/"+par.histFilePrefix()+"_"+fileNameId(par.inFileListMC(),false)+".root");


  // Set binning
  par.setConfigAdm("/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/resolutionFit/config/Analysis2012/Binning/BinningAdmin2012_v1.cfg");
  par.setConfigBin("/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/resolutionFit/config/Analysis2012/Binning/Binning2012_v1_skims.cfg");
  par.setConfigBinMerged("/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/resolutionFit/config/Analysis2012/Binning/Binning2012_v1_mergedPtBins.cfg");


  // Setup binning objects
  std::cout << "Defining binning and trigger thresholds" << std::endl;
  BinningAdmin binAdmin(par.configAdm(),par.configBin());
  binAdmin.print();
  binAdmin.printBinning();
  BinningAdmin binAdminMerged(par.configAdm(),par.configBinMerged());
  assert( binAdmin.nPtSoftBins() == binAdminMerged.nPtSoftBins() );
  assert( binAdmin.nEtaBins() == binAdminMerged.nEtaBins() );



  
  // ++++ PRODUCE HISTOGRAMS +++++++++++++++++++++++++++++++++++
  // First, run on data
  std::cout << "\n\n***** (1/3) RUNNING ON DATA *****************************\n\n" << std::endl;
  fillResTailInputHists(par,binAdmin,true,nEvts);
  // Then, run on MC
  std::cout << "\n\n\n\n***** (2/3) RUNNING ON MC *******************************\n\n" << std::endl;
  fillResTailInputHists(par,binAdmin,false,nEvts);
  // Finally, merge output and scale MC to data
  std::cout << "\n\n\n\n***** (3/3) MERGING BINS ********************************\n\n" << std::endl;
  binAdminMerged.printBinning();
  mergePtBins(par,binAdmin,binAdminMerged);
  std::cout << "\n\nDone." << std::endl;
}

