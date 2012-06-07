// $Id: CommanderCool.cc,v 1.13 2012/06/05 22:44:46 mschrode Exp $

#include <iomanip>
#include <iostream>
#include <set>

#include "CommanderCool.h"
#include "PtBin.h"

#include "../util/ConfigParser.h"


namespace resolutionFit {
  // -------------------------------------------------------------------
  CommanderCool::CommanderCool(const Parameters* par)
    : par_(par) {

    // Create bins
    if( par_->verbosity() > 0 ) std::cout << "CommanderCool::CommanderCool: Creating bins" << std::endl;
    for(unsigned int etaBin = 0; etaBin < par_->nEtaBins(); ++etaBin) {
      etaBins_.push_back(new EtaBin(etaBin,par_->nPtBins(etaBin),par_));
    }

    // Hand over reference to bins to PlotMaker
    plotMaker_ = new PlotMaker(par_,etaBins_);
  }

  

  // -------------------------------------------------------------------
  CommanderCool::~CommanderCool() {
    for(std::vector<EtaBin*>::iterator it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      delete *it;
    }
    delete plotMaker_;
  }



  // -------------------------------------------------------------------
  void CommanderCool::setMCTruthResolution(const TString &fileName, ResolutionFunction::Type type) {
    if( par_->verbosity() > 0 ) std::cout << "CommanderCool::setMCTruthResolution(): Setting MC truth resolution from config file '" << fileName << "'" << std::endl;

    // Parse file for parameters and create MCTruth resolution functions
    util::ConfigParser parser(fileName.Data());
    for(unsigned int etaBin = 0; etaBin < par_->nEtaBins(); ++etaBin) {
      TString key = JetProperties::toString(par_->jetType())+"Eta";
      key += etaBin;
      std::vector<double> param = parser.readDoubleVec(key.Data(),":");
      std::vector<double> paramErr(param.size(),0.);
      etaBins_.at(etaBin)->setMCTruthResolution(ResolutionFunction::createResolutionFunction(type,param,paramErr));
    }
  }



  //! \brief Fit the intrinsic particle level imbalance
  //!
  //! This will add a Sample of type MCTruth, which has the correct
  //! FitResult::Type PtGenAsym. The fitted ptGen resolution is used
  //! as PLI correction to compensate for jet hadronisation and
  //! fragmentation effects.
  // -------------------------------------------------------------------
  void CommanderCool::fitPLI(const TString &label, const TString &fileName, ResolutionFunction::Type type) {
    if( par_->verbosity() > 0 ) std::cout << "CommanderCool::fitPLI: Fitting PLI" << std::endl;

    if( !isConsistentInputName(fileName) ) {
      std::cerr << "WARNING in CommanderCool::fitPLI(): name of file contains a jet type string different than the current type '" << JetProperties::toString(par_->jetType()) << "'" << std::endl;
      exit(1);
    }

    std::cout << "Adding MCTruthSample '" << label << "'" << std::endl;
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      (*it)->fitPLI(label,fileName,type);
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::setPLI(const TString &fileName, ResolutionFunction::Type type) {
    if( par_->verbosity() > 0 ) std::cout << "CommanderCool::setPLI: Setting Particle Level Imbalance from config file '" << fileName << "'" << std::endl;

    // Parse file for parameters and create PLI resolution functions
    util::ConfigParser parser(fileName.Data());
    for(unsigned int etaBin = 0; etaBin < par_->nEtaBins(); ++etaBin) {
      TString key = "Eta";
      key += etaBin;
      std::vector<double> param = parser.readDoubleVec(key.Data(),":");
      std::vector<double> paramErr(param.size(),0.);
      etaBins_.at(etaBin)->setPLI(ResolutionFunction::createResolutionFunction(type,param,paramErr));
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addDataSample(const TString &label, const TString &fileName, const TString &printLabel) {
    if( !isConsistentInputName(fileName) ) {
      std::cerr << "WARNING in CommanderCool::addDataSample(): name of file contains a jet type string different than the current type '" << JetProperties::toString(par_->jetType()) << "'" << std::endl;
      exit(1);
    }
    std::cout << "Adding DataSample '" << label << "'" << std::endl;
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addDataSample(label,fileName,printLabel) ) ) {
	std::cerr << "ERROR adding DataSample '" << label << "'" << std::endl;
	exit(1);
      }
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addMCSample(const TString &label, const TString &fileName, const TString &printLabel) {
    if( !isConsistentInputName(fileName) ) {
      std::cerr << "WARNING in CommanderCool::addMCSample(): name of file contains a jet type string different than the current type '" << JetProperties::toString(par_->jetType()) << "'" << std::endl;
      exit(1);
    }

    std::cout << "Adding MCSample '" << label << "'" << std::endl;
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addMCSample(label,fileName,printLabel) ) ) {
	std::cerr << "ERROR adding MCSample '" << label << "'" << std::endl;
	exit(1);
      }
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addFitResult(FitResult::Type type) {
    std::cout << "Producing fit results of type '" << FitResult::toString(type) << "'" << std::endl;
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addFitResult(type)) ) {
	std::cerr << "ERROR adding FitResult of type '" << FitResult::toString(type) << "'" << std::endl;
	exit(1);
      }
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addExtrapolationUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addExtrapolationUncertainty(nominalSample,type,color)) ) {
	std::cerr << "ERROR adding extrapolation uncertainty" << std::endl;
	exit(1);
      }
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addPLIUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addPLIUncertainty(nominalSample,type,color)) ) {
	std::cerr << "ERROR adding particle level imbalance uncertainty" << std::endl;
	exit(1);
      }
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addMCClosureUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)-> addMCClosureUncertainty(nominalSample,type,color)) ) {
	std::cerr << "ERROR adding MC closure uncertainty" << std::endl;
	exit(1);
      }
    }
  }
   


  // -------------------------------------------------------------------
  void CommanderCool::addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSample, int color) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addUncertaintyFromVariedSample(uncertaintyLabel,fraction,nominalSample,type,variedSample,color)) ) {
	std::cerr << "ERROR adding '" << uncertaintyLabel << "' uncertainty" << std::endl;
	exit(1);
      }
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSampleDown, const TString &variedSampleUp, int color) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addUncertaintyFromVariedSample(uncertaintyLabel,fraction,nominalSample,type,variedSampleDown,variedSampleUp,color)) ) {
	std::cerr << "ERROR adding '" << uncertaintyLabel << "' uncertainty" << std::endl;
	exit(1);
      }
    }
  }



  //! \brief Marks the specified samples for direct comparison plots
  //! 
  //! Also, if only two samples are added (i.e. if this method gets
  //! only called once), k-values can be fitted from the resolution
  //! ratio of 'sample1' and 'sample2' (\sa CommanderCool::fitKValues)
  //! \todo Should be simplified: one vector of samples that are compared
  //! to one reference sample
  // -------------------------------------------------------------------
  void CommanderCool::compareSamples(const SampleLabel &label1, const SampleLabel &label2) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->compareSamples(label1,label2)) ) {
	std::cerr << "ERROR comparing samples" << std::endl;
	exit(1);
      }
    }
  }



  //! \brief Fit the ratio of the measured resolutions of the
  //! samples marked for comparison (\sa CommanderCool::compareSamples)
  // -------------------------------------------------------------------
  void CommanderCool::fitKValues(FitResult::Type type) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      (*it)->fitKValue(type);
    }
  }



  // -------------------------------------------------------------------
  void CommanderCool::printSetup() const {
    std::cout << "\n\n++++++ Setup of CommanderCool ++++++++++++++++++++++++" << std::endl;

    std::cout << "\nSamples:" << std::endl;
    for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
      std::cout << "  " << sTIt->first << " (" << Sample::toString(sTIt->second) << " sample)" << std::endl; 
    }

    std::cout << "\nBinnig:" << std::endl;
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      for(PtBinIt ptBinIt = (*etaBinIt)->ptBinsBegin(); ptBinIt != (*etaBinIt)->ptBinsEnd(); ++ptBinIt) {
	std::cout << "  " << (*ptBinIt)->toTString() << std::endl;
      }
    }

    std::cout << "\nOutput:" << std::endl;
    std::cout << "  " << std::flush;
    if( par_->outMode() == OutputManager::PSAllInOne ) std::cout << "One ps file containing all plots" << std::endl;
    else if( par_->outMode() == OutputManager::EPSSingleFiles ) std::cout << "One eps file per plot" << std::endl;
    else if( par_->outMode() == OutputManager::EPSSingleFilesPlusROOT ) std::cout << "One eps file per plot and all plots in one root file" << std::endl;
    
    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
  }



  // -------------------------------------------------------------------
  void CommanderCool::printResult() const {
    std::cout << "\n\n++++++ Fitted Resolution +++++++++++++++++++++++++++++\n" << std::endl;

    std::cout.setf(std::ios::fixed,std::ios::floatfield);
    for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin();
	rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
      FitResult::Type fitResType = *rIt;
      std::cout << "   Method: '" << FitResult::toString(fitResType) << "'\n\n";

      for(SampleTypeIt sTIt = etaBins_.front()->sampleTypesBegin();
	  sTIt != etaBins_.front()->sampleTypesEnd(); ++sTIt) {
	SampleLabel sampleLabel = sTIt->first;
	Sample::Type sampleType = sTIt->second;
	
	for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	  const EtaBin* etaBin = *etaBinIt;
	  std::cout << "  " << sampleLabel << " (" << Sample::toString(sampleType) << " sample)" << std::endl;
	  
	  for(PtBinIt ptBinIt = etaBin->ptBinsBegin(); ptBinIt != etaBin->ptBinsEnd(); ++ptBinIt) {
	    std::cout << std::setprecision(0) << "   " << (*ptBinIt)->toTString() << ": \t";
	    std::cout << std::setprecision(1) << etaBin->meanPt(sampleLabel,fitResType,(*ptBinIt)->ptBin()) << " \t";
	    std::cout << std::setprecision(4) << etaBin->correctedResolution(sampleLabel,fitResType,(*ptBinIt)->ptBin()) << " +/- ";
	    std::cout << std::setprecision(5) << etaBin->correctedResolutionStatUncert(sampleLabel,fitResType,(*ptBinIt)->ptBin()) << " \t";
	    std::cout << std::setprecision(4) << etaBin->mcTruthResolution(sampleLabel,fitResType,(*ptBinIt)->ptBin()) << std::endl;
	  }
	}
	std::cout << std::endl;
      }
    }

    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;

    int pre = 4;		// precision
    for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin();
	rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
      FitResult::Type fitResType = *rIt;
      
      for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	  sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	SampleLabel sLabel1 = (*sCIt)->label1();
	SampleLabel sLabel2 = (*sCIt)->label2();
	
	if( etaBins_.front()->hasKValue(sLabel1,sLabel2,fitResType) ) {
	  std::cout << "\n\n++++++ " << sLabel1 << " / " << sLabel2 << " ratio +++++++++++++++++++++++++++++\n" << std::endl;

	  std::cout << "   Method: '" << FitResult::toString(fitResType) << "'\n\n";
	  std::cout.setf(std::ios::fixed,std::ios::floatfield);
	  
	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
	    const EtaBin* etaBin = *etaBinIt;

	    std::cout << std::setprecision(0) << "   " << etaBin->toString() << ":   \t$ ";
	    std::cout << std::setprecision(1) << par_->etaMin(etaBin->etaBin()) << "$ -- $" << par_->etaMax(etaBin->etaBin()) << "$ & $";
	    std::cout << std::setprecision(pre) << util::round(etaBin->kValue(sLabel1,sLabel2,fitResType),pre) << " \\pm ";
	    std::cout << std::setprecision(pre) << util::round(etaBin->kStat(sLabel1,sLabel2,fitResType),pre) << " _{-";
	    std::cout << std::setprecision(pre) << util::round(etaBin->kSystDown(sLabel1,sLabel2,fitResType),pre) << "} ^{+";
	    std::cout << std::setprecision(pre) << util::round(etaBin->kSystUp(sLabel1,sLabel2,fitResType),pre) << "} (-";
	    std::cout << std::setprecision(pre) << util::round(etaBin->kTotalDown(sLabel1,sLabel2,fitResType),pre) << " +";
	    std::cout << std::setprecision(pre) << util::round(etaBin->kTotalUp(sLabel1,sLabel2,fitResType),pre) << ") $ \\\\" << std::endl;

	  }
	  std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
	}
      }
    }



//     for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin();
// 	rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
//       FitResult::Type fitResType = *rIt;
      
//       for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
// 	  sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
// 	SampleLabel sLabel1 = (*sCIt)->label1();
// 	SampleLabel sLabel2 = (*sCIt)->label2();
	
// 	if( etaBins_.front()->hasKValue(sLabel1,sLabel2,fitResType) ) {
// 	  std::cout.setf(std::ios::fixed,std::ios::floatfield);
	  
// 	  for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
// 	    const EtaBin* etaBin = *etaBinIt;

// 	    std::cout << std::setprecision(0) << "   " << etaBin->toString() << ":   \t ";
// 	    std::cout << std::setprecision(3) << etaBin->kValue(sLabel1,sLabel2,fitResType)-1. << " (-)";
// 	    std::cout << std::setprecision(4) << std::max(etaBin->kValue(sLabel1,sLabel2,fitResType) - etaBin->kTotalDown(sLabel1,sLabel2,fitResType)-1.,0.) << " (+)";
// 	    std::cout << std::setprecision(4) << std::max(etaBin->kValue(sLabel1,sLabel2,fitResType) + etaBin->kTotalUp(sLabel1,sLabel2,fitResType)-1.,0.) << std::endl;
// 	  }
// 	}
//       }
//     }

  }


  // -------------------------------------------------------------------
  void CommanderCool::printPLI() const {
    std::cout << "\n\n++++++ Particle-level imbalance ++++++++++++++++++++++\n" << std::endl;

    std::cout.setf(std::ios::fixed,std::ios::floatfield);
    for(EtaBinIt etaBinIt = etaBins_.begin(); etaBinIt != etaBins_.end(); ++etaBinIt) {
      const EtaBin* etaBin = *etaBinIt;
      TF1* pli = etaBin->pliFunc("PLIFuncTmp");
      std::cout << std::setprecision(0) << "   " << etaBin->toString() << ":   \t$ ";
      std::cout << std::setprecision(1) << par_->etaMin(etaBin->etaBin()) << "$ -- $" << par_->etaMax(etaBin->etaBin());
      for(int i = 0; i < pli->GetNpar(); ++i) {
	std::cout << std::setprecision(4) << "$ & $" << pli->GetParameter(i) << "$ & $";
	std::cout << std::setprecision(4) << pli->GetParError(i);
      }
      std::cout << "$ \\\\" << std::endl;
      delete pli;
    }

    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
  }



  //! \brief For combination with photon+jet data
  // -------------------------------------------------------------------
  void CommanderCool::printRatioToCombinationFormat() const {
    std::cout << "\n\n++++++ Combination Format +++++++++++++++++++++++++++++\n" << std::endl;

    for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin();
	rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
      FitResult::Type fitResType = *rIt;
      
      for(ComparedSamplesIt sCIt = etaBins_.front()->comparedSamplesBegin();
	  sCIt != etaBins_.front()->comparedSamplesEnd(); ++sCIt) {
	SampleLabel sLabel1 = (*sCIt)->label1();
	SampleLabel sLabel2 = (*sCIt)->label2();

	if( etaBins_.front()->hasKValue(sLabel1,sLabel2,fitResType) ) {
	  for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
	    const EtaBin* etaBin = *it;
	    TGraphAsymmErrors* g = etaBin->ratioGraph(sLabel1,sLabel2,fitResType);

	    std::cout << std::setprecision(1) << std::flush;
	    std::cout << "\n\n#Eta " << par_->etaMin(etaBin->etaBin()) << " -- " << par_->etaMax(etaBin->etaBin()) << std::endl;
	    std::cout << std::setprecision(5) << std::flush;
	    std::cout << "MeanPt:" << std::flush;
	    for(unsigned int i = 0; i < g->GetN(); ++i) {
	      std::cout << "  " << g->GetX()[i] << std::flush;
	    }
	    std::cout << "\nMeanPtError:" << std::flush;
	    for(unsigned int i = 0; i < g->GetN(); ++i) {
	      std::cout << "  " << g->GetEXhigh()[i] << std::flush;
	    }
	    std::cout << "\nRatio:" << std::flush;
	    for(unsigned int i = 0; i < g->GetN(); ++i) {
	      std::cout << "  " << g->GetY()[i] << std::flush;
	    }
	    std::cout << "\nRatioError:" << std::flush;
	    for(unsigned int i = 0; i < g->GetN(); ++i) {
	      std::cout << "  " << g->GetEYhigh()[i] << std::flush;
	    }
	    
	    const SystematicUncertainty* un = 0;
	    if( etaBin->findSystematicUncertainty(sLabel2,fitResType,un) ) {
	      std::cout << "\nSystJESUp:" << std::flush;
	      for(unsigned int i = 0; i < g->GetN(); ++i) {
		std::cout << "  " << g->GetY()[i]*un->relUncertUp(i) << std::flush;
	      }
	      std::cout << "\nSystJESDown:" << std::flush;
	      for(unsigned int i = 0; i < g->GetN(); ++i) {
		std::cout << "  " << g->GetY()[i]*un->relUncertDown(i) << std::flush;
	      }
	      std::cout << "\nSystOtherUp:" << std::flush;
	      for(unsigned int i = 0; i < g->GetN(); ++i) {
		std::cout << "  " << g->GetY()[i]*un->relUncertUp(i) << std::flush;
	      }
	      std::cout << "\nSystOtherDown:" << std::flush;
	      for(unsigned int i = 0; i < g->GetN(); ++i) {
		std::cout << "  " << g->GetY()[i]*un->relUncertDown(i) << std::flush;
	      }
	      std::cout << std::endl;
	    }
	  }
	}
      }
    }
  }



  //! \brief Just some means for sanity checks I wrote when I was bored...
  // -------------------------------------------------------------------
  bool CommanderCool::isConsistentInputName(const TString &name) const {
    // Set of all known jet types
    std::set<JetProperties::Type> types;
    types.insert(JetProperties::Calo);
    types.insert(JetProperties::JPT);
    types.insert(JetProperties::PF);
    
    // Remove current jet type
    types.erase(par_->jetType());
    
    // Check whether name contains label of other jet types
    bool result = true;
    for(std::set<JetProperties::Type>::const_iterator it = types.begin(); it != types.end(); ++it) {
      if( name.Contains(JetProperties::toString(*it)) ) {
	result = false;
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------
  void CommanderCool::makeAllPlots() const {
    if( !etaBins_.front()->hasDataSample() ) {
      // If only MC samples, set solid marker styles
      std::cout << "Only MC samples: setting only solid marker styles" << std::endl;
      Sample::setOnlySolidMarkerStyles();
    }
    plotMaker_->makeAllPlots();
  }
}
