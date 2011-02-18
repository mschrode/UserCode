// $Id: CommanderCool.cc,v 1.1 2011/02/15 18:22:25 mschrode Exp $

#include <iomanip>
#include <iostream>
#include <set>

#include "CommanderCool.h"
#include "PtBin.h"

#include "../util/ConfigParser.h"


namespace resolutionFit {
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

  
  CommanderCool::~CommanderCool() {
    for(std::vector<EtaBin*>::iterator it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      delete *it;
    }
    delete plotMaker_;
  }


  void CommanderCool::setMCTruthResolution(const TString &fileName, ResolutionFunction::Type type) {
    if( par_->verbosity() > 0 ) std::cout << "CommanderCool::setMCTruthResolution(): Setting MC truth resolution from config file '" << fileName << "'" << std::endl;

    // Parse file for parameters and create MCTruth resolution functions
    util::ConfigParser parser(fileName.Data());
    for(unsigned int etaBin = 0; etaBin < par_->nEtaBins(); ++etaBin) {
      TString key = JetProperties::toString(par_->jetType())+"Eta";
      key += etaBin;
      std::vector<double> param = parser.readDoubleVec(key.Data(),":");
      etaBins_.at(etaBin)->setMCTruthResolution(ResolutionFunction::createResolutionFunction(type,param));
    }
  }


  void CommanderCool::setPLI(const TString &fileName, ResolutionFunction::Type type) {
    if( par_->verbosity() > 0 ) std::cout << "CommanderCool::setPLI: Setting Particle Level Imbalance from config file '" << fileName << "'" << std::endl;

    // Parse file for parameters and create PLI resolution functions
    util::ConfigParser parser(fileName.Data());
    for(unsigned int etaBin = 0; etaBin < par_->nEtaBins(); ++etaBin) {
      TString key = "Eta";
      key += etaBin;
      std::vector<double> param = parser.readDoubleVec(key.Data(),":");
      etaBins_.at(etaBin)->setPLI(ResolutionFunction::createResolutionFunction(type,param));
    }
  }



  void CommanderCool::addDataSample(const TString &label, const TString &baseFileName) {
    if( !isConsistentInputName(baseFileName) ) {
      std::cerr << "WARNING in CommanderCool::addDataSample(): name of file contains a jet type string different than the current type '" << JetProperties::toString(par_->jetType()) << "'" << std::endl;
      exit(1);
    }

    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addDataSample(label,baseFileName)) ) {
	std::cerr << "ERROR adding DataSample '" << label << "'" << std::endl;
	exit(1);
      }
    }
  }



  void CommanderCool::addMCSample(const TString &label, const TString &baseFileName) {
    if( !isConsistentInputName(baseFileName) ) {
      std::cerr << "WARNING in CommanderCool::addDataSample(): name of file contains a jet type string different than the current type '" << JetProperties::toString(par_->jetType()) << "'" << std::endl;
      exit(1);
    }

    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addMCSample(label,baseFileName) ) ) {
	std::cerr << "ERROR adding MCSample '" << label << "'" << std::endl;
	exit(1);
      }
    }
  }


  void CommanderCool::addFitResult(FitResult::Type type) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addFitResult(type)) ) {
	std::cerr << "ERROR adding FitResult of type '" << FitResult::toString(type) << "'" << std::endl;
	exit(1);
      }
    }
  }


  void CommanderCool::compareSamples(const SampleLabel &label1, const SampleLabel &label2) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->compareSamples(label1,label2)) ) {
	std::cerr << "ERROR comparing samples" << std::endl;
	exit(1);
      }
    }
  }


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

    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
  }


  void CommanderCool::printResult() const {
    std::cout << "\n\n++++++ Fitted Resolution +++++++++++++++++++++++++++++\n" << std::endl;

    std::cout.setf(std::ios::fixed,std::ios::floatfield);
    for(FitResultTypeIt rIt = etaBins_.front()->fitResultTypesBegin();
	rIt != etaBins_.front()->fitResultTypesEnd(); ++rIt) {
      FitResult::Type fitResType = *rIt;
      std::cout << " " << FitResult::toString(fitResType) << std::endl;

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
  }


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
}
