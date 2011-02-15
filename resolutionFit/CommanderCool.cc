// $Id: $

#include <iomanip>
#include <iostream>

#include "CommanderCool.h"
#include "PtBin.h"

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


  // This shoudl take a config file holding the parameters as input
  void CommanderCool::setMCTruthResolution(ResolutionFunction::Type type) {
    if( type == ResolutionFunction::ModifiedNSC ) {
      std::cout << "WARNING in CommanderCool::setMCTruthResolution(): using hard-coded parameters!" << std::endl;
      etaBins_.at(0)->setMCTruthResolution(new ResolutionFunctionModifiedNSC(10.,1000.,3.8663,0.728714,0.,0.224013));
    } else {
      std::cerr << "ERROR in CommanderCool::setMCTruthResolution(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
  }


  // This shoudl take a config file holding the parameters as input
  void CommanderCool::setPLI(ResolutionFunction::Type type) {
    if( type == ResolutionFunction::NSC ) {
      std::cout << "WARNING in CommanderCool::setPLI(): using hard-coded parameters!" << std::endl;
      etaBins_.at(0)->setPLI(new ResolutionFunctionNSC(10.,1000.,2.54877,0.149045,0.0109168));
    } else {
      std::cerr << "ERROR in CommanderCool::setPLI(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
  }



  void CommanderCool::addDataSample(const TString &label, const TString &baseFileName) {
    for(EtaBinIt it = etaBins_.begin(); it != etaBins_.end(); ++it) {
      if( !((*it)->addDataSample(label,baseFileName)) ) {
	std::cerr << "ERROR adding DataSample '" << label << "'" << std::endl;
	exit(1);
      }
    }
  }



  void CommanderCool::addMCSample(const TString &label, const TString &baseFileName) {
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
      }
    }

    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
  }
}
