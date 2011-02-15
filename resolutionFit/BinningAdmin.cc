#include "BinningAdmin.h"

namespace sampleTools {



  // -------------------------------------------------------------------------------------
  BinningAdmin::BinningAdmin() {
    bins_ = new Binning();
  }


  // -------------------------------------------------------------------------------------
  BinningAdmin::BinningAdmin(const TString &fileName) {
    util::ConfigParser parser(fileName.Data());

    bins_ = new Binning(parser.readString("Binning config"));
			
			
    std::vector<std::string> triggerNames = parser.readStringVec("Trigger");
    for(std::vector<std::string>::const_iterator trigIt = triggerNames.begin();
	trigIt != triggerNames.end(); ++trigIt) {
      std::vector<double> info = parser.readDoubleVec(*trigIt);
      if( info.size() < 2 ) {
	std::cerr << "ERROR in BinningAdmin: too few arguments for trigger '" << *trigIt << "'\n";
	exit(1);
      }
      HltMaxInfo *hltInfo = new HltMaxInfo(info.at(0),info.at(1));
      for(unsigned int etaBin = 0; etaBin < nEtaBins(); ++etaBin) {
	unsigned int firstPtBin = 0; // first pt bin where trigger is fully efficient
	unsigned int lastPtBin = nPtBins(etaBin)-1;
	if( findPtBin(hltInfo->turnOn(),etaBin,firstPtBin) ||
	    hltInfo->turnOn() < ptMin(etaBin,firstPtBin) ) {
	  if( hltInfo->turnOn() > ptMin(etaBin,firstPtBin) ) {
	    // If turn-on is within this bin, move to next higher bin
	    if( firstPtBin < lastPtBin ) {
	      ++firstPtBin;
	    } else {
	      std::cerr << "WARNING in BinningAdmin: '" << *trigIt << "' turn-on " << hltInfo->turnOn() << " in last pt bin\n";
	    }
	  }
	  // Loop over triggers already read and adjust
	  // first / last pt bins
	  for(HltInfoIt hltIt = hltInfos_.begin(); hltIt != hltInfos_.end(); ++hltIt) {
	    if( hltIt->second->minPtBin(etaBin) <= lastPtBin &&
		hltIt->second->minPtBin(etaBin) > firstPtBin &&
		lastPtBin > firstPtBin )
	      lastPtBin = hltIt->second->minPtBin(etaBin)-1;
	    if( hltIt->second->maxPtBin(etaBin) >= firstPtBin &&
		hltIt->second->minPtBin(etaBin) <= firstPtBin-1 )
	      hltIt->second->setPtMaxBin(etaBin,firstPtBin-1);
	  }
	  hltInfo->addEtaPtBinRange(etaBin,firstPtBin,lastPtBin);
	} else {
	  std::cerr << "ERROR in BinningAdmin: '" << *trigIt << "' turn-on " << hltInfo->turnOn() << " out of binning\n";
	  exit(1);
	}
      } // End of loop over eta bins
      
      hltInfos_[*trigIt] = hltInfo;
    } // End of loop over trigger names

    // PtSoft binning:
    ptSoftMin_ = parser.readDoubleVec("Min ptSoft"); 
    ptSoftMax_ = parser.readDoubleVec("Max ptSoft"); 
    if( ptSoftMin_.size() != ptSoftMax_.size() ) {
      std::cerr << "\nERROR in BinningAdmin: 'ptSoftMin_.size() != ptSoftMax_.size()'\n";
      for(unsigned int i = 0; i < ptSoftMin_.size(); ++i) {
 	std::cout << "  ptSoftMin_[" << i << "] = " << ptSoftMin_[i] << std::endl;
      }
      for(unsigned int i = 0; i < ptSoftMax_.size(); ++i) {
 	std::cout << "  ptSoftMax_[" << i << "] = " << ptSoftMax_[i] << std::endl;
      }
      exit(1);
    }
  }


  // -------------------------------------------------------------------------------------
  BinningAdmin::~BinningAdmin() {
    delete bins_;
    for(std::map<TString,HltMaxInfo*>::iterator it = hltInfos_.begin();
	it != hltInfos_.end(); ++it) {
      delete it->second;
    }
    hltInfos_.clear();
  }


  // -------------------------------------------------------------------------------------
  unsigned int BinningAdmin::nPtBins(const TString &hltName, unsigned int etaBin) const {
    unsigned int min = hltMinPtBin(hltName,etaBin);
    unsigned int max = hltMaxPtBin(hltName,etaBin);

    return 1 + max - min;    
  }
  

  // -------------------------------------------------------------------------------------
  bool BinningAdmin::findPtBin(const TString &hltName, double pt, unsigned int etaBin, unsigned int &ptBin) const {
    bool inRange = false;
    if( findPtBin(pt,etaBin,ptBin) ) {
      unsigned int min = hltMinPtBin(hltName,etaBin);
      unsigned int max = hltMaxPtBin(hltName,etaBin);
      if( ptBin >= min && ptBin <= max ) {
	inRange = true;
      }
    }

    return inRange;
  }


  // -------------------------------------------------------------------------------------
  bool BinningAdmin::findHltMax(unsigned int etaBin, unsigned int ptBin, TString &hltName) const {
    bool inRange = false;
    hltName = "none";
    for(HltInfoIt hltIt = hltInfos_.begin(); hltIt != hltInfos_.end(); ++hltIt) {
      if( ptBin >= hltIt->second->minPtBin(etaBin) && ptBin <= hltIt->second->maxPtBin(etaBin) ) {
	hltName = hltIt->first;
	inRange = true;
	break;
      }
    }

    return inRange;
  }


  // -------------------------------------------------------------------------------------
  void BinningAdmin::print() const {
    std::cout << "\n\n ++++++ BINNING ++++++++++++++++++++++++++++\n\n";
    for(unsigned int etaBin = 0; etaBin < nEtaBins(); ++etaBin) {
      std::cout << "Eta " << etaBin << ": " << etaMin(etaBin) << " - " << etaMax(etaBin) << std::endl;
      for(HltInfoIt hltIt = hltInfos_.begin(); hltIt != hltInfos_.end(); ++hltIt) {
	std::cout << "  " << hltIt->first << ":  Pt " << hltIt->second->minPtBin(etaBin) << " - " << hltIt->second->maxPtBin(etaBin) << std::endl;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  void BinningAdmin::print(const TString &hlt) const {
    if( hltInfos_.find(hlt) != hltInfos_.end() ) {
      std::cout << "\nBinning for " << hlt << ":" << std::endl;
      for(unsigned int etaBin = 0; etaBin < nEtaBins(); ++etaBin) {
	unsigned int min = hltMinPtBin(hlt,etaBin);
	unsigned int max = hltMaxPtBin(hlt,etaBin);
	std::cout << "  Eta " << etaBin << " (" << etaMin(etaBin) << " - " << etaMax(etaBin) << "): Pt " << min << " (" << ptMin(etaBin,min) << " GeV) to " << max << " (" << ptMax(etaBin,max) << " GeV)" << std::endl;
      }
    } else {
      std::cerr << "BinningAdmin::print(): No trigger with name '" << hlt << "'\n";
    }
  }


  // -------------------------------------------------------------------------------------
  TString BinningAdmin::triggerName(double thres) const {
    TString hlt;
    if( thres == 50 ) hlt = "HLT_DiJetAve50U";
    else if( thres == 70 ) hlt = "HLT_DiJetAve70U";
    else if( thres == 100 ) hlt = "HLT_DiJetAve100U";
    else if( thres == 140 ) hlt = "HLT_DiJetAve140U";

    return hltInfos_.find(hlt) != hltInfos_.end() ? hlt : "none";
  }

}
