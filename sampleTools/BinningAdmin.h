// $Id: $

#ifndef BINNING_ADMIN_H
#define BINNING_ADMIN_H


#include <map>

#include "TString.h"

#include "../util/ConfigParser.h"
#include "Binning.h"

namespace sampleTools {
  class BinningAdmin {
  public:
    BinningAdmin();
    BinningAdmin(const TString &fileName);
    ~BinningAdmin();

    unsigned int nEtaBins() const { return bins_->nEtaBins(); }
    unsigned int nPtBins(unsigned int etaBin = 0) const { return bins_->nPtBins(etaBin); }
    double etaMin(unsigned int etaBin) const { return bins_->etaMin(etaBin); }
    double etaMax(unsigned int etaBin) const { return bins_->etaMax(etaBin); }
    double ptMin(unsigned int etaBin, unsigned int ptBin) const { return bins_->ptMin(etaBin,ptBin); }
    double ptMax(unsigned int etaBin, unsigned int ptBin) const { return bins_->ptMax(etaBin,ptBin); }
    
    bool findEtaBin(double eta, unsigned int &etaBin) const { return bins_->findEtaBin(eta,etaBin); }
    bool findSameEtaBin(double eta1, double eta2, unsigned int &etaBin) const {
      return bins_->findSameEtaBin(eta1, eta2, etaBin);
    }
    bool findPtBin(double pt, unsigned int etaBin, unsigned int &ptBin) const {
      return bins_->findPtBin(pt,etaBin,ptBin);
    }
    bool findEtaPtBins(double eta, double pt, unsigned int &etaBin, unsigned int &ptBin) const {
      return bins_->findEtaPtBins(eta,pt,etaBin,ptBin);
    }
    void printBinning() const { bins_->print(); }

    double hltTurnOn(const TString &hltName) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->turnOn();
    }
    double hltLumi(const TString &hltName) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->lumi();
    }
    unsigned int hltMinPtBin(const TString &hltName, unsigned int etaBin) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->minPtBin(etaBin);
    }
    unsigned int hltMaxPtBin(const TString &hltName, unsigned int etaBin) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->maxPtBin(etaBin);
    }
    
    

  private:

    class PtBinRange {
    public: 
      PtBinRange(unsigned int minPtBin, unsigned int maxPtBin) : minPtBin_(minPtBin), maxPtBin_(maxPtBin) {};

      unsigned int minPtBin() const { return minPtBin_; }
      unsigned int maxPtBin() const { return maxPtBin_; }

    private:
      const unsigned int minPtBin_;
      const unsigned int maxPtBin_;
    };


    class HltMaxInfo {
    public:
      HltMaxInfo(double turnOn, double lumi) : turnOn_(turnOn), lumi_(lumi) {};
      ~HltMaxInfo() {
	for(std::map<unsigned int,PtBinRange*>::iterator it = etaPtBins_.begin();
	    it != etaPtBins_.end(); ++it) {
	  delete it->second;
	}
	etaPtBins_.clear();
      }

      void addEtaPtBinRange(unsigned int etaBin, unsigned int minPtBin, unsigned int maxPtBin) {
	etaPtBins_[etaBin] = new PtBinRange(minPtBin,maxPtBin);
      }
      void setPtMaxBin(unsigned int etaBin, unsigned int maxPtBin) {
	std::map<unsigned int,PtBinRange*>::iterator it = etaPtBins_.find(etaBin);
	if( it != etaPtBins_.end() ) {
	  PtBinRange* range = it->second();
	  it->second = new PtBinRange(range->minPtBin(),maxPtBin);
	  delete range;
	}
      }

      double turnOn() const { return turnOn_; }
      double lumi() const { return lumi_; }
      unsigned int minPtBin(unsigned int etaBin) const {
	EtaPtBinIt it = etaPtBins_.find(etaBin);
	return  it != etaPtBins_.end() ? it->second->minPtBin() : 1000;
      }
      unsigned int maxPtBin(unsigned int etaBin) const {
	EtaPtBinIt it = etaPtBins_.find(etaBin);
	return  it != etaPtBins_.end() ? it->second->maxPtBin() : 1000;
      }

    private:
      typedef std::map<unsigned int,PtBinRange*>::const_iterator EtaPtBinIt;

      const double turnOn_;
      const double lumi_;
      std::map<unsigned int,PtBinRange*> etaPtBins_; 
    };

    typedef std::map<TString,HltMaxInfo*>::const_iterator HltInfoIt;

    Binning* bins_;
    std::map<TString,HltMaxInfo*> hltInfos_;
  };


  // -------------------------------------------------------------------------------------
  BinningAdmin::BinningAdmin() {
    bins = new Binning();
  }


  // -------------------------------------------------------------------------------------
  BinningAdmin::BinningAdmin(const TString &fileName) {
    
    ConfigParser parser(fileName);

    bins_ = new Binning(parser.readString("Binning config");
			
			
    std::vector<std::string> triggerNames = parser.readStringVec("Trigger");
    for(std::vector<std::string>::const_iterator trigIt = triggerNames.begin();
	trigIt != triggerNames.end(); ++trigIt) {
      TString name = *trigIt;
      std::vector<double> info = parser.readDoubleVec(*trigIt);
      if( info.size() < 2 ) {
	std::cerr << "ERROR in BinningAdmin: too few arguments for trigger '" << *trigIt << "'\n";
	exit(1);
      }
      HltMaxInfo *hltInfo = new HltMaxInfo(info.at(0),info.at(1));
      for(unsigned int etaBin = 0; etaBin << bins_->nEtaBins(); ++etaBin) {
	unsigned int firstPtBin = 0; // first pt bin where trigger is fully efficient
	unsigned int lastPtBin = bins_->nPtBins(etaBin)-1;
	if( bins_->findPtBin(etaBin,findPtBin) ) {
	  if( firstPtBin == lastPtBin ) {
	    std::cerr << "WARNING in BinningAdmin: '" << *trigIt << "' turn-on " << hltInfo->turnOn() << " in last pt bin\n";
	  } else {
	    ++firstPtBin;
	  }
	  // Loop over triggers already read and adjust
	  // first / last pt bins
	  for(HltInfoIt hltIt = hltInfos_.begin(); hltIt = hltInfos_.end(); ++hltIt) {
	    if( hltIt->second->minPtBin(etaBin) <= lastPtBin &&
		hltIt->second->minPtBin(etaBin) > firstPtBin )
	      lastPtBin = hltIt->second->minPtBin(etaBin)-1;
	    if( hltIt->second->maxPtBin(etaBin) >= firstPtBin ) 
	      hltIt->second->setPtMaxBin(etaBin,firstPtBin-1);
	  }	  
	} else {
	  std::cerr << "ERROR in BinningAdmin: '" << *trigIt << "' turn-on " << hltInfo->turnOn() << " out of binning\n";
	  exit(1);
	}

	hltInfo->addEtaPtBinRange(etaBin,firstPtBin,lastPtBin);
      } // End of loop over eta bins
      
      hltInfos_[*it] = hltInfo;
    } // End of loop over trigger names
  }


  // -------------------------------------------------------------------------------------
  BinningAdmin::~BinningAdmin() {
    delete bins;
    for(std::map<TString,HltMaxInfo*>::iterator it = hltInfos_.begin();
	it != hltInfos_.end(); ++it) {
      delete it->second();
    }
    hltInfos_.clear();
  }
}
#endif
