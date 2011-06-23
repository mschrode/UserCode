// $Id: BinningAdmin.h,v 1.6 2011/05/20 10:00:07 mschrode Exp $

#ifndef BINNING_ADMIN_H
#define BINNING_ADMIN_H

#include <cassert>
#include <map>

#include "TString.h"

#include "../util/ConfigParser.h"
#include "Binning.h"

namespace sampleTools {

  //! Administrate eta and pt binning and
  //! HLT thresholds
  // -------------------------------------------------------------------------------------
  class BinningAdmin {
  public:
    inline BinningAdmin();
    inline BinningAdmin(const TString &fileName);
    inline ~BinningAdmin();

    inline unsigned int nEtaBins() const { return bins_->nEtaBins(); }
    inline unsigned int nPtBins(unsigned int etaBin) const { return bins_->nPtBins(etaBin); }
    inline unsigned int nPtSoftBins() const { return ptSoftMin_.size(); }
    inline double etaMin(unsigned int etaBin) const { return bins_->etaMin(etaBin); }
    inline double etaMax(unsigned int etaBin) const { return bins_->etaMax(etaBin); }
    inline double ptMin(unsigned int etaBin, unsigned int ptBin) const { return bins_->ptMin(etaBin,ptBin); }
    inline double ptMax(unsigned int etaBin, unsigned int ptBin) const { return bins_->ptMax(etaBin,ptBin); }
    inline double ptMin(unsigned int etaBin) const { return bins_->ptMin(etaBin); }
    inline double ptMax(unsigned int etaBin) const { return bins_->ptMax(etaBin); }
    inline const std::vector<double> ptBinEdges(unsigned int etaBin) const { return bins_->ptBinEdges(etaBin); }
    inline const std::vector<double> ptBinEdgesInt(unsigned int etaBin) const { return bins_->ptBinEdgesInt(etaBin);}
    inline double ptSoftMin(unsigned int ptSoftBin) const { return ptSoftMin_.at(ptSoftBin); }
    inline double ptSoftMax(unsigned int ptSoftBin) const { return ptSoftMax_.at(ptSoftBin); }
    
    inline bool findEtaBin(double eta, unsigned int &etaBin) const { return bins_->findEtaBin(eta,etaBin); }
    inline bool findSameEtaBin(double eta1, double eta2, unsigned int &etaBin) const {
      return bins_->findSameEtaBin(eta1, eta2, etaBin);
    }
    inline bool findPtBin(double pt, unsigned int etaBin, unsigned int &ptBin) const {
      return bins_->findPtBin(pt,etaBin,ptBin);
    }
    inline bool findEtaPtBins(double eta, double pt, unsigned int &etaBin, unsigned int &ptBin) const {
      return bins_->findEtaPtBins(eta,pt,etaBin,ptBin);
    }
    inline unsigned int nPtBins(const TString &hltName, unsigned int etaBin) const;
    inline bool findPtBin(const TString &hltName, double pt, unsigned int etaBin, unsigned int &ptBin) const;
    inline bool findHltMax(unsigned int etaBin, unsigned int ptBin, TString &hltName) const;
    
    inline void printBinning() const { bins_->print(); }
    inline void printPtSoftBins() const;

    inline double hltTurnOn(const TString &hltName) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->turnOn() : 0.;
    }
    inline double hltLumi(const TString &hltName) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->lumi() : 0.;
    }
    inline double hltLumi(unsigned int etaBin, unsigned int ptBin) const {
      TString hltName;
      findHltMax(etaBin,ptBin,hltName);
      return hltLumi(hltName);
    }
    inline unsigned int hltMinPtBin(const TString &hltName, unsigned int etaBin) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->minPtBin(etaBin) : 0;
    }
    inline unsigned int hltMaxPtBin(const TString &hltName, unsigned int etaBin) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->maxPtBin(etaBin) : nPtBins(etaBin)-1;
    }
    inline TString triggerName(double thres) const;

    inline void print() const;
    inline void print(const TString &hlt) const;
    
    

  private:

    class PtBinRange {
    public: 
      inline PtBinRange(unsigned int minBin, unsigned int maxBin) : minPtBin_(minBin), maxPtBin_(maxBin) {};

      inline unsigned int minPtBin() const { return minPtBin_; }
      inline unsigned int maxPtBin() const { return maxPtBin_; }

    private:
      const unsigned int minPtBin_;
      const unsigned int maxPtBin_;
    };


    class HltMaxInfo {
    public:
      inline HltMaxInfo(double hltTurnOn, double hltLumi) : turnOn_(hltTurnOn), lumi_(hltLumi) {};
      inline ~HltMaxInfo() {
	for(std::map<unsigned int,PtBinRange*>::iterator it = etaPtBins_.begin();
	    it != etaPtBins_.end(); ++it) {
	  delete it->second;
	}
	etaPtBins_.clear();
      }

      inline void addEtaPtBinRange(unsigned int etaBinIdx, unsigned int minPtBinIdx, unsigned int maxPtBinIdx) {
	etaPtBins_[etaBinIdx] = new PtBinRange(minPtBinIdx,maxPtBinIdx);
      }
      inline void setPtMaxBin(unsigned int etaBinIdx, unsigned int maxPtBinIdx) {
	std::map<unsigned int,PtBinRange*>::iterator it = etaPtBins_.find(etaBinIdx);
	if( it != etaPtBins_.end() ) {
	  PtBinRange* range = it->second;
	  it->second = new PtBinRange(range->minPtBin(),maxPtBinIdx);
	  delete range;
	}
      }

      inline double turnOn() const { return turnOn_; }
      inline double lumi() const { return lumi_; }
      inline unsigned int minPtBin(unsigned int etaBin) const {
	EtaPtBinIt it = etaPtBins_.find(etaBin);
	return  it != etaPtBins_.end() ? it->second->minPtBin() : 1000;
      }
      inline unsigned int maxPtBin(unsigned int etaBin) const {
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
    std::vector<double> ptSoftMin_;
    std::vector<double> ptSoftMax_;
  };


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
	  std::cerr << "WARNING in BinningAdmin: '" << *trigIt << "' turn-on " << hltInfo->turnOn() << " out of binning\n";
	}
      } // End of loop over eta bins
      
      hltInfos_[*trigIt] = hltInfo;
    } // End of loop over trigger names

    // PtSoft binning:
    ptSoftMin_ = parser.readDoubleVec("Min ptSoft"); 
    ptSoftMax_ = parser.readDoubleVec("Max ptSoft"); 
    if( ptSoftMin_.size() != ptSoftMax_.size() ) {
      if( ptSoftMin_.size() == 0 ) {
	for(unsigned int i = 0; i < ptSoftMax_.size(); ++i) {
	  ptSoftMin_.push_back(0.);
	}
      } else {
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
  void BinningAdmin::printPtSoftBins() const {
    for(unsigned int i = 0; i < nPtSoftBins(); ++i) {
      if( i == 0 ) std::cout << "if" << std::flush;
      else std::cout << "} else if" << std::flush;
      std::cout << "( pt3Bin == " << i << " ) {" << std::endl;
      std::cout << "pt3Cut = " << ptSoftMax(i) << ";" << std::endl;
      std::cout << "pt3CutStr = \"PtSoft" << i << ".root\";" << std::endl;
    }
    std::cout << "}" << std::endl;
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
    
/*     // 2010  */
/*     if( thres == 15 ) hlt = "HLT_DiJetAve15U"; */
/*     else if( thres == 30 ) hlt = "HLT_DiJetAve30U"; */
/*     else if( thres == 50 ) hlt = "HLT_DiJetAve50U"; */
/*     else if( thres == 70 ) hlt = "HLT_DiJetAve70U"; */
/*     else if( thres == 100 ) hlt = "HLT_DiJetAve100U"; */
/*     else if( thres == 140 ) hlt = "HLT_DiJetAve140U"; */

    // 2011
    if( thres == 30 ) hlt = "HLT_DiJetAve30";
    else if( thres == 60 ) hlt = "HLT_DiJetAve60";
    else if( thres == 80 ) hlt = "HLT_DiJetAve80";
    else if( thres == 110 ) hlt = "HLT_DiJetAve110";
    else if( thres == 150 ) hlt = "HLT_DiJetAve150";
    else if( thres == 190 ) hlt = "HLT_DiJetAve190";
    else if( thres == 240 ) hlt = "HLT_DiJetAve240";
    else if( thres == 300 ) hlt = "HLT_DiJetAve300";
    else if( thres == 370 ) hlt = "HLT_DiJetAve370";

    return hltInfos_.find(hlt) != hltInfos_.end() ? hlt : "none";
  }
}
#endif
