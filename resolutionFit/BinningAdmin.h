// $Id: BinningAdmin.h,v 1.4 2010/11/28 22:48:26 mschrode Exp $

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
     BinningAdmin();
     BinningAdmin(const TString &fileName);
     ~BinningAdmin();

     unsigned int nEtaBins() const { return bins_->nEtaBins(); }
     unsigned int nPtBins(unsigned int etaBin) const { return bins_->nPtBins(etaBin); }
     unsigned int nPtSoftBins() const { return ptSoftMin_.size(); }
     double etaMin(unsigned int etaBin) const { return bins_->etaMin(etaBin); }
     double etaMax(unsigned int etaBin) const { return bins_->etaMax(etaBin); }
     double ptMin(unsigned int etaBin, unsigned int ptBin) const { return bins_->ptMin(etaBin,ptBin); }
     double ptMax(unsigned int etaBin, unsigned int ptBin) const { return bins_->ptMax(etaBin,ptBin); }
     double ptMin(unsigned int etaBin) const { return bins_->ptMin(etaBin); }
     double ptMax(unsigned int etaBin) const { return bins_->ptMax(etaBin); }
     const std::vector<double> ptBinEdges(unsigned int etaBin) const { return bins_->ptBinEdges(etaBin); }
     const std::vector<double> ptBinEdgesInt(unsigned int etaBin) const { return bins_->ptBinEdgesInt(etaBin);}
     const std::vector<double>& ptSoftMax() const { return ptSoftMax_; }
     double ptSoftMin(unsigned int ptSoftBin) const { return ptSoftMin_.at(ptSoftBin); }
     double ptSoftMax(unsigned int ptSoftBin) const { return ptSoftMax_.at(ptSoftBin); }
    
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
     unsigned int nPtBins(const TString &hltName, unsigned int etaBin) const;
     bool findPtBin(const TString &hltName, double pt, unsigned int etaBin, unsigned int &ptBin) const;
     bool findHltMax(unsigned int etaBin, unsigned int ptBin, TString &hltName) const;
    
     void printBinning() const { bins_->print(); }

     double hltTurnOn(const TString &hltName) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->turnOn() : 0.;
    }
     double hltLumi(const TString &hltName) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->lumi() : 0.;
    }
     double hltLumi(unsigned int etaBin, unsigned int ptBin) const {
      TString hltName;
      findHltMax(etaBin,ptBin,hltName);
      return hltLumi(hltName);
    }
     unsigned int hltMinPtBin(const TString &hltName, unsigned int etaBin) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->minPtBin(etaBin) : 0;
    }
     unsigned int hltMaxPtBin(const TString &hltName, unsigned int etaBin) const {
      HltInfoIt it = hltInfos_.find(hltName);
      return it != hltInfos_.end() ? it->second->maxPtBin(etaBin) : nPtBins(etaBin)-1;
    }
     TString triggerName(double thres) const;

     void print() const;
     void print(const TString &hlt) const;
    
    

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
	  PtBinRange* range = it->second;
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
    std::vector<double> ptSoftMin_;
    std::vector<double> ptSoftMax_;
  };
}
#endif
