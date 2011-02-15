// $Id: Binning.h,v 1.3 2010/11/28 22:48:26 mschrode Exp $

#ifndef BINNING_H
#define BINNING_H

#include <cassert>
#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>

#include "TString.h"

namespace sampleTools {

  //! An eta and pt binning
  class Binning {

  public:
     Binning();
     Binning(const std::vector<double> &etaBinEdges, const std::vector<double> &ptBinEdges);
     Binning(const std::vector<double> &etaBinEdges, const std::vector< std::vector<double> > &ptBinEdges);
     Binning(const TString &fileName);
    
     unsigned int nEtaBins() const { return etaBinEdges_.size()-1; }
     unsigned int nPtBins(unsigned int etaBin = 0) const { return ptBinEdges_.at(etaBin).size()-1; }

     double etaMin(unsigned int etaBin) const { return etaBinEdges_.at(etaBin); }
     double etaMax(unsigned int etaBin) const { return etaBinEdges_.at(etaBin+1); }
     double etaEdge(unsigned int etaBinEdge) const { return etaBinEdges_.at(etaBinEdge); }
     double ptMin(unsigned int etaBin) const { return ptBinEdges_.at(etaBin).front(); }
     double ptMax(unsigned int etaBin) const { return ptBinEdges_.at(etaBin).back(); }
     double ptMin(unsigned int etaBin, unsigned int ptBin) const { return ptBinEdges_.at(etaBin).at(ptBin); }
     double ptMax(unsigned int etaBin, unsigned int ptBin) const { return ptBinEdges_.at(etaBin).at(ptBin+1); }
     double ptEdge(unsigned int etaBin, unsigned int ptBinEdge) const { return ptBinEdges_.at(etaBin).at(ptBinEdge); }
     const std::vector<double> ptBinEdges(unsigned int etaBin) const { return ptBinEdges_.at(etaBin); }
     const std::vector<double> ptBinEdgesInt(unsigned int etaBin) const;
    
     bool findEtaBin(double eta, unsigned int &etaBin) const { return findBin(std::abs(eta),etaBinEdges_,etaBin); }
     bool findSameEtaBin(double eta1, double eta2, unsigned int &etaBin) const;
     bool findPtBin(double pt, unsigned int &ptBin) const { return findPtBin(pt,0,ptBin); }
     bool findPtBin(double pt, unsigned int etaBin, unsigned int &ptBin) const { return findBin(pt,ptBinEdges_.at(etaBin),ptBin); }
     bool findEtaPtBins(double eta, double pt, unsigned int &etaBin, unsigned int &ptBin) const;

     void print() const;


  private:
    std::vector<double> etaBinEdges_;
    std::vector< std::vector<double> > ptBinEdges_;
    
    bool findBin(double x, const std::vector<double> &binEdges, unsigned int &bin) const;
    bool isSane() const;
  };
}
#endif
