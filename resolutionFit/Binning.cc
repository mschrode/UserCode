#include "Binning.h"

namespace sampleTools {



  //! Default constructor; one eta and pt bin
  // -------------------------------------------------------------------------------------
  Binning::Binning() {
    etaBinEdges_.push_back(0.);
    etaBinEdges_.push_back(10.);

    std::vector<double> tmp;
    tmp.push_back(0.);
    tmp.push_back(10000.);
    ptBinEdges_.push_back(tmp);

    assert( isSane() );
  }



  //! Constructor for 'rectangular' binning i.e. same
  //! pt bins per eta bin
  // -------------------------------------------------------------------------------------
  Binning::Binning(const std::vector<double> &etaBinEdges, const std::vector<double> &ptBinEdges) {
    etaBinEdges_ = etaBinEdges;
    ptBinEdges_.push_back(ptBinEdges);

    assert( isSane() );
  }



  //! Constructor for general binning i.e. different
  //! pt bins per eta bin
  // -------------------------------------------------------------------------------------
  Binning::Binning(const std::vector<double> &etaBinEdges, const std::vector< std::vector<double> > &ptBinEdges) {
    etaBinEdges_ = etaBinEdges;
    ptBinEdges_ = ptBinEdges;

    assert( isSane() );
  }



  //! Constructor for general binning i.e. different
  //! pt bins per eta bin which are read from file
  // -------------------------------------------------------------------------------------
  Binning::Binning(const TString &fileName) {

    std::ifstream file;
    file.open(fileName.Data());
    if( file.is_open() ) {
      double val = 0.;
      std::vector<double> tmpPtBins;
      while( !file.eof() ) {
	tmpPtBins.clear();

	// Read first double
	file >> val;
	// Check for empty last line in file
	if( ptBinEdges_.size() ) {
	  if( val == ptBinEdges_.back().back() ) break;
	}
	// Read eta bin edges
	if( etaBinEdges_.size() == 0 ) { // First eta bin
	  etaBinEdges_.push_back(val);
	  file >> val;
	  etaBinEdges_.push_back(val);
	} else { // other eta bins
	  if( val != etaBinEdges_.back() ) {
	    std::cerr << "ERROR in Binning(): eta bin edges non-continous!" << std::endl;
	    exit(-1);
	  }
	  file >> val;
	  etaBinEdges_.push_back(val);
	}
	// Read pt bin edges
	file >> val;
	unsigned int nPtBinEdges = static_cast<unsigned int>(val);
	for(unsigned int i = 0; i < nPtBinEdges; ++i) {
	  file >> val;
	  tmpPtBins.push_back(val);
	}
	ptBinEdges_.push_back(tmpPtBins);
      }
    } else {
      std::cerr << "ERROR in Binning(): opening file '" << fileName << "'\n";
      exit(1);
    }
    file.close();

    assert( isSane() );
  }


  // -------------------------------------------------------------------------------------
  bool Binning::findSameEtaBin(double eta1, double eta2, unsigned int &etaBin) const { 
    bool sameEtaBin = false;
    unsigned int etaBin2 = 1000;
    if( findEtaBin(eta1,etaBin) && findEtaBin(eta2,etaBin2) ) {
      if( etaBin == etaBin2 ) {
	sameEtaBin = true;
      }
    }

    return sameEtaBin;
  }


  // -------------------------------------------------------------------------------------
  bool Binning::findEtaPtBins(double eta, double pt, unsigned int &etaBin, unsigned int &ptBin) const {
    bool inRange = findBin(std::abs(eta),etaBinEdges_,etaBin);
    if( inRange ) inRange = findBin(pt,ptBinEdges_.at(etaBin),ptBin);
    return inRange;
  }


  // -------------------------------------------------------------------------------------
  void Binning::print() const {
    std::cout << "\n\nETA-PT BINNING\n\n";
    for(unsigned int etaBin = 0; etaBin < nEtaBins(); ++etaBin) {
      std::cout << "  Eta (" << etaMin(etaBin) << ", " << etaMax(etaBin) << "):  \t" << std::flush;
      for(unsigned int ptBin = 0; ptBin < nPtBins(etaBin); ++ptBin) {
	std::cout << ptMin(etaBin,ptBin) << "  " << std::flush;
      }
      std::cout << ptMax(etaBin) << std::endl;
    }
  }


  // -------------------------------------------------------------------------------------
  bool Binning::findBin(double x, const std::vector<double> &binEdges, unsigned int &bin) const {
    bin = 0;
    bool inRange = false;
    if( x >= binEdges.front() && x <= binEdges.back() ) {
      inRange = true;
      for(unsigned int i = 0; i < (binEdges.size()-1); ++i) {
	if( x >= binEdges[i] ) bin = i;
	else break;
      }
    }
    
    return inRange;
  }


  // -------------------------------------------------------------------------------------
  const std::vector<double> Binning::ptBinEdgesInt(unsigned int etaBin) const {
    // Should be done in a smarter way
    std::vector<double> ptBinEdgesInt;
    ptBinEdgesInt.push_back(ptMin(etaBin));
    if( nPtBins(etaBin) > 4  && nPtBins(etaBin) < 10 ) {
      ptBinEdgesInt.push_back(ptMax(etaBin,1));
    } else {
      ptBinEdgesInt.push_back(ptMax(etaBin,2));
      ptBinEdgesInt.push_back(ptMax(etaBin,5));
    }    
    ptBinEdgesInt.push_back(ptMax(etaBin));
    
    return ptBinEdgesInt;
  }


  // -------------------------------------------------------------------------------------
  bool Binning::isSane() const {
    bool sane = true;

    if( etaBinEdges_.size() < 2 ) {
      sane = false;
    } else {
      for(unsigned int i = 1; i < etaBinEdges_.size(); ++i) {
	if( etaBinEdges_[i] <= etaBinEdges_[i-1] ) {
	  sane = false;
	  break;
	}
      }
    }

    if( sane ) {
      if( ptBinEdges_.size() < 1 ) {
 	sane = false;
      } else {
 	for(unsigned int i = 0; i < ptBinEdges_.size(); ++i) {
 	  if( ptBinEdges_[i].size() < 2 ) {
	    sane = false;
	    break;
 	  }
 	}
 	if( sane ) {
 	  for(unsigned int i = 0; i < ptBinEdges_.size(); ++i) {
 	    for(unsigned int j = 1; j < ptBinEdges_.at(i).size(); ++j) {
 	      if( ptBinEdges_[i][j] <= ptBinEdges_[i][j-1] ) {
 		sane = false;
 		break;
 	      }
 	    }
 	  }
 	}
      }
    }

    return sane;
  }
}
