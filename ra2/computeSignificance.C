#include <cassert>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"


class Parameters {
public:
  static unsigned int nBins() { return 36; }
  static unsigned int maxYield() { return 50000; }
};


class Observation {
public:
  unsigned int nBins() const { return yields_.size(); }
  std::vector<unsigned int>::const_iterator yields() const { return yields_.begin(); }
  void addBin(unsigned int yield) { yields_.push_back(yield); }

private:
  std::vector<unsigned int> yields_;
};


class Background {
public:
  Background(const TString &bkgName) : name_(bkgName) {};

  TString name() const { return name_; }
  unsigned int nBins() const { return meanPredictions_.size(); }
  std::vector<double>::const_iterator meanPredictions() const { return meanPredictions_.begin(); }
  std::vector<double>::const_iterator uncorrelatedUncerts() const { return uncorrelatedUncerts_.begin(); }
  std::vector<double>::const_iterator correlatedUncerts() const { return correlatedUncerts_.begin(); }

  void addBin(double meanPrediction, double uncorrelatedUncert, double correlatedUncert);
  
  
private:
  const TString name_;

  std::vector<double> meanPredictions_;
  std::vector<double> uncorrelatedUncerts_;
  std::vector<double> correlatedUncerts_;
};


void Background::addBin(double meanPrediction, double uncorrelatedUncert, double correlatedUncert) {
  meanPredictions_.push_back(     meanPrediction     );
  uncorrelatedUncerts_.push_back( uncorrelatedUncert );
  correlatedUncerts_.push_back(   correlatedUncert   );
}



class ToyExperiments {
public:
  ToyExperiments();
  ~ToyExperiments();

  void printSetup() const;
  void run(unsigned int nPredictions, unsigned int nExperiments) const;
  double localPValue(unsigned int bin, unsigned int obs) const;
  void printLocalPValue() const;
  void printGlobalPValueOfLocalFluctuation(unsigned int bin, unsigned int nExperiments) const;
  void printGlobalPValueOfObservation(unsigned int nExperiments) const;
  double testStatistic(const std::vector<unsigned int> &obs) const;

  void addBackground(const Background &bkg);
  void setObservation(const Observation &obs);
  
  
private:
  TRandom* rand_;
  std::vector<double> meanPredictions_;
  std::vector<double> uncorrelatedUncerts_;
  std::vector<double> correlatedUncerts_;
  std::vector<unsigned int> observedYields_;
  std::vector<TH1*> predictedYields_;

  double findMinValidRandomNumberForCorrelatedUncertainties() const;
  std::vector<unsigned int> yields(double localPVal) const;
};


ToyExperiments::ToyExperiments()
  : rand_(new TRandom3(0)) {
  meanPredictions_     = std::vector<double>(Parameters::nBins(),0.);
  uncorrelatedUncerts_ = std::vector<double>(Parameters::nBins(),0.);
  correlatedUncerts_   = std::vector<double>(Parameters::nBins(),0.);
  observedYields_      = std::vector<unsigned int>(Parameters::nBins(),0);
  for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
    TString name = "YieldsBin";
    name += bin;
    TString title = "Bin ";
    title += bin;
    TH1* h = new TH1D(name,title+";Toy Yield",Parameters::maxYield()+1,0,Parameters::maxYield()+1);
    predictedYields_.push_back(h);
  }
}


ToyExperiments::~ToyExperiments() {
  delete rand_;
  for(std::vector<TH1*>::iterator it = predictedYields_.begin();
      it != predictedYields_.end(); ++it) {
    delete *it;
  }
}


// Throw toy experiments to predict observed yields:
// - Throw nPredictions mean values (prediction) from the background estimates
//   considering their uncertainties
//   - Per prediction, throw nExperiments observed yields from a Poisson
//     distribution with mean value prediction
void ToyExperiments::run(unsigned int nPredictions, unsigned int nExperiments) const {
  std::cout << "Predicting event yields in " << Parameters::nBins() << " bins from " << nPredictions*nExperiments << " toy experiments ...  " << std::flush;

  // Minimal value (in standard deviations) allowed for
  //correlated fluctuation
  const double minCorr = findMinValidRandomNumberForCorrelatedUncertainties();

  // Throw mean values
  for(unsigned int p = 0; p < nPredictions; ++p) {
    // Throw one (normalized) random number for correlated
    // uncertainties that is valid in all bins
    double rCorr = rand_->Gaus(0.,1.);
    while( rCorr <= minCorr ) {
      rCorr = rand_->Gaus(0.,1.);
    }

    // Loop over all bins and get individual predictions
    for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
      double prediction = -1.;
      bool negativePrediction = true;
      while( negativePrediction ) {
	// Throw one (normalized) random number for uncorrelated
	// uncertainties that is valid in this bin only
	double rUncorr = rand_->Gaus(0.,1.);
	
	// Scale the normalized random numbers by the uncertainties' size
	// to obtain variation of yield
	double uncorrVar = rUncorr * uncorrelatedUncerts_.at(bin);
	double corrVar   = rCorr   * correlatedUncerts_.at(bin);
	
	// Add variations to yield
	prediction = meanPredictions_.at(bin) + uncorrVar + corrVar;
	
	// Check if prediction is positive
	if( prediction >= 0. ) {
	  negativePrediction = false;
	}
      }

      // Throw predicted yields from Poisson with
      // this mean
      for(unsigned int e = 0; e < nExperiments; ++e) {
	predictedYields_.at(bin)->Fill(rand_->Poisson(prediction));
      }
    }
  }
  std::cout << "ok" << std::endl;

  for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
    if( predictedYields_.at(bin)->GetBinContent(Parameters::maxYield()+1) > 0 ) {
      std::cerr << "\n\nWARNING: Overflows in yield histograms!" << std::endl;
      std::cerr << "This is probably safe, but better increase Parameters::maxYield.\n\n" << std::endl;
    }
  }
}


// Number of allowed standard deviations to fluctuate downward
// considering only the correlated uncertainties
double ToyExperiments::findMinValidRandomNumberForCorrelatedUncertainties() const {
  double min = -9999.;
  for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
    if( correlatedUncerts_.at(bin) > 0. ) {
      double minInBin = -1.*meanPredictions_.at(bin)/correlatedUncerts_.at(bin);
      if( minInBin > min ) min = minInBin;
    }
  }
  //  std::cout << "MIN(corr) = " << min << std::endl;

  return min;
}


void ToyExperiments::printSetup() const {
  char txt[100];
  std::vector< std::vector<TString> > table;
  std::vector<TString> headCells;
  headCells.push_back("bin");
  headCells.push_back("bkg");
  headCells.push_back("tot unc");
  headCells.push_back("uncorr unc");
  headCells.push_back("corr unc");
  headCells.push_back("observed");
  table.push_back(headCells);
  for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
    std::vector<TString> cells;
    sprintf(txt,"%d",static_cast<int>(bin));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",meanPredictions_.at(bin));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",sqrt( uncorrelatedUncerts_.at(bin)*uncorrelatedUncerts_.at(bin) + correlatedUncerts_.at(bin)*correlatedUncerts_.at(bin) ));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",uncorrelatedUncerts_.at(bin));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",correlatedUncerts_.at(bin));
    cells.push_back(txt);
    sprintf(txt,"%u",observedYields_.at(bin));
    cells.push_back(txt);

    table.push_back(cells);
  }

  std::vector<int> colWidth(table.front().size(),0);
  for(unsigned int row = 0; row < table.size(); ++row) {
    for(unsigned int col = 0; col < table.at(row).size(); ++col) {
      if( table.at(row).at(col).Length() > colWidth.at(col) ) {
	colWidth.at(col) = table.at(row).at(col).Length();
      }
    }
  }
  for(unsigned int row = 0; row < table.size(); ++row) {
    for(unsigned int col = 0; col < table.at(row).size(); ++col) {
      while( table.at(row).at(col).Length() < colWidth.at(col) ) {
	table.at(row).at(col) = " "+table.at(row).at(col);
      }
    }
  }

  std::cout << "\n\n----- Setup -----" << std::endl;
  for(unsigned int row = 0; row < table.size(); ++row) {
    std::cout << " | ";
    for(unsigned int col = 0; col < table.at(row).size(); ++col) {
      std::cout << table.at(row).at(col) << " | ";
    }
    std::cout << std::endl;
  }
}



void ToyExperiments::printLocalPValue() const {
  char txt[100];
  std::vector< std::vector<TString> > table;
  std::vector<TString> headCells;
  headCells.push_back("bin");
  headCells.push_back("bkg");
  headCells.push_back("tot unc");
  headCells.push_back("uncorr unc");
  headCells.push_back("corr unc");
  headCells.push_back("observed");
  headCells.push_back("p-value");
  headCells.push_back("significance");
  table.push_back(headCells);
  for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
    std::vector<TString> cells;
    sprintf(txt,"%d",static_cast<int>(bin));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",meanPredictions_.at(bin));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",sqrt( uncorrelatedUncerts_.at(bin)*uncorrelatedUncerts_.at(bin) + correlatedUncerts_.at(bin)*correlatedUncerts_.at(bin) ));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",uncorrelatedUncerts_.at(bin));
    cells.push_back(txt);
    sprintf(txt,"%.2lf",correlatedUncerts_.at(bin));
    cells.push_back(txt);
    sprintf(txt,"%u",observedYields_.at(bin));
    cells.push_back(txt);
    double lpv = localPValue(bin,observedYields_.at(bin));
    sprintf(txt,"%.5lf",lpv);
    cells.push_back(txt);
    sprintf(txt,"%.5lf",TMath::NormQuantile(1.-lpv));
    cells.push_back(txt);

    table.push_back(cells);
  }

  std::vector<int> colWidth(table.front().size(),0);
  for(unsigned int row = 0; row < table.size(); ++row) {
    for(unsigned int col = 0; col < table.at(row).size(); ++col) {
      if( table.at(row).at(col).Length() > colWidth.at(col) ) {
	colWidth.at(col) = table.at(row).at(col).Length();
      }
    }
  }
  for(unsigned int row = 0; row < table.size(); ++row) {
    for(unsigned int col = 0; col < table.at(row).size(); ++col) {
      while( table.at(row).at(col).Length() < colWidth.at(col) ) {
	table.at(row).at(col) = " "+table.at(row).at(col);
      }
    }
  }

  std::cout << "\n\n----- Observed local p-values -----" << std::endl;
  for(unsigned int row = 0; row < table.size(); ++row) {
    std::cout << " | ";
    for(unsigned int col = 0; col < table.at(row).size(); ++col) {
      std::cout << table.at(row).at(col) << " | ";
    }
    std::cout << std::endl;
  }
}


double ToyExperiments::localPValue(unsigned int bin, unsigned int obs) const {
  assert( bin < Parameters::nBins() );
  int binNObs = predictedYields_.at(bin)->FindBin(obs);
  int binTot  = predictedYields_.at(bin)->GetNbinsX()+1;
  
  return predictedYields_.at(bin)->Integral(binNObs,binTot)/predictedYields_.at(bin)->Integral(1,binTot);
}



// PValue for finding a local p-value as observed in 'bin' or worse 
void ToyExperiments::printGlobalPValueOfLocalFluctuation(unsigned int bin, unsigned int nExperiments) const {
  std::cout << "Determining global p-value for observed fluctuation in bin " << bin << " from " << nExperiments << " toy experiments ...  " << std::flush;

  // Find the predicted yields that correspond
  // to the local p-value 'localPValue'
  std::vector<unsigned int> limitYields = yields(localPValue(bin,observedYields_.at(bin)));

  TH1* hIsAbovePValue = new TH1D("hIsAbovePValue","",2,0,2);

  const double minCorr = findMinValidRandomNumberForCorrelatedUncertainties();

  for(unsigned int p = 0; p < nExperiments; ++p) {
    bool isAbovePValue = false;
    double rCorr = rand_->Gaus(0.,1.);
    while( rCorr <= minCorr ) {
      rCorr = rand_->Gaus(0.,1.);
    }
    for(unsigned int b = 0; b < Parameters::nBins(); ++b) {
      double prediction = -1.;
      bool negativePrediction = true;
      while( negativePrediction ) {
	double rUncorr = rand_->Gaus(0.,1.);
	double uncorrVar = rUncorr * uncorrelatedUncerts_.at(b);
	double corrVar   = rCorr   * correlatedUncerts_.at(b);
	prediction = meanPredictions_.at(b) + uncorrVar + corrVar;
	if( prediction >= 0. ) {
	  negativePrediction = false;
	}
      }
      double predictedYield = rand_->Poisson(prediction);
      if( predictedYield >= limitYields.at(b) ) {
	isAbovePValue = true;
	break;
      }      
    }
    if( isAbovePValue ) {
      hIsAbovePValue->Fill(1);
    } else {
      hIsAbovePValue->Fill(0);
    }
  }
  std::cout << "ok" << std::endl;

  double lpv      = localPValue(bin,observedYields_.at(bin));
  double gpUncorr = 1. - pow(1.-lpv,Parameters::nBins());
  double gpCorr   = hIsAbovePValue->Integral(2,2)/hIsAbovePValue->Integral(1,2);

  std::cout << "\n\n----- Global p-value for observed fluctuation in bin " << bin << " -----" << std::endl;
  std::cout << "  local p-value                           : " << lpv << " (" << TMath::NormQuantile(1.-lpv) << "sig)" << std::endl;
  std::cout << "  global p-value (without correlations)   : " << gpUncorr  << " (" << TMath::NormQuantile(1.-gpUncorr) << "sig)" << std::endl;
  std::cout << "  global p-value (including correlations) : " << gpCorr  << " (" << TMath::NormQuantile(1.-gpCorr) << "sig)" << std::endl;
}


// For each bin, find the yield that corresponds to a local
// p-value of 'localPValue'
std::vector<unsigned int> ToyExperiments::yields(double localPVal) const {
  std::vector<unsigned int> result;
  for(std::vector<TH1*>::const_iterator iY = predictedYields_.begin();
      iY != predictedYields_.end(); ++iY) {
    int maxHistBin   = (*iY)->GetNbinsX()+1;
    int yieldHistBin = maxHistBin;
    double norm = (*iY)->Integral(1,maxHistBin);
    if( norm > 0. ) {
      double alpha = (*iY)->Integral(yieldHistBin,maxHistBin)/norm;
      while( alpha < localPVal && yieldHistBin > 1 ) {
	--yieldHistBin;
	alpha = (*iY)->Integral(yieldHistBin,maxHistBin)/norm;
      }
      result.push_back(yieldHistBin-1);
    } else {
      result.push_back(0);
    } 
    //    std::cout << "yield for " << localPVal << ": " << result.back() << std::endl;
  }

  return result;
}


// A test statistic for a set of observed yields (one per bin)
// given the background expectation. The test statistic is
// defined as negative logarithm of the product of the local p-values.
double ToyExperiments::testStatistic(const std::vector<unsigned int> &obs) const {
  double testStat = 0.;
  for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
    testStat += -1.*log(localPValue(bin,obs.at(bin)));
  }

  return testStat;
}


// The probability to find a value of the test statistic equal or
// worse than the one obtained with the observed yields
void ToyExperiments::printGlobalPValueOfObservation(unsigned int nExperiments) const {
  std::cout << "Performing " << nExperiments << " toy experiments to compute global p-value of observation ...  " << std::flush;

  // The test statistic
  TH1* hTestStat = new TH1D("hTestStat","",10000*Parameters::nBins(),0,1000);
  std::vector<unsigned int> obs(Parameters::nBins(),0);

  // Minimal value (in standard deviations) allowed for
  //correlated fluctuation
  const double minCorr = findMinValidRandomNumberForCorrelatedUncertainties();

  // Throw mean values
  for(unsigned int p = 0; p < nExperiments; ++p) {
    // Throw one (normalized) random number for correlated
    // uncertainties that is valid in all bins
    double rCorr = rand_->Gaus(0.,1.);
    while( rCorr <= minCorr ) {
      rCorr = rand_->Gaus(0.,1.);
    }

    // Loop over all bins and get individual predictions
    for(unsigned int bin = 0; bin < Parameters::nBins(); ++bin) {
      double prediction = -1.;
      bool negativePrediction = true;
      while( negativePrediction ) {
	// Throw one (normalized) random number for uncorrelated
	// uncertainties that is valid in this bin only
	double rUncorr = rand_->Gaus(0.,1.);
	
	// Scale the normalized random numbers by the uncertainties' size
	// to obtain variation of yield
	double uncorrVar = rUncorr * uncorrelatedUncerts_.at(bin);
	double corrVar   = rCorr   * correlatedUncerts_.at(bin);
	
	// Add variations to yield
	prediction = meanPredictions_.at(bin) + uncorrVar + corrVar;
	
	// Check if prediction is positive
	if( prediction >= 0. ) {
	  negativePrediction = false;
	}
      }

      // Throw predicted yields from Poisson with
      // this mean
      obs.at(bin) = rand_->Poisson(prediction);
    } // End of loop over bins
    
    hTestStat->Fill(testStatistic(obs));
  }
  std::cout << "ok" << std::endl;
  
  std::cout << "\n\n----- Global p-value of observation -----" << std::endl;
  double tSObs = testStatistic(observedYields_);
  int binTSObs = hTestStat->FindBin(tSObs);
  int binTSTot = hTestStat->GetNbinsX()+1; // Include overflows
  double gpVal = hTestStat->Integral(binTSObs,binTSTot)/hTestStat->Integral(1,binTSTot);
  if( binTSObs == binTSTot ) {
    std::cerr << "  ERROR: observed value of test statistic out of histogram range" << std::endl;
  } else if( hTestStat->GetBinContent(binTSTot) > 0 ) {
    std::cerr << "  WARNING: distribution of test statistic has overflows (N = " << hTestStat->GetBinContent(binTSTot) << ")" << std::endl;
  }
  std::cout << "         observed value of test statistic : " << tSObs << std::endl;
  std::cout << "  global p-value (including correlations) : " << gpVal  << " (" << TMath::NormQuantile(1.-gpVal) << "sig)" << std::endl;

  hTestStat->Draw();
  gPad->SetLogy();
  gPad->SaveAs("TestStatistic.pdf");

  delete hTestStat;
}


void ToyExperiments::addBackground(const Background &bkg) {
  assert( bkg.nBins() >= Parameters::nBins() );

  std::vector<double>::iterator iTot;
  std::vector<double>::const_iterator iBkg;

  // Add yields to total background
  iTot = meanPredictions_.begin();
  iBkg = bkg.meanPredictions();
  for(; iTot != meanPredictions_.end(); ++iTot, ++iBkg) {
    *iTot += *iBkg;
  }

  // Add uncertainties in quadrature
  iTot = uncorrelatedUncerts_.begin();
  iBkg = bkg.uncorrelatedUncerts();
  for(; iTot != uncorrelatedUncerts_.end(); ++iTot, ++iBkg) {
    double t = *iTot;
    double b = *iBkg;
    *iTot = sqrt( t*t + b*b );
  }
  iTot = correlatedUncerts_.begin();
  iBkg = bkg.correlatedUncerts();
  for(; iTot != correlatedUncerts_.end(); ++iTot, ++iBkg) {
    double t = *iTot;
    double b = *iBkg;
    *iTot = sqrt( t*t + b*b );
  }
}


void ToyExperiments::setObservation(const Observation &obs) {
  assert( obs.nBins() >= Parameters::nBins() );

  std::vector<unsigned int>::iterator iOld;
  std::vector<unsigned int>::const_iterator iNew;

  // Set observed yields
  iOld = observedYields_.begin();
  iNew = obs.yields();
  for(; iOld != observedYields_.end(); ++iOld, ++iNew) {
    *iOld = *iNew;
  }
}


Background ra2BkgFullyCorrelated() {
  Background bkg("FullyCorrelated");
  bkg.addBin(6024.9,0.,620.4);
  bkg.addBin(2281.0,0.,232.3);
  bkg.addBin( 419.7,0., 55.0);
  bkg.addBin(  57.6,0.,  9.7);
  bkg.addBin( 776.9,0.,101.6);
  bkg.addBin( 330.2,0., 38.3);
  bkg.addBin( 107.5,0., 14.5);
  bkg.addBin(  55.8,0.,  9.0);
  bkg.addBin( 305.7,0., 40.1);
  bkg.addBin( 132.6,0., 17.3);
  bkg.addBin(  32.2,0.,  5.9);
  bkg.addBin(  22.8,0.,  4.7);
  bkg.addBin( 109.8,0., 17.5);
  bkg.addBin(  40.2,0.,  7.1);
  bkg.addBin(  17.6,0.,  4.0);
  bkg.addBin(  77.6,0., 16.1);
  bkg.addBin(  29.6,0.,  5.8);
  bkg.addBin( 300.5,0., 63.4);
  bkg.addBin(  52.4,0., 12.1);
  bkg.addBin(   0.8,0.,  1.7);
  bkg.addBin( 122.6,0., 27.7);
  bkg.addBin(  29.3,0.,  6.6);
  bkg.addBin(   6.1,0.,  2.7);
  bkg.addBin(  64.4,0., 14.6);
  bkg.addBin(  22.6,0.,  5.7);
  bkg.addBin(   2.3,0.,  2.1);
  bkg.addBin(  27.7,0.,  8.1);
  bkg.addBin(   8.8,0.,  3.4);
  bkg.addBin(   0.5,0.,  1.3);
  bkg.addBin(  21.5,0.,  8.1);
  bkg.addBin(   8.0,0.,  3.3);
  bkg.addBin(   4.8,0.,  2.1);
  bkg.addBin(   8.7,0.,  3.3);
  bkg.addBin(   5.8,0.,  2.2);
  bkg.addBin(   6.9,0.,  3.7);
  bkg.addBin(   2.4,0.,  2.8);

  return bkg;
}


Background ra2BkgZInv() {
  Background bkg("ZInv");
  bkg.addBin(	1821.28	,	26.31	,	322.71	);
  bkg.addBin(	993.58	,	19.78	,	175.21	);
  bkg.addBin(	273.21	,	10.31	,	49.51	);
  bkg.addBin(	42.01	,	4.04	,	7.59	);
  bkg.addBin(	215.84	,	9.04	,	38.65	);
  bkg.addBin(	124.08	,	6.96	,	22.48	);
  bkg.addBin(	46.86	,	4.26	,	8.71	);
  bkg.addBin(	35.29	,	3.70	,	6.46	);
  bkg.addBin(	76.28	,	5.37	,	13.71	);
  bkg.addBin(	39.34	,	3.92	,	7.15	);
  bkg.addBin(	18.08	,	2.64	,	3.47	);
  bkg.addBin(	17.81	,	2.63	,	3.31	);
  bkg.addBin(	25.34	,	3.10	,	4.54	);
  bkg.addBin(	16.71	,	2.55	,	3.06	);
  bkg.addBin(	12.35	,	2.18	,	2.36	);
  bkg.addBin(	10.51	,	1.99	,	1.91	);
  bkg.addBin(	10.93	,	2.07	,	1.96	);
  bkg.addBin(	22.72	,	2.78	,	5.40	);
  bkg.addBin(	9.92	,	1.88	,	2.42	);
  bkg.addBin(	0.70	,	0.50	,	0.23	);
  bkg.addBin(	9.09	,	1.75	,	2.22	);
  bkg.addBin(	4.23	,	1.22	,	1.07	);
  bkg.addBin(	1.77	,	0.79	,	0.55	);
  bkg.addBin(	4.37	,	1.21	,	1.06	);
  bkg.addBin(	3.52	,	1.11	,	0.88	);
  bkg.addBin(	1.41	,	0.70	,	0.45	);
  bkg.addBin(	3.33	,	1.05	,	0.84	);
  bkg.addBin(	1.42	,	0.71	,	0.35	);
  bkg.addBin(	0.36	,	0.36	,	0.11	);
  bkg.addBin(	1.33	,	0.67	,	0.33	);
  bkg.addBin(	1.06	,	0.61	,	0.32	);
  bkg.addBin(	0.00	,	0.58	,	0.48	);
  bkg.addBin(	0.64	,	0.45	,	0.25	);
  bkg.addBin(	0.60	,	0.43	,	0.22	);
  bkg.addBin(	0.00	,	0.61	,	0.55	);
  bkg.addBin(	0.00	,	0.55	,	0.45	);

  return bkg;
}

Background ra2BkgLLep() {
  Background bkg("LLep");
  bkg.addBin(	2210.65	,	51.93	,	421.59	);
  bkg.addBin(	660.14	,	27.19	,	125.56	);
  bkg.addBin(	77.34	,	9.21	,	14.83	);
  bkg.addBin(	9.52	,	3.39	,	2.04	);
  bkg.addBin(	277.49	,	18.93	,	53.01	);
  bkg.addBin(	112.80	,	11.99	,	21.74	);
  bkg.addBin(	36.11	,	6.51	,	7.27	);
  bkg.addBin(	9.01	,	3.20	,	1.84	);
  bkg.addBin(	103.54	,	12.15	,	20.33	);
  bkg.addBin(	52.45	,	8.42	,	10.42	);
  bkg.addBin(	6.95	,	2.84	,	1.41	);
  bkg.addBin(	2.35	,	1.68	,	0.50	);
  bkg.addBin(	31.02	,	6.49	,	6.00	);
  bkg.addBin(	10.13	,	3.84	,	2.03	);
  bkg.addBin(	2.28	,	1.62	,	0.49	);
  bkg.addBin(	16.71	,	4.65	,	3.30	);
  bkg.addBin(	9.74	,	3.69	,	2.03	);
  bkg.addBin(	132.48	,	11.86	,	56.83	);
  bkg.addBin(	22.01	,	4.59	,	9.40	);
  bkg.addBin(	0.00	,	0.00	,	1.60	);
  bkg.addBin(	55.76	,	7.82	,	24.01	);
  bkg.addBin(	10.40	,	3.01	,	4.50	);
  bkg.addBin(	2.87	,	2.04	,	1.42	);
  bkg.addBin(	24.05	,	5.26	,	10.37	);
  bkg.addBin(	7.97	,	3.02	,	3.56	);
  bkg.addBin(	0.00	,	0.00	,	1.80	);
  bkg.addBin(	11.53	,	4.08	,	5.05	);
  bkg.addBin(	3.52	,	2.04	,	1.63	);
  bkg.addBin(	0.00	,	0.00	,	1.20	);
  bkg.addBin(	10.04	,	5.02	,	4.66	);
  bkg.addBin(	3.21	,	2.28	,	1.54	);
  bkg.addBin(	1.86	,	1.31	,	0.78	);
  bkg.addBin(	4.79	,	1.96	,	2.17	);
  bkg.addBin(	1.41	,	1.41	,	0.62	);
  bkg.addBin(	5.08	,	2.54	,	2.36	);
  bkg.addBin(	0.00	,	0.00	,	2.10	);
 
  return bkg;
}


Background ra2BkgHadT() {
  Background bkg("HadT");
  bkg.addBin(	1683.72	,	36.26	,	164.11	);
  bkg.addBin(	591.85	,	22.31	,	57.16	);
  bkg.addBin(	67.60	,	6.71	,	6.59	);
  bkg.addBin(	5.98	,	1.74	,	0.64	);
  bkg.addBin(	191.64	,	12.42	,	19.25	);
  bkg.addBin(	83.29	,	7.33	,	8.295	);
  bkg.addBin(	23.65	,	3.06	,	2.32	);
  bkg.addBin(	11.39	,	2.97	,	1.12	);
  bkg.addBin(	66.79	,	6.82	,	7.155	);
  bkg.addBin(	35.70	,	4.96	,	3.62	);
  bkg.addBin(	6.60	,	2.03	,	0.66	);
  bkg.addBin(	2.54	,	0.93	,	0.245	);
  bkg.addBin(	22.21	,	3.24	,	2.07	);
  bkg.addBin(	11.06	,	3.37	,	1.19	);
  bkg.addBin(	2.76	,	1.45	,	0.315	);
  bkg.addBin(	15.22	,	2.91	,	1.725	);
  bkg.addBin(	6.45	,	1.84	,	0.635	);
  bkg.addBin(	127.08	,	7.85	,	19.85	);
  bkg.addBin(	18.56	,	3.12	,	2.895	);
  bkg.addBin(	0.11	,	0.25	,	0.0215	);
  bkg.addBin(	44.65	,	4.25	,	7	);
  bkg.addBin(	12.76	,	2.37	,	2	);
  bkg.addBin(	1.32	,	0.46	,	0.235	);
  bkg.addBin(	24.02	,	3.88	,	3.82	);
  bkg.addBin(	9.55	,	2.03	,	1.53	);
  bkg.addBin(	0.76	,	0.53	,	0.14	);
  bkg.addBin(	6.06	,	2.23	,	1.02	);
  bkg.addBin(	2.93	,	1.41	,	0.5	);
  bkg.addBin(	0.07	,	0.2	,	0.0225	);
  bkg.addBin(	2.30	,	1.23	,	0.455	);
  bkg.addBin(	2.87	,	1.13	,	0.48	);
  bkg.addBin(	2.76	,	1.15	,	0.61	);
  bkg.addBin(	2.68	,	0.9	,	0.59	);
  bkg.addBin(	3.08	,	1	,	0.675	);
  bkg.addBin(	1.32	,	0.71	,	0.345	);
  bkg.addBin(	1.53	,	0.95	,	0.435	);

  return bkg;
}

Background ra2BkgQCD() {
  Background bkg("QCD");
  bkg.addBin(	307.40	,	18.48	,	219.93	);
  bkg.addBin(	34.45	,	5.84	,	22.42	);
  bkg.addBin(	1.32	,	1.15	,	0.99	);
  bkg.addBin(	0.09	,	0.29	,	0.09	);
  bkg.addBin(	91.68	,	10.2	,	64.65	);
  bkg.addBin(	9.94	,	3.22	,	6.72	);
  bkg.addBin(	0.84	,	0.91	,	0.93	);
  bkg.addBin(	0.12	,	0.35	,	0.12	);
  bkg.addBin(	59.00	,	7.2	,	24.93	);
  bkg.addBin(	5.09	,	2.16	,	1.9	);
  bkg.addBin(	0.51	,	0.67	,	0.23	);
  bkg.addBin(	0.12	,	0.33	,	0.07	);
  bkg.addBin(	31.17	,	5.3	,	11.96	);
  bkg.addBin(	2.27	,	1.27	,	1	);
  bkg.addBin(	0.24	,	0.5	,	0.13	);
  bkg.addBin(	35.14	,	6.12	,	11.89	);
  bkg.addBin(	2.43	,	1.37	,	1.42	);
  bkg.addBin(	18.19	,	3.85	,	8.5	);
  bkg.addBin(	1.85	,	1.37	,	0.97	);
  bkg.addBin(	0.01	,	0.1	,	0.06	);
  bkg.addBin(	13.13	,	3.4	,	5.98	);
  bkg.addBin(	1.95	,	1.14	,	0.8	);
  bkg.addBin(	0.15	,	0.39	,	0.11	);
  bkg.addBin(	11.92	,	3.84	,	4.36	);
  bkg.addBin(	1.54	,	1.31	,	0.72	);
  bkg.addBin(	0.11	,	0.32	,	0.16	);
  bkg.addBin(	6.82	,	2.95	,	2.54	);
  bkg.addBin(	0.92	,	1.03	,	0.72	);
  bkg.addBin(	0.09	,	0.29	,	0.06	);
  bkg.addBin(	8.00	,	2.84	,	3.13	);
  bkg.addBin(	0.83	,	0.92	,	0.56	);
  bkg.addBin(	0.14	,	0.38	,	0.19	);
  bkg.addBin(	0.54	,	0.69	,	0.56	);
  bkg.addBin(	0.73	,	0.78	,	0.59	);
  bkg.addBin(	0.54	,	0.75	,	0.52	);
  bkg.addBin(	0.89	,	0.94	,	0.95	);

  return bkg;
}


Observation ra2Data() {
  Observation obs;
  obs.addBin(6159);
  obs.addBin(2305);
  obs.addBin(454);
  obs.addBin(62);
  obs.addBin(808);
  obs.addBin(305);
  obs.addBin(124);
  obs.addBin(52);
  obs.addBin(335);
  obs.addBin(129);
  obs.addBin(34);
  obs.addBin(32);
  obs.addBin(98);
  obs.addBin(38);
  obs.addBin(23);
  obs.addBin(94);
  obs.addBin(39);
  obs.addBin(266);
  obs.addBin(62);
  obs.addBin(9);
  obs.addBin(111);
  obs.addBin(35);
  obs.addBin(4);
  obs.addBin(67);
  obs.addBin(20);
  obs.addBin(4);
  obs.addBin(24);
  obs.addBin(5);
  obs.addBin(2);
  obs.addBin(18);
  obs.addBin(3);
  obs.addBin(8);
  obs.addBin(9);
  obs.addBin(8);
  obs.addBin(5);
  obs.addBin(2);

  return obs;
}



void sig(unsigned int n = 5000) {
  // Test
  //   Background bkg("TEST");
  //   bkg.addBin(6,0,1.1);
  //   bkg.addBin(6,0,1.1);
  //   bkg.addBin(6,0,1.1);
  //   Observation obs;
  //   obs.addBin(9);
  //   obs.addBin(9);
  //   obs.addBin(9);

  //   Background bkg("TEST");
  //   bkg.addBin(0.8,0,1.7);
  //   bkg.addBin(0.8,0,1.7);
  //   bkg.addBin(0.8,0,1.7);
  //   Observation obs;
  //   obs.addBin(9);
  //   obs.addBin(9);
  //   obs.addBin(9);


  ToyExperiments toys;
  toys.addBackground(ra2BkgZInv());
  toys.addBackground(ra2BkgLLep());
  toys.addBackground(ra2BkgHadT());
  toys.addBackground(ra2BkgQCD());
  toys.setObservation(ra2Data());
  toys.printSetup();
  toys.run(n,n);
  toys.printLocalPValue();
  //  toys.printGlobalPValueOfLocalFluctuation(19,n*n);
  toys.printGlobalPValueOfObservation(n*n);
}
