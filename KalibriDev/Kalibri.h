//  $Id: Kalibri.h,v 1.5 2010/04/24 14:10:52 mschrode Exp $

//!  \mainpage
//!
//!  \image html kalibriLogoSmall.jpg
//!  Package for data driven calibration using an unbinned fit (see also the related
//!  <A HREF="https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration">
//!  Twiki Page</A>).
//!
//!  \section label_sec_src Source Code
//!  The source code can be found
//!  <A HREF="http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Bromo/Calibration/CalibCore/">here</A>.
//!
//!  \section label_sec_workflow Workflow
//!  \image html kalibri_workflow.png
//!  (Graphic in <A HREF="../graphic/kalibri.eps">eps</A> format.)
//!
//!  \section label_sec_authors Authors
//!  - Christian Autermann
//!  - Ulla Gebbert
//!  - Robert Klanner
//!  - Bj&ouml;rn Kolodzey
//!  - Sebastian Naumann-Emme
//!  - Christian Sander
//!  - Matthias Schr&ouml;der
//!  - Torben Schum
//!  - Hartmut Stadie
//!  - Jan Thomsen
//!  - Roger Wolf

 

#ifndef caliber_h
#define caliber_h

#include <vector>
#include <string>

class TParameters;
class Controlplots;
class Event;
class Measurement;


//!  \brief Main program
//!  \note  For profiling:
//!         To prevent gprof from missing the threads: 
//!         wget http://sam.zoy.org/writings/programming/gprof-helper.c
//!         gcc -shared -fPIC gprof-helper.c -o gprof-helper.so -lpthread -ldl 
//!         LD_PRELOAD=./gprof-helper.so ./junk
//!  \authors Christian Autermann, Hartmut Stadie, Matthias Schroeder
//!  \date Wed Jul 18 13:54:50 CEST 2007
//!  $Id: Kalibri.h,v 1.5 2010/04/24 14:10:52 mschrode Exp $
// -----------------------------------------------------------------
class Kalibri {
public :
  //!  \brief Constructor
  //!  \param f Name of the configuration file
  // -----------------------------------------------------------------
  Kalibri(const std::string& f)
  : configFile_(f),
  par_(0),
  fitMethod_(1),
  nThreads_(1),
  nGammajetEvents_(0),
  nDijetEvents_(0),
  nTrijetEvents_(0),
  nTrackClusterEvents_(0),
  nZjetEvents_(0),
  nTopEvents_(0),
  printParNDeriv_(false),
  derivStep_(1e-03),
  mvec_(6),
  nIter_(100),
  eps_(1e-02),
  wlf1_(1e-04),
  wlf2_(0.9),
  calcCov_(false)
  {};

  ~Kalibri(){};

  void init();  //!< Read parameters from configfile, read data
  void run();    //!< Run the fit
  void done();   //!< Make control plots, clean up
  const char * getOutputFile() { return outputFile_.c_str(); }; //!< Get the ouputfile name

protected:  
  //internal functions
  void run_Lvmini();  //!< Run the fit

private:
  //internal variables
  std::string configFile_;    //!< The configuration file name
  std::string outputFile_;    //!< The output file name
  TParameters * par_;         //!< Fit parameters, depend on number of bins & geometry
  int fitMethod_;             //!< Running mode
  int nThreads_;              //!< Number of threads
  std::vector<Event*> data_;  //!< The data
  int nGammajetEvents_;       //!< Number of gamma-jet events
  int nDijetEvents_;          //!< Number of dijet events
  int nTrijetEvents_;         //!< Number of trijet events
  int nTrackClusterEvents_;   //!< Number of track-cluster events
  int nZjetEvents_;           //!< Number of Zjet events
  int nTopEvents_;            //!< Number of top events
  int mode_;

  // control parameters of fit
  bool printParNDeriv_;     //!< Control whether to print derivatives in each iteration
  std::vector<int> residualScalingScheme_;    //!< Iteration scheme of scaling of residuals
  double outlierChi2Cut_;                     //!< Cut on outlier when no scaling is chosen
  std::vector<int> fixedJetPars_;             //!< List of fixed jet parameters
  std::vector<int> fixedGlobalJetPars_;       //!< List of fixed global jet parameters

  // LVMINI parameters
  double derivStep_;        //!< Step width for derivative calculation
  int mvec_;                //!< Number of stored vector pairs in LVMINI
  int nIter_;               //!< Number of iterations in LVMINI
  float eps_;               //!< Convergence parameter in LVMINI
  float wlf1_;              //!< Parameter 1 of strong Wolfe condition in LVMINI
  float wlf2_;              //!< Parameter 2 of strong Wolfe condition in LVMINI
  bool calcCov_;            //!< If true, calculate covariance matrix of fitted parameters
};

#endif
