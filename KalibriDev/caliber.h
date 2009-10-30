//  $Id: caliber.h,v 1.50 2009/08/07 12:20:50 mschrode Exp $

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
//!  \section label_sec_milestones Milestones
//!  - July 2009: Determination of the Summer08 L2L3 MC truth corrections
//!    - <A HREF="http://indico.cern.ch/conferenceDisplay.py?confId=64137">
//!      Presentation in Jet Energy Corrections Meeting</A>
//!    - <A HREF="../results/L2L3MCTruthCorrections/KalibriL2L3fromSummer08DiJetMC.txt">
//!      Calibration constants</A>
//!    - Controlplots
//!      - <A HREF="../results/L2L3MCTruthCorrections/controlplots.root">
//!        controlplots.root</A>
//!      - <A HREF="../results/L2L3MCTruthCorrections/controlplotsJetTruthEventResponse.ps">
//!        controlplotsJetTruthEventResponse.ps</A>
//!      - <A HREF="../results/L2L3MCTruthCorrections/controlplotsJetTruthEventResolution.ps">
//!        controlplotsJetTruthEventResolution.ps</A>
//!      - <A HREF="../results/L2L3MCTruthCorrections/controlplotsBinnedResponse.ps">
//!        controlplotsBinnedResponse.ps</A>
//!    - <A HREF="../results/L2L3MCTruthCorrections/L2L3Summer08Dijets.cfg">
//!      Configuration file</A>
//!
//!  \section label_sec_geninfo General information about and results of the calibration method
//!
//!  - J. Naumann-Emme:
//!    <A HREF="http://indico.cern.ch/conferenceDisplay.py?confId=65022">
//!    JEC from top using "Kalibri" method</A>,
//!    CMS Jet Energy Corrections Meeting, 24th July, 2009
//!
//!  - M. Schr&ouml;der:
//!    <A HREF="http://indico.cern.ch/conferenceDisplay.py?confId=64137">
//!    An Unbinned Fit for Jet Energy Corrections</A>,
//!    CMS Jet Energy Corrections Meeting, 10th July, 2009
//!
//!  - J. Thomsen:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=1638">
//!    Inclusion of track information</A>,
//!    UHH Jetcalibration meeting, 8th May, 2009
//!
//!  - H. Stadie:
//!    <A HREF="http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=55161">
//!    Jet Calibraton Activities in Hamburg</A>,
//!    FSP-CMS Meeting, 28th April 2009
//!
//!  - S. Naumann-Emme:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=3&amp;resId=0&amp;materialId=slides&amp;confId=1683">
//!    Jet Energy Corrections from Top Quark Decays</A>,
//!    CMS Hamburg Meeting, 28th January 2009
//!
//!  - M. Schr&ouml;der:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=2&amp;resId=0&amp;materialId=slides&amp;confId=1683">
//!    Conceptual Studies for a Jet Energy Correction</A>,
//!    CMS Hamburg Meeting, 28th January 2009
//!
//!  - H. Stadie:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&amp;resId=0&amp;materialId=slides&amp;confId=1683">
//!    CMS Calorimeter and Jet Calibration</A>,
//!    CMS Hamburg Meeting, 28th January 2009
//!
//!  - S. Naumann-Emme: Top as a Calibration Tool,
//!    2nd "Physics at the Terascale" Workshop, 27th November 2008
//!
//!  - T. Schum:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&amp;resId=0&amp;materialId=slides&amp;confId=627">
//!    HCAL Calibration using a Global Fit Ansatz</A>,
//!    Hamburg CMS Meeting, 12th March 2008
//!
//!  - R. Wolf:
//!    <A HREF="http://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=29582">
//!    Data-Driven Calorimeter Calibration Exploiting a Global-Fit Ansatz</A>,
//!    JetMET Meeting, 26th February 2008
//!
//!  - C. Autermann:
//!    <A HREF="http://indico.cern.ch/getFile.py/access?contribId=7&resId=0&materialId=slides&confId=22705">
//!    A global fit approach to HCAL/jet calibration</A>,
//!    JetMET Meeting, 18th October 2007
//!
//!  \section label_sec_techinfo Technical information and tutorials
//!  - <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=1643">
//!    SmearData: Extension of framework to jet smearing method</A>
//!    UHH Jetcalibration meeting, 12th June, 2009
//!
//!  - <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=1623">
//!    Inversion technique including towers</A>
//!    UHH Jetcalibration meeting, 23rd January, 2009
//!    
//!  - <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=1622">
//!    Details of the inversion technique</A>
//!    UHH Jetcalibration meeting, 16th January, 2009
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
class TControlPlots;
class TData;
class TMeasurement;


//!  \brief Main program
//!  \note  For profiling:
//!         To prevent gprof from missing the threads: 
//!         wget http://sam.zoy.org/writings/programming/gprof-helper.c
//!         gcc -shared -fPIC gprof-helper.c -o gprof-helper.so -lpthread -ldl 
//!         LD_PRELOAD=./gprof-helper.so ./junk
//!  \author Christian Autermann
//!  \date Wed Jul 18 13:54:50 CEST 2007
//!  $Id: caliber.h,v 1.50 2009/08/07 12:20:50 mschrode Exp $
// -----------------------------------------------------------------
class TCaliber {
public :
  //!  \brief Constructor
  //!  \param f Name of the configuration file
  // -----------------------------------------------------------------
  TCaliber(const std::string& f)
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
  wlf2_(0.9)
  {};

  ~TCaliber(){};

  void init();   //!< Read parameters from configfile, read data
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
  std::vector<TData*> data_;  //!< The data
  int nGammajetEvents_;       //!< Number of gamma-jet events
  int nDijetEvents_;          //!< Number of dijet events
  int nTrijetEvents_;         //!< Number of trijet events
  int nTrackClusterEvents_;   //!< Number of track-cluster events
  int nZjetEvents_;           //!< Number of Zjet events
  int nTopEvents_;            //!< Number of top events

  // control parameters of fit
  bool printParNDeriv_;     //!< Control whether to print derivatives in each iteration
  std::vector<int> residualScalingScheme_;    //!< Iteration scheme of scaling of residuals
  double outlierChi2Cut_;                     //!< Cut on outlier when no scaling is chosen
  std::vector<int> globalJetPars_;            //!< List of global jet parameters
  std::vector<int> fixedJetPars_;             //!< List of fixed jet parameters
  std::vector<int> fixedGlobalJetPars_;       //!< List of fixed global jet parameters

  // LVMINI parameters
  double derivStep_;        //!< Step width for derivative calculation
  int mvec_;                //!< Number of stored vector pairs in LVMINI
  int nIter_;               //!< Number of iterations in LVMINI
  float eps_;               //!< Convergence parameter in LVMINI
  float wlf1_;              //!< Parameter 1 of strong Wolfe condition in LVMINI
  float wlf2_;              //!< Parameter 2 of strong Wolfe condition in LVMINI
};

#endif
