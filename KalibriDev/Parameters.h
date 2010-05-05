//
// Original Authors:  Christian Autermann, Hartmut Stadie
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: Parameters.h,v 1.60 2010/04/13 13:53:21 mschrode Exp $
//
#ifndef TParameters_h
#define TParameters_h

//C++ libs
#include <vector>
#include <map>
#include <string>
#include <utility> 

#include <iostream>
#include <cmath>
#include <cstring>

#include "ConfigFile.h"
#include "Parametrization.h"
#include "Function.h"
#include "SmearFunction.h"


//!  \brief Connection between detector geometry and fit parameters,
//!         interface to response and error parametrizations
//!  \author Christian Autermann
//!  \date   Wed Jul 18 13:54:50 CEST 2007
//!  $Id: Parameters.h,v 1.60 2010/04/13 13:53:21 mschrode Exp $
// -----------------------------------------------------------------
class TParameters {  
public :
  
  static TParameters* CreateParameters(const ConfigFile& config);

  std::string GetName() const;

  int GetEtaBin(int const eta_id) const { return GetEtaBin(eta_id, eta_granularity, phi_granularity, eta_symmetry);}
  int GetPhiBin(int const phi_id) const { return GetPhiBin(phi_id, phi_granularity);}
  int GetJetEtaBin(int const eta_id) const { return GetEtaBin(eta_id, eta_granularity_jet, phi_granularity_jet, eta_symmetry);}
  int GetJetPhiBin(int const phi_id) const { return GetPhiBin(phi_id, phi_granularity_jet);}
  int GetTrackEtaBin(int const eta_id) const { return GetEtaBin(eta_id, eta_granularity_track, phi_granularity_track, eta_symmetry);}
  int GetTrackPhiBin(int const phi_id) const { return GetPhiBin(phi_id, phi_granularity_jet);}
  int GetBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity + phibin;}
  int GetJetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return etabin*phi_granularity_jet + phibin;}
  int GetTrackBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity_track + phibin;}

  int GetNumberOfTowerParameters() const{return p->nTowerPars() *eta_granularity*phi_granularity;}
  int GetNumberOfJetParameters() const{return p->nJetPars()*eta_granularity_jet*phi_granularity_jet;}
  int GetNumberOfTrackParameters() const{return p->nTrackPars()*eta_granularity_track*phi_granularity_track;}
  int GetNumberOfGlobalJetParameters() const{return p->nGlobalJetPars();}
  int GetNumberOfFixedParameters() const {
    int n = 0;
    for(std::vector<bool>::const_iterator it = isFixedPar_.begin();
	it != isFixedPar_.end(); it++) {
      if( *it ) n++;
    }
    return n;
  }

  int GetNumberOfParameters() const{return GetNumberOfTowerParameters()+GetNumberOfJetParameters() + GetNumberOfTrackParameters()+GetNumberOfGlobalJetParameters();}
  int GetNumberOfTowerParametersPerBin() const {return p->nTowerPars();}
  int GetNumberOfJetParametersPerBin() const {return p->nJetPars();}
  int GetNumberOfTrackParametersPerBin() const {return p->nTrackPars();}
  int GetNumberOfCovCoeffs() const { 
    return (GetNumberOfParameters()*GetNumberOfParameters()+GetNumberOfParameters())/2;
  }

  int GetEtaGranularity() const { return eta_granularity;}
  int GetPhiGranularity() const { return phi_granularity;}
  int GetEtaGranularityJet() const { return eta_granularity_jet;}
  int GetPhiGranularityJet() const { return phi_granularity_jet;}
  int GetEtaGranularityTrack() const { return eta_granularity_track;}
  int GetPhiGranularityTrack() const { return phi_granularity_track;}

  void writeCalibrationCfi(const char* name); //!< write calibration constants to cfi file
  void writeCalibrationTxt(const char* name); //!< write calibration constants to txt file
  void writeCalibrationTex(const char* name, const ConfigFile& config); //!< write calibration constants and some paraemters of the fit to tex file

  double* GetTowerParRef(int bin) { return k + bin*p->nTowerPars(); }
  double* GetJetParRef(int jetbin)  { return k + GetNumberOfTowerParameters()+jetbin*p->nJetPars();}
  double* GetTrackParRef(int trackbin)  { return k + GetNumberOfTowerParameters() + GetNumberOfJetParameters() +trackbin*p->nTrackPars();}
  double* GetGlobalJetParRef()  { return k + GetNumberOfTowerParameters() + GetNumberOfJetParameters() + GetNumberOfTrackParameters();}

  double* GetTowerParErrorRef(int bin) { 
    return parErrors_ + bin*p->nTowerPars();
  }
  double* GetJetParErrorRef(int jetbin)  { 
    return parErrors_ + GetNumberOfTowerParameters()+jetbin*p->nJetPars();
  }
  double* GetTrackParErrorRef(int trackbin)  {
    return parErrors_ + GetNumberOfTowerParameters() + GetNumberOfJetParameters() +trackbin*p->nTrackPars();
  }
  double* GetGlobalJetParErrorRef()  { 
    return parErrors_ + GetNumberOfTowerParameters() + GetNumberOfJetParameters() + GetNumberOfTrackParameters();
  }

  double* GetTowerParGlobalCorrCoeffRef(int bin) { 
    return parGCorr_ + bin*p->nTowerPars();
  }
  double* GetJetParGlobalCorrCoeffRef(int jetbin)  { 
    return parGCorr_ + GetNumberOfTowerParameters()+jetbin*p->nJetPars();
  }
  double* GetTrackParGlobalCorrCoeffRef(int trackbin)  {
    return parGCorr_ + GetNumberOfTowerParameters() + GetNumberOfJetParameters() +trackbin*p->nTrackPars();
  }
  double* GetGlobalJetParGlobalCorrCoeffRef()  { 
    return parGCorr_ + GetNumberOfTowerParameters() + GetNumberOfJetParameters() + GetNumberOfTrackParameters();
  }

  bool isFixedPar(int i) const { 
    assert( i >= 0 && i < GetNumberOfParameters() );
    return isFixedPar_[i];
  }
  std::string parName(int i) const {
    assert( i >= 0 && i < GetNumberOfParameters() );
    return parNames_[i];
  }

  void SetParameters(double *np) {
    std::memcpy(k,np,GetNumberOfParameters()*sizeof(double));
  }
  void SetErrors(double *ne) {
    std::memcpy(parErrors_,ne,GetNumberOfParameters()*sizeof(double));
  }  
  void SetGlobalCorrCoeff(double *gcc) {
    std::memcpy(parGCorr_,gcc,GetNumberOfParameters()*sizeof(double));
  }  
  void SetCovCoeff(double *cov) {
    std::memcpy(parCov_,cov,GetNumberOfCovCoeffs()*sizeof(double));
  }
  void fixPar(int i) {
    assert( i >= 0 && i < GetNumberOfParameters() );
    isFixedPar_[i] = true;
  }
  void SetFitChi2(double chi2) { fitchi2 = chi2;}
  double GetFitChi2() const { return fitchi2;}
  void FillErrors(double* copy) const {
    std::memcpy(copy,parErrors_,GetNumberOfParameters()*sizeof(double));
  }
  double* GetPars() { return k; }
  double* GetErrors() { return parErrors_; }
  double* GetGlobalCorrCoeff() { return parGCorr_; }
  double* GetCovCoeff() { return parCov_; }
  double* GetEffMap() {return trackEff;}
  int GetTrackEffBin(double pt, double eta);

  void print() const;

  const char * name() const { return p->name(); }
  bool needsUpdate() const { return p->needsUpdate(); }
  void update() { p->update(GetPars()); }
  
  static double tower_parametrization(const Measurement* x, const double* par) {
    return instance->p->correctedTowerEt(x,par);
  }
  static double jet_parametrization(const Measurement* x, const double* par) {
    return instance->p->correctedJetEt(x,par);
  }  
  static double inv_jet_parametrization(const Measurement* x, const double* par) {
    return instance->p->inverseJetCorrection(x,par);
  }  
  static double track_parametrization(const Measurement* x, const double* par) {
    return instance->p->GetExpectedResponse(x,par);
  }
  static double global_jet_parametrization(const Measurement* x, const double* par) {
    return instance->p->correctedGlobalJetEt(x,par);
  }

  static double dummy_parametrization(const Measurement* x, const double* par) {
    return x->pt;
  }

  //Error parametrization functions:
  template<int Et> static double const_error(const double *x, const Measurement *xorig=0, double errorig=0) {
    return Et;
  }
  static double tower_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0) { 
    return (x[0]>0 ?  1.25 * sqrt( x[0])   :   1.25 * sqrt(-x[0]) );
  }
  static double jet_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0) {
    return (x[0]>0. ? 0.033*x[0] + 5.6   :   0.033*(-x[0]) + 5.6 ); 
  }

  static double track_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0) { 
    //for full error also see Grooms paper 0605164v4, p.25
    double error=0,error2=0;
    error =  (x[0]>0 ? x[0] *( 0.05 + 0.00015 * x[0])   : (-x[0]) *(  0.05 + 0.00015 * (-x[0]) )); //trackerror  to be checken and dependent on pt, eta, chi2, nohits, ....
    error2 = error * error;

    //Pi0 Fehler s.Clemens
    if(x[0] > 3)      error = x[0] * 0.15 + 3;              //p. 70
    else              error = x[0] * 1.15;
    error2 += error * error;

    //error2 += (1-1/1.48)*(1-1/1.48)*0.125*0.125*x[0]*x[0];   //*(x[0]/100)^(-0.076)          //1/1.48 = h/e
    //following term has to be checked!!!!
    //double a = 1/(1.48 * 1.48) * 1.25 * 1.25 * pow((fabs(x[0])* (xorig->E / xorig->pt) / 0.96),(0.816 - 1));  // 1- Pi0 * error(h)^2 (h/e)^2
    //error2 += (x[0]>0 ?  a * x[0]  : a * (-x[0]));    //intrinsic term (HCAL)
    error = sqrt(error2);
    return error;
  }


  static double jet_only_tower_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0) { 
    return 0;
  }



  //!  \brief Parameters from V. Chetluru's fit to L2L3 corrected jets
  //!
  //!  Use results from V. Chetluru's talk:
  //!  <A HREF="http://indico.cern.ch/getFile.py/access?contribId=1&resId=1&materialId=slides&confId=52598">
  //!  Jet energy resolution studies
  //!  </A>.
  //!
  //!  The absolute resolution is given by
  //!  \f[
  //!   \sigma^{2} = a^{2} + b^{2}p_{T} + c^{2}p^{2}_{T}
  //!  \f]
  //!  with the \f$ \eta \f$ dependent parameters
  //!  <TABLE>
  //!   <TR>
  //!    <TD>  </TD>
  //!    <TD> a </TD>
  //!    <TD> b </TD>
  //!    <TD> c </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 0 < \eta < 0.8 \f$ </TD>
  //!    <TD> 4.44 </TD>
  //!    <TD> 1.11 </TD>
  //!    <TD> 0.03 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 0.8 < \eta < 1.5 \f$ </TD>
  //!    <TD> 4.35 </TD>
  //!    <TD> 1.17 </TD>
  //!    <TD> 0.04 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 1.5 < \eta < 2.4 \f$ </TD>
  //!    <TD> 4.34 </TD>
  //!    <TD> 0.85 </TD>
  //!    <TD> 0.03 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 2.4 < \eta < 3.2 \f$ </TD>
  //!    <TD> 4.08 </TD>
  //!    <TD> 0.45 </TD>
  //!    <TD> 0.04 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 3.2 < \eta \f$ </TD>
  //!    <TD> 3.90 </TD>
  //!    <TD> 0.29 </TD>
  //!    <TD> 0.09 </TD>
  //!   </TR>
  //!  </TABLE>
  //!
  //!  \return The absolute resolution
  // -----------------------------------------------------
  static double jet_only_jet_error_parametrization_et(const double *x, const Measurement *xorig=0, double errorig=0) {
    const static double a[5] = { 4.44 * 4.44, 4.35 * 4.35, 4.34 * 4.34 , 4.08 * 4.08, 3.90 * 3.90 };
    const static double b[5] = { 1.11 * 1.11, 1.17 * 1.17, 0.85 * 0.85, 0.45 * 0.45, 0.29 * 0.29};
    const static double c[5] = { 0.03 * 0.03, 0.04 * 0.04, 0.03 * 0.03, 0.04 * 0.04, 0.09 * 0.09};

    double abseta = std::abs(xorig->eta);
    int i = (abseta < 0.8) ? 0 : ((abseta < 1.5) ? 1 : ((abseta < 2.4) ? 2 : (abseta < 3.2) ? 3 : 4));
    return sqrt(a[i] + (b[i] + c[i] *x[0]) * x[0]);
  }

  static double jet_only_jet_error_parametrization_energy(const double *x, const Measurement *xorig=0, double errorig=0) {
    /*
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    //constant before stochastic term is not properly knowen
    return (x[0]>0. ? 0.033*x[0] + 5.6 + 1.0 * sqrt(pmess)  :   0.033*(-x[0]) + 5.6 + 1.0 * sqrt(-pmess) ); 
    */
    double E = x[0] * xorig->E/xorig->pt;
    //double sqE = sqrt(E);
    return sqrt(1.3*1.3/E + 0.056 * 0.056) * x[0];
  }


  static double dummy_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0) {        
    return x[0];  
  }
  static double fast_error_parametrization(const double *x, const Measurement *xorig, double errorig)  {
    return (xorig->pt==0. ? errorig : errorig*x[0]/xorig->pt );  
  }
  static double jans_E_tower_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0)  {
    
    // E = x[0]*xorig[7];  x[0]=param. mess;    xorig == _mess
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    return (xorig->E!=0. ? tower_error_parametrization(&pmess,xorig,errorig) * xorig->pt / xorig->E : 0.0);
    
    //return 0;
  }

  static double toy_tower_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0);
  
  static double toy_jet_error_parametrization(const double *x, const Measurement *xorig=0, double errorig=0);

  static double const_error_parametrization(const double *x, const Measurement *xorig, double errorig)  {
    return errorig;  
  }

  //Plot paramterization stuff
  static double plot_parametrization(const Measurement* x, const double* par) {
    return tower_parametrization(x,par)/x->pt; 
  }

  static double jes_plot_parametrization(double * x,double * par)  {
    //return jet_parametrization(x,par)/x->pt;
    return ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x[0]));   

  }

  //Limiting parameters
  static double parameter_limit(const Measurement* x, const double *par) {
    double min = x->pt;
    double max = x->EMF;
    if(par[0] < min) return (min-par[0]);
    if(par[0] > max) return (par[0]-max);
    return 0;
    //return 1e4/(1+exp(k* (par[0] - min))) + 1e4/(1+exp(-k* (par[0] - max));
  }

  float etaEdge(int const etaBin, bool lowerEdge);
  //! return upper edge of bin in eta
  float etaUpperEdge(int const etaBin) { return etaEdge(etaBin, false); };
  //! return lower edge of bin in eta
  float etaLowerEdge(int const etaBin) { return etaEdge(etaBin, true ); };

  Function tower_function(int etaid, int phiid);
  Function jet_function(int etaid, int phiid);
  Function track_function(int etaid, int phiid);
  Function global_jet_function();
  SmearFunction resolutionFitPDF(int etaid, int phiid);

  void readCalibrationCfi(const std::string& file);
  void readCalibrationTxt(const std::string& file);
  void readCalibrationJetMET(const std::vector<std::string>& inputFileNames);
  void readCalibrationJetMETL2(const std::string& inputFileName);
  void readCalibrationJetMETL3(const std::string& inputFileName);


protected:
  TParameters(Parametrization* p) 
    : p(p),k(0),parErrors_(0),parGCorr_(0),parCov_(0),trackEff(0),fitchi2(0) {
  };
  virtual ~TParameters();


private:
  TParameters();
  TParameters(const TParameters&) {}
  int GetEtaBin(int phi_id, int etagranu, int phigranu, bool etasym) const;
  int GetPhiBin(int phi_id, int phigranu) const;
  //! Return one line of LaTeX tabular containing the name and value of a given parameter from config file
  template<class T> std::string texTabularLine(const ConfigFile& config, const std::string& fieldname) const;
  //! Return submatrix of covariance matrix for \p nPar parameters from \p firstPar
  std::vector<int> findCovIndices(int firstPar, int nPar) const;
  //! Return stati (is fixed?) for \p nPar parameters from \p firstPar
  std::vector<bool> findParStatus(int firstPar, int nPar) const;

  //Towers in Eta-, Phi- direction (according to PTDR Vol I, p.201)
  unsigned const static eta_ntwr=82, phi_ntwr=72;
  unsigned eta_ntwr_used;
  bool eta_symmetry;
  unsigned int eta_granularity, phi_granularity,eta_granularity_jet, phi_granularity_jet, eta_granularity_track, phi_granularity_track;
  std::vector<double> start_values, jet_start_values, track_start_values, global_jet_start_values;
  std::vector<std::string> parNames_;

  //The parametrization functions:
  Parametrization* p;

  double * k; //!< all fit-parameters
  std::vector<bool> isFixedPar_;
  double * parErrors_; //!< all fit-parameter errors
  double * parGCorr_; //!< Global correlation coefficients of parameters
  double * parCov_;
  double * trackEff; //!< track Efficiency 13eta X 13 ptbins;
  double fitchi2;

  /// ------------------------------------------------------
  /// private functions

  void Init(const ConfigFile& config);
  void readTrackEffTxt(const std::string& file);
  std::string trim(std::string const& source, char const* delims = " {}\t\r\n");

  static TParameters *instance; 

  static Parametrization* CreateParametrization(const std::string& name, const ConfigFile& config);
  
  class Cleaner
  {
  public:
    Cleaner() {}
    ~Cleaner()
    {
      if(TParameters::instance) { 
	delete TParameters::instance; 
	TParameters::instance = 0; 
      }
    }
  };
  friend class Cleaner;
};

#endif
