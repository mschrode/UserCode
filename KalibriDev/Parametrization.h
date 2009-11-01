//  $Id: Parametrization.h,v 1.2 2009/10/30 13:14:55 mschrode Exp $

#ifndef CALIBCORE_PARAMETRIZATION_H
#define CALIBCORE_PARAMETRIZATION_H

#include <cmath>
#include <vector>

#include "TH1D.h"
#include "TMath.h"

#include "CalibData.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
     


//!  \brief Abstract base class for parametrizations of
//!         correction functions
//!
//!  Interface to different parametrizations for the 
//!  underlying hypothesis of the GlobalFit. Allows
//!  to correct a tower or jet measurement.
//!  \author Hartmut Stadie
//!  \date Thu Apr 03 17:09:50 CEST 2008
//!  $Id: Parametrization.h,v 1.2 2009/10/30 13:14:55 mschrode Exp $
// -----------------------------------------------------------------
class Parametrization 
{
public:
  //!  \brief Constructor
  //!  \param ntowerpars Number of parameters for tower parametrization
  //!  \param njetpars Number of parameters for (eta,phi)-dependent jet parametrization
  //!  \param ntrackpars Number of parameters for track parametrization
  //!  \param nglobaljetpars Number of parameters for global jet parametrization
  Parametrization(unsigned int ntowerpars, unsigned int njetpars, 
		  unsigned int ntrackpars, unsigned int nglobaljetpars) : 
    ntowerpars_(ntowerpars), njetpars_(njetpars), ntrackpars_(ntrackpars), 
    nglobaljetpars_(nglobaljetpars) {}

  virtual ~Parametrization() {}


  //!  \brief Corrects the measured tower Et
  //!
  //!  The parameters of the tower correction function
  //!  may depend on the tower Et, eta, and phi.
  //!
  //!  \param x Tower Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also TMeasurement):
  //!           - x->pt   : et  of whole tower
  //!           - x->EMF  : et  of ECAL  part
  //!           - x->HadF : et  of HCAL  part
  //!           - x->OutF : et  of Outer part
  //!           - x->E    : en  of Outer part
  //!           - x->eta  : eta of tower
  //!           - x->phi  : phi of tower
  //!  \param par Parameters of the correction function of this tower
  //!  \return The corrected Et of a tower 
  // -----------------------------------------------------------------
  virtual double correctedTowerEt(const TMeasurement *x,const double *par) const = 0;


  //!  \brief Corrects the measured jet Et
  //!
  //!  The parameters of the jet correction function
  //!  may depend on the jet Et, eta, and phi.
  //!
  //!  \param x Jet Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also TMeasurement):
  //!           - x->pt   : et  of whole jet
  //!           - x->EMF  : et  of ECAL  part
  //!           - x->HadF : et  of HCAL  part
  //!           - x->OutF : et  of Outer part
  //!           - x->E    : en  of Outer part
  //!           - x->eta  : eta of jet
  //!           - x->phi  : phi of jet
  //!  \param par Parameters of the correction function of this jet
  //!  \return The corrected Et of a jet
  // -----------------------------------------------------------------
  virtual double correctedJetEt(const TMeasurement *x,const double *par) const = 0;

 
  //!  \brief Returns the expected signal of a track in the Calorimeter
  //!
  //!  \param x Track Et measurement from which the expected response
  //!           of the calorimeter is calculated (see also TTrack).
  //!  \param par Parameters of the response function of this track
  //!  \return The expected calorimeter response
  // -----------------------------------------------------------------
  virtual double GetExpectedResponse(const TMeasurement *x,const double *par) const { return x->pt;}


  //!  \brief Corrects the measured jet Et with global
  //!         correction function
  //!
  //!  The parameters of the global jet correction function
  //!  are independent of the jet Et, eta, and phi.
  //!
  //!  \param x Jet Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also TMeasurement):
  //!           - x->pt   : et  of whole jet
  //!           - x->EMF  : et  of ECAL  part
  //!           - x->HadF : et  of HCAL  part
  //!           - x->OutF : et  of Outer part
  //!           - x->E    : en  of Outer part
  //!           - x->eta  : eta of jet
  //!           - x->phi  : phi of jet
  //!  \param par Parameters of the global correction function
  //!  \return The corrected Et of a jet
  // -----------------------------------------------------------------
  virtual double correctedGlobalJetEt(const TMeasurement *x,const double *par) const { return x->pt;}


  //!  \brief Get the name of the parametrization class
  //!  \return Name of the parametrization class
  // -----------------------------------------------------------------
  virtual const char * name() const = 0;


  //!  \brief Get the number of parameters of the tower parametrization
  //!  \return Number of parameters of the tower parametrization
  // -----------------------------------------------------------------
  unsigned int nTowerPars() const { return ntowerpars_;}


  //!  \brief Get the number of parameters of the jet parametrization
  //!  \return Number of parameters of the jet parametrization
  // -----------------------------------------------------------------
  unsigned int nJetPars() const { return njetpars_;}


  //!  \brief Get the number of parameters of the track parametrization
  //!  \return Number of parameters of the track parametrization
  // -----------------------------------------------------------------
  unsigned int nTrackPars() const { return ntrackpars_;}


  //!  \brief Get the number of parameters of the global jet parametrization
  //!  \return Number of parameters of the global jet parametrization
  // -----------------------------------------------------------------
  unsigned int nGlobalJetPars() const { return nglobaljetpars_;}

  //!  \brief Returns the expected response for a jet with true x-> Et 
  //!         (this is the inverted jet correction
  //!
  //!  \param x TMeasurement describing the true jet properties
  //!  \param par Parameters of the jet response function
  //!  \return the expected calorimeter response
  // -----------------------------------------------------------------
  virtual double inverseJetCorrection(const TMeasurement *x,const double *par) const { return -1;}

  
  //!  \brief Get the number of parameters of the track p
  //!  \return Number of parameters of the track parametrization
  // -----------------------------------------------------------------
  virtual bool hasInvertedCorrection() const { return false;}

  virtual bool needsUpdate() const { return false; }
  virtual void update(const double * par) {;}

private: 
  Parametrization();
  unsigned int ntowerpars_;      //!< Number of parameters of the tower parametrization
  unsigned int njetpars_;        //!< Number of parameters of the jet parametrization
  unsigned int ntrackpars_;      //!< Number of parameters of the track parametrization
  unsigned int nglobaljetpars_;  //!< Number of parameters of the global jet parametrization
};



//!  \brief Parametrization that does not change a thing  
//!
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class ConstParametrization : public Parametrization { 
public:
  ConstParametrization() : Parametrization(0,0,0,0) {}
  const char* name() const { return "ConstParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return  x->pt;   
  }
};

//!  \brief Parametrization of the hadronic tower response 
//!         by a step function in Et
//!
//!  The total jet measurement remains unchanged.
//!
//!  This parametrization has 12 tower parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepParametrization : public Parametrization { 
public:
  StepParametrization() : Parametrization(12,0,0,0) {}
  const char* name() const { return "StepParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result = 0;
    
    if(x->HadF>=0.0  && x->HadF<=1.0)  result = x->EMF+x->OutF + par[0]*x->HadF; 
    else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[11]*x->HadF;
    return result;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return  x->pt;   
    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    //return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt)); 
    /*
      double result = 0;
      if(x->pt>=0.0  && x->pt<=5.0)          result =  par[0]*x->pt + par[1];
      else if (x->pt>5.0   && x->pt<=20.0)   result =  par[2]*x->pt + par[3];
      else if (x->pt>20.0  && x->pt<=80.0)   result =  par[4]*x->pt + par[5];
      else if (x->pt>80.0 )                 result =  par[6]*x->pt + par[7];
      return result;
    */
  }
};



//!  \brief Parametrization of the hadronic tower response 
//!         by a step function in E
//!
//!  Additionally, the total jet measurement is corrected
//!  by an Et-dependent function.
//!
//!  This parametrization has 12 tower parameters and 2
//!  jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepParametrizationEnergy : public Parametrization { 
public:
  StepParametrizationEnergy() : Parametrization(12,2,0,0) {}
  const char* name() const { return "StepParametrizationEnergy";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result = 0;
    // reweight from et to en
    double e =  x->HadF * x->E / x->pt;
    
    if(e>=0.0  && e<=1.0)  result = x->EMF+x->OutF + par[0]*x->HadF;
    else if (e>   1.0 && e<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if (e>   2.0 && e<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if (e>   5.0 && e<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if (e>  10.0 && e<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if (e>  20.0 && e<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if (e>  40.0 && e<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if (e>  80.0 && e<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if (e> 160.0 && e<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if (e> 300.0 && e<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if (e> 600.0 && e<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if (e>1000.0 )             result = x->EMF+x->OutF+par[11]*x->HadF;
    
    return result;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   
  }
};



//!  \brief Parametrization of the hadronic tower
//!         response by a step function in EMF
//!
//!  Parametrization of hadronic response of a tower
//!  by a step function with 3 sets of parametrizations for 
//!  different em fractions.
//!
//!  This parametrization has 36 tower parameters.
//!
//!  The total jet measurement remains unchanged.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepEfracParametrization : public Parametrization {
public:
  StepEfracParametrization() : Parametrization(36,0,0,0) {}  //(36,2) {}
  const char* name() const { return "StepEfracParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result=0;
    
    double Efrac = x->EMF/(x->HadF+x->OutF+x->EMF);
    if( Efrac < 0.2 ) {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[ 0]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[11]*x->HadF;
    } else if (Efrac < 0.5) {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[12]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[13]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[14]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[15]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[16]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[17]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[18]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[19]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[20]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[21]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[22]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[23]*x->HadF;
    } else {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[24]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[25]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[26]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[27]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[28]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[29]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[30]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[31]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[32]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[33]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[34]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[35]*x->HadF;
    }
    return result;
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return  x->pt;   
    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    //return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   
  }
};



//!  \brief Parametrization of the hadronic jet response 
//!         by a step function in EMF
//!
//!  Parametrization of hadronic response of a jet by a step 
//!  function with 3 sets of parametrizations for 
//!  different em fractions; the tower response remains
//!  unchanged.
//!
//!  This parametrization has 65 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepJetParametrization : public Parametrization { 
public:
  StepJetParametrization() : Parametrization(0,65,0,0) {}
  const char* name() const { return "StepJetParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double pt     = x->pt;
    double Efrac  = x->EMF/( x->EMF+x->HadF+x->OutF);
    double result = 0.;
    if(Efrac < 0.2) {
      if(pt>=0.0 && pt<=10.0)         result = x->EMF+x->OutF+par[ 0]*x->HadF;
      else if (pt> 10.0 && pt<= 20.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
      else if (pt> 20.0 && pt<= 30.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
      else if (pt> 30.0 && pt<= 40.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
      else if (pt> 40.0 && pt<= 60.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
      else if (pt> 60.0 && pt<= 80.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
      else if (pt> 80.0 && pt<=100.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
      else if (pt>100.0 && pt<=120.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
      else if (pt>120.0 && pt<=140.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
      else if (pt>140.0 && pt<=160.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF+par[10]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF+par[11]*x->HadF;
      else if (pt>200.0 && pt<=225.0) result = x->EMF+x->OutF+par[12]*x->HadF;
      else if (pt>225.0 && pt<=250.0) result = x->EMF+x->OutF+par[13]*x->HadF;
      else if (pt>250.0 && pt<=275.0) result = x->EMF+x->OutF+par[14]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF+par[15]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF+par[16]*x->HadF;
      else if (pt>350.0 && pt<=400.0) result = x->EMF+x->OutF+par[17]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF+par[18]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF+par[19]*x->HadF;
      else if (pt>700.0 )             result = x->EMF+x->OutF+par[20]*x->HadF;
    } else if (Efrac < 0.5) {
      if(pt>=0.0 && pt<=10.0)         result = x->EMF+x->OutF+par[21]*x->HadF;
      else if (pt> 10.0 && pt<= 20.0) result = x->EMF+x->OutF+par[22]*x->HadF;
      else if (pt> 20.0 && pt<= 30.0) result = x->EMF+x->OutF+par[23]*x->HadF;
      else if (pt> 30.0 && pt<= 40.0) result = x->EMF+x->OutF+par[24]*x->HadF;
      else if (pt> 40.0 && pt<= 60.0) result = x->EMF+x->OutF+par[25]*x->HadF;
      else if (pt> 60.0 && pt<= 80.0) result = x->EMF+x->OutF+par[26]*x->HadF;
      else if (pt> 80.0 && pt<=100.0) result = x->EMF+x->OutF+par[27]*x->HadF;
      else if (pt>100.0 && pt<=120.0) result = x->EMF+x->OutF+par[28]*x->HadF;
      else if (pt>120.0 && pt<=140.0) result = x->EMF+x->OutF+par[29]*x->HadF;
      else if (pt>140.0 && pt<=160.0) result = x->EMF+x->OutF+par[30]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF+par[31]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF+par[32]*x->HadF;
      else if (pt>200.0 && pt<=225.0) result = x->EMF+x->OutF+par[33]*x->HadF;
      else if (pt>225.0 && pt<=250.0) result = x->EMF+x->OutF+par[34]*x->HadF;
      else if (pt>250.0 && pt<=275.0) result = x->EMF+x->OutF+par[35]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF+par[36]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF+par[37]*x->HadF;
      else if (pt>350.0 && pt<=400.0) result = x->EMF+x->OutF+par[38]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF+par[39]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF+par[40]*x->HadF;
      else if (pt>700.0 )             result = x->EMF+x->OutF+par[41]*x->HadF;
    } else {
      if(pt>=0.0 && pt<=10.0)   result = x->EMF+x->OutF+par[42]*x->HadF;
      else if (pt> 10.0 && pt<= 20.0) result = x->EMF+x->OutF+par[43]*x->HadF;
      else if (pt> 20.0 && pt<= 30.0) result = x->EMF+x->OutF+par[44]*x->HadF;
      else if (pt> 30.0 && pt<= 40.0) result = x->EMF+x->OutF+par[45]*x->HadF;
      else if (pt> 40.0 && pt<= 60.0) result = x->EMF+x->OutF+par[46]*x->HadF;
      else if (pt> 60.0 && pt<= 80.0) result = x->EMF+x->OutF+par[47]*x->HadF;
      else if (pt> 80.0 && pt<=100.0) result = x->EMF+x->OutF+par[48]*x->HadF;
      else if (pt>100.0 && pt<=120.0) result = x->EMF+x->OutF+par[49]*x->HadF;
      else if (pt>120.0 && pt<=140.0) result = x->EMF+x->OutF+par[50]*x->HadF;
      else if (pt>140.0 && pt<=160.0) result = x->EMF+x->OutF+par[51]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF+par[52]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF+par[53]*x->HadF;
      else if (pt>200.0 && pt<=225.0) result = x->EMF+x->OutF+par[54]*x->HadF;
      else if (pt>225.0 && pt<=250.0) result = x->EMF+x->OutF+par[55]*x->HadF;
      else if (pt>250.0 && pt<=275.0) result = x->EMF+x->OutF+par[56]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF+par[57]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF+par[58]*x->HadF;
      else if (pt>350.0 && pt<=400.0) result = x->EMF+x->OutF+par[59]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF+par[60]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF+par[61]*x->HadF;
      else if (pt>700.0 )             result = x->EMF+x->OutF+par[62]*x->HadF;
    }

    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    result *= ( 1. + 0.295 * par[63] * exp(- 0.02566 * par[64] * result));   

    return  result;
  }
};



//!  \brief Parametrization of tower and jet response
//!         by some "clever" function
//!
//!
//!  This parametrization has 2 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class MyParametrization: public Parametrization {
 public:
  MyParametrization() : Parametrization(0,2,0,0) {}
  const char* name() const { return "MyParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->EMF + x->HadF + x->OutF;
    //return x->EMF + par[0]*x->HadF + par[1]*log(x->pt) + par[2];
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->EMF + (par[0] + x->HadF * par[1]/1000) * x->HadF + x->OutF;;
  }
};



//!  \brief Parametrization of tower and jet response
//!         with some ideas from the JetMET group
//!
//!  This parametrization has 3 tower parameters and 5
//!  jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class JetMETParametrization: public Parametrization {
public:
  JetMETParametrization() : Parametrization(3,5,0,0) {}
  const char* name() const { return "JetMETParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return par[1] * x->HadF + par[2] * x->EMF + x->OutF + par[0];
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double logx = log(x->HadF);
    if(logx < 0) logx = 0;
    return (par[0] - par[1]/(pow(logx,par[2]) + par[3]) + par[4]/x->HadF) * x->HadF;  
  }
};



//!  \brief Simple parametrization that uses a global scale factor
//!
//!  This parametrization has 0 tower parameters, 0 track parameters
//!  and 1 global jet parameter, no eta or phi dependence.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class GlobalScaleFactorParametrization: public Parametrization {
public:
  GlobalScaleFactorParametrization() : Parametrization(0,0,0,1) {}
  const char* name() const { return "GlobalScaleFactorParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedGlobalJetEt(const TMeasurement *x, const double *par) const {
    return  par[0] * x->pt;
  }

};



//!  \brief Simple tower and jet parametrization
//!
//!  This parametrization has 3 tower parameters and 3
//!  jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class SimpleParametrization: public Parametrization {
public:
  SimpleParametrization() : Parametrization(3,3,0,0) {}
  const char* name() const { return "SimpleParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return par[1] * x->EMF + par[2] * x->HadF + x->OutF + par[0];
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt * ( par[2] + par[0] * exp( -par[1] * x->pt ) );  
  }
};



//!  \brief Parametrization for toy MC with constant response
//!
//!  This is the parametrization of the correction for ToyMC
//!  events with a constant tower response i.e. when specifying
//!  only one tower constant. In this parametrization, the
//!  hadronic part of the tower Et is corrected.
//!
//!  This parametrization has 1 tower parameter.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyParametrization: public Parametrization {
public:
  ToyParametrization() : Parametrization(1,0,0,0) {}
  const char* name() const { return "ToyParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return par[0] * x->HadF + x->EMF + x->OutF;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt;  
  }
};



//!  \brief Parametrization for toy MC with constant response
//!
//!  This is the parametrization of the correction for ToyMC
//!  events with a constant tower response i.e. when specifying
//!  only one tower constant.  In this parametrization, the
//!  hadronic part of the jet Et is corrected.
//!
//!  This parametrization has 1 jet parameters.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyJetParametrization: public Parametrization {
public:
  ToyJetParametrization() : Parametrization(0,1,0,0) {}
  const char* name() const { return "ToyJetParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //return par[0] * x->HadF + x->EMF + x->OutF;
    return par[0] * x->pt;
  }
};



//!  \brief Parametrization for toy MC with step function of
//!         hadronic tower Et
//!
//!  In this parametrization, the hadronic part of
//!  tower Et is corrected by a step function. It is
//!  intended for a pt spectrum from 0 - 300 GeV
//!
//!  This parametrization has 15 tower parameters.
//!
//!  \note This parametrization is intended for studying
//!        cutoff and resolution effects on the minimization
//!        procedure.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyStepParametrization : public Parametrization { 
public:
  ToyStepParametrization() : Parametrization(15,0,0,0) {}
  const char* name() const { return "ToyStepParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result = 0.;
    double pt = x->HadF;
    
    if(pt < 2.)        result = x->EMF+x->OutF+par[ 0]*x->HadF;
    else if(pt <   5.) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if(pt <  10.) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if(pt <  20.) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if(pt <  30.) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if(pt <  40.) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if(pt <  50.) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if(pt <  60.) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if(pt <  70.) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if(pt <  80.) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if(pt <  90.) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if(pt < 100.) result = x->EMF+x->OutF+par[11]*x->HadF;
    else if(pt < 110.) result = x->EMF+x->OutF+par[12]*x->HadF;
    else if(pt < 120.) result = x->EMF+x->OutF+par[13]*x->HadF;
    else               result = x->EMF+x->OutF+par[14]*x->HadF;
    return result;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
};



//!  \brief Parametrization for toy MC with step function of
//!         hadronic jet Et
//!
//!  In this parametrization, the hadronic part of the
//!  jet Et is corrected by a step function. It is
//!  intended for a pt spectrum from 0 - 300 GeV
//!
//!  This parametrization has 15 jet parameters.
//!
//!  \note This parametrization is intended for studying
//!        cutoff and resolution effects on the minimization
//!        procedure.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyStepJetParametrization : public Parametrization { 
 public:
  ToyStepJetParametrization() : Parametrization(0,15,0,0) {}
  const char* name() const { return "ToyStepJetParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double pt = x->HadF;
    double result = 0.;
    
    if(pt < 25. )      result = x->EMF+x->OutF+par[ 0]*x->HadF;
    else if(pt <  50.) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if(pt <  75.) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if(pt < 100.) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if(pt < 125.) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if(pt < 150.) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if(pt < 175.) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if(pt < 200.) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if(pt < 225.) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if(pt < 250.) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if(pt < 275.) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if(pt < 300.) result = x->EMF+x->OutF+par[11]*x->HadF;
    else if(pt < 325.) result = x->EMF+x->OutF+par[12]*x->HadF;
    else if(pt < 350.) result = x->EMF+x->OutF+par[13]*x->HadF;
    else               result = x->EMF+x->OutF+par[14]*x->HadF;
    return  result;
  }
};



//!  \brief Complete Track Parametrization
//!
//!  Same parametrization as StepEfracParametrization,
//!  if track outside tracker or track errors too large
//!
//!  This parametrization has 12 tower parameters,
//!  3 jet parameters, and 6 track parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class TrackParametrization : public Parametrization {
 public:
  TrackParametrization() : Parametrization(12,3,6,0) {}  //(36,3,3,0) {}
  const char* name() const { return "TrackParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const 
    {
      double result=0;
      
      //double Efrac = x->EMF/(x->HadF+x->OutF+x->EMF);
      //if( Efrac < 0.2 ) {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[ 0]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[11]*x->HadF;

      if(result < x->EMF+x->OutF) result = x->EMF+x->OutF;
      return result;
    }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const 
    {
    double result;
    /*
    //result = x->pt * ( 1. + 0.295 * par[3] * exp(- 0.02566 * par[4] * x->pt));
    result = x->pt; 
    */ 
    if(x->E < -500) //set to -1000 or -800 for track jets! Positive for all others
      {
	if(x->E < 900)  //set to -1000 for Calo Rest
	  result = x->pt; //*par[0];
	else //set to -800 for complete Track Jet
	  result = x->pt * par[2];
      }
    else 
      result = par[1] * x->pt;
    if(result < 0)   result = 0; 
    return result;
  }
   
  double GetExpectedResponse(const TMeasurement *x,const double *par) const   
    {
    double result=0;
    double PiFrac;
    //Groom
    double eh = 1.48  * par[0]; 
    //Wigman
    //double eh = 1.39  * par[0]; 

    //double ehECAL = 1.6 ;//* par[1]; 
    double TrackEMF = 0;
    bool LS = false;
    
    TTrack* temp = (TTrack*)(x);
    
    if(temp->EM1 < 1.2) LS = true;   //late showering particle
    
    //this is Groom's Parametrization:
    TrackEMF = temp->EM1 / (temp->Had1 + temp->EM1);  //bei 1X1 EMF ueberschaetzt, da shower schmaler, bei 3X3 zu viele andere Tracks
    if(TrackEMF > 1) TrackEMF = 1; 
    if(TrackEMF < 0) TrackEMF = 0;
    //JPT alg (2004/15):
    //TrackEMF = 0.4;  //mean value from Test Beam

    //Groom
    PiFrac = 1 - pow(((x->E + par[4]) /(fabs(par[2]) * 0.96) ),(par[3] * 0.816) - 1 );  
    //Wigman
    //PiFrac = 0.11*par[2] * log(par[3] * x->E); 

    if(PiFrac > 1)  PiFrac = 1;   //do not allow unphysical results
    if(PiFrac < 0)  PiFrac = 0;  //happens e.g. for TrackPt<2GeV
    if(eh < 0.1) eh=0.1; 

    
    double responseH = (1 + (eh - 1) * PiFrac) / eh;
    
    /*
    //double responseE = (1 + (ehECAL - 1) * PiFrac) / ehECAL;
    double responseE = 1;
    
    double resultE = x->pt * TrackEMF * responseE;
    //if(LS) resultE = temp->EM1;
    if((temp->EM5 > 0) && (temp->EM5  < resultE) )  resultE = temp->EM5;
    
    double resultH = x->pt * (1 - TrackEMF) * responseH;
    if((temp->Had5 > 0) && (temp->Had5  < resultH) )  resultH = temp->Had5;
    result = resultE + resultH;
    */

    result = x->pt * responseH; 
 
    return result;
  }
};



//!  \brief L2L3 JetMET parametrization
//!
//!  This parametrization uses the correction functions
//!  from the JetMET group:
//!  - There is no tower correction function
//!  - The jet Et is corrected eta-dependent with the
//!    L2 correction function
//!  - The jet Et is corrected globally with the L3
//!    correction function
//!
//!  This parametrization has 6 jet parameters and
//!  4 global jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class L2L3JetParametrization : public Parametrization { 
public:
  L2L3JetParametrization() : Parametrization(0,3,0,4) {}
  const char* name() const { return "L2L3JetParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  //!  \brief Code from L2RelativeCorrector
  //!  \code
  //!  double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
  //!  double logpt = log10(pt);
  //!  double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
  //!  \endcode   
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double  pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt;
    double logpt = log10(pt);
    //    double c1 = par[0]+logpt*(0.1 * par[1]+logpt *(0.01* par[2]+logpt*(par[3]+logpt*(par[4]+logpt*par[5]))));
    //double c1 = par[0]+logpt*(0.1 * par[1]+logpt *(0.01* par[2]+logpt*(0.01*par[3]+logpt*(0.01*par[4]))));
    double c1 = par[0]+logpt*(0.1 * par[1]+logpt * 0.01* par[2]);
    return c1 * x->pt;
  }

  //!  \brief Code from SimpleL3AbsoluteCorrector
  //!  \code
  //!  double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
  //!  double log10pt = log10(pt);
  //!  double result = p[2]+p[3]/(pow(log10pt,p[4])+p[5]);
  //!  \endcode
  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt;
    double logpt = log10(pt);
    double c2 = par[0] + par[1]/(pow(logpt,par[2]) + par[3]);
    return  c2 * x->pt;
  }
};



//!  \brief L2L3 JetMET parametrization
//!
//!  This parametrization uses the correction functions
//!  from the JetMET group:
//!  - There is no tower correction function
//!  - The jet Et is corrected eta-dependent with the
//!    product of the L2 and L3 correction functions.
//!
//!  This parametrization has 7 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class L2L3JetParametrization2 : public Parametrization { 
public:
  L2L3JetParametrization2() : Parametrization(0,7,0,0) {}
  const char* name() const { return "L2L3JetParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //code from L2RelativeCorrector
    //double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
    //double logpt = log10(pt);
    //double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //double result = par[0]+logpt*(par[1]+logpt*(par[2]+logpt*(par[3]+logpt*(par[4]+logpt*par[5]))));
    double c1 = par[0]+logpt*(par[1]*0.01+logpt*0.001*par[2]);
    pt = (c1 * x->pt < 4.0) ? 4.0 : (c1 * x->pt > 2000.0) ? 2000.0 : c1 * x->pt; 
    logpt = log10(pt);
    //result = par[6] + par[7]/(pow(logpt,par[8]) + par[9]);
    double c2 = par[3] + par[4]/(pow(logpt,par[5]) + par[6]);
    //std::cout << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << ", " 
    //	      <<  par[4] << " c2:" << c2 << " Et:" << x->pt << '\n';
    return  c2 * c1 * x->pt; 
  }
};



//!  \brief Complete Track Parametrization
//!
//!  Same parametrization as L2L3JetParametrization,
//!  if track outside tracker or track errors too large
//!
//!  This parametrization has 3 jet parameters, 5 track
//!  parameters, and 4 global jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------

class L2L3JetTrackParametrization : public Parametrization { 
public:
  L2L3JetTrackParametrization() : Parametrization(0,3,5,4) {}
  const char* name() const { return "L2L3JetTrackParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //code from L2RelativeCorrector
    //double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
    //double logpt = log10(pt);
    //double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //double result = par[0]+logpt*(par[1]+logpt*(par[2]+logpt*(par[3]+logpt*(par[4]+logpt*par[5]))));
    double c1 = par[0]+logpt*(0.1 * par[1]+logpt * 0.1* par[2]);
    return c1 * x->pt;
  }

  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    // code from SimpleL3AbsoluteCorrector
    //double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
    //double log10pt = log10(pt);
    //double result = p[2]+p[3]/(pow(log10pt,p[4])+p[5]);
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //result = par[6] + par[7]/(pow(logpt,par[8]) + par[9]);
    double c2 = par[0] + par[1]/(pow(logpt,par[2]) + par[3]);
    if(c2 < par[0]) c2 = par[0];
    //std::cout << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << ", " 
    //	      << " c2:" << c2 << " Et:" << x->pt << '\n';
    return  c2 * x->pt; 
  }
  
  double GetExpectedResponse(const TMeasurement *x,const double *par) const {
    double result=0;
    //Groom
    double eh = 1.48 * par[0]; 
    //Wigman
    //double eh = 1.39  * par[0]; 

    //double ehECAL = 1.6 ;//* par[1]; 
    double TrackEMF = 0;
    bool LS = false;
    
    TTrack* temp = (TTrack*)(x);
    
    if(temp->EM1 < 1.2) LS = true;   //late showering particle
    
    //this is Groom's Parametrization:
    TrackEMF = temp->EM1 / (temp->Had1 + temp->EM1);  //bei 1X1 EMF ueberschaetzt, da shower schmaler, bei 3X3 zu viele andere Tracks
    if(TrackEMF > 1) TrackEMF = 1; 
    if(TrackEMF < 0) TrackEMF = 0;
    //JPT alg (2004/15):
    //TrackEMF = 0.4;  //mean value from Test Beam
    
    //Groom
    double Ecor = (x->E + par[4] * 0.01) /par[2] /0.96;
    double PiFrac = Ecor > 0 ? 1 - pow(Ecor,par[3] * 0.816 - 1 ) : 1.0;  
    //PiFrac = 1 - pow(Ecor,-0.16);
    //Wigman
    //PiFrac = 0.11*par[2] * log(par[3] * x->E); 
    assert(PiFrac == PiFrac);
    if(PiFrac > 1)  PiFrac = 1;   //do not allow unphysical results
    if(PiFrac < 0)  PiFrac = 0;  //happens e.g. for TrackPt<2GeV
    if(eh < 0.1) eh=0.1; 
    
    
    double responseH = (1 + (eh - 1) * PiFrac) / eh;
    
    /*
    //double responseE = (1 + (ehECAL - 1) * PiFrac) / ehECAL;
    double responseE = 1;
    
    double resultE = x->pt * TrackEMF * responseE;
    //if(LS) resultE = temp->EM1;
    if((temp->EM5 > 0) && (temp->EM5  < resultE) )  resultE = temp->EM5;
    
    double resultH = x->pt * (1 - TrackEMF) * responseH;
    if((temp->Had5 > 0) && (temp->Had5  < resultH) )  resultH = temp->Had5;
    result = resultE + resultH;
    */

    result = x->pt * responseH; 

    return result;
  }
};


//!  \brief Parametrization for Toy MC with ToyMC::Response SimpleInverse
//!
//!  This parametrization contains the correction function
//!  for the ToyMC::Response SimpleInverse
//!  \f[ C(E_{T}) = \sqrt{ A^{2} + A_{1}E_{T}} - A \f]
//!  where \f$ A = -0.5(A_{1} - A_{0} - E_{T}) \f$.
//!  This function was derived by analytical inversion.
//!  - There is no tower correction function
//!  - There is no jet Et correction function
//!  - The hadronic part of the jet Et is corrected
//!    globally with the correction function given above
//!
//!  This parametrization 2 global jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class ToySimpleInverseParametrization : public Parametrization { 
public:
  ToySimpleInverseParametrization() : Parametrization(0,1,0,2) {}
  const char* name() const { return "ToySimpleInverseParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return par[0]*x->pt;
  }

  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    double a = 0.5*(par[1] - par[0] - x->pt);
    return sqrt( a*a + par[1]*x->pt ) - a;
  }
};



//!  \brief Parametrization used for SmearFunction estimation
// -----------------------------------------------------------------
class SmearFermiTail : public Parametrization{
 public:
  SmearFermiTail() : Parametrization(0,3,0,0) {}
  const char * name() const { return "SmearFermiTail"; }

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return 0.;
  }

  //!  \brief Returns probability density of response
  //!  \param x   TMeasurement::E is response for which the 
  //!             probability density is returned
  //!  \param par Pointer to parameters
  //!  \return Probability density of response
  // ------------------------------------------------------------------------
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double c      = par[0];
    double mu     = 1.;
    double sigma  = par[1];
    double T      = par[2];

    if(c < 0.) c = 0.;
    if(c > 1.) c = 1.;
    if(T < 0.) T = 0.;

    double p = c / sqrt(2.* M_PI ) / sigma * exp(-pow((x->E - mu) / sigma, 2.) / 2.);
    p += (1. - c) / (T * log(1 + exp(mu / T))) / (exp((x->E - mu) / T) + 1.);

    return p;
  }

  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    return 0.;
  }
};



//!  \brief Parametrization used for SmearFunction estimation
//!
//!  This pdf consists of two Gaussians, a central one around
//!  1. and a second one describing the tail:
//!  \f[
//!   p(r) = c_{i} \cdot G_{0}(\mu_{0},\sigma_{0}) + (1-c_{i}) \cdot G_{1}(\mu_{1},\sigma_{1})
//!  \f]
//!  The mean of the Gaussian is parameterized by
//!  \f[
//!   \mu(t) = a_{1} + a_{2}\cdot t
//!  \f]
//!  where the \f$ a_{i} \f$ are specified via the config file.
//!  The width of the Gaussian is parameterized by
//!  \f[
//!   \frac{\sigma(t)}{p_{T}} = \sqrt{\frac{a^{2}_{1}}{t} + \frac{a^{2}_{2}}{t} + a^{2}_{3}},
//!  \f]
//!  where the \f$ a_{i} \f$ are parameters of the fit.
//!
//!  The parameters of the tail Gaussian are binned
//!  in \f$ N_{i} \f$ \p t bins.
//!
//!  The truth pdf is \f$ f(t) = N/t^{n} \f$. The normalization
//!  is such that the dijet probability
//!  \f[
//!    p = \int\;dt\;f(t)r(m_{1}/t)r(m_{2}/t)
//!  \f]
//!  is normalized to one. Therefore: \f$ n > 3 \f$
//!
//!  Number of free parameters:
//!   - towers: 0
//!   - jets:   5
//!     - 0: Normalization \f$ c \f$
//!     - 1: Mean \f$ \mu_{0} \f$ of main Gaussian
//!     - 2: Width \f$ \sigma_{0} \f$ of main Gaussian
//!     - 3: Mean \f$ \mu_{1} \f$ of tail Gaussian
//!     - 4: Width \f$ \sigma_{1} \f$ of tail Gaussian
//!   - tracks: 0
//!   - global: 1
// ------------------------------------------------------------------------
class SmearTwoGauss : public Parametrization { 
 public:
  //!  \brief Constructor
  //!  \param ptBinEdges  Edges of pt bins of tail Gaussian
  //!  \param meanRespPar Parameters of mean response parametrization
  //!  \param rParScales  Response parameter scales
  // ------------------------------------------------------------------------
  SmearTwoGauss(const std::vector<double>& ptBinEdges, const std::vector<double>& meanRespPar, const std::vector<double>& rParScales)
    : Parametrization(0,3+3*(ptBinEdges.size()-1),0,1),
    meanRespPar_(meanRespPar),
    ptBinEdges_(ptBinEdges),
    nPtBins_(static_cast<int>(ptBinEdges_.size())-1),
    ptMin_(ptBinEdges_.front()),
    ptMax_(ptBinEdges_.back()),
    respParScales_(rParScales)
    {
      // Duplicate jet parameter scales for pt binned
      // function part
      int nRespParScales = static_cast<int>(respParScales_.size());
      if( nRespParScales == 6 ) {
	for(int ptBin = 1; ptBin < nPtBins_; ptBin++) {
	  for(int i = 3; i < nRespParScales; i++) {
	    respParScales_.push_back( respParScales_.at(i) );
	  }
	}
      }

      print();

      assert( meanRespPar.size() == 2 );
      for(int ptBin = 0; ptBin < nPtBins_; ptBin++) {
	assert( ptBinEdges_.at(ptBin) < ptBinEdges_.at(ptBin+1) );
      }
      assert( 0.0 <= ptMin_ && ptMin_ < ptMax_ );
      assert( respParScales_.size() == nJetPars() );
    }

  const char* name() const { return "SmearTwoGauss";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return 0.;
  }

  //!  \brief Returns probability density of response
  //!  \param x   TMeasurement::E is response for which the 
  //!             probability density is returned
  //!  \param par Pointer to parameters
  //!  \return Probability density of response
  // ------------------------------------------------------------------------
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    // Central Gaussian
    double u0 = meanRespPar_.at(0) + meanRespPar_.at(1)*x->pt;
    double a1 = respParScales_.at(0)*par[0];
    double a2 = respParScales_.at(1)*par[1];
    double a3 = respParScales_.at(2)*par[2];
    double s0 = sqrt( a1*a1/x->pt/x->pt + a2*a2/x->pt + a3*a3 );

    // Offset to ptBin
    int idx = 3*findPtBin(x->pt);

    // Normalization
    double c  = respParScales_.at(3+idx)*par[3+idx];
    if(c < 0.) c = 0.;
    if(c > 1.) c = 1.;

    // Tail Gaussian
    double u1 = respParScales_.at(4+idx)*par[4+idx];
    double s1 = respParScales_.at(5+idx)*par[5+idx];

    // Combined probability density
    double p  =       c  / sqrt(2* M_PI) / s0 * exp(-pow((x->E - u0) / s0, 2) / 2);
    p        += (1. - c) / sqrt(2* M_PI) / s1 * exp(-pow((x->E - u1) / s1, 2) / 2);

    return p;
  }



  //!  \brief Returns probability density of true pt multiplied by normalization
  //!         of dijet probability (see also \p SmearDiJet::truthPDF(t)).
  //!  \param x   \p TMeasurement::pt is true pt for which the
  //!             probability density is returned
  //!  \param par Pointer to parameters
  //!  \return Probability density of truth times normalization of dijet probability
  // ------------------------------------------------------------------------
  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    double p = 0.;
    if( ptMin_ < x->pt && x->pt < ptMax_ ) {
      double n    = par[0];        // The exponent
      double k    = n - 3.;
      double norm = k / ( pow(ptMin_,-k) - pow(ptMax_,-k) );
      
      p = norm / pow( x->pt, n );
    }

    return p;
  }


 private:
  const std::vector<double> meanRespPar_;    //!< Parameters of mean response
  const std::vector<double> ptBinEdges_;     //!< Pt bin edges
  const int    nPtBins_;                     //!< Number of pt bins
  const double ptMin_;                       //!< Minimum of non-zero range of truth pdf
  const double ptMax_;                       //!< Maximum of non-zero range of truth pdf
  std::vector<double> respParScales_;  //!< Response parameter scales



  //!  \brief Return the pt bin for a given \p pt value
  //!
  //!  There are \p nPtBins_ from 0 to <tt> nPtBins_-1 </tt>.
  //!  The underflow bin is also 0, the overflow bin also
  //!  <tt> nPtBins_-1 </tt>.
  // ------------------------------------------------------------------------
  int findPtBin(double pt) const {
    int ptBin = 0;
    for(int bin = 0; bin < nPtBins_; bin++) {
      if( pt > ptBinEdges_.at(bin) )
	ptBin = bin;
      else
	break;
    }

    return ptBin;
  }



  //!  \brief Print some initialization details
  // ------------------------------------------------------------------------
  void print() const {
    std::cout << "Parametrization class '" << name() << "'\n";
    std::cout << "  Probability density of ptTrue spectrum:\n";
    std::cout << "    Powerlaw\n";
    std::cout << "    " << ptMin_ << " < ptTrue < " << ptMax_ << " GeV\n";
    std::cout << "  Probability density of response:\n";
    std::cout << "    Mean response '" << meanRespPar_.at(0) << " + " << meanRespPar_.at(1) << " * pt'\n";
    std::cout << "    " << nPtBins_ << " ptTrue bins:\n";
    for(int i = 0; i < nPtBins_; i++) {
      std::cout << "      " << i << ": " << ptBinEdges_.at(i) << " - " << ptBinEdges_.at(1+i) << " GeV\n";
    }
    std::cout << "    Response parameter scale factors:\n";
    for(size_t i = 0; i < respParScales_.size(); i++) {
      std::cout << "      " << i << ": " << respParScales_.at(i) << "\n";
    }
  }
};



//!  \brief Parametrization used for SmearFunction estimation
//!
//!  This pdf consists of a central Gaussians around 1 and a
//!  step function (Histogram) parametrizing the tail:
//!  \f[
//!   p(r) = c \cdot G(1,\sigma) + (1-c) \cdot H(b_{i})
//!  \f]
//!
//!  Number of free parameters:
//!   - towers: 0
//!   - jets:   12
//!     - 0:      Normalization \f$ c \f$
//!     - 1:      Width \f$ \sigma \f$ of main Gaussian
//!     - 2 - 11: Bin content of histogram
//!   - tracks: 0
//!   - global: 0
// ------------------------------------------------------------------------
class SmearStepGauss : public Parametrization
{ 
 public:
  SmearStepGauss(double min, double max, int nsteps, const std::vector<double>& scale)
    : Parametrization(0,nsteps+3,0,0),
    nGausPar_(3),
    min_(min),
    max_(max),
    nStepPar_(nsteps),
    stepWidth_((max_ - min_)/nStepPar_),
    scale_(scale)
    {
      for(int i = 0; i < nStepPar_+1; i++) {
	binEdges_.push_back( min_ + i*stepWidth_ );
      }
      assert( 0.0 <= min_ && min_ < max_ );
      assert( scale_.size() >= nJetPars() );
    }
  ~SmearStepGauss() { binEdges_.clear(); }

  const char* name() const { return "SmearStepGauss";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const { return 0.; }

  //!  \brief Returns probability density of response
  //!  \param x   TMeasurement::E is response for which the 
  //!             probability density is returned
  //!  \param par Pointer to parameters
  //!  \return Probability density of response
  // ------------------------------------------------------------------------
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    // Probability density from Gaussian part
    double c = scale_.at(0)*par[0];
    double u = scale_.at(1)*par[1];
    double s = scale_.at(2)*par[2];
    
    double p = c / sqrt(2* M_PI ) / s * exp(-pow((x->E - u) / s, 2) / 2);
    
    
    // Probability density from step function part
    // Copy parameter values to temp vector for normalization
    std::vector<double> stepPar = std::vector<double>(nStepPar_,1.);
    double norm = 0.;
    for(int i = 0; i < nStepPar_; i++) {
      stepPar.at(i) = scale_.at(nGausPar_+i)*par[nGausPar_+i];
      if( stepPar.at(i) < 0. ) stepPar.at(i) = 0.;
      norm += stepPar.at(i);
    }
    norm = (1.-c)/norm/stepWidth_;


    // Find parameter idx for given response x
    // and add content to p
    if( x->E >= binEdges_.front() && x->E < binEdges_.back() ) {
      int i = 0;
      while( x->E > binEdges_.at(1+i) ) i++;
      p += norm*stepPar.at(i);
    }

    return p;
  }

 private:
  const int    nGausPar_;              //!< Number of parameters of Gaussian part of pdf
  const double min_;                   //!< Minimum of non-zero range of pdf
  const double max_;                   //!< Maximum of non-zero range of pdf
  const int    nStepPar_;              //!< Number of parameters of step part of pdf
  const double stepWidth_;             //!< Step width
  const std::vector<double> scale_;    //!< Parameters scales
  std::vector<double> binEdges_;       //!< Edges of steps
};



//!  \brief Parametrization used for SmearFunction estimation
//!
//!  The response pdf consists of a central Gaussians around 1 and an
//!  interpolated step function (Histogram) parametrizing the tail:
//!  \f[
//!   p(r) = c \cdot G(1,\sigma) + (1-c) \cdot H_{inter}(b_{i})
//!  \f]
//!
//!  The truth pdf is \f$ f(t) = N/t^{n} \f$. The normalization
//!  is such that the dijet probability
//!  \f[
//!    p = \int\;dt\;f(t)r(m_{1}/t)r(m_{2}/t)
//!  \f]
//!  is normalized to one. Therefore: \f$ n > 3 \f$
//!
//!  Number of free parameters:
//!   - towers: 0
//!   - jets:   12
//!     - 0:      Normalization \f$ c \f$
//!     - 1:      Width \f$ \sigma \f$ of main Gaussian
//!     - 2 - 11: Bin content of histogram
//!   - tracks: 0
//!   - global: 1
//!     - The exponent of the truth spectrum
// ------------------------------------------------------------------------
class SmearStepGaussInter : public Parametrization
{ 
 public:
  //!  \brief Constructor
  //!  \param tMin          Minimum of non-zero part of truth pdf
  //!  \param tMax          Maximum of non-zero part of truth pdf
  //!  \param rMin          Minimum of binned part of response pdf \f$ H_{inter} \f$
  //!  \param rMax          Maximum of binned part of response pdf \f$ H_{inter} \f$
  //!  \param ptDijetMin    Minimum of dijet pt
  //!  \param ptDijetMax    Maximum of dijet pt
  //!  \param rNBins        Number of bins of binned part of response pdf \f$ H_{inter} \f$
  //!  \param rParScales    Jet parameter scales
  //!  \param meanRespPar   Parameters for mean of Gaussian
  // ------------------------------------------------------------------------
  SmearStepGaussInter(double tMin, double tMax, double rMin, double rMax, int rNBins, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales, const std::vector<double>& gaussPar)
    : Parametrization(0,rNBins+1,0,1),
    tMin_(tMin),
    tMax_(tMax),
    rMin_(rMin),
    rMax_(rMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    nStepPar_(rNBins),
    binWidth_((rMax_ - rMin_)/nStepPar_),
    respParScales_(rParScales),
    gaussPar_(gaussPar)
    {
      for(int i = 0; i < nStepPar_+1; i++) {
	binCenters_.push_back( rMin_ + (0.5+i)*binWidth_ );
      }
      assert( 0.0 <= tMin_ && tMin_ < tMax_ );
      assert( 0.0 <= rMin_ && rMin_ < rMax_ );
      assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
      assert( respParScales_.size() >= nJetPars() );
      assert( gaussPar_.size() >= 4 );

      print();

      // Integral over dijet resolution for truth pdf
      ptDijetCutInt_ = new TH1D("norm","",400,tMin_,tMax_);
    }

  ~SmearStepGaussInter() { binCenters_.clear(); delete ptDijetCutInt_; }

  const char* name() const { return "SmearStepGaussInter";}

  virtual bool needsUpdate() const { return true; }

  virtual void update(const double * par) {
    std::cout << "'" << name() << "': Updating ptDijet cut integral... ";
    ptDijetCutInt_->Reset();

    for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
      TMeasurement x;
      x.pt = ptDijetCutInt_->GetBinCenter(bin);
      double integral = 0.;
      int nSteps = 400;
      double dPtDijet = (ptDijetMax_ - ptDijetMin_) / nSteps;
      for(int i = 0; i < nSteps; i++) {
	double ptDijet = ptDijetMin_ + i*dPtDijet;
	x.E = ptDijet / x.pt;
	double prob = correctedJetEt(&x,par);
	prob *= x.pt;
	integral += prob;
      }
      integral /= sqrt(2.);
      integral *= dPtDijet;
      ptDijetCutInt_->SetBinContent(bin,integral);
    }
    std::cout << "ok\n";
  }


  double correctedTowerEt(const TMeasurement *x,const double *par) const { return 0.; }


  //!  \brief Returns probability density of response
  //!  \param x   TMeasurement::E is the response,
  //!  \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
  //!  \return Probability density of response
  // ------------------------------------------------------------------------
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    // Probability density from Gaussian part
    double c = respParScales_.at(0)*par[0];
    double u = gaussPar_.at(0);
    double a1 = gaussPar_.at(1);
    double a2 = gaussPar_.at(2);
    double a3 = gaussPar_.at(3);
    double s = sqrt( a1*a1/x->pt/x->pt + a2*a2/x->pt + a3*a3 );

    double p = c / sqrt(2* M_PI ) / s * exp(-pow((x->E - u) / s, 2) / 2);
        
    // Probability density from step function part
    // Copy parameter values to temp vector for normalization
    std::vector<double> stepPar(nStepPar_,1.);
    double norm = 0.;
    for(int i = 0; i < nStepPar_; i++) {
      stepPar.at(i) = respParScales_.at(1+i)*par[1+i];
      if( stepPar.at(i) < 0. ) stepPar.at(i) = 0.;
      norm += stepPar.at(i);
    }
    norm = (1.-c)/norm/binWidth_;
    
    p += norm*interpolate(x->E,stepPar);
    
    return p;
  }



  //!  \brief Returns probability density of true pt multiplied by normalization
  //!         of dijet probability (see also \p SmearDiJet::truthPDF(t)).
  //!  \param x   \p TMeasurement::pt is truth for which the 
  //!             probability density is returned
  //!  \param par Pointer to parameters (parameters are multiplied by the corresponding scale)
  //!  \return Probability density of truth times normalization of dijet probability
  // ------------------------------------------------------------------------
  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    // Norm of probability of dijet event configuration
    double norm = 0.;
    for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
      double pt = ptDijetCutInt_->GetBinCenter(bin);
      norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetBinContent(bin) * ptDijetCutInt_->GetXaxis()->GetBinWidth(1);
    }

    double p = 0.;
    if( norm ) {
      p = truthPDF(x->pt,par[0]);
      p /= norm;
    }
    
    return p;
  }


 private:
  const double tMin_;                   //!< Minimum of non-zero range of truth pdf
  const double tMax_;                   //!< Maximum of non-zero range of truth pdf
  const double rMin_;                   //!< Minimum of non-zero range of response pdf
  const double rMax_;                   //!< Maximum of non-zero range of response pdf
  const double ptDijetMin_;             //!< Minimum of pt dijet
  const double ptDijetMax_;             //!< Maximum of pt dijet
  const int    nStepPar_;               //!< Number of parameters of step part of pdf
  const double binWidth_;               //!< Bin width
  const std::vector<double> respParScales_;     //!< Parameter scales
  const std::vector<double> gaussPar_;  //!< Parameters of Gaussian part of response pdf
  std::vector<double> binCenters_;      //!< Centers of response bins
  TH1D * ptDijetCutInt_;	        //!< Integral over dijet resolution for truth pdf for different t



  //!  \brief Linear interpolation between two bins
  //!
  //!  Returns the linearly weighted mean value of two adjacent bin
  //!  contents, where the two bins are the ones whose bin centers
  //!  are closest to \p r. The contents are weighted with the distance
  //!  of r from the corresponding bin center. The interpolation
  //!  assumes \p nStepPar_ bins. At the left and right edges of the
  //!  binned range, the interpolation is done assuming 0 bin content
  //!  outside the binned range.
  //!
  //!  \param r   Response
  //!  \param par Parameters / bin content to be interpolated
  //!  \return    Interpolated bin content
  // ------------------------------------------------------------------------
  double interpolate(double r, const std::vector<double>& par) const {
    double p = 0.;     // The interpolated parameter

    // Check that r is in range of binning +/- 0.5*mBinWidth
    if( r > rMin_ - 0.5*binWidth_  &&  r < rMax_ + 0.5*binWidth_ ) {
      // Find index i of the next larger bin center
      unsigned int i = 0;                              
      while( r > binCenters_.at(i) ) i++;
      assert( i < binCenters_.size() );

      double dx  = binCenters_.at(i) - r;  // Distance to next larger bin center
      dx        /= binWidth_;              // dx relative to bin width
      double a   = 0.;
      double b   = 0.;

      // Find bin contents to be interpolated
      if( i == 0 ) { // Left edge
	a = 0.;
	b = par.front();
      }
      else if( i == binCenters_.size()-1 ) { // Right edge
	a = par.back();
	b = 0.;
      }
      else { // Central part
	a = par.at(i-1);
	b = par.at(i);
      }

      // Interpolate contents
      p = dx * a + (1-dx) * b;
    }

    return p;
  }



  //!  \brief Returns probability density (not normalised!) of truth
  //!         considering cuts on dijet pt
  //!
  //!  The probability density of \p pt if cuts
  //!  \f$ \texttt{ptDijetMin\_} < p^{\textrm{dijet}}_{T} < \texttt{ptDijetMax\_} \f$
  //!  are applied:
  //!  \f[ t^{-n} \int^{\texttt{ptDijetMax\_}}_{\texttt{ptDijetMin\_}}
  //!      dx\,\frac{r(x/t) \cdot t}{\sqrt{2}} \f]
  // ------------------------------------------------------------------------
  double truthPDF(double pt, double n) const {
    double prob = 0.;
    if( tMin_ < pt && pt < tMax_ ) {
      prob = 1. / pow( pt, n );
      int bin = ptDijetCutInt_->FindBin(pt);
      prob *= ptDijetCutInt_->GetBinContent(bin);
    }

    return prob;
  }



  //!  \brief Print some initialization details
  // ------------------------------------------------------------------------
  void print() const {
    std::cout << "Parametrization class '" << name() << "'\n";
    std::cout << "  Probability density of ptTrue spectrum:\n";
    std::cout << "    Powerlaw\n";
    std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
    std::cout << "    " << ptDijetMin_ << " < ptDijet < " << ptDijetMax_ << " GeV\n";
    std::cout << "  Probability density of response:\n";
    std::cout << "    Non-zero for " << rMin_ << " < R < " << rMax_ << "\n";
    std::cout << "    " << nStepPar_ << " step parameters\n";
    std::cout << "    " << gaussPar_.size() << " Gauss parameters: ";
    for(size_t i = 0; i < gaussPar_.size(); i++) {
      std::cout << gaussPar_.at(i);
      if( i != gaussPar_.size() - 1 ) std::cout << ", ";
    }
    std::cout << std::endl;
  }
};



//!  \brief Parametrization used for SmearFunction estimation
//!
//!  The response pdf consists of a central Gaussians and an
//!  interpolated step function (with \f$ N_{j} \f$ response bins)
//!  parametrizing the tail:
//!  \f[
//!   p(r) = c_{i} \cdot G(\mu(t),\sigma(t)) + (1-c_{i}) \cdot H_{inter}(b_{i,j})
//!  \f]
//!  The mean of the Gaussian is parameterized by
//!  \f[
//!   \mu(t) = 1.008 - 8E-06\cdot t
//!  \f]
//!  from a fit on the mean of L2L3 corrected jets.
//!  The width of the Gaussian is parameterized by
//!  \f[
//!   \frac{\sigma(t)}{p_{T}} = \sqrt{\frac{a^{2}_{1}}{t} + \frac{a^{2}_{2}}{t} + a^{2}_{3}},
//!  \f]
//!  where the \f$ a_{i} \f$ are parameters of the fit. (So far,
//!  \f$ a_{1} = 0 \f$.)
//!  The step function part is additionally binned in \f$ N_{i} \f$ \p t bins.
//!
//!  The truth pdf is \f$ f(t) = N/t^{n} \f$. The normalization
//!  is such that the dijet probability
//!  \f[
//!    p = \int\;dt\;f(t)r(m_{1}/t)r(m_{2}/t)
//!  \f]
//!  is normalized to one. Therefore: \f$ n > 3 \f$
//!
//!  Number of free parameters:
//!   - towers: 0
//!   - jets:   \f$ 2 + N_{i}\cdot(1+n_{j}) \f$
//!     - \f$ 0 - 1 \f$:  \f$ a_{1} \f$ and \f$ a_{2} \f$    Normalization \f$ c \f$
//!     - \f$ 2 + N_{i}\cdot(N_{j}+1) \f$: Normalization \f$ c_{i} \f$
//!     - \f$ 3\ldots + N_{i} \cdot (N_{j}+1) \f$: Parameters \f$ b_{i,j} \f$ of the step function
//!   - tracks: 0
//!   - global: 1
//!     - The exponent of the truth spectrum
// ------------------------------------------------------------------------
class SmearStepGaussInterPtBinned : public Parametrization
{
 public:
  //!  \brief Constructor
  //!  \param ptBinEdges  Edges of pt bins of step function
  //!  \param rMin        Minimum response of step function \f$ H_{inter} \f$
  //!  \param rMax        Maximum response of step function \f$ H_{inter} \f$
  //!  \param rNBins      Number of bins per pt bin of step function \f$ H_{inter} \f$
  //!  \param rParScales  Response parameter scales
  // ------------------------------------------------------------------------
  SmearStepGaussInterPtBinned(const std::vector<double>& ptBinEdges, double rMin, double rMax, int rNBins, const std::vector<double>& rParScales)
    : Parametrization( 0, 2 + (ptBinEdges.size()-1) * (1+rNBins), 0, 1 ),
    ptBinEdges_(ptBinEdges),
    ptMin_(ptBinEdges_.front()),
    ptMax_(ptBinEdges_.back()),
    nPtBins_(static_cast<int>(ptBinEdges_.size())-1),
    rMin_(rMin),
    rMax_(rMax),
    nGaussPar_(2),
    nStepParPerPtBin_(rNBins),
    binWidth_((rMax_ - rMin_)/nStepParPerPtBin_),
    respParScales_(rParScales)
    {
      // Find step function bin centers
      for(int i = 0; i < nStepParPerPtBin_+1; i++) {
	binCenters_.push_back( rMin_ + (0.5+i)*binWidth_ );
      }
      // Duplicate jet parameter scales for pt binned step
      // function part
      int nRespParScales = static_cast<int>(respParScales_.size());
      if( nRespParScales == nGaussPar_ + 1 + nStepParPerPtBin_ ) {
	for(int ptBin = 1; ptBin < nPtBins_; ptBin++) {
	  for(int i = nGaussPar_; i < nRespParScales; i++) {
	    respParScales_.push_back( respParScales_.at(i) );
	  }
	}
      }

      print();

      for(int ptBin = 0; ptBin < nPtBins_; ptBin++) {
	assert( ptBinEdges_.at(ptBin) < ptBinEdges_.at(ptBin+1) );
      }
      assert( 0.0 <= ptMin_ && ptMin_ < ptMax_ );
      assert( 0.0 <= rMin_ && rMin_ < rMax_ );
      assert( respParScales_.size() == nJetPars() );
    }

  ~SmearStepGaussInterPtBinned() {};

  const char* name() const { return "SmearStepGaussInterPtBinned";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const { return 0.; }



  //!  \brief Returns probability density of response
  //!  \param x \p TMeasurement::E is the response,
  //!           \p TMeasurement::pt is the true pt
  //!  \param par Pointer to parameters
  //!  \return Probability density of response
  // ------------------------------------------------------------------------
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    // Find index to start of step parameters of this pt bin
    // Per pt bin there are
    //  - 1 normalization constant,
    //  - nStepParPerPtBin_ step parameters
    int idx = nGaussPar_ + findPtBin(x->pt)*(1+nStepParPerPtBin_);

    // Probability density from Gaussian part
    double c = respParScales_.at(idx)*par[idx];
    double u = meanResponse(x->pt);

    double a1 = 0.;
    double a2 = respParScales_.at(0)*par[0];
    double a3 = respParScales_.at(1)*par[1];
    double s  = sqrt( a1*a1/x->pt/x->pt + a2*a2/x->pt + a3*a3 );

    double p = c / sqrt(2* M_PI ) / s * exp(-pow((x->E - u) / s, 2) / 2);
    
    
    // Probability density from step function part
    // Copy parameter values to temp vector for normalization
    std::vector<double> stepPar(nStepParPerPtBin_,1.);
    double norm = 0.;
    idx++;
    for(int i = 0; i < nStepParPerPtBin_; i++) {
      stepPar.at(i) = respParScales_.at(idx+i)*par[idx+i];
      if( stepPar.at(i) < 0. ) stepPar.at(i) = 0.;
      norm += stepPar.at(i);
    }
    norm = (1.-c)/norm/binWidth_;
    
    p += norm*interpolate(x->E,stepPar);
    
    return p;
  }



  //!  \brief Returns probability density of true pt multiplied by normalization
  //!         of dijet probability (see also \p SmearDiJet::truthPDF(t)).
  //!  \param x   \p TMeasurement::pt is true pt for which the
  //!             probability density is returned
  //!  \param par Pointer to parameters
  //!  \return Probability density of truth times normalization of dijet probability
  // ------------------------------------------------------------------------
  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    double p = 0.;
    if( ptMin_ < x->pt && x->pt < ptMax_ ) {
      double n    = par[0];        // The exponent
      double k    = n - 3.;
      double norm = k / ( pow(ptMin_,-k) - pow(ptMax_,-k) );
      
      p = norm / pow( x->pt, n );

/*       double norm = ( pow(ptMax_,3) + pow(ptMin_,3) ) / 3.; */
/*       p = 1./norm; */
    }

    return p;
  }


 private:
  const std::vector<double> ptBinEdges_;  //!< Pt bin edges
  const double ptMin_;                  //!< Minimum of non-zero range of truth pdf
  const double ptMax_;                  //!< Maximum of non-zero range of truth pdf
  const int    nPtBins_;                //!< Number of pt bins
  const double rMin_;                   //!< Minimum of non-zero range of response pdf
  const double rMax_;                   //!< Maximum of non-zero range of response pdf
  const int    nGaussPar_;              //!< Number of parameters of gauss part of pdf
  const int    nStepParPerPtBin_;       //!< Number of parameters of step part of pdf per pt bin
  const double binWidth_;               //!< Response bin width
  std::vector<double> respParScales_;   //!< Jet parameter scales
  std::vector<double> binCenters_;      //!< Centers of response bins



  //!  \brief Return the pt bin for a given \p pt value
  //!
  //!  There are \p nPtBins_ from 0 to <tt> nPtBins_-1 </tt>.
  //!  The underflow bin is also 0, the overflow bin also
  //!  <tt> nPtBins_-1 </tt>.
  // ------------------------------------------------------------------------
  int findPtBin(double pt) const {
    int ptBin = 0;
    for(int bin = 0; bin < nPtBins_; bin++) {
      if( pt > ptBinEdges_.at(bin) )
	ptBin = bin;
      else
	break;
    }

    return ptBin;
  }



  //!  \brief Linear interpolation between two bins
  //!
  //!  Returns the linearly weighted mean value of two adjacent bin
  //!  contents, where the two bins are the ones whose bin centers
  //!  are closest to \p r. The contents are weighted with the distance
  //!  of r from the corresponding bin center. The interpolation
  //!  assumes \p nStepParPerPtBin_ bins. At the left and right edges of the
  //!  binned range, the interpolation is done assuming 0 bin content
  //!  outside the binned range.
  //!
  //!  \param r   Response
  //!  \param par Parameters / bin content to be interpolated
  //!  \return    Interpolated bin content
  // ------------------------------------------------------------------------
  double interpolate(double r, const std::vector<double>& par) const {
    double p = 0.;     // The interpolated parameter

    // Check that r is in range of binning +/- 0.5*mBinWidth
    if( r > rMin_ - 0.5*binWidth_  &&  r < rMax_ + 0.5*binWidth_ ) {
      // Find index i of the next larger bin center
      unsigned int i = 0;
      while( r > binCenters_.at(i) ) i++;
      assert( i < binCenters_.size() );

      double dx  = binCenters_.at(i) - r;  // Distance to next larger bin center
      dx        /= binWidth_;              // dx relative to bin width
      double a   = 0.;
      double b   = 0.;

      // Find bin contents to be interpolated
      if( i == 0 ) { // Left edge
	a = 0.;
	b = par.front();
      }
      else if( i == binCenters_.size()-1 ) { // Right edge
	a = par.back();
	b = 0.;
      }
      else { // Central part
	a = par.at(i-1);
	b = par.at(i);
      }

      // Interpolate contents
      p = dx * a + (1-dx) * b;
    }

    return p;
  }


  //!  \brief Return the mean response for a given \p pt
  //!
  //!  The mean response \p R is linearly parametrized as
  //!  \f[ R = 1.008 - 8\cdot10^{-6}\cdot p_{T} \f],
  //!  where the parameter values come from a fit of 
  //!  \f$ <p^{jet}_{T} / p^{gen}_{T}> \f$ of Summer08
  //!  L2L3 corrected jets.
  // ------------------------------------------------------------------------
  double meanResponse(double pt) const {
    double a = 1.008;
    double b = -8E-06;

    return a + b*pt;
  }



  //!  \brief Print some initialization details
  // ------------------------------------------------------------------------
  void print() const {
    std::cout << "Parametrization class '" << name() << "'\n";
    std::cout << "  Probability density of ptTrue spectrum:\n";
    std::cout << "    Powerlaw\n";
    std::cout << "    " << ptMin_ << " < ptTrue < " << ptMax_ << " GeV\n";
    std::cout << "  Probability density of response:\n";
    std::cout << "    Non-zero for " << rMin_ << " < R < " << rMax_ << "\n";
    std::cout << "    " << nGaussPar_ << " gauss parameters\n";
    std::cout << "    " << nStepParPerPtBin_ << " step parameters per pt bin\n";
    std::cout << "    " << nPtBins_ << " ptTrue bins:\n";
    for(int i = 0; i < nPtBins_; i++) {
      std::cout << "      " << i << ": " << ptBinEdges_.at(i) << " - " << ptBinEdges_.at(1+i) << " GeV\n";
    }
    std::cout << "    Response parameter scale factors:\n";
    for(int i = 0; i < 1+nGaussPar_+nStepParPerPtBin_; i++) {
      std::cout << "      " << i << ": " << respParScales_.at(i) << "\n";
      if( i == (nGaussPar_-1) || i == nGaussPar_ ) 
	std::cout << "     ----------\n";
    }
  }
};




//!  \brief Jet parametrization using power law from Groom
//!
//!  This parametrization has 0 tower parameters,3
//!  jet parameters and 0 global jet parameters
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class GroomParametrization: public Parametrization {
public:
  GroomParametrization() : Parametrization(0,3,0,0) {
    gsl_set_error_handler_off();
  }
  const char* name() const { return "GroomParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //return 1/par[0] * x->pt; int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    f_par p(x->pt,par);
    
    F.function = &f;
    F.params = &p;

    double x_lo = 0.1 * x->pt;
    double x_hi = 5 * x->pt;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    if(gsl_root_fsolver_set(s,&F,x_lo,x_hi)) {
      //std::cout << "Warning: root not bracketed\n";
      return 0;
    }
    int status, iter = 0;
    double r;
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo,x_hi,0, eps_);
    }
    while (status == GSL_CONTINUE && iter < max_iter_);
    assert(status == GSL_SUCCESS);
    gsl_root_fsolver_free(s);
    return r;
  }

  double inverseJetCorrection(const TMeasurement *x,const double *par) const { 
    return (par[0] - par[1] * pow(0.01 * x->pt, -par[2])) * x->pt;
  }
 
  bool hasInvertedCorrection() const { return true;}

 private:
  static const double eps_ =  1.0e-8;
  static const int max_iter_ = 20;
  struct f_par {
    double y;
    const double *par;
    f_par(double y, const double *par) : y(y), par(par) {}
  };

  static double f(double x, void* params) {
    f_par* p = (f_par*)params;
    //return y - (par[0] - par[1] * pow(0.01 * x,-par[2])) * x;
    return p->y - (p->par[0] - p->par[1] * pow(0.01 * x, -p->par[2])) * x;
  };
};
#endif
