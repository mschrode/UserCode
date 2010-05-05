//
//  $Id: Parametrization.h,v 1.60 2010/04/27 16:03:41 stadie Exp $
//
#ifndef CALIBCORE_PARAMETRIZATION_H
#define CALIBCORE_PARAMETRIZATION_H

#include <cmath>
#include <vector>

#include "CalibData.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
     
class TH1;


//!  \brief Abstract base class for parametrizations of
//!         correction functions
//!
//!  Interface to different parametrizations for the 
//!  underlying hypothesis of the GlobalFit. Allows
//!  to correct a tower or jet measurement.
//!  \author Hartmut Stadie
//!  \date Thu Apr 03 17:09:50 CEST 2008
//!  $Id: Parametrization.h,v 1.60 2010/04/27 16:03:41 stadie Exp $
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
  //!           members (see also Measurement):
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
  virtual double correctedTowerEt(const Measurement *x,const double *par) const = 0;


  //!  \brief Corrects the measured jet Et
  //!
  //!  The parameters of the jet correction function
  //!  may depend on the jet Et, eta, and phi.
  //!
  //!  \param x Jet Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also Measurement):
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
  virtual double correctedJetEt(const Measurement *x,const double *par) const = 0;

 
  //!  \brief Returns the expected signal of a track in the Calorimeter
  //!
  //!  \param x Track Et measurement from which the expected response
  //!           of the calorimeter is calculated (see also TTrack).
  //!  \param par Parameters of the response function of this track
  //!  \return The expected calorimeter response
  // -----------------------------------------------------------------
  virtual double GetExpectedResponse(const Measurement *x,const double *par) const { return x->pt;}


  //!  \brief Corrects the measured jet Et with global
  //!         correction function
  //!
  //!  The parameters of the global jet correction function
  //!  are independent of the jet Et, eta, and phi.
  //!
  //!  \param x Jet Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also Measurement):
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
  virtual double correctedGlobalJetEt(const Measurement *x,const double *par) const { return x->pt;}


  //! Returns probability density of measured jet pt given a true jet pt
  virtual double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const { return 0.; }
  virtual double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const { return 0.; }
  //! Returns probability density of true jet pt
  virtual double pdfPtTrue(double ptTrue, const double *par) const { return 0.; }
  virtual double pdfPtTrueError(double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const { return 0.; }
  //! Returns probability density of response given a true jet pt
  virtual double pdfResponse(double r, double ptTrue, const double *par) const { return 0.; }
  virtual double pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const { return 0.; }
  virtual double pdfDijetAsym(double a, double ptTrue, const double *par) const { return 0.; }


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
  //!  \param x Measurement describing the true jet properties
  //!  \param par Parameters of the jet response function
  //!  \return the expected calorimeter response
  // -----------------------------------------------------------------
  virtual double inverseJetCorrection(const Measurement *x,const double *par) const { return -1;}

  
  //!  \brief Get the number of parameters of the track p
  //!  \return Number of parameters of the track parametrization
  // -----------------------------------------------------------------
  virtual bool hasInvertedCorrection() const { return false;}

  virtual bool needsUpdate() const { return false; }
  virtual void update(const double * par) {;}

 protected:
  //interpolotion code from JetMETObjects/Utilities
  static double quadraticInterpolation(double fZ, const double fX[3], const double fY[3])
  {
    // Quadratic interpolation through the points (x[i],y[i]). First find the parabola that
    // is defined by the points and then calculate the y(z).
    float D[4],a[3];
    D[0] = fX[0]*fX[1]*(fX[0]-fX[1])+fX[1]*fX[2]*(fX[1]-fX[2])+fX[2]*fX[0]*(fX[2]-fX[0]);
    D[3] = fY[0]*(fX[1]-fX[2])+fY[1]*(fX[2]-fX[0])+fY[2]*(fX[0]-fX[1]);
    D[2] = fY[0]*(pow(fX[2],2)-pow(fX[1],2))+fY[1]*(pow(fX[0],2)-pow(fX[2],2))+fY[2]*(pow(fX[1],2)-pow(fX[0],2));
    D[1] = fY[0]*fX[1]*fX[2]*(fX[1]-fX[2])+fY[1]*fX[0]*fX[2]*(fX[2]-fX[0])+fY[2]*fX[0]*fX[1]*(fX[0]-fX[1]);
    if (D[0] != 0)
      {
        a[0] = D[1]/D[0];
        a[1] = D[2]/D[0];
        a[2] = D[3]/D[0];
      }
    else
      {
        a[0] = 0.0;
        a[1] = 0.0;
        a[2] = 0.0;
      }
    return a[0]+fZ*(a[1]+fZ*a[2]);
  }
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
    
    double correctedTowerEt(const Measurement *x,const double *par) const {
      return x->pt;
    }
    
    double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
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
  
  double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
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
    
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
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
  
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->EMF + x->HadF + x->OutF;
    //return x->EMF + par[0]*x->HadF + par[1]*log(x->pt) + par[2];
  }
  
  double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return par[1] * x->HadF + par[2] * x->EMF + x->OutF + par[0];
  }
  
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedGlobalJetEt(const Measurement *x, const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
    return par[1] * x->EMF + par[2] * x->HadF + x->OutF + par[0];
  }

  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
    return par[0] * x->HadF + x->EMF + x->OutF;
  }

  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
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
    
  double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const 
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
  
  double correctedJetEt(const Measurement *x,const double *par) const 
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
   
  double GetExpectedResponse(const Measurement *x,const double *par) const   
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
  L2L3JetParametrization() : Parametrization(0,5,0,4) {} 
  const char* name() const { return "L2L3JetParametrization";}
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }

  //!  \brief Code from L2RelativeCorrector
  //!  \code
  //!  double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
  //!  double logpt = log10(pt);
  //!  double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
  //!  \endcode   
  double correctedJetEt(const Measurement *x,const double *par) const {
    double  pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt;
    double logpt = log10(pt);
    double c = par[0]+logpt*(0.1 * par[1]+logpt *(0.01* par[2]+logpt*(0.01*par[3]+logpt*(0.01*par[4]))));
    
    if(c < 0.2) c = 0.2;
    if(c > 5.0) c = 5.0;
    //assert(c > 0);
    return c * x->pt;
  }

  //!  \brief Code from SimpleL3AbsoluteCorrector
  //!  \code
  //!  double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
  //!  double log10pt = log10(pt);
  //!  double result = p[2]+p[3]/(pow(log10pt,p[4])+p[5]);
  //!  \endcode
  double correctedGlobalJetEt(const Measurement *x,const double *par) const {
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt;
    double logpt = log10(pt);
    double c = par[0] + par[1]/(pow(logpt,par[2]) + par[3]);
    if(c < 0.3) c = 0.3;
    if(c > 3.0) c = 3.0;
    assert(c > 0);
    return  c * x->pt;
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedGlobalJetEt(const Measurement *x,const double *par) const {
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
  
  double GetExpectedResponse(const Measurement *x,const double *par) const {
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const Measurement *x,const double *par) const {
    return par[0]*x->pt;
  }

  double correctedGlobalJetEt(const Measurement *x,const double *par) const {
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

  double correctedTowerEt(const Measurement *x,const double *par) const {
    return 0.;
  }

  //!  \brief Returns probability density of response
  //!  \param x   Measurement::E is response for which the 
  //!             probability density is returned
  //!  \param par Pointer to parameters
  //!  \return Probability density of response
  // ------------------------------------------------------------------------
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double correctedGlobalJetEt(const Measurement *x,const double *par) const {
    return 0.;
  }
};



//! \brief Parametrization used for resolution function estimation
//!        with a Gaussian function
// ------------------------------------------------------------------------
class SmearGauss : public Parametrization
{ 
 public:
  //! Constructor
  SmearGauss(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales);

  ~SmearGauss();

  const char* name() const { return "SmearGauss";}

  virtual bool needsUpdate() const { return false; }

  //! Returns 0
  double correctedTowerEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedJetEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedGlobalJetEt(const Measurement *x,const double *par) const { return 0.; }

  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtTrue(double ptTrue, const double *par) const;
  double pdfResponse(double r, double ptTrue, const double *par) const;
  double pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const;
  double pdfDijetAsym(double a, double ptTrue, const double *par) const;


 private:
  const double tMin_;                   //!< Minimum of non-zero range of truth pdf
  const double tMax_;                   //!< Maximum of non-zero range of truth pdf
  const double xMin_;             //!< Minimum of pt dijet
  const double xMax_;             //!< Maximum of pt dijet
  const std::vector<double> scale_;     //!< Parameter scales
  TH1 *hPdfPtTrue_;

  double sigma(double ptTrue, const double *par) const {
    double a[3] = { 0., 0., 0. };
    for(int i = 1; i < 2; i++) {
      if( par[i] > 1E-3 ) {
	a[i] = scale_[i]*par[i];
      } else {
	a[i] = scale_[i]*1E-3;
      }
    }

    return sqrt( a[0]*a[0] + a[1]*a[1]*ptTrue + a[2]*a[2]*ptTrue*ptTrue );
  }
  double exponent(const double *par) const {
    return par[3] > 1E-3 ? scale_[3]*par[3] : scale_[3]*1E-3;
  }
  double specPar(const double *par, int i) const {
    return scale_[3+i]*par[3+i];
  }

  double corr3rdJet(double ptTrue) const;
  double pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const;
  double pdfPtTrueNotNorm(double ptTrue, const double *par) const;

  //! Print some initialization details
  void print() const;
};



//! \brief Parametrization used for resolution function estimation
//!        with a Gaussian function with one sigma (not pt-dependent)
// ------------------------------------------------------------------------
class SmearGaussPtBin : public Parametrization
{ 
 public:
  //! Constructor
  SmearGaussPtBin(double tMin, double tMax, double xMin, double xMax, const std::vector<double> &parScales, const std::vector<double> &startPar);

  ~SmearGaussPtBin();

  const char* name() const { return "SmearGaussPtBin";}

  virtual bool needsUpdate() const { return false; }

  //! Returns 0
  double correctedTowerEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedJetEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedGlobalJetEt(const Measurement *x,const double *par) const { return 0.; }

  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtTrue(double ptTrue, const double *par) const;
  double pdfResponse(double r, double ptTrue, const double *par) const;
  double pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const;
  double pdfDijetAsym(double a, double ptTrue, const double *par) const;


 private:
  const double tMin_;                   //!< Minimum of non-zero range of truth pdf
  const double tMax_;                   //!< Maximum of non-zero range of truth pdf
  const double xMin_;                   //!< Minimum of pt dijet
  const double xMax_;                   //!< Maximum of pt dijet
  const std::vector<double> scale_;     //!< Parameter scales
  TH1 *hashTablePdfPtTrue_;
  TH1 *hPdfPtTrue_;

  double sigma(const double *par) const {
    double p = par[0] > 0 ? par[0] : 1E-3;
    return scale_[0]*p;
  }
  double specSigmaPar(const double *par, int i) const {
    return scale_[1+i]*par[1+i];
  }
  double specSlopePar(const double *par, int i) const {
    return scale_[4+i]*par[4+i];
  }

  void hashPdfPtTrue(const double *par) const;
  double pdfPtTrueNotNorm(double ptTrue, const double *par) const;
  double underlyingPdfPtTrue(double ptTrue, const double *par) const;

  //! Print some initialization details
  void print() const;
};


//! \brief Parametrization used for resolution function estimation
//!        with a Gaussian function
// ------------------------------------------------------------------------
class SmearGaussExtrapolation : public Parametrization
{ 
 public:
  //! Constructor
  SmearGaussExtrapolation(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales);
  ~SmearGaussExtrapolation();

  const char* name() const { return "SmearGaussExtrapolation";}

  //! Returns 0
  double correctedTowerEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedJetEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedGlobalJetEt(const Measurement *x,const double *par) const { return 0.; }

  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtTrue(double ptTrue, const double *par) const;
  double pdfResponse(double r, double ptTrue, const double *par) const;
  double pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const;
  double pdfDijetAsym(double a, double ptTrue, const double *par) const;


 private:
  const double tMin_;                   //!< Minimum of non-zero range of truth pdf
  const double tMax_;                   //!< Maximum of non-zero range of truth pdf
  const double xMin_;             //!< Minimum of pt dijet
  const double xMax_;             //!< Maximum of pt dijet
  const std::vector<double> scale_;     //!< Parameter scales
  TH1 *hPdfPtTrue_;

  double sigma(double ptTrue, const double *par) const {
    double a[3] = { 0., 0., 0. };
    for(int i = 1; i < 2; i++) {
      if( par[i] > 1E-3 ) {
	a[i] = scale_[i]*par[i];
      } else {
	a[i] = scale_[i]*1E-3;
      }
    }

    return sqrt( a[0]*a[0] + a[1]*a[1]*ptTrue + a[2]*a[2]*ptTrue*ptTrue );
  }
  double sigmaOffsetSlope(double ptTrue, const double *par) const {
    return par[3] > 1E-3 ? scale_[3]*par[3] : scale_[3]*1E-3;
  }
  double sigmaWithOffset(double ptTrue, double pt3Rel, const double *par) const {
    return sigma(ptTrue,par) + sigmaOffsetSlope(ptTrue,par)*pt3Rel;
  }
  double exponent(const double *par) const {
    return par[5] > 1E-3 ? scale_[5]*par[5] : scale_[5]*1E-3;
  }
  double specPar(const double *par, int i) const {
    return scale_[5+i]*par[5+i];
  }

  double pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const;

  //! Print some initialization details
  void print() const;
};


//! \brief Parametrization used for resolution function estimation
//!        with a Gaussian function
// ------------------------------------------------------------------------
class SmearGaussImbalance : public Parametrization
{ 
 public:
  //! Constructor
  SmearGaussImbalance(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales);
  ~SmearGaussImbalance();

  const char* name() const { return "SmearGaussImbalance";}

  //! Returns 0
  double correctedTowerEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedJetEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedGlobalJetEt(const Measurement *x,const double *par) const { return 0.; }

  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtTrue(double ptTrue, const double *par) const;
  double pdfResponse(double r, double ptTrue, const double *par) const;
  double pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const;
  double pdfDijetAsym(double a, double ptTrue, const double *par) const;


 private:
  const double tMin_;                   //!< Minimum of non-zero range of truth pdf
  const double tMax_;                   //!< Maximum of non-zero range of truth pdf
  const double xMin_;             //!< Minimum of pt dijet
  const double xMax_;             //!< Maximum of pt dijet
  const std::vector<double> scale_;     //!< Parameter scales
  const double dPtTrue_;
  std::vector<double> parMin_;
  TH1 *hPdfPtTrue_;

  double ptTrue2Min(double ptTrue) const {
    return ptTrue*(1.-dPtTrue_);
  }
  double ptTrue2Max(double ptTrue) const {
    return ptTrue*(1.+dPtTrue_);
  }

  double param(int i, const double *par) const {
    return par[i] > parMin_[i] ? scale_[i]*par[i] : scale_[i]*parMin_[i];
  }
  double sigma(double ptTrue, const double *par) const {
    double a[3] = { parMin_[0], parMin_[1], parMin_[2] };
    for(int i = 1; i < 3; i++) {
      if( par[i] > parMin_[i] ) {
	a[i] = scale_[i]*par[i];
      } else {
	a[i] = scale_[i]*parMin_[i];
      }
    }

    return sqrt( a[0]*a[0] + a[1]*a[1]*ptTrue + a[2]*a[2]*ptTrue*ptTrue );
  }
  double exponent(const double *par) const {
    return par[3] > 1E-3 ? scale_[3]*par[3] : scale_[3]*1E-3;
  }
  double specPar(const double *par, int i) const {
    return scale_[3+i]*par[3+i];
  }


  double pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const;

  //! Print some initialization details
  void print() const;
};


//! \brief Parametrization used for resolution function estimation
//!        with a Crystal Ball function with one sigma (not pt-dependent)
// ------------------------------------------------------------------------
class SmearCrystalBallPtBin : public Parametrization
{ 
 public:
  //! Constructor
  SmearCrystalBallPtBin(double tMin, double tMax, double xMin, double xMax, const std::vector<double> &parScales, const std::vector<double> &startPar);

  ~SmearCrystalBallPtBin();

  const char* name() const { return "SmearCrystalBallPtBin";}

  virtual bool needsUpdate() const { return false; }

  //! Returns 0
  double correctedTowerEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedJetEt(const Measurement *x,const double *par) const { return 0.; }
  double correctedGlobalJetEt(const Measurement *x,const double *par) const { return 0.; }

  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const;
  double pdfPtTrue(double ptTrue, const double *par) const;
  double pdfResponse(double r, double ptTrue, const double *par) const;
  double pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const { return 0.; }
  double pdfDijetAsym(double a, double ptTrue, const double *par) const;


 private:
  const double tMin_;                   //!< Minimum of non-zero range of truth pdf
  const double tMax_;                   //!< Maximum of non-zero range of truth pdf
  const double xMin_;                   //!< Minimum of pt dijet
  const double xMax_;                   //!< Maximum of pt dijet
  const std::vector<double> scale_;     //!< Parameter scales
  TH1 *hashTablePdfPtTrue_;

  double sigma(const double *par) const {
    double p = par[0] > 0 ? par[0] : 1E-3;
    return scale_[0]*p;
  }
  double alpha(const double *par) const {
    double p = par[1] > 0 ? par[1] : 1E-3;
    return scale_[1]*p;
  }
  double n(const double *par) const {
    double p = par[2] > 0 ? par[2] : 1E-3;
    return scale_[2]*p;
  }
  double specSigmaPar(const double *par, int i) const {
    return scale_[3+i]*par[3+i];
  }
  double specSlopePar(const double *par, int i) const {
    return scale_[6+i]*par[6+i];
  }

  double crystalBallFunc(double x, double mean, double sigma, double alpha, double n) const;
  double crystalBallNorm(double mean, double sigma, double alpha, double n) const;
  double crystalBallInt(double mean, double sigma, double alpha, double n, double min, double max) const;

  void hashPdfPtTrue(const double *par) const;
  double pdfPtTrueNotNorm(double ptTrue, const double *par) const;
  double underlyingPdfPtTrue(double ptTrue, const double *par) const;

  //! Print some initialization details
  void print() const;
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
  
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedJetEt(const Measurement *x,const double *par) const {
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

  double inverseJetCorrection(const Measurement *x,const double *par) const { 
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

//!  \brief Parametrization of jet response based to eta-eta moment
//!
//!  This parametrization has 14 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class EtaEtaParametrization: public Parametrization {
public:
  EtaEtaParametrization() : Parametrization(0,14,0,0) {}
  const char* name() const { return "EtaEtaParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    if(std::abs(x->eta) > 1.2) return x->pt;
    double y = x->etaeta > 0.25 ? 0.12 : x->etaeta - 0.13;
    
    double c = 1/(par[0] + exp((x->pt-par[1])/par[2])) + par[3];
    c += y * ( par[4] - par[5]/(pow(log(x->pt),par[6])+par[7]) + par[8]/x->pt);
    c += y *y *(par[9] * x->pt + par[10] * pow(x->pt,2.0/3.0) + par[11] * pow(x->pt,-1.0/3.0) + par[12] * (x->pt - par[13])/x->pt);
    
    return c * x->pt;  
  }
};

//!  \brief Parametrization of jet response based to eta-eta moment
//!
//!  This parametrization has 14 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class PhiPhiParametrization: public Parametrization {
public:
  PhiPhiParametrization() : Parametrization(0,14,0,0) {}
  const char* name() const { return "PhiPhiParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    if(std::abs(x->eta) > 1.2) return x->pt;
    double y = x->phiphi > 0.25 ? 0.12 : x->phiphi - 0.13;
    double pt =  (x->pt < 20) ? 20 : ((x->pt > 1000) ? 1000 : x->pt);
    double c = 1/(par[0] + exp((pt-par[1])/par[2])) + par[3];
    c += y * ( par[4] - par[5]/(pow(log(pt),par[6])+par[7]) + par[8]/pt);
    c += y *y *(par[9] * pt + par[10] * pow(pt,2.0/3.0) + par[11] * pow(pt,-1.0/3.0) + par[12] * (pt - par[13])/pt);
    // std::cout << x->pt << ":" << c << ", " << 1/c << std::endl;
    return 1/c * x->pt;  
  }
};


//!  \brief Parametrization of jet response based on emf
//!
//!  This parametrization has 40 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class BinnedEMFParametrization: public Parametrization {
public:
  BinnedEMFParametrization() : Parametrization(0,40,0,0) {}
  const char* name() const { return "BinnedEMFParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    if(std::abs(x->eta) > 1.2) return x->pt;
    // please note that Measurement::EMF is EmEt!!!
    double emf = x->EMF / (x->EMF+x->HadF);
    if(emf < 0) return x->pt;
    if(emf > 1) return x->pt;
    int bin = (int)(emf * 10);
    if(bin == 10) bin = 9;
    
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    double fX[3],fY[3];
    if((bin > 0) && (bin < 9)) {
      for(int i = 0 ; i < 3 ; ++i) {
	int id = 4 * (bin + i - 1);
	fX[i] = (bin + i - 1) / 10.0 + 0.05;
	fY[i] = par[id] + logpt *(0.1*par[id+1] + logpt * (0.001*par[id+2]+ 0.001*par[id+3]*logpt));
	//std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
      }
      double c = quadraticInterpolation(emf,fX,fY);
      if(c != c) { 
	for(int i = 0 ; i < 3 ; ++i) {
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
	std::cout << "interpolated:" << emf << ":" << c << '\n';
      }
      return c * x->pt;	
    } 
    int id = 4 * bin;
    double c = par[id] + logpt *(0.1*par[id+1] + logpt * (0.001*par[id+2]+ 0.001*par[id+3]*logpt));
    return c * x->pt;
  }
};

//!  \brief Parametrization of jet response based on width in phi
//!
//!  This parametrization has 40 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class BinnedPhiPhiParametrization: public Parametrization {
public:
  BinnedPhiPhiParametrization() : Parametrization(0,50,0,0) {}
  const char* name() const { return "BinnedPhiPhiParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    if(std::abs(x->eta) > 1.2) return x->pt;
    
    int bin = (int)(x->phiphi * 40 );
    if(bin > 9) bin = 9;
    if(bin < 0) {
      std::cout << "error: wrong bin " << bin << " for phiphi: " << x->phiphi 
		<< std::endl;
      assert(bin >= 0);
    }
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);

    double fX[3],fY[3];
    if((bin > 0) && (bin < 9)) {
      for(int i = 0 ; i < 3 ; ++i) {
	int id = 5 * (bin + i - 1);
	fX[i] = (bin + i - 1) / 40.0;
	fY[i] = par[id] + logpt *(0.1*par[id+1] + logpt * (0.001*par[id+2]+ logpt * (0.001*par[id+3] + logpt * 0.001*par[id+4])));
	//std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
      }
      double c = quadraticInterpolation(x->phiphi,fX,fY);
      if(c != c) { 
	for(int i = 0 ; i < 3 ; ++i) {
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
	std::cout << "interpolated:" << x->phiphi << ":" << c << '\n';
      }
      return c * x->pt;	
    } 
    int id = 5 * bin;
    double c = par[id] + logpt *(0.1*par[id+1] + logpt * (0.001*par[id+2]+ logpt * (0.001*par[id+3] + logpt * 0.001 * par[id+4])));
    return c * x->pt;
  }  
};


class BinnedPhiPhiParametrization2 : public Parametrization {
public:
  BinnedPhiPhiParametrization2() : Parametrization(0,40,0,0) {}
  const char* name() const { return "BinnedPhiPhiParametrization2";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    //if(std::abs(x->eta) > 3.0) return x->pt;
    
    int bin = (int)((x->phiphi-0.025) * 45 );
    if(bin > 9) bin = 9;
    if(bin < 0) bin = 0;

    double pt = (x->E < 15.0) ? 15.0 : (x->E > 1200.0) ? 1200.0 : x->pt; 
    double logpt = log10(pt);
    double fX[3],fY[3];
    
    if((bin > 0) && (bin < 9)) {
      for(int i = 0 ; i < 3 ; ++i) {
	int id = 4 * (bin + i - 1);
	fX[i] = (bin + i - 1) / 40.0;
	fY[i] = par[id] + 0.1 * logpt * ( par[id+1] + 0.1 * logpt * ( par[id+2] + 0.1 * par[id+3] * logpt));
	if(std::isnan(fY[i])) { 
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << ", " << pt << ", " << logpt << '\n';
	}
      }
      double c = quadraticInterpolation(x->phiphi,fX,fY);
      if(c != c) { 
	for(int i = 0 ; i < 3 ; ++i) {
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
	std::cout << "interpolated:" << x->phiphi << ":" << c << '\n';
      }
      return 1/c * x->pt;	
    } 
    int id = 4 * bin;
    double c = par[id] + 0.1 * logpt * ( par[id+1] + 0.1 * logpt * ( par[id+2] + 0.1 * par[id+3] * logpt));
    if(c <= 0)  {
      std::cout << "bad factor:" << c << " from " <<  par[id] << ", " 
		<< par[id+1] << ", " << par[id+2] << ", " << par[id+3] 
		<< " and pt " << pt << std::endl; 
    }
    return 1/c * x->pt;
  }  
};

//!  \brief Parametrization of jet response based on width in phi
//!
//!  This parametrization has 10 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class BinnedScaledPhiPhiParametrization: public Parametrization {
public:
  BinnedScaledPhiPhiParametrization() : Parametrization(0,45,0,0) {}
  const char* name() const { return "BinnedScaledPhiPhiParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    if(std::abs(x->eta) > 1.2) return x->pt;

    //hack use original Jet pt for scaling
    double scaled = x->phiphi /(0.2 - 0.02 * log(x->EMF + x->HadF));
    int bin = (int)(scaled * 10);
    if(bin > 14) bin = 14;
    if(bin < 0) bin = 0;
    
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 1200.0) ? 1200.0 : x->pt; 
    double logpt = log10(pt);
    double fX[3],fY[3];
    if((bin > 0) && (bin < 9)) {
      for(int i = 0 ; i < 3 ; ++i) {
	int id = 3 * (bin + i - 1);
	fX[i] = (bin + i - 1) / 20.0 + 0.2;
	fY[i] = par[id] + 0.1 * logpt * ( par[id+1] + 0.1 * logpt );
	if(std::isnan(fY[i])) { 
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
      }
      double c = quadraticInterpolation(scaled,fX,fY);
      if(c != c) { 
	for(int i = 0 ; i < 3 ; ++i) {
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
	std::cout << "interpolated:" << scaled << ":" << c << '\n';
      }
      return 1/c * x->pt;	
    } 
    int id = 3 * bin;
    double c = par[id] + 0.1 * logpt * ( par[id+1] + 0.1 * logpt);
    if(c <= 0)  {
      std::cout << "bad factor:" << c << " from " <<  par[bin] << std::endl; 
    }
    return 1/c * x->pt;
  }
};


//!  \brief Parametrization of jet response based on width in eta
//!
//!  This parametrization has 10 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class BinnedScaledEtaEtaParametrization: public Parametrization {
public:
  BinnedScaledEtaEtaParametrization() : Parametrization(0,10,0,0) {}
  const char* name() const { return "BinnedScaledEtaEtaParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    if(std::abs(x->eta) > 1.2) return x->pt;

    //hack use original Jet pt for scaling
    double scaled = x->etaeta * log(x->EMF + x->HadF);
    int bin = (int)((scaled - 0.2) * 20);
    if(bin > 9) bin = 9;
    if(bin < 0) bin = 0;
    
    double fX[3],fY[3];
    if((bin > 0) && (bin < 9)) {
      for(int i = 0 ; i < 3 ; ++i) {
	int id = bin + i - 1;
	fX[i] = (bin + i - 1) / 20.0 + 0.2;
	fY[i] = par[id];
	if(std::isnan(fY[i])) { 
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
      }
      double c = quadraticInterpolation(scaled,fX,fY);
      if(c != c) { 
	for(int i = 0 ; i < 3 ; ++i) {
	  std::cout << i << ":" << fX[i] << ", " << fY[i] << '\n';
	}
	std::cout << "interpolated:" << scaled << ":" << c << '\n';
      }
      return c * x->pt;	
    } 
    
    double c = par[bin];
    if(c <= 0)  {
      std::cout << "bad factor:" << c << " from " <<  par[bin] << std::endl; 
    }
    return c * x->pt;
  }
};


//!  \brief Parametrization of jet response based on width in phi
//!
//!  This parametrization has 5 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class SimplePhiPhiParametrization : public Parametrization {
public:
  SimplePhiPhiParametrization() : Parametrization(0,10,0,0) {}
  const char* name() const { return "SimplePhiPhiParametrization";}
    
  double correctedTowerEt(const Measurement *x,const double *par) const {
    return x->pt;
  }
 
  double correctedJetEt(const Measurement *x,const double *par) const {
    //if(std::abs(x->eta) > 1.2) return x->pt;
    
    double cut = 30 * x->pt / x->E;
    double pt = (x->pt < cut) ? cut : (x->pt > 800.0) ? 800.0 : x->pt; 
    double logpt = log10(pt);
    double phi = x->phiphi;
    //phi *= (1 + (par[0] + (0.1 * par[1] + (0.01 * par[2] + 0.001 * par[3] * logpt) * logpt) * logpt) * 
    double phimean =  0.1 * par[8] - 0.01 * par[9] * logpt;
    //std::cout << par[8] << ", " << par[9] << ", " <<  phimean << '\n';
    if(phimean > 0.2) phimean = 0.2;
    if(phimean < 0.05) phimean = 0.05;
    phi -= phimean;
    //if(phi < -0.15) phi = -0.15;
    //else if(phi > 0.3) phi = 0.3;
    //phi =- 0.213 - 0.0486 * logpt;
    double a = par[0] + logpt * ( par[1] + logpt * (par[2] + logpt * par[3]));
    double b = par[4] + logpt * ( par[5] + logpt * (par[6] + logpt * par[7])); 
    double c = 1 + (0.1 * a + 0.1 * b * phi) * phi;
    //if(c <= 0) {
    //  std::cout << pt << ", " << phi << '\n';
    //}
    if(c < 0.3) c = 0.3;
    if(c > 3.0) c = 3.0;
    assert(c > 0);
    return x->pt / c;
  }  
};

#endif
