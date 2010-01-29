//
// $Id: CalibData.h,v 1.2 2010/01/21 16:48:40 mschrode Exp $
//
#ifndef CalibData_h
#define CalibData_h

#include <iostream>
#include <vector> 
#include <cmath>
#include <cassert>


//!  \brief Type of data
//!
//!  \sa TAbstractData 
enum DataType {Default, TrackTower, GammaJet, TrackCluster, MessMess, PtBalance,
               InvMass, typeTowerConstraint, ParLimit, TypeSmearPhotonJet, TypeSmearDiJet, JetConstraint };

//!  \brief Base class of a measurement
//!
//!  A measurement can represent a tower, a track, or a jet.
//!
//!  \note: The parametrized (pt-)measurement, wich will be compared to the
//!      'truth' will remain a single 'double' value (i.e. the 
//!       same type as 'truth')!
//!
//!  \sa Jet, Tower, Track, JetWithTowers, JetWithTracks
//!
//!  \author Christian Autermann
//!  $Id: CalibData.h,v 1.2 2010/01/21 16:48:40 mschrode Exp $
class Measurement
{
public:
 Measurement() :
  pt(0.),EMF(0.),HadF(0.),OutF(0.),E(0.),eta(0.),phi(0.),etaeta(0.0)
    {
    }
 Measurement(double Et,double EmEt,double HadEt,double OutEt,double E,
	     double eta,double phi, double etaeta = 0)
   : pt(Et),EMF(EmEt),HadF(HadEt),OutF(OutEt),E(E),eta(eta),phi(phi),
    etaeta(etaeta) 
  {
  }
  virtual ~Measurement() {};
  //all common variables
  double pt;     //!< Total transverse momentum (pt = EMF + HadF + OutF)
  double EMF;    //!< Pt from the ECAL part of the tower(s)		
  double HadF;   //!< Pt from the HCAL part of the towers(s)		
  double OutF;   //!< Pt fromt the HO part of the tower(s)		
  double E;      //!< Total energy					
  double eta;    //!< Pseudorapidity eta				
  double phi;    //!< Polar angle phi  
  double etaeta; //!< Eta-Eta moment (width in eta) 
};




//!  \brief A track measurement
//!
//!  \sa Measurement, TJet, TTower, Jet, JetWithTowers
//!
//!  \todo Document members
//!
//!  \author Jan Thomsen
//!  $Id: CalibData.h,v 1.2 2010/01/21 16:48:40 mschrode Exp $
class TTrack : public Measurement
{
public:
  TTrack():Measurement(){};
  TTrack(double Et, double EmEt, double HadEt ,double OutEt, double E,double eta,
	 double phi,int TrackId, int TowerId, double DR, double DRout, 
	 double etaOut, double phiOut, double EM1, double EM5, double Had1, 
	 double Had5, double TrackChi2, int NValidHits, bool TrackQualityT, 
	 double MuDR, double MuDE, double Efficiency) 
    : Measurement(Et,EmEt,HadEt,OutEt,E,eta,phi),TrackId(TrackId),TowerId(TowerId),
    DR(DR),DRout(DRout),etaOut(etaOut),phiOut(phiOut),EM1(EM1),EM5(EM5),Had1(Had1),
    Had5(Had5),TrackChi2(TrackChi2),NValidHits(NValidHits),TrackQualityT(TrackQualityT),
    MuDR(MuDR),MuDE(MuDE),Efficiency(Efficiency) {}
  virtual ~TTrack() {}
//variables specific only to Tracks
  int TrackId;
  int TowerId;
  double DR;
  double DRout;
  double etaOut;
  double phiOut;
  double EM1;
  double EM5;
  double Had1;
  double Had5;
  double TrackChi2;
  int NValidHits;
  bool TrackQualityT;
  double MuDR;
  double MuDE;
  double Efficiency;
};


//!  \brief Interface to the data 
//!
//!  A Event object represents one event. It holds the measured
//!  quantities of that event (see Measurement) and allows
//!  access to the corrected measurement. Moreover, the normalized,
//!  weighted, squared, and squared residual \f$ z^{2} \f$ of this
//!  event, which enters the global \f$ \chi^{2} = \sum z^{2} \f$
//!  function, is calculated.
//!
//!  Event is a virtual base class. The derived interfaces are
//!  specific for a certain type of data.
//!
//!  There are currently two different calibration schemes, resulting
//!  in two different sets of data classes derived from Event:
//!  -# Original calibration scheme ("Correction of the measurement")
//!     There is a second base class for this scheme, TAbstractData,
//!     deriving from Event. All interfaces for specific data types
//!     derive from TAbstractData in this scheme.
//!  -# New calibration scheme ("Variation of the truth") 
//!     The available data types are:
//!  \author Christian Autermann
//!  \date Wed Jul 18 13:54:50 CEST 2007
//! $Id: CalibData.h,v 1.2 2010/01/21 16:48:40 mschrode Exp $
class Event
{
public:
  virtual ~Event() {}
  virtual Measurement *GetMess() const = 0;                           //!< Get Measurement object
  virtual double GetTruth() const = 0;                                 //!< Get truth of measurement
  virtual double GetParametrizedMess() const = 0;                      //!< Get corrected measurement
  virtual void ChangeParAddress(double* oldpar, double* newpar) = 0;   //!< Change adress of parameter array
  virtual DataType GetType() const = 0;                                //!< Get DataType
  virtual double GetWeight() const = 0;                                //!< Get weight
  virtual void setWeight(double w) = 0;                              //!< Set weight
  double ptHat() const { return ptHat_; }                      //!< Get event scale


  //!  \brief Get the normalized, squared residual \f$ z^{2} \f$ of this event
  //!
  //!  The normalized, squared residual \f$ z^{2} \f$ of
  //!  this event is calculated. It is weighted with
  //!  GetWeight(), and scaled with ScaleResidual.
  //!  It enters the global \f$ \chi^{2} = \sum z^{2} \f$
  //!  function.
  //!
  //!  \return The normalized, squared residual\f$ z^{2} \f$ of this event
  virtual double chi2() const = 0;


  //!  \brief Chi2 value from last iteration
  //!  \return Chi2 value from last iteration
  // ------------------------------------------
  virtual double chi2_plots() const = 0;


  //!  \brief Get the normalized, squared residual\f$ z^{2} \f$ of this event
  //!         and calculate the first and second derivatives
  //!
  //!  The normalized, squared residual \f$ z^{2} \f$ of
  //!  this event is calculated. It is weighted with
  //!  GetWeight(), and scaled with ScaleResidual.
  //!  It enters the global \f$ \chi^{2} = \sum z^{2} \f$
  //!  function.
  //!
  //!  Moreover, the contribution of this event to the 
  //!  first and second derivative ('temp_derivative1',
  //!  'temp_derivative2') of the global \f$ \chi^{2} \f$
  //!  function is calculated numerically and returned
  //!  by reference, where
  //!  \f[
  //!   \frac{\partial \chi^{2} }{\partial p}
  //!   = \sum \frac{\textrm{temp\_derivative1}}{2\epsilon}
  //!  \f]
  //!  and
  //!  \f[
  //!   \frac{\partial^{2} \chi^{2} }{\partial p^{2}}
  //!   = \sum \frac{\textrm{temp\_derivative2}}{\epsilon^{2}}
  //!  \f]
  //!
  //!  \param temp_derivative1 Pointer to first derivative contribution
  //!  \param temp_derivative2 Pointer to second derivative contribution
  //!  \param epsilon Step size \f$ \epsilon \f$  for derivative calculation
  //!  \return The normalized, squared residual\f$ z^{2} \f$ of this event
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;


  virtual void updateError() = 0;  //!< Update error terms using current corrected energies


  //!  \brief Scale residual for outlier treatment
  //!
  //!  Points to one of the following functions to
  //!  scale the squared, normalized, and weighted
  //!  residual
  //!  \f$ z^{2} = \chi^{2}/\textrm{weight} \f$:
  //!   - ScaleNone(double z2)
  //!   - ScaleCauchy(double z2)
  //!   - ScaleHuber(double z2)
  //!   - ScaleTukey(double z2)
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double (*ScaleResidual)(double z2);


  //!  \brief No scaling of residuals
  //!
  //!  \note This is the default
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double ScaleNone(double z2){ return z2; }


  static double ScaleCauchy(double z2);  //!< Scaling of residual with Cauchy function
  static double ScaleHuber(double z2);   //!< Scaling of residual with Huber function  
  
  //!  \brief Cut on residuals
  //!
  //!  discards events with \f$ |residual| > 1.5 \sigma \f$
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double ScaleTukey(double z2);  //!< Scaling of residual a la  Tukey

 protected:
  double ptHat_;
};


//!  \brief Interface to the data for the original calibration
//!         scheme ("Correction of the measurement")
//!
//!  For a description of the functionality see Event
//!
//!  TAbstractData is a virtual base class. The derived interfaces are
//!  specific for a certain type of data.
//!
//!  \author Hartmut Stadie
//!  \date Thu Dec 11 17:20:25 2008 UTC
//!  $Id: CalibData.h,v 1.2 2010/01/21 16:48:40 mschrode Exp $
class TAbstractData : public Event
{
public:
  //!  \brief Constructor (default)
  TAbstractData() : Event() {_par=0;};


  //!  \brief Constructor
  //!
  //!  \param index  Index of the first of the successive parameters
  //!                covered by this event, see TParameters
  //!  \param mess   Pointer to the measurement, see Measurement
  //!  \param truth  Truth
  //!  \param error  Error on measurement
  //!  \param weight Weight of event in \f$ \chi^{2} \f$ sum
  //!  \param par    Pointer to the first of the successive elements in
  //!                parameter array covered by this event, see TParameters
  //!  \param n_par  Number of succesive parameters covered by this event,
  //!                see TParameters
  //!  \param *func  Pointer to correction function, see Parametrization
  //!  \param *err   Pointer to error function, see TParameters
  TAbstractData(unsigned short int index, Measurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double (*func)(const Measurement*, const double*),
	double (*err)(const double *,const Measurement *,double))
  : _index(index), _mess(mess),_truth(truth),_error(error),_weight(weight),_par(par),_n_par(n_par),_func(func),_err(err){};


  //!  \brief Destructor
  virtual ~TAbstractData(){
    delete _mess;
  };

  Measurement *GetMess() const { return _mess;}; //!< Get Measurement object


  //!  \brief Get corrected measurement
  //!
  //!  Calculates the corrected measurement from the
  //!  original measurement as returned byGetMess().
  //!  The correction is given by the correction function (see
  //!  Parametrization, TAbstractData) using the parameter
  //!  values currently stored in the parameter
  //!  array in TParameters.
  //!
  //!  \return Corrected measurement
  //!
  //!  \sa GetParametrizedMess(double *const paramess)
  double GetParametrizedMess() const { return _func(_mess,_par); };  

  //!  \brief Correct a user defined measurement
  //!
  //!  Calculates the corrected measurement from a \e specified
  //!  measurement 'paramess'.
  //!  The correction is given by the correction function (see
  //!  Parametrization, TAbstractData) using the parameter
  //!  values currently stored in the parameter
  //!  array in TParameters.
  //!
  //!  \note This method is intended for  derivative calculation
  //!  (see chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) )
  //!  in a Event_TruthMultMess object. If the tower energies 
  //!  are corrected with different parameters, a different jet
  //!  energy (sum of tower energies) enters the jet correction 
  //!  function. This can be done using this method.
  //!
  //!  \param paramess Array of measurements which are to be corrected
  //!  \return Corrected measurement
  //!
  //!  \sa GetParametrizedMess()
  virtual double GetParametrizedMess(double *const paramess) const {
    Measurement m(*_mess);
    m.pt = paramess[0];
    return _func(&m,_par);
  };


  //!  \brief Get error of measurement 
  //!
  //!  Returns the error of the original measurement.
  //!
  //!  \note This method might have a different functionality
  //!  in derived classes. In that case, the error of the
  //!  corrected measurement is calculated using the error
  //!  parametrization function _err.
  //!
  //!  \return Error of original measurement
  virtual double GetParametrizedErr() const { return _error;};


  //!  \brief Get error of a user defined measurement 
  //!
  //!  Returns the error of a \e specified measurement.
  //!  The error of the specified measurement 'paramess'is
  //!  calculated using the error parametrization function
  //!  _err.
  //!
  //!  \param paramess Pointer to measurements for which the error is to be calculated
  //!  \return Error of the measurement
  virtual double GetParametrizedErr(double *const paramess) const { return _err(paramess,_mess,_error);};


  //!  \brief Get error square of a user defined measurement 
  //!
  //!  Returns the squared error of a \e specified measurement.
  //!  This is the same as
  //!  GetParametrizedErr(double *const paramess) * GetParametrizedErr(double *const paramess)
  //!
  //!  \param paramess Pointer to measurements for which the error is to be calculated
  //!  \return Squared error of the measurement
  //!  \sa GetParametrizedErr(double *const paramess)
  virtual double GetParametrizedErr2(double *const paramess){ 
    double error = GetParametrizedErr(paramess);
    return error *error;
  };

  double GetTruth() const { return _truth;};
  virtual double GetScale() const {return GetTruth();};//flatten spectrum w.r.t. this
  virtual void updateError(){};
  double GetError() const { return _error;};
  double GetWeight() const { return _weight;};
  virtual void setWeight(double w) { _weight = w;};
  virtual double ptHat() const { return 0.; }                                    //!< Dummy
  DataType GetType() const {return _type;};
  void SetType(DataType type) {_type=type;};
  unsigned short int GetIndex(){return _index;};
  virtual const std::vector<TAbstractData*>& GetRef() = 0;
  virtual const std::vector<TAbstractData*>& GetRefTrack() = 0;
  virtual double chi2() const {return 0.;};
  double chi2_plots() const { return chi2(); }
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  double * GetPar(){return _par;};
  unsigned short int GetNumberOfPars() const {return _n_par;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { _par += newpar - oldpar;}
  virtual bool GetTrackuse() {return true;}

  static unsigned int total_n_pars;

protected:
  unsigned short int _index; //limited from 0 to 65535
  Measurement *_mess;
  double _truth, _error, _weight;
  double *_par;
  double _invisError;
  unsigned short int _n_par; //limited from 0 to 65535
  double (*_func)(const Measurement *x, const double *par);
  double (*_err)(const double *x, const Measurement *x_original, double error);
  DataType _type;
};


//!  \brief Data class to limit a parameter
class TData_ParLimit : public TAbstractData
{
 public:
 TData_ParLimit(unsigned short int index, Measurement *mess, double error,
		double *par,double (*func)(const Measurement *,const double*))
   : TAbstractData(index,mess,0,error,1.0,par,1,func,0)
    { 
      _type=ParLimit;
      _invisError=0;
    };
  
  virtual const std::vector<TAbstractData*>& GetRef() { 
    _cache.clear();
    _cache.push_back(this);
    return _cache;
  };
  //only a dummy
  virtual const std::vector<TAbstractData*>& GetRefTrack() { 
    _cache.clear();
    _cache.push_back(this);
    return _cache;
  };
  
  virtual double GetParametrizedErr(double *const paramess) const{ 
    return _error;
  };
  
  virtual double chi2() const{  
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return new_mess * new_mess / (new_error * new_error);
  };
  
  virtual double chi2_fast(double* temp_derivative1, 
			   double* temp_derivative2, double const epsilon) const;
  
 private:
  static std::vector<TAbstractData*> _cache;
};

#endif
