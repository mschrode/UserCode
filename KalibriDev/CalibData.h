//  \author Christian Autermann
//  \date Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
//
#ifndef CalibData_h
#define CalibData_h

#include <iostream>
using namespace std;
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
//!  \sa TJet, TTower, TTrack, Jet, JetWithTowers
//!
//!  \author Christian Autermann
//!  $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
class TMeasurement
{
public:
 TMeasurement():pt(0.),EMF(0.),HadF(0.),OutF(0.),E(0.),eta(0.),phi(0.){};  //!< All quantities are initialized with 0
 TMeasurement(double Et,double EmEt,double HadEt,double OutEt,double E,
	       double eta,double phi)
    : pt(Et),EMF(EmEt),HadF(HadEt),OutF(OutEt),E(E),eta(eta),phi(phi) {}
 TMeasurement(TMeasurement* m):pt(m->pt),EMF(m->EMF),HadF(m->HadF),OutF(m->OutF),
                                E(m->E),eta(m->eta),phi(m->phi){};
  virtual ~TMeasurement() {};
  //all common variables
  double pt;     //!< Total transverse momentum (pt = EMF + HadF + OutF)
  double EMF;    //!< Pt from the ECAL part of the tower(s)		
  double HadF;   //!< Pt from the HCAL part of the towers(s)		
  double OutF;   //!< Pt fromt the HO part of the tower(s)		
  double E;      //!< Total energy					
  double eta;    //!< Pseudorapidity eta				
  double phi;    //!< Polar angle phi  
};




//!  \brief A tower measurement
//!
//!  \sa TMeasurement, TJet, TTrack, Jet, JetWithTowers
//!
//!  \author Christian Autermann
//!  $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
class TTower : public TMeasurement
{ 
public:
  TTower():TMeasurement(){}; 
  TTower(double Et,double EmEt,double HadEt,double OutEt,double E,
	 double eta,double phi)
    : TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi) {}
  TTower(TMeasurement* t):TMeasurement(t){};
  //TTower(TTower* t):TMeasurement(t){/*further initialization*/};
  virtual ~TTower() {}
//variables specific only to towers (i.e. # EM cells)
};




//!  \brief A jet measurement
//!
//!  \sa TMeasurement, TTower, TTrack, Jet, JetWithTowers
//!
//!  \author Christian Autermann
//!  $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
class TJet : public TMeasurement
{
public:
  //!  \brief   Container class for jet correction factors
  class CorFactors
  {
  public :
    CorFactors(double L1=1.0, double L2=1.0, double L3=1.0, double L4=1.0, double L5=1.0,
	       double JPT=1.0, double JPTL2L3=1.0) :
      l1_(L1), l2_(L2), l3_(L3), l4_(L4), l5_(L5), jpt_(JPT), jptL2L3_(JPTL2L3) {};
    double getL1()  const { return l1_; }    //!< Return L1 correction factor (zero-suppression)
    double getL2()  const { return l2_; }    //!< Return L2 correction factor (relative, in eta)
    double getL3()  const { return l3_; }    //!< Return L3 correction factor (absolute, in pt)
    double getL4()  const { return l4_; }    //!< Return L4 correction factor (electromagnetic fraction)
    double getL5()  const { return l5_; }    //!< Return L5 correction factor (flavor)
    double getJPT() const { return jpt_; }   //!< Return Jet+Track correction factor
    double getL2L3() const { return l2_*l3_; }   //!< Return product of L2 and L3 correction factors
    double getJPTL2L3() const { return jptL2L3_; }   //!< Return product of L2 and L3 correction factors for Jet+Track
    double getToL2() const { return l1_*l2_; }         //!< Return factor needed to get L2 corrected from raw jets: L1*L2
    double getToL3() const { return getToL2()*l3_; }   //!< Return factor needed to get L3 corrected from raw jets: L1*L2*L3
    double getToL4() const { return getToL3()*l4_; }   //!< Return factor needed to get L4 corrected from raw jets: L1*L2*L3*L4
    double getToL5() const { return getToL4()*l5_; }   //!< Return factor needed to get L5 corrected from raw jets: L1*L2*L3*L4*L5
    double getToJPTL3()
      const { return jpt_*l1_*jptL2L3_; }   //!< Return factor needed to get L3 corrected from raw jets for JPT: JPT*L1*JPTL2L3
  private :
    double l1_;      //!< Level 1 correction factor (zero-suppression)
    double l2_;      //!< Level 2 correction factor (relative, in eta)
    double l3_;      //!< Level 3 correction factor (absolute, in pt)
    double l4_;      //!< Level 4 correction factor (electromagnetic fraction)
    double l5_;      //!< Level 5 correction factor (flavor)
    double jpt_;     //!< Jet+Track correction factor
    double jptL2L3_; //!< Product of level 2 and level 3 correction factors for Jet+Track
  };
  //!  \brief Jet flavor
  //!
  //!  The possible flavors are
  //!  - 0: Gluon
  //!  - 1: u, d, or s quark
  //!  - 2: c quark
  //!  - 3: b quark
  enum Flavor{ gluon=0, uds=1, c=2, b=3 };
  TJet() : TMeasurement(), flavor(gluon), genPt(0.), dR(0.), ptHat(0.) {}; 
  TJet(double Et,double EmEt,double HadEt,double OutEt,double E,double eta,
       double phi, Flavor flavor, double genPt, double dR,
       CorFactors corFactors)
    : TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi),
    flavor(flavor), genPt(genPt), dR(dR), ptHat(0.), corFactors(corFactors) {};
  TJet(TMeasurement* j):TMeasurement(j){};
  virtual ~TJet() {}
  
  Flavor flavor;           //!< The jet's Flavor
  double genPt;            //!< The genjet pt
  double dR;               //!< \f$ \Delta R \f$ between jet and genjet
  double ptHat;            //!< \f$ \hat{p}_{T} \f$ of the event
  CorFactors corFactors;   //!< The correction factors
};



//!  \brief A track measurement
//!
//!  \sa TMeasurement, TJet, TTower, Jet, JetWithTowers
//!
//!  \todo Document members
//!
//!  \author Jan Thomsen
//!  $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
class TTrack : public TMeasurement
{
public:
  TTrack():TMeasurement(){};
  TTrack(double Et, double EmEt, double HadEt ,double OutEt, double E,double eta,
	 double phi,int TrackId, int TowerId, double DR, double DRout, 
	 double etaOut, double phiOut, double EM1, double EM5, double Had1, 
	 double Had5, double TrackChi2, int NValidHits, bool TrackQualityT, 
	 double MuDR, double MuDE, double Efficiency) 
    : TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi),TrackId(TrackId),TowerId(TowerId),
    DR(DR),DRout(DRout),etaOut(etaOut),phiOut(phiOut),EM1(EM1),EM5(EM5),Had1(Had1),
    Had5(Had5),TrackChi2(TrackChi2),NValidHits(NValidHits),TrackQualityT(TrackQualityT),
    MuDR(MuDR),MuDE(MuDE),Efficiency(Efficiency) {}
  TTrack(TMeasurement* tr):TMeasurement(tr){};
  //TTrack(TTrack* tr):TMeasurement(tr){/*further initialization*/};
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
//!  A TData object represents one event. It holds the measured
//!  quantities of that event (see TMeasurement) and allows
//!  access to the corrected measurement. Moreover, the normalized,
//!  weighted, squared, and squared residual \f$ z^{2} \f$ of this
//!  event, which enters the global \f$ \chi^{2} = \sum z^{2} \f$
//!  function, is calculated.
//!
//!  TData is a virtual base class. The derived interfaces are
//!  specific for a certain type of data.
//!
//!  There are currently two different calibration schemes, resulting
//!  in two different sets of data classes derived from TData:
//!  -# Original calibration scheme ("Correction of the measurement")
//!     There is a second base class for this scheme, TAbstractData,
//!     deriving from TData. All interfaces for specific data types
//!     derive from TAbstractData in this scheme.
//!  -# New calibration scheme ("Variation of the truth") 
//!     The available data types are:
//!  \author Christian Autermann
//!  \date Wed Jul 18 13:54:50 CEST 2007
//! $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
class TData
{
public:
  virtual ~TData() {}
  virtual TMeasurement *GetMess() const = 0;                           //!< Get TMeasurement object
  virtual double GetTruth() const = 0;                                 //!< Get truth of measurement
  virtual double GetParametrizedMess() const = 0;                      //!< Get corrected measurement
  virtual void ChangeParAddress(double* oldpar, double* newpar) = 0;   //!< Change adress of parameter array
  virtual DataType GetType() const = 0;                                //!< Get DataType
  virtual double GetWeight() const = 0;                                //!< Get weight
  virtual void setWeight(double w) = 0;                              //!< Set weight
  virtual double ptHat() const = 0;                                    //!< Get event scale


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
  //!  discards events with $|residual| > 1.5 \sigma$
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double ScaleTukey(double z2) { return ( z2 > 2.25) ? 0 : z2;}    //!< Scaling of residual a la  Tukey
};




//!  \brief Interface to the data for the original calibration
//!         scheme ("Correction of the measurement")
//!
//!  For a description of the functionality see TData
//!
//!  TAbstractData is a virtual base class. The derived interfaces are
//!  specific for a certain type of data.
//!
//!  \author Hartmut Stadie
//!  \date Thu Dec 11 17:20:25 2008 UTC
//!  $Id: CalibData.h,v 1.73 2009/10/30 08:14:24 mschrode Exp $
class TAbstractData : public TData
{
public:
  //!  \brief Constructor (default)
  TAbstractData() : TData() {_par=0;};


  //!  \brief Constructor
  //!
  //!  \param index  Index of the first of the successive parameters
  //!                covered by this event, see TParameters
  //!  \param mess   Pointer to the measurement, see TMeasurement
  //!  \param truth  Truth
  //!  \param error  Error on measurement
  //!  \param weight Weight of event in \f$ \chi^{2} \f$ sum
  //!  \param par    Pointer to the first of the successive elements in
  //!                parameter array covered by this event, see TParameters
  //!  \param n_par  Number of succesive parameters covered by this event,
  //!                see TParameters
  //!  \param *func  Pointer to correction function, see Parametrization
  //!  \param *err   Pointer to error function, see TParameters
  TAbstractData(unsigned short int index, TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double (*func)(const TMeasurement*, const double*),
	double (*err)(const double *,const TMeasurement *,double))
  : _index(index), _mess(mess),_truth(truth),_error(error),_weight(weight),_par(par),_n_par(n_par),_func(func),_err(err){};


  //!  \brief Destructor
  virtual ~TAbstractData(){
    delete _mess;
  };

  TMeasurement *GetMess() const { return _mess;}; //!< Get TMeasurement object


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
  virtual double GetParametrizedMess() const { return _func(_mess,_par); };  



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
  //!  in a TData_TruthMultMess object. If the tower energies 
  //!  are corrected with different parameters, a different jet
  //!  energy (sum of tower energies) enters the jet correction 
  //!  function. This can be done using this method.
  //!
  //!  \param paramess Array of measurements which are to be corrected
  //!  \return Corrected measurement
  //!
  //!  \sa GetParametrizedMess()
  virtual double GetParametrizedMess(double *const paramess) const {
    TMeasurement m(_mess);
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
  TMeasurement *_mess;
  double _truth, _error, _weight;
  double *_par;
  double _invisError;
  unsigned short int _n_par; //limited from 0 to 65535
  double (*_func)(const TMeasurement *x, const double *par);
  double (*_err)(const double *x, const TMeasurement *x_original, double error);
  DataType _type;
};

//data class for data providing one truth and one messurement, 
//e.g. track-tower
class TData_TruthMess : public TAbstractData
{
public:
  TData_TruthMess(unsigned short int index,  TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par, double (*func)(const TMeasurement *,const double*),  double (*err)(const double*,const TMeasurement*,double))
  : TAbstractData(index, mess, truth, error, weight, par, n_par, func, err ){_type=TrackTower;_invisError=0;};

  virtual const std::vector<TAbstractData*>& GetRef() { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };
  
  //only a dummy
  virtual const std::vector<TAbstractData*>& GetRefTrack()  { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };


  virtual double chi2() const{ 
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return GetWeight() * (*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  
private:
  static std::vector<TAbstractData*> resultcache;
};


//data class for data providing one truth and multiple messurements, 
//e.g. gamma-jet or track-cluster
class TData_TruthMultMess : public TData_TruthMess
{
public:
  TData_TruthMultMess(unsigned short int index, double truth, double error, 
		      double weight, double * par, unsigned short int n_par,
        	      double (*func)(const TMeasurement *,const double*),
		      double (*err)(const double*,const TMeasurement*,double),
		      TMeasurement *mess)
  : TData_TruthMess(index,  mess, truth, error, weight, par, n_par, func, err){
    _type=GammaJet; trackuse = false;_invisError=0;};
  virtual ~TData_TruthMultMess() {
    for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it)
      delete *it;
    _vecmess.clear();	
  };
  void AddMess(TData_TruthMess * m){ _vecmess.push_back(m);};
  void AddTrack(TData_TruthMess * m){ _vectrack.push_back(m);};

  void UseTracks(bool use){   //check if tracks shall be used in this jet:
    if((fabs(_mess->eta) < 2.1) &&  ( _vectrack.size() > 0) && use)      //@ eta > 2.1or eta < -2.1 parts of the cone are outside the tracker. 
      trackuse = true;
  };
  virtual bool GetTrackuse(){return trackuse;};

  virtual void updateError() {
    if(trackuse)
      {
	double CaloRest = _mess->pt;
	bool IsMuon;
	double ConeRadius = 0.5;
	for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	     it!=_vectrack.end(); ++it) {
	  //if( ((TTrack*)(*it)->GetMess())->TrackChi2 > 6  || ((TTrack*)(*it)->GetMess())->NValidHits < 9) // || TrackQuality != 1)
	  if( !((TTrack*)(*it)->GetMess())->TrackQualityT)
	    continue;  // qualityTracks = false;

	  if( (*it)->GetMess()->pt > 100) continue;

	  TTrack* temp = (TTrack*)(*it)->GetMess();
	  if(temp->TrackId == 13)  IsMuon = true;
	  else                     IsMuon = false;
	  if(temp->DR < ConeRadius) 
	    {
	      if(temp->DRout < ConeRadius)
		{
		  if(!IsMuon) {
		    CaloRest -= (*it)->GetParametrizedMess();   //this depends on start values. Error will be updated, but look for a more stable way  
		  }
		}
	    }
	}
	std::vector<TAbstractData*>::const_iterator tower=_vecmess.begin();
	if(CaloRest > 0)
	  _invisError = (*tower)->GetParametrizedErr(&CaloRest);  //Rest Error same as tower error, save this invisible energy error
      }
  }

  virtual double GetParametrizedMess() const{
    double result = 0;
    if(trackuse) 
      { 
	result = GetParametrizedTrackMess();
      }
    else
      {    //no tracks are used
	double tower_pt_sum=0.0;
	for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	     it!=_vecmess.end(); ++it){
	  tower_pt_sum += (*it)->GetParametrizedMess(); // Sum of tower Pt
	}
	TJet jet(_mess);
	jet.pt = tower_pt_sum;
	result = _func(&jet, _par);
      }
    return result;
  };

  virtual double GetParametrizedErr() const { //returns total jet error and not only the non-tower part like _err
    double error = 0, new_mess=0, sum_mess=0, new_error=0,sum_error2=0;
    if(trackuse) 
      {
	for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	     it!=_vectrack.end(); ++it) {

	  if( !((TTrack*)(*it)->GetMess())->TrackQualityT)  continue;
	  if( (*it)->GetMess()->pt > 100) continue;

	  new_mess =  (*it)->GetMess()->pt;
	  new_error   = (*it)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
	}
	sum_error2 += _invisError * _invisError;
      }
    else
      {
      for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	     it!=_vecmess.end(); ++it) {
	  new_mess    = (*it)->GetParametrizedMess();
	  sum_mess   += new_mess;
	  new_error   = (*it)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
	}
	TJet jet(_mess);
	jet.pt    = sum_mess;
	new_mess  = _func(&jet, _par);
	new_error =  _err( &new_mess, _mess, _error );
	sum_error2 += new_error * new_error;
	}
    error = sqrt(sum_error2);
    /*
    TJet* myTJet =  (TJet*)(GetMess());
    double pt =myTJet->genPt;
    double error = 5.6 + 1.25 * sqrt(pt) + 0.033 * pt;
    */
    return error;
  };

  virtual double GetParametrizedMess(double *const paramess) const { // For derivative calculation
    double result = 0;
    if(trackuse) 	result = GetParametrizedTrackMess();
    else {
      TJet jet(_mess);
      jet.pt = paramess[0];
      result =  _func(&jet,_par);
    }
    return result;
  };

  virtual double GetParametrizedTrackMess() const{
    double JetPt = 0, CaloRest, CaloTrackPt;
    int NoUsedTracks = 0;
    bool IsMuon;
    const double ConeRadius = 0.5;      //Jet Cone Radius should not be hard coded or must be changed if Radius != 0.5
    const double MIPsignal = 4;      
    CaloRest = _mess->pt;
    TJet* myTJet =  (TJet*)(GetMess());
    CaloRest *= myTJet->corFactors.getL1();   //ZSP correction

    /*
    //possible Tower correction here (instead of ZSP?)
    double tower_pt_sum=0.0;
    for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it){
      tower_pt_sum += (*it)->GetParametrizedMess(); // Sum of tower Pt
    }
    tower_pt_sum *= myTJet->ZSPcor;   //ZSP correction
    CaloRest =  tower_pt_sum;
    */


    for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	 it!=_vectrack.end(); ++it) {
      if( !((TTrack*)(*it)->GetMess())->TrackQualityT)
	continue; 

      if( (*it)->GetMess()->pt > 100) continue;

      NoUsedTracks++;
      TTrack* temp = (TTrack*)(*it)->GetMess();

      CaloTrackPt = (*it)->GetParametrizedMess();  //Expected Signal of track in Calorimeter

      if(temp->TrackId == 13)  IsMuon = true;
      else                     IsMuon = false;
      if(temp->DR < ConeRadius) 
	{
	  if(temp->DRout < ConeRadius)
	    {
	      if(!IsMuon) {
		JetPt +=  temp->pt;
		CaloRest -= CaloTrackPt;
		//JetPt += ((1 - temp->Efficiency)/temp->Efficiency) * (temp->pt - CaloTrackPt);  //track reco correction
	      }
	    }
	  else //Out of Cone
	    {
	      JetPt +=  temp->pt;
	    }
	}
      else
	{
	  if(IsMuon)     CaloRest -= MIPsignal;
	  else           CaloRest -= CaloTrackPt;
	}
    }

    if(CaloRest > 0) {
      TJet jet(_mess);
      jet.pt = CaloRest;
      jet.E = -1000;  //signal that this is only the calo rest of a track jet!
      JetPt += _func(&jet, _par);  //ein anderer Parameter als bei der normalen korrektur nehmen (da hier nur neutr. Had + other rest)
    }
    TJet jet(_mess);
    jet.pt = JetPt;
    jet.E = -800; //signal that this is the factor for a track jet!
    JetPt = _func(&jet, _par);  // factor like JetMET on top of JPT

    return JetPt;
  };


  virtual double chi2() const{ 
    double weight = GetWeight(); 
    double new_mess, new_error;
    double sum_mess = 0.;
    double sum_error2 = 0.;
    if(trackuse) {
      new_error = GetParametrizedErr();
      new_mess = GetParametrizedTrackMess();
    }
    else{
      for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	   it!=_vecmess.end(); ++it) {
	new_mess    = (*it)->GetParametrizedMess();
	sum_mess   += new_mess;
	new_error   = (*it)->GetParametrizedErr(&new_mess);
	sum_error2 += new_error * new_error;
      }
      TJet jet(_mess);
      jet.pt    = sum_mess;
      new_mess  = _func(&jet, _par);
      new_error =  _err( &new_mess, _mess, _error );
    }
    return (new_error!=0. ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) ):0.0);
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  virtual const std::vector<TAbstractData*>& GetRef() {return _vecmess;};
  virtual const std::vector<TAbstractData*>& GetRefTrack() {return _vectrack;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { 
    TAbstractData::ChangeParAddress(oldpar,newpar);
    for (std::vector<TAbstractData*>::iterator it=_vecmess.begin();
	 it !=_vecmess.end(); ++it) { (*it)->ChangeParAddress(oldpar,newpar);} 
    if(trackuse) {
      for (std::vector<TAbstractData*>::iterator it = _vectrack.begin() ;
	   it !=_vectrack.end() ; ++it) { (*it)->ChangeParAddress(oldpar,newpar);} 
    }
  }
protected:  
  std::vector<TAbstractData*> _vecmess; 
  std::vector<TAbstractData*> _vectrack; 
  bool trackuse;
  //double JetError2;
};

//virtual data class for data providing only multiple messurements
//NOT DIRECTLY USED, SINCE "combine" FUNCTION IS NOT YET DEFINED HERE!
//Base class for PT-balance and Inv-Mass !!!
class TData_MessMess : public TData_TruthMultMess
{
public:
  TData_MessMess(unsigned short int index, double * dir, double truth, double error, 
                 double weight, double * par, unsigned short int n_par,
       	         double (*func)(const TMeasurement *,const double*),
		 double (*err)(const double*,const TMeasurement*,double),
		 TMeasurement *mess)
  : TData_TruthMultMess(index, truth, error, weight, par, n_par, func, err,mess),
    _direction(dir){_type=MessMess; trackuse = false;_invisError=0;};
  virtual ~TData_MessMess() {
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it)
      delete *it;
    _m2.clear();
    delete [] _direction;	
  };
  virtual void AddNewMultMess(TData_MessMess * m2 ){assert(m2->_m2.empty());_m2.push_back(m2);};
  void ClearMultMess() { _m2.clear();}
  unsigned MultMessSize(){return _m2.size();};
  virtual double chi2() const{ 
    double sum_error2=0.0, new_error, new_mess;
    double weight = GetWeight();
    
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
         it!=_m2.end(); ++it) {
      new_error   = (*it)->GetParametrizedErr();
      sum_error2 += new_error * new_error;
    }
    new_error = GetParametrizedErr();
    new_mess  = GetMessCombination();

    return (sum_error2!=0 ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) ) : 0.0) ;
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  virtual void ChangeParAddress(double* oldpar, double* newpar) { 
    TData_TruthMultMess::ChangeParAddress(oldpar,newpar);
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin(); it!=_m2.end(); ++it)
      (*it)->ChangeParAddress(oldpar,newpar);
  }
  virtual double GetMessCombination()const{ return combine(); }; // for plotting
  virtual double * GetDirection() const{ return _direction; };
  virtual std::vector<TData_MessMess*>const * GetSecondaryJets() const{return &_m2;};
   
protected:
  virtual double combine() const {return 0.;};
  std::vector<TData_MessMess*> _m2;
  double * _direction;
};

//e.g. jet-jet, jet-jet-jet, etc..
class TData_PtBalance : public TData_MessMess
{
public:
  TData_PtBalance(unsigned short int index, double * dir, double truth, double error, 
		  double weight, double * par, unsigned short int n_par,
        	  double (*func)(const TMeasurement *,const double*),
		  double (*err)(const double*,const TMeasurement*,double),
		  TMeasurement *mess)
  : TData_MessMess(index, dir, truth, error, weight, par, n_par, func, err, mess){_type=PtBalance; trackuse = false;_invisError=0;};
  virtual ~TData_PtBalance(){};
  virtual double GetScale() const{
    double scale = GetMess()->pt;
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();it!=_m2.end(); ++it)
      scale += (*it)->GetMess()->pt;
    return scale/2.;
  }
  virtual void updateError() {
    if(trackuse) //update error on invisible energy 1st jet
      {
	double CaloRest = _mess->pt;
	bool IsMuon;
	double ConeRadius = 0.5;
	for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	     it!=_vectrack.end(); ++it) {
	  //if( ((TTrack*)(*it)->GetMess())->TrackChi2 > 6  || ((TTrack*)(*it)->GetMess())->NValidHits < 9) // || TrackQuality != 1)
	  if( !((TTrack*)(*it)->GetMess())->TrackQualityT)
	    continue;  // qualityTracks = false;

	  if( (*it)->GetMess()->pt > 100) continue;

	  TTrack* temp = (TTrack*)(*it)->GetMess();
	  if(temp->TrackId == 13)  IsMuon = true;
	  else                     IsMuon = false;
	  if(temp->DR < ConeRadius) 
	    {
	      if(temp->DRout < ConeRadius)
		{
		  if(!IsMuon) {
		    CaloRest -= (*it)->GetParametrizedMess();   //this depends on start values. Error will be updated, but look for a more stable way  
		  }
		}
	    }
	}
	std::vector<TAbstractData*>::const_iterator tower=_vecmess.begin();
	if(CaloRest > 0)
	  _invisError = (*tower)->GetParametrizedErr(&CaloRest);  //Rest Error same as tower error, save this invisible energy error
	else 
	  _invisError = 0;
      }

    double totalsum = GetParametrizedMess();
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      (*it)->updateError();  //update error on invisible energy other jets
      totalsum += (*it)->GetParametrizedMess();
    }
    
    double scale = totalsum / 2; //should be changed to ptsum (scalar) of the part projected to leading jet axis devided by 2 for n>2 n-jets     (also chi2() and chi2fast())

    double new_mess = GetParametrizedMess();
    double sum = totalsum - new_mess; 
    double new_error = GetParametrizedErr();
    new_error *= sum / (scale * scale);
    double sum_error2 = new_error * new_error;
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      new_mess    = (*it)->GetParametrizedMess();
      new_error   = (*it)->GetParametrizedErr();
      sum = totalsum - new_mess;
      new_error *= sum / (scale * scale);
      sum_error2 += new_error * new_error;
    }
    if( (GetParametrizedMess()>0) && (totalsum>GetParametrizedMess()))   //update only if new value is reasonable (Pt > 0)
    _error = sqrt(sum_error2);  
  };

  virtual double chi2() const{ 
    double new_mess;
    double weight = GetWeight();
      
    double parascale = GetParametrizedMess();
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      parascale += (*it)->GetParametrizedMess();
    }
    
    //should be changed to ptsum (scalar) of the part projected to leading jet axis devided by 2 for n>2 n-jets           same for parascale and in ChiFast...
    parascale /= 2;

    /*
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
         it!=_m2.end(); ++it) {
      //new_mess    = (*it)->GetParametrizedMess();
      new_mess    = (*it)->GetMess()->pt;	 //const error
      new_error   = (*it)->GetParametrizedErr(&new_mess);
      sum = totalsum -  (*it)->GetMess()->pt;   //new_mess?
      new_error *= sum / (scale * scale);
      sum_error2 += new_error * new_error;
    }
    new_mess  = GetMess()->pt;          //const error
    sum = totalsum - GetMess()->pt;   //new_mess?
    new_error = GetParametrizedErr( &new_mess );
    new_error *= sum / (scale * scale);   //
    */
    new_mess  = combine() / parascale;

    return (_error!=0 ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(_error * _error) ) : 0.0) ;
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
protected:
  virtual double combine() const;  
};

//e.g. top-Mass, W-mass, etc..
class TData_InvMass2 : public TData_MessMess
{
public:
  TData_InvMass2(unsigned short int index, double * dir, double truth, double error, double weight,
		double * par, unsigned short int n_par,
        	double (*func)(const TMeasurement *,const double*),
		double (*err)(const double*,const TMeasurement*,double),
		TMeasurement *mess)
  : TData_MessMess(index, dir, truth, error, weight, par, n_par, func, err, mess) { _type=InvMass; trackuse = false;_invisError=0; };
  virtual ~TData_InvMass2(){};
 protected:
  virtual double combine() const;  
};

//!  \brief Data class to limit a parameter
class TData_ParLimit : public TAbstractData
{
 public:
  TData_ParLimit(unsigned short int index, TMeasurement *mess, double error,
                 double *par,double (*func)(const TMeasurement *,const double*))
  : TAbstractData(index,mess,0,error,1.0,par,1,func,0){ _type=ParLimit;_invisError=0;};
    
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
