#ifndef JET_H
#define JET_H

#include"CalibData.h"
#include "Function.h"


//!    \brief Class for basic jets 
//!
//!    \author Hartmut Stadie
//!
//!    \date 2008/12/14
//!
//!    $Id: Jet.h,v 1.23 2009/10/26 21:00:36 mschrode Exp $


#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
     
class Jet : public TJet
{
 public:
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor, const Function& f, 
      double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
      const Function& gf, double Etmin = 0);
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor, double genPt, double dR,
      TJet::CorFactors corFactors, const Function& f,
      double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
      const Function& gf, double Etmin = 0); 
  virtual ~Jet() {};

  double Et()     const {return pt;}                 //!< Return transverse energy Et
  double EmEt()   const {return EMF;}                //!< Return Et from the ECAL part of the towers
  double HadEt()  const {return HadF;}               //!< Return Et from the HCAL part of the towers
  double OutEt()  const {return OutF;}               //!< Return Et from the HOut part of the towers
  double E()      const {return TMeasurement::E;}    //!< Return energy
  double eta()    const {return TMeasurement::eta;}  //!< Return pseudorapidity
  double phi()    const {return TMeasurement::phi;}  //!< Return azimuthal angle
  Flavor flavor() const {return TJet::flavor;}       //!< Return jet flavor
  double GenPt()  const {return TJet::genPt;}        //!< Return Pt for corresponding GenJet

  //!  \brief Change address of parameters covered by this jet
  //!  \sa Parameters
  // ---------------------------------------------------------
  virtual void ChangeParAddress(double* oldpar, double* newpar) {
    f.changeParBase(oldpar,newpar);
    gf.changeParBase(oldpar,newpar);
  }
  virtual double correctedEt() const { return correctedEt(Et()); }
  virtual double correctedEt(double Et, bool fast = false) const;
  double expectedEt(double truth, double start, double& error,
		    bool fast = false);

  //!  \brief Calculate error from original measurement
  //!
  //!  The error is calculated using the error function
  //!  errf and TMeasurement::pt, the pt of the original
  //!  measurement.
  //!  
  //!  \return Error of original measurement
  // ---------------------------------------------------------
  virtual double Error() const {return errf(&(TMeasurement::pt),this,0);}

  //!  \brief Calculate error from given Et
  //!
  //!  The error is calculated using the error function
  //!  errf and a given Et.
  //!
  //!  \param et Jet Et
  //!  \return Error from jet Et
  // ---------------------------------------------------------
  virtual double expectedError(double et) const { return errf(&et,this,0);}

  //!  \brief Get number of parameters of this jet
  //!
  //!  Returns the number of parameters of the
  //!  local and global jet correction functions
  //!  (see Parametrization).
  //!
  //!  \return Number of parameters of this jet
  // ---------------------------------------------------------
  virtual int nPar() const {return f.nPars() + gf.nPars();}

  //!  \brief Represents a corrected jet if one parameter of
  //!         the correction function is varied (for derivative
  //!         calculation)
  //!
  //!  Stores the corrected Et and error if a parameter is
  //!  varied by +/- eps for derivative calculation.
  // ---------------------------------------------------------
  struct ParameterVariation {
    int    parid;        //!< Id of varied parameter
    double upperEt;      //!< Expected Et if parameter is varied by +eps
    double lowerEt;      //!< Expected Et if parameter is varied by -eps
    double upperError;   //!< Expected error if parameter is varied by +eps
    double lowerError;   //!< Expected error if parameter is varied by -eps
    double upperEtDeriv; //!< Derivative of Et if parameter is  varied by +eps
    double lowerEtDeriv; //!< Derivative of Et if parameter is  varied by +eps
    bool operator==(int b) const { return parid == b;} //!< Two ParameterVariation are the same if they have the same parid
  };
  typedef std::vector<ParameterVariation> VariationColl;
  typedef std::vector<ParameterVariation>::const_iterator VariationCollIter;
  virtual const VariationColl& varyPars(double eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(double eps);

  void print();                       //!< Print some jet members
  static void printInversionStats();  //!< Print some info on inversion

  int parIndex() const { return f.parIndex(); }


 protected:
  mutable VariationColl varcoll;
  virtual double expectedEt(double truth, double start, bool fast = false);

 private:
  double    error;                //!< Stores error for constant error mode
  Function  f;                    //!< Jet correction function
  Function  gf;                   //!< Global jet correction function
  double    (*errf)(const double *x, const TMeasurement *xorig, double err);   //!< Error function
  double    etmin;                //!< Lower cut on measured Et

  bool      secant(double truth, double& x1, double& x2, double eps);
  bool      falseposition(double truth, double& x1, double& x2, double eps);

  static long long ncalls;        //!< Number of calls of inversion methods secant and falseposition
  static long long ntries;        //!< Number of tries in iteration during inversion
  static long long nfails;        //!< Number of failed tries during inversion
  static long long nwarns;        //!< Number of warnings during inversion

  mutable TMeasurement temp;
  const double EoverPt;
  class GslImplementation {
    struct rf_par {
      double y;
      const Jet* jet;
      rf_par(double y, const Jet *jet) : y(y), jet(jet) {}
    } par;
    gsl_root_fsolver *s;
    gsl_function F;
    static double rf(double x, void* params) {
      rf_par* p = (rf_par*)params;
      return p->y - p->jet->correctedEt(x,true);
    };
  public:
    GslImplementation(const Jet* jet);
    ~GslImplementation();
    bool root(double truth, double& x1, double& x2, double eps);
  } gsl_impl;
};

#endif
