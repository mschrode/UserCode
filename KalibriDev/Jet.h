//!    \brief Class for basic jets 
//!
//!    \author Hartmut Stadie
//!
//!    \date 2008/12/14
//!
//!    $Id: Jet.h,v 1.35 2010/04/13 13:44:10 mschrode Exp $
#ifndef JET_H
#define JET_H

#include"CalibData.h"
#include "Function.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
     
class CorFactors;

class Jet : public Measurement
{
 public:
  //! For sorting jets in calo pt
  static bool caloPtGreaterThan(const Jet *j1, const Jet *j2) {
    // check for 0
    if (j1 == 0) {
      return j2 != 0;
    } else if (j2 == 0) {
      return false;
    } else {
      return j1->pt() > j2->pt();
    }
  }


  //!  \brief Jet flavor
  //!
  //!  The possible flavors are
  //!  - 0: Gluon
  //!  - 1: u, d, or s quark
  //!  - 2: c quark
  //!  - 3: b quark
  enum Flavor{ unknown = -1, gluon=0, uds=1, c=2, b=3 };

  //! return flavor for pdg id
  static Flavor flavorFromPDG(int pdg) {
    if(pdg == 21) return gluon;
    if(pdg == 0) return unknown;
    unsigned int id = std::abs(pdg);
    if(id < 4) return uds;
    if(id == 4) return c;
    if(id == 5) return b;
    return unknown;
  }

 public:
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, double phiphi, double etaeta, Flavor flavor, 
      double genPt, double dR, CorFactors* corFactors, const Function& f,
      double (*errfunc)(const double *x, const Measurement *xorig, double err), 
      const Function& gf, double Etmin = 0); 
  virtual ~Jet();

  double Et()     const {return Measurement::pt;}                 //!< Return transverse energy Et
  double pt()     const {return Measurement::pt;}                 //!< Return transverse energy Et
  double EmEt()   const {return EMF;}                //!< Return Et from the ECAL part of the towers
  double HadEt()  const {return HadF;}               //!< Return Et from the HCAL part of the towers
  double OutEt()  const {return OutF;}               //!< Return Et from the HOut part of the towers
  double E()      const {return Measurement::E;}    //!< Return energy
  double eta()    const {return Measurement::eta;}  //!< Return pseudorapidity
  double phi()    const {return Measurement::phi;}  //!< Return azimuthal angle
  double momentPhiPhi() const {return Measurement::phiphi;}  //!< Return phi-phi moment (width of jet in phi)
  double momentEtaEta() const {return Measurement::etaeta;}  //!< Return eta-eta moment (width of jet in eta)
  double scaledPhiWidth() const { return (sphi_ > 0 ? sphi_ : sphi_ = momentPhiPhi() /( 0.2 - 0.02 * log(Et())));}
  double scaledEtaWidth() const { return (seta_ > 0 ? seta_ : seta_ = momentEtaEta() /(0.2 -  0.02 * log(Et())));}  
  Flavor flavor() const {return flavor_;}       //!< Return jet flavor
  double genPt()  const {return genPt_;}        //!< Return Pt for corresponding GenJet 
  double dR() const {return dR_;}               //!< \f$ \Delta R \f$ between jet and genjet
  const CorFactors& corFactors() const { return *corFactors_;}
  void updateCorFactors(CorFactors *cor);
  //! Correct measurement by product \p L1*L2*L3
  void correctToL3();
  //! Correct measurement by product \p L2*L3
  void correctL2L3();

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
  //!  errf and Measurement::pt, the pt of the original
  //!  measurement.
  //!  
  //!  \return Error of original measurement
  // ---------------------------------------------------------
  virtual double Error() const {return errf(&(Measurement::pt),this,0);}

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
  virtual const VariationColl& varyParsDirectly(double eps, bool computeDeriv = true);
  
  void print();                       //!< Print some jet members
  static void printInversionStats();  //!< Print some info on inversion

  int parIndex() const { return f.parIndex(); }

  virtual Jet* clone() const { return new Jet(*this);} //!< Clone this jet
  void setGlobalFunction(const Function& ngf) { gf = ngf;} //!< Set global correction function, needed for constraints

 protected:
  mutable VariationColl varcoll;
  virtual double expectedEt(double truth, double start, bool fast = false);
  Jet(const Jet&j); //!< disallow copies!

 private: 
  Flavor flavor_;           //!< The jet's Flavor
  double genPt_;            //!< The genjet pt
  double dR_;               //!< \f$ \Delta R \f$ between jet and genjet
  const CorFactors* corFactors_;   //!< The correction factors
  double    error;                //!< Stores error for constant error mode
  Function  f;                    //!< Jet correction function
  Function  gf;                   //!< Global jet correction function
  double    (*errf)(const double *x, const Measurement *xorig, double err);   //!< Error function
  double    etmin;                //!< Lower cut on measured Et

  bool      secant(double truth, double& x1, double& x2, double eps);
  bool      falseposition(double truth, double& x1, double& x2, double eps);

  static long long ncalls;        //!< Number of calls of inversion methods secant and falseposition
  static long long ntries;        //!< Number of tries in iteration during inversion
  static long long nfails;        //!< Number of failed tries during inversion
  static long long nwarns;        //!< Number of warnings during inversion

  mutable Measurement temp;
  mutable double root;
  mutable double sphi_,seta_;
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
