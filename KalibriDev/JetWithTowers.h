#ifndef JETWITHTOWERS_H
#define JETWITHTOWERS_H

#include"Jet.h"

#include <vector>
#include <map>


//!
//!    \brief Class for jets with towers 
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/25
//!    $Id: JetWithTowers.h,v 1.21 2010/02/15 12:40:18 stadie Exp $
// ----------------------------------------------------------------   
class JetWithTowers : public Jet
{
 public:
  JetWithTowers(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi, double phiphi, double etaeta, 
		Flavor flavor, double genPt, double dR, CorFactors* corFactors,
		const Function& f,
		double (*errfunc)(const double *x, const Measurement *xorig, double err), 
		const Function& gf, double Etmin = 0); 
  virtual ~JetWithTowers(); 
  virtual int nPar() const {return Jet::nPar() + towerpars.size() * ntowerpars;}
  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double correctedEt(double Et,bool fast = false) const; 
  virtual double Error() const;
  virtual double expectedError(double et) const;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(double eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(double eps, bool computeDeriv = false);

  void addTower(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi,const Function& f,
		double (*errfunc)(const double *x, const Measurement *xorig, double err));
  virtual Jet* clone() const { return new JetWithTowers(*this);} //!< Clone this jet
 private:
  JetWithTowers(const JetWithTowers& j); //!< disallow copies!
  class Tower : public Measurement {
  public:
    Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double eta,double phi, double alpha, const Function& func,
	  double (*errfunc)(const double *x, const Measurement *xorig, double err));
    Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double EmEttrue, double HadEttrue, double OutEttrue,
	  double eta,double phi, double alpha, const Function& func,
	  double (*errfunc)(const double *x, const Measurement *xorig, double err));
    virtual ~Tower() {}
    double Et()     const {return pt;}
    double EmEt()   const {return EMF;}
    double HadEt()  const {return HadF;}
    double OutEt()  const {return OutF;}
    double E()      const {return Measurement::E;}
    double eta()    const {return Measurement::eta;}
    double phi()    const {return Measurement::phi;}
    double projectionToJetAxis() const {return alpha;}
    double fractionOfJetHadEt() const { return fraction;}
    void setFractionOfJetHadEt(double frac) const { fraction = frac;}
    void ChangeParAddress(double* oldpar, double* newpar) {f.changeParBase(oldpar,newpar);}
    double correctedHadEt(double HadEt) const;
    double lastCorrectedHadEt() const { return lastCorHadEt;}  
    double Error() const {return errf(&(Measurement::pt),this,0);}
    double expectedError(double et) const { return  errf(&et,this,0);}
    int nPar() const {return f.nPars();}
    int FirstPar() const {return f.parIndex();}
    double *Par() const {return f.firstPar();}
  private:
    double alpha;                //!< Projection factor onto jet axis
    double error;                //!< Error for constant error mode
    double mEttrue;              //!< True total transverse energy
    double mEmEttrue;            //!< True Et from the ECAL part of the tower		
    double mHadEttrue;           //!< True Et from the HCAL part of the towers
    double mOutEttrue;           //!< True Et from the HO part of the tower
    mutable Measurement temp;
    mutable double lastCorHadEt;
    mutable double fraction;
    Function f;
    double (*errf)(const double *x, const Measurement *xorig, double err);

    friend class JetWithTowers;
  };
  typedef std::vector<Tower*> TowerColl;
  typedef TowerColl::iterator TowerCollIter;
  typedef TowerColl::const_iterator TowerCollConstIter;
  TowerColl towers;
  int ntowerpars;
  std::map<int,double*> towerpars;
};

#endif
