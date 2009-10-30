#ifndef JETWITHTRACKS_H
#define JETWITHTRACKS_H

#include"Jet.h"

#include <vector>
#include <map>


//!
//!    \brief Class for jets with towers 
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/25
//!    $Id: JetWithTracks.h,v 1.4 2009/07/13 12:04:39 snaumann Exp $
// ---------------------------------------------------------------   
class JetWithTracks : public Jet
{
 public:
  JetWithTracks(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi, Flavor flavor,const Function& f, 
		double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
		const Function& gf, double Etmin = 0); 
  JetWithTracks(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi, Flavor flavor,double genPt, double dR,
		TJet::CorFactors corFactors, const Function& f,
		double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
		const Function& gf, double Etmin = 0); 
  virtual ~JetWithTracks(); 
  virtual int nPar() const {return Jet::nPar() + trackpars.size() * ntrackpars;}
  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double correctedEt(double Et,bool fast = false) const; 
  virtual double Error() const;
  virtual double expectedError(double et) const;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(double eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(double eps);

  void addTrack(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi,int TrackId, int TowerId, double DR, double DRout,
		double etaOut, double phiOut, double EM1, double EM5, double Had1, 
		double Had5, double TrackChi2, int NValidHits, bool TrackQualityT, 
		double MuDR, double MuDE, double Efficiency, const Function& f,
		double (*errfunc)(const double *x, const TMeasurement *xorig, double err));
 protected:
  virtual double expectedEt(double truth, double start, bool fast = false);
 private:
  class Track : public TTrack {
  public:
    Track(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double eta,double phi, int TrackId, int TowerId, double DR, double DRout,
	  double etaOut, double phiOut, double EM1, double EM5, double Had1, double Had5,
	  double TrackChi2, int NValidHits, bool TrackQualityT, double MuDR, double MuDE,
	  double Efficiency, const Function& func,
	  double (*errfunc)(const double *x, const TMeasurement *xorig, double err));
    virtual ~Track() {}
    double Et()     const {return pt;}
    double EmEt()   const {return EMF;}
    double HadEt()  const {return HadF;}
    double OutEt()  const {return OutF;}
    double E()      const {return TMeasurement::E;}
    double eta()    const {return TMeasurement::eta;}
    double phi()    const {return TMeasurement::phi;}
    int trackId() const { return TrackId;}
    int towerId() const { return TowerId;}
    double dR() const { return DR;}
    double dRout() const { return DRout;}
    bool goodTrack() const { return TrackQualityT;}
    void ChangeParAddress(double* oldpar, double* newpar) {f.changeParBase(oldpar,newpar);}
    double expectedEt() const;
    double Error() const {return 0;}
    int nPar() const {return f.nPars();}
    int FirstPar() const {return f.parIndex();}
    double *Par() const {return f.firstPar();}
  private:
    Function f;
    double (*errf)(const double *x, const TMeasurement *xorig, double err);
  };
  typedef std::vector<Track*> TrackColl;
  typedef TrackColl::iterator TrackCollIter;
  typedef TrackColl::const_iterator TrackCollConstIter;
  TrackColl tracks;
  int ntrackpars;
  std::map<int,double*> trackpars;
  mutable double expectedCaloEt;
  mutable double trackPt;
  
};

#endif
