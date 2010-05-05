//
//    Class for jets with tracks 
//
//    first version: Hartmut Stadie 2009/04/08
//    $Id: JetWithTracks.cc,v 1.11 2010/02/15 12:40:18 stadie Exp $
//   
#include"JetWithTracks.h"

#include "TLorentzVector.h"

JetWithTracks::JetWithTracks(double Et, double EmEt, double HadEt ,double OutEt, double E,
			     double eta,double phi, double phiphi, double etaeta, Flavor flavor, double genPt, double dR,
			     CorFactors* corFactors, const Function& func, 
			     double (*errfunc)(const double *x, const Measurement *xorig, double err), 
			     const Function& gfunc, double Etmin) 
  :  Jet(Et,EmEt,HadEt,OutEt,E,eta,phi,phiphi,etaeta,flavor,genPt,dR,corFactors,
	 func,errfunc,gfunc,Etmin),
     ntrackpars(0), expectedCaloEt(0), trackPt(0)
{
}

JetWithTracks::JetWithTracks(const JetWithTracks& j) 
  :  Jet(j),ntrackpars(0), expectedCaloEt(0), trackPt(0)
{
  for(TrackCollConstIter i = j.tracks.begin() ; i != j.tracks.end() ; ++i) {
    const Track *t = *i;
    addTrack(t->Et(),t->EmEt(),t->HadEt(),t->OutEt(),t->E(),t->eta(),t->phi(),
	     t->TrackId,t->TowerId,t->DR,t->DRout,t->etaOut,t->phiOut,
	     t->EM1,t->EM5,t->Had1,t->Had5,t->TrackChi2,t->NValidHits,
	     t->TrackQualityT,t->MuDR,t->MuDE,t->Efficiency,t->f,t->errf);
  }
}

JetWithTracks::~JetWithTracks() 
{
  for(TrackCollIter i = tracks.begin() ; i != tracks.end() ; ++i) {
    delete *i;
  }
}  

void JetWithTracks::ChangeParAddress(double* oldpar, double* newpar) 
{
  Jet::ChangeParAddress(oldpar,newpar);
  for(TrackCollIter i = tracks.begin() ; i != tracks.end() ; ++i) {
    (*i)->ChangeParAddress(oldpar,newpar);
  }
  for(std::map<int,double*>::iterator iter = trackpars.begin() ;
      iter != trackpars.end() ; ++iter) {
    iter->second += newpar - oldpar;
  }
}

double JetWithTracks::correctedEt(double Et,bool fast) const
{
  const double ConeRadius = 0.5; 
  const double MIPsignal = 4;   
  int NoUsedTracks = 0;
  bool isMuon;
  trackPt = 0;
  expectedCaloEt = 0;
  for(TrackCollConstIter i = tracks.begin() ; i != tracks.end() ; ++i) {
    Track *t = (*i);
    if(! t->goodTrack()) continue; 
    if( t->pt > 100) continue;

    ++NoUsedTracks;
    isMuon = ( t->trackId() == 13);
  
    if(t->dR() < ConeRadius) {
      if(t->dRout() < ConeRadius) {
	if(!isMuon) {
	  trackPt += t->pt;
	  expectedCaloEt += t->expectedEt();
	}
      } else { // out-of-cone
	trackPt += t->pt;
      }
    } else { // subtract curl-ins
      if(isMuon) expectedCaloEt -= MIPsignal;
      else expectedCaloEt -= t->expectedEt();
    }
  }
  //correct neutral rest
  double res = (Et > expectedCaloEt) ? Jet::correctedEt(Et -  expectedCaloEt) : 0;
  //std::cout << "Et:" << Et << " track corrected:" << ccet << " jet cor:" << res << '\n';
  return res + trackPt;
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTracks::varyPars(double eps, double Et, double start)
{
  Jet::varyPars(eps,Et,start);
  int i = Jet::nPar();

  for(std::map<int,double*>::const_iterator iter = trackpars.begin() ;
      iter != trackpars.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " track[0]:" << tracks[0]->Par() << '\n';
    for(int trkpar = 0 ; trkpar < ntrackpars ; ++trkpar) {
      //std::cout <<  "truth: " << Et << " alternating par:" << id + trkpar << "  = " << p[trkpar] << std::endl;
      double orig = p[trkpar]; 
      p[trkpar] += eps;
      varcoll[i].upperEt = Jet::expectedEt(Et,start,varcoll[i].upperError);
      p[trkpar] = orig - eps;
      varcoll[i].lowerEt = Jet::expectedEt(Et,start,varcoll[i].lowerError); 
      p[trkpar] = orig;
      varcoll[i].parid = id + trkpar;
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTracks::varyParsDirectly(double eps, bool computeDeriv)
{
  Jet::varyParsDirectly(eps,computeDeriv);
  int i = Jet::nPar();

  const double deltaE = eps * 100.0;
  for(std::map<int,double*>::const_iterator iter = trackpars.begin() ;
      iter != trackpars.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " track[0]:" << tracks[0]->Par() << '\n';
    for(int trkpar = 0 ; trkpar < ntrackpars ; ++trkpar) {
      //std::cout <<  " alternating par:" << id + trkpar << "  = " << p[trkpar] 
      //		<< ", " << p << std::endl;
      double orig = p[trkpar]; 
      p[trkpar] += eps;
      varcoll[i].upperEt = correctedEt(Measurement::pt);
      varcoll[i].upperError = expectedError(varcoll[i].upperEt);
      if(computeDeriv) {
	varcoll[i].upperEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }
      p[trkpar] = orig - eps;
      varcoll[i].lowerEt = correctedEt(Measurement::pt); 
      varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
      if(computeDeriv) {
	varcoll[i].lowerEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }      
      p[trkpar] = orig;
      varcoll[i].parid = id + trkpar;
      //std::cout << "up:" << varcoll[i].upperEt << " low:" << varcoll[i].lowerEt << '\n'; 
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}

double JetWithTracks::Error() const {
  double var = 0, err;
  for(TrackCollConstIter i = tracks.begin() ; i != tracks.end() ; ++i) {
    err = (*i)->Error();
    var += err * err;
  }   
  double jeterr = Jet::Error();
  return sqrt(var + jeterr * jeterr);
}

double JetWithTracks::expectedError(double et) const
{
  double var = 0, err;
  for(TrackCollConstIter i = tracks.begin() ; i != tracks.end() ; ++i) {
    err = (*i)->Error();
    var += err * err;
  }
  double jeterr = Jet::expectedError(et);
  return sqrt(var + jeterr * jeterr);
}
  
void JetWithTracks::addTrack(double Et, double EmEt, double HadEt ,
			     double OutEt, double E,double eta,double phi,
			     int TrackId, int TowerId, double DR, double DRout,
			     double etaOut, double phiOut, double EM1, double EM5, 
			     double Had1, double Had5, double TrackChi2, 
			     int NValidHits, bool TrackQualityT, double MuDR, 
			     double MuDE, double Efficiency, const Function& func,
			     double (*errfunc)(const double *x, const Measurement *xorig, double err))
{
  tracks.push_back(new Track(Et,EmEt,HadEt,OutEt,E,eta,phi,TrackId,TowerId,DR,DRout,etaOut,phiOut,
			     EM1,EM5,Had1,Had5,TrackChi2,NValidHits,TrackQualityT,MuDR,MuDE,Efficiency,
			     func,errfunc)); 
  ntrackpars = func.nPars();
  trackpars[func.parIndex()] = func.firstPar();
  varcoll.resize(Jet::nPar() + trackpars.size() * ntrackpars);
}



JetWithTracks::Track::Track(double Et, double EmEt, double HadEt ,
			    double OutEt, double E,double eta,double phi,
			    int TrackId, int TowerId, double DR, double DRout,
			    double etaOut, double phiOut, double EM1, double EM5, 
			    double Had1, double Had5, double TrackChi2, 
			    int NValidHits, bool TrackQualityT, double MuDR, 
			    double MuDE, double Efficiency, const Function& func,
			    double (*errfunc)(const double *x, const Measurement *xorig, double err))
  :  TTrack(Et,EmEt,HadEt,OutEt,E,eta,phi,TrackId,TowerId,DR,DRout,etaOut,phiOut,EM1,EM5,Had1,Had5,
	    TrackChi2,NValidHits,TrackQualityT,MuDR,MuDE,Efficiency), f(func), errf(errfunc)
{ 
}

double JetWithTracks::Track::expectedEt() const
{
  double et = f(this);
  assert(et == et);
  assert(et >= 0);
  return et;
}

  
double JetWithTracks::expectedEt(double truth, double start, bool fast)
{
  correctedEt(start);
  //std::cout << truth << ", " << pt << ", " << cet << ", " 
  //	    <<  expectedCaloEt << ", " << trackPt << "\n";
  if(truth < trackPt) return (expectedCaloEt > 0) ? expectedCaloEt : 1.0;
  if(expectedCaloEt >=  Measurement::pt) return expectedCaloEt;
  double m = Jet::expectedEt(truth, start, fast);
  //std::cout << "expected ET:" << m << '\n';
  //assert(m > 0);
  return m;
}
