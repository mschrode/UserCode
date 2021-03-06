//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.cc,v 1.43 2010/04/13 13:44:10 mschrode Exp $
//   
#include "Jet.h"  

#include "CorFactors.h"

#include <iostream>
#include <iomanip>

Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
         double eta,double phi, double phiphi, double etaeta, Flavor flavor, 
	 double genPt, double dR, CorFactors* corFactors, const Function& f, 
	 double (*errfunc)(const double *x, const Measurement *xorig, double err), 
	 const Function& gf, double Etmin) 
  : Measurement(Et,EmEt,HadEt,OutEt,E,eta,phi,phiphi,etaeta),flavor_(flavor), 
    genPt_(genPt),dR_(dR),corFactors_(corFactors),f(f),gf(gf),
    errf(errfunc),etmin(Etmin),root(Et),EoverPt(E/Et),gsl_impl(this)
{
  temp = *this;
  varcoll.resize(f.nPars() + gf.nPars());
  assert(Measurement::phiphi == Measurement::phiphi);
}

Jet::Jet(const Jet& j) 
  : Measurement(j), flavor_(j.flavor_), genPt_(j.genPt_),dR_(j.dR_), 
    corFactors_(new CorFactors(*(j.corFactors_))),f(j.f),
    gf(j.gf),errf(j.errf),etmin(j.etmin),root(j.root),EoverPt(j.EoverPt),
    gsl_impl(this)
{
  temp = *this;
  varcoll.resize(f.nPars() + gf.nPars());
}

Jet::~Jet()
{
  delete corFactors_;
}

void Jet::updateCorFactors(CorFactors *cor)
{
  delete corFactors_;
  corFactors_ = cor;
}

void Jet::correctToL3()
{
  Measurement::pt   *= corFactors_->getToL3();
  Measurement::E    *= corFactors_->getToL3();
  Measurement::EMF  *= corFactors_->getToL3();
  Measurement::HadF *= corFactors_->getToL3();
  Measurement::OutF *= corFactors_->getToL3();
  temp = *this;
}

void Jet::correctL2L3()
{
  Measurement::pt   *= corFactors_->getL2L3();
  Measurement::E    *= corFactors_->getL2L3();
  Measurement::EMF  *= corFactors_->getL2L3();
  Measurement::HadF *= corFactors_->getL2L3();
  Measurement::OutF *= corFactors_->getL2L3();
  temp = *this;
}

//!  \brief Varies all parameters for this jet by +/-eps
//!
//!  The corrected Et and errors obtained by applying the
//!  correction function with the varied parameters are
//!  stored in a VariationColl.
//!
//!  \param eps Amount by which the parameters are varied
//!  \param Et Et which is to be corrected after parameter variation
//!  \param start Start value for expectedEt(double truth, double start, bool fast)
//!  \return Vector of ParameterVariation
//!  \sa varyParsDirectly
// -------------------------------------------------------
const Jet::VariationColl& Jet::varyPars(double eps, double Et, double start)
{
  //start = Et;
  double oldroot = root;
  for(int i = 0 ; i < f.nPars() ; ++i) {
    double orig = f.firstPar()[i];
    f.firstPar()[i] += eps;
    varcoll[i].upperEt = expectedEt(Et,start,varcoll[i].upperError);
    root = oldroot;
    //if( varcoll[i].upperEt < 0) return varcoll;
    //varcoll[i].upperEt = expectedEt(Et,s,false);
    f.firstPar()[i] = orig - eps;;
    varcoll[i].lowerEt = expectedEt(Et,start,varcoll[i].lowerError); 
    root = oldroot;
    //if( varcoll[i].lowerEt < 0) return varcoll;
    //varcoll[i].lowerEt = expectedEt(Et,s,false);
    f.firstPar()[i] = orig;
    varcoll[i].parid = f.parIndex() + i;
  }
  for(int i = 0,j =  f.nPars(); i < gf.nPars() ; ++i,++j) {
    double orig = gf.firstPar()[i];
    gf.firstPar()[i] += eps;
    varcoll[j].upperEt = expectedEt(Et,start,varcoll[j].upperError);
    root = oldroot;
    //if( varcoll[j].upperEt < 0) return varcoll;
    //varcoll[j].upperEt = expectedEt(Et,s,false);
    gf.firstPar()[i] = orig - eps;;
    varcoll[j].lowerEt = expectedEt(Et,start,varcoll[j].lowerError);
    root = oldroot;
    //if( varcoll[j].lowerEt < 0) return varcoll;
    //varcoll[j].lowerEt = expectedEt(Et,s,false);
    gf.firstPar()[i] = orig;
    varcoll[j].parid = gf.parIndex() + i;
  }
  return varcoll;
}

//!  \brief Varies all parameters for this jet by +/-eps
//!
//!  The corrected original Et and errors obtained by applying the
//!  correction function with the varied parameters are
//!  stored in a VariationColl.
//!
//!  \note In contrast to varyPars, the originally measured
//!        Et is corrected.
//!  \param eps Amount by which the parameters are varied
//!  \return Vector of ParameterVariation
//!  \sa varyPars
// -------------------------------------------------------
const Jet::VariationColl& Jet::varyParsDirectly(double eps, bool computeDeriv)
{
  const double deltaE = eps * 100.0;
  for(int i = 0 ; i < f.nPars() ; ++i) {
    double orig = f.firstPar()[i];
    f.firstPar()[i] += eps;
    varcoll[i].upperEt = correctedEt(Measurement::pt);
    varcoll[i].upperError = expectedError(varcoll[i].upperEt);
    if(computeDeriv) {
      varcoll[i].upperEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
    }
    f.firstPar()[i] = orig - eps;
    varcoll[i].lowerEt = correctedEt(Measurement::pt); 
    varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
    if(computeDeriv) {
      varcoll[i].lowerEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
    }
    f.firstPar()[i] = orig;
    varcoll[i].parid = f.parIndex() + i;
  }  
  for(int i = 0, j =  f.nPars(); i < gf.nPars() ; ++i,++j) {
    double orig = gf.firstPar()[i];
    gf.firstPar()[i] += eps;
    varcoll[j].upperEt = correctedEt(Measurement::pt);
    varcoll[j].upperError = expectedError(varcoll[j].upperEt);
    if(computeDeriv) {
      varcoll[j].upperEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
    }
    gf.firstPar()[i] = orig - eps;
    varcoll[j].lowerEt = correctedEt(Measurement::pt); 
    varcoll[j].lowerError = expectedError(varcoll[j].lowerEt);
    if(computeDeriv) {
      varcoll[j].lowerEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
    }
    gf.firstPar()[i] = orig;
    varcoll[j].parid = gf.parIndex() + i;
  }
  return varcoll;
}


//!  \brief Correct a given jet Et
//!
//!  The given Et is corrected applying successively
//!  the local and the global jet correction functions
//!  f and gf (see also Parametrization).
//!
//!  \note Modifies only hadronic part of tower Et
//!  \param Et Jet Et which gets corrected
//!  \param fast No functionality yet
//!  \return Corrected jet Et
// -------------------------------------------------------
double Jet::correctedEt(double Et, bool fast) const {
  
  //std::cout << "Pars:" << f.firstPar()[0] << ", " << f.firstPar()[1] << ", " << f.firstPar()[2]
  //	    << ", " << gf.firstPar()[0] << ", " << gf.firstPar()[1] << ", " << gf.firstPar()[2]
  //	    << ", " <<  gf.firstPar()[3] << '\n';
  // 
  //assume that only the hadronic energy gets modified!
  temp.pt   = Et;  
  //temp.HadF = Et - OutF - EMF;
  // if(temp.HadF < 0) temp.HadF = 0;
  temp.E    = EoverPt * Et;
  temp.pt = f(&temp);
  /*
    if(corEt != corEt) 
    std::cout << "Et:" << Et << "  orig Et:" << pt << " cor Et:" << corEt << "\n";
    assert(corEt == corEt);
    //if(corEt <  OutF + EMF) corEt = OutF + EMF;
  */
  if(temp.pt <= 0.0) {
    //std::cout << "WARNING: jet cor. Et <= 0.0 GeV:" << corEt << " at eta " << TJet::eta << '\n';
    temp.pt = 1.0;
    temp.E = EoverPt;
  }
  //temp.HadF = corEt - OutF - EMF;
  //if(temp.HadF < 0) temp.HadF = 0;
  temp.E = EoverPt * temp.pt;  
  temp.pt = gf(&temp);
  //if(corEt != corEt) std::cout << "Et:" << Et << "  orig Et:" << pt << " cor Et:" << corEt << "\n";
  //assert(corEt == corEt);
  //if(corEt <  OutF + EMF) corEt = OutF + EMF;
  if(temp.pt <= 1.0) {
    //std::cout << "WARNING: global jet cor. Et <= 1.0 GeV:" << corEt << " at eta " << TJet::eta << '\n';
    temp.pt = 1.0;
    temp.E = EoverPt;
  }
  return temp.pt;
}


//!  \brief Find mean measured Et from correction function for a given truth
//!
//!  Finds the mean measured Et \f$ \bar{E}_{T} \f$ corresponding to a
//!  given truth from the correction function \f$ g_{p} \f$ by solving
//!  \f[
//!    0 = E^{\textrm{true}}_{T} - g_{p}(\bar{E}_{T})
//!  \f]
//!  Here, \f$ p \f$ is the set of current parameters.
//!  The mean Et is defined by the (unknown) response \f$ R \f$ by
//!  \f[
//!   \bar{E}_{T} = R(E^{\textrm{true}}_{T}) \cdot E^{\textrm{true}}_{T}.
//!  \f]
//!  The solution is found numerically using 
//!  secant(double truth, double& x1, double& x2, double eps) or
//!  falseposition(double truth, double& x1, double& x2, double eps)
//!  
//!  \param truth The true Et
//!  \param start Start value for inversion procedure
//!  \param fast No functionality yet
//!  \return Mean measured Et for given truth
// -------------------------------------------------------
double Jet::expectedEt(double truth, double start, bool fast)
{
  if(f.hasInverse()) { 
    temp.pt   = truth;  
    temp.HadF = truth - OutF - EMF;
    if(temp.HadF < 0) temp.HadF = 0;
    temp.E    = Measurement::E * truth/Measurement::pt;
    return f.inverse(&temp);
  }
  static const double eps = 1.0e-12;
  double x1 = root;
  //find root of truth - jet->correctedEt(expectedEt)
  // x: expectedEt
  // y: truth -  jet->correctedEt(expectedEt)
  // f: jet->correctedEt(expectedEt)
  double f1 = correctedEt(x1,false);

  if(std::abs(f1 - truth) < eps * truth) {
    //std::cout << "this really happens!?\n";
    return x1;
  }
 
  //bracket the root
  double x2 = x1 + 0.2;
  
  double f2 = correctedEt(x2,false);
  for( int i = 0 ; ; ++i) {
    if(f1 >= f2) {
      if((f1 > truth) && (f2 < truth)) break; 
      double step = 0.5 * i;
      x1 -= 2 * step;
      x2 += 2 * step;
      f1 = correctedEt(x1,false);
      f2 = correctedEt(x2,false);
    }
    else {
      if(f1 > truth) {
	double step = (f1 - truth) * (x2-x1)/(f2-f1);
	if(step < 0.1) step = 0.1;
	x1 -=  2 * step;
	f1 = correctedEt(x1,false);
	++ntries;
      } else if(f2 < truth) {
	double step = (truth - f2) * (x2-x1)/(f2-f1);
	if(step < 0.1) step = 0.1;
	x2 += 2 * step;
	f2 = correctedEt(x2,false);
	++ntries;
      } else break;
    }
    if(i > 20) {
      ++nfails;
      std::cout << "Warning failed to bag: " << x1 << ", " << x2 << ":" << f1 << " < " << truth << " < " << f2 << std::endl;
      //assert(i < 10);
      return -1;
    }
  }
  //assert((x1 == x1) && (x2 == x2));
  if(! gsl_impl.root(truth,x1,x2,eps)) return -1;
  return root = x2;
  
  //get second point assuming a constant correction factor
  //x2 = (truth - EMF - OutF) * (x1 - EMF - OutF)/(f1 - EMF - OutF) + EMF;
  //x2 = truth * x1 / f1;
  //if(! secant(truth,x2,x1,eps)) return -1;
  //assert(std::abs(correctedEt(x2)-truth) < eps * truth); 
  //root = x2;
  //return x2;
}


//!  \brief Find mean measured Et and error from correction function for a given truth
//!
//!  Finds the mean measured Et and the corresponding error
//!  corresponding to a given truth from the correction function
//!  (see expectedEt(double truth, double start, bool fast)).
//!  The error is calculated from the obtained mean Et using
//!  expectedError(double et) and returned by reference.
//!
//!  Both Et and the error are corrected for effects due to
//!  a cut on the Et spectrum (etmin).
//!  
//!  \param truth The true Et
//!  \param start Start value for inversion procedure
//!  \param error Will be filled with the error
//!  \param fast No functionality yet
//!  \return Mean measured Et for given truth
// -------------------------------------------------------
double Jet::expectedEt(double truth, double start, double& error,bool fast)
{
  //truncate mean for jet min Et-cut
  double m = expectedEt(truth,start,fast);
  if(m < 0) return m;
  double s = expectedError(m);

  // hack
  error = s;
  return m;

//   double x = (etmin - m)/s;
//   if(x < -10) {
//     error = s;
//     return m;
//   }
//   //truncated mean:
//   //m + (E^(-((a - m)^2/(2 s^2))) Sqrt[2/\[Pi]] s)/Erfc[(a - m)/(Sq[2] s)]
//   //truncated RMS
//   //m^2 + s^2 + (E^(-((a - m)^2/(2 s^2))) (a + m) Sqrt[2/\[Pi]] s)/Erfc[(a - m)/(Sqrt[2] s)]
//   double l = exp(-x*x/2) * sqrt(2/M_PI) * s/TMath::Erfc(x/sqrt(2));
//   m += l;
//   s =  expectedError(m);
//   error = sqrt(l*(etmin - m) + s * s);
//   return m;
}


bool Jet::falseposition(double truth, double& x1, double& x2,double eps)
{
  //x2 is the best estimate!
  double f1 = correctedEt(x1,true);
  double f2 = correctedEt(x2,true);
  double temp;
  double step = 0.1 * truth;
  ++ncalls;
  if(x1 > x2) {
    temp = x1;
    x1 = x2;
    x2 = temp;
    temp = f1;
    f1 = f2;
    f2 = temp;
  }  
  double y2 = truth - f2;
  double y1 = truth - f1;
  //std::cout << "x1,2:" << x1 << "," << x2 << " f1,2:" << f1 << ", " << f2 
  //	    << " y1,2:" << y1 << ", " << y2 << std::endl;
  int i = 0;
  while( 1 ) {
    //std::cout << "x1,2:" << x1 << "," << x2 << " f1,2:" << f1 << ", " << f2 
    //      << " y1,2:" << y1 << ", " << y2 << std::endl;
    if(f1 > truth) {
      x1 -= step;
      f1 = correctedEt(x1,true);
      y1 = truth - f1;
      ++ntries;
    } else if(f2 < truth) {
      x2 += step;
      f2 = correctedEt(x2,true);
      y2 = truth - f2;
      ++ntries;
    }
    else break;
    ++i;
    //if(i > 5) break;
    assert(i < 100);
  } 
  i = 0;
  while(std::abs(y2)/truth > eps) {
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << std::endl;
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    double f3 = correctedEt(x3,true);
    double y3 = truth - f3;
    //std::cout << i << ":" << x3 << ":" << y3 << std::endl;
    ++i;
    if(std::abs(y3) / truth < eps) {
      x2 = x3;
      x1 = x1;
      return true;
    }
    if(y3 > 0) {
      x1 = x3;
      f1 = f3;
      y1 = y3;
    } else {
      x2 = x3;
      f2 = f3;
      y2 = y3;
    }
    if(i > 100) {
      ++nfails;
      ntries += i;
      return false;
    }
  } 
  ntries += i;
  return true;
}


bool Jet::secant(double truth, double& x2, double& x1,double eps)
{
  //x2 is the best estimate!
  const double up = 5 * x2;
  const double low = 1.0; 
 
  double f1 = correctedEt(x1,true);
  double f2 = correctedEt(x2,true);
  double y2 = truth - f2;
  double y1 = truth - f1;
  ++ncalls;
  int i = 0;
  double dx = std::abs(x1-x2), dx1 = truth, dx2 = truth;
  if(dx < 1e-12) {
    x2 = 1.0001 * x1;
    dx = std::abs(x1-x2);
  }
  //std::cout << "first intervall size:" << dx/x1 << '\n';
  while((std::abs(y2) > eps * truth)&&(i < 1000)) {
    if(f1 == f2) {
      std::cout << "Warning: no difference in corrected Et:" << f1 << "," << f2 << '\n';
      //print();
      x2 = 0.5 * (x1 + x2);
      ++nfails;
      return false;
    }
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << ", " << x3 << std::endl;
    if(x3 < low) x3 = low;
    if(x3 > up) x3 = up;

    dx2 = dx1;
    dx1 = dx;
    dx = std::abs(x2 - x3);
    if(dx2 < dx) {
      //std::cout << "Warning: fit alternating!\n";
      //std::cout << i << ": last three intervall sizes " << dx << ", " 
      //		<< dx1 << ", " << dx2 << std::endl;
      ++nwarns;
      x3 = 0.5 * (x2+x1);
      x1 = x2;
      x2 = x3;
      dx = truth;
      dx1 = truth;
      dx2 = truth;
      ++i;
      continue;
    }
    double f3 = correctedEt(x3,true);
    //std::cout << x1 << ":" << f1 << " " << x2 << ":" << f2 << " " << x3 << ":" << f3 << '\n';
    double y3 = truth - f3;
    //use false position if root is bracketed
    if(((y3 < 0) && (y1 > 0)) || ((y3 > 0) && (y1 < 0))) {
      x2 = x3;
      f2 = f3;
      y2 = y3;
    } else {
      x1 = x2;
      f1 = f2;
      y1 = y2;
      x2 = x3;
      f2 = f3;
      y2 = y3;
    }
    ++i;
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f1 << "; " << truth - f2 << "\n";
    //     if(i > 100) {
    //       //std::cout << "failed to find good root\n";
    //       //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
    //       x2 = 0.5 * (x2+x1);
    //       ++nfails;
    //       ntries += i;
    //       return false;
    //     } 
  } 
  ntries += i;
  if( i >= 1000) {
    //std::cout << "failed to find good root\n";
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << " truth:" << truth << "fs:" << f1 << "," << f2 << "\n";
    ++nfails;
    return false;
  }
  
  if(x2 != x2) {
    std::cout << "failed to find good root\n";
    std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
    ++nfails;
    return false;
  }
  
  return true;
}

bool Jet::GslImplementation::root(double truth, double& x1, double& x2, double eps) {
  par.y = truth;
  ++ncalls;
 
  if(gsl_root_fsolver_set(s,&F,x1,x2)) {
    //std::cout << "Warning: root not bracketed\n";
    ++nfails;
    return false;
  }
  int status, iter = 0;
  double r;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x1 = gsl_root_fsolver_x_lower(s);
    x2 = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x1,x2,0, eps);
    }
  while(status == GSL_CONTINUE && iter < 100);
  ntries += iter;
  if(status != GSL_SUCCESS) {
    std::cout << "inversion failed:" << x1 << ", " << x2 << " truth:" 
	      << truth << std::endl;
    ++nfails;
    return false;
  }
  x2 = r;
  return true;
}

// ------------------------------------------------------------
void Jet::print()
{
  std::cout << "Jet  Et: " << Et() << " GeV, eta: " << eta() << std::endl;
  std::cout << "Jet par: ";
  for(int i = 0; i < f.nPars(); i++) {
    std::cout << f.firstPar()[i] << "  ";
  }
  std::cout << "\nGlobal par: ";
  for(int i = 0; i < gf.nPars(); i++) {
    std::cout << gf.firstPar()[i] << "  ";
  }
  std::cout << "\n";
}

Jet::GslImplementation::GslImplementation(const Jet* jet) 
  : par(0,jet),s(gsl_root_fsolver_alloc(gsl_root_fsolver_brent))
{
  F.function = &rf;
  F.params = &par;
  gsl_set_error_handler_off();
}

Jet::GslImplementation::~GslImplementation() {
  gsl_root_fsolver_free(s);
} 

// ------------------------------------------------------------
long long Jet::ncalls = 0;
long long Jet::ntries = 0;
long long Jet::nfails = 0;
long long Jet::nwarns = 0;
void Jet::printInversionStats()
{
  if(ncalls) {
    std::cout << "Inversion statistics for Jet::expectedEt:\n";
    std::cout << "calls: " << ncalls << " average number of iterations:"
	      << (double)ntries/ncalls << " failures:" << (double)nfails/ncalls*100
	      << "% warnings:" << (double)nwarns/ntries*100 << "%" <<std::endl;
  }
}
