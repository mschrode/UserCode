//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.cc,v 1.18 2009/11/24 16:52:58 stadie Exp $
//   

#include "JetTruthEvent.h"

int JetTruthEvent::nflagged = 0;

JetTruthEvent::~JetTruthEvent() 
{ 
  delete jet;
}

double JetTruthEvent::chi2() const
{
  double diff = (jet->correctedEt(jet->Et()) - truth)/jet->Error();
  return weight * Event::ScaleResidual(diff*diff);
}

double JetTruthEvent::chi2_fast_blobel(double * temp_derivative1, 
				       double * temp_derivative2, 
				       double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double err2inv = jet->Error();
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * Event::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    //err2inv = jet->expectedError(i->lowerEt);;
    //err2inv *= err2inv;
    //err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * Event::ScaleResidual(temp1);
    temp2 = truth - i->upperEt; 
    //err2inv = jet->expectedError(i->upperEt);
    //err2inv *= err2inv;
    //err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * Event::ScaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double JetTruthEvent::chi2_fast_simple_scaled(double * temp_derivative1, 
					      double * temp_derivative2, 
					      double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double c = et / jet->Et();
  //if(c <= 0) c = 1.0;
  double err2inv = c * jet->expectedError(truth/c);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * Event::ScaleResidual(chi2);
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << et << ", " <<  jet->Et() << ", " << c << ", " << chi2 << '\n';
  }
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    c = i->lowerEt / jet->Et();  
    //if(c <= 0) c = 1.0;
    err2inv = c*jet->expectedError(truth/c );
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * Event::ScaleResidual(temp1);
    assert(temp1 == temp1);
    temp2 = truth - i->upperEt;
    c = i->upperEt / jet->Et();  
    //if(c <= 0) c = 1.0;
    err2inv = c * jet->expectedError(truth/c); 
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * Event::ScaleResidual(temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative   
  }
  assert(chi2 == chi2);
  return chi2;
}


double JetTruthEvent::chi2_fast_scaled(double * temp_derivative1, 
				       double * temp_derivative2, 
				       double const epsilon) const
{
  if(flagged_bad) return 0;
  double et = jet->correctedEt(jet->Et());
  double c = et / jet->Et();
  const double deltaE = epsilon;
  double etprime  = (jet->correctedEt(jet->Et() + deltaE) - 
		     jet->correctedEt(jet->Et() - deltaE))/2/deltaE;
  if((etprime == 0) || (c <= 0)) {
    std::cout << "warning: deriv zero!\n";
    flagged_bad = true;
    ++nflagged;
    return 0;
  }
  double err2 = etprime * jet->expectedError(jet->Et());
  err2 *= err2;
  double chi2 = truth - et;
  chi2 *= chi2 / err2;
  chi2 = weight * Event::ScaleResidual(log(err2) + chi2);
  //chi2 = weight *  Event::ScaleResidual(chi2);
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << et << ", " <<  jet->Et() << ", " << c << ", " << chi2 << '\n';
    assert(chi2 == chi2);
  }
  assert(! std::isinf(chi2));
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    c = i->lowerEt / jet->Et();
    //if(c <= 0) c = 1.0;
    err2 = i->lowerEtDeriv * jet->expectedError(jet->Et());
    if((i->lowerEtDeriv == 0)|| (c <= 0)) {
      std::cout << "warning: deriv zero!\n"; 
      flagged_bad = true;
      ++nflagged;
      return 0;
      //err2 = c * jet->expectedError(jet->Et());
    }
    err2 *= err2;
    temp1 *= temp1 / err2; 
    //temp1 = weight * (log(err2) + temp1);
    temp1 = weight * Event::ScaleResidual(log(err2) + temp1);
    //temp1 = weight * Event::ScaleResidual(temp1);
    assert(temp1 == temp1);
    assert(! std::isinf(temp1));
    temp2 = truth - i->upperEt;
    c = i->upperEt / jet->Et();  
    //if(c <= 0) c = 1.0;
    err2 = i->upperEtDeriv * jet->expectedError(jet->Et());
    if((i->upperEtDeriv == 0) || (c <= 0)) {
      std::cout << "warning: deriv zero!\n";
      flagged_bad = true; 
      ++nflagged;
      return 0;
      //err2 = c * jet->expectedError(jet->Et());
    }
    err2 *= err2;
    temp2 *= temp2 / err2;
    temp2 = weight * Event::ScaleResidual(log(err2) + temp2);
    //temp2 = weight * Event::ScaleResidual(temp2);
    assert(temp2 == temp2);
    assert(! std::isinf(temp2));
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative   
  }
  assert(chi2 == chi2);
  return chi2;
}

double JetTruthEvent::chi2_fast_simple(double * temp_derivative1, 
				       double * temp_derivative2, 
				       double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double c = et/jet->Et(); 
  double err2inv = jet->expectedError(truth/c);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << et << ", " <<  jet->Et() << ", " << c << ", " << chi2 << '\n';
  }
  chi2 = weight * Event::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    c = i->lowerEt/jet->Et(); 
    err2inv = jet->expectedError(truth/c);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * Event::ScaleResidual(temp1);
    temp2 = truth - i->upperEt;
    c = i->upperEt/jet->Et(); 
    err2inv = jet->expectedError(truth/c); 
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * Event::ScaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double JetTruthEvent::chi2_fast_invert(double * temp_derivative1, 
				double * temp_derivative2, 
				double const epsilon) const
{
  if(flagged_bad) return 0;
  //find expected measurement of jet Et 
  double err2inv;
  double expectedEt = jet->expectedEt(truth,truth,err2inv);
  if(expectedEt < 0) {
    flagged_bad = true;
    ++nflagged;
    //std::cout << "Inversion failed: flag event bad.\n";
    return 0;
    //return chi2_fast_simple_scaled(temp_derivative1,temp_derivative2,epsilon);
  }
  //calculate chi2
  //double err2inv = jet->expectedError(expectedEt);
  //std::cout << "Jet Et:" << expectedEt << "," << jet->Et() << " error:" << err2inv 
  //  	    << "  true Et:" << truth << '\n';
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = jet->Et() - expectedEt;
  chi2 *= chi2 * err2inv;
  chi2 = weight * Event::ScaleResidual(chi2);
  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << expectedEt << ", " << chi2 << '\n';
    assert(false);
  }
  
  //calculate chi2 for derivatives
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyPars(epsilon,truth,truth);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if(( i->lowerEt < 0 ) || ( i->upperEt < 0 )) {
      flagged_bad = true;
      ++nflagged;
      return 0;
    }
  }

  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    //std::cout << i->parid << ":" << i->lowerEt << ";" << i->upperEt << " diff:" 
    //	      << i->upperEt - i->lowerEt << ", " << i->upperError <<
    //  ", " << i->lowerError << '\n';
    if((std::abs((i->lowerEt - expectedEt)/expectedEt) > 0.01) || 
       (std::abs((i->upperEt - expectedEt)/expectedEt) > 0.01)) {
      std::cout << "strange extrapolation result modifying par:" << i->parid << ":" 
		<< expectedEt << "  - " << i->lowerEt << "  + " << i->upperEt 
		<< "  uncor  jet Et:" << jet->Et() << " truth:" << truth << std::endl;
      continue;
    }   
    temp1 = i->lowerEt - jet->Et();
    err2inv = i->lowerError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * Event::ScaleResidual(temp1);
    assert(temp1 == temp1);
    temp2 = i->upperEt - jet->Et(); 
    err2inv = i->upperError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * Event::ScaleResidual(temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    assert( temp_derivative1[i->parid] ==  temp_derivative1[i->parid] );
    assert( temp_derivative2[i->parid] ==  temp_derivative2[i->parid] );
  }
  assert(chi2 == chi2);
  if(chi2/weight > 1000) {
    std::cout << "from invert: Jet Et:" << jet->Et() << "  expected Et:" << expectedEt << " error:" << sqrt(1/err2inv)  
	      << "  true Et:" << truth << "  chi2:" << chi2 << '\n';
  }
  return chi2;
}



//!  \brief Use correct likelihood
//!
//!  Calculates the summand of the negative log-likelihood
//!  for non-constant error:
//!  \f[
//!   -2\ln L = c +
//!   \sum_{i}\left(\ln\sigma^{2}_{i} + \left(\frac{x_{i}-\mu_{i}}{\sigma_{i}}\right)^{2}\right)
//!  \f]
//!
//!  \param temp_derivative1 Summand of this event to first derivative
//!  \param temp_derivative2 Summand of this event to second derivative
//!  \param epsilon Step width in numerical derivative calculation
//!  \return Summand of this event to negative log-likelihood
// ------------------------------------------------------------
double JetTruthEvent::chi2_log_fast_invert(double * temp_derivative1, 
					   double * temp_derivative2, 
					   double const epsilon) const
{
  if(flagged_bad) return 0;
  //find expected measurement of jet Et 
  double err2;
  double expectedEt = jet->expectedEt(truth,truth,err2);
  if(expectedEt < 0) {
    flagged_bad = true;
    ++nflagged;
    return 0;
  }
  err2 *= err2;
  double chi2 = jet->Et() - expectedEt;
  chi2 *= chi2 / err2;
  chi2 = weight * Event::ScaleResidual(log(err2) + chi2);
  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << expectedEt << ", " << chi2 << '\n';
  }

  //calculate chi2 for derivatives
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyPars(epsilon,truth,expectedEt);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if(( i->lowerEt < 0 ) || ( i->upperEt < 0 )) {
      flagged_bad = true;
      ++nflagged;
      return 0;
    }
  }

  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if((std::abs((i->lowerEt - expectedEt)/expectedEt) > 0.01) || 
       (std::abs((i->upperEt - expectedEt)/expectedEt) > 0.01)) {
      std::cout << "strange extrapolation result modifying par:" << i->parid << ":" 
		<< expectedEt << "  - " << i->lowerEt << "  + " << i->upperEt 
		<< "  uncor  jet Et:" << jet->Et() << " truth:" << truth << std::endl;
      continue;
    }   
    temp1 = i->lowerEt - jet->Et();
    err2 = i->lowerError;
    err2 *= err2;
    temp1 *= temp1 / err2;
    temp1 = weight * Event::ScaleResidual(log(err2) + temp1);
    assert(temp1 == temp1);
    temp2 = i->upperEt - jet->Et(); 
    err2 = i->upperError;
    err2 *= err2;
    temp2 *= temp2 / err2;
    temp2 = weight * Event::ScaleResidual(log(err2) + temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    assert( temp_derivative1[i->parid] ==  temp_derivative1[i->parid] );
    assert( temp_derivative2[i->parid] ==  temp_derivative2[i->parid] );
  }
  assert(chi2 == chi2);
  if(chi2/weight > 1000) {
    std::cout << "from invert: Jet Et:" << jet->Et() << "  expected Et:" << expectedEt << " error:" << sqrt(err2)  
	      << "  true Et:" << truth << "  chi2:" << chi2 << '\n';
  }
  return chi2;
}
 

void JetTruthEvent::printStats()
{
  std::cout << "JetTruthEvent: " << nflagged << " events flagged bad.\n";
}
