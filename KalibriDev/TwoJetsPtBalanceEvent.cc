// $Id: TwoJetsPtBalanceEvent.cc,v 1.7 2009/11/24 16:52:59 stadie Exp $

#include "TwoJetsPtBalanceEvent.h"

#include <iomanip>
#include <algorithm>

//!  \brief Calculates \f$ \chi^{2} \f$ from pt difference
// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_simple(double * temp_derivative1, 
					       double * temp_derivative2, double const epsilon) const
{
  // Corrected jet pt
  double pt1 = GetParametrizedMess();
  double pt2 = GetParametrizedMess2();

  // Residual
  double res = chi2_fast_simple_res(pt1,pt2);

  // Squared error on residual
  double dRes2 = chi2_fast_simple_dRes2(pt1, pt2);

  // Likelihood
  double chi2 = res * res / dRes2;
  if(chi2 != chi2) {//check for NAN
    std::cout <<pt1 << ", " << pt2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = GetWeight() * Event::ScaleResidual(chi2);

  if(!temp_derivative1) return chi2;

  // Derivative calculation
  const Jet::VariationColl& varColl1 = jet1_->varyParsDirectly(epsilon);
  const Jet::VariationColl& varColl2 = jet2_->varyParsDirectly(epsilon);

  // Variation of parameters of first jet
  for(Jet::VariationCollIter i1 = varColl1.begin() ; i1 != varColl1.end() ; ++i1) {
    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = i1->lowerEt;
    double pt2tmp = 0.;
    Jet::VariationCollIter i2 = find(varColl2.begin(),varColl2.end(),i1->parid);
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->lowerEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = GetWeight() * Event::ScaleResidual(temp1);


    // Corrected pt and derivative in case of 
    // upper parameter variation
    pt1tmp = i1->upperEt;
    pt2tmp = 0.;
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->upperEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = GetWeight() * Event::ScaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of first jet

  // Variation of parameters of second jet
  for(Jet::VariationCollIter i2 = varColl2.begin() ; i2 != varColl2.end() ; ++i2) {
    Jet::VariationCollIter i1 = find(varColl1.begin(),varColl1.end(),i2->parid);
    if(i1 != varColl1.end()) continue; // Parameter already covered by jet 1 par variation

    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = pt1;
    double pt2tmp = i2->lowerEt;

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = GetWeight() * Event::ScaleResidual(temp1);


    // Corrected pt in case of 
    // upper parameter variation
    pt2tmp = i2->upperEt;

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = GetWeight() * Event::ScaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of second jet

  return chi2;
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_simple_res(double pt1, double pt2) const {
  return pt1 - pt2;
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_simple_dRes2(double pt1, double pt2) const {
  // Derivativ d(corr pt) / d(pt); for simplification
  // assuming constant correction function
  double dPt1 = pt1 / jet1_->Et();
  double dPt2 = pt2 / jet2_->Et();

  double dRes2 = dPt1 * dPt1 * error1_ * error1_;
  dRes2       += dPt2 * dPt2 * error2_ * error2_;
  
  return dRes2;
}



//!  \brief Calculates \f$ \chi^{2} \f$ from pt balance
// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_balance(double * temp_derivative1, 
						double * temp_derivative2, double const epsilon) const
{
  // Corrected jet pt
  double pt1 = GetParametrizedMess();
  double pt2 = GetParametrizedMess2();

  // Residual
  double res = chi2_fast_balance_res(pt1,pt2);

  // Squared error on residual
  double dRes2 = chi2_fast_balance_dRes2(pt1, pt2);

  // Likelihood
  double chi2 = res * res / dRes2;
  if(chi2 != chi2) {//check for NAN
    std::cout <<pt1 << ", " << pt2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = GetWeight() * Event::ScaleResidual(chi2);

  if(!temp_derivative1) return chi2;

  // Derivative calculation
  const Jet::VariationColl& varColl1 = jet1_->varyParsDirectly(epsilon);
  const Jet::VariationColl& varColl2 = jet2_->varyParsDirectly(epsilon);

  // Variation of parameters of first jet
  for(Jet::VariationCollIter i1 = varColl1.begin() ; i1 != varColl1.end() ; ++i1) {
    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = i1->lowerEt;
    double pt2tmp = 0.;
    Jet::VariationCollIter i2 = find(varColl2.begin(),varColl2.end(),i1->parid);
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->lowerEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = GetWeight() * Event::ScaleResidual(temp1);


    // Corrected pt and derivative in case of 
    // upper parameter variation
    pt1tmp = i1->upperEt;
    pt2tmp = 0.;
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->upperEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = GetWeight() * Event::ScaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of first jet

  // Variation of parameters of second jet
  for(Jet::VariationCollIter i2 = varColl2.begin() ; i2 != varColl2.end() ; ++i2) {
    Jet::VariationCollIter i1 = find(varColl1.begin(),varColl1.end(),i2->parid);
    if(i1 != varColl1.end()) continue; // Parameter already covered by jet 1 par variation

    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = pt1;
    double pt2tmp = i2->lowerEt;

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = GetWeight() * Event::ScaleResidual(temp1);


    // Corrected pt in case of 
    // upper parameter variation
    pt2tmp = i2->upperEt;

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = GetWeight() * Event::ScaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of second jet

  return chi2;
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_balance_res(double pt1, double pt2) const {
  return 2. * (pt1 - pt2) / (pt1 + pt2);
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_balance_dRes2(double pt1, double pt2) const {
  double ptDijet = (pt1 + pt2) / 2.;
  double ptBal = (pt1 - pt2) / (pt1 + pt2);

  // Derivativ d(corr pt) / d(pt); for simplification
  // assuming constant correction function
  double dPt1 = pt1 / jet1_->Et();
  double dPt2 = pt2 / jet2_->Et();

  double dRes1 = (1. - ptBal) * dPt1 * error1_;
  double dRes2 = (1. + ptBal) * dPt2 * error2_;
  
  return ( dRes1*dRes1 + dRes2*dRes2 ) / ptDijet / ptDijet;
}



//!  \brief Return absolute uncorrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbs(double pt1, double pt2) const {
  double x = pt1 * cos(jet1_->phi()) + pt2 * cos(jet2_->phi());
  double y = pt1 * sin(jet1_->phi()) + pt2 * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute uncorrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbs() const {
  double x = jet1_->Et() * cos(jet1_->phi()) + jet2_->Et() * cos(jet2_->phi());
  double y = jet1_->Et() * sin(jet1_->phi()) + jet2_->Et() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute generated pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbsGen() const {
  double x = jet1_->genPt() * cos(jet1_->phi()) + jet2_->genPt() * cos(jet2_->phi());
  double y = jet1_->genPt() * sin(jet1_->phi()) + jet2_->genPt() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute corrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbsCorr() const {
  double x = jet1_->correctedEt() * cos(jet1_->phi()) + jet2_->correctedEt() * cos(jet2_->phi());
  double y = jet1_->correctedEt() * sin(jet1_->phi()) + jet2_->correctedEt() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute L2L3 corrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbsCorrL2L3() const {
  double x = jet1_->corFactors().getL2L3() * jet1_->Et() * cos(jet1_->phi())
    + jet2_->corFactors().getL2L3() * jet2_->Et() * cos(jet2_->phi());
  double y = jet1_->corFactors().getL2L3() * jet1_->Et() * sin(jet1_->phi())
    + jet2_->corFactors().getL2L3() * jet2_->Et() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}
