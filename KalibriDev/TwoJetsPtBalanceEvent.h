#ifndef TWO_JETS_PT_BALANCE_EVENT_H
#define TWO_JETS_PT_BALANCE_EVENT_H

#include <cmath>

#include "CalibData.h"
#include "Jet.h"

//!
//!  \brief Class for relative calibration in pseudorapidity
//!         using dijet events
//!
//!  \author Matthias Schroeder
//!  \date Mon Oct 26 21:03:43 CET 2009 
//!  $Id: TwoJetsPtBalanceEvent.h,v 1.3 2009/11/06 14:14:18 mschrode Exp $
// --------------------------------------------------
class TwoJetsPtBalanceEvent : public TData {
 public:
  TwoJetsPtBalanceEvent(Jet *j1, Jet *j2, Jet *j3, double ptHat, double w) 
    : ptHat_(ptHat),
    jet1_(j1),
    jet2_(j2),
    jet3_(j3),
    weight_(w),
    flaggedBad_(false),
    chi2Plots_(1000.) {
    error1_ = jet1_->Error();
    error2_ = jet2_->Error();
  }
  ~TwoJetsPtBalanceEvent() { delete jet1_; delete jet2_; if( hasJet3() ) delete jet3_; }

  virtual TMeasurement *GetMess() const { return jet1_; }
  virtual double GetParametrizedMess() const { return jet1_->correctedEt();}
  
  TMeasurement *GetMess2() const { return jet2_; }
  double GetParametrizedMess2() const { return jet2_->correctedEt();}

  TMeasurement *GetMess3() const { return jet3_; }
  double GetParametrizedMess3() const { return hasJet3() ? jet3_->correctedEt() : 0.;}

  Jet *getJet1() const { return jet1_; }
  Jet *getJet2() const { return jet2_; }
  Jet *getJet3() const { return jet3_; }

  bool hasJet3() const { return jet3_ != 0 ? true : false; }

  virtual void ChangeParAddress(double* oldpar, double* newpar) {
    jet1_->ChangeParAddress(oldpar,newpar);
    jet2_->ChangeParAddress(oldpar,newpar);
    if( hasJet3() ) jet3_->ChangeParAddress(oldpar,newpar);
  }

  virtual double GetTruth() const { return 0.; }
  virtual DataType GetType() const { return PtBalance; } 
  virtual double GetWeight() const { return weight_; }
  virtual void setWeight(double w) { weight_ = w; }
  virtual double ptHat() const { return ptHat_; }

  double correctedMass() const;  
  virtual double chi2() const { return chi2_fast(0, 0, 0); }
  virtual double chi2_plots() const { return chi2Plots_; }
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const { 
    chi2Plots_ = chi2_fast_balance(temp_derivative1,temp_derivative2,epsilon);
    return chi2Plots_;
  }
  virtual void updateError() {
    error1_ = jet1_->expectedError(ptDijetCorr());
    error2_ = jet2_->expectedError(ptDijetCorr());
  }

  double ptDijet() const { return 0.5*(jet1_->Et()+jet2_->Et()); }
  double ptDijetGen() const { return 0.5*(jet1_->GenPt()+jet2_->GenPt()); }
  double ptDijetCorr() const { return 0.5*(GetParametrizedMess()+GetParametrizedMess2()); }
  double ptDijetCorrL2L3() const { return 0.5*( jet1_->corFactors.getL2L3() * jet1_->Et() + jet2_->corFactors.getL2L3() * jet2_->Et() ); }

  double ptBalance() const { return (jet1_->Et() - jet2_->Et()) / ptDijet(); }
  double ptBalanceGen() const { return (jet1_->GenPt()-jet2_->GenPt()) / ptDijetGen(); }
  double ptBalanceCorr() const { return (GetParametrizedMess()-GetParametrizedMess2()) / ptDijetCorr(); }
  double ptBalanceCorrL2L3() const { return ( jet1_->corFactors.getL2L3() * jet1_->Et() - jet2_->corFactors.getL2L3() * jet2_->Et() ) / ptDijetCorrL2L3(); }

  double ptSumAbs(double pt1, double pt2) const;
  double ptSumAbs() const;
  double ptSumAbsCorr() const;
  double ptSumAbsGen() const;
  double ptSumAbsCorrL2L3() const;

  double relPtJet3() const { return hasJet3() ? jet3_->pt / ptDijet() : 0.; }

  bool flaggedBad() const { return flaggedBad_; }  //!< Status


 private:
  const double ptHat_;

  double error1_;		//!< Store jet1_ error during major iteration
  double error2_;		//!< Store jet2_ error during major iteration
  Jet *jet1_;
  Jet *jet2_;
  Jet *jet3_;
  double weight_;

  mutable bool flaggedBad_;
  mutable double chi2Plots_;   //!< Store chi2 value from last iteration for plots

  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_simple_res(double pt1, double pt2) const;
  double chi2_fast_simple_dRes2(double pt1, double pt2) const;

  double chi2_fast_balance(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_balance_res(double pt1, double pt2) const;
  double chi2_fast_balance_dRes2(double pt1, double pt2) const;
};

#endif
