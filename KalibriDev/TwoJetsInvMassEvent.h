//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//
//    $Id: TwoJetsInvMassEvent.h,v 1.9 2009/11/24 16:52:59 stadie Exp $
//   

#ifndef TWOJETSINVMASSEVENT_H
#define TWOJETSINVMASSEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class TwoJetsInvMassEvent : public Event
{
public:
  TwoJetsInvMassEvent(Jet *j1, Jet *j2, double t, double w, double* p) 
    : jet1(j1), jet2(j2),truth(t),weight(w),flagged_bad(false),chi2plots(1000.),par(p) {}
  ~TwoJetsInvMassEvent() { delete jet1; delete jet2;}

  //interface from Event
  Measurement *GetMess() const {return jet1;}
  double GetTruth() const { return truth;}
  double GetParametrizedMess() const { return jet1->correctedEt(jet1->Et());}
  
  Measurement *GetMess2() const {return jet2;}
  double GetParametrizedMess2() const { return jet2->correctedEt(jet2->Et());}

  Jet *GetJet1() const {return jet1;}
  Jet *GetJet2() const {return jet2;}

  void ChangeParAddress(double* oldpar, double* newpar) {
    par = newpar;
    jet1->ChangeParAddress(oldpar,newpar);
    jet2->ChangeParAddress(oldpar,newpar);
  }
  DataType GetType() const { return InvMass;} 
  double GetWeight() const { return weight;}
  virtual void setWeight(double w) { weight = w; }
  virtual double ptHat() const { return 0.; }             //!< Dummy

  double correctedMass() const;
  
  double chi2() const;
  double chi2_plots() const { return chi2plots; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const { 
    //chi2plots = chi2_fast_simple(temp_derivative1,temp_derivative2,epsilon);
    //chi2plots = chi2_fast_scaled(temp_derivative1,temp_derivative2,epsilon);
    //chi2plots = chi2_fast_const_error(temp_derivative1,temp_derivative2,epsilon);
    chi2plots = chi2_fast_inv(temp_derivative1,temp_derivative2,epsilon);
    return chi2plots;
  }
  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_const_error(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_scaled(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_inv(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  void updateError() {}

 private:
  Jet *jet1,*jet2;
  double truth;
  double weight;
  mutable bool flagged_bad;
  mutable double chi2plots;   //!< Store chi2 value from last iteration for plots
  double *par;
};

#endif
