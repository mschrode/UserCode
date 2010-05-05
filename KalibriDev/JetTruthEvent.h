//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.h,v 1.14 2010/04/13 13:44:10 mschrode Exp $
//   
#ifndef JETTRUTHEVENT_H
#define JETTRUTHEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class JetTruthEvent : public Event
{
public:
  JetTruthEvent(Jet *j, double t, double w) : jet_(j),truth_(t),weight_(w),chi2plots_(1000.),flagged_bad_(false) {}
  ~JetTruthEvent();

  //interface from Event
  Measurement *GetMess() const {return jet_;}
  double GetTruth() const { return truth_;}
  double GetParametrizedMess() const { return jet_->correctedEt(jet_->Et());}

  void ChangeParAddress(double* oldpar, double* newpar) { jet_->ChangeParAddress(oldpar,newpar);}
  DataType GetType() const { return GammaJet;} 
  double GetWeight() const { return weight_;}
  void setWeight(double w) { weight_ = w; }
  Jet* jet() const {return jet_;}
  
  double chi2() const;
  double chi2_plots() const { return chi2plots_; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const { 
    //chi2plots = chi2_fast_invert(temp_derivative1,temp_derivative2,epsilon);
    chi2plots_ = chi2_log_fast_invert(temp_derivative1,temp_derivative2,epsilon);
    //chi2plots = chi2_fast_scaled(temp_derivative1,temp_derivative2,epsilon);
    //chi2plots = chi2_fast_simple_scaled(temp_derivative1,temp_derivative2,epsilon);
    return chi2plots_;
  }
  double chi2_fast_blobel(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_scaled(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_simple_scaled(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_invert(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_log_fast_invert(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  void updateError() { }
  bool FlaggedBad() const { return flagged_bad_; }  //!< Status from inversion procedure

  static void printStats();
 private:
  Jet* jet_;
  double truth_;
  double weight_;
  mutable double chi2plots_;   //!< Store chi2 value from last iteration for plots
  mutable bool flagged_bad_;
  static int nflagged_;
};

#endif
