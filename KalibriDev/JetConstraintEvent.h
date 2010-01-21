//
//    Class for constraints on the jet correction
//
//    first version: Hartmut Stadie 2009/07/23
//    $Id: JetConstraintEvent.h,v 1.4 2009/11/24 16:52:59 stadie Exp $
//   
#ifndef JETCONSTRAINTEVENT_H
#define JETCONSTRAINTEVENT_H

#include <map>
#include <vector>

#include"CalibData.h"
#include "Jet.h"




//interface to Data
class JetConstraintEvent : public Event
{
 public:
  JetConstraintEvent(double t, double w) : truth(t), weight(w) {}
  ~JetConstraintEvent();
    
  void addJet(Jet* j);
  

  //interface from TData
  Measurement *GetMess() const {return jets[0];}
  double GetTruth() const { return truth;}
  double GetParametrizedMess() const { return jets[0]->correctedEt(jets[0]->Et());}

  void ChangeParAddress(double* oldpar, double* newpar);
  DataType GetType() const { return JetConstraint;} 
  double GetWeight() const { return weight;}
  
  double chi2() const;
  double chi2_plots() const { return chi2plots; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  void updateError() { } 
  void setWeight(double w) {weight = w;}
  double ptHat() const { return truth;}
 private:
  struct Variation{
    double uppersum;
    double lowersum;
    Variation() : uppersum(0), lowersum(0) {}
  };
  typedef std::map<int, Variation> VarMap;
  mutable VarMap varmap;
  std::vector<Jet*> jets;
  double truth;
  double weight;
  mutable double chi2plots;   //!< Store chi2 value from last iteration for plots
};

#endif
