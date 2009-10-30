//
//    Reader for Jet Constraints
//
//    This class add user defined jet constraints
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: JetConstraintsReader.h,v 1.2 2009/10/09 16:40:44 stadie Exp $
//   
#ifndef JETCONSTRAINTSREADER_H
#define JETCONSTRAINTSREADER_H

#include "EventReader.h"

#include <string>
#include <vector>

class Parametrization;

class JetConstraintsReader : public EventReader{
 public:
  JetConstraintsReader(const std::string& configfile, TParameters *p);
  virtual ~JetConstraintsReader();
  int readEvents(std::vector<TData*>& data);
 private:
  class JetConstraint {
  public:
    int mineta;
    int maxeta;
    double Et;
    double weight;
    JetConstraint(int mineta,int maxeta,double Et, double weight) :
      mineta(mineta),maxeta(maxeta),Et(Et),weight(weight) {}
  };
  std::vector<JetConstraint> jet_constraints;
  Parametrization* cp;
};


#endif
