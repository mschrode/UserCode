//
//    Reader for Tower Constraints
//
//    This class add user defined tower constraints
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TowerConstraintsReader.h,v 1.1 2008/12/12 17:06:00 stadie Exp $
//   
#ifndef TOWERCONSTRAINTSREADER_H
#define TOWERCONSTRAINTSREADER_H

#include "EventReader.h"

#include <string>

class TowerConstraintsReader : public EventReader{
 public:
  TowerConstraintsReader(const std::string& configfile, TParameters *p);
  virtual ~TowerConstraintsReader();
  int readEvents(std::vector<TData*>& data);
 private:
  class TowerConstraint {
  public:
    int mineta;
    int maxeta;
    double hadEt;
    double emEt;
    double weight;
    TowerConstraint(int mineta,int maxeta,double hadEt, double emEt, double weight) :
      mineta(mineta),maxeta(maxeta),hadEt(hadEt),emEt(emEt),weight(weight) {}
  };
  std::vector<TowerConstraint> tower_constraints;
};


#endif
