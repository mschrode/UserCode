//
//    Reader for Parameter Limits
//
//    This class add user defined parameter limits
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ParameterLimitsReader.h,v 1.1 2008/12/12 17:06:00 stadie Exp $
//   
#ifndef PARAMETERLIMITSREADER_H
#define PARAMETERLIMITSREADER_H

#include "EventReader.h"

#include <string>

class ParameterLimitsReader : public EventReader{
 public:
  ParameterLimitsReader(const std::string& configfile, TParameters *p);
  virtual ~ParameterLimitsReader();
  int readEvents(std::vector<TData*>& data);
 private:
  class ParameterLimit {
  public:
    int index;
    double min;
    double max;
    double k;
    ParameterLimit(int index, double min, double max, double k) 
      : index(index), min(min), max(max), k(k) {}
  };
  std::vector<ParameterLimit> par_limits;
};


#endif
