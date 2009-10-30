#ifndef EVENTREADER_H
#define EVENTREADER_H

class TData;
class TParameters;
class ConfigFile;
class TMeasurement;

#include <vector>
#include <string>


//!
//!    \brief Abstract base class for all event readers
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/12
//!    $Id: EventReader.h,v 1.4 2009/10/26 20:56:29 mschrode Exp $
// ----------------------------------------------------------------   
class EventReader
{
 public:
  static unsigned int numberOfEventReaders_;   //!< Number of initialized event readers

  EventReader(const std::string& configfile, TParameters* p);
  virtual ~EventReader();
  virtual int readEvents(std::vector<TData*>& data) = 0;

 protected:
  ConfigFile* config_;   //!< The configfile
  TParameters* par_;     //!< The parametrization
  bool useTracks_;       //!< True, if tracks are used in calibration

  double (*tower_error_param)(const double *x, const TMeasurement *xorig, double err);
  double (*jet_error_param)  (const double *x, const TMeasurement *xorig, double err);
  double (*track_error_param)(const double *x, const TMeasurement *xorig, double err);
};


#endif
