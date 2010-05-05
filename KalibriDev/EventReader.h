//
// $Id: EventReader.h,v 1.11 2010/04/13 13:44:09 mschrode Exp $
//
#ifndef EVENTREADER_H
#define EVENTREADER_H

class Event;
class TParameters;
class ConfigFile;
class Measurement;
class CorFactors;
class CorFactorsFactory;
class TTree;
class JetConstraintEvent;
class Parametrization;

#include <vector>
#include <string>


class EventReader
{
 public:
  static unsigned int numberOfEventReaders_;   //!< Number of initialized event readers

  EventReader(const std::string& configfile, TParameters* p);
  virtual ~EventReader();
  virtual int readEvents(std::vector<Event*>& data) = 0;

  static int addConstraints(std::vector<Event*>& data);

 protected:
  //! Read CorFactors from Ntuple
  virtual CorFactors* createCorFactors(int jetid) const { return 0;}
  //! Create TTree with data files
  TTree * createTree(const std::string &dataType) const;

  ConfigFile* config_;   //!< The configfile
  TParameters* par_;     //!< The parametrization
  bool useTracks_;       //!< True, if tracks are used in calibration
  CorFactorsFactory* corFactorsFactory_; //! Factory class for external source of CorFactors;
  //! Correct jets to L3 i.e. with L1*L2*L3
  bool correctToL3_;
  //!< Correct jets with L2*L3 corrections
  bool correctL2L3_;


  double (*tower_error_param)(const double *x, const Measurement *xorig, double err);
  double (*jet_error_param)  (const double *x, const Measurement *xorig, double err);
  double (*track_error_param)(const double *x, const Measurement *xorig, double err);  

  static std::vector<JetConstraintEvent*> constraints_;
  Parametrization *cp_;
};


#endif
