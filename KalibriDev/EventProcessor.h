//
// $Id: EventProcessor.h,v 1.4 2009/11/24 16:52:59 stadie Exp $
//
#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

class Event;
class TParameters;
class TH1;
class TF1;

#include <vector>
#include <string>



// -----------------------------------------------------------------
class EventProcessor
{
 public:
  EventProcessor(const std::string& configfile, TParameters* param);
  virtual ~EventProcessor();
  virtual int process(std::vector<Event*>& data);

 private:
  static void fitWithoutBottom(TH1 * hist, TF1 * func, double bottom=0.33);
  static double gaussStep(double *x, double *par);

  int flattenSpectra(std::vector<Event*>& data);
  void balanceSpectra(std::vector<Event*>& data);
  int getSpectraBin(double m1, double m2, double m3);

  TParameters* par_;
  double etCutOnGamma_;
  double etCutOnJet_;
  bool flattenSpectra_;
  double relWeight_[7];//@@ Replace 7 by something meaningful
  //  std::vector<int> residualScalingScheme_;          // Iteration scheme of scaling of residuals

  typedef std::vector<Event*>::iterator DataIter;
  typedef std::vector<Event*>::const_iterator DataConstIter;
  //"Not-Balanced" Rejection: Make average-fitting equal to peak-fitting
  struct NotBalancedRejection {
  NotBalancedRejection(double *cut, double min, double max):
    _cut(cut),_min(min),_max(max){};
    bool operator()(Event *d);
    double *_cut;
    double _min, _max;
  };
};


#endif
