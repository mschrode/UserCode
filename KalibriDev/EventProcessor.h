//!
//!    \brief Possible base class for all event processors
//!    
//!    Right now we only have the spectrum correctors.
//!    Thus they are implemented directly in this class
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/14
//!    $Id: EventProcessor.h,v 1.3 2009/07/23 13:46:20 mschrode Exp $
// -----------------------------------------------------------------
#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

class TData;
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
  virtual int process(std::vector<TData*>& data);

 private:
  static void fitWithoutBottom(TH1 * hist, TF1 * func, double bottom=0.33);
  static double gaussStep(double *x, double *par);

  int flattenSpectra(std::vector<TData*>& data);
  void balanceSpectra(std::vector<TData*>& data);
  int getSpectraBin(double m1, double m2, double m3);

  TParameters* par_;
  double etCutOnGamma_;
  double etCutOnJet_;
  bool flattenSpectra_;
  double relWeight_[7];//@@ Replace 7 by something meaningful
  //  std::vector<int> residualScalingScheme_;          // Iteration scheme of scaling of residuals

  typedef std::vector<TData*>::iterator DataIter;
  typedef std::vector<TData*>::const_iterator DataConstIter;
  //"Not-Balanced" Rejection: Make average-fitting equal to peak-fitting
  struct NotBalancedRejection {
  NotBalancedRejection(double *cut, double min, double max):
    _cut(cut),_min(min),_max(max){};
    bool operator()(TData *d);
    double *_cut;
    double _min, _max;
  };
};


#endif
