// $Id :$

#ifndef RESOLUTION_TAILS_FIT_PARAMETERS
#define RESOLUTION_TAILS_FIT_PARAMETERS


// Collection of fit parameters

namespace resolutionTails {
  class FitParameters {
  public:
    FitParameters() : nSigTailStart_(0.), nSigTailEnd_(0.), nSigCore_(0.), fixDataShape_(false), minPt3Data_(0.) {}
    
    double nSigTailStart() const { return nSigTailStart_; }
    double nSigTailEnd() const { return nSigTailEnd_; }
    double nSigCore() const { return nSigCore_; }
    bool fixDataShape() const { return fixDataShape_; }
    double minPt3Data() const { return minPt3Data_; }

    void setNSigTailStart(double val) { nSigTailStart_ = val; }
    void setNSigTailEnd(double val) { nSigTailEnd_ = val; }
    void setNSigCore(double val) { nSigCore_ = val; }
    void setFixDataShape(bool val) { fixDataShape_ = val; }
    void setMinPt3Data(double val) { minPt3Data_ = val; }


  private:
    double nSigTailStart_;
    double nSigTailEnd_;
    double nSigCore_;
    bool fixDataShape_;
    double minPt3Data_;
  };
}
#endif
