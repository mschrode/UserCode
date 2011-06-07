// $Id: FitResult.h,v 1.7 2011/03/02 16:06:01 mschrode Exp $

#ifndef FIT_RESULT_H
#define FIT_RESULT_H

#include <map>
#include <set>
#include <vector>

#include "TF1.h"
#include "TString.h"

#include "Extrapolation.h"
#include "Measurement.h"

class TH1;

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class FitResult {
  public:
    enum Type { MaxLikeKSoftRel, FullMaxLikeRel, FullMaxLikeAbs, SimpleMaxLike, PtAsym, PtGenAsym };

    static bool validType(Type type);
    static FitResult* createFitResult(Type type, const std::vector<Measurement*> &meas, unsigned int verbosity);
    static TString toString(FitResult::Type type);

    FitResult(const std::vector<Measurement*> &meas, unsigned int verbosity);
    virtual ~FitResult();

    virtual FitResult::Type fitResultType() const = 0;

    double minPt() const { return meas_.front()->ptMin(); }
    double maxPt() const { return meas_.front()->ptMax(); }
    double meanPt() const { return meanPt_; }
    double meanPtUncert() const { return meanPtUncert_; }
    double ptSoft(unsigned int ptSoftBin) const { return ptSoft_.at(ptSoftBin); }
    virtual TString labelQuantityInExtrapolation() const {
      return "#sigma / p^{ref}_{T}";
    }
    virtual double valueInExtrapolation(unsigned int ptSoftBin) const { 
      return values_.at(ptSoftBin);
    }
    virtual double statUncertInExtrapolation(unsigned int ptSoftBin) const {
      return values_.at(ptSoftBin);
    }
    double value(unsigned int ptSoftBin) const { return values_.at(ptSoftBin); }
    double statUncert(unsigned int ptSoftBin) const { return statUncerts_.at(ptSoftBin); }
    TF1* extrapolationFunction(const TString &name) const {
      return static_cast<TF1*>(extrapolation_->Clone(name)); }
    virtual double extrapolatedValue() const { return extrapolatedValue_; }
    double extrapolatedStatUncert() const { return extrapolatedStatUncert_; }
    double extrapolatedSystUncert() const { return extrapolatedSystUncert_; }
    virtual TH1* spectrum() const;
    double kSoftSlope() const;
    double kSoftSlopeStatUncert() const;
    TF1* kSoftFit(const TString &name) const { return static_cast<TF1*>(kSoftFit_->Clone(name)); }

    virtual void setKSoftFit(const TF1* fit);

    virtual bool init() = 0;


    
  protected:
    static unsigned int HIST_COUNT;

    const std::vector<Measurement*> meas_;
    const unsigned int verbosity_;
    const unsigned int workingPointBin_;

    double meanPt_;
    double meanPtUncert_;
    std::vector<double> values_;
    std::vector<double> statUncerts_;
    std::vector<double> ptSoft_;
    double extrapolatedValue_;
    double extrapolatedStatUncert_;
    double extrapolatedSystUncert_;
    TF1* extrapolation_;
    TF1* kSoftFit_;

    bool extrapolate();
  };

  typedef std::set<FitResult::Type> FitResultTypes;
  typedef std::set<FitResult::Type>::const_iterator FitResultTypeIt;



  // -------------------------------------------------------------------------------------
  class FitResultMaxLikeKSoftRel : public FitResult {
  public:
    FitResultMaxLikeKSoftRel(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::MaxLikeKSoftRel; }

    virtual double extrapolatedValue() const { return value(workingPointBin_)*kSoftFit_->Eval(meanPt()); }
    virtual void setKSoftFit(const TF1* fit);

    bool init();
  };


  

  // -------------------------------------------------------------------------------------
  class FitResultFullMaxLikeRel : public FitResult {
  public:
    FitResultFullMaxLikeRel(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::FullMaxLikeRel; }

    bool init();
  };



  // -------------------------------------------------------------------------------------
  class FitResultFullMaxLikeAbs : public FitResult {
  public:
    FitResultFullMaxLikeAbs(const std::vector<Measurement*> meas, unsigned int verbosity);
    ~FitResultFullMaxLikeAbs();

    FitResult::Type fitResultType() const { return FitResult::FullMaxLikeAbs; }

    bool init();

    virtual TString labelQuantityInExtrapolation() const {
      return "#sigma (GeV)";
    }
    double valueInExtrapolation(unsigned int ptSoftBin) const { 
      return meanPt()*values_.at(ptSoftBin);
    }
    double statUncertInExtrapolation(unsigned int ptSoftBin) const {
      return meanPt()*statUncerts_.at(ptSoftBin);
    }
    TH1* spectrum() const;

    
  private:
    TH1* spectrum_;
  };



  // -------------------------------------------------------------------------------------
  class FitResultPtAsym : public FitResult {
  public:
    FitResultPtAsym(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::PtAsym; }

    bool init();
  };



  // -------------------------------------------------------------------------------------
  class FitResultPtGenAsym : public FitResult {
  public:
    FitResultPtGenAsym(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::PtGenAsym; }

    bool init();
  };
}
#endif
