// $Id: ResolutionFunction.h,v 1.5 2011/11/21 17:18:05 mschrode Exp $

#ifndef RESOLUTION_FUNCTION_H
#define RESOLUTION_FUNCTION_H

#include <vector>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class ResolutionFunction {
  public:
    enum Type { NSC, ModifiedNSC, ScaledModifiedNSC };

    static ResolutionFunction* createScaledResolutionFunction(const ResolutionFunction* rf, double scale, double scaleErr, double relUncertDown = 0., double relUncertUp = 0.);
    static ResolutionFunction* createResolutionFunction(Type type, const TF1* f, double relUncertDown = 0., double relUncertUp = 0.);
    static ResolutionFunction* createResolutionFunction(Type type, std::vector<double> &param, std::vector<double> &paramErr, double relUncertDown = 0., double relUncertUp = 0.);
    static ResolutionFunction* fitTGraph(const TGraphAsymmErrors* g, ResolutionFunction::Type type);


    virtual ~ResolutionFunction();

    virtual Type type() const = 0;

    TF1* func(const TString &name) const { return static_cast<TF1*>(func_->Clone(name)); }
    unsigned int nPars() const { return func_->GetNpar(); }
    double par(unsigned int parIdx) const { return func_->GetParameter(parIdx); }
    double parErr(unsigned int parIdx) const { return func_->GetParError(parIdx); }
    double val(double pt) const { return func_->Eval(pt); }
    double relUncertDown(double pt) const { return relUncertDown_; }
    double relUncertUp(double pt) const { return relUncertUp_; }
    double uncertDown(double pt) const { return val(pt)*relUncertDown(pt); }
    double uncertUp(double pt) const { return val(pt)*relUncertUp(pt); }
    double ptMin() const { return min_; }
    double ptMinErr() const { return minErr_; }
    double ptMax() const { return max_; }
    double ptMaxErr() const { return maxErr_; }


  protected:
    TF1* func_;

    ResolutionFunction(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double relUncertDown, double relUncertUp);


  private:
    const double min_;
    const double minErr_;
    const double max_;
    const double maxErr_;
    const double relUncertDown_;
    const double relUncertUp_;
  };



  // -------------------------------------------------------------------------------------
  class ResolutionFunctionNSC : public ResolutionFunction {
  public:
    ResolutionFunctionNSC();
    ResolutionFunctionNSC(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double n, double ne, double s, double se, double c, double ce, double relUncertDown, double relUncertUp);

    Type type() const { return NSC; }
  };



  // -------------------------------------------------------------------------------------
  class ResolutionFunctionModifiedNSC : public ResolutionFunction {
  public:
    ResolutionFunctionModifiedNSC();
    ResolutionFunctionModifiedNSC(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double n, double ne, double s, double se, double c, double ce, double m, double me, double relUncertDown, double relUncertUp);

    Type type() const { return ModifiedNSC; }
  };


  // -------------------------------------------------------------------------------------
  class ResolutionFunctionScaledModifiedNSC : public ResolutionFunction {
  public:
    ResolutionFunctionScaledModifiedNSC();
    ResolutionFunctionScaledModifiedNSC(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double n, double ne, double s, double se, double c, double ce, double m, double me, double scale, double scaleErr, double relUncertDown, double relUncertUp);

    Type type() const { return ScaledModifiedNSC; }
  };

}
#endif
