// $Id: ResolutionFunction.h,v 1.4 2011/03/01 16:52:41 mschrode Exp $

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

    static ResolutionFunction* createScaledResolutionFunction(const ResolutionFunction* rf, double scale, double relUncertDown = 0., double relUncertUp = 0.);
    static ResolutionFunction* createResolutionFunction(Type type, const TF1* f, double relUncertDown = 0., double relUncertUp = 0.);
    static ResolutionFunction* createResolutionFunction(Type type, std::vector<double> &param, double relUncertDown = 0., double relUncertUp = 0.);
    static ResolutionFunction* fitTGraph(const TGraphAsymmErrors* g, ResolutionFunction::Type type);


    virtual ~ResolutionFunction();

    virtual Type type() const = 0;

    TF1* func(const TString &name) const { return static_cast<TF1*>(func_->Clone(name)); }
    unsigned int nPars() const { return func_->GetNpar(); }
    double par(unsigned int parIdx) const { return func_->GetParameter(parIdx); }
    double val(double pt) const { return func_->Eval(pt); }
    double relUncertDown(double pt) const { return relUncertDown_; }
    double relUncertUp(double pt) const { return relUncertUp_; }
    double uncertDown(double pt) const { return val(pt)*relUncertDown(pt); }
    double uncertUp(double pt) const { return val(pt)*relUncertUp(pt); }
    double ptMin() const { return min_; }
    double ptMax() const { return max_; }


  protected:
    TF1* func_;

    ResolutionFunction(double ptMin, double ptMax, double relUncertDown, double relUncertUp);


  private:
    const double min_;
    const double max_;
    const double relUncertDown_;
    const double relUncertUp_;
  };



  // -------------------------------------------------------------------------------------
  class ResolutionFunctionNSC : public ResolutionFunction {
  public:
    ResolutionFunctionNSC(double ptMin, double ptMax, double n, double s, double c, double relUncertDown, double relUncertUp);

    Type type() const { return NSC; }
  };



  // -------------------------------------------------------------------------------------
  class ResolutionFunctionModifiedNSC : public ResolutionFunction {
  public:
    ResolutionFunctionModifiedNSC(double ptMin, double ptMax, double n, double s, double c, double m, double relUncertDown, double relUncertUp);

    Type type() const { return ModifiedNSC; }
  };


  // -------------------------------------------------------------------------------------
  class ResolutionFunctionScaledModifiedNSC : public ResolutionFunction {
  public:
    ResolutionFunctionScaledModifiedNSC(double ptMin, double ptMax, double n, double s, double c, double m, double scale, double relUncertDown, double relUncertUp);

    Type type() const { return ScaledModifiedNSC; }
  };

}
#endif
