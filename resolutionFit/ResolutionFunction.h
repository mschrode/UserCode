// $Id: $

#ifndef RESOLUTION_FUNCTION_H
#define RESOLUTION_FUNCTION_H

#include "TF1.h"
#include "TString.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class ResolutionFunction {
  public:
    enum Type { NSC, ModifiedNSC };
    
    virtual ~ResolutionFunction();

    virtual Type type() const = 0;

    TF1* func(const TString &name) const { return static_cast<TF1*>(func_->Clone(name)); }
    unsigned int nPars() const { return func_->GetNpar(); }
    double par(unsigned int parIdx) const { return func_->GetParameter(parIdx); }
    double val(double pt) const { return func_->Eval(pt); }


  protected:
    TF1* func_;
  };



  // -------------------------------------------------------------------------------------
  class ResolutionFunctionNSC : public ResolutionFunction {
  public:
    ResolutionFunctionNSC(double ptMin, double ptMax, double n, double s, double c);

    Type type() const { return NSC; }
  };



  // -------------------------------------------------------------------------------------
  class ResolutionFunctionModifiedNSC : public ResolutionFunction {
  public:
    ResolutionFunctionModifiedNSC(double ptMin, double ptMax, double n, double s, double c, double m);

    Type type() const { return ModifiedNSC; }
  };
}
#endif
