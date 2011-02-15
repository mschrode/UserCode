// $Id: $

#include "ResolutionFunction.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  ResolutionFunction::~ResolutionFunction() {
    delete func_;
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunctionNSC::ResolutionFunctionNSC(double ptMin, double ptMax, double n, double s, double c) {
    func_ = new TF1("ResolutionFunctionNSC:TF1","sqrt( sq([0]/x) + sq([1])/x + sq([2]) )",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParameter(1,s);
    func_->SetParameter(2,c);
    func_->SetLineWidth(1);
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunctionModifiedNSC::ResolutionFunctionModifiedNSC(double ptMin, double ptMax, double n, double s, double c, double m) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParameter(1,s);
    func_->SetParameter(2,c);
    func_->SetParameter(3,m);      
    func_->SetLineWidth(1);
  }
}
