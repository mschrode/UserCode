// $Id: ResolutionFunction.cc,v 1.1 2011/02/15 18:22:25 mschrode Exp $

#include "ResolutionFunction.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createResolutionFunction(Type type, std::vector<double> &param) {
    ResolutionFunction* func = 0;
    if( type == NSC ) {
      if( param.size() >= 5 ) {
	func = new ResolutionFunctionNSC(param.at(0),param.at(1),param.at(2),param.at(3),param.at(4));
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'NSC'" << std::endl;
      exit(1);
      }
    } else if( type == ModifiedNSC ) {
      if( param.size() >= 6 ) {
	func = new ResolutionFunctionModifiedNSC(param.at(0),param.at(1),param.at(2),param.at(3),param.at(4),param.at(5));
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'ModifiedNSC'" << std::endl;
      exit(1);
      }
    } else {
      std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }

    return func;
  }



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
