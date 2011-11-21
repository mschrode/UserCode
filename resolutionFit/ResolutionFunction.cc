// $Id: ResolutionFunction.cc,v 1.4 2011/03/01 16:52:42 mschrode Exp $

#include <algorithm>

#include "ResolutionFunction.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createScaledResolutionFunction(const ResolutionFunction* rf, double scale, double relUncertDown, double relUncertUp) {
    std::vector<double> param;
    param.push_back(rf->ptMin());
    param.push_back(rf->ptMax());
    for(unsigned int i = 0; i < rf->nPars(); ++i) {
      param.push_back(rf->par(i));
    }
    param.push_back(scale);

    Type type;
    if( rf->type() == ModifiedNSC ) type = ScaledModifiedNSC;

    return createResolutionFunction(type,param,relUncertDown,relUncertUp);
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createResolutionFunction(Type type, const TF1* f, double relUncertDown, double relUncertUp) {
    std::vector<double> param;
    param.push_back(f->GetXmin());
    param.push_back(f->GetXmax());
    for(unsigned int i = 0; i < f->GetNpar(); ++i) {
      param.push_back(f->GetParameter(i));
    }

    return createResolutionFunction(type,param,relUncertDown,relUncertUp);
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createResolutionFunction(Type type, std::vector<double> &param, double relUncertDown, double relUncertUp) {
    ResolutionFunction* func = 0;
    if( type == NSC ) {
      if( param.size() >= 5 ) {
	func = new ResolutionFunctionNSC(param.at(0),param.at(1),param.at(2),param.at(3),param.at(4),relUncertDown,relUncertUp);
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'NSC'" << std::endl;
      exit(1);
      }
    } else if( type == ModifiedNSC ) {
      if( param.size() >= 6 ) {
	func = new ResolutionFunctionModifiedNSC(param.at(0),param.at(1),param.at(2),param.at(3),param.at(4),param.at(5),relUncertDown,relUncertUp);
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'ModifiedNSC'" << std::endl;
      exit(1);
      }
    } else if( type == ScaledModifiedNSC ) {
      if( param.size() >= 7 ) {
	func = new ResolutionFunctionScaledModifiedNSC(param.at(0),param.at(1),param.at(2),param.at(3),param.at(4),param.at(5),param.at(6),relUncertDown,relUncertUp);
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'ScaledModifiedNSC'" << std::endl;
      exit(1);
      }
    } else {
      std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }

    return func;
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::fitTGraph(const TGraphAsymmErrors* g, ResolutionFunction::Type type) {

    TGraphAsymmErrors* gc = static_cast<TGraphAsymmErrors*>(g->Clone());

    std::vector<double> param;
    // Min and max pt
    param.push_back(*std::min_element(gc->GetX(),gc->GetX()+gc->GetN()));
    param.push_back(*std::max_element(gc->GetX(),gc->GetX()+gc->GetN()));
    // Choose typical start values
    if( type == ResolutionFunction::NSC ) {
      param.push_back(4.);
      param.push_back(1.);
      param.push_back(0.03);
    } else if( type == ResolutionFunction::ModifiedNSC ) {
      param.push_back(4.);
      param.push_back(0.7);
      param.push_back(0.);
      param.push_back(0.2);
    }
    ResolutionFunction* f = ResolutionFunction::createResolutionFunction(type,param);

    // Get TF1 to fit
    TF1* fit = f->func("tmp");
    if( type == ResolutionFunction::ModifiedNSC ) fit->FixParameter(2,0.);
    gc->Fit(fit,"0QBR");
    std::cout << "Parameters of fitted resolution:\n";
    for(int i = 0; i < fit->GetNpar(); ++i) {
      std::cout << "  " << fit->GetParameter(i) << std::flush;
    }
    std::cout << std::endl;
    ResolutionFunction* result = ResolutionFunction::createResolutionFunction(type,fit);
    delete f;
    delete fit;
    delete gc;
    
    return result;
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunction::ResolutionFunction(double ptMin, double ptMax, double relUncertDown, double relUncertUp)
    : min_(ptMin), max_(ptMax), relUncertDown_(relUncertDown), relUncertUp_(relUncertUp) {}



  // -------------------------------------------------------------------------------------
  ResolutionFunction::~ResolutionFunction() {
    delete func_;
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunctionNSC::ResolutionFunctionNSC(double ptMin, double ptMax, double n, double s, double c, double relUncertDown, double relUncertUp)
 : ResolutionFunction(ptMin,ptMax,relUncertDown,relUncertUp) {
    func_ = new TF1("ResolutionFunctionNSC:TF1","sqrt( sq([0]/x) + sq([1])/x + sq([2]) )",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParameter(1,s);
    func_->SetParameter(2,c);
    func_->SetLineWidth(1);
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunctionModifiedNSC::ResolutionFunctionModifiedNSC(double ptMin, double ptMax, double n, double s, double c, double m, double relUncertDown, double relUncertUp)
 : ResolutionFunction(ptMin,ptMax,relUncertDown,relUncertUp) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParameter(1,s);
    func_->SetParameter(2,c);
    func_->SetParameter(3,m);      
    func_->SetLineWidth(1);
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunctionScaledModifiedNSC::ResolutionFunctionScaledModifiedNSC(double ptMin, double ptMax, double n, double s, double c, double m, double scale, double relUncertDown, double relUncertUp)
    : ResolutionFunction(ptMin,ptMax,relUncertDown,relUncertUp) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","[4]*sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParameter(1,s);
    func_->SetParameter(2,c);
    func_->SetParameter(3,m);      
    func_->SetParameter(4,scale);
    func_->SetLineWidth(1);
  }
}
