// $Id: ResolutionFunction.cc,v 1.7 2012/05/31 20:17:43 mschrode Exp $

#include <algorithm>

#include "ResolutionFunction.h"

namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createScaledResolutionFunction(const ResolutionFunction* rf, double scale, double scaleErr, double relUncertDown, double relUncertUp) {
    std::vector<double> param;
    std::vector<double> paramErr;
    param.push_back(rf->ptMin());
    paramErr.push_back(rf->ptMinErr());
    param.push_back(rf->ptMax());
    paramErr.push_back(rf->ptMaxErr());
    for(unsigned int i = 0; i < rf->nPars(); ++i) {
      param.push_back(rf->par(i));
      paramErr.push_back(rf->parErr(i));
    }
    param.push_back(scale);
    paramErr.push_back(scaleErr);

    Type type;
    if( rf->type() == ModifiedNSC ) type = ScaledModifiedNSC;

    return createResolutionFunction(type,param,paramErr,relUncertDown,relUncertUp);
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createResolutionFunction(Type type, const TF1* f, double relUncertDown, double relUncertUp) {
    std::vector<double> param;
    std::vector<double> paramErr;
    param.push_back(f->GetXmin());
    paramErr.push_back(0.);
    param.push_back(f->GetXmax());
    paramErr.push_back(0.);
    for(unsigned int i = 0; i < f->GetNpar(); ++i) {
      param.push_back(f->GetParameter(i));
      paramErr.push_back(f->GetParError(i));
    }

    return createResolutionFunction(type,param,paramErr,relUncertDown,relUncertUp);
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunction* ResolutionFunction::createResolutionFunction(Type type, std::vector<double> &param, std::vector<double> &paramErr, double relUncertDown, double relUncertUp) {
    ResolutionFunction* func = 0;
    if( type == NSC ) {
      if( param.size() >= 5 ) {
	func = new ResolutionFunctionNSC(param.at(0),paramErr.at(0),param.at(1),paramErr.at(1),param.at(2),paramErr.at(2),param.at(3),paramErr.at(3),param.at(4),paramErr.at(4),relUncertDown,relUncertUp);
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'NSC'" << std::endl;
      exit(1);
      }
    } else if( type == ModifiedNSC ) {
      if( param.size() >= 6 ) {
	func = new ResolutionFunctionModifiedNSC(param.at(0),paramErr.at(0),param.at(1),paramErr.at(1),param.at(2),paramErr.at(2),param.at(3),paramErr.at(3),param.at(4),paramErr.at(4),param.at(5),paramErr.at(5),relUncertDown,relUncertUp);
      } else {
	std::cerr << "ERROR in ResolutionFunction::createResolutionFunction(): Too few parameters for type 'ModifiedNSC'" << std::endl;
      exit(1);
      }
    } else if( type == ScaledModifiedNSC ) {
      if( param.size() >= 7 ) {
	func = new ResolutionFunctionScaledModifiedNSC(param.at(0),paramErr.at(0),param.at(1),paramErr.at(1),param.at(2),paramErr.at(2),param.at(3),paramErr.at(3),param.at(4),paramErr.at(4),param.at(5),paramErr.at(5),param.at(6),paramErr.at(6),relUncertDown,relUncertUp);
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
  ResolutionFunction* ResolutionFunction::fitTGraph(const TGraphAsymmErrors* g, ResolutionFunction::Type type, bool setPar2to0) {

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
      param.push_back(0.01);
      param.push_back(0.2);
    }
    std::vector<double> paramErr(param.size(),0.);
    ResolutionFunction* f = ResolutionFunction::createResolutionFunction(type,param,paramErr);

    // Get TF1 to fit
    TF1* fit = f->func("tmp");
    if( setPar2to0 && type == ResolutionFunction::ModifiedNSC ) fit->FixParameter(2,0.);
    gc->Fit(fit,"0QR");
    std::cout << "Parameters of fitted resolution:\n";
    for(int i = 0; i < fit->GetNpar(); ++i) {
      std::cout << "  " << i << ": " << fit->GetParameter(i) << " +/- " << fit->GetParError(i) << std::endl;
    }
    std::cout << std::endl;
    ResolutionFunction* result = ResolutionFunction::createResolutionFunction(type,fit);
    delete f;
    delete fit;
    delete gc;
    
    return result;
  }



  // -------------------------------------------------------------------------------------
  ResolutionFunction::ResolutionFunction(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double relUncertDown, double relUncertUp)
    : min_(ptMin), minErr_(ptMinErr), max_(ptMax), maxErr_(ptMaxErr), relUncertDown_(relUncertDown), relUncertUp_(relUncertUp) {}



  // -------------------------------------------------------------------------------------
  ResolutionFunction::~ResolutionFunction() {
    delete func_;
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunctionNSC::ResolutionFunctionNSC()
    : ResolutionFunction(0.,0.,1.,0.,0.,0.) {
    func_ = new TF1("ResolutionFunctionNSC:TF1","sqrt( sq([0]/x) + sq([1])/x + sq([2]) )",0.,1.);
    for(int i = 0; i < func_->GetNpar(); ++i) {
      func_->SetParameter(i,0.);
      func_->SetParError(i,0.);
    }
  }

  // -------------------------------------------------------------------------------------
  ResolutionFunctionNSC::ResolutionFunctionNSC(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double n, double ne, double s, double se, double c, double ce, double relUncertDown, double relUncertUp)
    : ResolutionFunction(ptMin,ptMinErr,ptMax,ptMaxErr,relUncertDown,relUncertUp) {
    func_ = new TF1("ResolutionFunctionNSC:TF1","sqrt( sq([0]/x) + sq([1])/x + sq([2]) )",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParError(0,ne);
    func_->SetParameter(1,s);
    func_->SetParError(1,se);
    func_->SetParameter(2,c);
    func_->SetParError(2,ce);
    func_->SetLineWidth(1);
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunctionModifiedNSC::ResolutionFunctionModifiedNSC()
    : ResolutionFunction(0,1.,0.,0.,0.,0.) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",0.,1.);
    for(int i = 0; i < func_->GetNpar(); ++i) {
      func_->SetParameter(i,0.);
      func_->SetParError(i,0.);
    }
  }

  // -------------------------------------------------------------------------------------
  ResolutionFunctionModifiedNSC::ResolutionFunctionModifiedNSC(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double n, double ne, double s, double se, double c, double ce, double m, double me, double relUncertDown, double relUncertUp)
    : ResolutionFunction(ptMin,ptMinErr,ptMax,ptMaxErr,relUncertDown,relUncertUp) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParError(0,ne);
    func_->SetParameter(1,s);
    func_->SetParError(1,se);
    func_->SetParameter(2,c);
    func_->SetParError(2,ce);
    func_->SetParameter(3,m);      
    func_->SetParError(3,me);
    func_->SetLineWidth(1);
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunctionScaledModifiedNSC::ResolutionFunctionScaledModifiedNSC()
    : ResolutionFunction(0.,1.,0.,0.,0.,0.) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","[4]*sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",0.,1.);
    for(int i = 0; i < func_->GetNpar(); ++i) {
      func_->SetParameter(i,0.);
      func_->SetParError(i,0.);
    }
  }

  // -------------------------------------------------------------------------------------
  ResolutionFunctionScaledModifiedNSC::ResolutionFunctionScaledModifiedNSC(double ptMin, double ptMinErr, double ptMax, double ptMaxErr, double n, double ne, double s, double se, double c, double ce, double m, double me, double scale, double scaleErr, double relUncertDown, double relUncertUp) 
    : ResolutionFunction(ptMin,ptMinErr,ptMax,ptMaxErr,relUncertDown,relUncertUp) {
    func_ = new TF1("ResolutionFunctionModifiedNSC:TF1","[4]*sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",ptMin,ptMax);
    func_->SetParameter(0,n);
    func_->SetParError(0,ne);
    func_->SetParameter(1,s);
    func_->SetParError(1,se);
    func_->SetParameter(2,c);
    func_->SetParError(2,ce);
    func_->SetParameter(3,m);      
    func_->SetParError(3,me);
    func_->SetParameter(4,scale);
    func_->SetParError(4,scaleErr);
    func_->SetLineWidth(1);
  }
}
