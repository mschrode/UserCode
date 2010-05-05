//
//    Class representing a correction function
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Function.h,v 1.6 2010/02/15 12:40:18 stadie Exp $
//   
#ifndef FUNCTION_H
#define FUNCTION_H
 
class Measurement;
class Parametrization;


class Function {
  
  typedef double (Parametrization::*ParametrizationFunction)(const Measurement*, const double*) const; 

 public:
  Function(ParametrizationFunction func, ParametrizationFunction invfunc,
	   double *firstpar, int parindex, int npars, const Parametrization* p)
    : func_(func),invfunc_(invfunc),firstpar_(firstpar),parindex_(parindex),
    npars_(npars),param_(p)
    {}
  double* firstPar() const { return firstpar_;}
  int parIndex() const { return parindex_;}
  int nPars() const { return npars_;}
  double operator()(const Measurement* x) const { return (param_->*func_)(x,firstpar_);}
  void changeParBase(double* oldpar, double* newpar) { firstpar_ += newpar - oldpar;}
  
  bool hasInverse() { return invfunc_;}
  double inverse(const Measurement* x) const { return invfunc_ ? (param_->*invfunc_)(x,firstpar_) : 0;}
 private:
  ParametrizationFunction func_;
  ParametrizationFunction invfunc_;
  double *firstpar_;
  int parindex_;
  int npars_;
  const Parametrization* param_;
};
#endif
