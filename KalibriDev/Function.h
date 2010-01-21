//
//    Class representing a correction function
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Function.h,v 1.4 2009/11/24 16:52:59 stadie Exp $
//   
#ifndef FUNCTION_H
#define FUNCTION_H
 
class Measurement;
class Parametrization;


class Function {
  
  typedef double (Parametrization::*ParametrizationFunction)(const Measurement*, const double*) const; 

 public:
  Function(ParametrizationFunction func, ParametrizationFunction invfunc,
	   int parindex, int npars, double *firstpar, const Parametrization* p) 
    : func_(func),invfunc_(invfunc),firstpar_(firstpar),firstError_(0),
    cov_(0),parindex_(parindex),npars_(npars),param_(p)
    {}
  Function(ParametrizationFunction func, ParametrizationFunction invfunc,
	   int parindex, int npars, double *firstpar, double *firstError, double *cov,
	   const std::vector<int> &covIdx, const Parametrization* p) 
    : func_(func),invfunc_(invfunc),firstpar_(firstpar),firstError_(firstError),
    cov_(cov),covIdx_(covIdx),parindex_(parindex),npars_(npars),param_(p)
    {}
  double* firstPar() const { return firstpar_;}
  double* firstError() const { return firstError_; }
  double* cov() const { return cov_; }
  int parIndex() const { return parindex_;}
  int nPars() const { return npars_;}
  double operator()(const Measurement* x) const { return (param_->*func_)(x,firstpar_);}
  void changeParBase(double* oldpar, double* newpar) { firstpar_ += newpar - oldpar;}

  double par(int i) const { assert( i >= 0 && i < nPars() ); return firstPar()[i]; }
  double parError(int i) const { assert( i >= 0 && i < nPars() ); return firstError()[i]; }
  double parCov(int i, int j) const;
  
  bool hasInverse() { return invfunc_;}
  double inverse(const Measurement* x) const { return invfunc_ ? (param_->*invfunc_)(x,firstpar_) : 0;}

 private:
  const ParametrizationFunction func_;
  const ParametrizationFunction invfunc_;
  double *firstpar_;
  double *firstError_;
  double *cov_;
  std::vector<int> covIdx_;
  int parindex_;
  int npars_;
  const Parametrization* param_;
};
#endif
