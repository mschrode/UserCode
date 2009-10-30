#ifndef FUNCTION_H
#define FUNCTION_H
 
class TMeasurement;
class Parametrization;


class Function {
  
  typedef double (Parametrization::*ParametrizationFunction)(const TMeasurement*, const double*) const; 

 public:
  Function(ParametrizationFunction func, ParametrizationFunction invfunc,
	   double *firstpar, int parindex, int npars, const Parametrization* p) 
    : func(func),invfunc(invfunc),firstpar(firstpar),parindex(parindex),npars(npars),param(p)
    {}
  double* firstPar() const { return firstpar;}
  int parIndex() const { return parindex;}
  int nPars() const { return npars;}
  double operator()(const TMeasurement* x) const { return (param->*func)(x,firstpar);}
  void changeParBase(double* oldpar, double* newpar) { firstpar += newpar - oldpar;}
  
  bool hasInverse() { return invfunc;}
  double inverse(const TMeasurement* x) const { return invfunc ? (param->*invfunc)(x,firstpar) : 0;}
 private:
  const ParametrizationFunction func;
  const ParametrizationFunction invfunc;
  double *firstpar;
  int parindex;
  int npars;
  const Parametrization* param;
};
#endif
