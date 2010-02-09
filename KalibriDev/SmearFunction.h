#ifndef SMEAR_FUNCTION_H
#define SMEAR_FUNCTION_H

#include <cassert>
 
class Parametrization;


class SmearFunction {
  typedef double (Parametrization::*RespPDF)(double r, double pt, const double*) const; 
  typedef double (Parametrization::*RespError)(double r, double pt, const double*, const double*, const std::vector<int>&) const; 
  typedef double (Parametrization::*TruthPDF)(double pt, const double*, const double*) const; 
  typedef double (Parametrization::*TruthError)(double pt, const double*, const double*, const double*, const std::vector<int>&) const; 

 public:
  SmearFunction(RespPDF respPDF, RespError respErr, int respParIdx, int nRespPars, double *firstRespPar,
		const std::vector<int>& respCovIdx, const std::vector<bool>& isFixedRespPar,
		TruthPDF truthPDF, TruthError truthErr, int truthParIdx, int nTruthPars, double *firstTruthPar,
		const std::vector<int>& truthCovIdx, const std::vector<bool>& isFixedTruthPar,
		double *firstCov, const Parametrization *p)
    : respPDF_(respPDF), respError_(respErr), respParIdx_(respParIdx),
    nRespPars_(nRespPars), firstRespPar_(firstRespPar),
    respCovIdx_(respCovIdx), isFixedRespPar_(isFixedRespPar),
    truthPDF_(truthPDF), truthError_(truthErr), truthParIdx_(truthParIdx),
    nTruthPars_(nTruthPars), firstTruthPar_(firstTruthPar),
    truthCovIdx_(truthCovIdx), isFixedTruthPar_(isFixedTruthPar),
    cov_(firstCov), param_(p) {}

  int nRespPars() const { return nRespPars_;}
  int nTruthPars() const { return nTruthPars_;}
  int respParIdx() const { return respParIdx_; }
  int truthParIdx() const { return truthParIdx_; }
  double respPar(int i) const { assert( i >= 0 && i < nRespPars() ); return firstRespPar_[i]; }
  double truthPar(int i) const { assert( i >= 0 && i < nTruthPars() ); return firstTruthPar_[i]; }
  bool isFixedRespPar(int i) const { assert( i >= 0 && i < nRespPars() ); return isFixedRespPar_[i]; }
  bool isFixedTruthPar(int i) const { assert( i >= 0 && i < nTruthPars() ); return isFixedTruthPar_[i]; }
  double respPDF(double r, double pt) const { return (param_->*respPDF_)(r,pt,firstRespPar_); }
  double respError(double r, double pt) const {
    return (param_->*respError_)(r,pt,firstRespPar_,cov_,respCovIdx_);
  }
  double truthPDF(double pt) const { return (param_->*truthPDF_)(pt,firstTruthPar_,firstRespPar_); }
  double truthError(double pt) const {
    return (param_->*truthError_)(pt,firstTruthPar_,firstRespPar_,cov_,truthCovIdx_);
  }

  void changeParBase(double* oldpar, double* newpar) {
    firstRespPar_ += newpar - oldpar;
    firstTruthPar_ += newpar - oldpar;
  }

  double *respPar() const { return firstRespPar_; }
  double *truthPar() const { return firstTruthPar_; }


 private:
  const RespPDF respPDF_;
  const RespError respError_;
  int respParIdx_;
  int nRespPars_;
  double *firstRespPar_;
  std::vector<int> respCovIdx_;
  std::vector<bool> isFixedRespPar_;

  const TruthPDF truthPDF_;
  const TruthError truthError_;
  int truthParIdx_;
  int nTruthPars_;
  double *firstTruthPar_;
  std::vector<int> truthCovIdx_;
  std::vector<bool> isFixedTruthPar_;

  double *cov_;
  bool *parStatus_;

  const Parametrization* param_;
};
#endif
