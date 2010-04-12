#ifndef SMEAR_FUNCTION_H
#define SMEAR_FUNCTION_H

#include <cassert>
 
class Parametrization;


class SmearFunction {
  typedef double (Parametrization::*PdfPtMeas)(double ptMeas, double ptTrue, const double*) const;
  typedef double (Parametrization::*PdfPtTrue)(double ptTrue, const double*) const;
  typedef double (Parametrization::*PdfPtTrueError)(double ptTrue, const double*, const double*, const std::vector<int>&) const;
  typedef double (Parametrization::*PdfResp)(double r, double ptTrue, const double*) const;
  typedef double (Parametrization::*PdfRespError)(double r, double ptTrue, const double*, const double*, const std::vector<int>&) const;
  typedef double (Parametrization::*PdfDijetAsym)(double a, double ptTrue, const double*) const;

 public:
  SmearFunction(PdfPtMeas pdfPtMeasJet1, PdfPtMeas pdfPtMeasJet2,
		PdfPtTrue pdfPtTrue, PdfPtTrueError pdfPtTrueError,
		PdfResp pdfResp, PdfRespError pdfRespError,
		PdfDijetAsym pdfDijetAsym,
		int parIdx, int nPars, double *firstPar, const std::vector<bool>& isFixedPar,
		const std::vector<int>& covIdx, double *firstCov, const Parametrization *p)
    : pdfPtMeasJet1_(pdfPtMeasJet1),
    pdfPtMeasJet2_(pdfPtMeasJet2),
    pdfPtTrue_(pdfPtTrue),
    pdfPtTrueError_(pdfPtTrueError),
    pdfResp_(pdfResp),
    pdfRespError_(pdfRespError),
    pdfDijetAsym_(pdfDijetAsym),
    parIdx_(parIdx),
    nPars_(nPars),
    firstPar_(firstPar),
    isFixedPar_(isFixedPar),
    covIdx_(covIdx),
    firstCov_(firstCov),
    param_(p) {};
    
    
    int nPars() const { return nPars_;}
    int parIdx() const { return parIdx_; }
    double par(int i) const { assert( i >= 0 && i < nPars() ); return firstPar_[i]; }
    bool isFixedPar(int i) const { assert( i >= 0 && i < nPars() ); return isFixedPar_[i]; }

    double pdfPtMeasJet1(double ptMeas, double ptTrue) const {
      return (param_->*pdfPtMeasJet1_)(ptMeas,ptTrue,firstPar_);
    }
    double pdfPtMeasJet2(double ptMeas, double ptTrue) const {
      return (param_->*pdfPtMeasJet2_)(ptMeas,ptTrue,firstPar_);
    }
    double pdfPtTrue(double ptTrue) const {
      return (param_->*pdfPtTrue_)(ptTrue,firstPar_);
    }
    double pdfPtTrueError(double ptTrue) const {
      return (param_->*pdfPtTrueError_)(ptTrue,firstPar_,firstCov_,covIdx_);
    }
    double pdfResp(double r, double ptTrue) const {
      return (param_->*pdfResp_)(r,ptTrue,firstPar_);
    }
    double pdfRespError(double r, double ptTrue) const {
      return (param_->*pdfRespError_)(r,ptTrue,firstPar_,firstCov_,covIdx_);
    }
    double pdfDijetAsym(double a, double ptTrue) const {
      return (param_->*pdfDijetAsym_)(a,ptTrue,firstPar_);
    }

    void changeParBase(double* oldpar, double* newpar) {
      firstPar_ += newpar - oldpar;
    }

    double *par() const { return firstPar_; }


 private:
    const PdfPtMeas pdfPtMeasJet1_;
    const PdfPtMeas pdfPtMeasJet2_;
    const PdfPtTrue pdfPtTrue_;
    const PdfPtTrueError pdfPtTrueError_;
    const PdfResp pdfResp_;
    const PdfRespError pdfRespError_;
    const PdfDijetAsym pdfDijetAsym_;

    const int parIdx_;
    const int nPars_;
    double *firstPar_;
    const std::vector<bool> isFixedPar_;

    const std::vector<int> covIdx_;
    const double *firstCov_;

    const Parametrization* param_;
};
#endif
