#ifndef RESPONSE_FUNCTION_H
#define RESPONSE_FUNCTION_H

#include <vector>

class TRandom3;
class TString;
namespace resolutionFit {
  class ResponseFunction {
  public: 
    enum Type { Gauss, CrystalBall, TruncCrystalBall };

    ResponseFunction(Type type);
    ~ResponseFunction();

    double operator()(double pt, const std::vector<double> &par) const {
      return pdf(pt,par); 
    }
    double pdf(double pt, const std::vector<double> &par) const {
      return pdf(pt,type(),par);
    }
    Type type() const { return type_; }
    TString typeLabel() const;
    int nPars() const;
    double pdfAsymmetry(double a, const std::vector<double> &par) const {
      return pdfAsymmetry(a,type(),par);
    }

    double pdfGauss(double pt, const std::vector<double> &par) const;
    double pdfCrystalBall(double pt, const std::vector<double> &par) const {
      return fCB_->pdf(pt,par[0],par[1],par[2],par[3]);
    }
    double pdfTruncCrystalBall(double pt, const std::vector<double> &par) const {
      return fCB_->truncPdf(pt,par[0],par[1],par[2],par[3],par[4]);
    }

    double random(const std::vector<double> &par) const {
      return random(type(),par);
    }
    double randomGauss(const std::vector<double> &par) const;
    double randomCrystalBall(const std::vector<double> &par) const {
      return fCB_->random(par[0],par[1],par[2],par[3]);
    }
    double randomTruncCrystalBall(const std::vector<double> &par) const {
      return fCB_->truncRandom(par[0],par[1],par[2],par[3],par[4]);
    }

    double pdfAsymmetryGauss(double a, const std::vector<double> &par) const {
      return pdfAsymmetry(a,Gauss,par);
    }
    double pdfAsymmetryCrystalBall(double a, const std::vector<double> &par) const {
      return pdfAsymmetry(a,CrystalBall,par);
    }
    double pdfAsymmetryTruncCrystalBall(double a, const std::vector<double> &par) const {
      return pdfAsymmetry(a,TruncCrystalBall,par);
    }

    double sigmaGauss(double pt, const std::vector<double> &par) const;


  private:
    class CrystalBallFunction {
    public:
      CrystalBallFunction();
      CrystalBallFunction(const std::vector<double> &par);
      ~CrystalBallFunction();
	  
      double integral(double min, double max) const {
	return integral(min,max,par_[0],par_[1],par_[2],par_[3]);
      }
      double integral(double min, double max, double mean, double sigma, double alpha, double n) const;
      double norm() const {
	return norm(par_[0],par_[1],par_[2],par_[3]);
      }
      double norm(double mean, double sigma, double alpha, double n) const;
      double pdf(double x) const {
	return pdf(x,par_[0],par_[1],par_[2],par_[3]);
      }
      double pdf(double x, double mean, double sigma, double alpha, double n) const {
	return norm(mean,sigma,alpha,n)*value(x,mean,sigma,alpha,n);
      }
      double random() const {
	return random(par_[0],par_[1],par_[2],par_[3]);
      }
      double random(double mean, double sigma, double alpha, double n) const;
      double value(double x) const {
	return value(x,par_[0],par_[1],par_[2],par_[3]);
      }
      double value(double x, double mean, double sigma, double alpha, double n) const;

 
      double truncPdf(double x, double mean, double sigma, double alpha, double n, double min) const;
      double truncRandom(double mean, double sigma, double alpha, double n, double min) const;
      double truncValue(double x, double mean, double sigma, double alpha, double n, double min) const;
      
	  
    private:
      TRandom3 *rand_;
      std::vector<double> par_;
    };


    const int nBinsIntegration_;
    const double zMinIntegration_;
    const double zMaxIntegration_;
    const double deltaZIntegration_;
    const double epsIntegration_;
    const int maxNIterIntegration_;

    Type type_;
    TRandom3 *rand_;
    CrystalBallFunction *fCB_;

    double pdf(double pt, Type type, const std::vector<double> &par) const;
    double pdfAsymmetry(double a, Type type, const std::vector<double> &par) const;
    double random(Type type, const std::vector<double> &par) const;
  };
}
#endif
