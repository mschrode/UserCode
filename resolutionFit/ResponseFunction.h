#ifndef RESPONSE_FUNCTION_H
#define RESPONSE_FUNCTION_H

#include <vector>

namespace resolutionFit {
  class ResponseFunction {
  public: 
    enum Type { Gauss, CrystalBall };

    ResponseFunction(Type type)
      : type_(type) {};

    double operator()(double x, const std::vector<double> &par) const;
    Type type() const { return type_; }
    int nPars() const;

    double pdfGauss(double x, double mean, double sigma) const;
    double pdfCrystalBall(double x, double mean, double sigma, double alpha, double n) const;

    double sigmaGauss(double pt, const std::vector<double> &par) const;


  private:
    Type type_;

    double crystalBallInt(double mean, double sigma, double alpha, double n, double min, double max) const;
  };
}
#endif
