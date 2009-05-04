// $Id: Fitter.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $

#ifndef JS_FITTER_H
#define JS_FITTER_H

#include <string>
#include <vector>

#include "TF1.h"

#include "Event.h"
#include "NJetEvent.h"
#include "external.h"

namespace js
{
  class Fitter
  {
  public:
    static double FermiTail(double * x, double * par);
    static double TwoGauss(double * x, double * par);
    static double ThreeGauss(double * x, double * par);
    static double ExpTail(double * x, double * par);

  Fitter(const Data& data, double min, double max,
	 const std::string& model, const std::vector<double>& par)
    : mData(data), mMin(min), mMax(max), mModel(model), mPar(par)
    {
      assert( GetNPar() == GetNPar(model) );
    }
    ~Fitter() {};

    void Fit();
    int GetNPar() const  { return static_cast<int>(mPar.size()); }
    int GetNPar(const std::string& model) const;
    double GetPar(int i) { assert( i >=0 && i < GetNPar() ); return mPar.at(i); }
    void SetPar(int i, double par) { assert( i >=0 && i < GetNPar() ); mPar.at(i) = par; }
    TF1 * GetTF1() const;


  private:
    const Data mData;
    double mMin;
    double mMax;
    std::string mModel;
    std::vector<double> mPar;

    double NLogL(std::vector<double>& grad);
    double ProbDiJet(const DiJetEvent * dijet);
    double Prob(double x) const;
    double ProbFermiTail(double x, const std::vector<double>& par) const;
    double ProbTwoGauss(double x, const std::vector<double>& par) const;
    double ProbThreeGauss(double x, const std::vector<double>& par) const;
    double ProbExpTail(double x, const std::vector<double>& par) const;

  };
}
#endif
