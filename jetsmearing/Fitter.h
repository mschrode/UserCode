// $Id: Fitter.h,v 1.4 2009/05/05 13:58:37 mschrode Exp $

#ifndef JS_FITTER_H
#define JS_FITTER_H

#include <string>
#include <vector>

#include "TF1.h"

#include "Event.h"
#include "NJetEvent.h"
#include "PhotonJetEvent.h"
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

  Fitter(Data& data, double min, double max,
	 const std::string& model, const std::vector<double>& par)
    : mData(data), mMin(min), mMax(max), mModel(model), mPar(par), mNBadEvts(0)
    {
      assert( GetNPar() == GetNPar(model) );
    }
    ~Fitter() {};

    void Fit();
    int GetNPar() const  { return static_cast<int>(mPar.size()); }
    int GetNPar(const std::string& model) const;
    double GetPar(int i) { assert( i >=0 && i < GetNPar() ); return mPar.at(i); }
    void SetPar(int i, double par) { assert( i >=0 && i < GetNPar() ); mPar.at(i) = par; }
    void SetPar(const std::vector<double>& par) { assert( static_cast<int>(par.size()) == GetNPar(mModel) ); mPar = par; }
    TF1 * GetTF1() const;


  private:
    Data mData;
    double mMin;
    double mMax;
    std::string mModel;
    std::vector<double> mPar;
    mutable int mNBadEvts;

    double NLogLSum(std::vector<double>& grad);
    double NLogL(Event * evt) const;
    double NLogLDiJet(DiJetEvent * dijet) const;
    double NLogLPhotonJet(PhotonJetEvent * photonjet) const;

    double PDF(double x) const;
    double PDFFermiTail(double x, const std::vector<double>& par) const;
    double PDFTwoGauss(double x, const std::vector<double>& par) const;
    double PDFThreeGauss(double x, const std::vector<double>& par) const;
    double PDFExpTail(double x, const std::vector<double>& par) const;
  };
}
#endif
