// $Id: Fitter.h,v 1.5 2009/05/07 15:22:23 mschrode Exp $

#ifndef JS_FITTER_H
#define JS_FITTER_H

#include <string>
#include <vector>

#include "TH1D.h"
#include "TH1F.h"
#include "TF1.h"

#include "Event.h"
#include "NJetEvent.h"
#include "PhotonJetEvent.h"
#include "external.h"


//!  \brief Test program for jet smearing method to estimate
//!         MET contribution from mismeasured QCD
//!
//!  See the talk 
//!  <A HREF="https://indico.desy.de/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=1984">
//!  "Data driven determination of smearing function"</A>
//!  by C. Sander and his
//!  <A HREF="http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/csander/ToyMC/toyMC_smearing.cpp?revision=1.3&view=markup">
//!  source code</A>.
// --------------------------------------------------
namespace js
{
  class Fitter
  {
  public:
    Fitter(Data& data, double min, double max, const std::string& model, const std::vector<double>& par);
    ~Fitter();

    void Fit();
    int GetNPar() const  { return static_cast<int>(mPar.size()); }
    int GetNPar(const std::string& model) const;
    double GetPar(int i) { assert( i >=0 && i < GetNPar() ); return mPar.at(i); }
    void SetPar(int i, double par) { assert( i >=0 && i < GetNPar() ); mPar.at(i) = par; }
    void SetPar(const std::vector<double>& par) { assert( static_cast<int>(par.size()) == GetNPar(mModel) ); mPar = par; }
    TF1 * GetTF1pdf() const;
    TH1F * GetTH1Fpdf() const;


  private:
    Data mData;                 //!< Input data for fit
    double mMin;                //!< Minimum of considered truth spectrum
    double mMax;                //!< Maximum of considered truth spectrum
    std::string mModel;         //!< Assumed response model
    std::vector<double> mPar;   //!< Fitted parameters
    double mPDFHistMin;         //!< Minimum of pdf histogram
    double mPDFHistMax;         //!< Maximum of pdf histogram
    mutable int mNBadEvts;      //!< Number of events flagged as bad
    mutable TH1D * mPDFHist;    //!< Histogrammed pdf

    double NLogLSum(std::vector<double>& grad);
    double NLogL(Event * evt) const;
    double NLogLDiJet(DiJetEvent * dijet) const;
    double NLogLPhotonJet(PhotonJetEvent * photonjet) const;

    double PDF(double x) const;
    double PDFFermiTail(double x, const std::vector<double>& par) const;
    double PDFTwoGauss(double x, const std::vector<double>& par) const;
    double PDFThreeGauss(double x, const std::vector<double>& par) const;
    double PDFExpTail(double x, const std::vector<double>& par) const;
    double PDFHist(double x, const std::vector<double>& par) const;

    static double FermiTail(double * x, double * par);
    static double TwoGauss(double * x, double * par);
    static double ThreeGauss(double * x, double * par);
    static double ExpTail(double * x, double * par);
  };
}
#endif
