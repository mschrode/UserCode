// $Id: Fitter.h,v 1.8 2009/05/12 11:54:04 mschrode Exp $

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
  //!  \brief Fit a certain response model to the data
  //!
  //!
  //!  The Fitte finds the most probable probability density
  //!  function (pdf) \f$ p \f$ of the response \f$ x \f$
  //!  compatible with the data.
  //!  The data is assumed to be correctly calibrated i.e.
  //!  the mean response is assumed to be 1.
  //!  Possible response pdf models are:
  //!
  //!   - "FermiTail":
  //!     Gaussian \f$ G \f$ plus constant contribution
  //!     from a Fermi function \f$ F \f$ which drops to
  //!     0 for \f$ x > 1 \f$:
  //!     \f[
  //!      p(x) = c \cdot G(1,\sigma) + (1-c) \cdot F(1,T)
  //!     \f]
  //!     This model has 3 parameters:
  //!      - 0: Normalization \f$ c \f$
  //!      - 1: Width \f$ \sigma \f$ of Gaussian
  //!      - 2: Temperature \f$ T \f$ of Fermi-function    
  //!
  //!   - "TwoGauss": See NGauss
  //!
  //!   - "ThreeGauss": See NGauss
  //!
  //!   - NGauss (The name of the actual models is "TwoGauss" and
  //!     "ThreeGauss"):
  //!     One central and \f$ N \f$ side Gaussians \f$ G \f$:
  //!     \f[
  //!      p(x) = (1-(c_{1}+c_{2})) \cdot G(1,\sigma_{0})
  //!             + \sum^{N}_{i=1} c_{i} \cdot G(\mu_{i},\sigma_{i})
  //!     \f]
  //!     This model has \f$ 1+3N \f$ parameters:
  //!      - 0: Width \f$ \sigma_{0} \f$ of central Gaussian
  //!      - 1: Normalization \f$ c_{1} \f$ of first side Gaussian
  //!      - 2: Mean \f$ \mu_{1} \f$ of first side Gaussian
  //!      - 3: Width \f$ \sigma_{1} \f$ of first side Gaussian
  //!      - ...
  //!
  //!   - "Hist":
  //!     The response pdf is given by a histogram.
  //!     This model has a specified number of parameters
  //!     corresponding to the number of bins in a specified
  //!     response range. Outside this range, the pdf is 0.
  //!  
  //!   - "HistGauss":
  //!     The response pdf is given by the sum of a central
  //!     Gaussian and a histogram describing the tails:
  //!     \f[
  //!      p(x) = c \cdot G(1,\sigma) + (1-c) \cdot H(b_{i})
  //!     \f]
  //!     This model has 12 parameters:
  //!      - 0: Normalization \f$ c \f$
  //!      - 1: Width \f$ \sigma \f$ of central Gaussian
  //!      - 2 - 11: Bin content of the histogram between 0 and 2
  //!
  //!  The Fitter uses the
  //!  <A HREF="http://www.desy.de/~blobel/largesc.html">
  //!  LVMINI</A> program by V. Blobel.
  // --------------------------------------------------
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
    double PDFHistGauss(double x, const std::vector<double>& par) const;

    static double FermiTail(double * x, double * par);
    static double TwoGauss(double * x, double * par);
    static double ThreeGauss(double * x, double * par);
    static double ExpTail(double * x, double * par);
  };
}
#endif
