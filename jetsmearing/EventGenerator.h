// $Id: EventGenerator.h,v 1.3 2009/05/05 13:58:37 mschrode Exp $

#ifndef JS_EVENTGENERATOR_H
#define JS_EVENTGENERATOR_H

#include <string>
#include <vector>

#include "TLorentzVector.h"

#include "TH1F.h"
#include "TRandom3.h"

#include "Event.h"

namespace js
{
  //!  \brief An event generator
  //!  
  //!  The transverse momentum pt and the polar angle phi
  //!  are smeared by a response function. Its parameters
  //!  can be specified.
  //!
  //!  Possible response models are
  //!   - "Gauss+Uniform":
  //!     Pt response is given by a Gaussian with mean 1 and
  //!     width as specified. For a specified fraction of
  //!     events, the response is uniformly distributed
  //!     between 0 and 1 simulating a low energy tail.
  //!     The phi response is an exponential
  //!     distribution with "decay constant" tau as specified.
  //!     The parameters are:
  //!      - 0: Width sigma of Gaussian
  //!      - 1: Fraction of uniformly distributed response
  //!           (Low energy tail)
  //!      - 2: Scale tau of exponential
  //!     (This is the default model.)
  //!   - "Histogram":
  //!     The pt response is distributed according to
  //!     a specified histogram. That needs to be
  //!     set via SetRespHist(TH1F * h).
  //!     The phi response is an exponential
  //!     distribution with "decay constant" tau as specified.
  //!     The parameters are:
  //!      - 0: Scale tau of exponential
  //!   - "TwoGauss":
  //!     The pt response is distributed according to two Gaussians.
  //!     \f[
  //!      p(x) = c \cdot G_{0}(1,\sigma_{0}) + (1-c) \cdot G_{1}(\mu_{1},\sigma_{1})
  //!     \f]
  //!     (Internally a histogram is filled from the smooth, two
  //!     Gaussian function, and the random number is generated from
  //!     the histogram.)  
  //!     The phi response is an exponential
  //!     distribution with "decay constant" tau as specified.
  //!     The parameters are:
  //!      - 0: Normalization \f$ c \f$
  //!      - 1: Width \f$ \sigma_{0} \f$ of central Gaussian
  //!      - 2: Mean \f$ \mu_{1} \f$ of second Gaussian
  //!      - 3: Width \f$ \sigma_{0} \f$ of second Gaussian
  //!     
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 15:22:42 CEST 2009
  //!  $Id: EventGenerator.h,v 1.3 2009/05/05 13:58:37 mschrode Exp $
  // --------------------------------------------------
  class EventGenerator
  {
  public:
    EventGenerator();
    EventGenerator(const std::string& modelTruth, const std::vector<double>& parTruth,
		   const std::string& modelResp, const std::vector<double>& parResp);
    ~EventGenerator();

    int GetNParResp() const { return static_cast<int>(mParResp.size()); }
    double GetParResp(int i) const { assert( i >= 0 && i < GetNParResp() ); return mParResp.at(i); }
    int GetNParTruth() const { return static_cast<int>(mParTruth.size()); }
    double GetParTruth(int i) const { assert( i >= 0 && i < GetNParTruth() ); return mParTruth.at(i); }
    void SetRespHist(TH1F * h) { mHistResp = h; }
    Data GenerateDijetEvents(int n) const;
    Data GeneratePhotonJetEvents(int n) const;
    void WriteResponseHist(const std::string& name) const;
    


  private:
    TRandom3 *             mRandom;      //!< Random number generator
    std::string            mModelTruth;  //!< Truth model (so far uniform)
    std::vector<double>    mParTruth;    //!< Parameter for truth model (so min and max)
    std::string            mModelResp;   //!< Response model
    std::vector<double>    mParResp;     //!< Parameter for response model
    TH1F *                 mHistResp;    //!< Holds histogram in case of histogramed response model

    TLorentzVector GenerateTruth() const;
    TLorentzVector SmearTruth(const TLorentzVector& pTrue) const;
    void SmearTruthGaussUniform(TLorentzVector& p) const;
    void SmearTruthHistogram(TLorentzVector& p) const;
    void SmearTruthTwoGauss(TLorentzVector& p) const;
  };
}
#endif
