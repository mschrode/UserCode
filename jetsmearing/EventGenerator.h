// $Id: EventGenerator.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $

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
  //!  Possible response models are
  //!   - Gauss+Uniform:
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
  //!   - Histogram:
  //!     The pt response is distributed according to
  //!     a specified histogram. That needs to be
  //!     set via SetRespHist(TH1F * h).
  //!     The phi response is an exponential
  //!     distribution with "decay constant" tau as specified.
  //!     The parameters are:
  //!      - 0: Scale tau of exponential
  //!     
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 15:22:42 CEST 2009
  //!  $Id: EventGenerator.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $
  // --------------------------------------------------
  class EventGenerator
  {
  public:
    EventGenerator();
    EventGenerator(double min, double max, 
		   const std::string& modelTruth, const std::vector<double>& parTruth,
		   const std::string& modelResp, const std::vector<double>& parResp);
    ~EventGenerator();

    int GetNParResp() const { return static_cast<int>(mParResp.size()); }
    double GetParResp(int i) const { assert( i >= 0 && i < GetNParResp() ); return mParResp.at(i); }
    int GetNParTruth() const { return static_cast<int>(mParTruth.size()); }
    double GetParTruth(int i) const { assert( i >= 0 && i < GetNParTruth() ); return mParTruth.at(i); }
    void SetRespHist(TH1F * h) { mHistResp = h; }
    Data GenerateDijetEvents(int n) const;
    Data GeneratePhotonJetEvents(int n) const;
    void SetModelTruth();


  private:
    TRandom3 * mRandom;
    std::string mModelTruth;
    std::vector<double> mParTruth;
    std::string mModelResp;
    std::vector<double> mParResp;
    TH1F * mHistResp;
    double mMin;
    double mMax;

    TLorentzVector GenerateTruth() const;
    TLorentzVector SmearTruth(const TLorentzVector& pTrue) const;
    void SmearTruthGaussUniform(TLorentzVector& p) const;
    void SmearTruthHistogram(TLorentzVector& p) const;
  };
}
#endif
