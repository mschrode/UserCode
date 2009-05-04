#ifndef JS_JET_H
#define JS_JET_H

#include <vector>

#include "TLorentzVector.h"

namespace js
{
  //!  \brief A jet
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 13:22:36 CEST 2009
  // --------------------------------------------------
  class Jet
  {
  public:
    Jet();
    Jet(const Jet& jet);
    Jet(const TLorentzVector& pTrue, const TLorentzVector& pMeas);
    ~Jet();

    Jet& operator=(const Jet& jet);

    double PxMeas() const { return mPMeas->Px(); }
    double PyMeas() const { return mPMeas->Py(); }
    double PzMeas() const { return mPMeas->Pz(); }
    double PtMeas() const { return mPMeas->Perp(); }
    double PhiMeas() const { return mPMeas->Phi(); }
    double ThetaMeas() const { return mPMeas->Theta(); }
    double EtaMeas() const { return mPMeas->Eta(); }

    double PxTrue() const { return mPTrue->Px(); }
    double PyTrue() const { return mPTrue->Py(); }
    double PzTrue() const { return mPTrue->Pz(); }
    double PtTrue() const { return mPTrue->Perp(); }
    double PhiTrue() const { return mPTrue->Phi(); }
    double ThetaTrue() const { return mPTrue->Theta(); }
    double EtaTrue() const { return mPTrue->Eta(); }

    
  private:
    TLorentzVector * mPMeas;
    TLorentzVector * mPTrue;
  };


  typedef std::vector<js::Jet*>::const_iterator JetIt;
}
#endif
