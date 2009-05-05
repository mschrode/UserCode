// $Id: $

#ifndef JS_PHOTONJET_EVENT_H
#define JS_PHOTONJET_EVENT_H

#include <string>
#include <vector>

#include "TLorentzVector.h"

#include "Event.h"
#include "Jet.h"

namespace js
{
  //!  \brief An event with one photon balancing one jet
  //!  \author Matthias Schroeder
  //!  \date Tue May  5 08:32:31 CEST 2009
  //!  $Id: $
  // --------------------------------------------------
  class PhotonJetEvent : public Event
  {
  public:
    PhotonJetEvent();
    PhotonJetEvent(const TLorentzVector& photonPt, const Jet& jet);
    ~PhotonJetEvent();

    std::string Type() const { return "PhotonJetEvent"; }

    double PxMeas() const { return mJet->PxMeas(); }
    double PyMeas() const { return mJet->PyMeas(); }
    double PtMeas() const { return mJet->PtMeas(); }
    double PhiMeas() const { return mJet->PhiMeas(); }
    double PxTrue() const { return mJet->PxTrue(); }
    double PyTrue() const { return mJet->PyTrue(); }
    double PtTrue() const { return mJet->PtTrue(); }
    double PhiTrue() const { return mJet->PhiTrue(); }
    double PxPhoton() const { return mPPhoton->Px(); }
    double PyPhoton() const { return mPPhoton->Py(); }
    double PtPhoton() const { return mPPhoton->Perp(); }
    double PhiPhoton() const { return mPPhoton->Phi(); }

    void Print() const;


  private:
    Jet * mJet;                 //!< The jet
    TLorentzVector * mPPhoton;  //!< The photon 4-momentum
  };


  typedef std::vector<PhotonJetEvent*> PhotonJetData;
  typedef std::vector<PhotonJetEvent*>::const_iterator PhotonJetDataIt;
}
#endif
