// $Id: PhotonJetEvent.cc,v 1.1 2009/05/05 13:58:37 mschrode Exp $

#include "PhotonJetEvent.h"

#include <iostream>

namespace js
{
  //!  \brief Default constructor
  //!
  //!  Photon and jet with 4-momenta 0
  //!  are initialized.
  // --------------------------------------------------
  PhotonJetEvent::PhotonJetEvent()
    : Event(),
      mJet(new Jet()),
      mPPhoton(new TLorentzVector()) {}



  //!  \brief Constructor
  //!
  //!  Photon and jet with 4-momenta 0
  //!  are initialized.
  // --------------------------------------------------
  PhotonJetEvent::PhotonJetEvent(const TLorentzVector& photonPt, const Jet& jet)
    : mJet(new Jet(jet)),
      mPPhoton(new TLorentzVector(photonPt)) {}



  //!  \brief Destructor
  // --------------------------------------------------
  PhotonJetEvent::~PhotonJetEvent()
  {
    delete mJet;
    delete mPPhoton;
  }



  // --------------------------------------------------
  void PhotonJetEvent::Print() const
  {
    std::cout << "\n";
    std::cout << "Photon:" << std::flush;
    std::cout << "  (PtPhoton) " << PtPhoton() << std::flush;
    std::cout << "  (PhiPhoton) " << PhiPhoton() << std::endl;
    std::cout << "Jet:" << std::flush;
    std::cout << "  (PtTrue) " << PtTrue() << std::flush;
    std::cout << "  (PtMeas) " << PtMeas() << std::flush;
    std::cout << "  (PhiTrue) " << PhiTrue() << std::flush;
    std::cout << "  (PhiMeas) " << PhiMeas() << std::endl;
  }
}
