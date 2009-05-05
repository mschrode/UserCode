// $Id: Jet.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $

#ifndef JS_JET_H
#define JS_JET_H

#include <vector>

#include "TLorentzVector.h"

namespace js
{
  //!  \brief A jet
  //!
  //!  A jet consists of two Lorentz vectors, one
  //!  holding the truth information and one the
  //!  measured information.
  //!
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 13:22:36 CEST 2009
  //!  $Id: Jet.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $
  // --------------------------------------------------
  class Jet
  {
  public:
    //!  \brief Default constructor
    //!
    //!  All quantities are initialized with 0.
    // --------------------------------------------------
    Jet();

    //!  \brief Copy constructor
    // --------------------------------------------------
    Jet(const Jet& jet);

    //!  \brief Constructor
    //!
    //!  Lorentz vector of truth and measurement are
    //!  specified.
    //!  
    //!  \param pTrue Lorentz vector with truth
    //!  \param pTrue Lorentz vector with measurement
    // --------------------------------------------------
    Jet(const TLorentzVector& pTrue, const TLorentzVector& pMeas);
    ~Jet();

    Jet& operator=(const Jet& jet);

    //!  \brief Get x component of measured momentum
    //!  \return x component of measured momentum
    // --------------------------------------------------
    double PxMeas() const { return mPMeas->Px(); }

    //!  \brief Get y component of measured momentum
    //!  \return y component of measured momentum
    // --------------------------------------------------
    double PyMeas() const { return mPMeas->Py(); }

    //!  \brief Get z component of measured momentum
    //!  \return z component of measured momentum
    // --------------------------------------------------
    double PzMeas() const { return mPMeas->Pz(); }

    //!  \brief Get transverse component of measured momentum
    //!  \return Transverse component of measured momentum
    // --------------------------------------------------
    double PtMeas() const { return mPMeas->Perp(); }

    //!  \brief Get polar angle of measured momentum
    //!  \return Polar angle of measured momentum
    // --------------------------------------------------
    double PhiMeas() const { return mPMeas->Phi(); }

    //!  \brief Get azimuthal angle of measured momentum
    //!  \return Azimuthal angle of measured momentum
    // --------------------------------------------------
    double ThetaMeas() const { return mPMeas->Theta(); }

    //!  \brief Get pseudorapidity of measured momentum
    //!  \return Pseudorapidity of measured momentum
    // --------------------------------------------------
    double EtaMeas() const { return mPMeas->Eta(); }

    //!  \brief Get x component of true momentum
    //!  \return x component of true momentum
    // --------------------------------------------------
    double PxTrue() const { return mPTrue->Px(); }

    //!  \brief Get y component of true momentum
    //!  \return y component of true momentum
    // --------------------------------------------------
    double PyTrue() const { return mPTrue->Py(); }

    //!  \brief Get z component of true momentum
    //!  \return z component of true momentum
    // --------------------------------------------------
    double PzTrue() const { return mPTrue->Pz(); }

    //!  \brief Get transverse component of true momentum
    //!  \return Transverse component of true momentum
    // --------------------------------------------------
    double PtTrue() const { return mPTrue->Perp(); }

    //!  \brief Get polar angle of true momentum
    //!  \return Polar angle of true momentum
    // --------------------------------------------------
    double PhiTrue() const { return mPTrue->Phi(); }

    //!  \brief Get azimuthal angle of true momentum
    //!  \return Azimuthal angle of true momentum
    // --------------------------------------------------
    double ThetaTrue() const { return mPTrue->Theta(); }

    //!  \brief Get pseudorapidity of true momentum
    //!  \return Pseudorapidity of true momentum
    // --------------------------------------------------
    double EtaTrue() const { return mPTrue->Eta(); }

    
  private:
    TLorentzVector * mPMeas;  //!< The measured 4-momentum
    TLorentzVector * mPTrue;  //!< The true 4-momentum
  };


  typedef std::vector<js::Jet*>::const_iterator JetIt; //!< Iterator of vector of jets
}
#endif
