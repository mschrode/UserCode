// $Id: SmearPhotonJet.h,v 1.5 2010/04/13 13:38:24 mschrode Exp $

#ifndef SmearPhotonJet_h
#define SmearPhotonJet_h

#include "Jet.h"
#include "SmearData.h"
#include "SmearFunction.h"




//!  \brief Photon-jet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearPhotonJet.h,v 1.5 2010/04/13 13:38:24 mschrode Exp $
// --------------------------------------------------
class SmearPhotonJet : public SmearData {
 public:
  SmearPhotonJet(Measurement * mess, double photonpt, double ptHat, double weight, const SmearFunction& pdf)
    : SmearData(TypeSmearPhotonJet,mess,photonpt,ptHat,weight,pdf) {};
  ~SmearPhotonJet() {};

  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  virtual void printInitStats() const {};  //!< No functionality yet
  virtual double ptHat() const { return GetTruth(); }

  const Jet * jet() const { return static_cast<Jet*>(mess_); }
};

#endif
