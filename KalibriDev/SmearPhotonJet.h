// $Id: SmearPhotonJet.h,v 1.3 2009/07/23 13:44:41 mschrode Exp $

#ifndef SmearPhotonJet_h
#define SmearPhotonJet_h

#include "CalibData.h"
#include "SmearData.h"


//!  \brief Photon-jet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearPhotonJet.h,v 1.3 2009/07/23 13:44:41 mschrode Exp $
// --------------------------------------------------
class SmearPhotonJet : public SmearData {
 public:
  SmearPhotonJet(TMeasurement * mess, double photonpt, double weight, const Function& respPDF)
    : SmearData(TypeSmearPhotonJet,mess,photonpt,weight,respPDF) {};
  ~SmearPhotonJet() {};

  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  virtual void printInitStats() const {};  //!< No functionality yet
  virtual double ptHat() const { return GetTruth(); }
};

#endif
