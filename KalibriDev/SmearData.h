// $Id: SmearData.h,v 1.4 2010/01/29 20:54:22 mschrode Exp $

#ifndef SmearData_h
#define SmearData_h

#include "CalibData.h"
#include "SmearFunction.h"


//!  \brief Abstract base class for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 15:24:49 CEST 2009
//!  $Id: SmearData.h,v 1.4 2010/01/29 20:54:22 mschrode Exp $
// --------------------------------------------------
class SmearData : public Event {
 public:
  SmearData(DataType type, Measurement * mess, double truth, double ptHat, double weight, const SmearFunction& pdf);
  virtual ~SmearData() { delete mess_; }

  //!  \brief Get the negative log-likelihood of this event
  //!  \return The negative log-likelihood of this event
  // --------------------------------------------------
  virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  virtual void printInitStats() const = 0;

  virtual void ChangeParAddress(double* oldpar, double* newpar) { pdf_.changeParBase(oldpar,newpar); }
  virtual Measurement * GetMess() const { return mess_; }
  virtual double GetTruth() const { return kTruth_; }
  virtual DataType GetType() const { return kType_; }
  virtual double GetWeight() const { return weight_; }

  double respPar(int i) { return pdf_.respPar(i); }
  double respPDF(double r, double pt) const { return pdf_.respPDF(r,pt); }
  double respPDFError(double r, double pt) const { return pdf_.respError(r,pt); }
  double truthPDF(double pt) const { return pdf_.truthPDF(pt); }
  double truthPDFError(double pt) const { return pdf_.truthError(pt); }

  virtual void setWeight(double w) { weight_ = w; } 

  virtual double chi2_plots() const { return 0.; }                 //!< Dummy, no functionality
  virtual double GetParametrizedMess() const { return 0.; }        //!< Dummy, no functionality
  virtual void updateError() { }                                   //!< Dummy, no functionality


 protected:
  SmearFunction pdf_;                    //!< Response pdf
  Measurement * mess_;                       //!< The jet
  //  mutable Measurement tmpMeas_;


 private:
  const double    kTruth_;                     //!< Truth
  const DataType  kType_;                      //!< Event type
  double          weight_;                     //!< Event weight
};

#endif
