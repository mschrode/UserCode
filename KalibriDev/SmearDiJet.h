// $Id: SmearDiJet.h,v 1.7 2010/04/13 13:38:24 mschrode Exp $

#ifndef SmearDiJet_h
#define SmearDiJet_h

#include "SmearData.h"
#include "SmearFunction.h"
#include "Jet.h"


//!  \brief Dijet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearDiJet.h,v 1.7 2010/04/13 13:38:24 mschrode Exp $
// --------------------------------------------------
class SmearDiJet : public SmearData {
 public:
  SmearDiJet(Jet * jet1,
	     Jet * jet2,
	     Jet * jet3,
	     double ptHat,
	     double weight,
	     const SmearFunction& pdf,
	     double min,
	     double max,
	     double eps,
	     int niter);
  ~SmearDiJet();

  //  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  virtual void printInitStats() const;

  const Jet * jet1() const { return static_cast<Jet*>(mess_); }
  const Jet * jet2() const { return jet2_; }
  const Jet * jet3() const { return jet3_; }

  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel) const { return pdf_.pdfPtMeasJet1(ptMeas,ptTrue,pt3Rel); }
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel) const { return pdf_.pdfPtMeasJet2(ptMeas,ptTrue,pt3Rel); }

  //! Get dijet pt \f$ \frac{1}{2} (p^{1}_{T} + p^{2}_{T}) \f$
  double dijetPt() const { return 0.5 * (jet1()->pt() + jet2()->pt()); } 
  double relJet3Pt() const { return jet3()->pt() / dijetPt(); }


 private:
  const int    kMaxNIter_;   //!< Max number of iterations in integration
  const double kEps_;        //!< Integration precision for convergence
  const double kMin_;        //!< Minimum of truth pdf
  const double kMax_;        //!< Maximum of truth pdf

  Jet * jet2_; //!< Second jet
  Jet * jet3_; //!< Third jet
};
#endif
