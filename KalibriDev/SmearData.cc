// $Id: SmearData.cc,v 1.4 2010/01/29 20:54:22 mschrode Exp $

#include "SmearData.h"

#include "CalibData.h"


//!  \brief Constructor
//!  \param type Data type
//!  \param mess The measurement
//!  \param truth The truth of the measurement
//!  \param weight Event weight
//!  \param respPDF Response probability density
// --------------------------------------------------
SmearData::SmearData(DataType type, Measurement * mess, double truth, double ptHat, double weight, const SmearFunction& pdf)
  : Event(),
    pdf_(pdf),
    mess_(mess),
    kTruth_(truth),
    kType_(type),
    weight_(weight) {
  ptHat_ = ptHat;
};



// //!  \brief Response pdf
// //!  \param r Response
// //!  \param pt pt
// //!  \return The probability density of the response \p r 
// //!          for an event with true pt \p pt
// // --------------------------------------------------
// double SmearData::respPDF(double r, double pt) const {
//   tmpMeas_.E  = r;
//   tmpMeas_.pt = pt;
//   return respPDF_(&tmpMeas_);
// }



// // --------------------------------------------------
// double SmearData::respPDFSigma(double r, double pt) const {
//   tmpMeas_.E  = r;
//   tmpMeas_.pt = pt;
//   return respPDF_.sigma(&tmpMeas_);
// }



