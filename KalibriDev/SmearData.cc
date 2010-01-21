// $Id: SmearData.cc,v 1.8 2010/01/08 18:16:02 mschrode Exp $

#include "SmearData.h"

#include "CalibData.h"


//!  \brief Constructor
//!  \param type Data type
//!  \param mess The measurement
//!  \param truth The truth of the measurement
//!  \param weight Event weight
//!  \param respPDF Response probability density
// --------------------------------------------------
SmearData::SmearData(DataType type, Measurement * mess, double truth, double weight, const Function& respPDF)
  : Event(),
    respPDF_(respPDF),
    mess_(mess),
    kTruth_(truth),
    kType_(type),
    weight_(weight) {};



//!  \brief Response pdf
//!  \param r Response
//!  \param pt pt
//!  \return The probability density of the response \p r 
//!          for an event with true pt \p pt
// --------------------------------------------------
double SmearData::respPDF(double r, double pt) const {
  Measurement meas;
  meas.E  = r;
  meas.pt = pt;
  return respPDF_(&meas);
}
