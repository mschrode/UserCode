// $Id: SmearData.cc,v 1.6 2009/09/02 13:52:26 mschrode Exp $

#include "SmearData.h"



//!  \brief Constructor
//!  \param type Data type
//!  \param mess The measurement
//!  \param truth The truth of the measurement
//!  \param weight Event weight
//!  \param respPDF Response probability density
// --------------------------------------------------
SmearData::SmearData(DataType type, TMeasurement * mess, double truth, double weight, const Function& respPDF)
  : TData(),
    respPDF_(respPDF),
    kTruth_(truth),
    kType_(type),
    weight_(weight),
    mess_(mess) {};



//!  \brief Response pdf
//!  \param r Response
//!  \param pt pt
//!  \return The probability density of the response \p r 
//!          for an event with true pt \p pt
// --------------------------------------------------
double SmearData::respPDF(double r, double pt) const {
  TMeasurement meas;
  meas.E  = r;
  meas.pt = pt;
  return respPDF_(&meas);
}
