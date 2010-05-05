// $Id: SmearData.cc,v 1.9 2010/04/13 13:38:23 mschrode Exp $

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



