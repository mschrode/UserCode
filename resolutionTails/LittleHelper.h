// $Id: $

#ifndef RESOLUTION_TAILS_LITTLE_HELPER
#define RESOLUTION_TAILS_LITTLE_HELPER

#define UTILS_AS_HEADER_FILE
#include "../util/HistOps.h"



namespace resolutionTails {
  class LittleHelper {
  public:
    static void correctAsymmetryWidth(TH1* hOrig, double nSigCore, double coreScale, TH1* &hSmeared, double &width, double &smearedWidth);
  };


  // ------------------------------------------------------------------------------------
  void LittleHelper::correctAsymmetryWidth(TH1* hOrig, double nSigCore, double coreScale, TH1* &hSmeared, double &width, double &smearedWidth) {
    // Fit width of original distribution
    width = 0.;
    double widthErr = 1000.;
    if( !util::HistOps::fitCoreWidth(hOrig,nSigCore,width,widthErr) ) {
      width = hOrig->GetRMS();
    }
    if( width > 2.*hOrig->GetRMS() ) {
      width = hOrig->GetRMS();
    }
    // Smear original distribution
    util::HistOps::smearHistogram(hOrig,hSmeared,width,coreScale);
    // Width of correced distribution
    if( !util::HistOps::fitCoreWidth(hSmeared,nSigCore,smearedWidth,widthErr) ) {
      smearedWidth = hSmeared->GetRMS();
    }
    if( smearedWidth > 2.*hSmeared->GetRMS() ) {
      smearedWidth = hSmeared->GetRMS();
    }
  }
}
#endif
