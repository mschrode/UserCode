// $Id: ControlPlots.h,v 1.4 2009/05/06 12:17:32 mschrode Exp $

#ifndef JS_CONTROLPLOTS_H
#define JS_CONTROLPLOTS_H

#include <string>

#include "TF1.h"
#include "TObject.h"

#include "Jet.h"
#include "Event.h"

namespace js
{
  //!  \brief Generates validation plots
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 19:02:46 CEST 2009
  //!  $Id: ControlPlots.h,v 1.4 2009/05/06 12:17:32 mschrode Exp $
  // --------------------------------------------------
  class ControlPlots
  {
  public:
    ControlPlots(const Data& data);
    ~ControlPlots();

    void PlotDijets() const;
    void PlotPhotonJets() const;
    void PlotResponse(TObject * pdf) const;
    void SetFileNameSuffix(std::string suffix) { mFileNameSuffix = suffix; }
    void SetRespBinning(int nbins, double min, double max) { mRespNBins = nbins; mRespMin = min; mRespMax = max;}


  private:
    Data         mData;
    int          mDijetNBins;
    double       mDijetMin;
    double       mDijetMax;
    int          mPhotonJetNBins;
    double       mPhotonJetMin;
    double       mPhotonJetMax;
    int          mRespNBins;
    double       mRespMin;
    double       mRespMax;
    std::string  mRootFileName;
    std::string  mFileNameSuffix;
    std::string  mDir;

    void SetGStyle() const;
  };
}
#endif
