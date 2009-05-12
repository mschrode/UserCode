// $Id: ControlPlots.h,v 1.5 2009/05/08 12:13:28 mschrode Exp $

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
  //!  $Id: ControlPlots.h,v 1.5 2009/05/08 12:13:28 mschrode Exp $
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
    void SetLog(bool setlog) { mSetLog = setlog; }
    void SetGrid(bool setgrid) { mSetGrid = setgrid; }


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
    bool         mSetLog;
    bool         mSetGrid;

    void SetGStyle() const;
  };
}
#endif
