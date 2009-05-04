#ifndef JS_CONTROLPLOTS_H
#define JS_CONTROLPLOTS_H

#include <string>

#include "TF1.h"

#include "EventGenerator.h"
#include "Jet.h"
#include "Event.h"
#include "NJetEvent.h"

namespace js
{
  //!  \brief Generates validation plots
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 19:02:46 CEST 2009
  // --------------------------------------------------
  class ControlPlots
  {
  public:
    ControlPlots(const Data& data);
    ~ControlPlots();

    void PlotDijets() const;
    void PlotResponse(TF1 * pdf) const;


  private:
    Data         mData;
    int          mDijetNBins;
    double       mDijetMin;
    double       mDijetMax;
    int          mRespNBins;
    double       mRespMin;
    double       mRespMax;
    std::string  mRootFileName;

    void SetGStyle() const;
  };
}
#endif
