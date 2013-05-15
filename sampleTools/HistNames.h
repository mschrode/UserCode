// $Id: $

#ifndef HISTNAMES
#define HISTNAMES

#include "TString.h"


namespace sampleTools {
  class HistNames {
  public:
    TString ptAve() const { return "hPtAve"; }
    TString ptAveGenBin() const { return "hPtAveGenBin"; }
    TString ptAsym() const { return "hPtAsym"; }
    TString ptAsymGenBin() const { return "hPtAsymGenBin"; }
    TString ptAsymAbs() const { return "hPtAsymAbs"; }
    TString ptAsymAbsGenBin() const { return "hPtAsymAbsGenBin"; }
    TString respGenBin() const { return "hRespGenBin"; }
    TString deltaPhi() const { return "hDeltaPhi"; }
    TString pt1() const { return "hPtJet1"; }
    TString pt2() const { return "hPtJet2"; }    
    TString eta1() const { return "hEtaJet1"; }
    TString eta2() const { return "hEtaJet2"; }
    TString pt3Rel() const { return "hPt3Rel"; }
    TString respVsPtGen() const { return "hRespVsPtGen"; }
    TString respVsEtaGen() const { return "hRespVsEtaGen"; }
  };
}
#endif
