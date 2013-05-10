// $Id: Uncertainty.h,v 1.1 2013/05/10 13:22:18 mschrode Exp $

#ifndef RESOLUTION_TAILS_UNCERTAINTY
#define RESOLUTION_TAILS_UNCERTAINTY

#include <cstdlib>

#include "TString.h"

namespace resolutionTails {

  class Uncertainty {
  public:
    enum SystematicVariation { Nominal, CoreUp, CoreDn, Extrapolation, Closure, PUUp, PUDn };

    static TString id(Uncertainty::SystematicVariation var);
    static TString name(Uncertainty::SystematicVariation var);
  };


  // Unique id string to be used in file names
  TString Uncertainty::id(Uncertainty::SystematicVariation var) {
    TString res = "";
    if( var == Nominal ) res = "";
    else if( var == CoreUp ) res = "VarCoreUp";
    else if( var == CoreDn ) res = "VarCoreDn";
    else if( var == Extrapolation ) res = "VarExtra";
    else if( var == Closure ) res = "VarClosure";
    else if( var == PUUp ) res = "VarPUUp";
    else if( var == PUDn ) res = "VarPUDn";
    else {
      std::cerr << "ERROR: undefined systematic variation " << var << std::endl;
      exit(-1);
    }

    return res;
  }


  // Name for printout
  TString Uncertainty::name(Uncertainty::SystematicVariation var) {
    TString res = "";
    if( var == Nominal ) res = "nominal";
    else if( var == CoreUp ) res = "core up";
    else if( var == CoreDn ) res = "core down";
    else if( var == Extrapolation ) res = "extrapolation";
    else if( var == Closure ) res = "closure";
    else if( var == PUUp ) res = "PU up";
    else if( var == PUDn ) res = "PU down";
    else {
      std::cerr << "ERROR: undefined systematic variation " << var << std::endl;
      exit(-1);
    }

    return res;
  }
}
#endif
