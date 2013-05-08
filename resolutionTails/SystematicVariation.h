// $Id: $

#ifndef RESOLUTION_TAILS_SYSTEMATIC_VARIATION
#define RESOLUTION_TAILS_SYSTEMATIC_VARIATION

#include <vector>

#include "TString.h"

namespace resolutionTails {

  class SystematicVariation {
  public:
    static std::vector<SystematicVariation*> create();

    virtual ~SystematicVariation() {};

    virtual TString name() const { return ""; }
    virtual TString label() const { return ""; }

    virtual bool isNominal() const { return false; }
    virtual bool isCoreUp() const { return false; }
    virtual bool isCoreDn() const { return false; }
    virtual bool isExtrapolation() const { return false; }
    virtual bool isClosure() const { return false; }
    virtual bool isPUUp() const { return false; }
    virtual bool isPUDn() const { return false; }
  };


  class VariationNominal : public SystematicVariation {
  public:
    VariationNominal() : SystematicVariation() {};

    TString name() const { return "nominal"; }
    TString label() const { return ""; }
    bool isNominal() const { return true; }
  };

  class VariationCoreUp : public SystematicVariation {
  public:
    VariationCoreUp() : SystematicVariation() {};

    TString name() const { return "core up"; }
    TString label() const { return "VarCoreUp"; }
    bool isCoreUp() const { return true; }
  };

  class VariationCoreDn : public SystematicVariation {
  public:
    VariationCoreDn() : SystematicVariation() {};

    TString name() const { return "core down"; }
    TString label() const { return "VarCoreDn"; }
    bool isCoreDn() const { return true; }
  };

  class VariationExtrapolation : public SystematicVariation {
  public:
    VariationExtrapolation() : SystematicVariation() {};

    TString name() const { return "extrapolation"; }
    TString label() const { return "VarExtra"; }
    bool isExtrapolation() const { return true; }
  };

  class VariationClosure : public SystematicVariation {
  public:
    VariationClosure() : SystematicVariation() {};

    TString name() const { return "closure"; }
    TString label() const { return "VarClosure"; }
    bool isClosure() const { return true; }
  };

  class VariationPUUp : public SystematicVariation {
  public:
    VariationPUUp() : SystematicVariation() {};

    TString name() const { return "PU up"; }
    TString label() const { return "VarPUUp"; }
    bool isPUUp() const { return true; }
  };

  class VariationPUDn : public SystematicVariation {
  public:
    VariationPUDn() : SystematicVariation() {};

    TString name() const { return "PU down"; }
    TString label() const { return "VarPUDn"; }
    bool isPUDn() const { return true; }
  };





  std::vector<resolutionTails::SystematicVariation*> SystematicVariation::create() {
    std::vector<resolutionTails::SystematicVariation*> vars;
    vars.push_back(new VariationNominal());
    vars.push_back(new VariationCoreUp());
    vars.push_back(new VariationCoreDn());
    vars.push_back(new VariationExtrapolation());
    vars.push_back(new VariationClosure());
    //    vars.push_back(new VariationPUUp());
    //    vars.push_back(new VariationPUDn());
    
    return vars;
  }
}
#endif
