//!  \brief   Container class for jet correction factors
//
//    $Id: JetMETCorFactorsFactory.h,v 1.4 2009/11/26 13:23:40 stadie Exp $
//   
#ifndef JETMETCORFACTORSFACTORY_H
#define JETMETCORFACTORSFACTORY_H

#include "CorFactorsFactory.h"

class FactorizedJetCorrector;

class JetMETCorFactorsFactory : public CorFactorsFactory
{
 public:
  JetMETCorFactorsFactory(const std::string& name, const std::string& files);
  ~JetMETCorFactorsFactory();

  CorFactors* create(const Jet* j);
 private:
  FactorizedJetCorrector* cor_;
  
  class Register {
  public:
    Register();
  private:
    JetMETCorFactorsFactory* create(const std::string& name, const std::string& files) const;
  };
  static Register register_;
};



#endif
