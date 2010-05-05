//
//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactorsFactory.cc,v 1.3 2010/04/29 13:29:41 stadie Exp $
// 
#include "CorFactorsFactory.h"  
#include "CorFactors.h"

#include <iostream>

std::map<std::string,CorFactorsFactory*> CorFactorsFactory::map;


CorFactorsFactory::CorFactorsFactory(const std::string& name) : name_(name)
{
  static Cleaner cleaner;
  //std::cout << "creating CorFactorsFactory: " << name_ << '\n';

  map[name] = this;
}

CorFactorsFactory::~CorFactorsFactory()
{
  //std::cout << "deleting CorFactorsFactory: " << name_ << '\n';
  map[name_] = 0;
}
