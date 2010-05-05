//
//    $Id: JetMETCorFactorsFactory.cc,v 1.7 2010/04/29 13:29:41 stadie Exp $
//   
#include "JetMETCorFactorsFactory.h"
#include "CorFactors.h"
#include "Jet.h"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"

#include <iostream>

JetMETCorFactorsFactory::JetMETCorFactorsFactory(const std::string& name,
						 const std::string& files)
  : CorFactorsFactory(name)
{
  cor_ = new FactorizedJetCorrector("L2Relative:L3Absolute",files);
  //std::cout << "created JetMETCorFactorsFactory: " << name << '\n';
}

JetMETCorFactorsFactory::~JetMETCorFactorsFactory()
{
  delete cor_;
}


CorFactors* JetMETCorFactorsFactory::create(const Jet* j)
{
  cor_->setJetEta(j->eta());
  cor_->setJetPt(j->pt()); 
  cor_->setJetE(j->E());
  cor_->setJetPhi(j->phi());
  cor_->setJetEMF(j->EmEt()/j->pt()); 
  
  std::vector<float> levels = cor_->getSubCorrections();

  return new CorFactors(1.0,
			levels[0],
			levels[1]/levels[0],
			1.0,1.0,0.0,0.0);			
}
JetMETCorFactorsFactory::Register JetMETCorFactorsFactory::register_;

JetMETCorFactorsFactory::Register::Register() 
{
  //create("Summer09_7TeV_AK5Calo","JetMETObjects/data/Summer09_7TeV_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_L3Absolute_AK5Calo.txt");
  //create("Summer09_AK5Calo","JetMETObjects/data/Summer09_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_L3Absolute_AK5Calo.txt");
  create("Summer09_7TeV_ReReco332_AK5Calo","JetMETObjects/data/Summer09_7TeV_ReReco332_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_ReReco332_L3Absolute_AK5Calo.txt");  
  create("Spring10_AK5Calo","JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt");
}

JetMETCorFactorsFactory* JetMETCorFactorsFactory::Register::create(const std::string& name, const std::string& files) const
{
  JetMETCorFactorsFactory* jmcff = 0;
  try {
    jmcff = new JetMETCorFactorsFactory(name,files);
  } 
  catch(std::exception& e) {
    std::cout << "...failed to create " << name << ":\n";
    std::cout << "     " << e.what() << std::endl;
    jmcff = 0;
  } 
  return jmcff;
}
