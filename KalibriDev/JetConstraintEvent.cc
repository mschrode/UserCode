//
//    Class for constraints on the jet correction
//
//    first version: Hartmut Stadie 2009/07/23
//    $Id: JetConstraintEvent.cc,v 1.2 2009/07/24 09:38:34 stadie Exp $
//   


#include "JetConstraintEvent.h"


//interface to Data
JetConstraintEvent::~JetConstraintEvent()
{
  for(unsigned int i = 0, njets = jets.size() ; i < njets ; ++i) {
    delete jets[i];
  }
  jets.clear();
}

void JetConstraintEvent::addJet(Jet* j) 
{
  jets.push_back(j); 
  //add parameter ids to variation  map
  const Jet::VariationColl& varcoll = j->varyParsDirectly(0.001);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    varmap[i->parid] = Variation();
  }
}


void JetConstraintEvent::ChangeParAddress(double* oldpar, double* newpar) { 
  for(unsigned int i = 0, njets = jets.size() ; i < njets ; ++i) {
    jets[i]->ChangeParAddress(oldpar,newpar);
  }
}
 
double JetConstraintEvent::chi2() const
{
  double sumet = 0;
  unsigned int njets = jets.size();
  for(unsigned int i = 0 ; i < njets ; ++i) {
    sumet += jets[i]->correctedEt(jets[i]->Et());
  }
  double chi2 = sumet / njets - truth;
  return weight * chi2 * chi2;
}

double JetConstraintEvent::chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const
{
  double sumet = 0;
  unsigned int njets = jets.size();
  for(unsigned int i = 0 ; i < njets ; ++i) {
    sumet += jets[i]->correctedEt(jets[i]->Et());
  }
  double chi2 = sumet / njets - truth;
  chi2 *= chi2;
  chi2 *= weight;
  //derivative....
  for(VarMap::iterator i = varmap.begin() ; i != varmap.end() ; ++i) {
    i->second.uppersum = sumet;
    i->second.lowersum = sumet;
  }
  for(unsigned int i = 0 ; i < njets ; ++i) {
    const Jet::VariationColl& varcoll = jets[i]->varyParsDirectly(epsilon);
    double et = jets[i]->correctedEt(jets[i]->Et());
    for(Jet::VariationCollIter j = varcoll.begin() ; j != varcoll.end() ; ++j) {
      VarMap::iterator k = varmap.find(j->parid);
      assert(k != varmap.end());
      k->second.uppersum += j->upperEt - et;
      k->second.lowersum += j->lowerEt - et;
    }
  }
  for(VarMap::iterator i = varmap.begin() ; i != varmap.end() ; ++i) {
    double temp1 = i->second.lowersum / njets - truth;
    temp1 *= temp1;
    temp1 *= weight;
    double temp2 = i->second.uppersum / njets - truth;
    temp2 *= temp2;
    temp2 *= weight;
    temp_derivative1[i->first] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->first] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative   
  }
  chi2plots = chi2;
  return chi2;
}


