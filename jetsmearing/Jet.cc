#include "Jet.h"


namespace js
{
  // --------------------------------------------------
  Jet::Jet()
    : mPMeas(new TLorentzVector(0.,0.,0.,0.)),
      mPTrue(new TLorentzVector(0.,0.,0.,0.)) {}



  // --------------------------------------------------
  Jet::Jet(const Jet& jet)
    : mPMeas( new TLorentzVector( *(jet.mPMeas) ) ),
      mPTrue( new TLorentzVector( *(jet.mPTrue) ) ) {}



  // --------------------------------------------------
  Jet::Jet(const TLorentzVector& pTrue, const TLorentzVector& pMeas)
    : mPMeas(new TLorentzVector(pMeas)),
      mPTrue(new TLorentzVector(pTrue)) {}



  // --------------------------------------------------
  Jet::~Jet()
  {
    delete mPMeas;
    delete mPTrue;
  }



  // --------------------------------------------------
  Jet& Jet::operator=(const Jet& jet)
  {
    if(this != &jet)
      {
	delete mPMeas;
	delete mPTrue;
	
	mPMeas = new TLorentzVector( *(jet.mPMeas) );
	mPTrue = new TLorentzVector( *(jet.mPTrue) );
      }

    return *this;
  }
}

