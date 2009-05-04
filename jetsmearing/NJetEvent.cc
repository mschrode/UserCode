// $Id: $

#include "NJetEvent.h"

#include <iostream>

namespace js
{
  //!  \brief Destructor
  // --------------------------------------------------
  NJetEvent::~NJetEvent()
  {
    for(std::vector<js::Jet*>::iterator it = mJets.begin();
	it != mJets.end(); it++)
      {
	delete *it;
      }
    mJets.clear();
  }



  // --------------------------------------------------
  void NJetEvent::Print() const
  {
    std::cout << "\n";
    for(int i = 0; i < GetNJets(); i++)
      {
	std::cout << "Jet " << i << ":" << std::flush;
	std::cout << "  (PtTrue) " << PtTrue(i) << std::flush;
	std::cout << "  (PtMeas) " << PtMeas(i) << std::flush;
	std::cout << "  (PhiTrue) " << PhiTrue(i) << std::flush;
	std::cout << "  (PhiMeas) " << PhiMeas(i) << std::endl;
      }
  }
}
