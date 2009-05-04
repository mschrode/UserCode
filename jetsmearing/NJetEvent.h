// $Id: $

#ifndef JS_NJET_EVENT_H
#define JS_NJET_EVENT_H

#include <string>
#include <vector>

#include "Event.h"
#include "Jet.h"

namespace js
{
  //!  \brief An event with n jets
  //!  \author Matthias Schroeder
  //!  \date Tue Apr 28 13:31:38 CEST 2009
  //!  $Id: $
  // --------------------------------------------------
  class NJetEvent : public Event
  {
  public:
    NJetEvent() {};
    ~NJetEvent();

    virtual std::string Type() const { return "NJetEvent"; }

    virtual void AddJet(const Jet& jet) { mJets.push_back(new Jet(jet)); }
    virtual int GetNJets() const { return static_cast<int>(mJets.size()); }

    JetIt Begin() const { return mJets.begin(); }
    JetIt End() const { return mJets.end(); }

    double PxMeas(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PxMeas(); }
    double PyMeas(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PyMeas(); }
    double PtMeas(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PtMeas(); }
    double PhiMeas(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PhiMeas(); }
    double PxTrue(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PxTrue(); }
    double PyTrue(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PyTrue(); }
    double PtTrue(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PtTrue(); }
    double PhiTrue(int i) const { assert( i >= 0 && i < GetNJets() ); return mJets.at(i)->PhiTrue(); }

    void Print() const;


  protected:
    std::vector<js::Jet*> mJets;
  };


  typedef std::vector<NJetEvent*> NJetData;
  typedef std::vector<NJetEvent*>::const_iterator NJetDataIt;



  //!  \brief An event with 2 jets
  //!  \author Matthias Schroeder
  //!  \date Wed Apr 29 14:43:45 CEST 2009
  //!  $Id: $
  // --------------------------------------------------
  class DiJetEvent : public NJetEvent
  {
  public:
    DiJetEvent() {};
    ~DiJetEvent() {};

    void AddJet(const Jet& jet) { assert( mJets.size() < 2 ); mJets.push_back(new Jet(jet)); }
    std::string Type() const { return "DiJetEvent"; }
    int GetNJets() const { return static_cast<int>(std::min(1.*(mJets.size()),2.)); }
  };


  typedef std::vector<DiJetEvent*> DiJetData;
  typedef std::vector<DiJetEvent*>::const_iterator DiJetDataIt;
}
#endif
