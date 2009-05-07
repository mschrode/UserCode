// $Id: Event.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $

#ifndef JS_EVENT_H
#define JS_EVENT_H

#include <string>
#include <vector>

namespace js
{
  //!  \brief An abstract event class
  //!  \author Matthias Schroeder
  //!  \date Wed Apr 29 11:13:55 CEST 2009
  //!  $Id: Event.h,v 1.2 2009/05/04 14:35:04 mschrode Exp $
  // --------------------------------------------------
  class Event
  {
  public:
    Event();
    ~Event();

    virtual std::string Type() const = 0;
    virtual void Print() const = 0;
    int Status() const { return mStatus; }           //!< 0 is good
    void SetStatus(int status) { mStatus = status; }

  private:
    bool mStatus;
  };


  typedef std::vector<Event*> Data;
  typedef std::vector<Event*>::const_iterator DataIt;
}
#endif
