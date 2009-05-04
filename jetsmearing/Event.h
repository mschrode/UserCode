#ifndef JS_EVENT_H
#define JS_EVENT_H

#include <string>
#include <vector>

namespace js
{
  //!  \brief An abstract event class
  //!  \author Matthias Schroeder
  //!  \date Wed Apr 29 11:13:55 CEST 2009
  // --------------------------------------------------
  class Event
  {
  public:
    Event();
    ~Event();

    virtual std::string Type() const = 0;
    virtual void Print() const = 0;
  };


  typedef std::vector<Event*> Data;
  typedef std::vector<Event*>::const_iterator DataIt;
}
#endif
