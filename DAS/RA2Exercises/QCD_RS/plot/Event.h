#ifndef EVENT_H
#define EVENT_H

#include "NtupleSelector.h"

#include <string>
#include <vector>

class Event : public NtupleSelector {
 public:
   void CopyEvent(const Event& o);
   float HT() const;
   float MHT() const;

 private:
};

void ReadEvents( std::string file, std::vector<Event*>& evts);
void MakePseudoEvents( std::string file, std::vector<Event*>& evts);
void WriteEvents(std::string file, std::vector<Event*>& evts);


#endif
