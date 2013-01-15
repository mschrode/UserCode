#ifndef SELECTION_H
#define SELECTION_H

#include "Event.h"

#include <vector>


void Cut(std::vector<Event*>& evts, float(*var)(const Event*), const char op, const float val);


#endif
