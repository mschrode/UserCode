//
//    Reader for Tri-Jet Events
//
//    This class reads events according to the NJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TriJetReader.h,v 1.2 2009/11/24 16:52:59 stadie Exp $
//   
#ifndef TRIJETREADER_H
#define TRIJETREADER_H

#include "EventReader.h"

#include <string>

class NJetSel;

class TriJetReader : public EventReader{
 public:
  TriJetReader(const std::string& configfile, TParameters *p);
  virtual ~TriJetReader();
  int readEvents(std::vector<Event*>& data);
 private:
  NJetSel* njet;
  double Et_cut_nplus1Jet,Rel_cut_on_nJet;
  int n_trijet_events;
};


#endif
