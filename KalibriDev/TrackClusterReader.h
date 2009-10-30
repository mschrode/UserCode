//
//    Reader for Track Cluster Events
//
//    This class reads events according to the TrackClusterSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TrackClusterReader.h,v 1.1 2008/12/12 17:06:00 stadie Exp $
//   
#ifndef TRACKCLUSTEREADER_H
#define TRACKCLUSTERREADER_H

#include "EventReader.h"

#include <string>

#include "TrackClusterSel.h"

class TrackClusterReader : public EventReader{
 public:
  TrackClusterReader(const std::string& configfile, TParameters *p);
  virtual ~TrackClusterReader();
  int readEvents(std::vector<TData*>& data);
 private:
  TrackClusterSel trackcluster;
  double Et_cut_on_track, Et_cut_on_cluster;
  int n_trackcluster_events;
};


#endif
