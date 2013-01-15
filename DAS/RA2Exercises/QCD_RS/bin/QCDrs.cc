#include "QCDrs.h"
#include "Event.h"
#include "Selection.h"
#include "Rebalance.h"
#include "Smear.h"


//global functions, used by the 'Cut' function declared in 'Selection.h'
float HT(const Event * evt){ return evt->HT(); }
float MHT(const Event * evt){ return evt->MHT(); }
float Jet1Pt(const Event * evt){ return evt->recoJetPt[0]; }
float Jet2Pt(const Event * evt){ return evt->recoJetPt[1]; }
float Jet3Pt(const Event * evt){ return evt->recoJetPt[2]; }
float Weight(const Event * evt){ return evt->EvtWgt; }


//The main function for Rebalance and Smearing
//
// Put your rebalancing code in function 
// 'bool Rebalance(const Event*, Event*,JetResolution*)' in file 'Rebalance.cc'
// 
// Put your smearing code in function
// 'void Smear(const Event*, Event*, JetResolution*)' in file 'Smear.cc'
//
int main(int argc, char* argv[])
{
  //STL container for all events
  std::vector<Event*> events, rebalanced_events, rs_events;
  
  //Read the control sample
  std::string filename="QCDcontrol_data.root";
  ReadEvents( filename, events);
 
  //Rebalance events 'events', resulting events are 'rebalanced_events'
  Rebalance( events,  rebalanced_events );
  WriteEvents("QCDcontrol_rebalanced.root", rebalanced_events);
  
  //Smear events 'rebalanced_events' by jet energy resolution, result is 'rs_events'  
  Smear(rebalanced_events, rs_events);

  //Apply the same selection cuts as for the oroginal control sample, e.g.
  Cut(rs_events, Jet1Pt, '<', 50); 
  Cut(rs_events, Jet2Pt, '<', 50); 
  Cut(rs_events, Jet3Pt, '<', 50); 
  
  //Write the 'rs_events' to file
  WriteEvents("QCDcontrol_rs.root", rs_events);

  return 0;
}
