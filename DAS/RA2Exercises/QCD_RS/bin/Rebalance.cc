#include "Rebalance.h"
#include "JetResolution.h"

#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"

#include <vector>
#include <map>
#include <iostream>

//this 'struct Jet' can be used to sort the jets in an event by Pt
struct Jet{
  Jet(float Pt,float Px,float Py,float Pz,float E,float Phi,float Eta):
      Pt(Pt),Px(Px),Py(Py),Pz(Pz),E(E),Phi(Phi),Eta(Eta){}
  bool operator<(const Jet& o) const {return Pt>o.Pt;}    
  float Pt,Px,Py,Pz,E,Phi,Eta;
};

// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//
//  your code goes here:
//
bool Rebalance( const Event* evt, Event* rebalanced,  JetResolution * JetRes=0)
{
   //copy event 'evt' to the resulting event 'rebalanced'
   rebalanced->CopyEvent( *evt );

   //used to sort the jets according to Pt after rebalancing:
   std::vector<Jet> rebjets;

   /// ... rebalancing

   // At this point we have the rebalanced jets. To sort them, fill them into the 'rebjets' vector:
   for (int i = 0; i<evt->NrecoJet; ++i) {
          //example: use the unbalanced jets from the original event 'evt':
	  rebjets.push_back(Jet(evt->recoJetPt[i],evt->recoJetPx[i],evt->recoJetPy[i],
	       evt->recoJetPz[i],evt->recoJetE[i],evt->recoJetPhi[i],evt->recoJetEta[i]));
   }
   
   //sort them by Pt:
   std::sort(rebjets.begin(), rebjets.end());

   //write the sorted jets into the result event 'rebalanced'
   for (std::vector<Jet>::const_iterator jet=rebjets.begin(); jet!=rebjets.end(); ++jet){
      rebalanced->recoJetPt[jet-rebjets.begin()] = jet->Pt;
      rebalanced->recoJetPx[jet-rebjets.begin()] = jet->Px;
      rebalanced->recoJetPy[jet-rebjets.begin()] = jet->Py;
      rebalanced->recoJetPz[jet-rebjets.begin()] = jet->Pz;
      rebalanced->recoJetE[ jet-rebjets.begin()] = jet->E;
      rebalanced->recoJetPhi[jet-rebjets.begin()] = jet->Phi;
      rebalanced->recoJetEta[jet-rebjets.begin()] = jet->Eta;
   }  

   //return true if successful, false otherwise
   return true;
}
//
//
// ======================================================================================



//Get the jet resolution and call Rebalance() defined above for each event:
void Rebalance( const std::vector<Event*>& evts, std::vector<Event*>& rebalanced_events )
{
  int not_converged=0;
  std::cout<<"...reading jet energy resolutions for rebalancing"<<std::endl;
  JetResolution * JetRes = new JetResolution();
  std::cerr<<"...rebalancing: ";
  for (std::vector<Event*>::const_iterator it=evts.begin(); it!=evts.end(); ++it){
    Event * rebalanced = new Event;
    if ( Rebalance( *it, rebalanced, JetRes ) ) {
      rebalanced_events.push_back( rebalanced );
    } else {
      ++not_converged;
      delete rebalanced;
    } 
    if ((it-evts.begin())%(evts.size()/10)==0)std::cerr<<"->"<<(it-evts.begin())/(evts.size()/10)*10<<"%"; 
  }
  delete JetRes;
  std::cout << "\nSuccessfully rebalanced "<<rebalanced_events.size()<<" out of "<<evts.size()
            <<" events. Fit failed for "<<not_converged<<" events."<<std::endl;
}
