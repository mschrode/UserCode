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


bool Rebalance_KinFitter( const Event* evt, Event* rebalanced,  JetResolution * JetRes=0)
{
   //// Interface to KinFitter
   TKinFitter* myFit = new TKinFitter();
   std::vector<TLorentzVector*> lvec_m;
   std::vector<TMatrixD*> covMat_m;
   std::vector<TFitParticleEtEtaPhi*> fitted;
   std::vector<TFitParticleEtEtaPhi*> measured;
   
   bool result = true;
   double MHTx_high = 0;
   double MHTy_high = 0;
   std::vector<Jet> rebjets;

   //// Fill measured particles to vector
   for (int i = 0; i<evt->NrecoJet; ++i) {

      if (evt->recoJetPt[i] >= 10. ) {  //use all jets >= 20. GeV
         MHTx_high -= evt->recoJetPx[i];
         MHTy_high -= evt->recoJetPy[i];

         // The particles before fitting
         double tmppx, tmppy, tmppz, tmpe;
         tmppx = evt->recoJetPx[i];
         tmppy = evt->recoJetPy[i];
         tmppz = evt->recoJetPz[i];
         tmpe  = evt->recoJetE[i];

         TLorentzVector* lv = new TLorentzVector(tmppx, tmppy, tmppz, tmpe);
         lvec_m.push_back(lv);
         TMatrixD* cM = new TMatrixD(3, 3);
         (*cM)(0, 0) = JetRes->JetResolution_Pt2(evt->recoJetPt[i], evt->recoJetEta[i], i);
         (*cM)(1, 1) = JetRes->JetResolution_Eta2(evt->recoJetE[i], evt->recoJetEta[i]);
         (*cM)(2, 2) = JetRes->JetResolution_Phi2(evt->recoJetE[i], evt->recoJetEta[i]);
         covMat_m.push_back(cM);
         char name[10];
         sprintf(name, "jet%i", i);
         TFitParticleEtEtaPhi* jet1 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
         measured.push_back(jet1);
         TFitParticleEtEtaPhi* jet2 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
         fitted.push_back(jet2);
         myFit->addMeasParticle(fitted.back());
      }
      else {
	  rebjets.push_back(Jet( 
	       evt->recoJetPt[i],
	       evt->recoJetPx[i],evt->recoJetPy[i],evt->recoJetPz[i],
	       evt->recoJetE[i],evt->recoJetPhi[i],evt->recoJetEta[i]));
      }
   }

   //// Add momentum constraints
   double MET_constraint_x = 0.;
   double MET_constraint_y = 0.;
   TFitConstraintEp* momentumConstr1 = new TFitConstraintEp("px", "px", 0, TFitConstraintEp::pX, MET_constraint_x);
   TFitConstraintEp* momentumConstr2 = new TFitConstraintEp("py", "py", 0, TFitConstraintEp::pY, MET_constraint_y);
   for (unsigned int i = 0; i < fitted.size(); ++i) {
      momentumConstr1->addParticle(fitted.at(i));
      momentumConstr2->addParticle(fitted.at(i));
   }
   myFit->addConstraint(momentumConstr1);
   myFit->addConstraint(momentumConstr2);

   //// Set fit parameters
   myFit->setVerbosity(0);
   myFit->setMaxNbIter(100);
   myFit->setMaxF(0.01 * 2);
   myFit->setMaxDeltaS(1.e-3);
   myFit->fit();
   int status = myFit->getStatus();

   double chi2 = 0;
   double F = 0;
   double prob = 0;
   if (status == 0) {
      chi2 = myFit->getS();
      F = myFit->getF();
      int dof = myFit->getNDF();
      prob = TMath::Prob(chi2, dof);
      if (prob < 1.e-2) result = false;
   } else {
      chi2 = 99999;
      prob = 0;
      F = 99999;
      result = false;
   }
   //   h_fitProb->Fill(prob);  /// plot fit probability

   //// Get the output of KinFitter
   rebalanced->CopyEvent( *evt );
   //rebalanced->NrecoJet = measured.size();
   for (unsigned int i = 0; i < measured.size(); ++i) {
          rebjets.push_back(Jet( 
	         fitted.at(i)->getCurr4Vec()->Pt(),
	         fitted.at(i)->getCurr4Vec()->Px(),fitted.at(i)->getCurr4Vec()->Py(),
		 fitted.at(i)->getCurr4Vec()->Pz(),fitted.at(i)->getCurr4Vec()->E(),
		 fitted.at(i)->getCurr4Vec()->Phi(),fitted.at(i)->getCurr4Vec()->Eta()));
   }
   
   std::sort(rebjets.begin(), rebjets.end());
   for (std::vector<Jet>::const_iterator jet=rebjets.begin(); jet!=rebjets.end(); ++jet){
      rebalanced->recoJetPt[jet-rebjets.begin()] = jet->Pt;
      rebalanced->recoJetPx[jet-rebjets.begin()] = jet->Px;
      rebalanced->recoJetPy[jet-rebjets.begin()] = jet->Py;
      rebalanced->recoJetPz[jet-rebjets.begin()] = jet->Pz;
      rebalanced->recoJetE[ jet-rebjets.begin()] = jet->E;
      rebalanced->recoJetPhi[jet-rebjets.begin()] = jet->Phi;
      rebalanced->recoJetEta[jet-rebjets.begin()] = jet->Eta;
   }  


   delete myFit;
   for (unsigned int i = 0; i < measured.size(); ++i) {
      delete lvec_m.at(i);
      delete covMat_m.at(i);
      delete measured.at(i);
      delete fitted.at(i);
   }
   delete momentumConstr1;
   delete momentumConstr2;

   return result;
}

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
//    if ( Rebalance( *it, rebalanced, JetRes ) ) {
    if ( Rebalance_KinFitter( *it, rebalanced, JetRes ) ) {
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
