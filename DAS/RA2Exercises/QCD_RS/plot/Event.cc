//system
#include <iostream>
#include <cmath>
#include <cstdlib>

//root includes
#include "TFile.h"
#include "TTree.h"

//own includes
#include "Event.h"

float Event::MHT() const
{
  double mhtx=0., mhty=0;
  for (unsigned i=0; i<(unsigned)NrecoJet; ++i)
    if (recoJetPt[i]>30. && fabs(recoJetEta[i])<5.0 ) {
      mhtx += recoJetPx[i];
      mhty += recoJetPy[i];
    }  
  return sqrt(mhtx*mhtx+mhty*mhty);      
}

float Event::HT() const
{
  double ht=0.;
  for (unsigned i=0; i<(unsigned)NrecoJet; ++i)
    if (recoJetPt[i]>50. && fabs(recoJetEta[i])<2.5 )
      ht += recoJetPt[i];
  return ht;      
}


/// Event IO ==============================================================================

void Event::CopyEvent(const Event& o){
  fChain =        o.fChain;
  RunNr =         o.RunNr;
  EvtNr =	  o.EvtNr;
  LumiB =	  o.LumiB;
  EvtWgt =	  o.EvtWgt;
  //NrecoJetGen =	  o.NrecoJetGen;
  //recoMetGen =	  o.recoMetGen;
  //recoMetGenPhi = o.recoMetGenPhi;
  NrecoJet =	  o.NrecoJet;
  recoMetCal =	  o.recoMetCal;
  recoMetCalPhi = o.recoMetCalPhi;
  NVtx = 	  o.NVtx;
  NrecoMu =	  o.NrecoMu;
  NrecoEle =	  o.NrecoEle;
  NrecoPho =	  o.NrecoPho;

  int f= sizeof(float);
/*
  memcpy(&recoJetGenPx, &o.recoJetGenPx,  f*sizeof(&recoJetGenPx)); 
  memcpy(&recoJetGenPy, &o.recoJetGenPy,  f*sizeof(&recoJetGenPy)); 
  memcpy(&recoJetGenPz, &o.recoJetGenPz,  f*sizeof(&recoJetGenPz)); 
  memcpy(&recoJetGenE,  &o.recoJetGenE,   f*sizeof(&recoJetGenE));  
  memcpy(&recoJetGenPt, &o.recoJetGenPt,  f*sizeof(&recoJetGenPt)); 
  memcpy(&recoJetGenPhi,&o.recoJetGenPhi, f*sizeof(&recoJetGenPhi));
  memcpy(&recoJetGenEta,&o.recoJetGenEta, f*sizeof(&recoJetGenEta));
*/
  memcpy(&recoJetPx,    &o.recoJetPx,     f*sizeof(&recoJetPx));    
  memcpy(&recoJetPy,    &o.recoJetPy,     f*sizeof(&recoJetPy));    
  memcpy(&recoJetPz,    &o.recoJetPz,     f*sizeof(&recoJetPz));    
  memcpy(&recoJetE,     &o. recoJetE,     f*sizeof(&recoJetE));	  
  memcpy(&recoJetPt,    &o.recoJetPt,     f*sizeof(&recoJetPt));    
  memcpy(&recoJetPhi,   &o.recoJetPhi,    f*sizeof(&recoJetPhi));   
  memcpy(&recoJetEta,   &o.recoJetEta,    f*sizeof(&recoJetEta));   
  memcpy(&VtxZ,         &o.VtxZ,          f*sizeof(&VtxZ));	  
  memcpy(&recoMuQ,      &o.  recoMuQ,     f*sizeof(&recoMuQ));	  
  memcpy(&recoMuPx,     &o. recoMuPx,     f*sizeof(&recoMuPx));	
  memcpy(&recoMuPy,     &o. recoMuPy,     f*sizeof(&recoMuPy));	
  memcpy(&recoMuPz,     &o. recoMuPz,     f*sizeof(&recoMuPz));	
  memcpy(&recoMuEn,     &o. recoMuEn,     f*sizeof(&recoMuEn));	
  memcpy(&recoMuPt,     &o. recoMuPt,     f*sizeof(&recoMuPt));	
  memcpy(&recoMuPhi,    &o.recoMuPhi,     f*sizeof(&recoMuPhi));  
  memcpy(&recoMuEta,    &o.recoMuEta,     f*sizeof(&recoMuEta));  
  memcpy(&recoEleQ,     &o.recoEleQ,      f*sizeof(&recoEleQ));	
  memcpy(&recoElePx,    &o.recoElePx,     f*sizeof(&recoElePx));  
  memcpy(&recoElePy,    &o.recoElePy,     f*sizeof(&recoElePy));  
  memcpy(&recoElePz,    &o.recoElePz,     f*sizeof(&recoElePz));  
  memcpy(&recoEleEn,    &o.recoEleEn,     f*sizeof(&recoEleEn));  
  memcpy(&recoElePt,    &o.recoElePt,     f*sizeof(&recoElePt));  
  memcpy(&recoElePhi,   &o.recoElePhi,    f*sizeof(&recoElePhi)); 
  memcpy(&recoEleEta,   &o.recoEleEta,    f*sizeof(&recoEleEta)); 
  memcpy(&recoPhoPx,    &o.recoPhoPx,     f*sizeof(&recoPhoPx));  
  memcpy(&recoPhoPy,    &o.recoPhoPy,     f*sizeof(&recoPhoPy));  
  memcpy(&recoPhoPz,    &o.recoPhoPz,     f*sizeof(&recoPhoPz));  
  memcpy(&recoPhoEn,    &o.recoPhoEn,     f*sizeof(&recoPhoEn));  
  memcpy(&recoPhoPt,    &o.recoPhoPt,     f*sizeof(&recoPhoPt));  
  memcpy(&recoPhoPhi,   &o.recoPhoPhi,    f*sizeof(&recoPhoPhi)); 
  memcpy(&recoPhoEta,   &o.recoPhoEta,    f*sizeof(&recoPhoEta)); 
}


void MakePseudoEvents( std::string file, std::vector<Event*>& evts)
{
  
  //Make pseudo events for testing //@@TODO: Replace with ntuple reading code
  srand(1234567890);
  for (int i=0; i<1000; ++i) {
    Event * evt = new Event;
    evt->EvtNr = i;
    evt->RunNr = 1234567890;
    evt->recoJetPx[ 0]=400.*fabs((float) rand() / ((float)RAND_MAX+1)) - 200.;
    evt->recoJetPy[ 0]=400.*fabs((float) rand() / ((float)RAND_MAX+1)) - 200.;
    evt->recoJetPz[ 0]=1000.*fabs((float) rand() / ((float)RAND_MAX+1)) - 500.;
    evt->recoJetPt[ 0]=sqrt( evt->recoJetPx[0]*evt->recoJetPx[0] + evt->recoJetPy[0]*evt->recoJetPy[0] );
    evt->recoJetE[  0]=sqrt( evt->recoJetPx[0]*evt->recoJetPx[0] + evt->recoJetPy[0]*evt->recoJetPy[0]  + evt->recoJetPz[0]*evt->recoJetPz[0] );
    //evt->recoJetPhi[0]=3.1415*fabs((float) rand() / ((float)RAND_MAX+1)) ;
    //evt->recoJetEta[0]=4.0*fabs((float) rand() / ((float)RAND_MAX+1)) - 2.0;

    evt->recoJetPx[ 1]=300.*fabs((float) rand() / ((float)RAND_MAX+1)) - 150.;
    evt->recoJetPy[ 1]=300.*fabs((float) rand() / ((float)RAND_MAX+1)) - 150.;
    evt->recoJetPz[ 1]=800.*fabs((float) rand() / ((float)RAND_MAX+1)) - 400.;
    evt->recoJetPt[ 1]=sqrt( evt->recoJetPx[1]*evt->recoJetPx[1] + evt->recoJetPy[1]*evt->recoJetPy[1] );
    evt->recoJetE[  1]=sqrt( evt->recoJetPx[1]*evt->recoJetPx[1] + evt->recoJetPy[1]*evt->recoJetPy[1]  + evt->recoJetPz[1]*evt->recoJetPz[1] );
    //evt->recoJetPhi[1]=3.1415*fabs((float) rand() / ((float)RAND_MAX+1)) ;
    //evt->recoJetEta[1]=4.0*fabs((float) rand() / ((float)RAND_MAX+1)) - 2.0;

    evt->recoJetPx[ 2]=200.*fabs((float) rand() / ((float)RAND_MAX+1)) -100.;
    evt->recoJetPy[ 2]=200.*fabs((float) rand() / ((float)RAND_MAX+1)) -100.;
    evt->recoJetPz[ 2]=500.*fabs((float) rand() / ((float)RAND_MAX+1)) -250.;
    evt->recoJetPt[ 2]=sqrt( evt->recoJetPx[2]*evt->recoJetPx[2] + evt->recoJetPy[2]*evt->recoJetPy[2] );
    evt->recoJetE[  2]=sqrt( evt->recoJetPx[2]*evt->recoJetPx[2] + evt->recoJetPy[2]*evt->recoJetPy[2]  + evt->recoJetPz[3]*evt->recoJetPz[3] );
    //evt->recoJetPhi[2]=3.1415*fabs((float) rand() / ((float)RAND_MAX+1)) ;
    //evt->recoJetEta[2]=4.0*fabs((float) rand() / ((float)RAND_MAX+1)) - 2.0;

    evt->recoJetPx[ 3]=100.*fabs((float) rand() / ((float)RAND_MAX+1)) -50.;
    evt->recoJetPy[ 3]=100.*fabs((float) rand() / ((float)RAND_MAX+1)) -50.;
    evt->recoJetPz[ 3]=300.*fabs((float) rand() / ((float)RAND_MAX+1)) -150.;
    evt->recoJetPt[ 3]=sqrt( evt->recoJetPx[3]*evt->recoJetPx[3] + evt->recoJetPy[3]*evt->recoJetPy[3] );
    evt->recoJetE[  3]=sqrt( evt->recoJetPx[3]*evt->recoJetPx[3] + evt->recoJetPy[3]*evt->recoJetPy[3]  + evt->recoJetPz[3]*evt->recoJetPz[3] );
    //evt->recoJetPhi[3]=3.1415*fabs((float) rand() / ((float)RAND_MAX+1)) ;
    //evt->recoJetEta[3]=4.0*fabs((float) rand() / ((float)RAND_MAX+1)) - 2.0;

    evt->fChain = 0;
    evt->NrecoJet = 4;
    //evt->NrecoJetGen =
    evt->NVtx=evt->NrecoMu=evt->NrecoEle=evt->NrecoPho=0;
    evts.push_back( evt );
    
  }
  std::cout << "Created "<<evts.size()<<" pseudo events '"<<file<<"'."<<std::endl;
}

void ReadEvents( std::string file, std::vector<Event*>& evts)
{
  TChain * tchain = new TChain( "AnaTree" );
  tchain->Add( file.c_str() );
  Event *evt = new Event;
  evt->Init( tchain );
  Long64_t nentries = evt->fChain->GetEntries();
  for (Long64_t i = 0; i < nentries; i++){
     evt->fChain->GetEntry(i);
     Event *copy = new Event;
     copy->CopyEvent( *evt );
     evts.push_back( copy );
  }
  std::cout << "Read "<<evts.size()<<" events from file '"<<file<<"'."<<std::endl;
  delete evt;
}


void WriteEvents(std::string file, std::vector<Event*>& evts)
{
  TFile* o_file = new TFile(file.c_str(),"RECREATE");
  TTree* myTree = new TTree("AnaTree","");
  Event *evt = new Event;
    
  //Failed to use the NtupleReader->Init(TTree*) function for this; need to define the branches by hand:
  //--------Tree Branches---------
  //----event info
  myTree->Branch("RunNr",&evt->RunNr,"RunNr/I");
  myTree->Branch("EvtNr",&evt->EvtNr,"EvtNr/I");
  myTree->Branch("LumiB",&evt->LumiB,"LumiB/I");
  myTree->Branch("EvtWgt",&evt->EvtWgt,"EvtWgt/F");

/*
  //---generator level jets
  myTree->Branch("NrecoJetGen",&evt->NrecoJetGen,"NrecoJetGen/I");
  myTree->Branch("recoJetGenPx",evt->recoJetGenPx,"recoJetGenPx[NrecoJetGen]/F");
  myTree->Branch("recoJetGenPy",evt->recoJetGenPy,"recoJetGenPy[NrecoJetGen]/F");
  myTree->Branch("recoJetGenPz",evt->recoJetGenPy,"recoJetGenPz[NrecoJetGen]/F");
  myTree->Branch("recoJetGenE",evt->recoJetGenE,"recoJetGenE[NrecoJetGen]/F");
  myTree->Branch("recoJetGenPt",evt->recoJetGenPt,"recoJetGenPt[NrecoJetGen]/F");
  myTree->Branch("recoJetGenPhi",evt->recoJetGenPhi,"recoJetGenPhi[NrecoJetGen]/F");
  myTree->Branch("recoJetGenEta",evt->recoJetGenEta,"recoJetGenEta[NrecoJetGen]/F");
  //---generator level MET
  myTree->Branch("recoMetGen",&evt->recoMetGen,"recoMetGen/F");
  myTree->Branch("recoMetGenPhi",&evt->recoMetGenPhi,"recoMetGenPhi/F");
*/

  //---reco level jets
  myTree->Branch("NrecoJet",&evt->NrecoJet,"NrecoJet/I");
  myTree->Branch("recoJetPx",evt->recoJetPx,"recoJetPx[NrecoJet]/F");
  myTree->Branch("recoJetPy",evt->recoJetPy,"recoJetPy[NrecoJet]/F");
  myTree->Branch("recoJetPz",evt->recoJetPy,"recoJetPz[NrecoJet]/F");
  myTree->Branch("recoJetE",evt->recoJetE,"recoJetE[NrecoJet]/F");
  myTree->Branch("recoJetPt",evt->recoJetPt,"recoJetPt[NrecoJet]/F");
  myTree->Branch("recoJetPhi",evt->recoJetPhi,"recoJetPhi[NrecoJet]/F");
  myTree->Branch("recoJetEta",evt->recoJetEta,"recoJetEta[NrecoJet]/F");

  //---MET
  myTree->Branch("recoMetCal",&evt->recoMetCal,"recoMetCal/F");
  myTree->Branch("recoMetCalPhi",&evt->recoMetCalPhi,"recoMetCalPhi/F");
  
  //---vertex
  myTree->Branch("NVtx",&evt->NVtx,"NVtx/I");
  myTree->Branch("VtxZ",evt->VtxZ,"VtxZ[NVtx]/F");
  
  //----muons
  myTree->Branch("NrecoMu",&evt->NrecoMu,"NrecoMu/I");    
  myTree->Branch("recoMuQ",evt->recoMuQ,"recoMuQ[NrecoMu]/F");   
  myTree->Branch("recoMuPx",evt->recoMuPx,"recoMuPx[NrecoMu]/F");
  myTree->Branch("recoMuPy",evt->recoMuPy,"recoMuPy[NrecoMu]/F");
  myTree->Branch("recoMuPz",evt->recoMuPz,"recoMuPz[NrecoMu]/F");
  myTree->Branch("recoMuEn",evt->recoMuEn,"recoMuEn[NrecoMu]/F");
  myTree->Branch("recoMuPt",evt->recoMuPt,"recoMuPt[NrecoMu]/F");
  myTree->Branch("recoMuPhi",evt->recoMuPhi,"recoMuPhi[NrecoMu]/F");
  myTree->Branch("recoMuEta",evt->recoMuEta,"recoMuEta[NrecoMu]/F");
  //---electrons
  myTree->Branch("NrecoEle",&evt->NrecoEle,"NrecoEle/I");    
  myTree->Branch("recoEleQ",evt->recoEleQ,"recoEleQ[NrecoEle]/F");   
  myTree->Branch("recoElePx",evt->recoElePx,"recoElePx[NrecoEle]/F");
  myTree->Branch("recoElePy",evt->recoElePy,"recoElePy[NrecoEle]/F");
  myTree->Branch("recoElePz",evt->recoElePz,"recoElePz[NrecoEle]/F");
  myTree->Branch("recoEleEn",evt->recoEleEn,"recoEleEn[NrecoEle]/F");
  myTree->Branch("recoElePt",evt->recoElePt,"recoElePt[NrecoEle]/F");
  myTree->Branch("recoElePhi",evt->recoElePhi,"recoElePhi[NrecoEle]/F");
  myTree->Branch("recoEleEta",evt->recoEleEta,"recoEleEta[NrecoEle]/F");
  //---photons
  myTree->Branch("NrecoPho",&evt->NrecoPho,"NrecoPho/I");    
  myTree->Branch("recoPhoPx",evt->recoPhoPx,"recoPhoPx[NrecoPho]/F");
  myTree->Branch("recoPhoPy",evt->recoPhoPy,"recoPhoPy[NrecoPho]/F");
  myTree->Branch("recoPhoPz",evt->recoPhoPz,"recoPhoPz[NrecoPho]/F");
  myTree->Branch("recoPhoEn",evt->recoPhoEn,"recoPhoEn[NrecoPho]/F");
  myTree->Branch("recoPhoPt",evt->recoPhoPt,"recoPhoPt[NrecoPho]/F");
  myTree->Branch("recoPhoPhi",evt->recoPhoPhi,"recoPhoPhi[NrecoPho]/F");
  myTree->Branch("recoPhoEta",evt->recoPhoEta,"recoPhoEta[NrecoPho]/F");

  //fill TTree
  for (std::vector<Event*>::const_iterator it=evts.begin(); it!=evts.end(); ++it){
    //std::cout<< evt->recoJetPt;
    //std::cout << "Calling CopyEvent for # " << (*it)->EvtNr << std::endl;
    evt->CopyEvent( **it );
    //std::cout<<" <> "<<evt->recoJetPt
    //         <<", jetpt[0]="<<evt->recoJetPt[0]
    //         <<", jetpt[1]="<<evt->recoJetPt[1]
    //         <<", vec-jetpt[1]="<<(*it)->recoJetPt[1]
    //	     <<std::endl;
    myTree->Fill();
  }
  
  //writing stuff to file
  o_file->cd();
  myTree->Write();
  delete myTree;

  if (o_file!=0){
    o_file->Write();
    delete o_file;
  }
  delete evt;
  std::cout << "Wrote "<<evts.size()<<" events to file '"<<file<<"'."<<std::endl;
}

