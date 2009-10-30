//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 14 17:47:13 2007 by ROOT version 5.14/00b
// from TTree CalibTree/
// found on file: calib.root
//////////////////////////////////////////////////////////

#ifndef GammaJetSel_h
#define GammaJetSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <iostream>
class GammaJetSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leave types
   Int_t           NobjTowCal;
   Int_t           NobjTrack;
   Int_t           TowId[200];   //[NobjTowCal]
   Int_t           TowId_phi[200];   //[NobjTowCal]
   Int_t           TowId_eta[200];   //[NobjTowCal]
   Float_t         TowEt[200];   //[NobjTowCal]
   Float_t         TowEta[200];   //[NobjTowCal]
   Float_t         TowPhi[200];   //[NobjTowCal]
   Float_t         TowE[200];   //[NobjTowCal]
   Float_t         TowEm[200];   //[NobjTowCal]
   Float_t         TowHad[200];   //[NobjTowCal]
   Float_t         TowOE[200];   //[NobjTowCal]
   Int_t           TrackId[200];   //[NobjTrack]
   Int_t           TrackTowId[200];   //[NobjTrack]
   Int_t           TrackTowIdPhi[200];   //[NobjTrack]
   Int_t           TrackTowIdEta[200];   //[NobjTrack]
   Float_t         TrackPt[200];   //[NobjTrack]
   Float_t         TrackEta[200];   //[NobjTrack]
   Float_t         TrackEtaOut[200];   //[NobjTrack]
   Float_t         TrackPhi[200];   //[NobjTrack]
   Float_t         TrackPhiOut[200];   //[NobjTrack]
   Float_t         TrackDR[200];   //[NobjTrack]
   Float_t         TrackDROut[200];   //[NobjTrack]
   Float_t         TrackP[200];   //[NobjTrack]
   Float_t         TrackEMC1[200];   //[NobjTrack]
   Float_t         TrackEMC3[200];   //[NobjTrack]
   Float_t         TrackEMC5[200];   //[NobjTrack]
   Float_t         TrackHAC1[200];   //[NobjTrack]
   Float_t         TrackHAC3[200];   //[NobjTrack]
   Float_t         TrackHAC5[200];   //[NobjTrack]
   Float_t         TrackChi2[200];   //[NobjTrack]
   Int_t           TrackNHits[200];   //[NobjTrack]
   Bool_t          TrackQualityL[200];   //[NobjTrack]
   Bool_t          TrackQualityT[200];   //[NobjTrack]
   Bool_t          TrackQualityHP[200];   //[NobjTrack]
   Float_t         MuDE[200];   //[NobjTrack]
   Float_t         MuDR[200];   //[NobjTrack]

   Float_t         JetCalPt;
   Float_t         JetCalPhi;
   Float_t         JetCalEta;
   Float_t         JetCalEt;
   Float_t         JetCalE; 
   Float_t         JetGenPt;
   Float_t         JetGenPhi;
   Float_t         JetGenEta;
   Float_t         JetGenEt;
   Float_t         JetGenE;
   Float_t         JetCorrZSP;
   Float_t         JetCorrL2;
   Float_t         JetCorrL3;
   Float_t         JetCorrJPT;
   Float_t         JetCorrL2L3;
   Float_t         JetCorrL2L3JPT;
   Float_t         MetCal;
   Float_t         MetCalPhi;
   Float_t         MetCalSum;
   Float_t         NonLeadingJetPt;
   Float_t         PhotonPt;
   Float_t         PhotonPhi;
   Float_t         PhotonEta;
   Float_t         PhotonEt;
   Float_t         PhotonE;
   Float_t         GenPhotonPt;
   Float_t         GenPhotonPhi;
   Float_t         GenPhotonEta;
   Float_t         GenPhotonEt;
   Float_t         GenPhotonE;
   //Int_t           ProcessID;
   Float_t         EventWeight;   

   // List of branches
   TBranch        *b_NobjTowCal;   //!
   TBranch        *b_NobjTrack;   //!
   TBranch        *b_TowId;   //!
   TBranch        *b_TowId_phi;   //!
   TBranch        *b_TowId_eta;   //!
   TBranch        *b_TowEt;   //!
   TBranch        *b_TowEta;   //!
   TBranch        *b_TowPhi;   //!
   TBranch        *b_TowE;   //!
   TBranch        *b_TowEm;   //!
   TBranch        *b_TowHad;   //!
   TBranch        *b_TowOE;   //!
   TBranch        *b_TrackId;   //
   TBranch        *b_TrackTowId;   //
   TBranch        *b_TrackTowIdPhi;   //
   TBranch        *b_TrackTowIdEta;   //
   TBranch        *b_TrackPt;   //
   TBranch        *b_TrackEta;   //
   TBranch        *b_TrackEtaOut;   //
   TBranch        *b_TrackPhi;   //
   TBranch        *b_TrackPhiOut;   //
   TBranch        *b_TrackDR;   //
   TBranch        *b_TrackDROut;   //
   TBranch        *b_TrackP;   //
   TBranch        *b_TrackEMC1;   //
   TBranch        *b_TrackEMC3;   //
   TBranch        *b_TrackEMC5;   //
   TBranch        *b_TrackHAC1;   //
   TBranch        *b_TrackHAC3;   //
   TBranch        *b_TrackHAC5;   //
   TBranch        *b_TrackChi2;   //
   TBranch        *b_TrackNHits;   //
   TBranch        *b_TrackQualityL;   //!
   TBranch        *b_TrackQualityT;   //!
   TBranch        *b_TrackQualityHP;   //!
   TBranch        *b_MuDE;   //
   TBranch        *b_MuDR;   //
   TBranch        *b_JetCalPt;   //!
   TBranch        *b_JetCalPhi;   //!
   TBranch        *b_JetCalEta;   //!
   TBranch        *b_JetCalEt;   //!
   TBranch        *b_JetCalE;   //! 
   TBranch        *b_JetCorrZSP;   //!
   TBranch        *b_JetCorrL2;   //!
   TBranch        *b_JetCorrL3;   //!
   TBranch        *b_JetCorrJPT;   //! 
   TBranch        *b_JetCorrL2L3;   //!
   TBranch        *b_JetCorrL2L3JPT;   //!
   TBranch        *b_JetGenPt;   //!
   TBranch        *b_JetGenPhi;   //!
   TBranch        *b_JetGenEta;   //!
   TBranch        *b_JetGenEt;   //!
   TBranch        *b_JetGenE;   //!
   TBranch        *b_MetCal;   //!
   TBranch        *b_MetCalPhi;   //!
   TBranch        *b_MetCalSum;   //!
   TBranch        *b_NonLeadingJetPt;  //!
   TBranch        *b_PhotonPt;   //!
   TBranch        *b_PhotonPhi;   //!
   TBranch        *b_PhotonEta;   //!
   TBranch        *b_PhtonEt;   //!
   TBranch        *b_PhotonE;   //!
   TBranch        *b_GenPhotonPt;   //!
   TBranch        *b_GenPhotonPhi;   //!
   TBranch        *b_GenPhotonEta;   //!
   TBranch        *b_GenPhtonEt;   //!
   TBranch        *b_GenPhotonE;   //!
   //TBranch        *b_ProcessID;   //!
   TBranch        *b_EventWeight;   //!

   //GammaJetSel(TTree * /*tree*/ =0) : ProcessID(-1), EventWeight(1) { }
   GammaJetSel(TTree * /*tree*/ =0) : EventWeight(1) { }
   virtual ~GammaJetSel() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   //   ClassDef(GammaJetSel,0);
};

#endif

#ifdef GammaJetSel_cxx
void GammaJetSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   // Set branch addresses and branch pointers
   if (!tree) { 
     std::cerr<<"tree is bad!"<<std::endl;
     return;
   }
   fChain = tree;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("NobjTowCal", &NobjTowCal, &b_NobjTowCal);
   fChain->SetBranchAddress("NobjTrack", &NobjTrack, &b_NobjTrack);
   fChain->SetBranchAddress("TowId", TowId, &b_TowId);
   fChain->SetBranchAddress("TowId_phi", TowId_phi, &b_TowId_phi);
   fChain->SetBranchAddress("TowId_eta", TowId_eta, &b_TowId_eta);
   fChain->SetBranchAddress("TowEt", TowEt, &b_TowEt);
   fChain->SetBranchAddress("TowEta", TowEta, &b_TowEta);
   fChain->SetBranchAddress("TowPhi", TowPhi, &b_TowPhi);
   fChain->SetBranchAddress("TowE", TowE, &b_TowE);
   fChain->SetBranchAddress("TowEm", TowEm, &b_TowEm);
   fChain->SetBranchAddress("TowHad", TowHad, &b_TowHad);
   fChain->SetBranchAddress("TowOE", TowOE, &b_TowOE);
   fChain->SetBranchAddress("TrackId", TrackId, &b_TrackId);
   fChain->SetBranchAddress("TrackTowId", TrackTowId, &b_TrackTowId);
   fChain->SetBranchAddress("TrackTowIdPhi", TrackTowIdPhi, &b_TrackTowIdPhi);
   fChain->SetBranchAddress("TrackTowIdEta", TrackTowIdEta, &b_TrackTowIdEta);
   fChain->SetBranchAddress("TrackPt", TrackPt, &b_TrackPt);
   fChain->SetBranchAddress("TrackEta", TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackEtaOut", TrackEtaOut, &b_TrackEtaOut);
   fChain->SetBranchAddress("TrackPhi", TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackPhiOut", TrackPhiOut, &b_TrackPhiOut);
   fChain->SetBranchAddress("TrackDR", TrackDR, &b_TrackDR);
   fChain->SetBranchAddress("TrackDROut", TrackDROut, &b_TrackDROut);
   fChain->SetBranchAddress("TrackP", TrackP, &b_TrackP);
   fChain->SetBranchAddress("TrackEMC1", TrackEMC1, &b_TrackEMC1);
   fChain->SetBranchAddress("TrackEMC3", TrackEMC3, &b_TrackEMC3);
   fChain->SetBranchAddress("TrackEMC5", TrackEMC5, &b_TrackEMC5);
   fChain->SetBranchAddress("TrackHAC1", TrackHAC1, &b_TrackHAC1);
   fChain->SetBranchAddress("TrackHAC3", TrackHAC3, &b_TrackHAC3);
   fChain->SetBranchAddress("TrackHAC5", TrackHAC5, &b_TrackHAC5);
   fChain->SetBranchAddress("TrackChi2", TrackChi2, &b_TrackChi2);
   fChain->SetBranchAddress("TrackNHits", TrackNHits, &b_TrackNHits);
   fChain->SetBranchAddress("TrackQualityL", TrackQualityL, &b_TrackQualityL);
   fChain->SetBranchAddress("TrackQualityT", TrackQualityT, &b_TrackQualityT);
   fChain->SetBranchAddress("TrackQualityHP", TrackQualityHP, &b_TrackQualityHP);
   fChain->SetBranchAddress("MuDR", MuDR, &b_MuDR);
   fChain->SetBranchAddress("MuDE", MuDE, &b_MuDE);
   fChain->SetBranchAddress("JetCalPt", &JetCalPt, &b_JetCalPt);
   fChain->SetBranchAddress("JetCalPhi", &JetCalPhi, &b_JetCalPhi);
   fChain->SetBranchAddress("JetCalEta", &JetCalEta, &b_JetCalEta);
   fChain->SetBranchAddress("JetCalEt", &JetCalEt, &b_JetCalEt);
   fChain->SetBranchAddress("JetCalE", &JetCalE, &b_JetCalE);
   fChain->SetBranchAddress("JetCorrZSP", &JetCorrZSP, &b_JetCorrZSP);
   fChain->SetBranchAddress("JetCorrL2", &JetCorrL2, &b_JetCorrL2);
   fChain->SetBranchAddress("JetCorrL3", &JetCorrL3, &b_JetCorrL3);
   fChain->SetBranchAddress("JetCorrJPT", &JetCorrJPT, &b_JetCorrJPT);
   fChain->SetBranchAddress("JetCorrL2L3", &JetCorrL2L3, &b_JetCorrL2L3); 
   fChain->SetBranchAddress("JetCorrL2L3JPT", &JetCorrL2L3JPT, &b_JetCorrL2L3JPT); 
   fChain->SetBranchAddress("JetGenPt", &JetGenPt, &b_JetGenPt);
   fChain->SetBranchAddress("JetGenPhi", &JetGenPhi, &b_JetGenPhi);
   fChain->SetBranchAddress("JetGenEta", &JetGenEta, &b_JetGenEta);
   fChain->SetBranchAddress("JetGenEt", &JetGenEt, &b_JetGenEt);
   fChain->SetBranchAddress("JetGenE", &JetGenE, &b_JetGenE);
   fChain->SetBranchAddress("MetCal", &MetCal, &b_MetCal);
   fChain->SetBranchAddress("MetCalPhi", &MetCalPhi, &b_MetCalPhi);
   fChain->SetBranchAddress("MetCalSum", &MetCalSum, &b_MetCalSum);
   fChain->SetBranchAddress("NonLeadingJetPt", &NonLeadingJetPt, &b_NonLeadingJetPt);
   fChain->SetBranchAddress("PhotonPt", &PhotonPt, &b_PhotonPt);
   fChain->SetBranchAddress("PhotonPhi", &PhotonPhi, &b_PhotonPhi);
   fChain->SetBranchAddress("PhotonEta", &PhotonEta, &b_PhotonEta);
   fChain->SetBranchAddress("PhotonEt", &PhotonEt, &b_PhtonEt);
   fChain->SetBranchAddress("PhotonE", &PhotonE, &b_PhotonE);
   fChain->SetBranchAddress("GenPhotonPt", &GenPhotonPt, &b_GenPhotonPt);
   fChain->SetBranchAddress("GenPhotonPhi", &GenPhotonPhi, &b_GenPhotonPhi);
   fChain->SetBranchAddress("GenPhotonEta", &GenPhotonEta, &b_GenPhotonEta);
   fChain->SetBranchAddress("GenPhotonEt", &GenPhotonEt, &b_GenPhtonEt);
   fChain->SetBranchAddress("GenPhotonE", &GenPhotonE, &b_GenPhotonE);
   //fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
}

Bool_t GammaJetSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef GammaJetSel_cxx
