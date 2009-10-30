//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 27 18:11:17 2009 by ROOT version 5.20/00
// from TTree DiJetTree/
// found on file: /scratch/current/cms/user/stadie/DiJet_Track_230_300_rereco_incomplete.root
//////////////////////////////////////////////////////////

#ifndef NJetSel_h
#define NJetSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class NJetSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           NobjTow;
   Int_t           TowId[10000];   //[NobjTow]
   Int_t           TowId_phi[10000];   //[NobjTow]
   Int_t           TowId_eta[10000];   //[NobjTow]
   Float_t         TowEt[10000];   //[NobjTow]
   Float_t         TowEta[10000];   //[NobjTow]
   Float_t         TowPhi[10000];   //[NobjTow]
   Float_t         TowE[10000];   //[NobjTow]
   Float_t         TowEm[10000];   //[NobjTow]
   Float_t         TowHad[10000];   //[NobjTow]
   Float_t         TowOE[10000];   //[NobjTow]
   Int_t           Tow_jetidx[10000];   //[NobjTow]
   Int_t           NobjTrack;
   Int_t           TrackTowId[10000];   //[NobjTrack]
   Int_t           TrackTowIdPhi[10000];   //[NobjTrack]
   Int_t           TrackTowIdEta[10000];   //[NobjTrack]
   Int_t           TrackId[10000];   //[NobjTrack]
   Int_t           TrackNHits[10000];   //[NobjTrack]
   Bool_t          TrackQualityL[10000];   //[NobjTrack]
   Bool_t          TrackQualityT[10000];   //[NobjTrack]
   Bool_t          TrackQualityHP[10000];   //[NobjTrack]
   Float_t         TrackChi2[10000];   //[NobjTrack]
   Float_t         TrackPt[10000];   //[NobjTrack]
   Float_t         TrackEta[10000];   //[NobjTrack]
   Float_t         TrackPhi[10000];   //[NobjTrack]
   Float_t         TrackP[10000];   //[NobjTrack]
   Float_t         TrackDR[10000];   //[NobjTrack]
   Float_t         TrackPhiOut[10000];   //[NobjTrack]
   Float_t         TrackEtaOut[10000];   //[NobjTrack]
   Float_t         TrackDROut[10000];   //[NobjTrack]
   Float_t         TrackEMC1[10000];   //[NobjTrack]
   Float_t         TrackEMC3[10000];   //[NobjTrack]
   Float_t         TrackEMC5[10000];   //[NobjTrack]
   Float_t         TrackHAC1[10000];   //[NobjTrack]
   Float_t         TrackHAC3[10000];   //[NobjTrack]
   Float_t         TrackHAC5[10000];   //[NobjTrack]
   Int_t           Track_jetidx[10000];   //[NobjTrack]
   Float_t         MuDR[10000];   //[NobjTrack]
   Float_t         MuDE[10000];   //[NobjTrack]
   Int_t           NobjJet;
   Float_t         JetPt[100];   //[NobjJet]
   Float_t         JetPhi[100];   //[NobjJet]
   Float_t         JetEta[100];   //[NobjJet]
   Float_t         JetEt[100];   //[NobjJet]
   Float_t         JetE[100];   //[NobjJet]
   Float_t         JetCorrZSP[100];   //[NobjJet]
   Float_t         JetCorrL2[100];   //[NobjJet]
   Float_t         JetCorrL3[100];   //[NobjJet]
   Float_t         JetCorrJPT[100];   //[NobjJet]
   Float_t         JetCorrL2L3[100];   //[NobjJet]
   Float_t         JetCorrL2L3JPT[100];   //[NobjJet]
   Float_t         GenJetPt[100];   //[NobjJet]
   Float_t         GenJetPhi[100];   //[NobjJet]
   Float_t         GenJetEta[100];   //[NobjJet]
   Float_t         GenJetEt[100];   //[NobjJet]
   Float_t         GenJetE[100];   //[NobjJet]
   Int_t           NobjGenJet;
   Float_t         GenJetColPt[100];   //[NobjJet]
   Float_t         GenJetColPhi[100];   //[NobjJet]
   Float_t         GenJetColEta[100];   //[NobjJet]
   Float_t         GenJetColEt[100];   //[NobjJet]
   Float_t         GenJetColE[100];   //[NobjJet]
   Int_t           GenJetColJetIdx[100];
   Float_t         GenEvtScale;
   Float_t         Met;
   Float_t         MetPhi;
   Float_t         MetSum;
   Float_t         Weight;

   // List of branches
   TBranch        *b_NobjTow;   //!
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
   TBranch        *b_Tow_jetidx;   //!
   TBranch        *b_NobjTrack;   //!
   TBranch        *b_TrackTowId;   //!
   TBranch        *b_TrackTowIdPhi;   //!
   TBranch        *b_TrackTowIdEta;   //!
   TBranch        *b_TrackId;   //!
   TBranch        *b_TrackNHits;   //!
   TBranch        *b_TrackQualityL;   //!
   TBranch        *b_TrackQualityT;   //!
   TBranch        *b_TrackQualityHP;   //!
   TBranch        *b_TrackChi2;   //!
   TBranch        *b_TrackPt;   //!
   TBranch        *b_TrackEta;   //!
   TBranch        *b_TrackPhi;   //!
   TBranch        *b_TrackP;   //!
   TBranch        *b_TrackDR;   //!
   TBranch        *b_TrackPhiOut;   //!
   TBranch        *b_TrackEtaOut;   //!
   TBranch        *b_TrackDROut;   //!
   TBranch        *b_TrackEMC1;   //!
   TBranch        *b_TrackEMC3;   //!
   TBranch        *b_TrackEMC5;   //!
   TBranch        *b_TrackHAC1;   //!
   TBranch        *b_TrackHAC3;   //!
   TBranch        *b_TrackHAC5;   //!
   TBranch        *b_Track_jetidx;   //!
   TBranch        *b_MuDR;   //!
   TBranch        *b_MuDE;   //!
   TBranch        *b_NobjJet;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetEt;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_JetCorrZSP;   //!
   TBranch        *b_JetCorrL2;   //!
   TBranch        *b_JetCorrL3;   //!
   TBranch        *b_JetCorrJPT;   //!
   TBranch        *b_JetCorrL2L3;   //!
   TBranch        *b_JetCorrL2L3JPT;   //!
   TBranch        *b_GenJetPt;   //!
   TBranch        *b_GenJetPhi;   //!
   TBranch        *b_GenJetEta;   //!
   TBranch        *b_GenJetEt;   //!
   TBranch        *b_GenJetE;   //!
   TBranch        *b_NobjGenJet;
   TBranch        *b_GenJetColPt;   //!
   TBranch        *b_GenJetColPhi;   //!
   TBranch        *b_GenJetColEta;   //!
   TBranch        *b_GenJetColEt;   //!
   TBranch        *b_GenJetColE;   //!
   TBranch        *b_GenJetColJetIdx; //!
   TBranch        *b_GenEvtScale; //!
   TBranch        *b_Met;   //!
   TBranch        *b_MetPhi;   //!
   TBranch        *b_MetSum;   //!
   TBranch        *b_Weight;   //!

   NJetSel(TTree * /*tree*/ =0) { }
   virtual ~NJetSel() { }
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

   //ClassDef(NJetSel,0);
};

#endif

#ifdef NJetSel_cxx
void NJetSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NobjTow", &NobjTow, &b_NobjTow);
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
   fChain->SetBranchAddress("Tow_jetidx", Tow_jetidx, &b_Tow_jetidx);
   fChain->SetBranchAddress("NobjTrack", &NobjTrack, &b_NobjTrack);
   fChain->SetBranchAddress("TrackTowId", TrackTowId, &b_TrackTowId);
   fChain->SetBranchAddress("TrackTowIdPhi", TrackTowIdPhi, &b_TrackTowIdPhi);
   fChain->SetBranchAddress("TrackTowIdEta", TrackTowIdEta, &b_TrackTowIdEta);
   fChain->SetBranchAddress("TrackId", TrackId, &b_TrackId);
   fChain->SetBranchAddress("TrackNHits", TrackNHits, &b_TrackNHits);
   fChain->SetBranchAddress("TrackQualityL", TrackQualityL, &b_TrackQualityL);
   fChain->SetBranchAddress("TrackQualityT", TrackQualityT, &b_TrackQualityT);
   fChain->SetBranchAddress("TrackQualityHP", TrackQualityHP, &b_TrackQualityHP);
   fChain->SetBranchAddress("TrackChi2", TrackChi2, &b_TrackChi2);
   fChain->SetBranchAddress("TrackPt", TrackPt, &b_TrackPt);
   fChain->SetBranchAddress("TrackEta", TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackPhi", TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackP", TrackP, &b_TrackP);
   fChain->SetBranchAddress("TrackDR", TrackDR, &b_TrackDR);
   fChain->SetBranchAddress("TrackPhiOut", TrackPhiOut, &b_TrackPhiOut);
   fChain->SetBranchAddress("TrackEtaOut", TrackEtaOut, &b_TrackEtaOut);
   fChain->SetBranchAddress("TrackDROut", TrackDROut, &b_TrackDROut);
   fChain->SetBranchAddress("TrackEMC1", TrackEMC1, &b_TrackEMC1);
   fChain->SetBranchAddress("TrackEMC3", TrackEMC3, &b_TrackEMC3);
   fChain->SetBranchAddress("TrackEMC5", TrackEMC5, &b_TrackEMC5);
   fChain->SetBranchAddress("TrackHAC1", TrackHAC1, &b_TrackHAC1);
   fChain->SetBranchAddress("TrackHAC3", TrackHAC3, &b_TrackHAC3);
   fChain->SetBranchAddress("TrackHAC5", TrackHAC5, &b_TrackHAC5);
   fChain->SetBranchAddress("Track_jetidx", Track_jetidx, &b_Track_jetidx);
   fChain->SetBranchAddress("MuDR", MuDR, &b_MuDR);
   fChain->SetBranchAddress("MuDE", MuDE, &b_MuDE);
   fChain->SetBranchAddress("NobjJet", &NobjJet, &b_NobjJet);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetEt", JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("JetCorrZSP", JetCorrZSP, &b_JetCorrZSP);
   fChain->SetBranchAddress("JetCorrL2", JetCorrL2, &b_JetCorrL2);
   fChain->SetBranchAddress("JetCorrL3", JetCorrL3, &b_JetCorrL3);
   fChain->SetBranchAddress("JetCorrJPT", JetCorrJPT, &b_JetCorrJPT);
   fChain->SetBranchAddress("JetCorrL2L3", JetCorrL2L3, &b_JetCorrL2L3);
   fChain->SetBranchAddress("JetCorrL2L3JPT", JetCorrL2L3JPT, &b_JetCorrL2L3JPT);
   fChain->SetBranchAddress("GenJetPt", GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetPhi", GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetEta", GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetEt", GenJetEt, &b_GenJetEt);
   fChain->SetBranchAddress("GenJetE", GenJetE, &b_GenJetE);
   fChain->SetBranchAddress("NobjGenJet", &NobjGenJet, &b_NobjGenJet);
   fChain->SetBranchAddress("GenJetColPt", GenJetColPt, &b_GenJetColPt);
   fChain->SetBranchAddress("GenJetColPhi", GenJetColPhi, &b_GenJetColPhi);
   fChain->SetBranchAddress("GenJetColEta", GenJetColEta, &b_GenJetColEta);
   fChain->SetBranchAddress("GenJetColEt", GenJetColEt, &b_GenJetColEt);
   fChain->SetBranchAddress("GenJetColE", GenJetColE, &b_GenJetColE);
   fChain->SetBranchAddress("GenJetColJetIdx", GenJetColJetIdx, &b_GenJetColJetIdx);
   fChain->SetBranchAddress("GenEvtScale", &GenEvtScale, &b_GenEvtScale);
   fChain->SetBranchAddress("Met", &Met, &b_Met);
   fChain->SetBranchAddress("MetPhi", &MetPhi, &b_MetPhi);
   fChain->SetBranchAddress("MetSum", &MetSum, &b_MetSum);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
}

Bool_t NJetSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NJetSel_cxx
