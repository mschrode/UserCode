//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 23 22:39:24 100007 by ROOT version 5.14/00b
// from TTree TrackTowerTree/
// found on file: calib.root
//////////////////////////////////////////////////////////

#ifndef TrackTowerSel_h
#define TrackTowerSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class TrackTowerSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leave types
   Int_t           NobjTowCal;
   Int_t           TowId[10000];   //[NobjTowCal]
   Int_t           TowId_phi[10000];   //[NobjTowCal]
   Int_t           TowId_eta[10000];   //[NobjTowCal]
   Float_t         TowEt[10000];   //[NobjTowCal]
   Float_t         TowEta[10000];   //[NobjTowCal]
   Float_t         TowPhi[10000];   //[NobjTowCal]
   Float_t         TowE[10000];   //[NobjTowCal]
   Float_t         TowEm[10000];   //[NobjTowCal]
   Float_t         TowHad[10000];   //[NobjTowCal]
   Float_t         TowOE[10000];   //[NobjTowCal]
   Int_t           NobjTrackCal;
   Float_t         TrackEt[10000];   //[NobjTrackCal]
   Float_t         TrackEterr[10000];   //[NobjTrackCal]
   Float_t         TrackEta[10000];   //[NobjTrackCal]
   Float_t         TrackPhi[10000];   //[NobjTrackCal]
   Float_t         TrackE[10000];   //[NobjTrackCal]
   Int_t           ProcessID;
   Float_t         EventWeight;

   // List of branches
   TBranch        *b_NobjTowCal;   //!
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
   TBranch        *b_NobjTrackCal;   //!
   TBranch        *b_TrackEt;   //!
   TBranch        *b_TrackEterr;   //!
   TBranch        *b_TrackEta;   //!
   TBranch        *b_TrackPhi;   //!
   TBranch        *b_TrackE;   //!
   TBranch        *b_ProcessID;
   TBranch        *b_EventWeight;

   TrackTowerSel(TTree * /*tree*/ =0) { }
   virtual ~TrackTowerSel() { }
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

   //ClassDef(TrackTowerSel,0);
};

#endif

#ifdef TrackTowerSel_cxx
void TrackTowerSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NobjTowCal", &NobjTowCal, &b_NobjTowCal);
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
   fChain->SetBranchAddress("NobjTrackCal", &NobjTrackCal, &b_NobjTrackCal);
   fChain->SetBranchAddress("TrackEt", TrackEt, &b_TrackEt);
   fChain->SetBranchAddress("TrackEterr", TrackEterr, &b_TrackEterr);
   fChain->SetBranchAddress("TrackEta", TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackPhi", TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackE", TrackE, &b_TrackE);
   fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
}

Bool_t TrackTowerSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TrackTowerSel_cxx
