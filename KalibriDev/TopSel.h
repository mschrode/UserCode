//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 13 17:16:28 2008 by ROOT version 5.18/00a
// from TTree TopTree/
// found on file: Top_Calib.root
//////////////////////////////////////////////////////////

#ifndef TopSel_h
#define TopSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class TopSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           NobjTow;
   Int_t           TowId[1000];   //[NobjBTow]
   Int_t           TowId_phi[1000];   //[NobjBTow]
   Int_t           TowId_eta[1000];   //[NobjBTow]
   Float_t         TowEt[1000];   //[NobjBTow]
   Float_t         TowEta[1000];   //[NobjBTow]
   Float_t         TowPhi[1000];   //[NobjBTow]
   Float_t         TowE[1000];   //[NobjBTow]
   Float_t         TowEm[1000];   //[NobjBTow]
   Float_t         TowHad[1000];   //[NobjBTow]
   Float_t         TowOE[1000];   //[NobjBTow]
   Int_t           Tow_jetidx[1000];   //[NobjBTow]
   Int_t           NobjJet;
   Float_t         JetPt[8];   //[NobjBJet]
   Float_t         JetPhi[8];   //[NobjBJet]
   Float_t         JetEta[8];   //[NobjBJet]
   Float_t         JetEt[8];   //[NobjBJet]
   Float_t         JetE[8];   //[NobjBJet]
   Float_t         GenJetPt[8];   //[NobjJet]
   Float_t         GenJetPhi[8];   //[NobjJet]
   Float_t         GenJetEta[8];   //[NobjJet]
   Float_t         GenJetEt[8];   //[NobjJet]
   Float_t         GenJetE[8];   //[NobjJet]
   Int_t           JetFlavor[8];   //[NobjBJet]
   Int_t           JetTopID[8];   //[NobjBJet]
   Float_t         JetCorrL1[8];   //[NobjJet]
   Float_t         JetCorrL2[8];   //[NobjJet]
   Float_t         JetCorrL3[8];   //[NobjJet]
   Float_t         JetCorrL4[8];   //[NobjJet]
   Float_t         JetCorrL5[8];   //[NobjJet]
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
   TBranch        *b_NobjJet;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetEt;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_GenJetPt;   //!
   TBranch        *b_GenJetPhi;   //!
   TBranch        *b_GenJetEta;   //!
   TBranch        *b_GenJetEt;   //!
   TBranch        *b_GenJetE;   //!
   TBranch        *b_JetFlavor;   //!
   TBranch        *b_JetTopID;   //!
   TBranch        *b_JetCorrL1;   //!
   TBranch        *b_JetCorrL2;   //!
   TBranch        *b_JetCorrL3;   //!
   TBranch        *b_JetCorrL4;   //!
   TBranch        *b_JetCorrL5;   //!
   TBranch        *b_Weight;   //!

   TopSel(TTree * /*tree*/ =0) { }
   virtual ~TopSel() { }
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

   //ClassDef(TopSel,0);
};

#endif

#ifdef TopSel_cxx
void TopSel::Init(TTree *tree)
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
   fChain->SetBranchAddress("NobjJet", &NobjJet, &b_NobjJet);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetEt", JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("GenJetPt", GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetPhi", GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetEta", GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetEt", GenJetEt, &b_GenJetEt);
   fChain->SetBranchAddress("GenJetE", GenJetE, &b_GenJetE);
   fChain->SetBranchAddress("JetFlavor", JetFlavor, &b_JetFlavor);
   fChain->SetBranchAddress("JetTopID", JetTopID, &b_JetTopID);
   fChain->SetBranchAddress("JetCorrL1", JetCorrL1, &b_JetCorrL1);
   fChain->SetBranchAddress("JetCorrL2", JetCorrL2, &b_JetCorrL2);
   fChain->SetBranchAddress("JetCorrL3", JetCorrL3, &b_JetCorrL3);
   fChain->SetBranchAddress("JetCorrL4", JetCorrL4, &b_JetCorrL4);
   fChain->SetBranchAddress("JetCorrL5", JetCorrL5, &b_JetCorrL5);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
}

Bool_t TopSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TopSel_cxx
