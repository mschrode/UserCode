// $ Id: $

#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TPad.h"
#include "TString.h"
#include "TTree.h"

#include "BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


unsigned int COUNT = 0;

class EvtInfo {
public:
  double pt1_;
  double pt2_;
  double r1_;
  double r2_;
  int run_;
  int lumi_;
  int evt_;
  int nPU_;
};


void plotResponse(const TString &fileName, double maxResp, TH1* &hRespPreSel,  TH1* &hRespFullSel, TH1* &hNPU, TH1* &hNPUAll, std::vector<EvtInfo> &info) {

  std::cout << "Processing file '" << fileName << "'" << std::endl;
  hRespPreSel = new TH1D("hRespPreSel"+util::toTString(COUNT),";Response",101,0.,15.);
  hRespPreSel->SetLineWidth(2);
  hRespFullSel = static_cast<TH1*>(hRespPreSel->Clone("hRespFullSel"+util::toTString(COUNT)));
  hRespFullSel->SetLineColor(kRed);
  hNPU = new TH1D("hNPU"+util::toTString(COUNT),";N(PU Vertices)",30,0.,30.);
  hNPU->SetLineWidth(2);
  hNPUAll = static_cast<TH1*>(hNPU->Clone("hNPUAll"+util::toTString(COUNT)));
  hNPU->SetLineColor(kBlue);


  TFile file(fileName,"READ");
  TTree* djt = 0;
  file.GetObject("DiJetTree",djt);
  if( djt == 0 ) {
    std::cerr << "  ERROR reading tree from file" << std::endl;
    exit(1);
  }

  const int maxNJet = 50;
  Int_t NObjJet = 0;
  Int_t Idx[maxNJet];
  Float_t JetPt[maxNJet];
  Float_t GenJetPt[maxNJet];
  Float_t JetCorrL1[maxNJet];
  Float_t JetCorrL2L3[maxNJet];
  UInt_t RunNumber = 0;
  UInt_t LumiBlockNumber = 0;
  UInt_t EventNumber = 0;
  Int_t PUMCNumVtx = 0;
  djt->SetBranchAddress("NobjJet",&NObjJet);
  djt->SetBranchAddress("RunNumber", &RunNumber);
  djt->SetBranchAddress("LumiBlockNumber", &LumiBlockNumber);
  djt->SetBranchAddress("EventNumber", &EventNumber);
  djt->SetBranchAddress("PUMCNumVtx",&PUMCNumVtx);
  djt->SetBranchAddress("JetPt",JetPt);
  djt->SetBranchAddress("GenJetPt", GenJetPt);
  djt->SetBranchAddress("JetCorrL1",JetCorrL1);
  djt->SetBranchAddress("JetCorrL2L3",JetCorrL2L3);
  djt->SetBranchAddress("L2L3CorrJetColJetIdx",Idx);

  for(int i = 0; i < djt->GetEntries(); ++i) {
    djt->GetEntry(i);
    
    double p1 = JetCorrL1[Idx[0]]*JetCorrL2L3[Idx[0]]*JetPt[Idx[0]];
    double p2 = JetCorrL1[Idx[1]]*JetCorrL2L3[Idx[1]]*JetPt[Idx[1]];
    double pg1 = GenJetPt[Idx[0]];
    double pg2 = GenJetPt[Idx[1]];
    double r1 = 0.;
    double r2 = 0.;
    if( pg1 && pg2 ) {
      r1 = p1/pg1;
      r2 = p2/pg2;
      hRespPreSel->Fill(r1);
      hRespPreSel->Fill(r2);
      if( !( NObjJet > 2 && JetCorrL1[Idx[2]]*JetCorrL2L3[Idx[2]]*JetPt[Idx[2]] > 0.3*0.5*(p1+p2) ) ) {
	hRespFullSel->Fill(r1);
	hRespFullSel->Fill(r2);
      }	
      hNPUAll->Fill(PUMCNumVtx);
      if( r1 > maxResp || r2 > maxResp ) {
	EvtInfo inf;
	inf.pt1_ = p1;
	inf.pt2_ = p2;
	inf.r1_ = r1;
	inf.r2_ = r2;
	inf.run_ = RunNumber;
	inf.lumi_ = LumiBlockNumber;
	inf.evt_ = EventNumber;
	inf.nPU_ = PUMCNumVtx;
	info.push_back(inf);
	hNPU->Fill(PUMCNumVtx);
      }
    }
  }
  file.Close();

  if( hNPU->Integral() ) hNPU->Scale(1./hNPU->Integral("width"));
  if( hNPUAll->Integral() ) hNPUAll->Scale(1./hNPUAll->Integral("width"));

  ++COUNT;
}



void investigateMCResponse(const TString &binAdmCfg, double maxResp = 10.) {
  util::StyleSettings::presentationNoTitle();
  gErrorIgnoreLevel = 1001;
  
  sampleTools::BinningAdmin binAdmin(binAdmCfg);
  TString baseName = "~/lustre/Analysis2011/KalibriDiJetSkims/KalibriSkim_PF_L1FastJet_MCSummer11_Eta0";
  unsigned int etaBin = 0;
  for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
    //for(unsigned int ptBin = 0; ptBin < 1; ++ptBin) {
    TString fileName = baseName+"_Pt"+util::toTString(ptBin)+".root";
    std::vector<EvtInfo> info;
    TH1* hRespPreSel = 0;
    TH1* hRespFullSel = 0;
    TH1* hNPU = 0;
    TH1* hNPUAll = 0;
    plotResponse(fileName,maxResp,hRespPreSel,hRespFullSel,hNPU,hNPUAll,info);
    util::HistOps::setAxisTitles(hRespPreSel,"Response","","events");
    util::HistOps::setAxisTitles(hRespFullSel,"Response","","events");
    util::HistOps::setAxisTitles(hNPU,"N(PU Vertices)","","events",true);
    util::HistOps::setAxisTitles(hNPUAll,"N(PU Vertices)","","events",true);
    util::HistOps::setYRange(hRespPreSel,40,0.3);
    util::HistOps::setYRange(hNPU,4);
    util::HistOps::setYRange(hNPUAll,4);

    TPaveText* label = util::LabelFactory::createPaveText(2);
    label->AddText("Summer11, "+util::LabelFactory::labelJetAlgo(baseName));
    label->AddText(util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", "+util::toTString(binAdmin.ptMin(etaBin,ptBin))+" < p^{ave}_{T} < "+util::toTString(binAdmin.ptMax(etaBin,ptBin))+" GeV");
    TLegend* legRes = util::LabelFactory::createLegendWithOffset(2,2);
    legRes->AddEntry(hRespPreSel,"Dijet Preselection","L");
    legRes->AddEntry(hRespFullSel,"Dijet Preselection + p_{T,3} < 0.3 #upoint p^{ave}_{T}","L");
    TLegend* legPU = util::LabelFactory::createLegendWithOffset(2,2);
    legPU->AddEntry(hNPUAll,"All : N = "+util::toTString(hNPUAll->GetEntries(),0)+", #LTN(PU)#GT = "+util::toTString(hNPUAll->GetMean(),1),"L");
    legPU->AddEntry(hNPU,"R > "+util::toTString(maxResp)+" : N = "+util::toTString(hNPU->GetEntries(),0)+", #LTN(PU)#GT = "+util::toTString(hNPU->GetMean(),1),"L");

    TCanvas* cResp = new TCanvas("cResp"+util::toTString(COUNT),"Resp PtBin"+util::toTString(ptBin),500,500);
    cResp->cd();
    hRespPreSel->Draw();
    hRespFullSel->Draw("same");
    label->Draw("same");
    legRes->Draw("same");
    cResp->SetLogy();
    gPad->RedrawAxis();
    cResp->SaveAs("QCDFlat_Summer11_DijetResponse_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+".eps","eps");

    TCanvas* cNPU = new TCanvas("cNPU"+util::toTString(COUNT),"N(PU) PtBin"+util::toTString(ptBin),500,500);
    cNPU->cd();
    hNPU->Draw("HIST");
    hNPUAll->Draw("HISTsame");
    label->Draw("same");
    legPU->Draw("same");
    cNPU->SaveAs("QCDFlat_Summer11_DijetResponseNPUDependence_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+".eps","eps");
  }
}
