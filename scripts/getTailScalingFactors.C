// $Id: getTailScalingFactors.C,v 1.12 2011/01/03 18:20:23 mschrode Exp $

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"

#include "globalFunctions.h"
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



const double gBINLABEL_WIDTH = -0.48;
const double gLEG_WIDTH = 0.48;


///////////////////////// TYPE DEFINITIONS /////////////////////////////////////////////////


void setStyleDataMarker(TH1* h);
void setStyleMCFilled(TH1* h);
void setStyleData(TGraphAsymmErrors* g);
void setStyleMC(TGraphAsymmErrors* g);


// ------------------------------------------------------------------------------------
class Pt3Bin {
public:
  Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, double tailWindowDataMin, double tailWindowDataMax, double tailWindowMCMin, double tailWindowMCMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel);
  Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel);
  ~Pt3Bin();

  double pt3Thres() const { return pt3Thres_; }

  double coreScalingFactor() const { return coreScalingFactor_; }
  double fTailData() const { return fNData_; }
  double fTailDataErr() const { return fNDataErr_; }
  double fTailMC() const { return fNMC_; }
  double fTailMCErr() const { return fNMCErr_; }
  double fTailMCSmeared() const { return fNMCSmeared_; }
  double fTailMCSmearedErr() const { return fNMCSmearedErr_; }
  double fTailMCSmearedNonGauss() const { return fNMCSmearedNonGauss_; }
  double fTailMCSmearedNonGaussErr() const { return fNMCSmearedNonGaussErr_; }

  void plotAsymmetryDataMC(const TString &outNameIdw) const;
  void plotAsymmetryDataMCSmeared(const TString &outNameId) const;
  void plotMCTruthResponse(const TString &outNameId) const;


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const unsigned int pt3Bin_;
  const double pt3Thres_;

  double coreScalingFactor_;
  int tailStartBinData_;
  int tailEndBinData_;
  int tailStartBinMC_;
  int tailEndBinMC_;
  double fNMCSmeared_;
  double fNMCSmearedErr_;
  double fNMCSmearedNonGauss_;
  double fNMCSmearedNonGaussErr_;
  double fNMC_;
  double fNMCErr_;
  double fNData_;
  double fNDataErr_;

  TH1* hAsymData_;
  TH1* hAsymMC_;
  TH1* hAsymMCSmeared_;
  TH1* hResp_;
  TF1* fGaussMCTruth_;

  TPaveText* binLabel_;

  TH1* readHist(const TString &fileName, const TString &id) const;  
  TH1* readMCTruthResponse(const TString &fileName, const TString &id) const;
  void init(double nSigCore, double coreScalingFactor, double tailWindowDataMin, double tailWindowDataMax, double tailWindowMCMin, double tailWindowMCMax);
  void initBinLabel(const sampleTools::BinningAdmin* adm, const TString &jetLabel, bool isMCTruth = false);
  void getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr) const;
};



// ------------------------------------------------------------------------------------
class EtaPtBin {
public:
  EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor,const sampleTools::BinningAdmin* adm, const TString &jetLabel);
  ~EtaPtBin();


  unsigned int nPt3Bins() const { return pt3Bins_.size(); }
  unsigned int etaBin() const { return etaBin_; }
  unsigned int ptBin() const { return ptBin_; }

  double extraMC() const { return extraMC_; }
  double extraMCErr() const { return extraMCErr_; }
  double extraData() const { return extraData_; }
  double extraDataErr() const { return extraDataErr_; }
  double deltaExtra() const { return deltaEx_; }
  double deltaExtraErr() const { return deltaExErr_; }
  double scalingFactor() const { return scalingFactor_; }
  double scalingFactorErr() const { return scalingFactorErr_; }
  
  void plotAsymmetryDistributions(const TString &outNameId) const;
  void plotExtrapolation(const TString &outNameId) const;
  void plotMCTruthResponse(const TString &outNameId) const;

  bool hasMCTruthResponse() const { return hasMCTruthResponse_; }

  void addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC);
  void addMCTruthResponse(const TString &fileName);
  void findWindow(const TString &fileNameData, const TString &fileNameMC, double nSigTailWindowMin, double nSigTailWindowMax);
  void extrapolate(double minPt3Data, bool fixDataShape, bool mcTruthRef);


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const double nSigCore_;
  const double coreScalingFactor_;
  const double exMin_;
  const double exMax_;
  const sampleTools::BinningAdmin* binAdm_;
  const TString jetLabel_;

  double tailWindowDataMin_;
  double tailWindowDataMax_;
  double tailWindowMCMin_;
  double tailWindowMCMax_;
  TGraphAsymmErrors* gFTailMC_;
  TGraphAsymmErrors* gFTailData_;
  TGraphAsymmErrors* gFTailMCTruth_;
  TGraphAsymmErrors* gFTailMCTruthNonGauss_;
  TGraphAsymmErrors* gFTailSpreadData_;
  TGraphAsymmErrors* gFTailSpreadMC_;
  TF1* fExMC_;
  TF1* fExData_;
  double extraMC_;
  double extraMCErr_;
  double extraData_;
  double extraDataErr_;
  double deltaEx_;
  double deltaExErr_;
  double scalingFactor_;
  double scalingFactorErr_;

  std::vector<Pt3Bin*> pt3Bins_;
  Pt3Bin* mcTruthResponse_;
  bool hasMCTruthResponse_;

  TPaveText* binLabel_;

  TString binId() const { return "_EtaBin"+util::toTString(etaBin_)+"_PtBin"+util::toTString(ptBin_); }
  TString binId(unsigned int pt3BinIdx) const { return binId()+"_Pt3Bin"+util::toTString(pt3BinIdx)+"_"; }
};





//////////////////////////////// MAIN ROUTINE ///////////////////////////////////////////


// ------------------------------------------------------------------------------------
void getTailScalingFactors() {

  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting up parameters" << std::endl;

  const TString outLabel = "Tail";
  const double nSigCore = 2.;

  const double nSigTailStart = 3.5;
  const double nSigTailEnd = 6.5;

  const double minPt3Data = 0.;
  const bool fixDataShape = false;
  const bool mcTruthRef = true;

  const TString jetAlgo = "PF";

  util::StyleSettings::paperNoTitle();
  gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("../sampleTools/BinningAdminTailsRebinned.cfg");

  TString id = outLabel+"_"+jetAlgo;
  TString jetLabel = "Anti-k_{T} (d=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";




  // +++++ Bins +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating bins" << std::endl;

  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    //etaBin = 0;
    double coreScaling = 0.;
    if( etaBin == 0 )      coreScaling = 0.041;
    else if( etaBin == 1 ) coreScaling = 0.020;
    else if( etaBin == 2 ) coreScaling = 0.;
    else if( etaBin == 3 ) coreScaling = 0.041;

    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {

      TString fileNameData = "~/results/ResolutionFit/TailScalingNote/Tails_PF_DataRebinned_Eta";
      fileNameData += etaBin;
      TString fileNameMC = "~/results/ResolutionFit/TailScalingNote/Tails_PF_MCRebinned_Eta";
      fileNameMC += etaBin;

      // Create eta-pt bin
      EtaPtBin* bin = new EtaPtBin(etaBin,ptBin,nSigCore,coreScaling,binAdm,jetLabel);

      // Define window
      bin->findWindow(fileNameMC+"_PtSoft0.root",fileNameMC+"_PtSoft0.root",nSigTailStart,nSigTailEnd);

      // Add pt3 bins
      bin->addPt3Bin(0.050,fileNameData+"_PtSoft0.root",fileNameMC+"_PtSoft0.root");
      bin->addPt3Bin(0.075,fileNameData+"_PtSoft1.root",fileNameMC+"_PtSoft1.root");
      bin->addPt3Bin(0.100,fileNameData+"_PtSoft2.root",fileNameMC+"_PtSoft2.root");
      bin->addPt3Bin(0.125,fileNameData+"_PtSoft3.root",fileNameMC+"_PtSoft3.root");
      bin->addPt3Bin(0.150,fileNameData+"_PtSoft4.root",fileNameMC+"_PtSoft4.root");
      bin->addPt3Bin(0.175,fileNameData+"_PtSoft5.root",fileNameMC+"_PtSoft5.root");
      bin->addPt3Bin(0.200,fileNameData+"_PtSoft6.root",fileNameMC+"_PtSoft6.root");
      bin->addPt3Bin(0.225,fileNameData+"_PtSoft7.root",fileNameMC+"_PtSoft7.root");
      bin->addPt3Bin(0.250,fileNameData+"_PtSoft8.root",fileNameMC+"_PtSoft8.root");
      bin->addPt3Bin(0.275,fileNameData+"_PtSoft9.root",fileNameMC+"_PtSoft9.root");
      bin->addPt3Bin(0.300,fileNameData+"_PtSoft10.root",fileNameMC+"_PtSoft10.root");

      // Add mc truth response
      bin->addMCTruthResponse("results/SymmetrizedMCTruth_Rebinned_Eta"+util::toTString(etaBin)+".root");

      etaPtBins.push_back(bin);
    }
  }



  // +++++ Asymmetry and extrapolation plots ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating asymmetry and extrapolation plots" << std::endl;

  for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
      it != etaPtBins.end(); ++it ) {

    // Asymmetry
    (*it)->plotAsymmetryDistributions(id);
    (*it)->plotMCTruthResponse(id);
    
    // Extrapolation
    (*it)->extrapolate(minPt3Data,fixDataShape,mcTruthRef);  
    (*it)->plotExtrapolation(id);

    // Print fractional number of events  
//     std::cout << std::setprecision(2) << " & $" << (*it)->extraData() << " \\pm " << (*it)->extraDataErr() << "$" << std::flush;
//     std::cout << " & $" << (*it)->extraMC() << " \\pm " << (*it)->extraMCErr() << "$ " << std::endl;
  }



  // +++++ Scaling factors ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating scaling factor plots" << std::endl;

  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    TH1* hDelta = new TH1D("hDelta_Eta"+util::toTString(etaBin),
			   ";p^{ave}_{T} (GeV);#Delta = f^{Data}_{Asym}(0) - f^{MC}_{Asym}(0)",
			   binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hDelta->SetMarkerStyle(20);

    TH1* hScale = new TH1D("hScale_Eta"+util::toTString(etaBin),
			   ";p^{ave}_{T} (GeV);Scaling Factor",
			   binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hScale->SetMarkerStyle(20);

    for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
	it != etaPtBins.end(); ++it) {
      if( (*it)->etaBin() == etaBin ) {
	int bin = 1+(*it)->ptBin();
	hDelta->SetBinContent(bin,(*it)->deltaExtra());
	hDelta->SetBinError(bin,(*it)->deltaExtraErr());
	hScale->SetBinContent(bin,(*it)->scalingFactor());
	hScale->SetBinError(bin,(*it)->scalingFactorErr());
      }
    }


    // Plots
    TCanvas* can = new TCanvas("can","",500,500);

    TH1* hDeltaFrame = static_cast<TH1D*>(hDelta->Clone("hDeltaFrame_Eta"+util::toTString(etaBin)));
    hDeltaFrame->SetLineColor(kBlack);
    hDeltaFrame->SetLineStyle(2);
    hDeltaFrame->SetMarkerStyle(1);
    for(int i = 1; i <= hDeltaFrame->GetNbinsX(); ++i) {
      hDeltaFrame->SetBinContent(i,0.);
      hDeltaFrame->SetBinError(i,0.);
    }
    hDeltaFrame->GetYaxis()->SetRangeUser(-0.0045,0.013);
    hDeltaFrame->GetXaxis()->SetMoreLogLabels();
    can->cd();
    hDeltaFrame->Draw("HIST");
    hDelta->Draw("PE1same");
    can->SetLogx();
    can->SaveAs(id+"_EtaBin"+util::toTString(etaBin)+"_Delta.eps","eps");

    TH1* hScaleFrame = static_cast<TH1D*>(hScale->Clone("hScaleFrame_Eta"+util::toTString(etaBin)));
    hScaleFrame->SetLineColor(kBlack);
    hScaleFrame->SetLineStyle(2);
    hScaleFrame->SetMarkerStyle(1);
    for(int i = 1; i <= hScaleFrame->GetNbinsX(); ++i) {
      hScaleFrame->SetBinContent(i,1.);
      hScaleFrame->SetBinError(i,1.);
    }
    hScaleFrame->GetYaxis()->SetRangeUser(0.69,2.49);
    hScaleFrame->GetXaxis()->SetMoreLogLabels();
    can->cd();
    hScaleFrame->Draw("HIST");
    hScale->Draw("PE1same");
    can->SetLogx();
    can->SaveAs(id+"_EtaBin"+util::toTString(etaBin)+"_ScalingFactors.eps","eps");


    // To ROOT file
    TString fileMode = "UPDATE";
    if( etaBin == 0 ) fileMode == "RECREATE";
    TFile outFile(id+".root",fileMode);
    outFile.WriteTObject(hDelta);
    outFile.WriteTObject(hDeltaFrame);
    outFile.WriteTObject(hScale);    
    outFile.WriteTObject(hScaleFrame);   
    outFile.Close();
    

    delete hDelta;
    delete hScale;
    delete hDeltaFrame;
    delete hScaleFrame;
    delete can;
  }
}










//////////////////////////////// IMPLEMENTATIONS ///////////////////////////////////////


// ------------------------------------------------------------------------------------
EtaPtBin::EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, const sampleTools::BinningAdmin* binAdm, const TString &jetLabel)
  : etaBin_(etaBin), ptBin_(ptBin),
    nSigCore_(nSigCore), coreScalingFactor_(coreScalingFactor),
    exMin_(0.), exMax_(0.32),
    binAdm_(binAdm), jetLabel_(jetLabel) {

  tailWindowDataMin_ = 0.;
  tailWindowDataMax_ = 0.;
  tailWindowMCMin_ = 0.;
  tailWindowMCMax_ = 0.;

  gFTailData_ = new TGraphAsymmErrors(0);
  gFTailMC_ = new TGraphAsymmErrors(0);
  gFTailMCTruth_ = new TGraphAsymmErrors(0);
  gFTailSpreadMC_ = new TGraphAsymmErrors(0);
  gFTailSpreadData_  = new TGraphAsymmErrors(0);
  gFTailMCTruthNonGauss_ = new TGraphAsymmErrors(0);

  fExMC_ = new TF1("fExMC"+binId(),"[0] + sq([1])*x + [2]*sq(x) + [3]*x*x*x",exMin_,exMax_);
  fExMC_->SetParameter(0,0.005);
  fExMC_->SetParameter(1,0.);
  fExMC_->SetParameter(2,0.05);
  fExMC_->SetParameter(3,0.);
  fExMC_->SetLineWidth(1);

  fExData_ = static_cast<TF1*>(fExMC_->Clone("fExData"+binId()));
  fExData_->SetLineStyle(2);

  deltaEx_ = 0.;
  deltaExErr_ = 0.;
  scalingFactor_ = 1.;
  scalingFactorErr_ = 0.;

  mcTruthResponse_ = 0;
  hasMCTruthResponse_ = false;

  binLabel_ = util::LabelFactory::createPaveText(4,gBINLABEL_WIDTH);
  binLabel_->AddText(jetLabel_);
  binLabel_->AddText(util::toTString(binAdm_->etaMin(etaBin_))+" < |#eta| < "+util::toTString(binAdm_->etaMax(etaBin_)));
  binLabel_->AddText(util::toTString(binAdm_->ptMin(etaBin_,ptBin_))+" < p^{ave}_{T} < "+util::toTString(binAdm_->ptMax(etaBin_,ptBin_))+" GeV");
}


// ------------------------------------------------------------------------------------
EtaPtBin::~EtaPtBin() {
  for(std::vector<Pt3Bin*>::iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    delete *it;
  }
  delete gFTailMC_;
  delete gFTailData_;
  delete gFTailMCTruth_;
  delete gFTailMCTruthNonGauss_;
  delete gFTailSpreadData_;
  delete gFTailSpreadMC_;
  delete fExMC_;
  delete fExData_;
  if( hasMCTruthResponse() ) delete mcTruthResponse_;
  if( binLabel_ ) delete binLabel_;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotAsymmetryDistributions(const TString &outNameId) const {
  unsigned int pt3Bin = 0;
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it, ++pt3Bin) {
    (*it)->plotAsymmetryDataMC(outNameId+binId(pt3Bin));
    (*it)->plotAsymmetryDataMCSmeared(outNameId+binId(pt3Bin));
  }
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotMCTruthResponse(const TString &outNameId) const {
  if( hasMCTruthResponse() ) mcTruthResponse_->plotMCTruthResponse(outNameId+binId());
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotExtrapolation(const TString &outNameId) const {

  // Extrapolation MC Closure
  TLegend* leg = util::LabelFactory::createLegendCol(3,gLEG_WIDTH);
  leg->AddEntry(gFTailMC_,"Asymmetry","P");
  leg->AddEntry(fExMC_,"Extrapolation","L");
  if( hasMCTruthResponse() ) leg->AddEntry(gFTailMCTruth_,"Response","P");

  TCanvas* can = new TCanvas("can","Number of events",500,500);
  can->cd();
  TH1* hFrame = new TH1D("hFrame",";p^{rel}_{T,3} Threshold;f_{Asym}",1000,0.,0.44);
  hFrame->GetYaxis()->SetRangeUser(0.,2.3*gFTailMC_->GetY()[gFTailMC_->GetN()-1]);
  hFrame->Draw();
  fExMC_->Draw("same");
  if( hasMCTruthResponse() ) {
    gFTailMCTruth_->Draw("PE1same");
  }
  gFTailMC_->Draw("PE1same");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_ExtrapolationMCClosure.eps","eps");
  delete leg;

  // Extrapolation MC + Data
  leg = util::LabelFactory::createLegendCol(4,gLEG_WIDTH);
  leg->AddEntry(gFTailData_,"Asymmetry Data","P");
  leg->AddEntry(fExData_,"Extrapolation Data","L");
  leg->AddEntry(gFTailMC_,"Asymmetry MC","P");
  leg->AddEntry(fExMC_,"Extrapolation MC","L");

  can->cd();
  hFrame->Draw();
  fExMC_->Draw("same");
  fExData_->Draw("same");
  gFTailMC_->Draw("PE1same");
  gFTailData_->Draw("PE1same");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_Extrapolation.eps","eps");
  delete leg;

  // Extrapolation MC Closure + Data + Shifted Extrapolation
  leg = util::LabelFactory::createLegendCol(5,gLEG_WIDTH);
  leg->AddEntry(gFTailData_,"Asymmetry Data","P");
  leg->AddEntry(fExData_,"Extrapolation Data","L");
  leg->AddEntry(gFTailMC_,"Asymmetry MC","P");
  leg->AddEntry(fExMC_,"Extrapolation MC","L");
  if( hasMCTruthResponse() ) leg->AddEntry(gFTailMCTruth_,"Response","P");

  can->cd();
  hFrame->Draw();
  fExMC_->Draw("same");
  fExData_->Draw("same");
  gFTailMC_->Draw("PE1same");
  gFTailData_->Draw("PE1same");
  if( hasMCTruthResponse() ) {
    gFTailMCTruth_->Draw("PE1same");
    //gFTailMCTruthNonGauss_->Draw("PE1same");
  }
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_Extrapolation2.eps","eps");
  delete leg;

  // Spread of ftail for mc
  double relErr = extraMCErr_/extraMC();
  hFrame->GetYaxis()->SetTitle("( f_{Asym} - Fit ) / Fit");
  hFrame->GetYaxis()->SetRangeUser(std::min(-7.*relErr,0.),10.*relErr);
  for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
    hFrame->SetBinContent(i,0.);
  }
  hFrame->SetLineStyle(2);
  TH1* hBand = static_cast<TH1D*>(hFrame->Clone("hBand"));
  for(int i = 1; i <= hBand->GetNbinsX(); ++i) {
    hBand->SetBinError(i,relErr);
  }
  hBand->SetFillStyle(1001);
  hBand->SetFillColor(42);

  leg = util::LabelFactory::createLegendCol(2,gLEG_WIDTH);
  leg->AddEntry(gFTailSpreadMC_,"Asymmetry MC","P");
  leg->AddEntry(hBand,"Standard Dev.","F");

  can->cd();
  hBand->Draw("E3");
  hFrame->Draw("same");
  gFTailSpreadMC_->Draw("Psame");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_SpreadMC.eps","eps");  
  delete leg;

  // Spread of ftail for data
  relErr = extraDataErr_/extraData();
  hFrame->GetYaxis()->SetRangeUser(std::min(-7.*relErr,0.),10.*relErr);

  hBand->GetYaxis()->SetRangeUser(std::min(-7.*relErr,0.),10.*relErr);
  for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
    hFrame->SetBinContent(i,0.);
  }
  for(int i = 1; i <= hBand->GetNbinsX(); ++i) {
    hBand->SetBinError(i,relErr);
  }

  leg = util::LabelFactory::createLegendCol(2,gLEG_WIDTH);
  leg->AddEntry(gFTailSpreadData_,"Asymmetry Data","P");
  leg->AddEntry(hBand,"Standard Dev.","F");

  can->cd();
  hBand->Draw("E3");
  hFrame->Draw("same");
  gFTailSpreadData_->Draw("Psame");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_SpreadData.eps","eps");  
  delete leg;

  delete hFrame;
  delete hBand;
  delete can;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC) {
  pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bins_.size(),thres,nSigCore_,coreScalingFactor_,tailWindowDataMin_,tailWindowDataMax_,tailWindowMCMin_,tailWindowMCMax_,binAdm_,jetLabel_));
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addMCTruthResponse(const TString &fileName) {
  if( hasMCTruthResponse() ) delete mcTruthResponse_;
  double winMin = 1.+sqrt(2.)*tailWindowMCMin_;
  double winMax = 1.+sqrt(2.)*tailWindowMCMax_;
  mcTruthResponse_ = new Pt3Bin(fileName,etaBin_,ptBin_,nSigCore_,coreScalingFactor_,winMin,winMax,binAdm_,jetLabel_);
  hasMCTruthResponse_ = true;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::findWindow(const TString &fileNameData, const TString &fileNameMC, double nSigTailWindowMin, double nSigTailWindowMax) {
  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);

  TH1* h = util::FileOps::readTH1(fileNameMC,histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
    tailWindowMCMin_ = std::abs(nSigTailWindowMin*width);
    tailWindowMCMax_ = std::abs(nSigTailWindowMax*width);
    //std::cout << "SIG(" << etaBin_ << ", " << ptBin_ << "): " << width << " \\pm " << widthErr << std::endl;
  }
  delete h;

  h = util::FileOps::readTH1(fileNameData,histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  width = 0.;
  widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
    tailWindowDataMin_ = std::abs(nSigTailWindowMin*width);
    tailWindowDataMax_ = std::abs(nSigTailWindowMax*width);
  }
  delete h;

  binLabel_->AddText(util::toTString(nSigTailWindowMin)+" - "+util::toTString(nSigTailWindowMax)+" #sigma  Window");
}


// ------------------------------------------------------------------------------------
void EtaPtBin::extrapolate(double minPt3Data, bool fixDataShape, bool mcTruthRef) {

  // Fill graphs of ftail vs pt3 for mc
  std::vector<double> pt3;
  std::vector<double> pt3e;
  std::vector<double> n;
  std::vector<double> ne;
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    double thres = (*it)->pt3Thres();
    pt3.push_back(thres);
    pt3e.push_back(0.);
    n.push_back((*it)->fTailMCSmeared());
    ne.push_back((*it)->fTailMCSmearedErr());
  }
  delete gFTailMC_;
  gFTailMC_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
				    &(pt3e.front()),&(pt3e.front()),
				    &(ne.front()),&(ne.front()));
  setStyleMC(gFTailMC_);

  // Fill graphs of ftail mc truth
  if( hasMCTruthResponse() ) {
    pt3.clear();
    pt3e.clear();
    n.clear();
    ne.clear();
    pt3.push_back(0.);
    pt3e.push_back(0.);
    n.push_back(mcTruthResponse_->fTailMCSmeared());
    ne.push_back(mcTruthResponse_->fTailMCSmearedErr());
    delete gFTailMCTruth_;
    gFTailMCTruth_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					   &(pt3e.front()),&(pt3e.front()),
					   &(ne.front()),&(ne.front()));
    gFTailMCTruth_->SetMarkerStyle(21);
    gFTailMCTruth_->SetMarkerColor(kRed);
    gFTailMCTruth_->SetLineColor(gFTailMCTruth_->GetMarkerColor());

    // After Gaussian subtraction
    pt3.clear();
    pt3e.clear();
    n.clear();
    ne.clear();
    pt3.push_back(0.);
    pt3e.push_back(0.);
    n.push_back(mcTruthResponse_->fTailMCSmearedNonGauss());
    ne.push_back(mcTruthResponse_->fTailMCSmearedNonGaussErr());
    delete gFTailMCTruthNonGauss_;
    gFTailMCTruthNonGauss_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
						   &(pt3e.front()),&(pt3e.front()),
						   &(ne.front()),&(ne.front()));
    gFTailMCTruthNonGauss_->SetMarkerStyle(27);
    gFTailMCTruthNonGauss_->SetMarkerColor(kGreen);
    gFTailMCTruthNonGauss_->SetLineColor(gFTailMCTruth_->GetMarkerColor());

    //    std::cout << ">>>> " << (mcTruthResponse_->fTailMCSmearedNonGauss()/mcTruthResponse_->fTailMCSmeared()) << std::endl;
  }

  // Extrapolate
  gFTailMC_->Fit(fExMC_,"0QR");
  fExMC_->SetLineColor(gFTailMC_->GetLineColor());

  // Fill graphs of spread of ftail vs pt3 for mc
  TH1* hSpread = new TH1D("hSpread","",10000,-0.1,0.1);
  pt3.clear();
  pt3e.clear();
  n.clear();
  ne.clear();
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    double thres = (*it)->pt3Thres();
    double spread = ((*it)->fTailMCSmeared()-fExMC_->Eval(thres))/fExMC_->Eval(thres);
    pt3.push_back(thres);
    pt3e.push_back(0.);
    n.push_back(spread);
    ne.push_back(0.);
    hSpread->Fill(spread);
  }
  delete gFTailSpreadMC_;
  gFTailSpreadMC_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					  &(pt3e.front()),&(pt3e.front()),
					  &(ne.front()),&(ne.front()));
  setStyleMC(gFTailSpreadMC_);

  // Extrapolation of ftail
  extraMC_ = fExMC_->GetParameter(0);
  extraMCErr_ = hSpread->GetRMS()*extraMC_;


  // Fill graphs of ftail vs pt3 for data
  pt3.clear();
  pt3e.clear();
  n.clear();
  ne.clear();
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    double thres = (*it)->pt3Thres();
    if( thres > minPt3Data && (*it)->fTailData() > 0. ) {
      pt3.push_back(thres);
      pt3e.push_back(0.);
      n.push_back((*it)->fTailData());
      ne.push_back((*it)->fTailDataErr());
    }
  }
  delete gFTailData_;
  gFTailData_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
				      &(pt3e.front()),&(pt3e.front()),
				      &(ne.front()),&(ne.front()));
  setStyleData(gFTailData_);


  // Fit extrapolation to data
  fExData_->SetLineColor(gFTailData_->GetLineColor());
  if( fixDataShape ) {
    fExData_->FixParameter(1,fExMC_->GetParameter(1));
    fExData_->FixParameter(2,fExMC_->GetParameter(2));
    fExData_->FixParameter(3,fExMC_->GetParameter(3));
  }
  gFTailData_->Fit(fExData_,"0QRB");

  // Fill graphs of spread of ftail vs pt3 for data
  hSpread->Reset();
  pt3.clear();
  pt3e.clear();
  n.clear();
  ne.clear();
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    double thres = (*it)->pt3Thres();
    if( thres > minPt3Data && (*it)->fTailData() > 0. ) {
      double spread = ((*it)->fTailData()-fExData_->Eval(thres))/fExData_->Eval(thres);
      pt3.push_back(thres);
      pt3e.push_back(0.);
      n.push_back(spread);
      ne.push_back(0.);
      hSpread->Fill(spread);
    }
  }
  delete gFTailSpreadData_;
  gFTailSpreadData_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					    &(pt3e.front()),&(pt3e.front()),
					    &(ne.front()),&(ne.front()));
  setStyleData(gFTailSpreadData_);

  // Extrapolated ftail
  extraData_ = fExData_->GetParameter(0);
  extraDataErr_ = hSpread->GetRMS()*extraData_;
  delete hSpread;


  // Set delta: absolute difference of extrapolated values
  deltaEx_ = extraData_ - extraMC_;
  // Uncertainty from fitted parameters (no fluctuations with fixed window)
  deltaExErr_ = sqrt( pow(extraDataErr_,2.) + pow(extraMCErr_,2.) );


  // Scaling factor (fmc(0)+delta)/fmc(0)
  double ref = extraMC_;
  double refE = extraMCErr_;
  if( hasMCTruthResponse() && mcTruthRef ) {
    ref = mcTruthResponse_->fTailMCSmeared();
    refE = mcTruthResponse_->fTailMCSmearedErr();
  }
  
  //  std::cout << "(" << etaBin_ << "," << ptBin_ << ")   $" << ref << " \\pm " << refE << "$" << std::endl;

  scalingFactor_ = (ref+deltaEx_)/ref;
  scalingFactorErr_ = sqrt( pow(deltaExErr_/ref,2.) + pow(deltaEx_*refE/ref/ref,2.) );

  std::cout << "PtBin " << ptBin_ << ": " << deltaEx_ << " +/- " << deltaExErr_ << " --> " << scalingFactor_ << " +/- " << scalingFactorErr_ << std::endl;
}




// Constructor for asymmetry
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, double tailWindowDataMin, double tailWindowDataMax, double tailWindowMCMin, double tailWindowMCMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel)
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(pt3Bin), pt3Thres_(pt3Thres),
    coreScalingFactor_(coreScalingFactor) {

  hAsymData_ = readHist(fileNameData,"Data");
  hAsymMC_ = readHist(fileNameMC,"MC");
  hResp_ = 0;
  fGaussMCTruth_ = 0;

  hAsymData_->UseCurrentStyle();
  hAsymMC_->UseCurrentStyle();

  // Smear MC asymmetry
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(hAsymMC_,nSigCore,width,widthErr) ) {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),width,coreScalingFactor);
  } else {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),0.,0.);
  }

  // Get relative number of entries in tail
  tailStartBinMC_ = hAsymMCSmeared_->FindBin(tailWindowMCMin);
  tailEndBinMC_ = hAsymMCSmeared_->FindBin(tailWindowMCMax);
  getFTail(hAsymMC_,hAsymMC_->GetEntries(),tailStartBinMC_,tailEndBinMC_,fNMC_,fNMCErr_);
  getFTail(hAsymMCSmeared_,hAsymMC_->GetEntries(),tailStartBinMC_,tailEndBinMC_,fNMCSmeared_,fNMCSmearedErr_);

  tailStartBinData_ = hAsymData_->FindBin(tailWindowDataMin);
  tailEndBinData_ = hAsymData_->FindBin(tailWindowDataMax);
  getFTail(hAsymData_,hAsymData_->GetEntries(),tailStartBinData_,tailEndBinData_,fNData_,fNDataErr_);

  setStyleDataMarker(hAsymData_);
  setStyleMCFilled(hAsymMC_);
  setStyleMCFilled(hAsymMCSmeared_);
  initBinLabel(adm,jetLabel);
}


// Constructor for MC truth response
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel) 
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(1000), pt3Thres_(1000.),
    coreScalingFactor_(coreScalingFactor) {

  hAsymData_ = 0;
  hAsymMC_ = 0;
  hAsymMCSmeared_ = 0;
  fGaussMCTruth_ = 0;

  
  TH1* h = readMCTruthResponse(fileName,"MCTruthResponse");
  double entries = h->GetEntries();
  util::HistOps::setAxisTitles(h,"Symmetrised Response","GeV","events");
  h->UseCurrentStyle();

  // Smear core
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore,width,widthErr) ) {
    func::smearHistogram(h,hResp_,h->GetEntries(),width,coreScalingFactor);
  } else {
    func::smearHistogram(h,hResp_,h->GetEntries(),0.,0.);
  }
  //  std::cout << "(" << etaBin_ << "," << ptBin_ << ") Entries " << h->GetEntries() << " (" << hResp_->GetEntries() << ")\n";

  delete h;


  // Get relative number of entries in tail
  tailStartBinMC_ = hResp_->FindBin(tailWindowMin);
  tailEndBinMC_ = hResp_->FindBin(tailWindowMax);
  getFTail(hResp_,entries,tailStartBinMC_,tailEndBinMC_,fNMCSmeared_,fNMCSmearedErr_);

  //  std::cout << "  " << tailWindowMin << " (" << tailStartBinMC_ << ")" << " - " << tailWindowMax << " (" << tailEndBinMC_ << ")\n";
  //  std::cout << "  " << fNMCSmeared_ << std::endl;

  // Get tail after Gaussian subtraction (note: this is a Gaussian to the smeared histogram!)
  if( util::HistOps::fitCoreWidth(hResp_,nSigCore,fGaussMCTruth_,width,widthErr) ) {
    // Subtract Gaussian
    double winMin = hResp_->GetXaxis()->GetBinLowEdge(tailStartBinMC_);
    double winMax = hResp_->GetXaxis()->GetBinUpEdge(tailEndBinMC_);
    double gaussContr = fGaussMCTruth_->Integral(winMin,winMax);
    fNMCSmearedNonGauss_ = std::max(fNMCSmeared_-gaussContr,0.);
    fNMCSmearedNonGaussErr_ = fNMCSmearedErr_;
    // Style
    fGaussMCTruth_->SetLineColor(kRed);
    fGaussMCTruth_->SetLineWidth(1);
  }

  setStyleDataMarker(hResp_);
  initBinLabel(adm,jetLabel,true);
}


// ------------------------------------------------------------------------------------
void Pt3Bin::initBinLabel(const sampleTools::BinningAdmin* adm, const TString &jetLabel, bool isMCTruth) {
  if( isMCTruth ) binLabel_ = util::LabelFactory::createPaveText(3,gBINLABEL_WIDTH);
  else binLabel_ = util::LabelFactory::createPaveText(4,gBINLABEL_WIDTH);
  binLabel_->AddText(jetLabel);
  binLabel_->AddText(util::toTString(adm->etaMin(etaBin_))+" < |#eta| < "+util::toTString(adm->etaMax(etaBin_)));
  binLabel_->AddText(util::toTString(adm->ptMin(etaBin_,ptBin_))+" < p^{ave}_{T} < "+util::toTString(adm->ptMax(etaBin_,ptBin_))+" GeV");
  if( !isMCTruth ) binLabel_->AddText(util::toTString(adm->ptSoftMin(pt3Bin_))+" < p^{rel}_{T,3} < "+util::toTString(adm->ptSoftMax(pt3Bin_)));
}


// ------------------------------------------------------------------------------------
Pt3Bin::~Pt3Bin() {
  if( hAsymData_ ) delete hAsymData_;
  if( hAsymMC_ ) delete hAsymMC_;
  if( hAsymMCSmeared_ ) delete hAsymMCSmeared_;
  if( hResp_ ) delete hResp_;
  if( fGaussMCTruth_ ) delete fGaussMCTruth_;
  if( binLabel_ ) delete binLabel_;
}


// Get asymmetry histograms from Kalibri input files
// ------------------------------------------------------------------------------------
TH1* Pt3Bin::readHist(const TString &fileName, const TString &id) const {

  //  TString histName = "hPtAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,id+histName+"_Pt3Bin"+util::toTString(pt3Bin_));
  h->GetXaxis()->SetRangeUser(-1.,1);
  h->Scale(1./h->Integral("width"));
  util::HistOps::setAxisTitles(h,"Asymmetry","","events",true);      
  h->SetTitle("");

  return h;
}


// Get symmetrized MC truth response histogram from Kalibri input files
// ------------------------------------------------------------------------------------
TH1* Pt3Bin::readMCTruthResponse(const TString &fileName, const TString &id) const {

  TString histName = "hRespSymAbs_"+util::toTString(ptBin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,id+histName);
  h->GetXaxis()->SetRangeUser(0.,2.);
  h->Scale(1./h->Integral("width"));
  util::HistOps::setAxisTitles(h,"Response","","events",true);      
  h->SetTitle("");

  return h;
}


// Gaussian approximation for binomial error
// ------------------------------------------------------------------------------------
void Pt3Bin::getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr) const {
  // Relative number of tail events
  // useage of Integral() because
  // a) possiblity to get #evts in certain interval
  // b) get total #evts consistently taking into account normalizations
  double nTail = h->Integral(start,end);
  double nTotal = h->Integral();
  fTail = nTail/nTotal;
  
  // Total number of events for uncertainty calculation
  // from number of entries (assume no overflow)
  nTotal = entries;
  fTailErr = sqrt( fTail*(1.-fTail)/nTotal );
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMC(const TString &outNameId) const {

  TLegend* leg = util::LabelFactory::createLegendCol(2,gLEG_WIDTH);
  leg->AddEntry(hAsymData_,"Data","P");
  leg->AddEntry(hAsymMC_,"MC","F");

  // Log scale
  util::HistOps::setYRange(hAsymMC_,4,3E-5);
  util::HistOps::setYRange(hAsymData_,4,3E-5);
  hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);

  TString canName = outNameId+"PtAsym";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hAsymMC_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtAsymLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hAsymMC_,4);
  hAsymMC_->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hAsymMC_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"PtAsym.eps","eps");
  
//   // Superimpose fit
//   TF1 *tailMC = mc_->fTail();
//   TF1 *tailData = data_->fTail();
//   can->cd();
//   hAsymMC_->Draw("HISTE");
//   hAsymData_->Draw("PE1same");
//   tailMC->Draw("same");
//   tailData->Draw("same");
//   label_->info(etaBin_,ptBin_)->Draw("same");
//   label_->eta(etaBin_)->Draw("same");
//   label_->legend(hAsymData_,hAsymMC_,tailData,tailMC)->Draw("same");
//   can->SetLogy(0);
//   can->SaveAs(outNameId+"PtAsymFit.eps","eps");

//   // Ratio data / MC
//   TH1 *hRatio = util::HistOps::createRatioPlot(hAsymData_,hAsymMC_);
//   TH1 *hRatioFrame = util::HistOps::createRatioFrame(hAsymData_,"Data / MC",0.,3.);
//   can->cd();
//   hRatioFrame->Draw();
//   hRatio->Draw("PE1same");
//   label_->info(etaBin_,ptBin_)->Draw("same");
//   label_->eta(etaBin_)->Draw("same");
//   can->SetLogy(0);
//   can->SaveAs(outNameId+"PtAsymLogRatio.eps","eps");

//   hRatioFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
//   can->cd();
//   hRatioFrame->Draw();
//   hRatio->Draw("PE1same");
//   label_->info(etaBin_,ptBin_)->Draw("same");
//   label_->eta(etaBin_)->Draw("same");
//   can->SetLogy(0);
//   can->SaveAs(outNameId+"PtAsymRatio.eps","eps");

//   // Bottom ratio plot
//   delete can;
//   can = util::HistOps::createRatioTopCanvas();
//   TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
//   TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hAsymMC_);
//   TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,"Asymmetry","",0.61,1.39);
//   can->cd();
//   bRatioTopFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
//   bRatioTopFrame->GetYaxis()->SetRangeUser(0.,9.5);
//   bRatioTopFrame->Draw("HISTE");
//   hAsymData_->Draw("PE1same");
//   label_->info(etaBin_,ptBin_)->Draw("same");
//   label_->eta(etaBin_)->Draw("same");
//   label_->legend(hAsymData_,hAsymMC_)->Draw("same");
//   bRatioBottomPad->Draw();
//   bRatioBottomPad->cd();
//   bRatioBottomFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
//   bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
//   bRatioBottomFrame->Draw();
//   hRatio->Draw("PE1same");
//   can->SaveAs(outNameId+"PtAsymBottomRatio.eps","eps");

  delete leg;
  delete can;

  hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);
  hAsymData_->GetXaxis()->SetRangeUser(-1.,1.);
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMCSmeared(const TString &outNameId) const {

  TLegend* leg = util::LabelFactory::createLegendCol(2,gLEG_WIDTH);
  leg->AddEntry(hAsymData_,"Data","P");
  leg->AddEntry(hAsymMCSmeared_,"Reweighted MC","F");

  // Log scale
  util::HistOps::setYRange(hAsymMCSmeared_,4,3E-5);
  util::HistOps::setYRange(hAsymData_,4,3E-5);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);

  TString canName = outNameId+"PtAsym";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtSmearAsymLog.eps","eps");

  // Including tail window
  TLine* winMCMin = new TLine(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBinMC_),0.,
			      hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBinMC_),10.);
  winMCMin->SetLineWidth(1);
  winMCMin->SetLineStyle(1);
  winMCMin->SetLineColor(kRed);
  TLine* winMCMax = new TLine(hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBinMC_),0.,
			      hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBinMC_),10.);
  winMCMax->SetLineWidth(1);
  winMCMax->SetLineStyle(1);
  winMCMax->SetLineColor(kRed);

  TLine* winDataMin = new TLine(hAsymData_->GetXaxis()->GetBinLowEdge(tailStartBinData_),0.,
				hAsymData_->GetXaxis()->GetBinLowEdge(tailStartBinData_),10.);
  winDataMin->SetLineWidth(1);
  winDataMin->SetLineStyle(1);
  winDataMin->SetLineColor(kRed);
  TLine* winDataMax = new TLine(hAsymData_->GetXaxis()->GetBinUpEdge(tailEndBinData_),0.,
				hAsymData_->GetXaxis()->GetBinUpEdge(tailEndBinData_),10.);
  winDataMax->SetLineWidth(1);
  winDataMax->SetLineStyle(1);
  winDataMax->SetLineColor(kRed);

  TBox* winMC = new TBox(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBinMC_),0.,
			 hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBinMC_),10.);
  winMC->SetLineWidth(1);
  winMC->SetFillStyle(3444);
  winMC->SetLineColor(kRed);
  winMC->SetFillColor(winMC->GetLineColor());

  TLegend* legWin = util::LabelFactory::createLegendCol(3,gLEG_WIDTH);
  legWin->AddEntry(hAsymData_,"Data","P");
  legWin->AddEntry(hAsymMCSmeared_,"Reweighted MC","F");
  legWin->AddEntry(winMC,"Window","F");

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  winMC->Draw("same");
  //winMCMin->Draw("same");
  //winMCMax->Draw("same");
  //winDataMin->Draw("same");
  //winDataMax->Draw("same");
  binLabel_->Draw("same");
  legWin->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtSmearAsymTail.eps","eps");
  delete legWin;
  delete winMC;
  delete winMCMin;
  delete winMCMax;

  // Linear scale
  util::HistOps::setYRange(hAsymMCSmeared_,4);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"PtSmearAsym.eps","eps");

  delete leg;
  delete can;

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);
  hAsymData_->GetXaxis()->SetRangeUser(-1.,1.);
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotMCTruthResponse(const TString &outNameId) const {

  // Log scale
  TBox* winMC = new TBox(hResp_->GetXaxis()->GetBinLowEdge(tailStartBinMC_),0.,
			 hResp_->GetXaxis()->GetBinUpEdge(tailEndBinMC_),10.);
  winMC->SetLineWidth(1);
  winMC->SetFillStyle(3444);
  winMC->SetLineColor(kRed);
  winMC->SetFillColor(winMC->GetLineColor());

  TLine* winMCMin = new TLine(hResp_->GetXaxis()->GetBinLowEdge(tailStartBinMC_),0.,
			      hResp_->GetXaxis()->GetBinLowEdge(tailStartBinMC_),10.);
  winMCMin->SetLineWidth(1);
  winMCMin->SetLineStyle(2);
  winMCMin->SetLineColor(kRed);
  TLine* winMCMax = new TLine(hResp_->GetXaxis()->GetBinUpEdge(tailEndBinMC_),0.,
			      hResp_->GetXaxis()->GetBinUpEdge(tailEndBinMC_),10.);
  winMCMax->SetLineWidth(1);
  winMCMax->SetLineStyle(2);
  winMCMax->SetLineColor(kRed);

  TLegend* leg = util::LabelFactory::createLegendCol(2,gLEG_WIDTH);
  leg->AddEntry(hResp_,"Reweighted MC","P");
  leg->AddEntry(winMC,"Window","F");

  util::HistOps::setYRange(hResp_,3,3E-5);
  hResp_->GetXaxis()->SetRangeUser(0.,2.);

  TString canName = outNameId+"_MCTruthResponse";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hResp_->Draw("PE1");
  //  fGaussMCTruth_->Draw("same");
  winMC->Draw("same");
//   winMCMin->Draw("same");
//   winMCMax->Draw("same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"_MCTruthResponseLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hResp_,3);
  hResp_->GetXaxis()->SetRangeUser(0.4,1.6);
  can->cd();
  hResp_->Draw("PE1");
  //  fGaussMCTruth_->Draw("same");
  binLabel_->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"_MCTruthResponse.eps","eps");

  hResp_->GetXaxis()->SetRangeUser(0.,2.);

  delete winMC;
  delete winMCMin;
  delete winMCMax;
  delete leg;
  delete can;
}


void setStyleDataMarker(TH1* h) {
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
}

void setStyleMCFilled(TH1* h) {
  h->SetLineWidth(1);
  h->SetMarkerStyle(1);
  h->SetFillColor(38);
}

void setStyleData(TGraphAsymmErrors* g) {
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
}

void setStyleMC(TGraphAsymmErrors* g) {
  g->SetMarkerStyle(25);
  g->SetMarkerColor(kBlue);
  g->SetLineColor(kBlue);
}




////////////////////// COMBINE FINAL RESULTS /////////////////////////////////////


TGraphAsymmErrors* relUncertainty(const TH1* hNom, int color, double weight, const TH1* hVarUp, const TH1* hVarDown = 0) {

  std::vector<double> x;
  std::vector<double> xe;
  std::vector<double> y;
  std::vector<double> ye;
  for(int bin = 1; bin <= hNom->GetNbinsX(); ++bin) {
    x.push_back(hNom->GetBinCenter(bin));
    xe.push_back(0.5*hNom->GetBinWidth(bin));
    y.push_back(0.);
    double nom = hNom->GetBinContent(bin);
    double var = hVarUp->GetBinContent(bin);
    double err = std::abs((var-nom)/nom);
    if( hVarDown ) {
      var = hVarDown->GetBinContent(bin);
      err += std::abs((var-nom)/nom);
      err /= 2.;
    }
    err *= weight;
    ye.push_back(err);
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
					       &(xe.front()),&(xe.front()),&(ye.front()),&(ye.front()));
  g->SetFillStyle(1001);
  g->SetFillColor(color);

  return g;
}


TGraphAsymmErrors* totalUncertainty(const std::vector<TGraphAsymmErrors*> &uncerts, int color) {

  TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(uncerts.at(0)->Clone());
  for(int i = 0; i < g->GetN(); ++i) {
    double eu2 = 0.;
    double ed2 = 0.;
    for(unsigned int j = 0; j < uncerts.size(); ++j) {
      eu2 += pow(uncerts.at(j)->GetEYhigh()[i],2); 
      ed2 += pow(uncerts.at(j)->GetEYlow()[i],2);
    }
    g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],sqrt(ed2),sqrt(eu2));
  }
  g->SetFillStyle(1001);
  g->SetFillColor(color);

  return g;
}


TGraphAsymmErrors* uncertaintyBand(const TH1* hNom, const TGraphAsymmErrors* gRelUncert) {
  TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(gRelUncert->Clone());
  for(int i = 0; i < g->GetN(); ++i) {
    double y = hNom->GetBinContent(1+i);
    g->SetPoint(i,g->GetX()[i],y);
    double eu = (g->GetEYhigh()[i])*y;
    double ed = (g->GetEYlow()[i])*y;
    g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],ed,eu);
  }

  return g;
}
  


void plotFinalResult() {
  util::StyleSettings::paper();
  gErrorIgnoreLevel = 1001;

  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("../sampleTools/BinningAdminTailsRebinned.cfg");  

  TString fileNamePrefix = "results/Tail";
  TString outNamePrefix = "TailScalingFactors_PF";

  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    // Read nominal scaling factors
    TString histName = "hScaleFrame_Eta"+util::toTString(etaBin);
    TH1* hScaleFrame = util::FileOps::readTH1(fileNamePrefix+"_PF.root",histName,histName);
    histName = "hScale_Eta"+util::toTString(etaBin);
    TH1* hScaleNom = util::FileOps::readTH1(fileNamePrefix+"_PF.root",histName,histName+"_Nominal");

    // Variations
    std::vector<TGraphAsymmErrors*> uncerts;
    TH1* hScaleVarClosure = util::FileOps::readTH1(fileNamePrefix+"VarClosure_PF.root",histName,histName+"_VarClosure");
    uncerts.push_back(relUncertainty(hScaleNom,14,0.5,hScaleVarClosure));
    TH1* hScaleVarExtra = util::FileOps::readTH1(fileNamePrefix+"VarExtrapolation_PF.root",histName,histName+"_VarExtrapolation");
    uncerts.push_back(relUncertainty(hScaleNom,5,0.5,hScaleVarExtra));

    TGraphAsymmErrors* gTotal = totalUncertainty(uncerts,38);
    TGraphAsymmErrors* gAbs = uncertaintyBand(hScaleNom,gTotal);

    // Label
    TPaveText* label = util::LabelFactory::createPaveText(2,gBINLABEL_WIDTH);
    label->AddText("Anti-k_{T} (d=0.5) PF Jets");
    label->AddText(util::toTString(binAdm->etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdm->etaMax(etaBin)));

    TLegend* legScale = util::LabelFactory::createLegendCol(2,gLEG_WIDTH);
    legScale->AddEntry(hScaleNom,"Scaling Factors","P");
    legScale->AddEntry(gAbs,"Syst. Uncertainty","F");

    TLegend* legUncert = util::LabelFactory::createLegendCol(uncerts.size()+1,gLEG_WIDTH);
    legUncert->AddEntry(uncerts.at(0),"Non-Closure","F");
    legUncert->AddEntry(uncerts.at(1),"Extrapolation","F");
    //    legUncert->AddEntry(uncerts.at(2),"Reweighting","F");
    legUncert->AddEntry(gTotal,"Total","F");


    // Plot scaling factors and total uncertainty    
    hScaleFrame->GetYaxis()->SetRangeUser(0.01,2.99);
    TCanvas* canScale = new TCanvas("canScale"+util::toTString(etaBin),"Scaling Factors Eta "+util::toTString(etaBin),500,500);
    canScale->cd();
    hScaleFrame->Draw("HIST");
    gAbs->Draw("E2same");
    hScaleFrame->Draw("HISTsame");
    hScaleNom->Draw("PE1same");
    canScale->SetLogx();
    label->Draw("same");
    legScale->Draw("same");
    canScale->SaveAs(outNamePrefix+"_Eta"+util::toTString(etaBin)+".eps","eps");


    // Plot relative uncertainties
    TH1* hUncertsFrame = static_cast<TH1D*>(hScaleFrame->Clone("hUncertsFrame"+util::toTString(etaBin)));
    for(int bin = 1; bin <= hUncertsFrame->GetNbinsX(); ++bin) {
      hUncertsFrame->SetBinContent(bin,0.);
    }
    hUncertsFrame->GetYaxis()->SetRangeUser(-0.99,1.49);
    hUncertsFrame->GetYaxis()->SetTitle("Relative Uncertainties");
    TCanvas* canRelUncerts = new TCanvas("canRelUncerts"+util::toTString(etaBin),"Relative Uncertainties Eta "+util::toTString(etaBin),500,500);
    canRelUncerts->cd();
    hUncertsFrame->Draw("HIST");
    gTotal->Draw("E2same");
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      uncerts.at(i)->Draw("E2same");
    }
    hUncertsFrame->Draw("HISTsame");
    label->Draw("same");
    legUncert->Draw("same");
    canRelUncerts->SetLogx();
    canRelUncerts->SaveAs(outNamePrefix+"Uncertainties_Eta"+util::toTString(etaBin)+".eps","eps");


    // Print factors with total uncertainty (stat + syst)
    if( etaBin == 0 ) {
      std::cout << "\\hline\n";
      std::cout << "\\eta & \\ptave & Factor & Uncertainty down & Uncertainty up \\\\\n";
      std::cout << "\\hline\n";
    }
    for(int bin = 1; bin <= hScaleNom->GetNbinsX(); ++bin) {
      std::cout << std::setprecision(1) << util::toTString(binAdm->etaMin(etaBin)) << " -- " << util::toTString(binAdm->etaMax(etaBin)) << " & ";
      std::cout << std::setprecision(0) << util::toTString(binAdm->ptMin(etaBin,bin-1)) << " -- " << util::toTString(binAdm->ptMax(etaBin,bin-1)) << " & ";
      std::cout << std::setprecision(4) << hScaleNom->GetBinContent(bin) << " & ";
      double estat = hScaleNom->GetBinError(bin);
      double esystd = gAbs->GetEYlow()[bin-1];
      double esystu = gAbs->GetEYhigh()[bin-1];
      std::cout << sqrt( estat*estat + esystd*esystd ) << " & ";
      std::cout << sqrt( estat*estat + esystu*esystu ) << " \\\\\n";
    }    
    std::cout << "\\hline\n";
    
  }

  delete binAdm;
}
