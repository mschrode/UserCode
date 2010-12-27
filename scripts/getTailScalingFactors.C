// $Id: getTailScalingFactors.C,v 1.2 2010/12/27 16:45:29 mschrode Exp $

#include <vector>

#include <cmath>
#include <iostream>
#include <vector>

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



///////////////////////// TYPE DEFINITIONS /////////////////////////////////////////////////


void setStyleDataMarker(TH1* h);
void setStyleMCFilled(TH1* h);
void setStyleData(TGraphAsymmErrors* g);
void setStyleMC(TGraphAsymmErrors* g);


// ------------------------------------------------------------------------------------
class Pt3Bin {
public:
  Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax);
  ~Pt3Bin();

  double pt3Thres() const { return pt3Thres_; }

  double coreScalingFactor() const { return coreScalingFactor_; }
  double fTailData() const { return fNData_; }
  double fTailDataErr() const { return fNDataErr_; }
  double fTailMC() const { return fNMC_; }
  double fTailMCErr() const { return fNMCErr_; }
  double fTailMCSmeared() const { return fNMCSmeared_; }
  double fTailMCSmearedErr() const { return fNMCSmearedErr_; }

  void plotAsymmetryDataMC(const TString &outNameId) const;
  void plotAsymmetryDataMCSmeared(const TString &outNameId) const;


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const unsigned int pt3Bin_;
  const double pt3Thres_;

  double coreScalingFactor_;
  int tailStartBin_;
  int tailEndBin_;
  double fNMCSmeared_;
  double fNMCSmearedErr_;
  double fNMC_;
  double fNMCErr_;
  double fNData_;
  double fNDataErr_;

  TH1* hAsymData_;
  TH1* hAsymMC_;
  TH1* hAsymMCSmeared_;

  TH1* readHist(const TString &fileName, const TString &id);  
  void init(double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax);
  void getFTail(const TH1* h, int start, int end, double &fTail, double &fTailErr) const;
};



// ------------------------------------------------------------------------------------
class EtaPtBin {
public:
  EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double nSigTailWindowMin, double nSigTailWindowMax);
  ~EtaPtBin();


  unsigned int nPt3Bins() const { return pt3Bins_.size(); }
  
  void plotAsymmetryDistributions(const TString &outNameId) const;
  void plotExtrapolation(const TString &outNameId) const;

  void addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC);
  void extrapolate(double minPt3Data);


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const double nSigCore_;
  const double coreScalingFactor_;
  const double nSigTailWindowMin_;
  const double nSigTailWindowMax_;
  const double exMin_;
  const double exMax_;
  
  double tailWindowMin_;
  double tailWindowMax_;
  TGraphAsymmErrors* gFTailMC_;
  TGraphAsymmErrors* gFTailData_;
  TF1* fExMC_;
  TF1* fExData_;

  std::vector<Pt3Bin*> pt3Bins_;

  TString binId() const { return "_EtaBin"+util::toTString(etaBin_)+"_PtBin"+util::toTString(ptBin_); }
  TString binId(unsigned int pt3BinIdx) const { return binId()+"_Pt3Bin"+util::toTString(pt3BinIdx)+"_"; }
};





//////////////////////////////// MAIN ROUTINE ///////////////////////////////////////////


// ------------------------------------------------------------------------------------
void getTailScalingFactors() {

  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting up parameters" << std::endl;

  const TString id = "Test";
  const double coreScaling = 0.041;
  const double nSigCore = 2.;
  const double nSigTailStart = 4.;
  const double nSigTailEnd = 7.;
  const double minPt3Data = 0.;

  util::StyleSettings::paper();
  gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("input/BinningAdminTailsMC.cfg");



  // +++++ Bins +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating bins" << std::endl;

  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = 1;//binAdm->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {

      // Create eta-pt bin
      EtaPtBin* bin = new EtaPtBin(etaBin,ptBin,nSigCore,coreScaling,nSigTailStart,nSigTailEnd);

      // Add pt3 bins
      TString fileNameData = "~/results/ResolutionFit/TailScaling/Tails_PF_DataFineBins_Eta";
      fileNameData += etaBin;
      TString fileNameMC = "~/results/ResolutionFit/TailScaling/Tails_PF_MCFineBins_Eta";
      fileNameMC += etaBin;
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

      etaPtBins.push_back(bin);
    }
  }



  // +++++ Asymmetry and extrapolation plots ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating asymmetry and extrapolation plots" << std::endl;

  for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
      it != etaPtBins.end(); ++it ) {

    // Asymmetry
    (*it)->plotAsymmetryDistributions(id);
    
    // Extrapolation
    (*it)->extrapolate(minPt3Data);  
    (*it)->plotExtrapolation(id);  
  }
}










//////////////////////////////// IMPLEMENTATIONS ///////////////////////////////////////


// ------------------------------------------------------------------------------------
EtaPtBin::EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double nSigTailWindowMin, double nSigTailWindowMax)
  : etaBin_(etaBin), ptBin_(ptBin),
    nSigCore_(nSigCore), coreScalingFactor_(coreScalingFactor),
    nSigTailWindowMin_(nSigTailWindowMin), nSigTailWindowMax_(nSigTailWindowMax),
    exMin_(0.), exMax_(0.32) {

  tailWindowMin_ = 0.;
  tailWindowMax_ = 0.;

  gFTailData_ = new TGraphAsymmErrors(0);
  gFTailMC_ = new TGraphAsymmErrors(0);

  fExMC_ = new TF1("fExMC"+binId(),"[0] + sq([1])*x + [2]*sq(x)",exMin_,exMax_);
  fExMC_->SetParameter(0,0.005);
  fExMC_->SetParameter(1,0.);
  fExMC_->SetParameter(2,0.05);
  fExMC_->SetLineWidth(1);

  fExData_ = static_cast<TF1*>(fExMC_->Clone("fExData"+binId()));
  fExData_->SetLineStyle(2);
}


// ------------------------------------------------------------------------------------
EtaPtBin::~EtaPtBin() {
  for(std::vector<Pt3Bin*>::iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    delete *it;
  }
  delete gFTailMC_;
  delete gFTailData_;
  delete fExMC_;
  delete fExData_;
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
void EtaPtBin::plotExtrapolation(const TString &outNameId) const {

  // Extrapolation MC Closure
  TCanvas* can = new TCanvas("can","Number of events",500,500);
  can->cd();
  TH1* hFrame = new TH1D("hFrame",";p_{T,3} Threshold;f_{Asym}",1000,0.,0.44);
  hFrame->GetYaxis()->SetRangeUser(0.,2.3*gFTailMC_->GetY()[gFTailMC_->GetN()-1]);
  hFrame->Draw();
  fExMC_->Draw("same");
  gFTailMC_->Draw("PE1same");
  can->SaveAs(outNameId+binId()+"_ExtrapolationMCClosure.eps","eps");

  // Extrapolation MC Closure + Data
  can->cd();
  hFrame->Draw();
  fExMC_->Draw("same");
  gFTailMC_->Draw("PE1same");
  gFTailData_->Draw("PE1same");
  can->SaveAs(outNameId+binId()+"_Extrapolation.eps","eps");

  // Extrapolation MC Closure + Data + Shifted Extrapolation
  can->cd();
  hFrame->Draw();
  fExMC_->Draw("same");
  fExData_->Draw("same");
  gFTailMC_->Draw("PE1same");
  gFTailData_->Draw("PE1same");
  can->SaveAs(outNameId+binId()+"_Extrapolation2.eps","eps");
  
  delete hFrame;
  delete can;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC) {
  // If has not been done yet, find asymmetry corresponding
  // to tail window borders from Gaussian fit to MC asymmetry
  // in this pt3 bin
  if( tailWindowMin_ == tailWindowMax_ ) {
    TString histName = "hPtAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
    TH1* h = util::FileOps::readTH1(fileNameMC,histName);
    h->GetXaxis()->SetRangeUser(-1.,1);
    double width = 0.;
    double widthErr = 1000.;
    if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
      tailWindowMin_ = std::abs(nSigTailWindowMin_*width);
      tailWindowMax_ = std::abs(nSigTailWindowMax_*width);
    }
    delete h;
  }

  // Add pt3 bin
  pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bins_.size(),thres,nSigCore_,coreScalingFactor_,tailWindowMin_,tailWindowMax_));
}


// ------------------------------------------------------------------------------------
void EtaPtBin::extrapolate(double minPt3Data) {

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

  // Extrapolate
  gFTailMC_->Fit(fExMC_,"0QR");
  fExMC_->SetLineColor(gFTailMC_->GetLineColor());


  // Fill graphs of ftail vs pt3 for data
  pt3.clear();
  pt3e.clear();
  n.clear();
  ne.clear();
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    double thres = (*it)->pt3Thres();
    if( thres > minPt3Data ) {
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
  fExData_->SetLineColor(gFTailMC_->GetLineColor());
  fExData_->FixParameter(1,fExMC_->GetParameter(1));
  fExData_->FixParameter(2,fExMC_->GetParameter(2));
  gFTailData_->Fit(fExData_,"0QRB");
}




// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax)
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(pt3Bin), pt3Thres_(pt3Thres),
    coreScalingFactor_(coreScalingFactor) {

  hAsymData_ = readHist(fileNameData,"Data");
  hAsymMC_ = readHist(fileNameMC,"MC");

  init(nSigCore,coreScalingFactor_,tailWindowMin,tailWindowMax);

  setStyleDataMarker(hAsymData_);
  setStyleMCFilled(hAsymMC_);
  setStyleMCFilled(hAsymMCSmeared_);
}


// ------------------------------------------------------------------------------------
Pt3Bin::~Pt3Bin() {
  delete hAsymData_;
  delete hAsymMC_;
  delete hAsymMCSmeared_;
}


// Get asymmetry histograms from Kalibri input files
// ------------------------------------------------------------------------------------
TH1* Pt3Bin::readHist(const TString &fileName, const TString &id) {

  TString histName = "hPtAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,id+histName+"_Pt3Bin"+util::toTString(pt3Bin_));
  h->GetXaxis()->SetRangeUser(-1.,1);
  util::HistOps::setAxisTitles(h,"Asymmetry","","events",true);      
  h->SetTitle("");

  return h;
}


// ------------------------------------------------------------------------------------
void Pt3Bin::init(double nSigCore, double coreScalingFactor, double tailWindowMin, double tailWindowMax) {

  // Smear MC asymmetry
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(hAsymMC_,nSigCore,width,widthErr) ) {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),width,coreScalingFactor);
  } else {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),0.,0.);
  }

  // Get relative number of entries in tail
  tailStartBin_ = hAsymMCSmeared_->FindBin(tailWindowMin);
  tailEndBin_ = hAsymMCSmeared_->FindBin(tailWindowMax);
  getFTail(hAsymMCSmeared_,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
  getFTail(hAsymMC_,tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
  getFTail(hAsymData_,tailStartBin_,tailEndBin_,fNData_,fNDataErr_);
}


// Gaussian approximation for binomial error
// ------------------------------------------------------------------------------------
void Pt3Bin::getFTail(const TH1* h, int start, int end, double &fTail, double &fTailErr) const {
  // Relative number of tail events
  // useage of Integral() because
  // a) possiblity to get #evts in certain interval
  // b) get total #evts consistently taking into account normalizations
  double nTail = h->Integral(start,end);
  double nTotal = h->Integral();
  fTail = nTail/nTotal;
  
  // Total number of events for uncertainty calculation
  // from number of entries (assume no overflow)
  nTotal = h->GetEntries();
  fTailErr = sqrt( fTail*(1.-fTail)/nTotal );
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMC(const TString &outNameId) const {

  // Log scale
  util::HistOps::setYRange(hAsymMC_,3,3E-5);
  util::HistOps::setYRange(hAsymData_,3,3E-5);
  hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);

  TString canName = outNameId+"PtAsym";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hAsymMC_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtAsymLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hAsymMC_,3);
  hAsymMC_->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hAsymMC_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
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

  delete can;

  hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);
  hAsymData_->GetXaxis()->SetRangeUser(-1.,1.);
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMCSmeared(const TString &outNameId) const {

  // Log scale
  util::HistOps::setYRange(hAsymMCSmeared_,3,3E-5);
  util::HistOps::setYRange(hAsymData_,3,3E-5);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);

  TString canName = outNameId+"PtAsym";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtSmearAsymLog.eps","eps");

  // Including tail window
  TLine* winMin = new TLine(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
			    hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),10.);
  winMin->SetLineWidth(1);
  winMin->SetLineStyle(2);
  winMin->SetLineColor(kRed);
  TLine* winMax = new TLine(hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_),0.,
			    hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_),10.);
  winMax->SetLineWidth(1);
  winMax->SetLineStyle(2);
  winMax->SetLineColor(kRed);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  winMin->Draw("same");
  winMax->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtSmearAsymTail.eps","eps");
  delete winMin;
  delete winMax;

  // Linear scale
  util::HistOps::setYRange(hAsymMCSmeared_,3);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"PtSmearAsym.eps","eps");

  delete can;

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);
  hAsymData_->GetXaxis()->SetRangeUser(-1.,1.);
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
