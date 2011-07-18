// $Id: getTailScalingFactors.C,v 1.25 2011/03/08 15:21:32 mschrode Exp $

#define UTILS_AS_HEADER_FILE

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TArrow.h"
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



const bool DEBUG = false;
const double BINLABEL_WIDTH = -0.48;
const double LEG_WIDTH = 0.48;
const TString PT3RELVAR = "p^{rel}_{T,3} Threshold";
const double PT3PLOTMAX = 0.29;
const TString FASYM = "f_{asym}";
const TString FASYMMC = "f^{mc}_{asym}";
const TString FASYMDATA = "f^{data}_{asym}";
const bool SHOW_HEADER = true;
const TString HEADER = SHOW_HEADER ? "CMS 2010,  #sqrt{s} = 7 TeV,  L = 36 pb^{ -1}" : "";
const TString LUMI_LABEL = "#sqrt{s} = 7 TeV,  L = 36 pb^{ -1}";


///////////////////////// TYPE DEFINITIONS /////////////////////////////////////////////////

class Pt3Bin;
class EtaPtBin;
typedef std::vector<EtaPtBin*>::const_iterator EtaPtBinConstIt;

void setStyleDataMarker(TH1* h);
void setStyleMCFilled(TH1* h);
void setStyleData(TGraphAsymmErrors* g);
void setStyleMC(TGraphAsymmErrors* g);

void printBin(unsigned int etaBin, unsigned int ptBin, const sampleTools::BinningAdmin* adm);
void printWindowBorders(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm);
void printMCClosure(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm);
void printExtrapolation(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm);


// ------------------------------------------------------------------------------------
class Pt3Bin {
public:
  Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel);
  Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double windowRespMin, double windowRespMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel); 
  Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel);
 ~Pt3Bin();

  double pt3Thres() const { return pt3Thres_; }
  int nPtAsymBins() const { return hAsymMC_->GetNbinsX(); }

  double ptAveMeanData() const { return ptAveMeanData_; }
  double ptAveMeanDataErr() const { return ptAveMeanDataErr_; }

  double coreScalingFactor() const { return coreScalingFactor_; }
  double sigma() const { return sig_; }
  double sigmaSmeared() const { return sigSmeared_; }

  double fTailData() const { return fNData_; }
  double fTailDataErr() const { return fNDataErr_; }
  double fTailMC() const { return fNMC_; }
  double fTailMCErr() const { return fNMCErr_; }
  double fTailMCSmeared() const { return fNMCSmeared_; }
  double fTailMCSmearedErr() const { return fNMCSmearedErr_; }
  double fTailMCSmearedNonGauss() const { return fNMCSmearedNonGauss_; }
  double fTailMCSmearedNonGaussErr() const { return fNMCSmearedNonGaussErr_; }
  double fTailMCSmearedGauss() const { return fNMCSmearedGauss_; }

  void plotAsymmetryDataMC(const TString &outNameIdw) const;
  void plotAsymmetryDataMCSmeared(const TString &outNameId, double nSigTailStart, double nSigTailEnd) const;
  void plotToyAsymmetry(const TString &outNameId) const;
  void plotSymMCTruthResponse(const TString &outNameId) const;
  void plotSpectra(const TString &outNameId) const;


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const unsigned int pt3Bin_;
  const double pt3Thres_;

  double coreScalingFactor_;
  double sig_;
  double sigSmeared_;
  int tailStartBin_;
  int tailEndBin_;
  double ptAveMeanData_;
  double ptAveMeanDataErr_;
  double fNMCSmeared_;
  double fNMCSmearedErr_;
  double fNMCSmearedNonGauss_;
  double fNMCSmearedNonGaussErr_;
  double fNMCSmearedGauss_;
  double fNMC_;
  double fNMCErr_;
  double fNData_;
  double fNDataErr_;

  TH1* hAsymData_;
  TH1* hAsymMC_;
  TH1* hAsymMCSmeared_;
  TH1* hResp_;
  TH1* hSymResp_;
  TF1* fGaussMCTruth_;
  TH1* hPtAveSpecData_;
  TH1* hPtAveSpecMC_;

  TPaveText* binLabel_;

  TH1* readHist(const TString &fileName, const TString &id) const;  
  TH1* readMCTruthResponse(const TString &fileName, const TString &type) const;
  //  void init(double nSigCore, double coreScalingFactor, double tailWindowDataMin, double tailWindowDataMax, double tailWindowMCMin, double tailWindowMCMax);
  void initBinLabel(const sampleTools::BinningAdmin* adm, const TString &jetLabel, bool isMCTruth = false);
  void getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr, double nExtraTailEvts = 0) const;
};



// ------------------------------------------------------------------------------------
class EtaPtBin {
public:
  EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor,const sampleTools::BinningAdmin* adm, const TString &jetLabel);
  ~EtaPtBin();


  unsigned int nPt3Bins() const { return pt3Bins_.size(); }
  unsigned int etaBin() const { return etaBin_; }
  unsigned int ptBin() const { return ptBin_; }

  double ptAveMeanData(int i) const { return pt3Bins_.at(i)->ptAveMeanData(); }
  double ptAveMeanDataErr(int i) const { return pt3Bins_.at(i)->ptAveMeanDataErr(); }

  double tailWindowMin() const { return tailWindowMin_; }
  double tailWindowMax() const { return tailWindowMax_; }
  double tailWindowEffMin() const { return tailWindowEffMin_; }
  double tailWindowEffMax() const { return tailWindowEffMax_; }
  double sigma(int i) const { return pt3Bins_.at(i)->sigma(); }
  double sigmaSmeared(int i) const { return pt3Bins_.at(i)->sigmaSmeared(); }

  double fTailMCSmeared(int i) const { return pt3Bins_.at(i)->fTailMCSmeared(); }
  double fTailMCSmearedErr(int i) const { return pt3Bins_.at(i)->fTailMCSmearedErr(); }
  double extraMC() const { return extraMC_; }
  double extraMCErr() const { return extraMCErr_; }
  double extraData() const { return extraData_; }
  double extraDataErr() const { return extraDataErr_; }
  double toyMC() const { return hasToyMC() ? toyMC_->fTailMCSmeared() : 0.; }
  double toyMCErr() const { return hasToyMC() ? toyMC_->fTailMCSmearedErr() : 0.; }
  double symResp() const { return hasSymMCTruth() ? symMCTruth_->fTailMCSmeared() : 0.; }
  double symRespErr() const { return hasSymMCTruth() ? symMCTruth_->fTailMCSmearedErr() : 0.; }
  double deltaExtra() const { return deltaEx_; }
  double deltaExtraErr() const { return deltaExErr_; }
  double scalingFactor() const { return scalingFactor_; }
  double scalingFactorErr() const { return scalingFactorErr_; }
  
  void plotAsymmetryDistributions(const TString &outNameId, double nSigTailStart, double nSigTailEnd) const;
  void plotSpectra(const TString &outNameId) const;
  void plotExtrapolation(const TString &outNameId) const;
  void plotMCTruth(const TString &outNameId) const;

  bool hasToyMC() const { return hasToyMC_; }
  bool hasSymMCTruth() const { return hasSymMCTruth_; }

  void addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC);
  void addMCTruthForToyAsym(const TString &fileName);
  void addSymMCTruthResponse(const TString &fileName);
  void findWindow(const TString &fileName, double nSigTailWindowMin, double nSigTailWindowMax);
  void findWindowFEvts(const TString &fileName, double nSigTailWindowMin, double fWin);
  void setWindow(const TString &fileName, unsigned int bin);
  void extrapolate(double minPt3Data, bool fixDataShape, bool useExtrapolatedValue, bool mcTruthRef);


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const double nSigCore_;
  const double coreScalingFactor_;
  const double exMin_;
  const double exMax_;
  const sampleTools::BinningAdmin* binAdm_;
  const TString jetLabel_;

  double tailWindowMin_;
  double tailWindowMax_;
  double tailWindowEffMin_;
  double tailWindowEffMax_;
  int tailMinBin_;
  int tailMaxBin_;
  TGraphAsymmErrors* gFTailMC_;
  TGraphAsymmErrors* gFTailData_;
  TGraphAsymmErrors* gFTailToyMC_; // Asymmetry from toy MC (MC truth --> asymmetry)
  TGraphAsymmErrors* gFTailMCTruth_; // Symmetrized MC truth
  TGraphAsymmErrors* gFTailMCTruthNonGauss_; // Symmetrized MC truth minus Gaussian component
  TGraphAsymmErrors* gFTailMCGauss_;
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
  Pt3Bin* symMCTruth_;
  bool hasSymMCTruth_;
  Pt3Bin* toyMC_;
  bool hasToyMC_;

  TPaveText* binLabel_;

  TString binId() const { return "_EtaBin"+util::toTString(etaBin_)+"_PtBin"+util::toTString(ptBin_); }
  TString binId(unsigned int pt3BinIdx) const { return binId()+"_Pt3Bin"+util::toTString(pt3BinIdx)+"_"; }
  int nPtAsymBins() const;
};





//////////////////////////////// MAIN ROUTINE ///////////////////////////////////////////


// ------------------------------------------------------------------------------------
void getTailScalingFactors() {

  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting up parameters" << std::endl;

  const TString outLabel = "TailSig25-Inf";
  //const TString outLabel = "TailSig25-InfVarExtrapolation";
  //const TString outLabel = "TailSig25-InfVarCoreDown";
  //const TString outLabel = "TailSig25-InfVarCoreUp";
  //const TString outLabel = "TailSig25-InfVarClosure";
  //const TString outLabel = "TailMCClosureSig15-25";

  const double nSigCore = 2.;

  const double nSigTailStart = 2.5;
  const double nSigTailEnd = 1000.;

  const double minPt3Data = 0.;
  const bool useExtrapolatedValue = outLabel.Contains("VarExtrapolation") ? false : true;
  const bool fixDataShape = false;
  const bool mcTruthRef = outLabel.Contains("VarClosure") ? true : false;


  const TString jetAlgo = "PF";

  if( SHOW_HEADER ) {
    util::StyleSettings::paper();
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleX(0.7);
    gStyle->SetTitleH(0.038);
  } else {
    util::StyleSettings::paperNoTitle();
  }

  gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("../sampleTools/BinningAdminTailsRebinned.cfg");

  TString id = outLabel+"_"+jetAlgo;
  TString jetLabel = "Anti-k_{T} (R=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";




  // +++++ Bins +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating bins" << std::endl;

  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  unsigned int globalBin = 0;
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    //etaBin = 0;
    double coreScaling = 0.;
    // From combined photon+jet and dijet results
    if( outLabel.Contains("VarCoreDown") ) {
      // Uncertainties down
      if( etaBin == 0 )      coreScaling = 0.011;
      else if( etaBin == 1 ) coreScaling = 0.012;
      else if( etaBin == 2 ) coreScaling = 0.;
      else if( etaBin == 3 ) coreScaling = 0.010;
    } else if( outLabel.Contains("VarCoreUp") ) {
      // Uncertainties up
      if( etaBin == 0 )      coreScaling = 0.081;
      else if( etaBin == 1 ) coreScaling = 0.097;
      else if( etaBin == 2 ) coreScaling = 0.045;
      else if( etaBin == 3 ) coreScaling = 0.190;
    } else {
      // Nominal
      if( etaBin == 0 )      coreScaling = 0.048;
      else if( etaBin == 1 ) coreScaling = 0.058;
      else if( etaBin == 2 ) coreScaling = 0.;
 	else if( etaBin == 3 ) coreScaling = 0.110;
    }
    
    
    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
      if( DEBUG ) std::cout << "Setting up bin (" << etaBin << ", " << ptBin << ")" << std::endl;
      
      TString path = "~/results/ResolutionFit/TailScalingRA2CleaningV1/";
      TString fileNameData = path+"Tails2010Run_PF_Data_Rebinned_Eta";
      TString fileNameMC = path+"Tails2010Run_PF_MC_Rebinned_Eta";

      fileNameData += etaBin;
      fileNameMC += etaBin;

      // Create eta-pt bin
      EtaPtBin* bin = new EtaPtBin(etaBin,ptBin,nSigCore,coreScaling,binAdm,jetLabel);

      // Define window
      bin->findWindow(fileNameMC+"_PtSoft0.root",nSigTailStart,nSigTailEnd);
      //bin->setWindow(fileNameMC+"_PtSoft0.root",globalBin);
      ++globalBin;

      // Add pt3 bins
      bin->addPt3Bin(0.050,fileNameData+"_PtSoft0.root",fileNameMC+"_PtSoft0.root");
      bin->addPt3Bin(0.075,fileNameData+"_PtSoft1.root",fileNameMC+"_PtSoft1.root");
      bin->addPt3Bin(0.100,fileNameData+"_PtSoft2.root",fileNameMC+"_PtSoft2.root");
      bin->addPt3Bin(0.125,fileNameData+"_PtSoft3.root",fileNameMC+"_PtSoft3.root");
      bin->addPt3Bin(0.150,fileNameData+"_PtSoft4.root",fileNameMC+"_PtSoft4.root");
      bin->addPt3Bin(0.175,fileNameData+"_PtSoft5.root",fileNameMC+"_PtSoft5.root");
      bin->addPt3Bin(0.200,fileNameData+"_PtSoft6.root",fileNameMC+"_PtSoft6.root");
      bin->addPt3Bin(0.225,fileNameData+"_PtSoft7.root",fileNameMC+"_PtSoft7.root");

//       bin->addPt3Bin(0.250,fileNameData+"_PtSoft8.root",fileNameMC+"_PtSoft8.root");
//       bin->addPt3Bin(0.275,fileNameData+"_PtSoft9.root",fileNameMC+"_PtSoft9.root");
//       bin->addPt3Bin(0.300,fileNameData+"_PtSoft10.root",fileNameMC+"_PtSoft10.root");

      // Add symmetrized mc truth response
      //bin->addSymMCTruthResponse(fileNameMC+"_PtSoft0.root");

      // Add mc truth response for toy asymmetry
      //bin->addMCTruthForToyAsym(fileNameMC+"_PtSoft2.root");
      bin->addMCTruthForToyAsym(path+"Tails2010Run_PF_PtGenSel_PtSoftCuts_Rebinned_Eta"+util::toTString(etaBin)+"_PtSoft0.root");

      etaPtBins.push_back(bin);
      if( DEBUG ) std::cout << "Done setting up bin" << std::endl;
    }
  }



  // +++++ Asymmetry and extrapolation plots ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating asymmetry and extrapolation plots" << std::endl;

  for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
      it != etaPtBins.end(); ++it ) {

    if( DEBUG ) std::cout << "Plotting and extrapolation" << std::endl;

//     // Asymmetry
    (*it)->plotAsymmetryDistributions(id,nSigTailStart,nSigTailEnd);
    (*it)->plotMCTruth(id);
    (*it)->plotSpectra(id);
    
    // Extrapolation
    (*it)->extrapolate(minPt3Data,fixDataShape,useExtrapolatedValue,mcTruthRef);  
    (*it)->plotExtrapolation(id);

    // Print fractional number of events  
//     std::cout << std::setprecision(2) << " & $" << (*it)->extraData() << " \\pm " << (*it)->extraDataErr() << "$" << std::flush;
//     std::cout << " & $" << (*it)->extraMC() << " \\pm " << (*it)->extraMCErr() << "$ " << std::endl;
  }



  // +++++ Scaling factors ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating scaling factor plots" << std::endl;

  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    TH1* hPtAveMeanData = new TH1D("hPtAveMeanData_Eta"+util::toTString(etaBin),
				   HEADER+";p^{ave}_{T} Bin Idx;<p^{ave}_{T}> (GeV/c)",
				   binAdm->nPtBins(etaBin),-0.5,binAdm->nPtBins(etaBin)-0.5);
    hPtAveMeanData->SetMarkerStyle(20);

    TH1* hDelta = new TH1D("hDelta_Eta"+util::toTString(etaBin),
			   HEADER+";p^{ave}_{T} (GeV/c);#Delta = f^{Data}_{Asym}(0) - f^{MC}_{Asym}(0)",
			   binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hDelta->SetMarkerStyle(20);

    TH1* hScale = new TH1D("hScale_Eta"+util::toTString(etaBin),
			   HEADER+";p^{ave}_{T} (GeV/c);(f^{Data}_{Asym}(0) - f^{MC}_{Asym}(0)) / f^{MC}_{Asym}(0)",
			   binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hScale->SetMarkerStyle(20);

    for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
	it != etaPtBins.end(); ++it) {
      if( (*it)->etaBin() == etaBin ) {
	int bin = 1+(*it)->ptBin();
	hPtAveMeanData->SetBinContent(bin,(*it)->ptAveMeanData(0));
	hPtAveMeanData->SetBinError(bin,(*it)->ptAveMeanDataErr(0));
	hDelta->SetBinContent(bin,(*it)->deltaExtra());
	hDelta->SetBinError(bin,(*it)->deltaExtraErr());
	hScale->SetBinContent(bin,(*it)->scalingFactor());
	hScale->SetBinError(bin,(*it)->scalingFactorErr());
      }
    }


    // Plots
    TCanvas* can = new TCanvas("can","",500,500);

    TH1* hDeltaFrame = static_cast<TH1D*>(hDelta->Clone("hDeltaFrame_Eta"+util::toTString(etaBin)));
    hDeltaFrame->SetTitle(HEADER);
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
      hScaleFrame->SetBinError(i,0.);
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
    if( etaBin == 0 ) fileMode = "RECREATE";
    TFile outFile(id+".root",fileMode);
    outFile.WriteTObject(hDelta);
    outFile.WriteTObject(hDeltaFrame);
    outFile.WriteTObject(hScale);    
    outFile.WriteTObject(hScaleFrame);   
    outFile.WriteTObject(hPtAveMeanData);
    outFile.Close();
    

    delete hDelta;
    delete hScale;
    delete hDeltaFrame;
    delete hScaleFrame;
    delete hPtAveMeanData;
    delete can;
  }

  printWindowBorders(etaPtBins,binAdm);
  printMCClosure(etaPtBins,binAdm);
  printExtrapolation(etaPtBins,binAdm);

  delete binAdm;
}










//////////////////////////////// IMPLEMENTATIONS ///////////////////////////////////////


// ------------------------------------------------------------------------------------
EtaPtBin::EtaPtBin(unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, const sampleTools::BinningAdmin* binAdm, const TString &jetLabel)
  : etaBin_(etaBin), ptBin_(ptBin),
    nSigCore_(nSigCore), coreScalingFactor_(coreScalingFactor),
    exMin_(0.), exMax_(0.25),
    binAdm_(binAdm), jetLabel_(jetLabel) {

  tailWindowMin_ = 0.;
  tailWindowMax_ = 0.;
  tailMinBin_ = 0;
  tailMaxBin_ = 0;

  gFTailData_ = new TGraphAsymmErrors(0);
  gFTailMC_ = new TGraphAsymmErrors(0);
  gFTailToyMC_ = new TGraphAsymmErrors(0);
  gFTailMCTruth_ = new TGraphAsymmErrors(0);
  gFTailSpreadMC_ = new TGraphAsymmErrors(0);
  gFTailSpreadData_  = new TGraphAsymmErrors(0);
  gFTailMCTruthNonGauss_ = new TGraphAsymmErrors(0);
  gFTailMCGauss_ = new TGraphAsymmErrors(0);

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

  symMCTruth_ = 0;
  hasSymMCTruth_ = false;
  toyMC_ = 0;
  hasToyMC_ = false;

  binLabel_ = util::LabelFactory::createPaveText(5,BINLABEL_WIDTH);
  binLabel_->AddText(LUMI_LABEL);
  binLabel_->AddText(jetLabel_);
  binLabel_->AddText(util::toTString(binAdm_->etaMin(etaBin_))+" < |#eta| < "+util::toTString(binAdm_->etaMax(etaBin_)));
  binLabel_->AddText(util::toTString(binAdm_->ptMin(etaBin_,ptBin_))+" < p^{ave}_{T} < "+util::toTString(binAdm_->ptMax(etaBin_,ptBin_))+" GeV/c");
}


// ------------------------------------------------------------------------------------
EtaPtBin::~EtaPtBin() {
  for(std::vector<Pt3Bin*>::iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    delete *it;
  }
  delete gFTailMC_;
  delete gFTailToyMC_;
  delete gFTailData_;
  delete gFTailToyMC_;
  delete gFTailMCTruth_;
  delete gFTailMCTruthNonGauss_;
  delete gFTailMCGauss_;
  delete gFTailSpreadData_;
  delete gFTailSpreadMC_;
  delete fExMC_;
  delete fExData_;
  if( hasSymMCTruth() ) delete symMCTruth_;
  if( hasToyMC() ) delete toyMC_;
  if( binLabel_ ) delete binLabel_;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotAsymmetryDistributions(const TString &outNameId, double nSigTailStart, double nSigTailEnd) const {
  unsigned int pt3Bin = 0;
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it, ++pt3Bin) {
    (*it)->plotAsymmetryDataMC(outNameId+binId(pt3Bin));
    (*it)->plotAsymmetryDataMCSmeared(outNameId+binId(pt3Bin), nSigTailStart, nSigTailEnd);
  }
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotMCTruth(const TString &outNameId) const {
  if( hasToyMC() ) toyMC_->plotToyAsymmetry(outNameId+binId());
  if( hasSymMCTruth() ) symMCTruth_->plotSymMCTruthResponse(outNameId+binId());
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotSpectra(const TString &outNameId) const {
  unsigned int pt3Bin = 0;
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it, ++pt3Bin) {
    (*it)->plotSpectra(outNameId+binId(pt3Bin));
  }
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotExtrapolation(const TString &outNameId) const {
  if( DEBUG ) std::cout << "Entering EtaPtBin::plotExtrapolation()" << std::endl;

  // Extrapolation MC Closure
  TLegend* leg = util::LabelFactory::createLegendCol(3,LEG_WIDTH);
  leg->AddEntry(gFTailMC_,"Asymmetry","P");
  leg->AddEntry(fExMC_,"Extrapolation","L");
  leg->AddEntry(gFTailToyMC_,"Toy MC","P");

  TCanvas* can = new TCanvas("can","Number of events",500,500);
  can->cd();
  TH1* hFrame = new TH1D("hFrame",HEADER+";"+PT3RELVAR,1000,0.,PT3PLOTMAX);
  hFrame->GetYaxis()->SetRangeUser(0.,2.3*gFTailMC_->GetY()[gFTailMC_->GetN()-1]);

  if( hasToyMC() ) {
    hFrame->GetYaxis()->SetTitle(FASYMMC);
    hFrame->Draw();
    fExMC_->Draw("same");
    if( hasSymMCTruth() ) gFTailMCTruth_->Draw("PE1same");
    gFTailToyMC_->Draw("PE1same");
    gFTailMC_->Draw("PE1same");
    binLabel_->DrawClone("same");
    leg->Draw("same");
    can->SaveAs(outNameId+binId()+"_ExtrapolationMCClosure.eps","eps");
  }
  delete leg;

  // Extrapolation MC + Data
  hFrame->GetYaxis()->SetTitle(FASYM);
  leg = util::LabelFactory::createLegendCol(5,LEG_WIDTH);
  leg->AddEntry(gFTailData_,"Asymmetry Data","P");
  leg->AddEntry(fExData_,"Extrapolation Data","L");
  leg->AddEntry(gFTailMC_,"Asymmetry MC","P");
  leg->AddEntry(fExMC_,"Extrapolation MC","L");
  leg->AddEntry(gFTailMCGauss_,"Asymmetry Gauss","P");

  can->cd();
  hFrame->Draw();
  fExMC_->Draw("same");
  fExData_->Draw("same");
  gFTailMC_->Draw("PE1same");
  gFTailData_->Draw("PE1same");
  gFTailMCGauss_->Draw("PE1same");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_Extrapolation.eps","eps");
  delete leg;


  // Extrapolation MC Closure + Data + Shifted Extrapolation
  if( hasToyMC() ) {
    hFrame->GetYaxis()->SetTitle("( "+FASYMMC+" - Fit ) / Fit");
    leg = util::LabelFactory::createLegendCol(5,LEG_WIDTH);
    leg->AddEntry(gFTailData_,"Asymmetry Data","P");
    leg->AddEntry(fExData_,"Extrapolation Data","L");
    leg->AddEntry(gFTailMC_,"Asymmetry MC","P");
    leg->AddEntry(fExMC_,"Extrapolation MC","L");
    leg->AddEntry(gFTailToyMC_,"Toy MC","P");
    
    can->cd();
    hFrame->Draw();
    fExMC_->Draw("same");
    fExData_->Draw("same");
    gFTailMC_->Draw("PE1same");
    gFTailData_->Draw("PE1same");
    if( hasSymMCTruth() ) {
      gFTailMCTruth_->Draw("PE1same");
      //gFTailMCTruthNonGauss_->Draw("PE1same");
    }
    gFTailToyMC_->Draw("PE1same");
    binLabel_->DrawClone("same");
    leg->Draw("same");
    can->SaveAs(outNameId+binId()+"_Extrapolation2.eps","eps");
    delete leg;
  }

  // Spread of ftail for mc
  double relErr = extraMCErr_/extraMC();
  hFrame->GetYaxis()->SetTitle("( "+FASYMMC+" - Fit ) / Fit");
  hFrame->GetYaxis()->SetRangeUser(std::min(-8.*relErr,0.),15.*relErr);
  for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
    hFrame->SetBinContent(i,0.);
  }
  hFrame->SetLineStyle(2);
  leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(gFTailSpreadMC_,"Asymmetry MC","P");
  leg->AddEntry(hFrame,"Fit","L");

  can->cd();
  hFrame->Draw();
  gFTailSpreadMC_->Draw("PE1same");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_SpreadMC.eps","eps");  
  delete leg;

  // Spread of ftail for data
  relErr = extraDataErr_/extraData();
  hFrame->GetYaxis()->SetRangeUser(std::min(-8.*relErr,0.),15.*relErr);
  hFrame->GetYaxis()->SetTitle("( "+FASYMDATA+" - Fit ) / Fit");
  for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
    hFrame->SetBinContent(i,0.);
  }
  leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(gFTailSpreadData_,"Asymmetry Data","P");
  leg->AddEntry(hFrame,"Fit","L");

  can->cd();
  hFrame->Draw();
  gFTailSpreadData_->Draw("PE1same");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  can->SaveAs(outNameId+binId()+"_SpreadData.eps","eps");  
  delete leg;

  delete hFrame;
  delete can;

  if( DEBUG ) std::cout << "Leaving EtaPtBin::plotExtrapolation()" << std::endl;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addPt3Bin(double thres, const TString &fileNameData, const TString &fileNameMC) {
  pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bins_.size(),thres,nSigCore_,coreScalingFactor_,tailMinBin_,tailMaxBin_,binAdm_,jetLabel_));
}


// ------------------------------------------------------------------------------------
int EtaPtBin::nPtAsymBins() const {
  if( pt3Bins_.size() == 0 ) {
    std::cerr << "ERROR in EtaPtBin::nPtAsymBins(): no asymmetry distributions available." << std::endl;
    std::cerr << "  Use addPt3Bin() at least once before calling this function" << std::endl;
    exit(0);
  }
  return pt3Bins_.at(0)->nPtAsymBins();
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addSymMCTruthResponse(const TString &fileName) {
  if( hasSymMCTruth() ) delete symMCTruth_;
  double winMin = 1.+sqrt(2.)*tailWindowMin_;
  double winMax = 1.+sqrt(2.)*tailWindowMax_;
  symMCTruth_ = new Pt3Bin(fileName,etaBin_,ptBin_,nSigCore_,coreScalingFactor_,winMin,winMax,binAdm_,jetLabel_);
  hasSymMCTruth_ = true;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addMCTruthForToyAsym(const TString &fileName) {
  if( hasToyMC() ) delete toyMC_;
  toyMC_ = new Pt3Bin(fileName,etaBin_,ptBin_,nSigCore_,coreScalingFactor_,nPtAsymBins(),tailMinBin_,tailMaxBin_,binAdm_,jetLabel_);
  hasToyMC_ = true;
}


// Find min and max asymmetry
// The borders are specified in number of sigmas, where sigma is the Gaussian
// width of the *unsmeared* asymmetry distribution from fileName
// ------------------------------------------------------------------------------------
void EtaPtBin::findWindow(const TString &fileName, double nSigTailWindowMin, double nSigTailWindowMax) {
  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);

  TH1* h = util::FileOps::readTH1(fileName,histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
    tailMinBin_ = h->FindBin(std::abs(nSigTailWindowMin*width));
    tailMaxBin_ = h->FindBin(std::abs(nSigTailWindowMax*width));
    tailWindowMin_ = h->GetXaxis()->GetBinLowEdge(tailMinBin_);
    tailWindowMax_ = h->GetXaxis()->GetBinUpEdge(tailMaxBin_);
    tailWindowEffMin_ = tailWindowMin_/width;
    tailWindowEffMax_ = tailWindowMax_/width;
    if( DEBUG ) std::cout << "EtaPtBin::findWindow(): sig(" << etaBin_ << ", " << ptBin_ << "): " << width << " \\pm " << widthErr << std::endl;
  } else {
    std::cerr << "ERROR in EtaPtBin::findWindow(): fit of core region did not converge" << std::endl;
    std::cerr << "  Window is not defined properly" << std::endl;
    tailWindowMin_ = 0.;
    tailWindowMax_ = 0.;
    tailWindowEffMin_ = 0.;
    tailWindowEffMax_ = 0.;
    tailMinBin_ = 0;
    tailMaxBin_ = 0;
  }
  delete h;
  
  if( nSigTailWindowMax < 50. ) binLabel_->AddText(util::toTString(nSigTailWindowMin)+" - "+util::toTString(nSigTailWindowMax)+" #sigma Window");
  else binLabel_->AddText("> "+util::toTString(nSigTailWindowMin)+" #sigma Window");
}


// Find min and max asymmetry
// The min border is defined in number of sigmas, where sigma is the Gaussian
// width of the *unsmeared* asymmetry distribution from fileName. The max border
// is computed such that the window contains a fractional number fWin of events
// ------------------------------------------------------------------------------------
void EtaPtBin::findWindowFEvts(const TString &fileName, double nSigTailWindowMin, double fWin) {
  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);

  TH1* h = util::FileOps::readTH1(fileName,histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
    tailMinBin_ = h->FindBin(std::abs(nSigTailWindowMin*width));
    tailMaxBin_ = tailMinBin_;
    double nTot = h->Integral();
    if( nTot > 0. ) {
      double fTmp = h->Integral(tailMinBin_,tailMaxBin_)/nTot;
      while( tailMaxBin_ < h->GetNbinsX() && fTmp < fWin ) {
	++tailMaxBin_;
	fTmp = h->Integral(tailMinBin_,tailMaxBin_)/nTot;
      }
    }
    tailWindowMin_ = h->GetXaxis()->GetBinLowEdge(tailMinBin_);
    tailWindowMax_ = h->GetXaxis()->GetBinUpEdge(tailMaxBin_);
    tailWindowEffMin_ = tailWindowMin_/width;
    tailWindowEffMax_ = tailWindowMax_/width;
  } else {
    std::cerr << "ERROR in EtaPtBin::findWindowFEvts(): fit of core region did not converge" << std::endl;
    std::cerr << "  Window is not defined properly" << std::endl;
    tailWindowMin_ = 0.;
    tailWindowMax_ = 0.;
    tailWindowEffMin_ = 0.;
    tailWindowEffMax_ = 0.;
    tailMinBin_ = 0;
    tailMaxBin_ = 0;
  }
  delete h;

  binLabel_->AddText(util::toTString(nSigTailWindowMin)+" #sigma + "+util::toTString(100.*fWin)+"% N(evts) Window");
}


// Set min and max asymmetry
// Hack: Set from fixed edges 
// ------------------------------------------------------------------------------------
void EtaPtBin::setWindow(const TString &fileName, unsigned int bin) {
  std::vector<double> min;
  min.push_back(0.128713);
  min.push_back(0.0693069);
  min.push_back(0.0693069);
  min.push_back(0.148515);
  min.push_back(0.0693069);
  min.push_back(0.128713);
  min.push_back(0.128713);

  std::vector<double> max;
  max.push_back(0.247525);
  max.push_back(0.148515);
  max.push_back(0.128713);
  max.push_back(0.267327);
  max.push_back(0.148515);
  max.push_back(0.227723);
  max.push_back(0.247525);

  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  TH1* h = util::FileOps::readTH1(fileName,histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  if( bin < min.size() || bin < max.size() ) {
    tailWindowMin_ = min.at(bin);
    tailWindowMax_ = max.at(bin);
    binLabel_->AddText(util::toTString(tailWindowMin_)+" < A < "+util::toTString(tailWindowMax_));
    tailMinBin_ = h->FindBin(tailWindowMin_+0.001);
    tailMaxBin_ = h->FindBin(tailWindowMax_-0.001);
    tailWindowEffMin_ = 0.;
    tailWindowEffMax_ = 0.;
  } else {
    std::cerr << "ERROR in EtaPtBin::setWindow(): borders not specified properly" << std::endl;
    tailWindowMin_ = 0.;
    tailWindowMax_ = 0.;
    tailWindowEffMin_ = 0.;
    tailWindowEffMax_ = 0.;
    tailMinBin_ = 0;
    tailMaxBin_ = 0;
  }
  delete h;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::extrapolate(double minPt3Data, bool fixDataShape, bool useExtrapolatedValue, bool mcTruthRef) {
  if( DEBUG ) std::cout << "Entering EtaPtBin::extrapolate()" << std::endl;

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

  // Purely Gaussian asymmetry
  n.clear();
  ne.clear();
  for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
      it != pt3Bins_.end(); ++it) {
    n.push_back((*it)->fTailMCSmearedGauss());
    ne.push_back(0.);
  }
  delete gFTailMCGauss_;
  gFTailMCGauss_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					 &(pt3e.front()),&(pt3e.front()),
					 &(ne.front()),&(ne.front()));
  gFTailMCGauss_->SetMarkerStyle(27);
  gFTailMCGauss_->SetMarkerColor(8);
  gFTailMCGauss_->SetLineColor(gFTailMCGauss_->GetMarkerColor());


  // Fill graphs of ftail from symmetrized mc truth
  if( hasSymMCTruth() ) {
    pt3.clear();
    pt3e.clear();
    n.clear();
    ne.clear();
    pt3.push_back(0.);
    pt3e.push_back(0.);
//     n.push_back(symMCTruth_->fTailMCSmeared());
//     ne.push_back(symMCTruth_->fTailMCSmearedErr());
    n.push_back(symResp());
    ne.push_back(symRespErr());
    delete gFTailMCTruth_;
    gFTailMCTruth_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					   &(pt3e.front()),&(pt3e.front()),
					   &(ne.front()),&(ne.front()));
    gFTailMCTruth_->SetMarkerStyle(24);
    gFTailMCTruth_->SetMarkerColor(kGreen+2);
    gFTailMCTruth_->SetLineColor(gFTailMCTruth_->GetMarkerColor());

    // After Gaussian subtraction
    pt3.clear();
    pt3e.clear();
    n.clear();
    ne.clear();
    pt3.push_back(0.);
    pt3e.push_back(0.);
    n.push_back(symMCTruth_->fTailMCSmearedNonGauss());
    ne.push_back(symMCTruth_->fTailMCSmearedNonGaussErr());
    delete gFTailMCTruthNonGauss_;
    gFTailMCTruthNonGauss_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
						   &(pt3e.front()),&(pt3e.front()),
						   &(ne.front()),&(ne.front()));
    gFTailMCTruthNonGauss_->SetMarkerStyle(27);
    gFTailMCTruthNonGauss_->SetMarkerColor(kGreen);
    gFTailMCTruthNonGauss_->SetLineColor(gFTailMCTruth_->GetMarkerColor());

    //    std::cout << ">>>> " << (symMCTruth_->fTailMCSmearedNonGauss()/symMCTruth_->fTailMCSmeared()) << std::endl;
  }

  // Fill graphs of ftail of toy asymmetry from mc truth
  if( hasToyMC() ) {
    pt3.clear();
    pt3e.clear();
    n.clear();
    ne.clear();
    pt3.push_back(0.);
    pt3e.push_back(0.);
//     n.push_back(toyMC_->fTailMCSmeared());
//     ne.push_back(toyMC_->fTailMCSmearedErr());
    n.push_back(toyMC());
    ne.push_back(toyMCErr());
    delete gFTailToyMC_;
    gFTailToyMC_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					 &(pt3e.front()),&(pt3e.front()),
					 &(ne.front()),&(ne.front()));
    gFTailToyMC_->SetMarkerStyle(21);
    gFTailToyMC_->SetMarkerColor(kRed);
    gFTailToyMC_->SetLineColor(gFTailToyMC_->GetMarkerColor());
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
    double spreade = (*it)->fTailMCSmearedErr()/fExMC_->Eval(thres);
    pt3.push_back(thres);
    pt3e.push_back(0.);
    n.push_back(spread);
    ne.push_back(spreade);
    hSpread->Fill(spread);
  }
  delete gFTailSpreadMC_;
  gFTailSpreadMC_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					  &(pt3e.front()),&(pt3e.front()),
					  &(ne.front()),&(ne.front()));
  setStyleMC(gFTailSpreadMC_);

  // Extrapolation of ftail
  extraMC_ = useExtrapolatedValue ? fExMC_->GetParameter(0) : fExMC_->Eval(gFTailMC_->GetX()[0]);
  extraMCErr_ = fExMC_->GetParError(0);//hSpread->GetRMS()*extraMC_;


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
      double spreade = (*it)->fTailDataErr()/fExData_->Eval(thres);
      pt3.push_back(thres);
      pt3e.push_back(0.);
      n.push_back(spread);
      ne.push_back(spreade);
      hSpread->Fill(spread);
    }
  }
  delete gFTailSpreadData_;
  gFTailSpreadData_ = new TGraphAsymmErrors(pt3.size(),&(pt3.front()),&(n.front()),
					    &(pt3e.front()),&(pt3e.front()),
					    &(ne.front()),&(ne.front()));
  setStyleData(gFTailSpreadData_);

  // Extrapolated ftail
  extraData_ = useExtrapolatedValue ? fExData_->GetParameter(0) : fExData_->Eval(gFTailData_->GetX()[0]);
  extraDataErr_ = fExData_->GetParError(0);//hSpread->GetRMS()*extraData_;
  delete hSpread;


  // Set delta: absolute difference of extrapolated values
  deltaEx_ = extraData_ - extraMC_;
  // Uncertainty from fitted parameters
  deltaExErr_ = sqrt( pow(extraDataErr_,2.) + pow(extraMCErr_,2.) );


  // Scaling factor (fmc(0)+delta)/fmc(0)
  double ref = extraMC_;
  double refE = extraMCErr_;
  if( hasToyMC() && mcTruthRef ) {
    ref = toyMC_->fTailMCSmeared();
    refE = toyMC_->fTailMCSmearedErr();
  }
  
  //  std::cout << "(" << etaBin_ << "," << ptBin_ << ")   $" << ref << " \\pm " << refE << "$" << std::endl;

  scalingFactor_ = (ref+deltaEx_)/ref;
  scalingFactorErr_ = sqrt( pow(deltaExErr_/ref,2.) + pow(deltaEx_*refE/ref/ref,2.) );

  std::cout << "PtBin " << ptBin_ << ": " << deltaEx_ << " +/- " << deltaExErr_ << " --> " << scalingFactor_ << " +/- " << scalingFactorErr_ << std::endl;

  if( DEBUG ) std::cout << "Leaving EtaPtBin::extrapolate()" << std::endl;
}




// Constructor for asymmetry
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Thres, double nSigCore, double coreScalingFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel)
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(pt3Bin), pt3Thres_(pt3Thres),
    coreScalingFactor_(coreScalingFactor), 
    tailStartBin_(windowMinBin), tailEndBin_(windowMaxBin) {

  // Non-read distributions
  hResp_ = 0;
  hSymResp_ = 0;
  fGaussMCTruth_ = 0;

  // Get spectra
  TString histName = "hPtAveCombined_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  hPtAveSpecData_ = util::FileOps::readTH1(fileNameData,histName,"Data"+histName+"_Pt3Bin"+util::toTString(pt3Bin_));
  util::HistOps::setAxisTitles(hPtAveSpecData_,"p^{ave}_{T}","GeV","events");      
  hPtAveSpecData_->SetTitle(HEADER);
  setStyleDataMarker(hPtAveSpecData_);
  ptAveMeanData_ = hPtAveSpecData_->GetMean();
  ptAveMeanDataErr_ = hPtAveSpecData_->GetMeanError();
  hPtAveSpecMC_ = util::FileOps::readTH1(fileNameMC,histName,"MC"+histName+"_Pt3Bin"+util::toTString(pt3Bin_));
  util::HistOps::setAxisTitles(hPtAveSpecMC_,"p^{ave}_{T}","GeV","events");      
  hPtAveSpecMC_->SetTitle(HEADER);
  setStyleMCFilled(hPtAveSpecMC_);

  // Get asymmetry distributions
  hAsymData_ = readHist(fileNameData,"Data");
  hAsymMC_ = readHist(fileNameMC,"MC");
  hAsymData_->UseCurrentStyle();
  hAsymData_->SetTitle(HEADER);
  hAsymMC_->UseCurrentStyle();
  hAsymMC_->SetTitle(HEADER);
  

  // Smear MC asymmetry
  sig_ = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(hAsymMC_,nSigCore,sig_,widthErr) ) {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),sig_,coreScalingFactor);
  } else {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,hAsymMC_->GetEntries(),0.,0.);
  }

  // Width of smeared asymmetry
  if( !util::HistOps::fitCoreWidth(hAsymMCSmeared_,nSigCore,sigSmeared_,widthErr) ) sigSmeared_ = 0.;

  // Get relative number of entries in tail
  getFTail(hAsymMC_,hAsymMC_->GetEntries(),tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
  getFTail(hAsymMCSmeared_,hAsymMC_->GetEntries(),tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
  double minEff = hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_);
  double maxEff = min(1.,hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_));
  fNMCSmearedGauss_ = (erf(maxEff/sqrt(2.)/sigSmeared_) - erf(minEff/sqrt(2.)/sigSmeared_))/2.;

  getFTail(hAsymData_,hAsymData_->GetEntries(),tailStartBin_,tailEndBin_,fNData_,fNDataErr_);

  setStyleDataMarker(hAsymData_);
  setStyleMCFilled(hAsymMC_);
  setStyleMCFilled(hAsymMCSmeared_);
  initBinLabel(adm,jetLabel);
}


// Constructor for symmetrized MC truth response
//  fGaussMCTruth from symmetrized and smeared MC truth response
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, double windowRespMin, double windowRespMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel) 
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(1000), pt3Thres_(1000.),
    coreScalingFactor_(coreScalingFactor) {

  hAsymData_ = 0;
  hAsymMC_ = 0;
  hAsymMCSmeared_ = 0;
  fGaussMCTruth_ = 0;
  hResp_ = 0;
  hPtAveSpecData_ = 0;
  hPtAveSpecMC_ = 0;
  sig_ = 0.;
  sigSmeared_ = 0.;
  ptAveMeanData_ = 0.;
  ptAveMeanDataErr_ = 0.;
  fNMCSmearedGauss_ = 0.;

  
  TH1* h = readMCTruthResponse(fileName,"SymmetrizedMCTruth");
  h->Scale(1./h->Integral("width"));
  double entries = h->GetEntries();

  // Smear core
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore,width,widthErr) ) {
    func::smearHistogram(h,hSymResp_,h->GetEntries(),width,coreScalingFactor);
  } else {
    func::smearHistogram(h,hSymResp_,h->GetEntries(),0.,0.);
  }
  //  std::cout << "(" << etaBin_ << "," << ptBin_ << ") Entries " << h->GetEntries() << " (" << hSymResp_->GetEntries() << ")\n";

  // Window borders from unsmeared distribution!
  tailStartBin_ = h->FindBin(windowRespMin);
  tailEndBin_ = h->FindBin(windowRespMax);

  delete h;


  // Get relative number of entries in tail
  getFTail(hSymResp_,entries,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);

  //  std::cout << "  " << tailWindowMin << " (" << tailStartBinMC_ << ")" << " - " << tailWindowMax << " (" << tailEndBinMC_ << ")\n";
  //  std::cout << "  " << fNMCSmeared_ << std::endl;

  // Get tail after Gaussian subtraction (note: this is a Gaussian to the smeared histogram!)
  if( util::HistOps::fitCoreWidth(hSymResp_,nSigCore,fGaussMCTruth_,width,widthErr) ) {
    // Subtract Gaussian
    double winMin = hSymResp_->GetXaxis()->GetBinLowEdge(tailStartBin_);
    double winMax = hSymResp_->GetXaxis()->GetBinUpEdge(tailEndBin_);
    double gaussContr = fGaussMCTruth_->Integral(winMin,winMax);
    fNMCSmearedNonGauss_ = std::max(fNMCSmeared_-gaussContr,0.);
    fNMCSmearedNonGaussErr_ = fNMCSmearedErr_;
    // Style
    fGaussMCTruth_->SetLineColor(kRed);
    fGaussMCTruth_->SetLineWidth(1);
  }

  setStyleDataMarker(hSymResp_);
  initBinLabel(adm,jetLabel,true);
}


// Constructor for toy asymmetry from MC truth response
//  fAsymMCSmeared from toy MC (MC truth --> asymmetry --> smearing)
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScalingFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel) 
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(1000), pt3Thres_(1000.),
    coreScalingFactor_(coreScalingFactor),
    tailStartBin_(windowMinBin), tailEndBin_(windowMaxBin) {

  if( DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Entering constructor for toy asymmetry" << std::endl;

  hAsymData_ = 0;
  hAsymMC_ = 0;
  hAsymMCSmeared_ = 0;
  fGaussMCTruth_ = 0;
  hSymResp_ = 0;
  hPtAveSpecData_ = 0;
  hPtAveSpecMC_ = 0;
  sig_ = 0.;
  sigSmeared_ = 0.;
  ptAveMeanData_ = 0.;
  ptAveMeanDataErr_ = 0.;
  fNMCSmearedGauss_ = 0.;


  // Get response from file  
  hResp_ = readMCTruthResponse(fileName,"MCTruth");
  double entries = hResp_->GetEntries(); // Use entries of original distribution for statistical precision
  hResp_->UseCurrentStyle();
  hResp_->SetMarkerStyle(20);
  util::HistOps::setAxisTitles(hResp_,"Response","","jets");
  hResp_->SetTitle(HEADER);

  // Generate asymmetry distribution from response
  if( DEBUG ) std::cout << "  Generating toy asymmetry distribution  . . .  " << std::flush;
  hAsymMC_ = new TH1D("ToyAsym_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin),HEADER,nAsymBins,-1.,1.);
  util::HistOps::setAxisTitles(hAsymMC_,"Asymmetry","","events",true);
  hAsymMC_->Sumw2();
  for(int i = 0; i < 1000000; ++i) {
    double x1 = hResp_->GetRandom();
    double x2 = hResp_->GetRandom();
    double sum = x1+x2;
    if( sum > 0. ) {
      hAsymMC_->Fill((x1-x2)/sum);
      hAsymMC_->Fill((x2-x1)/sum);
    }
  }
  if( hAsymMC_->Integral() ) hAsymMC_->Scale(1./hAsymMC_->Integral("width"));
  if( DEBUG ) std::cout << "ok" << std::endl;
  
  // Smear toy asymmetry
  if( DEBUG ) std::cout << "  Smearing toy asymmetry distribution  . . .  " << std::flush;
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(hAsymMC_,nSigCore,width,widthErr) ) {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,entries,width,coreScalingFactor);
  } else {
    func::smearHistogram(hAsymMC_,hAsymMCSmeared_,entries,0.,0.);
  }
  if( DEBUG ) std::cout << "ok" << std::endl;

  if( DEBUG ) std::cout << "  Getting relative number of tail events in toy asymmetry distribution  . . .  " << std::flush;
  // Get relative number of entries in tail
  getFTail(hAsymMC_,entries,tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
  getFTail(hAsymMCSmeared_,entries,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
  if( DEBUG ) std::cout << "ok" << std::endl;

  // Set style
  //setStyleMCFilled(hAsymMC_);
  //setStyleMCFilled(hAsymMCSmeared_);
  initBinLabel(adm,jetLabel,true);

  if( DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Leaving constructor for toy asymmetry" << std::endl;
}


// ------------------------------------------------------------------------------------
void Pt3Bin::initBinLabel(const sampleTools::BinningAdmin* adm, const TString &jetLabel, bool isMCTruth) {
  if( isMCTruth ) binLabel_ = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
  else binLabel_ = util::LabelFactory::createPaveText(5,BINLABEL_WIDTH);
  binLabel_->AddText(LUMI_LABEL);
  binLabel_->AddText(jetLabel);
  binLabel_->AddText(util::toTString(adm->etaMin(etaBin_))+" < |#eta| < "+util::toTString(adm->etaMax(etaBin_)));
  binLabel_->AddText(util::toTString(adm->ptMin(etaBin_,ptBin_))+" < p^{ave}_{T} < "+util::toTString(adm->ptMax(etaBin_,ptBin_))+" GeV/c");
  if( !isMCTruth ) binLabel_->AddText(util::toTString(adm->ptSoftMin(pt3Bin_))+" < p^{rel}_{T,3} < "+util::toTString(adm->ptSoftMax(pt3Bin_)));
}


// ------------------------------------------------------------------------------------
Pt3Bin::~Pt3Bin() {
  if( hAsymData_ ) delete hAsymData_;
  if( hAsymMC_ ) delete hAsymMC_;
  if( hAsymMCSmeared_ ) delete hAsymMCSmeared_;
  if( hResp_ ) delete hResp_;
  if( hSymResp_ ) delete hSymResp_;
  if( fGaussMCTruth_ ) delete fGaussMCTruth_;
  if( binLabel_ ) delete binLabel_;
  if( hPtAveSpecData_ ) delete hPtAveSpecData_;
  if( hPtAveSpecMC_ ) delete hPtAveSpecMC_;
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


// Get MC truth response histogram from Kalibri input files
// ------------------------------------------------------------------------------------
TH1* Pt3Bin::readMCTruthResponse(const TString &fileName, const TString &type) const {

  TString histName = "";
  TString title;
  if( type == "SymmetrizedMCTruth") {
    //histName = "hRespSymAbs_"+util::toTString(ptBin_);
    histName = "hRespSymAbs";
    title = "Symmetrised Response";
  } else if( type == "MCTruth") {
    //histName = "hRespMeasAbs_"+util::toTString(ptBin_);
    histName = "hRespMeasAbs";
    title = "Response";
  } else {
    std::cerr << "Pt3Bin::readMCTruthResponse(): Unknown type '" << type << "'" << std::endl;
    exit(0);
  }
  histName += "_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,type+"_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_));
  h->GetXaxis()->SetRangeUser(0.,2.);
  h->SetTitle("");

  return h;
}


// Gaussian approximation for binomial error
// ------------------------------------------------------------------------------------
void Pt3Bin::getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr, double nExtraTailEvts) const {
  // Relative number of tail events
  // useage of Integral() because
  // a) possiblity to get #evts in certain interval
  // b) get total #evts consistently taking into account normalizations
  double nTotal = h->Integral();
  double nTail = h->Integral(start,end);// + nExtraTailEvts*nTotal/entries;
  fTail = nTail/nTotal;

  //if( nExtraTailEvts > 1. ) std::cout << "nTail = " << nTail << " (" << nTotal << ")" << std::endl;
  
  // Total number of events for uncertainty calculation
  // from number of entries (assume no overflow)
  nTotal = entries;
  fTailErr = sqrt( fTail*(1.-fTail)/nTotal );
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMC(const TString &outNameId) const {

  TLegend* leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hAsymData_,"Data","P");
  leg->AddEntry(hAsymMC_,"MC","F");

  // Log scale
  util::HistOps::setYRange(hAsymMC_,5,3E-5);
  util::HistOps::setYRange(hAsymData_,5,3E-5);
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
  util::HistOps::setYRange(hAsymMC_,5);
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
void Pt3Bin::plotAsymmetryDataMCSmeared(const TString &outNameId, double nSigTailStart, double nSigTailEnd) const {

  TLegend* leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hAsymData_,"Data","P");
  leg->AddEntry(hAsymMCSmeared_,"Reweighted MC","F");

  // Log scale
  util::HistOps::setYRange(hAsymMCSmeared_,5,3E-5);
  util::HistOps::setYRange(hAsymData_,5,3E-5);
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
//   TBox* win = new TBox(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
//  		       min(1.,hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_)),10.);
//   win->SetLineWidth(1);
//   win->SetFillStyle(3444);
//   win->SetLineColor(kRed);
//   win->SetFillColor(win->GetLineColor());

  TLine* win = new TLine(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
			 hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),10.);
  win->SetLineWidth(2);
  win->SetLineColor(kRed);

  TArrow* arr = new TArrow(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),3.,
			   hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_)+0.2,3.);
  arr->SetLineWidth(2);
  arr->SetLineColor(kRed);
  arr->SetArrowSize(0.05);
  arr->SetAngle(30);

  // Gaussian contribution
  TF1* gauss = new TF1("gauss","gaus",0.,1.);
  gauss->SetParameter(0,1./sqrt(2.*M_PI)/sigmaSmeared());
  gauss->SetParameter(1,0.);
  gauss->SetParameter(2,sigmaSmeared());
  gauss->SetLineWidth(2);
  gauss->SetLineColor(8);
  gauss->SetFillStyle(3004);
  gauss->SetFillColor(gauss->GetLineColor());
  gauss->SetRange(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),gauss->GetX(3E-5,0.,1.));
  
  TLegend* legWin = util::LabelFactory::createLegendCol(4,LEG_WIDTH);
  legWin->AddEntry(hAsymData_,"Data","P");
  legWin->AddEntry(hAsymMCSmeared_,"Reweighted MC","F");
  if( nSigTailEnd < 50. ) legWin->AddEntry(win,util::toTString(nSigTailStart)+" - "+util::toTString(nSigTailEnd)+" #sigma Window","F");
  else legWin->AddEntry(win,"> "+util::toTString(nSigTailStart)+" #sigma Window","L");
  legWin->AddEntry(gauss,"Asymmetry Gauss","F");

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  gauss->Draw("same");
  hAsymData_->Draw("PE1same");
  win->Draw("same");
  arr->Draw();			// Arrow not drawn if option "same" called!
  binLabel_->Draw("same");
  legWin->Draw("same");
  can->SetLogy(1);
  gPad->RedrawAxis();
  can->SaveAs(outNameId+"PtSmearAsymTail.eps","eps");
  delete legWin;
  delete arr;

  // Without data to illustrate window definition
  if( pt3Bin_ == 0 ) {
    legWin = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
    legWin->AddEntry(hAsymMCSmeared_,"Reweighted MC","F");
    legWin->AddEntry(win,"Window","F");

    can->cd();
    hAsymMCSmeared_->Draw("HISTE");
    win->Draw("same");
    binLabel_->Draw("same");
    legWin->Draw("same");
    can->SetLogy(1);
    can->SaveAs(outNameId+"WindowDef.eps","eps");
    delete legWin;
  }
  delete win;

  // Linear scale
  util::HistOps::setYRange(hAsymMCSmeared_,5);
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
void Pt3Bin::plotSpectra(const TString &outNameId) const {

  TString canName = outNameId+"_PtAveSpectrum";
  TCanvas *can = new TCanvas(canName,canName,500,500);

  // Populated x-region
  int xMin = 1;
  int xMax = 1000;
  util::HistOps::findXRange(hPtAveSpecMC_,xMin,xMax);
  xMin = max(1,xMin-10);
  xMax = min(xMax+10,hPtAveSpecMC_->GetNbinsX());

  // Absolute data spectrum
  TLegend* leg = util::LabelFactory::createLegendCol(1,LEG_WIDTH);
  leg->AddEntry(hPtAveSpecData_,"Data","P");
  hPtAveSpecData_->GetXaxis()->SetRange(xMin,xMax);
  util::HistOps::setYRange(hPtAveSpecData_,6,3E-1);
  can->cd();
  hPtAveSpecData_->Draw("PE1");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtAveSpectrumData.eps","eps");

  // Absolute MC spectrum
  delete leg;
  leg = util::LabelFactory::createLegendCol(1,LEG_WIDTH);
  leg->AddEntry(hPtAveSpecMC_,"MC","F");
  util::HistOps::setYRange(hPtAveSpecMC_,6,3E-5);
  can->cd();
  hPtAveSpecMC_->GetXaxis()->SetRange(xMin,xMax);
  hPtAveSpecMC_->Draw("HIST");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"PtAveSpectrumMC.eps","eps");

  // Comparison (normalised distributions)
  delete leg;
  leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hPtAveSpecData_,"Data","P");
  leg->AddEntry(hPtAveSpecMC_,"MC","F");
  double scaleData = hPtAveSpecData_->Integral("width");
  double scaleMC = hPtAveSpecMC_->Integral("width");
  if( scaleData > 0. && scaleMC > 0. ) {
    hPtAveSpecData_->Scale(1./scaleData);
    hPtAveSpecMC_->Scale(1./scaleMC);
    hPtAveSpecMC_->GetYaxis()->SetRangeUser(3E-8,90.);
    can->cd();
    hPtAveSpecMC_->Draw("HISTE");
    hPtAveSpecData_->Draw("PE1same");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    can->SaveAs(outNameId+"PtAveSpectra.eps","eps");
    hPtAveSpecData_->Scale(scaleData);
    hPtAveSpecMC_->Scale(scaleMC);
  }  
  delete leg;
  delete can;
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotSymMCTruthResponse(const TString &outNameId) const {
  if( DEBUG ) std::cout << "Entering Pt3Bin::plotSymMCTruthResponse()" << std::endl;

  if( DEBUG ) std::cout << "  Creating window and legend  . . .  " << std::flush;
  // Log scale
  TBox* win = new TBox(hSymResp_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
		       hSymResp_->GetXaxis()->GetBinUpEdge(tailEndBin_),10.);
  win->SetLineWidth(1);
  win->SetFillStyle(3444);
  win->SetLineColor(kRed);
  win->SetFillColor(win->GetLineColor());

  TLegend* leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hSymResp_,"Reweighted MC","P");
  leg->AddEntry(win,"Window","F");
  if( DEBUG ) std::cout << "ok" << std::endl;

  if( DEBUG ) std::cout << "  Creating plots  . . .  " << std::flush;
  util::HistOps::setYRange(hSymResp_,4,3E-5);
  hSymResp_->GetXaxis()->SetRangeUser(0.,2.);

  TString canName = outNameId+"_SymMCTruthResponse";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hSymResp_->Draw("HIST");
  //  fGaussMCTruth_->Draw("same");
  win->Draw("same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"_SymMCTruthResponseLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hSymResp_,4);
  hSymResp_->GetXaxis()->SetRangeUser(0.4,1.6);
  can->cd();
  hSymResp_->Draw("HIST");
  //  fGaussMCTruth_->Draw("same");
  binLabel_->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"_SymMCTruthResponse.eps","eps");

  hSymResp_->GetXaxis()->SetRangeUser(0.,2.);

  if( DEBUG ) std::cout << "ok" << std::endl;

  delete win;
  delete leg;
  delete can;

  if( DEBUG ) std::cout << "Leaving Pt3Bin::plotSymMCTruthResponse()" << std::endl;
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotToyAsymmetry(const TString &outNameId) const {
  if( DEBUG ) std::cout << "Entering Pt3Bin::plotToyAsymmetry()" << std::endl;

  // Log scale
  TBox* win = new TBox(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
		       hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_),10.);
  win->SetLineWidth(1);
  win->SetFillStyle(3444);
  win->SetLineColor(kRed);
  win->SetFillColor(win->GetLineColor());

  TLegend* leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hAsymMCSmeared_,"Reweighted Toy MC","L");
//   leg->AddEntry(winMC,"Window","F");

  util::HistOps::setYRange(hAsymMCSmeared_,4,3E-5);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);

  TString canName = outNameId+"_ToyMC";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hAsymMCSmeared_->Draw("HIST");
  win->Draw("same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"_ToyMCTail.eps","eps");

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);
  can->cd();
  hAsymMCSmeared_->Draw("HIST");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"_ToyMCLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hAsymMCSmeared_,4);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hAsymMCSmeared_->Draw("HIST");
  binLabel_->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"_ToyMC.eps","eps");

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);

  // MC truth response
  delete leg;
  leg = util::LabelFactory::createLegendCol(1,LEG_WIDTH);
  leg->AddEntry(hAsymMCSmeared_,"MC Truth","P");

  util::HistOps::setYRange(hResp_,4,3E-5);
  hResp_->GetXaxis()->SetRangeUser(0.,2.);

  canName = outNameId+"_MCTruthResponse";
  can->SetName(canName);
  can->cd();
  hResp_->Draw("HIST");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  can->SaveAs(outNameId+"_MCTruthResponseLog.eps","eps");

  // Linear scale
  util::HistOps::setYRange(hResp_,4);
  hResp_->GetXaxis()->SetRangeUser(0.4,1.6);
  can->cd();
  hResp_->Draw("HIST");
  binLabel_->Draw("same");
  can->SetLogy(0);
  can->SaveAs(outNameId+"_MCTruthResponse.eps","eps");

  hResp_->GetXaxis()->SetRangeUser(0.,2.);


  delete win;
  delete leg;
  delete can;

  if( DEBUG ) std::cout << "Leaving Pt3Bin::plotToyAsymmetry()" << std::endl;
}





////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////

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


// ------------------------------------------------------------------------------------
void printBin(unsigned int etaBin, unsigned int ptBin, const sampleTools::BinningAdmin* adm) {
  cout.setf(ios::fixed,ios::floatfield);
  std::cout << setprecision(1) << "    $" << adm->etaMin(etaBin) << " - " << adm->etaMax(etaBin) << "$ & $";
  std::cout << setprecision(0) << adm->ptMin(etaBin,ptBin) << " - " << adm->ptMax(etaBin,ptBin) << "$";
  std::cout << setprecision(5);
}


// ------------------------------------------------------------------------------------
void printWindowBorders(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm) {
  std::cout << "\n\n\n***  Window Borders  ***\n\n";

  std::cout << "\\begin{tabular}[ht]{cccccc}\n\\hline\n";
  std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $A_{0} - A_{1}$ & $A_{0}/\\sigma - A_{1}/\\sigma$ & $\\int^{A_{1}}_{A_{0}}\\mathcal{G}$ & $\\fasymmc(\\pti{3} = 0.05)$ \\\\ \n\\hline\n";
  for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
    EtaPtBin* bin = *it;
    double min = bin->tailWindowMin();
    double max = bin->tailWindowMax();
    double sig = bin->sigma(0);
    double sigSmear = bin->sigmaSmeared(0);
    double fGauss = (erf(max/sqrt(2.)/sigSmear) - erf(min/sqrt(2.)/sigSmear))/2.;
    double fAsym = bin->fTailMCSmeared(0);
    double fAsymErr = bin->fTailMCSmearedErr(0);

    printBin(bin->etaBin(),bin->ptBin(),adm);
    std::cout << std::setprecision(3) << " & $" << min << " - " << max << "$ & $";
    std::cout << min/sigSmear << " - " << max/sigSmear << "$ & $";
    std::cout << std::setprecision(4) << fGauss << "$ & $" << fAsym << " \\pm " << fAsymErr << "$ \\\\" << std::endl;
    if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\hline\n";
  }
  std::cout << "  \\end{tabular}\n\n";

//   // hack
//   for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
//     EtaPtBin* bin = *it;
//     std::cout << "min.push_back(" << bin->tailWindowMin() << ");" << std::endl;
//   }
//   for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
//     EtaPtBin* bin = *it;
//     std::cout << "max.push_back(" << bin->tailWindowMax() << ");" << std::endl;
//   }
}


// ------------------------------------------------------------------------------------
void printMCClosure(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm) {
  if( bins.front()->hasToyMC() || bins.front()->hasSymMCTruth()  ) {
    std::cout << "\n\n\n***  MC Closure  ***\n";
    std::cout << "\\begin{tabular}[ht]{ccccc}\n\\hline\n";
    std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\fasymmc(0)$ & \\fasymtoy & $\\fasymmc(0))/\\fasymtoy$ \\\\ \n\\hline \n";
    for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
      EtaPtBin* bin = *it;
      printBin(bin->etaBin(),bin->ptBin(),adm);
      std::cout << " & $" << std::setprecision(3) << bin->extraMC() << " \\pm " << bin->extraMCErr();
      std::cout << "$ & $" << bin->toyMC() << " \\pm " << bin->toyMCErr();
      std::cout << "$ & $" << bin->extraMC()/bin->toyMC() << "$ \\\\ \n";
      if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\hline\n";
    }
    std::cout << "  \\end{tabular}\n\n";
  }
}


// ------------------------------------------------------------------------------------
void printExtrapolation(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm) {
  std::cout << "\n\n\n***  Extrapolation  ***\n";
  std::cout << "\\begin{tabular}[ht]{ccccc}\n\\hline\n";
  std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\mean{\\ptave} \\,(\\gevnospace)$ & $\\fasymdata(0)$ & $\\fasymmc(0)$ \\\\ \n\\hline \n";
  for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
    EtaPtBin* bin = *it;
    printBin(bin->etaBin(),bin->ptBin(),adm);
    std::cout << setprecision(1);
    std::cout << " & $"  << bin->ptAveMeanData(0) << " \\pm " << bin->ptAveMeanDataErr(0);
    std::cout << setprecision(3);
    std::cout << "$ & $"  << bin->extraData() << " \\pm " << bin->extraDataErr();  
    std::cout << "$ & $" << bin->extraMC() << " \\pm " << bin->extraMCErr() << "$ \\\\ \n";
    if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\hline\n";
  }
  std::cout << "  \\end{tabular}\n\n";
}






////////////////////// COMBINE TO FINAL RESULTS //////////////////////////////////


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* nomRatio(const TH1* hNom, const TH1* hPtMean) {

  std::vector<double> x;
  std::vector<double> xe;
  std::vector<double> y;
  std::vector<double> ye;
  for(int bin = 1; bin <= hNom->GetNbinsX(); ++bin) {
    x.push_back(hPtMean->GetBinContent(bin));
    xe.push_back(hPtMean->GetBinError(bin));
    y.push_back(hNom->GetBinContent(bin));
    ye.push_back(hNom->GetBinError(bin));
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
					       &(xe.front()),&(xe.front()),&(ye.front()),&(ye.front()));
  g->SetMarkerStyle(20);

  return g;
}


// ------------------------------------------------------------------------------------
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
  g->SetLineColor(color);

  return g;
}


// ------------------------------------------------------------------------------------
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
  g->SetFillStyle(0);
  g->SetFillColor(0);
  g->SetLineColor(color);

  return g;
}


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* uncertaintyBand(const TH1* hNom, const TGraphAsymmErrors* gRelUncert, int color) {
  TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(gRelUncert->Clone());
  for(int i = 0; i < g->GetN(); ++i) {
    double y = hNom->GetBinContent(1+i);
    g->SetPoint(i,g->GetX()[i],y);
    double eu = (g->GetEYhigh()[i])*y;
    double ed = (g->GetEYlow()[i])*y;
    g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],ed,eu);
  }
  g->SetFillStyle(1001);
  g->SetFillColor(color);
  g->SetLineColor(color);

  return g;
}
  

// ------------------------------------------------------------------------------------
void plotFinalResult(int windowIdx = 0) {
  if( SHOW_HEADER ) {
    util::StyleSettings::paper();
//     gStyle->SetTitleAlign(13);
//     gStyle->SetTitleFontSize(0.1);
//     gStyle->SetTitleX(0.7);
//     gStyle->SetTitleH(0.038);
  } else {
    util::StyleSettings::paperNoTitle();
  }
  gErrorIgnoreLevel = 1001;

  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("../resolutionFit/config/Analysis2010/Binning2010/BinningAdminTailsRebinned.cfg");  

  TString fileNamePrefix = "results/Analysis2010/RA2CleaningV1/TailSig15-25";
  TString outNamePrefix = "TailScalingFactorsSig15-25_PF";
  if( windowIdx == 1 ) {
    fileNamePrefix = "results/Analysis2010/RA2CleaningV1/TailSig25-50";
    outNamePrefix = "TailScalingFactorsSig25-50_PF";
  } else if( windowIdx == 2 ) {
    fileNamePrefix = "results/Analysis2010/RA2CleaningV1/TailSig15-Inf";
    outNamePrefix = "TailScalingFactorsSig15-Inf_PF";
  } else if( windowIdx == 3 ) {
    fileNamePrefix = "results/Analysis2010/RA2CleaningV1/TailSig25-Inf";
    outNamePrefix = "TailScalingFactorsSig25-Inf_PF";
  }

  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    // ****** Nominal scaling factors *****************************************************

    // Read mean pt
    //TString histName = "hPtAveMeanData_Eta"+util::toTString(etaBin);
    TString histName = "hPtAveMeanData_Eta"+util::toTString(etaBin);
    TH1* hPtAveMeanData = util::FileOps::readTH1(fileNamePrefix+"_PF.root",histName,histName);

    // Read nominal scaling factors
    histName = "hScaleFrame_Eta"+util::toTString(etaBin);
    TH1* hScaleFrame = util::FileOps::readTH1(fileNamePrefix+"_PF.root",histName,histName);
    hScaleFrame->GetXaxis()->SetTitle("p^{ave}_{T} (GeV/c)");
    histName = "hScale_Eta"+util::toTString(etaBin);
    TH1* hScaleNom = util::FileOps::readTH1(fileNamePrefix+"_PF.root",histName,histName+"_Nominal");

    // Graph of nominal scaling factors with ptAveMean as x
    TGraphAsymmErrors* gNomRatio = nomRatio(hScaleNom,hPtAveMeanData);



    // ****** Variations for systematic uncertainties *************************************

    // Read histograms with variations
    TH1* hScaleVarExtra = util::FileOps::readTH1(fileNamePrefix+"VarExtrapolation_PF.root",histName,histName+"_VarExtrapolation");
    TH1* hScaleVarCoreUp = util::FileOps::readTH1(fileNamePrefix+"VarCoreUp_PF.root",histName,histName+"_VarCoreUp");
    TH1* hScaleVarCoreDown = util::FileOps::readTH1(fileNamePrefix+"VarCoreDown_PF.root",histName,histName+"_VarCoreDown");
    TH1* hScaleVarClosure = util::FileOps::readTH1(fileNamePrefix+"VarClosure_PF.root",histName,histName+"_VarClosure");

    // Get relative uncertainties
    std::vector<TGraphAsymmErrors*> uncerts;
    uncerts.push_back(relUncertainty(hScaleNom,38,1.,hScaleVarExtra));
    uncerts.push_back(relUncertainty(hScaleNom,46,1.,hScaleVarCoreUp,hScaleVarCoreDown));
    uncerts.push_back(relUncertainty(hScaleNom,8,0.5,hScaleVarClosure));

    // Define labels
    std::vector<TString> uncertLabels;
    uncertLabels.push_back("Extrapolation");
    uncertLabels.push_back("Reweighting");
    uncertLabels.push_back("MC Closure");

    // Add up uncertainties
    TGraphAsymmErrors* gUncertRelTotal = totalUncertainty(uncerts,kBlack);
    TGraphAsymmErrors* gUncertAbs = uncertaintyBand(hScaleNom,gUncertRelTotal,5);
    TGraphAsymmErrors* gUncertAbsPos = static_cast<TGraphAsymmErrors*>(gUncertAbs->Clone());
    for(int i = 0; i < gUncertAbsPos->GetN(); ++i) {
      gUncertAbsPos->SetPointError(i,gUncertAbsPos->GetEXlow()[i],gUncertAbsPos->GetEXhigh()[i],
				   gUncertAbsPos->GetEYlow()[i]>gUncertAbsPos->GetY()[i] ? gUncertAbsPos->GetY()[i] : gUncertAbsPos->GetEYlow()[i],gUncertAbsPos->GetEYhigh()[i]);
    }

    // ****** Plotting ********************************************************************

    // Label
    TPaveText* label = util::LabelFactory::createPaveText(3,BINLABEL_WIDTH);
    //label->AddText(LUMI_LABEL);
    label->AddText("Anti-k_{T} (R=0.5) PF Jets");
    label->AddText(util::toTString(binAdm->etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdm->etaMax(etaBin)));
    if( windowIdx == 0 ) {
      label->AddText("1.5 - 2.5 #sigma Window");
    } else if( windowIdx == 1 ) {
      label->AddText("2.5 - 5.0 #sigma Window");
    } else if( windowIdx == 2 ) {
      label->AddText("> 1.5 #sigma Window");
    } else if( windowIdx == 3 ) {
      label->AddText("> 2.5 #sigma Window");
    }

    TLegend* legScale = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
    legScale->AddEntry(gNomRatio,"Scaling Factor","P");
    legScale->AddEntry(gUncertAbs,"Syst. Uncertainty","F");

    TLegend* legUncert = util::LabelFactory::createLegendCol(uncerts.size()+1,LEG_WIDTH);
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      legUncert->AddEntry(uncerts.at(i),uncertLabels.at(i),"F");
    }
    legUncert->AddEntry(gUncertRelTotal,"Total","F");


    // Plot scaling factors and total uncertainty    
    hScaleFrame->SetTitle(HEADER);
    hScaleFrame->GetYaxis()->SetRangeUser(0.01,2.99);
    if( gNomRatio->GetY()[0] > 2. ) hScaleFrame->GetYaxis()->SetRangeUser(0.01,6.99);
    hScaleFrame->GetYaxis()->SetTitle("Scaling Factor");
    hScaleFrame->GetXaxis()->SetMoreLogLabels();
    hScaleFrame->GetXaxis()->SetNoExponent();
    TCanvas* canScale = new TCanvas("canScale"+util::toTString(etaBin),"Scaling Factors Eta "+util::toTString(etaBin),500,500);
    canScale->cd();
    hScaleFrame->Draw("HIST");
    if( windowIdx == 0 && etaBin == 3 ) {
      gUncertAbsPos->Draw("E2same");
    } else {
      gUncertAbs->Draw("E2same");
    }
    hScaleFrame->Draw("HISTsame");
    gNomRatio->Draw("PE1same");
    label->Draw("same");
    legScale->Draw("same");
    canScale->SetLogx();
    gPad->RedrawAxis();
    canScale->SaveAs(outNamePrefix+"_Eta"+util::toTString(etaBin)+".eps","eps");
    canScale->SaveAs(outNamePrefix+"_Eta"+util::toTString(etaBin)+".root");


    // Plot relative uncertainties
    TH1* hUncertsFrame = static_cast<TH1D*>(hScaleFrame->Clone("hUncertsFrame"+util::toTString(etaBin)));
    for(int bin = 1; bin <= hUncertsFrame->GetNbinsX(); ++bin) {
      hUncertsFrame->SetBinContent(bin,0.);
    }
    hUncertsFrame->GetYaxis()->SetRangeUser(-0.99,1.49);
    hUncertsFrame->GetXaxis()->SetMoreLogLabels();
    hUncertsFrame->GetXaxis()->SetNoExponent();
    if( windowIdx == 0 && etaBin == 3 ) hUncertsFrame->GetYaxis()->SetRangeUser(-5.1,5.);
    hUncertsFrame->GetYaxis()->SetTitle("Relative Uncertainties");
    TCanvas* canRelUncerts = new TCanvas("canRelUncerts"+util::toTString(etaBin),"Relative Uncertainties Eta "+util::toTString(etaBin),500,500);
    canRelUncerts->cd();
    hUncertsFrame->Draw("HIST");
    gUncertRelTotal->Draw("E2same");
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      uncerts.at(i)->Draw("E2same");
    }
    hUncertsFrame->Draw("HISTsame");
    label->Draw("same");
    legUncert->Draw("same");
    canRelUncerts->SetLogx();
    gPad->RedrawAxis();
    canRelUncerts->SaveAs(outNamePrefix+"Uncertainties_Eta"+util::toTString(etaBin)+".eps","eps");




    // ****** Print results ***************************************************************

    // Print factors with total uncertainty (stat + syst)
    if( etaBin == 0 ) {
      std::cout << "\\begin{tabular}{ccccc}\n\\hline\n";
      std::cout << "$|\\eta|$ & $\\ptave\\,(\\gevnospace)$ & $\\mean{\\ptave}\\,(\\gevnospace)$ & Scaling Factor ($\\pm(\\text{stat}){}^{+\\text{syst}}_{-\\text{syst}}$) & Scaling Factor ($\\pm(\\text{total})$) \\\\\n";
      std::cout << "\\hline\n";
    }
    for(int bin = 1; bin <= hScaleNom->GetNbinsX(); ++bin) {
      cout.setf(ios::fixed,ios::floatfield);
      std::cout << std::setprecision(1) << "  $" << util::toTString(binAdm->etaMin(etaBin)) << " - " << util::toTString(binAdm->etaMax(etaBin)) << "$ & $";
      std::cout << std::setprecision(0) << util::toTString(binAdm->ptMin(etaBin,bin-1)) << " - " << util::toTString(binAdm->ptMax(etaBin,bin-1)) << "$ & $";
      std::cout << std::setprecision(1) << hPtAveMeanData->GetBinContent(bin) << " \\pm " << hPtAveMeanData->GetBinError(bin) << "$ & $";
      std::cout << std::setprecision(3) << hScaleNom->GetBinContent(bin);
      double estat = hScaleNom->GetBinError(bin);
      double esystd = gUncertAbs->GetEYlow()[bin-1];
      double esystu = gUncertAbs->GetEYhigh()[bin-1];
      double etotd = sqrt( estat*estat + esystd*esystd );
      double etotu = sqrt( estat*estat + esystu*esystu );
      std::cout << " \\pm " << estat << "^{+" << esystu << "}_{-" << esystd << "}$ & $";
      std::cout << hScaleNom->GetBinContent(bin) << "^{+" << etotu << "}_{-" << etotd << "} $ \\\\\n";
    }    
    std::cout << "\\hline\n";
    if( etaBin == 3 ) std::cout << "\\end{tabular}\n";
  }

  delete binAdm;
}



