// $Id: getTailScalingFactors.C,v 1.21 2012/06/15 23:08:26 mschrode Exp $

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
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
#include "TROOT.h"
#include "TString.h"
#include "TText.h"

#include "../sampleTools/BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"
#include "../util/MultiCanvas.h"


const bool DEBUG = false;
const bool SHOW_HEADER = false;
const double lumi = 4900.;

const TString LUMI = util::StyleSettings::luminosity(lumi);
const double BINLABEL_WIDTH = -0.52;
const double LEG_WIDTH = 0.48;
const double PT3PLOTMAX = 0.23;
const TString FASYM = "f_{asym}";
const TString FASYMMC = "f^{mc}_{asym}";
const TString FASYMDATA = "f^{data}_{asym}";
const TString FASYMGAUSS = "f^{gauss}_{asym}";
const TString FASYMTOY = "f^{toy}_{asym}";
const TString SCALE = "#rho_{tail}";
const TString ASYM_TAIL_START = "A_{tail}";
const TString ASYM_TAIL_START_EFF = "#hat{A}_{tail}";
const TString COMMON_SIGMA = "#sigma_{c}";
const TString DATA = util::LabelFactory::data(LUMI);
const TString MC = util::LabelFactory::mc();
const TString MCSMEAR = MC+" (Corrected #sigma_{A})^{#color[10]{A}}";
const TString MCGAUSS = "Gaussian Asymmetry";
const TString LABEL_NVTX = "N_{Vtx}";

const TString LUMI_LABEL = SHOW_HEADER ? "CMS preliminary, L = "+LUMI+",  #sqrt{s} = 7 TeV" : "#sqrt{s} = 7 TeV,  L = "+LUMI;
const TString HEADER = SHOW_HEADER ? LUMI_LABEL : "";

const int COLOR_GAUSS = 46;
const int COLOR_FILLED_ASYM = 38;
const int COLOR_FILLED_ASYM_SMEAR = 29;
const int COLOR_LINE_ASYM_SMEAR = 30;
const double MARKER_SIZE = 1.4;
const int LINE_WIDTH = 2;
const int HATCH_STYLE = 3354;



///////////////////////// TYPE DEFINITIONS /////////////////////////////////////////////////


const bool ROOT_OUTPUT = true;
bool EPS_OUTPUT = true;
TFile* ROOT_OUT_FILE = 0;

class Pt3Bin;
class EtaPtBin;
typedef std::vector<EtaPtBin*>::const_iterator EtaPtBinConstIt;

double coreScaleFactor(unsigned int etaBin, int variation);
void setStyleDataMarker(TH1* h);
void setStyleMCFilled(TH1* h);
void setStyleData(TGraphAsymmErrors* g);
void setStyleMC(TGraphAsymmErrors* g);

void printBin(unsigned int etaBin, unsigned int ptBin, const sampleTools::BinningAdmin* adm);
void printWindowBorders(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm, double nSigTailStart);
void printMCClosure(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm);
void printExtrapolation(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm);

bool completePicName(int padIdx, const sampleTools::BinningAdmin* adm, TString &picName);
void writeLaTeXSlides(const TString &outNameId, const sampleTools::BinningAdmin* adm);
TString labelWindow(double nSigMin, double nSigMax);
double getTailStart(const TString &name);
void toFiles(TNamed* obj, const TString &name = "");
void toFiles(TCanvas* obj, const TString &name = "");
void correctAsymmetryWidth(TH1* hOrig, double nSigCore, double coreScale, TH1* &hSmeared, double &width, double &smearedWidth);



// ------------------------------------------------------------------------------------
class Pt3Bin {
public:
  Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Threshold, double nSigCore, double coreScaleFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm);
  Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm);
 ~Pt3Bin();

  double pt3Thres() const { return pt3Thres_; }
  int nPtAsymBins() const { return hAsymMC_->GetNbinsX(); }

  double ptAveMeanData() const { return ptAveMeanData_; }
  double ptAveMeanDataErr() const { return ptAveMeanDataErr_; }

  double coreScalingFactor() const { return coreScalingFactor_; }
  double sigma() const { return sig_; }
  double sigmaSmeared() const { return sigSmeared_; }
  int tailStartBin() const { return tailStartBin_; }
  int tailEndBin() const { return tailEndBin_; }

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

  TH1* readHistAsym(const TString &fileName, const TString &id) const;  
  TH1* readMCTruthResponse(const TString &fileName, const TString &type) const;
  //  void init(double nSigCore, double coreScalingFactor, double tailWindowDataMin, double tailWindowDataMax, double tailWindowMCMin, double tailWindowMCMax);
  void initBinLabel(const sampleTools::BinningAdmin* adm, bool isMCTruth = false);
  void getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr) const;
};



// ------------------------------------------------------------------------------------
class EtaPtBin {
public:
  EtaPtBin(unsigned int etaBinIdx, unsigned int ptBinIdx, double nSigCore, double coreScaleFactor,const sampleTools::BinningAdmin* adm);
  ~EtaPtBin();


  unsigned int nPt3Bins() const { return pt3Bins_.size(); }
  unsigned int etaBin() const { return etaBin_; }
  unsigned int ptBin() const { return ptBin_; }

  double ptAveMeanData(int i) const { return pt3Bins_.at(i)->ptAveMeanData(); }
  double ptAveMeanDataErr(int i) const { return pt3Bins_.at(i)->ptAveMeanDataErr(); }

  double tailWindowMin() const { return tailWindowMin_; }
  double tailWindowMax() const { return tailWindowMax_; }
  double tailWindowMinEff() const { return tailWindowMinEff_; }
  double tailWindowMaxEff() const { return tailWindowMaxEff_; }
  double sigma(int i) const { return pt3Bins_.at(i)->sigma(); }
  double sigmaSmeared(int i) const { return pt3Bins_.at(i)->sigmaSmeared(); }

  double fTailData(int i) const { return pt3Bins_.at(i)->fTailData(); }
  double fTailDataErr(int i) const { return pt3Bins_.at(i)->fTailDataErr(); }
  double fTailMCSmeared(int i) const { return pt3Bins_.at(i)->fTailMCSmeared(); }
  double fTailMCSmearedErr(int i) const { return pt3Bins_.at(i)->fTailMCSmearedErr(); }
  double fTailMCSmearedGauss(int i) const { return pt3Bins_.at(i)->fTailMCSmearedGauss(); }
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
  void plotTailStart() const;

  bool hasToyMC() const { return hasToyMC_; }
  bool hasSymMCTruth() const { return hasSymMCTruth_; }

  void addPt3Bin(unsigned int pt3Bin, double thres, const TString &fileNameData, const TString &fileNameMC);
  void addMCTruthForToyAsym(const TString &fileName);
  void findWindow(const TString &fileName, unsigned int pt3Bin, double nSigTailWindowMin, double nSigTailWindowMax);
  void setWindow(const TString &fileName, double nSigTailWindowMin, double nSigTailWindowMax);
  void extrapolate(double minPt3Data, bool fixDataShape, bool useExtrapolatedValue, bool mcTruthRef);


private:
  const unsigned int etaBin_;
  const unsigned int ptBin_;
  const double nSigCore_;
  const double coreScalingFactor_;
  const double exMin_;
  const double exMax_;
  const sampleTools::BinningAdmin* binAdm_;

  double tailWindowMinEff_;
  double tailWindowMaxEff_;
  double tailWindowMin_;
  double tailWindowMax_;
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





//////////////////////////////// MAIN ROUTINES ///////////////////////////////////////////



// ------------------------------------------------------------------------------------
void getTailScalingFactors(double  nSigTailStart,
			   double  nSigTailEnd = 10000.,
			   bool    variationCoreUp = false,
			   bool    variationCoreDn = false,
			   bool    variationExtrapolation = false,
			   bool    variationClosure = false,
			   bool    variationPUDn = false,
			   bool    variationPUUp = false            ) {

  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting up parameters" << std::endl;

  // User specified parameters
  const TString uid                    = "163337-180252";

  const double  minPt3Data   = 0.;
  const bool    fixDataShape = false;
  const double  nSigCore     = 2.;
  const TString jetAlgo      = "PF";
  const bool    archivePlots = false;


  // +++++ Input Files +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const TString srcData = "../results/Analysis2011/Run2011_163337-180252_V10/";
  const TString srcMC   = "../results/Analysis2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/";
  
  // Nominal analysis
  TString fileNameData    = srcData+"ResTails_PtAveBins_Data2011_PF_L1FastJet_V10_REBINNED.root";
  TString fileNameMC      = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";
  if( variationPUDn )
    fileNameMC            = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_PUDn_REBINNED.root";
  else if( variationPUUp )
    fileNameMC            = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_PUUp_REBINNED.root";
  TString fileNameMCTruth = srcMC  +"ResTails_PtGenAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";
  TString fileNameTailWindow = "";

//    // PileUp selection: NVtx 0 - 6
//    fileNameTailWindow = "ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_Sig30-Inf_PF.root";
//    fileNameData            = srcData+"ResTails_PtAveBins_Data2011_PF_L1FastJet_V10_NVtx00-06_REBINNED.root";
//    fileNameMC              = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtx00-06_REBINNED.root";
//    if( variationPUDn )
//      fileNameMC            = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtx00-06_PUDn_REBINNED.root";
//    else if( variationPUUp )
//      fileNameMC            = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtx00-06_PUUp_REBINNED.root";

//      // PileUp selection: NVtx > 6
//      fileNameData            = srcData+"ResTails_PtAveBins_Data2011_PF_L1FastJet_V10_NVtx07-99_REBINNED.root";
//      fileNameMC              = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtx07-99_REBINNED.root";
//      if( variationPUDn )
//        fileNameMC            = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtx07-99_PUDn_REBINNED.root";
//      else if( variationPUUp )
//        fileNameMC            = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtx07-99_PUUp_REBINNED.root";


  // +++++ Setup +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // Sanity checks
  assert( nSigTailStart < nSigTailEnd );
  if( variationCoreUp )
    assert( !variationCoreDn && !variationExtrapolation && !variationClosure && !variationPUDn && !variationPUUp );
  else if( variationCoreDn )
    assert( !variationCoreUp && !variationExtrapolation && !variationClosure && !variationPUDn && !variationPUUp );
  else if( variationExtrapolation )
    assert( !variationCoreUp && !variationCoreDn && !variationClosure && !variationPUDn && !variationPUUp );
  else if( variationClosure )
    assert( !variationCoreUp && !variationCoreDn && !variationExtrapolation && !variationPUDn && !variationPUUp );
  else if( variationPUDn )
    assert( !variationCoreUp && !variationCoreDn && !variationExtrapolation && !variationClosure && !variationPUUp );
  else if( variationPUUp )
    assert( !variationCoreUp && !variationCoreDn && !variationExtrapolation && !variationPUDn && !variationClosure );

  if( variationCoreDn ) 
    std::cout << "\n***** Variation core down ***********************************\n" << std::endl;
  else if( variationCoreUp )
    std::cout << "\n***** Variation core up *************************************\n" << std::endl;
  else if( variationExtrapolation )
    std::cout << "\n***** Variation extrapolation *******************************\n" << std::endl;
  else if( variationClosure )
    std::cout << "\n***** Variation closure *************************************\n" << std::endl;
  else if( variationPUUp )
    std::cout << "\n***** Variation pile-up up **********************************\n" << std::endl;
  else if( variationPUDn )
    std::cout << "\n***** Variation pile-up down ********************************\n" << std::endl;


  // Create output directories, file names, and labels
  TString tmpLabelWindow = util::toTString(nSigTailStart);
  if( tmpLabelWindow.Length() == 1 ) tmpLabelWindow += "0";
  tmpLabelWindow += "-";
  if( nSigTailEnd > 10. ) tmpLabelWindow += "Inf";
  else tmpLabelWindow += util::toTString(nSigTailEnd);
  tmpLabelWindow.ReplaceAll(".","");
  TString outLabel = "Tail_"+uid;
  if( fileNameData.Contains("NVtx") ) {
    TString tmp = fileNameData;
    outLabel += "_"+tmp(tmp.First('NVtx')-3,9);
  }
  outLabel += "_Sig"+tmpLabelWindow+"_"+jetAlgo;
  if( variationCoreDn )
    outLabel += "_VarCoreDn";
  else if( variationCoreUp )
    outLabel += "_VarCoreUp";
  else if( variationExtrapolation )
    outLabel += "_VarExtrapolation";
  else if( variationClosure )
    outLabel += "_VarClosure";
  else if( variationPUUp )
    outLabel += "_VarPUUp";
  else if( variationPUDn )
    outLabel += "_VarPUDn";

  

  if( SHOW_HEADER ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStyleNoteNoTitle();
  }
  gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");

  // +++++ Bins +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating bins" << std::endl;

  std::vector<EtaPtBin*> etaPtBins;
  const unsigned int nEtaBins = binAdm->nEtaBins();
  unsigned int globalBin = 0;
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    double coreScaling = 0.;
    if     ( variationCoreDn ) coreScaling = coreScaleFactor(etaBin,-1);
    else if( variationCoreUp ) coreScaling = coreScaleFactor(etaBin,1);
    else                       coreScaling = coreScaleFactor(etaBin,0);
    std::cout << "  Eta " << etaBin << ": Using core scale factor " << coreScaling << std::endl;

    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
      if( DEBUG ) std::cout << "Setting up bin (" << etaBin << ", " << ptBin << ")" << std::endl;

      // Create eta-pt bin
      EtaPtBin* bin = new EtaPtBin(etaBin,ptBin,nSigCore,coreScaling,binAdm);

      // Define window
      if( fileNameTailWindow == "" ) {
	bin->findWindow(fileNameMC,0,nSigTailStart,nSigTailEnd);
      } else {
	bin->setWindow(fileNameTailWindow,nSigTailStart,nSigTailEnd);
      }

      ++globalBin;

      // Add pt3 bins
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdm->nPtSoftBins(); ++ptSoftBin) {
	bin->addPt3Bin(ptSoftBin,binAdm->ptSoftMax(ptSoftBin),fileNameData,fileNameMC);
      }

      // Add mc truth response for toy asymmetry
      bin->addMCTruthForToyAsym(fileNameMCTruth);

      etaPtBins.push_back(bin);
      if( DEBUG ) std::cout << "Done setting up bin" << std::endl;
    }
  }



  // +++++ Asymmetry and extrapolation plots ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating asymmetry and extrapolation plots" << std::endl;

  // ROOT output
  ROOT_OUT_FILE = new TFile(outLabel+".root","RECREATE");
  
  for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
      it != etaPtBins.end(); ++it ) {
    
    if( DEBUG ) std::cout << "Plotting and extrapolation" << std::endl;

    // Asymmetry
    (*it)->plotAsymmetryDistributions(outLabel,nSigTailStart,nSigTailEnd);
    (*it)->plotMCTruth(outLabel);
    (*it)->plotSpectra(outLabel);
    
    // Extrapolation
    (*it)->extrapolate(minPt3Data,fixDataShape,!variationExtrapolation,variationClosure);  
    (*it)->plotExtrapolation(outLabel);

    // Store tail-start bins
    (*it)->plotTailStart();
  }
  writeLaTeXSlides(outLabel,binAdm);



  // +++++ Scaling factors ++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating scaling factor plots" << std::endl;

  // Find min and max of fasym per pt3Bin
  std::vector<double> minFAsym(etaPtBins.front()->nPt3Bins(),1000.);
  std::vector<double> maxFAsym(etaPtBins.front()->nPt3Bins(),0.);
  for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin(); it != etaPtBins.end(); ++it) {
    for(unsigned int i = 0; i < (*it)->nPt3Bins(); ++i) {
      if( (*it)->fTailMCSmearedGauss(i) < minFAsym.at(i) ) 
	minFAsym.at(i) = (*it)->fTailMCSmearedGauss(i);
      if( (*it)->fTailData(i) > maxFAsym.at(i) )
	maxFAsym.at(i) = (*it)->fTailData(i);
    }
  }

  // For MultiCanvas figure
  std::vector< std::vector<TGraphAsymmErrors*> > vgFAsymData;
  std::vector< std::vector<TH1*> >               vhFAsymMCSmearedGauss;
  std::vector< std::vector<TH1*> >               vhFAsymMCSmeared;
  std::vector< std::vector<TGraphAsymmErrors*> > vgFAsymMCSmeared;
  std::vector< std::vector<TGraphAsymmErrors*> > vgFAsymDataRelToGauss;
  std::vector< std::vector<TH1*> >               vhFAsymMCSmearedRelToGauss;
  std::vector< std::vector<TGraphAsymmErrors*> > vgFAsymMCSmearedRelToGauss;

  // Loop over eta bins
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    TH1* hPtAveMeanData = new TH1D("hPtAveMeanData_Eta"+util::toTString(etaBin),
				   HEADER+";p^{ave}_{T} Bin Idx;<p^{ave}_{T}> (GeV)",
				   binAdm->nPtBins(etaBin),-0.5,binAdm->nPtBins(etaBin)-0.5);
    hPtAveMeanData->SetMarkerStyle(20);

    TH1* hDelta = new TH1D("hDelta_Eta"+util::toTString(etaBin),
			   HEADER+";p^{ave}_{T} (GeV);#Delta = "+FASYMDATA+"(0) - "+FASYMMC+"(0)",
			   binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hDelta->SetMarkerStyle(20);

    TH1* hScale = new TH1D("hScale_Eta"+util::toTString(etaBin),
			   HEADER+";p^{ave}_{T} (GeV);"+SCALE,
			   binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hScale->SetMarkerStyle(20);
    hScale->GetXaxis()->SetMoreLogLabels();
    hScale->GetXaxis()->SetNoExponent();

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


    // Start of tail region in asymmetry (same A for all pt3 bins)
    TH1* hAsymTailStart = new TH1D("hAsymTailStart_EtaBin"+util::toTString(etaBin),HEADER+";"+util::LabelFactory::ptAve()+" (GeV);Asymmetry at Tail Start",binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
    hAsymTailStart->SetLineWidth(LINE_WIDTH);
    hAsymTailStart->GetXaxis()->SetMoreLogLabels();
    hAsymTailStart->GetXaxis()->SetNoExponent();

    TH1* hAsymTailStartEff = static_cast<TH1*>(hAsymTailStart->Clone("hAsymTailStartEff_EtaBin"+util::toTString(etaBin)));
    hAsymTailStartEff->SetLineStyle(2);
    hAsymTailStartEff->SetLineColor(kBlue);

    for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
	it != etaPtBins.end(); ++it) {
      if( (*it)->etaBin() == etaBin ) {
	int bin = 1+(*it)->ptBin();
	hAsymTailStart->SetBinContent(bin,(*it)->tailWindowMin());
	hAsymTailStartEff->SetBinContent(bin,(*it)->tailWindowMinEff());
      }
    }

    // Labels
    TPaveText* label = util::LabelFactory::createPaveText(2);
    label->AddText(MCSMEAR);
    label->AddText(util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin))+",  "+labelWindow(nSigTailStart,nSigTailEnd));
    
    TLegend* leg = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
    leg->AddEntry(hAsymTailStart,ASYM_TAIL_START,"L");
    leg->AddEntry(hAsymTailStartEff,ASYM_TAIL_START_EFF,"L");

    // Plot
    TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
    can->cd();
    hAsymTailStart->GetYaxis()->SetRangeUser(0.071,0.169);
    hAsymTailStart->Draw("HIST");
    hAsymTailStartEff->Draw("HISTsame");
    hAsymTailStart->Draw("HISTsame");
    leg->Draw("same");
    label->Draw("same");
    gPad->RedrawAxis();
    can->SetLogx();
    toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_AsymTailStart");
    toFiles(hAsymTailStart);
    toFiles(hAsymTailStartEff);

    delete hAsymTailStart;
    delete hAsymTailStartEff;
    delete label;
    delete leg;
    delete can;
    
      
    // Loop over pt3 bins
    std::vector<TGraphAsymmErrors*> vgFAsymDataTmp;
    std::vector<TH1*>               vhFAsymMCSmearedGaussTmp;
    std::vector<TH1*>               vhFAsymMCSmearedTmp;
    std::vector<TGraphAsymmErrors*> vgFAsymMCSmearedTmp;
    std::vector<TGraphAsymmErrors*> vgFAsymDataRelToGaussTmp;
    std::vector<TH1*>               vhFAsymMCSmearedRelToGaussTmp;
    std::vector<TGraphAsymmErrors*> vgFAsymMCSmearedRelToGaussTmp;
    for(unsigned int pt3Bin = 0; pt3Bin < binAdm->nPtSoftBins(); ++pt3Bin) {
      
      TH1* hFAsymMCSmeared = new TH1D("hFAsymMCSmeared_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin),
				      HEADER+";p^{ave}_{T} (GeV);"+FASYM+"  (%)",
				      binAdm->nPtBins(etaBin),
				      &(binAdm->ptBinEdges(etaBin).front()));
      hFAsymMCSmeared->GetXaxis()->SetMoreLogLabels();
      hFAsymMCSmeared->GetYaxis()->SetMoreLogLabels();
      hFAsymMCSmeared->GetXaxis()->SetNoExponent();
      hFAsymMCSmeared->GetYaxis()->SetNoExponent();
      hFAsymMCSmeared->GetYaxis()->SetNdivisions(505);
      hFAsymMCSmeared->SetLineWidth(LINE_WIDTH);
      hFAsymMCSmeared->SetLineColor(COLOR_FILLED_ASYM_SMEAR);
      hFAsymMCSmeared->SetLineStyle(1);

      TH1* hFAsymData = static_cast<TH1*>(hFAsymMCSmeared->Clone("hFAsymData_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)));
      hFAsymData->SetMarkerStyle(20);
      hFAsymData->SetMarkerSize(MARKER_SIZE);
      hFAsymData->SetLineColor(kBlack);
      
      TH1* hFAsymMCSmearedGauss = static_cast<TH1*>(hFAsymMCSmeared->Clone("hFAsymMCSmearedGauss_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)));
      hFAsymMCSmearedGauss->SetLineColor(COLOR_GAUSS);
      hFAsymMCSmearedGauss->SetLineStyle(2);
      
      // Loop over pt bins and fill histograms
      for(std::vector<EtaPtBin*>::const_iterator it = etaPtBins.begin();
	    it != etaPtBins.end(); ++it) {
	if( (*it)->etaBin() == etaBin ) {
	  int bin = 1+(*it)->ptBin();
	  hFAsymData->SetBinContent(bin,100.*(*it)->fTailData(pt3Bin));
	  hFAsymData->SetBinError(bin,100.*(*it)->fTailDataErr(pt3Bin));
	  hFAsymMCSmeared->SetBinContent(bin,100.*(*it)->fTailMCSmeared(pt3Bin));
	  hFAsymMCSmeared->SetBinError(bin,100.*(*it)->fTailMCSmearedErr(pt3Bin));
	  hFAsymMCSmearedGauss->SetBinContent(bin,100.*(*it)->fTailMCSmearedGauss(pt3Bin));
	}
      } // End of loop over pt bins
      
      // Write basic histograms to ROOT file for further processing
      toFiles(hFAsymData);
      toFiles(hFAsymMCSmeared);
      toFiles(hFAsymMCSmearedGauss);

      // Convert fAsymData to graph and set x value to mean
      // of ptAve spectrum in each bin
      TGraphAsymmErrors* gFAsymData = util::HistOps::createTGraph(hFAsymData);
      for(int n = 0; n < gFAsymData->GetN(); ++n) {
	gFAsymData->GetX()[n] = hPtAveMeanData->GetBinContent(1+n);
	gFAsymData->SetPointEXhigh(n,hPtAveMeanData->GetBinError(1+n));
	gFAsymData->SetPointEXlow(n,hPtAveMeanData->GetBinError(1+n));
      }
      
      // Convert fAsymMC to graph to be able to plot error band
      TGraphAsymmErrors* gFAsymMCSmeared = util::HistOps::getUncertaintyBand(hFAsymMCSmeared,COLOR_FILLED_ASYM_SMEAR);
      gFAsymMCSmeared->SetLineColor(gFAsymMCSmeared->GetFillColor());
      gFAsymMCSmeared->SetLineStyle(1);
      gFAsymMCSmeared->SetLineWidth(LINE_WIDTH);
      gFAsymMCSmeared->SetFillStyle(HATCH_STYLE);

      // Ratios to Gauss
      TGraphAsymmErrors* gFAsymDataRelToGauss = util::HistOps::createRatioGraph(gFAsymData,hFAsymMCSmearedGauss);
      TGraphAsymmErrors* gFAsymMCSmearedRelToGauss = util::HistOps::createRatioGraph(gFAsymMCSmeared,hFAsymMCSmearedGauss);
      gFAsymMCSmearedRelToGauss->SetFillStyle(gFAsymMCSmeared->GetFillStyle());
      gFAsymMCSmearedRelToGauss->SetLineColor(gFAsymMCSmeared->GetLineColor());
      gFAsymMCSmearedRelToGauss->SetLineWidth(LINE_WIDTH);
      TH1* hFAsymMCSmearedRelToGauss = util::HistOps::createRatioPlot(hFAsymMCSmeared,hFAsymMCSmearedGauss);
      
      // Labels
      label = util::LabelFactory::createPaveText(2);
      label->AddText(util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin))+",  "+labelWindow(nSigTailStart,nSigTailEnd));
      label->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+util::LabelFactory::pt3RelCut(binAdm->ptSoftMax(pt3Bin)));
      
      leg = util::LabelFactory::createLegendWithOffset(3,label->GetSize());
      leg->AddEntry(gFAsymData,DATA,"P");
      leg->AddEntry(gFAsymMCSmeared,MCSMEAR,"LF");
      leg->AddEntry(hFAsymMCSmearedGauss,MCGAUSS,"L");
      
      // Plot
      can = util::HistOps::createTCanvas("can","",500,500);
      can->cd();
      double min = 100.*minFAsym.at(pt3Bin);
      double max = 100.*maxFAsym.at(pt3Bin);
      double delta = max - min;
      hFAsymMCSmeared->GetYaxis()->SetRangeUser(0.9*min,max+1.2*delta);
      hFAsymMCSmeared->Draw("HIST");
      hFAsymMCSmearedGauss->Draw("HISTsame");
      gFAsymMCSmeared->Draw("E2same");
      hFAsymMCSmeared->Draw("HISTsame");
      gFAsymData->Draw("PE1same");
      leg->Draw("same");
      label->Draw("same");
      gPad->RedrawAxis();
      can->SetLogx();
      toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsym");

      // Plot with ratio relative to Gauss
      delete can;
      can = util::HistOps::createRatioTopCanvas();
      TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
      //TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hFAsymData);
      TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hFAsymMCSmeared);
      TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hFAsymMCSmeared,util::LabelFactory::ptAve(),"GeV","#frac{Data, Sim.}{Gauss}",0.81,2.99);
      bRatioBottomFrame->GetYaxis()->SetRangeUser(0.95,1.1*(*std::max_element(gFAsymDataRelToGauss->GetY(),gFAsymDataRelToGauss->GetY()+gFAsymDataRelToGauss->GetN())));
      bRatioBottomFrame->SetLineStyle(hFAsymMCSmearedGauss->GetLineStyle());
      bRatioBottomFrame->SetLineColor(hFAsymMCSmearedGauss->GetLineColor());
      bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
      bRatioBottomFrame->GetXaxis()->SetNoExponent();
      bRatioTopFrame->GetYaxis()->SetRangeUser(0.9*min,max+1.2*delta);
      can->cd();
      bRatioTopFrame->Draw("HIST");
      hFAsymMCSmearedGauss->Draw("HISTsame");
      gFAsymMCSmeared->Draw("E2same");
      hFAsymMCSmeared->Draw("HISTsame");
      gFAsymData->Draw("PE1same");
      leg->Draw("same");
      label->Draw("same");
      gPad->RedrawAxis();
      can->SetLogx();
      bRatioBottomPad->Draw();
      bRatioBottomPad->cd();
      bRatioBottomFrame->Draw("HIST");
      gFAsymMCSmearedRelToGauss->Draw("E2same");
      hFAsymMCSmearedRelToGauss->Draw("HISTsame");
      gFAsymDataRelToGauss->Draw("PE1same");
      bRatioBottomPad->SetLogx();
      bRatioBottomPad->RedrawAxis();
      toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsymRelToGaussBottom");

      // Store for MultiCanvas figure
      vgFAsymDataTmp.push_back(gFAsymData);
      vhFAsymMCSmearedGaussTmp.push_back(hFAsymMCSmearedGauss);
      vhFAsymMCSmearedTmp.push_back(hFAsymMCSmeared);
      vgFAsymMCSmearedTmp.push_back(gFAsymMCSmeared);
      vgFAsymDataRelToGaussTmp.push_back(gFAsymDataRelToGauss);
      vhFAsymMCSmearedRelToGaussTmp.push_back(hFAsymMCSmearedRelToGauss);
      vgFAsymMCSmearedRelToGaussTmp.push_back(gFAsymMCSmearedRelToGauss);

//       delete hFAsymMCSmeared;
//       delete hFAsymMCSmearedGauss;
//       delete gFAsymMCSmeared;
//       delete gFAsymDataRelToGauss;
//       delete gFAsymMCSmearedRelToGauss;
//       delete hFAsymMCSmearedRelToGauss;
      //      delete gFAsymData;
      delete hFAsymData;
      delete bRatioBottomPad;
      delete bRatioTopFrame;
      delete bRatioBottomFrame;
      delete label;
      delete leg;
      delete can;
    } // End of loop over pt3 bins
    vgFAsymData.push_back(vgFAsymDataTmp);
    vhFAsymMCSmearedGauss.push_back(vhFAsymMCSmearedGaussTmp);
    vhFAsymMCSmeared.push_back(vhFAsymMCSmearedTmp);
    vgFAsymMCSmeared.push_back(vgFAsymMCSmearedTmp);
    vgFAsymDataRelToGauss.push_back(vgFAsymDataRelToGaussTmp);
    vhFAsymMCSmearedRelToGauss.push_back(vhFAsymMCSmearedRelToGaussTmp);
    vgFAsymMCSmearedRelToGauss.push_back(vgFAsymMCSmearedRelToGaussTmp);

    
    // Plots
    can = util::HistOps::createTCanvas("can","",500,500);
    
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
    toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_Delta");

    TH1* hScaleFrame = static_cast<TH1D*>(hScale->Clone("hScaleFrame_Eta"+util::toTString(etaBin)));
    hScaleFrame->SetLineColor(kBlack);
    hScaleFrame->SetLineStyle(2);
    hScaleFrame->SetMarkerStyle(1);
    for(int i = 1; i <= hScaleFrame->GetNbinsX(); ++i) {
      hScaleFrame->SetBinContent(i,1.);
      hScaleFrame->SetBinError(i,0.);
    }
    hScaleFrame->GetYaxis()->SetRangeUser(0.01,2.99);
    hScaleFrame->GetXaxis()->SetMoreLogLabels();
    can->cd();
    hScaleFrame->Draw("HIST");
    hScale->Draw("PE1same");
    can->SetLogx();
    toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_ScaleFactors");

    toFiles(hDelta);    
    toFiles(hDeltaFrame);   
    toFiles(hScale);    
    toFiles(hScaleFrame);   
    toFiles(hPtAveMeanData);

    delete hDelta;
    delete hScale;
    delete hDeltaFrame;
    delete hScaleFrame;
    delete hPtAveMeanData;
    delete can;
  } // End of loop over eta bins

  // MultiCanvas figure

  // Loop over pt3 bins
  for(unsigned int pt3Bin = 0; pt3Bin < binAdm->nPtSoftBins(); ++pt3Bin) {

    util::MultiCanvas* mc = new util::MultiCanvas(outLabel+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsymRelToGaussBottom",3,2,5,true);
	
    TLegend* leg = util::LabelFactory::createLegendWithOffset(3,0.2+util::LabelFactory::lineHeightMultiCan(),util::LabelFactory::lineHeightMultiCan());
    mc->adjustLegend(leg);
    mc->markForDeletion(leg);
      
    // Frames
    TH1* hFrame = static_cast<TH1*>(vhFAsymMCSmeared.at(0).at(pt3Bin)->Clone("FrameForMultiCanvas"));
    hFrame->Reset();
    hFrame->GetYaxis()->SetTitle(FASYM+"  (%)    ");
    double min = 100.*minFAsym.at(pt3Bin);
    double max = 100.*maxFAsym.at(pt3Bin);
    double delta = max - min;
    hFrame->GetYaxis()->SetRangeUser(0.8*min,max+0.6*delta);
    mc->markForDeletion(hFrame);

    for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

      // Labels
      TPaveText* label = util::LabelFactory::createPaveText(1,0.45);
      label->AddText(util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
      mc->markForDeletion(label);
      mc->adjustPaveText(label);
      if( etaBin == 0 ) {
	leg->AddEntry(vgFAsymData.at(etaBin).at(pt3Bin),DATA,"P");
	leg->AddEntry(vgFAsymMCSmeared.at(etaBin).at(pt3Bin),MCSMEAR,"LF");
	leg->AddEntry(vhFAsymMCSmearedGauss.at(etaBin).at(pt3Bin),MCGAUSS,"L");
      }

      // Plot
      TH1* mFrame = mc->mainFrame(hFrame,etaBin);
      TH1* rFrame = mc->ratioFrame(hFrame,"#frac{Dat, Sim}{Gauss}",0.81,2.99,etaBin);
      min = 0.95;
      max = 1.66;
      if( etaBin == 4 ) max = 2.4;
      if( nSigTailStart > 2. ) {
	min = 0.95;
	max = 4.1;
	if( etaBin == 4 ) max = 9.5;
      }
      rFrame->GetYaxis()->SetRangeUser(min,max);
      rFrame->SetLineStyle(vhFAsymMCSmearedGauss.at(0).at(pt3Bin)->GetLineStyle());
      rFrame->SetLineColor(vhFAsymMCSmearedGauss.at(0).at(pt3Bin)->GetLineColor());
      mc->canvas()->cd();
      TPad* padt = mc->mainPad(etaBin);
      padt->Draw();
      padt->cd();
      mFrame->Draw();
      vhFAsymMCSmearedGauss.at(etaBin).at(pt3Bin)->Draw("HISTsame");
      vgFAsymMCSmeared.at(etaBin).at(pt3Bin)->Draw("E2same");
      vhFAsymMCSmeared.at(etaBin).at(pt3Bin)->Draw("HISTsame");
      vgFAsymData.at(etaBin).at(pt3Bin)->Draw("PE1same");
      label->Draw("same");
      gPad->RedrawAxis();
      TPad* padr = mc->ratioPad(etaBin);
      padr->Draw();
      padr->cd();
      rFrame->Draw("HIST");
      vgFAsymMCSmearedRelToGauss.at(etaBin).at(pt3Bin)->Draw("E2same");
      vhFAsymMCSmearedRelToGauss.at(etaBin).at(pt3Bin)->Draw("HISTsame");
      vgFAsymDataRelToGauss.at(etaBin).at(pt3Bin)->Draw("PE1same");
      gPad->RedrawAxis();

    } // end of loop over eta bins
    TPaveText* label = util::LabelFactory::createPaveTextWithOffset(1,1.,0.2,util::LabelFactory::lineHeightMultiCan());
    label->AddText(labelWindow(nSigTailStart,nSigTailEnd)+",  "+util::LabelFactory::pt3RelCut(binAdm->ptSoftMax(pt3Bin)));
    mc->adjustPaveText(label);
    mc->markForDeletion(label);

    mc->canvas()->cd();
    TPad* padt = mc->mainPad(5);
    padt->Draw();
    padt->cd();
    label->Draw();
    leg->Draw();      
    mc->moreLogLabelsX();
    mc->noExponentX();
    mc->setLogx();
    toFiles(mc->canvas(),outLabel+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsymRelToGaussBottom");
    delete mc;
  } // End of loop over pt3 bins

  for(unsigned int i = 0; i < vgFAsymData.size(); ++i) {
    for(unsigned int j = 0; j < vgFAsymData.at(i).size(); ++j) {
      delete vgFAsymData.at(i).at(j);
      delete vhFAsymMCSmearedGauss.at(i).at(j);
      delete vhFAsymMCSmeared.at(i).at(j);
      delete vgFAsymMCSmeared.at(i).at(j);
      delete vgFAsymDataRelToGauss.at(i).at(j);
      delete vhFAsymMCSmearedRelToGauss.at(i).at(j);
      delete vgFAsymMCSmearedRelToGauss.at(i).at(j);
    }
  }

  ROOT_OUT_FILE->Close();
  delete ROOT_OUT_FILE;

  printWindowBorders(etaPtBins,binAdm,nSigTailStart);
//   printMCClosure(etaPtBins,binAdm);
//   printExtrapolation(etaPtBins,binAdm);

  for(std::vector<EtaPtBin*>::iterator it = etaPtBins.begin();
      it != etaPtBins.end(); ++it ) {
    delete *it;
  }
  delete binAdm;

  // Clean up working directory
  if( archivePlots ) {
    std::cout << "Cleaning up working directory" << std::endl;
    if( EPS_OUTPUT ) {
      TString filesInTar = outLabel+"*.tex "+outLabel+"*.eps ";
      gROOT->ProcessLine(".! tar -zcf "+outLabel+".tar.gz "+filesInTar);
      gROOT->ProcessLine(".! rm "+filesInTar);
      std::cout << "  Plots in eps format: "+outLabel+".tar.gz" << std::endl;
    } else {
      gROOT->ProcessLine(".! rm "+outLabel+"*.tex");
    }
    std::cout << "  Plots in ROOT format: "+outLabel+".root" << std::endl;
  }
}


void getTailScalingFactors() {

  std::vector<double> nSigTailStart;
  nSigTailStart.push_back(2.0);
  //nSigTailStart.push_back(2.5);
  nSigTailStart.push_back(3.0);

  double  nSigTailEnd = 1000.;
  for(std::vector<double>::const_iterator it = nSigTailStart.begin();
      it != nSigTailStart.end(); ++it) {
    EPS_OUTPUT = true;
    getTailScalingFactors(*it,nSigTailEnd,false,false,false,false,false,false);
    EPS_OUTPUT = false;
    getTailScalingFactors(*it,nSigTailEnd,true,false,false,false,false,false);
    getTailScalingFactors(*it,nSigTailEnd,false,true,false,false,false,false);
    getTailScalingFactors(*it,nSigTailEnd,false,false,true,false,false,false);
    getTailScalingFactors(*it,nSigTailEnd,false,false,false,true,false,false);
    getTailScalingFactors(*it,nSigTailEnd,false,false,false,false,true,false);
    getTailScalingFactors(*it,nSigTailEnd,false,false,false,false,false,true);
  }
}









//////////////////////////////// IMPLEMENTATIONS ///////////////////////////////////////


// ------------------------------------------------------------------------------------
EtaPtBin::EtaPtBin(unsigned int etaBinIdx, unsigned int ptBinIdx, double nSigCore, double coreScaleFactor, const sampleTools::BinningAdmin* binAdm)
  : etaBin_(etaBinIdx), ptBin_(ptBinIdx),
    nSigCore_(nSigCore), coreScalingFactor_(coreScaleFactor),
    exMin_(0.), exMax_(0.19),
    binAdm_(binAdm) {

  tailWindowMinEff_ = 0.;
  tailWindowMaxEff_ = 0.;
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

  fExMC_ = new TF1("fExMC"+binId(),"expo",exMin_,exMax_);
  fExMC_->SetLineWidth(LINE_WIDTH);


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

  binLabel_ = util::LabelFactory::createPaveText(2);
  binLabel_->AddText(util::LabelFactory::etaCut(binAdm_->etaMin(etaBin_),binAdm_->etaMax(etaBin_))+",  "+util::LabelFactory::ptAveCut(binAdm_->ptMin(etaBin_,ptBin_),binAdm_->ptMax(etaBin_,ptBin_)));
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
  TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
  can->cd();
  TH1* hFrame = new TH1D("hFrame",HEADER+";Threshold "+util::LabelFactory::pt3RelMax(),1000,0.,PT3PLOTMAX);
  hFrame->SetNdivisions(505);
  hFrame->GetYaxis()->SetNoExponent();
  hFrame->GetYaxis()->SetNdivisions(550);
  //hFrame->GetYaxis()->SetRangeUser(0.,2.3*gFTailMC_->GetY()[gFTailMC_->GetN()-1]);
  double minY = *std::min_element(gFTailMCGauss_->GetY(),gFTailMCGauss_->GetY()+gFTailMCGauss_->GetN());
  minY = std::min(minY,std::min(gFTailData_->GetY()[0],fExMC_->Eval(0.)));
  minY *= 0.8;
  double maxY = 1.2;

  if( hasToyMC() ) {
    TLegend* leg = util::LabelFactory::createLegendWithOffset(3,binLabel_->GetSize());
    leg->AddEntry(gFTailMC_,MCSMEAR,"PL");
    leg->AddEntry(gFTailMC_,FASYMMC+"("+util::LabelFactory::pt3RelMax()+"#rightarrow0) = "+util::toTString(extraMC(),4)+" #pm "+util::toTString(extraMCErr(),4),"");
    leg->AddEntry(gFTailToyMC_,FASYMTOY+" = "+util::toTString(toyMC(),4)+" #pm "+util::toTString(toyMCErr(),4),"P");

    hFrame->GetYaxis()->SetTitle(FASYMMC);
    hFrame->GetYaxis()->SetRangeUser(minY,0.7);
    hFrame->Draw();
    fExMC_->Draw("same");
    if( hasSymMCTruth() ) gFTailMCTruth_->Draw("PE1same");
    gFTailToyMC_->Draw("PE1same");
    gFTailMC_->Draw("PE1same");
    binLabel_->DrawClone("same");
    leg->Draw("same");
    gPad->SetLogy(1);
    toFiles(can,outNameId+binId()+"_ExtrapolationMCClosure");

    delete leg;
  }

  // Extrapolation MC + Data
  hFrame->GetYaxis()->SetRangeUser(minY,maxY);
  hFrame->GetYaxis()->SetTitle(FASYM);
  TLegend* leg = util::LabelFactory::createLegendWithOffset(5,binLabel_->GetSize());
  leg->AddEntry(gFTailData_,DATA,"PL");
  leg->AddEntry(gFTailData_,FASYMDATA+"("+util::LabelFactory::pt3RelMax()+"#rightarrow0) = "+util::toTString(extraData(),4)+" #pm "+util::toTString(extraDataErr(),4),"");
  leg->AddEntry(gFTailMC_,MCSMEAR,"PL");
  leg->AddEntry(gFTailMC_,FASYMMC+"("+util::LabelFactory::pt3RelMax()+"#rightarrow0) = "+util::toTString(extraMC(),4)+" #pm "+util::toTString(extraMCErr(),4),"");
  leg->AddEntry(gFTailMCGauss_,MCGAUSS,"P");

  can->cd();
  hFrame->Draw();
  fExMC_->Draw("same");
  fExData_->Draw("same");
  gFTailMC_->Draw("PE1same");
  gFTailData_->Draw("PE1same");
  gFTailMCGauss_->Draw("PE1same");
  binLabel_->DrawClone("same");
  leg->Draw("same");
  gPad->SetLogy(1);
  toFiles(can,outNameId+binId()+"_Extrapolation");
  delete leg;
  gPad->SetLogy(0);

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
  toFiles(can,outNameId+binId()+"_SpreadMC");  
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
  toFiles(can,outNameId+binId()+"_SpreadData");  
  delete leg;

  delete hFrame;
  delete can;

  if( DEBUG ) std::cout << "Leaving EtaPtBin::plotExtrapolation()" << std::endl;
}


// ------------------------------------------------------------------------------------
void EtaPtBin::plotTailStart() const {
  TH1* hTailStartBins = new TH1I("hTailStartBins_EtaBin"+util::toTString(etaBin())+"_PtBin"+util::toTString(ptBin()),HEADER+";Pt3Bin",nPt3Bins(),0,nPt3Bins());
  TH1* hTailEndBins = static_cast<TH1*>(hTailStartBins->Clone("hTailEndBins_EtaBin"+util::toTString(etaBin())+"_PtBin"+util::toTString(ptBin())));
  for(unsigned int i = 0; i < nPt3Bins(); ++i) {
    hTailStartBins->SetBinContent(1+i,pt3Bins_.at(i)->tailStartBin());
    hTailEndBins->SetBinContent(1+i,pt3Bins_.at(i)->tailEndBin());
  }
  toFiles(hTailStartBins);
  toFiles(hTailEndBins);
}


// ------------------------------------------------------------------------------------
void EtaPtBin::addPt3Bin(unsigned int pt3Bin, double thres, const TString &fileNameData, const TString &fileNameMC) {
  pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bin,thres,nSigCore_,coreScalingFactor_,tailMinBin_,tailMaxBin_,binAdm_));
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
void EtaPtBin::addMCTruthForToyAsym(const TString &fileName) {
  if( hasToyMC() ) delete toyMC_;
  toyMC_ = new Pt3Bin(fileName,etaBin_,ptBin_,nSigCore_,coreScalingFactor_,nPtAsymBins(),tailMinBin_,tailMaxBin_,binAdm_);
  hasToyMC_ = true;
}


// Find min and max asymmetry
// The borders are specified in number of sigmas, where sigma is the Gaussian
// width of the *smeared* asymmetry distribution from fileName
// ------------------------------------------------------------------------------------
void EtaPtBin::findWindow(const TString &fileName, unsigned int pt3Bin, double nSigTailWindowMin, double nSigTailWindowMax) {

  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft"+util::toTString(pt3Bin);

  TH1* hOrig = util::FileOps::readTH1(fileName,histName);
  hOrig->GetXaxis()->SetRangeUser(-1.,1);
  TH1* h = 0;
  double width = 0.;
  double widthOrig = 0.;
  correctAsymmetryWidth(hOrig,nSigCore_,coreScalingFactor_,h,widthOrig,width);
  delete hOrig;

  tailWindowMin_ = nSigTailWindowMin*width;
  tailWindowMax_ = nSigTailWindowMax*width;
  tailMinBin_ = h->FindBin(std::abs(tailWindowMin_));
  tailMaxBin_ = h->FindBin(std::abs(tailWindowMax_));
  tailWindowMinEff_ = h->GetXaxis()->GetBinLowEdge(tailMinBin_);
  tailWindowMaxEff_ = h->GetXaxis()->GetBinUpEdge(tailMaxBin_);

  delete h;
  
  binLabel_->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+labelWindow(nSigTailWindowMin,nSigTailWindowMax));
}


// Set window borders as specified in previous run and stored
// in histograms; nSigTailWindow??? are only needed for the
// labels
// ------------------------------------------------------------------------------------
void EtaPtBin::setWindow(const TString &fileName, double nSigTailWindowMin, double nSigTailWindowMax) {
  TH1* hMin = util::FileOps::readTH1(fileName,"hTailStartBins_EtaBin"+util::toTString(etaBin())+"_PtBin"+util::toTString(ptBin()));
  TH1* hMax = util::FileOps::readTH1(fileName,"hTailEndBins_EtaBin"+util::toTString(etaBin())+"_PtBin"+util::toTString(ptBin()));
  tailMinBin_ = static_cast<int>(hMin->GetBinContent(1));
  tailMaxBin_ = static_cast<int>(hMax->GetBinContent(1));
  delete hMin;
  delete hMax;

  hMin = util::FileOps::readTH1(fileName,"hAsymTailStart_EtaBin"+util::toTString(etaBin()));
  tailWindowMin_ = hMin->GetBinContent(1+ptBin());
  tailWindowMax_ = 10000.;
  delete hMin;

  hMin = util::FileOps::readTH1(fileName,"hAsymTailStartEff_EtaBin"+util::toTString(etaBin()));
  tailWindowMinEff_ = hMin->GetBinContent(1+ptBin());
  tailWindowMaxEff_ = 10000.;
  delete hMin;
  
  binLabel_->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+labelWindow(nSigTailWindowMin,nSigTailWindowMax));
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
  gFTailMCGauss_->SetMarkerColor(COLOR_GAUSS);
  gFTailMCGauss_->SetMarkerSize(MARKER_SIZE);
  gFTailMCGauss_->SetLineWidth(LINE_WIDTH);
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
    gFTailToyMC_->SetMarkerStyle(29);
    gFTailToyMC_->SetMarkerSize(1.8);
    gFTailToyMC_->SetMarkerColor(kRed);
    gFTailToyMC_->SetLineColor(gFTailToyMC_->GetMarkerColor());
  }

  // Extrapolate
  gFTailMC_->Fit(fExMC_,"0QR");
  fExMC_->SetLineColor(COLOR_LINE_ASYM_SMEAR);

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
  //extraMC_ = useExtrapolatedValue ? fExMC_->GetParameter(0) : fExMC_->Eval(gFTailMC_->GetX()[0]);
  extraMC_ = useExtrapolatedValue ? fExMC_->Eval(0) : fExMC_->Eval(gFTailMC_->GetX()[0]);
  extraMCErr_ = fExMC_->GetParError(0)*fExMC_->Eval(0);	// Parameter correlations do not need to be taken into account because the correlation term vanishes at x = 0


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
  extraData_ = useExtrapolatedValue ? fExData_->Eval(0) : fExData_->Eval(gFTailData_->GetX()[0]);
  extraDataErr_ = fExData_->GetParError(0)*fExData_->Eval(0);	// Parameter correlations do not need to be taken into account because the correlation term vanishes at x = 0
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

  std::cout << "  Eta " << etaBin() << ", Pt " << ptBin() << ": " << deltaEx_ << " +/- " << deltaExErr_ << " --> " << scalingFactor_ << " +/- " << scalingFactorErr_ << std::endl;

  if( DEBUG ) std::cout << "Leaving EtaPtBin::extrapolate()" << std::endl;
}




// Constructor for asymmetry
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Threshold, double nSigCore, double coreScaleFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm)
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(pt3Bin), pt3Thres_(pt3Threshold),
    coreScalingFactor_(coreScaleFactor), 
    tailStartBin_(windowMinBin), tailEndBin_(windowMaxBin) {

  // Non-read distributions
  hResp_ = 0;
  hSymResp_ = 0;
  fGaussMCTruth_ = 0;

  // Get spectra
  TString histName = "hPtAveCombined_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft"+util::toTString(pt3Bin);
  hPtAveSpecData_ = util::FileOps::readTH1(fileNameData,histName,"Data"+histName);
  util::HistOps::setAxisTitles(hPtAveSpecData_,"p^{ave}_{T}","GeV","events");      
  hPtAveSpecData_->SetTitle(HEADER);
  setStyleDataMarker(hPtAveSpecData_);
  ptAveMeanData_ = hPtAveSpecData_->GetMean();
  ptAveMeanDataErr_ = hPtAveSpecData_->GetMeanError();
  hPtAveSpecMC_ = util::FileOps::readTH1(fileNameMC,histName,"MC"+histName);
  util::HistOps::setAxisTitles(hPtAveSpecMC_,"p^{ave}_{T}","GeV","events");      
  hPtAveSpecMC_->SetTitle(HEADER);
  setStyleMCFilled(hPtAveSpecMC_);

  // Get asymmetry distributions
  hAsymData_ = readHistAsym(fileNameData,"Data");
  hAsymMC_ = readHistAsym(fileNameMC,"MC");
  hAsymData_->UseCurrentStyle();
  hAsymData_->SetTitle(HEADER);
  hAsymMC_->UseCurrentStyle();
  hAsymMC_->SetTitle(HEADER);
  
  // Smear MC asymmetry
  sig_ = 0.;
  correctAsymmetryWidth(hAsymMC_,nSigCore,coreScalingFactor_,hAsymMCSmeared_,sig_,sigSmeared_);

  // Total entries in asymmetry histograms
  // (to compute statistical uncertainty on fasym)
  // factor 1/2 bc each entry was filled twice
  double nTotalEntriesMC = hAsymMC_->GetEntries()/2.;
  double nTotalEntriesData = hAsymData_->GetEntries()/2.;

  // Get relative number of entries in tail
  getFTail(hAsymMC_,nTotalEntriesMC,tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
  getFTail(hAsymMCSmeared_,nTotalEntriesMC,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
  double minEff = hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_);
  double maxEff = min(1.,hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_));
  // 2*OneSidedGaussianTail/TotalGauss
  fNMCSmearedGauss_ = (erf(maxEff/sqrt(2.)/sigSmeared_) - erf(minEff/sqrt(2.)/sigSmeared_));

  getFTail(hAsymData_,nTotalEntriesData,tailStartBin_,tailEndBin_,fNData_,fNDataErr_);

  setStyleDataMarker(hAsymData_);
  setStyleMCFilled(hAsymMC_);
  setStyleMCFilled(hAsymMCSmeared_);
  initBinLabel(adm);
}


// Constructor for toy asymmetry from MC truth response
//  fAsymMCSmeared from toy MC (MC truth --> asymmetry --> smearing)
// ------------------------------------------------------------------------------------
Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm) 
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(1000), pt3Thres_(1000.),
    coreScalingFactor_(coreScaleFactor),
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
  util::HistOps::setAxisTitles(hAsymMC_,"|Asymmetry|","","events",true);
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
  sig_ = 0.;
  correctAsymmetryWidth(hAsymMC_,nSigCore,coreScalingFactor_,hAsymMCSmeared_,sig_,sigSmeared_);
  if( DEBUG ) std::cout << "ok" << std::endl;

  if( DEBUG ) std::cout << "  Getting relative number of tail events in toy asymmetry distribution  . . .  " << std::flush;
  // Get relative number of entries in tail
  getFTail(hAsymMC_,entries,tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
  getFTail(hAsymMCSmeared_,entries,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
  if( DEBUG ) std::cout << "ok" << std::endl;

  // Set style
  initBinLabel(adm,true);

  if( DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Leaving constructor for toy asymmetry" << std::endl;
}


// ------------------------------------------------------------------------------------
void Pt3Bin::initBinLabel(const sampleTools::BinningAdmin* adm, bool isMCTruth) {
  if( DEBUG ) std::cout << "Pt3Bin::initBinLabel(): Creating bin label . . . " << std::flush;
  if( isMCTruth ) binLabel_ = util::LabelFactory::createPaveText(1);
  else binLabel_ = util::LabelFactory::createPaveText(2);
  binLabel_->AddText(util::LabelFactory::etaCut(adm->etaMin(etaBin_),adm->etaMax(etaBin_))+",  "+util::LabelFactory::ptAveCut(adm->ptMin(etaBin_,ptBin_),adm->ptMax(etaBin_,ptBin_)));
  if( !isMCTruth ) binLabel_->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+util::LabelFactory::pt3RelCut(adm->ptSoftMax(pt3Bin_)));
  if( DEBUG ) std::cout << "ok" << std::endl;
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
// Histograms are assumed to have been symmetrised, i.e.
// for each event, +/-A has been filled into the histogram
// ------------------------------------------------------------------------------------
TH1* Pt3Bin::readHistAsym(const TString &fileName, const TString &id) const {
  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft"+util::toTString(pt3Bin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,id+histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  h->Scale(2./h->Integral("width"));
  util::HistOps::setAxisTitles(h,"|Asymmetry|","","events",true);      
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
  histName += "_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft0";
  TH1* h = util::FileOps::readTH1(fileName,histName,type+histName);
  h->GetXaxis()->SetRangeUser(0.,2.);
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
  double nTotal = h->Integral();
  double nTail = 2.*h->Integral(start,end);
  fTail = nTail/nTotal;
  
  // Total number of events for uncertainty calculation
  // from number of entries (assume no overflow)
  nTotal = entries;
  fTailErr = sqrt( fTail*(1.-fTail)/nTotal );
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMC(const TString &outNameId) const {

  TLegend* leg = util::LabelFactory::createLegendWithOffset(2,2);
  leg->AddEntry(hAsymData_,DATA,"P");
  leg->AddEntry(hAsymMC_,MC,"F");

  // Log scale
  util::HistOps::setYRange(hAsymMC_,5,3E-5);
  util::HistOps::setYRange(hAsymData_,5,3E-5);
  hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);

  TString canName = outNameId+"PtAsym";
  TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);
  can->cd();
  hAsymMC_->GetXaxis()->SetRangeUser(0.,1.);
  hAsymMC_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"PtAsymLog");

  // Linear scale
  util::HistOps::setYRange(hAsymMC_,5);
  hAsymMC_->GetXaxis()->SetRangeUser(0.,0.4);
  can->cd();
  hAsymMC_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(0);
  toFiles(can,outNameId+"PtAsym");
  
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
//   toFiles(can,outNameId+"PtAsymFit");

  // Ratio data / MC
  TH1 *hRatio = util::HistOps::createRatioPlot(hAsymData_,hAsymMC_);
//   TH1 *hRatioFrame = util::HistOps::createRatioFrame(hAsymData_,"Data / MC",0.,3.);
//   can->cd();
//   hRatioFrame->Draw();
//   hRatio->Draw("PE1same");
//   label_->info(etaBin_,ptBin_)->Draw("same");
//   label_->eta(etaBin_)->Draw("same");
//   can->SetLogy(0);
//   toFiles(can,outNameId+"PtAsymLogRatio");

//   hRatioFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
//   can->cd();
//   hRatioFrame->Draw();
//   hRatio->Draw("PE1same");
//   label_->info(etaBin_,ptBin_)->Draw("same");
//   label_->eta(etaBin_)->Draw("same");
//   can->SetLogy(0);
//   toFiles(can,outNameId+"PtAsymRatio");

  // Bottom ratio plot
  delete can;
  can = util::HistOps::createRatioTopCanvas();
  TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
  TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hAsymMC_);
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,static_cast<TString>(hAsymMC_->GetXaxis()->GetTitle()),"","#frac{Data}{Sim.}",0.61,1.39);
  bRatioBottomFrame->SetFillColor(0);
  bRatioBottomFrame->SetLineWidth(LINE_WIDTH);
  can->cd();
  bRatioTopFrame->GetXaxis()->SetRangeUser(0.,0.2);
  bRatioTopFrame->GetXaxis()->SetNdivisions(505);
  util::HistOps::setYRange(bRatioTopFrame,2+leg->GetNRows());
  bRatioTopFrame->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetRangeUser(0.,0.2);
  bRatioBottomFrame->GetXaxis()->SetNdivisions(505);
  bRatioBottomFrame->Draw();
  hRatio->Draw("PE1same");
  toFiles(can,outNameId+"PtAsymBottomRatio");

  delete bRatioTopFrame;
  delete bRatioBottomFrame;
  delete bRatioBottomPad;
  delete leg;
  delete can;

  hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);
  hAsymData_->GetXaxis()->SetRangeUser(-1.,1.);
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotAsymmetryDataMCSmeared(const TString &outNameId, double nSigTailStart, double nSigTailEnd) const {

  TLegend* leg = util::LabelFactory::createLegendWithOffset(2,2);
  leg->AddEntry(hAsymData_,DATA,"P");
  leg->AddEntry(hAsymMCSmeared_,MCSMEAR,"F");

  // Log scale
  util::HistOps::setYRange(hAsymMCSmeared_,5,3E-5);
  util::HistOps::setYRange(hAsymData_,5,3E-5);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);

  TString canName = outNameId+"PtAsym";
  TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"PtSmearAsymLog");

  // Including tail window
  TLine* win = new TLine(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
			 hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),10.);
  win->SetLineWidth(3);
  win->SetLineColor(kBlue);

  TArrow* arr = new TArrow(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),3.,
			   hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_)+0.2,3.);
  arr->SetLineWidth(3);
  arr->SetLineColor(win->GetLineColor());
  arr->SetArrowSize(0.05);
  arr->SetAngle(30);

  // Gaussian contribution
  TF1* gauss = new TF1("gauss","gaus",0.,1.);
  gauss->SetParameter(0,2./sqrt(2.*M_PI)/sigmaSmeared());
  gauss->SetParameter(1,0.);
  gauss->SetParameter(2,sigmaSmeared());
  gauss->SetLineWidth(2);
  gauss->SetLineColor(COLOR_GAUSS);
  gauss->SetFillStyle(HATCH_STYLE);
  gauss->SetFillColor(gauss->GetLineColor());
  gauss->SetRange(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),gauss->GetX(3E-5,0.,1.));
  
  TLegend* legWin = util::LabelFactory::createLegendWithOffset(4,2);
  legWin->AddEntry(hAsymData_,DATA,"P");
  legWin->AddEntry(hAsymMCSmeared_,MCSMEAR,"F");
  legWin->AddEntry(gauss,MCGAUSS,"F");
  legWin->AddEntry(arr,labelWindow(nSigTailStart,nSigTailEnd),"L");

  TPaveText* labelFAsym = util::LabelFactory::createPaveTextWithOffset(4,0.6,2+legWin->GetNRows());
  int nDigits = 4;
  labelFAsym->AddText(FASYMDATA+" = "+util::toTString(fTailData(),nDigits)+" #pm "+util::toTString(fTailDataErr(),nDigits));
  labelFAsym->AddText(FASYMMC+" = "+util::toTString(fTailMCSmeared(),nDigits)+" #pm "+util::toTString(fTailMCSmearedErr(),nDigits));
  labelFAsym->AddText(FASYMGAUSS+" = "+util::toTString(fTailMCSmearedGauss(),nDigits));

  hAsymMCSmeared_->GetXaxis()->SetNdivisions(505);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  gauss->Draw("same");
  hAsymData_->Draw("PE1same");
  win->Draw("same");
  arr->Draw();			// Arrow not drawn if option "same" called!
  binLabel_->Draw("same");
  legWin->Draw("same");
  labelFAsym->Draw("same");
  can->SetLogy(1);
  gPad->RedrawAxis();
  toFiles(can,outNameId+"PtSmearAsymTail");
  delete legWin;
  delete arr;

  // Without data to illustrate window definition
  if( pt3Bin_ == 0 ) {
    legWin = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
    legWin->AddEntry(hAsymMCSmeared_,MCSMEAR,"F");
    legWin->AddEntry(win,"Tail","F");

    can->cd();
    hAsymMCSmeared_->Draw("HISTE");
    win->Draw("same");
    binLabel_->Draw("same");
    legWin->Draw("same");
    can->SetLogy(1);
    toFiles(can,outNameId+"WindowDef");
    delete legWin;
  }
  delete win;

  // Linear scale
  util::HistOps::setYRange(hAsymMCSmeared_,5);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,0.4);
  can->cd();
  hAsymMCSmeared_->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(0);
  toFiles(can,outNameId+"PtSmearAsym");


  // Bottom ratio plot
  TH1 *hRatio = util::HistOps::createRatioPlot(hAsymData_,hAsymMCSmeared_);
  delete can;
  can = util::HistOps::createRatioTopCanvas();
  TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
  TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hAsymMCSmeared_);
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,static_cast<TString>(hAsymMC_->GetXaxis()->GetTitle()),"","#frac{Data}{Sim.}",0.61,1.39);
  bRatioBottomFrame->SetFillColor(0);
  bRatioBottomFrame->SetLineWidth(LINE_WIDTH);
  can->cd();
  bRatioTopFrame->GetXaxis()->SetRangeUser(0.,0.2);
  bRatioTopFrame->GetXaxis()->SetNdivisions(505);
  util::HistOps::setYRange(bRatioTopFrame,5);
  bRatioTopFrame->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetRangeUser(0.,0.2);
  bRatioBottomFrame->GetXaxis()->SetNdivisions(505);
  bRatioBottomFrame->Draw();
  hRatio->Draw("PE1same");
  toFiles(can,outNameId+"PtSmearAsymBottomRatio");

  delete bRatioTopFrame;
  delete bRatioBottomFrame;
  delete bRatioBottomPad;
  delete leg;
  delete can;

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);
  hAsymData_->GetXaxis()->SetRangeUser(-1.,1.);
}


// ------------------------------------------------------------------------------------
void Pt3Bin::plotSpectra(const TString &outNameId) const {

  TString canName = outNameId+"_PtAveSpectrum";
  TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);

  // Populated x-region
  int xMinBin = 1;
  int xMaxBin = 1000;
  util::HistOps::findXRange(hPtAveSpecData_,xMinBin,xMaxBin);
  xMinBin++;
  xMaxBin--;
  double yMin = 0.01;
  if( hPtAveSpecData_->GetBinCenter(xMaxBin) > 900. )
    yMin = std::min(0.3*hPtAveSpecData_->GetBinContent(xMinBin),0.2*hPtAveSpecData_->GetBinContent(xMaxBin));
  else
    yMin = std::min(0.4*hPtAveSpecData_->GetBinContent(xMinBin),0.7*hPtAveSpecData_->GetBinContent(xMaxBin));
  yMin = std::max(yMin,0.01);	// Enable log scale
  xMinBin = std::max(1,xMinBin-10);
  xMaxBin = std::min(xMaxBin+10,hPtAveSpecMC_->GetNbinsX());

  // Absolute data spectrum
  TLegend* leg = util::LabelFactory::createLegendWithOffset(1,binLabel_->GetSize());
  leg->AddEntry(hPtAveSpecData_,DATA+",  N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
  hPtAveSpecData_->GetXaxis()->SetRange(xMinBin,xMaxBin);
  util::HistOps::setYRange(hPtAveSpecData_,6,3E-1);
  can->cd();
  hPtAveSpecData_->Draw("PE1");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"PtAveSpectrumData");

  // Absolute MC spectrum
  delete leg;
  leg = util::LabelFactory::createLegendWithOffset(1,binLabel_->GetSize());
  leg->AddEntry(hPtAveSpecMC_,MC+",  N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
  util::HistOps::setYRange(hPtAveSpecMC_,6,3E-5);
  can->cd();
  hPtAveSpecMC_->GetXaxis()->SetRange(xMinBin,xMaxBin);
  hPtAveSpecMC_->Draw("HIST");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"PtAveSpectrumMC");

  // Comparison (MC entries scaled to data)
  delete leg;
  leg = util::LabelFactory::createLegendWithOffset(2,binLabel_->GetSize());
  leg->AddEntry(hPtAveSpecData_,DATA+",  N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
  leg->AddEntry(hPtAveSpecMC_,MC+",  N = "+util::toTString(hPtAveSpecMC_->GetEntries()),"F");
  double scaleData = hPtAveSpecData_->Integral("width");
  double scaleMC = hPtAveSpecMC_->Integral("width");
  if( scaleData > 0. && scaleMC > 0. ) {
    hPtAveSpecMC_->Scale(scaleData/scaleMC); // Scale mc to data
    util::HistOps::setYRange(hPtAveSpecMC_,binLabel_->GetSize()+leg->GetNRows(),yMin);
    can->cd();
    hPtAveSpecMC_->Draw("HISTE");
    hPtAveSpecData_->Draw("PE1same");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    toFiles(can,outNameId+"PtAveSpectra");
    // set back MC scale
    hPtAveSpecMC_->Scale(scaleMC/scaleData);
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
  leg->AddEntry(hSymResp_,MCSMEAR,"P");
  leg->AddEntry(win,"Tail","F");
  if( DEBUG ) std::cout << "ok" << std::endl;

  if( DEBUG ) std::cout << "  Creating plots  . . .  " << std::flush;
  util::HistOps::setYRange(hSymResp_,4,3E-5);
  hSymResp_->GetXaxis()->SetRangeUser(0.,2.);

  TString canName = outNameId+"_SymMCTruthResponse";
  TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);
  can->cd();
  hSymResp_->Draw("HIST");
  //  fGaussMCTruth_->Draw("same");
  win->Draw("same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"_SymMCTruthResponseLog");

  // Linear scale
  util::HistOps::setYRange(hSymResp_,4);
  hSymResp_->GetXaxis()->SetRangeUser(0.4,1.6);
  can->cd();
  hSymResp_->Draw("HIST");
  //  fGaussMCTruth_->Draw("same");
  binLabel_->Draw("same");
  can->SetLogy(0);
  toFiles(can,outNameId+"_SymMCTruthResponse");

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
  TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);
  can->cd();
  hAsymMCSmeared_->Draw("HIST");
  win->Draw("same");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"_ToyMCTail");

  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);
  can->cd();
  hAsymMCSmeared_->Draw("HIST");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"_ToyMCLog");

  // Linear scale
  util::HistOps::setYRange(hAsymMCSmeared_,4);
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-0.4,0.4);
  can->cd();
  hAsymMCSmeared_->Draw("HIST");
  binLabel_->Draw("same");
  can->SetLogy(0);
  toFiles(can,outNameId+"_ToyMC");

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
  toFiles(can,outNameId+"_MCTruthResponseLog");

  // Linear scale
  util::HistOps::setYRange(hResp_,4);
  hResp_->GetXaxis()->SetRangeUser(0.4,1.6);
  can->cd();
  hResp_->Draw("HIST");
  binLabel_->Draw("same");
  can->SetLogy(0);
  toFiles(can,outNameId+"_MCTruthResponse");

  hResp_->GetXaxis()->SetRangeUser(0.,2.);


  delete win;
  delete leg;
  delete can;

  if( DEBUG ) std::cout << "Leaving Pt3Bin::plotToyAsymmetry()" << std::endl;
}





////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////


double coreScaleFactor(unsigned int etaBin, int variation) {
  double coreScaling = 0.;

//   // Core scale factors from runs 163337-167151
//   // Version 2011/07/20: closure uncertainty removed
//   if( variation == -1 ) {	// core down 1 sigma
//     // Uncertainties down
//     if( etaBin == 0 )      coreScaling = 0.000;
//     else if( etaBin == 1 ) coreScaling = 0.001;
//     else if( etaBin == 2 ) coreScaling = 0.032;
//     else if( etaBin == 3 ) coreScaling = 0.043;
//     else if( etaBin == 4 ) coreScaling = 0.090;
//   } else if( variation == 1 ) {	// core up 1 sigma
//     // Uncertainties up
//     if( etaBin == 0 )      coreScaling = 0.116;
//     else if( etaBin == 1 ) coreScaling = 0.115;
//     else if( etaBin == 2 ) coreScaling = 0.162;
//     else if( etaBin == 3 ) coreScaling = 0.228;
//     else if( etaBin == 4 ) coreScaling = 0.489;
//   } else if( variation == 0 ) {	// nominal core scale
//     // Nominal
//     if( etaBin == 0 )      coreScaling = 0.052;
//     else if( etaBin == 1 ) coreScaling = 0.057;
//     else if( etaBin == 2 ) coreScaling = 0.096;
//     else if( etaBin == 3 ) coreScaling = 0.134;
//     else if( etaBin == 4 ) coreScaling = 0.288;
//   } else {
//     std::cerr << "ERROR: unknown core-scale variation '" << variation << "'" << std::endl;
//     exit(1);
//   }


  // Core scale factors from runs 163337-167151
  // Version 2012/06/09 (thesis):
  // - corrected MC statistical uncertainties
  // - corrected extrapolation uncertainty determination
  // - corrected uncertainty propagation to ratio
  if( variation == -1 ) {	// core down 1 sigma
    // Uncertainties down
    if( etaBin == 0 )      coreScaling = 0.0000;
    else if( etaBin == 1 ) coreScaling = 0.0000;
    else if( etaBin == 2 ) coreScaling = 0.0192;
    else if( etaBin == 3 ) coreScaling = 0.0283;
    else if( etaBin == 4 ) coreScaling = 0.0737;
  } else if( variation == 1 ) {	// core up 1 sigma
    // Uncertainties up
    if( etaBin == 0 )      coreScaling = 0.1232;
    else if( etaBin == 1 ) coreScaling = 0.1254;
    else if( etaBin == 2 ) coreScaling = 0.1659;
    else if( etaBin == 3 ) coreScaling = 0.2567;
    else if( etaBin == 4 ) coreScaling = 0.5174;
  } else if( variation == 0 ) {	// nominal core scale
    // Nominal
    if( etaBin == 0 )      coreScaling = 0.0541;
    else if( etaBin == 1 ) coreScaling = 0.0599;
    else if( etaBin == 2 ) coreScaling = 0.0920;
    else if( etaBin == 3 ) coreScaling = 0.1409;
    else if( etaBin == 4 ) coreScaling = 0.2942;
  } else {
    std::cerr << "ERROR: unknown core-scale variation '" << variation << "'" << std::endl;
    exit(1);
  }

  return coreScaling;
}


void setStyleDataMarker(TH1* h) {
  h->SetMarkerSize(MARKER_SIZE);
  h->SetLineWidth(LINE_WIDTH);
  h->SetMarkerStyle(20);
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
}

void setStyleMCFilled(TH1* h) {
  h->SetLineWidth(1);
  h->SetMarkerStyle(1);
  TString name = h->GetName();
  if( name.Contains("Smear") )
    h->SetFillColor(COLOR_FILLED_ASYM_SMEAR);
  else 
    h->SetFillColor(COLOR_FILLED_ASYM);
}

void setStyleData(TGraphAsymmErrors* g) {
  g->SetMarkerSize(MARKER_SIZE);
  g->SetLineWidth(LINE_WIDTH);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(g->GetMarkerColor());
  g->SetLineStyle(2);
}

void setStyleMC(TGraphAsymmErrors* g) {
  g->SetMarkerSize(MARKER_SIZE);
  g->SetLineWidth(LINE_WIDTH);
  g->SetMarkerStyle(21);
  g->SetMarkerColor(COLOR_FILLED_ASYM_SMEAR);
  g->SetLineColor(g->GetMarkerColor());
}


// ------------------------------------------------------------------------------------
void printBin(unsigned int etaBin, unsigned int ptBin, const sampleTools::BinningAdmin* adm) {
  cout.setf(ios::fixed,ios::floatfield);
  std::cout << setprecision(1) << "    $" << adm->etaMin(etaBin) << "$ -- $" << adm->etaMax(etaBin) << "$ & $";
  double min = adm->ptMin(etaBin,ptBin);
  double max = adm->ptMax(etaBin,ptBin);
  if( min < 1000 ) std::cout << " ";
  if( min <  100 ) std::cout << " ";
  std::cout << setprecision(0) << min << " - ";
  if( max < 1000 ) std::cout << " ";
  if( max <  100 ) std::cout << " ";
  std::cout << adm->ptMax(etaBin,ptBin) << "$";
  std::cout << setprecision(5);
}


// ------------------------------------------------------------------------------------
void printWindowBorders(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm, double nSigTailStart) {

  std::cout << "\n\n\n***  Window Borders  ***\n\n";

  std::cout << "\\begin{tabular}{cr@{ -- }rcc}\n\\toprule\n";
  std::cout << "\\multicolumn{3}{c}{Interval} & \\multicolumn{2}{c}{\\tailborder{" << util::toTString(nSigTailStart,1,true) << "}} \\\\\n";
  std::cout << "$|\\eta|$ & \\multicolumn{2}{c}{$\\ptave \\,(\\gevnospace)$} & \\asymtaileff & $\\asymtaileff/\\sigma_{\\asym}$ \\\\\n\\midrule\n";

  for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
    EtaPtBin* bin = *it;
    double min = bin->tailWindowMinEff();
    double sigSmear = bin->sigmaSmeared(0);
    
    printBin(bin->etaBin(),bin->ptBin(),adm);
    std::cout << std::setprecision(3) << " & $" << min << "$ & $" << min/sigSmear << "$ \\\\\n";
    if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\midrule\n";
  }
  std::cout << "\\bottomrule\n\\end{tabular}\n\n";


  //  std::cout << "\n\n\n***  Tail vs Gaussian  ***\n\n";

//   std::cout << "\\begin{tabular}[ht]{cccc}\n\\toprule\n";
//   std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\int^{";
//   if( exclRegion ) std::cout << "A_{1}";
//   else std::cout << "\\infty";
//   std::cout << "}_{A_{0}}\\mathcal{G}$ & $\\fasymmc(\\pti{3}/\\ptave = 0.05)$ \\\\ \n\\midrule\n";
//   for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
//     EtaPtBin* bin = *it;
//     double min = bin->tailWindowMin();
//     double max = bin->tailWindowMax();
//     double sigSmear = bin->sigmaSmeared(0);
//     double fGauss = (erf(max/sqrt(2.)/sigSmear) - erf(min/sqrt(2.)/sigSmear))/2.;
//     double fAsym = bin->fTailMCSmeared(0);
//     double fAsymErr = bin->fTailMCSmearedErr(0);

//     printBin(bin->etaBin(),bin->ptBin(),adm);
//     std::cout << std::setprecision(3) << " & $";
//     std::cout << std::setprecision(4) << fGauss << "$ & $" << fAsym << " \\pm " << fAsymErr << "$ \\\\" << std::endl;
//     if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\midrule\n";
//   }
//   std::cout << "  \\end{tabular}\n\n";
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
  std::cout << "\\begin{tabular}[ht]{ccccc}\n\\toprule\n";
  std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\mean{\\ptave} \\,(\\gevnospace)$ & $\\fasymdata(0)$ & $\\fasymmc(0)$ \\\\ \n\\midrule\n";
  for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
    EtaPtBin* bin = *it;
    printBin(bin->etaBin(),bin->ptBin(),adm);
    std::cout << setprecision(1);
    std::cout << " & $"  << bin->ptAveMeanData(0) << " \\pm " << bin->ptAveMeanDataErr(0);
    std::cout << setprecision(3);
    std::cout << "$ & $"  << bin->extraData() << " \\pm " << bin->extraDataErr();  
    std::cout << "$ & $" << bin->extraMC() << " \\pm " << bin->extraMCErr() << "$ \\\\ \n";
    if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\midrule\n";
  }
  std::cout << "  \\end{tabular}\n\n";
}


// ------------------------------------------------------------------------------------
bool completePicName(int padIdx, const sampleTools::BinningAdmin* adm, TString &picName) {
  bool picExists = true;
  if( padIdx == 0 )
    picName += "_Extrapolation.pdf";
  else if( padIdx-1 < static_cast<int>(adm->nPtSoftBins()) )
    picName += "_Pt3Bin"+util::toTString(padIdx-1)+"_PtSmearAsymTail.pdf";
  else if( padIdx-1-adm->nPtSoftBins() < adm->nPtSoftBins() )
    picName += "_Pt3Bin"+util::toTString(padIdx-1-adm->nPtSoftBins())+"_PtAveSpectra.pdf";
  else
    picExists = false;

  return picExists;
}



// ------------------------------------------------------------------------------------
void writeLaTeXSlides(const TString &outNameId, const sampleTools::BinningAdmin* adm) {
  ofstream outFile(outNameId+".tex");
  if( outFile.is_open() ) {
    for(unsigned int etaBin = 0; etaBin < adm->nEtaBins(); ++etaBin) {
      for(unsigned int ptBin = 0; ptBin < adm->nPtBins(etaBin); ++ptBin) {
	outFile << "\n\n\n% ==== Eta " << etaBin << ", Pt " << ptBin << " ================================================\n";

	int nPages = 1 + (1+2*adm->nPtSoftBins())%12;
	for(int page = 0; page < nPages; ++page) {
	  outFile << "\n% -----------------------------------------------------------------\n";
	  outFile << "\\begin{frame}\\frametitle{Asymmetry Tails $" << adm->etaMin(etaBin) << " < |\\eta| < " << adm->etaMax(etaBin) << "$}\n";
	  outFile << "  \\begin{columns}[T]\n";
	  for(int colIdx = 0; colIdx < 4; ++colIdx) {
	    outFile << "    \\begin{column}{0.25\\textwidth}\n";
	    outFile << "    \\centering\n";
	    for(int rowIdx = 0; rowIdx < 3; ++rowIdx) {
	      int padIdx = 12*page + 4*rowIdx + colIdx;
	      TString picName = outNameId+"_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
	      if( completePicName(padIdx,adm,picName) ) {
		outFile << "      \\includegraphics[width=\\textwidth]{figures/" << picName << "}\\\\\n";
	      } else {
		continue;
	      }
	    }
	    outFile << "    \\end{column}\n";
	  }
	  outFile << "  \\end{columns}\n";
	  outFile << "\\end{frame}\n";
	}
      }
    }
    outFile.close();
  } else {
    std::cerr << "ERROR writing LaTeX slides to file '" << outNameId << ".tex'" << std::endl;
  }
}


// ------------------------------------------------------------------------------------
TString labelWindow(double nSigMin, double nSigMax) {
  TString label = "";
  if( nSigMax < 50. ) 
    label = util::toTString(nSigMin)+" - "+util::toTString(nSigMax);
  else
    label = ASYM_TAIL_START+" = "+util::toTString(nSigMin);
  label += " "+COMMON_SIGMA;

  return label;
}


// ------------------------------------------------------------------------------------
double getTailStart(const TString &name) {
  double tailStart = 0.;
  if( name.Contains("Sig15") ) {
    tailStart = 1.5;
  } else if( name.Contains("Sig20") ) {
    tailStart = 2.0;
  } else if( name.Contains("Sig25") ) {
    tailStart = 2.5;
  } else if( name.Contains("Sig30") ) {
    tailStart = 3.0;
  } else if( name.Contains("Sig35") ) {
    tailStart = 3.5;
  }

  return tailStart;
}


// ------------------------------------------------------------------------------------
void correctAsymmetryWidth(TH1* hOrig, double nSigCore, double coreScale, TH1* &hSmeared, double &width, double &smearedWidth) {
  // Fit width of original distribution
  width = 0.;
  double widthErr = 1000.;
  if( !util::HistOps::fitCoreWidth(hOrig,nSigCore,width,widthErr) ) {
    width = hOrig->GetRMS();
  }
  if( width > 2.*hOrig->GetRMS() ) {
    width = hOrig->GetRMS();
  }
  // Smear original distribution
  util::HistOps::smearHistogram(hOrig,hSmeared,width,coreScale);
  // Width of correced distribution
  if( !util::HistOps::fitCoreWidth(hSmeared,nSigCore,smearedWidth,widthErr) ) {
    smearedWidth = hSmeared->GetRMS();
  }
  if( smearedWidth > 2.*hSmeared->GetRMS() ) {
    smearedWidth = hSmeared->GetRMS();
  }
}


// ------------------------------------------------------------------------------------
bool hugeWidth(const TH1* h, double &width) {
  bool isHuge = false;
  if( width > 2.*h->GetRMS() ) {
    std::cerr << ">>>> HUGE WIDTH" << std::endl;
    width = h->GetRMS();
    isHuge = true;
  }

  return isHuge;
}



// ------------------------------------------------------------------------------------
void toFiles(TNamed* obj, const TString &name) {
  if( ROOT_OUTPUT && ROOT_OUT_FILE ) {
    TString objName = obj->GetName();
    TString objTitle = obj->GetTitle();
    if( name != "" ) {
      obj->SetName(name);
      obj->SetTitle(name);
    }
    ROOT_OUT_FILE->WriteTObject(obj);
    obj->SetName(objName);
    obj->SetTitle(objTitle);
  }
}

void toFiles(TCanvas* obj, const TString &name) {
  if( EPS_OUTPUT ) {
    TString objName = obj->GetName();
    TString objTitle = obj->GetTitle();
    if( name != "" ) {
      obj->SetName(name);
      obj->SetTitle(name);
    }
    obj->SaveAs(name+".eps","eps");
    if( ROOT_OUTPUT && ROOT_OUT_FILE ) ROOT_OUT_FILE->WriteTObject(obj);
    obj->SetName(objName);
    obj->SetTitle(objTitle);
  }
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
  g->SetMarkerSize(MARKER_SIZE);
  g->SetLineWidth(LINE_WIDTH);

  return g;
}


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* relUncertainty(const TH1* hNom, int color, double weight, const TH1* hVarUp, const TH1* hVarDn = 0) {

  std::vector<double> x;
  std::vector<double> xe;
  std::vector<double> y;
  std::vector<double> yeu;
  std::vector<double> yed;
  for(int bin = 1; bin <= hNom->GetNbinsX(); ++bin) {
    x.push_back(hNom->GetBinCenter(bin));
    xe.push_back(0.5*hNom->GetBinWidth(bin));
    y.push_back(0.);
    double nom = hNom->GetBinContent(bin);
    double var = hVarUp->GetBinContent(bin);
    double err = std::abs((var-nom)/nom);
    if( hVarDn ) {
      var = hVarDn->GetBinContent(bin);
      err += std::abs((var-nom)/nom);
      err /= 2.;
    }
    err *= weight;
    yeu.push_back(err);
    yed.push_back(err);
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
					       &(xe.front()),&(xe.front()),&(yed.front()),&(yeu.front()));
  g->SetFillStyle(1001);
  g->SetFillColor(color);
  g->SetLineColor(color);

  return g;
}


// ------------------------------------------------------------------------------------
TGraphAsymmErrors* totalUncertainty(const std::vector<TGraphAsymmErrors*> &uncerts) {

  TGraphAsymmErrors* g = static_cast<TGraphAsymmErrors*>(uncerts.at(0)->Clone());
  for(int i = 0; i < g->GetN(); ++i) {
    double eu2 = 0.;
    double ed2 = 0.;
    for(unsigned int j = 0; j < uncerts.size(); ++j) {
      eu2 += pow(uncerts.at(j)->GetEYhigh()[i],2); 
      ed2 += pow(uncerts.at(j)->GetEYlow()[i],2);
    }
    if( ed2 == 0. ) ed2 = eu2;
    g->SetPointError(i,g->GetEXlow()[i],g->GetEXhigh()[i],sqrt(ed2),sqrt(eu2));
  }
  g->SetFillStyle(3444);
  g->SetFillColor(kBlack);
  g->SetLineColor(kBlack);


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
void plotFinalResult(const TString &fileName) {
  if( SHOW_HEADER ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStyleNoteNoTitle();
  }
  gErrorIgnoreLevel = 1001;

  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");  
  const TString fileNamePrefix = util::baseName(fileName);
  const TString outNamePrefix = util::fileName(fileNamePrefix+"_ScaleFactors");


  // ROOT output
  ROOT_OUT_FILE = new TFile(outNamePrefix+".root","RECREATE");

  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();


  TH1* hFrame = 0;
//   std::vector<TGraphAsymmErrors*> vgUncertAbs;
//   std::vector<TGraphAsymmErrors*> vgNomRatio;
//   std::vector< std::vector<TGraphAsymmErrors*> >


  util::MultiCanvas* mcScales = new util::MultiCanvas(outNamePrefix,3,2,5,false);
  util::MultiCanvas* mcUncert = new util::MultiCanvas(outNamePrefix+"_Uncertainties",3,2,5,false);
	
//   TLegend* leg = util::LabelFactory::createLegendWithOffset(3,0.2+util::LabelFactory::lineHeightMultiCan(),util::LabelFactory::lineHeightMultiCan());
//   mc->adjustLegend(leg);
//   mc->markForDeletion(leg);
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    // ****** Nominal scaling factors *****************************************************
    
    // Read mean pt
    TString histName = "hPtAveMeanData_Eta"+util::toTString(etaBin);
    TH1* hPtAveMeanData = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName);

    // Read nominal scaling factors
    histName = "hScaleFrame_Eta"+util::toTString(etaBin);
    TH1* hScaleFrame = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName);
    hScaleFrame->GetXaxis()->SetTitle("p^{ave}_{T} (GeV)");
    hScaleFrame->SetLineWidth(LINE_WIDTH);
    hScaleFrame->SetLineStyle(2);
    histName = "hScale_Eta"+util::toTString(etaBin);
    TH1* hScaleNom = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName+"_Nominal");

    // Graph of nominal scaling factors with ptAveMean as x
    TGraphAsymmErrors* gNomRatio = nomRatio(hScaleNom,hPtAveMeanData);


    // ****** Variations for systematic uncertainties *************************************

    // Read histograms with variations
    TH1* hScaleVarCoreUp = util::FileOps::readTH1(fileNamePrefix+"_VarCoreUp.root",histName,histName+"_VarCoreUp");
    TH1* hScaleVarCoreDn = util::FileOps::readTH1(fileNamePrefix+"_VarCoreDn.root",histName,histName+"_VarCoreDn");
    TH1* hScaleVarClosure = util::FileOps::readTH1(fileNamePrefix+"_VarClosure.root",histName,histName+"_VarClosure");
    TH1* hScaleVarPUUp = util::FileOps::readTH1(fileNamePrefix+"_VarPUUp.root",histName,histName+"_VarPUUp");
    TH1* hScaleVarPUDn = util::FileOps::readTH1(fileNamePrefix+"_VarPUDn.root",histName,histName+"_VarPUDn");
    TH1* hScaleVarExtra = util::FileOps::readTH1(fileNamePrefix+"_VarExtrapolation.root",histName,histName+"_VarExtrapolation");

    // Get relative uncertainties
    std::vector<TGraphAsymmErrors*> uncerts;
    uncerts.push_back(relUncertainty(hScaleNom,11,1.,hScaleVarPUUp,hScaleVarPUDn));
    uncerts.push_back(relUncertainty(hScaleNom,38,1.,hScaleVarExtra));
    uncerts.push_back(relUncertainty(hScaleNom,8,0.5,hScaleVarClosure));
    uncerts.push_back(relUncertainty(hScaleNom,46,1.,hScaleVarCoreUp,hScaleVarCoreDn));

    // Define labels
    std::vector<TString> uncertLabels;
    uncertLabels.push_back("PU");
    uncertLabels.push_back("Additional Jets");
    uncertLabels.push_back("Bias");
    uncertLabels.push_back("#sigma_{A} Correction");

    // Add up uncertainties
    TGraphAsymmErrors* gUncertRelTotal = totalUncertainty(uncerts);
    TGraphAsymmErrors* gUncertAbs = uncertaintyBand(hScaleNom,gUncertRelTotal,5);
    TGraphAsymmErrors* gUncertAbsPos = static_cast<TGraphAsymmErrors*>(gUncertAbs->Clone());
    for(int i = 0; i < gUncertAbsPos->GetN(); ++i) {
      gUncertAbsPos->SetPointError(i,gUncertAbsPos->GetEXlow()[i],gUncertAbsPos->GetEXhigh()[i],
				   gUncertAbsPos->GetEYlow()[i]>gUncertAbsPos->GetY()[i] ? gUncertAbsPos->GetY()[i] : gUncertAbsPos->GetEYlow()[i],gUncertAbsPos->GetEYhigh()[i]);
    }

    // Stacked (quadratic addition) uncertainties for plot of relative systematics
    // Only positive
    std::vector<TGraphAsymmErrors*> gUncertStack;
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      std::vector<double> x;
      std::vector<double> xed;
      std::vector<double> xeu;
      std::vector<double> y;
      std::vector<double> yed;
      std::vector<double> yeu;
      for(int n = 0; n < uncerts.at(i)->GetN(); ++n) {
	x.push_back(uncerts.at(i)->GetX()[n]);
	xed.push_back(uncerts.at(i)->GetEXlow()[n]);
	xeu.push_back(uncerts.at(i)->GetEXhigh()[n]);
	y.push_back(uncerts.at(i)->GetY()[n]);
	yed.push_back(0.);
	double err = 100.*uncerts.at(i)->GetEYhigh()[n];
	if( i > 0 ) {
	  double prevErr = gUncertStack.back()->GetEYhigh()[n];
	  err = sqrt( err*err + prevErr*prevErr );
	}
	yeu.push_back(err);
      }
      x.insert(x.begin(),binAdm->ptMin(etaBin,0));
      xed.insert(xed.begin(),0.);
      xeu.insert(xeu.begin(),0.);
      y.insert(y.begin(),y.front());
      yed.insert(yed.begin(),0.);
      yeu.insert(yeu.begin(),yeu.front());

      x.push_back(binAdm->ptMax(etaBin,binAdm->nPtBins(etaBin)-1));
      xed.push_back(0.);
      xeu.push_back(0.);
      y.push_back(y.back());
      yed.push_back(0.);
      yeu.push_back(yeu.back());

      TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),&(xed.front()),&(xeu.front()),&(yed.front()),&(yeu.front()));
      g->SetLineWidth(0);
      g->SetFillStyle(1001);
      g->SetFillColor(uncerts.at(i)->GetFillColor());
      g->SetLineColor(g->GetFillColor());
      gUncertStack.push_back(g);
    }


    // ****** Plotting ********************************************************************

    // Label
    TPaveText* label = util::LabelFactory::createPaveText(1);
    label->AddText(LUMI+",  "+labelWindow(getTailStart(fileNamePrefix),1000.)+",  "+util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));

    TLegend* legScale = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
    legScale->AddEntry(gNomRatio,"Extrapolation Uncertainty #deltaf_{ex}","L");
    legScale->AddEntry(gUncertAbs,"Systematic Uncertainty","F");

    // Plot scaling factors and total uncertainty    
    hScaleFrame->SetTitle(HEADER);
    hScaleFrame->GetYaxis()->SetRangeUser(0.01,3.3);
    hScaleFrame->GetYaxis()->SetTitle(SCALE);
    hScaleFrame->GetXaxis()->SetMoreLogLabels();
    hScaleFrame->GetXaxis()->SetNoExponent();
    TCanvas* canScale = util::HistOps::createTCanvas(outNamePrefix+"_canScale_"+util::toTString(etaBin),"Scale Factors Eta "+util::toTString(etaBin),500,500);
    canScale->cd();
    hScaleFrame->Draw("HIST");
    gUncertAbs->Draw("E2same");
    hScaleFrame->Draw("HISTsame");
    gNomRatio->Draw("PE1same");
    label->Draw("same");
    legScale->Draw("same");
    canScale->SetLogx();
    gPad->RedrawAxis();
    toFiles(canScale,outNamePrefix+"_Eta"+util::toTString(etaBin));

    // Multicanvas
    if( etaBin == 0 ) {
      hFrame = static_cast<TH1*>(hScaleFrame->Clone("FrameForMultiCan"));
      mcScales->markForDeletion(hFrame);
    }
    TH1* mFrame = mcScales->mainFrame(hFrame,etaBin);
    mcScales->canvas()->cd();
    TPad* padt = mcScales->mainPad(etaBin);
    padt->Draw();
    padt->cd();
    mFrame->Draw("HIST");
    gUncertAbs->Draw("E2same");
    mFrame->Draw("HISTsame");
    gNomRatio->Draw("PE1same");
    gPad->RedrawAxis();

    // Plot relative uncertainties
    TPaveText* labelU = util::LabelFactory::createPaveText(2);
    labelU->AddText(MCSMEAR);
    labelU->AddText(labelWindow(getTailStart(fileNamePrefix),1000.)+",   "+util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));

    TLegend* legUncert1 = util::LabelFactory::createLegendColWithOffset(2,-0.5,labelU->GetSize());
    TLegend* legUncert2 = util::LabelFactory::createLegendColWithOffset(2,0.5,labelU->GetSize());
    for(int i = static_cast<int>(gUncertStack.size())-1; i >= 0; --i) {
      if( i > 1 ) legUncert1->AddEntry(gUncertStack.at(i),uncertLabels.at(i),"F");
      else        legUncert2->AddEntry(gUncertStack.at(i),uncertLabels.at(i),"F");
    }

    TText* line = label->GetLine(0);
    line->SetText(line->GetX(),line->GetY(),MCSMEAR);
    TH1* hUncertsFrame = static_cast<TH1D*>(hScaleFrame->Clone("hUncertsFrame"+util::toTString(etaBin)));
    for(int bin = 1; bin <= hUncertsFrame->GetNbinsX(); ++bin) {
      hUncertsFrame->SetBinContent(bin,-1.);
    }
    hUncertsFrame->GetYaxis()->SetRangeUser(0.,95);
    hUncertsFrame->GetXaxis()->SetMoreLogLabels();
    hUncertsFrame->GetXaxis()->SetNoExponent();
    hUncertsFrame->GetYaxis()->SetTitle("Relative Uncertainty on "+SCALE+"  (%)");
    TCanvas* canRelUncerts = util::HistOps::createTCanvas(outNamePrefix+"_canRelUncerts_"+util::toTString(etaBin),"Relative Uncertainties Eta "+util::toTString(etaBin),500,500);
    canRelUncerts->cd();
    hUncertsFrame->Draw("HIST");
    for(std::vector<TGraphAsymmErrors*>::reverse_iterator it = gUncertStack.rbegin();
	it != gUncertStack.rend(); ++it) {
      (*it)->Draw("E3same");
    }
    labelU->Draw("same");
    legUncert1->Draw("same");
    legUncert2->Draw("same");
    canRelUncerts->SetLogx();
    gPad->RedrawAxis();
    toFiles(canRelUncerts,outNamePrefix+"_Uncertainties_Eta"+util::toTString(etaBin));

    
    // store factors and total uncertainties for evolution plots
    TH1* hOutStat = new TH1D("ScaleFactors_StatUncert_Eta"+util::toTString(etaBin),";PtAveBin;ScaleFactor",gNomRatio->GetN(),0,gNomRatio->GetN());
    TH1* hOutSyst = new TH1D("ScaleFactors_SystUncert_Eta"+util::toTString(etaBin),";PtAveBin;ScaleFactor",gNomRatio->GetN(),0,gNomRatio->GetN());
    for(int i = 0; i < gNomRatio->GetN(); ++i) {
      hOutStat->SetBinContent(1+i,gNomRatio->GetY()[i]);
      hOutStat->SetBinError(1+i,gNomRatio->GetEYhigh()[i]);
      hOutSyst->SetBinContent(1+i,gNomRatio->GetY()[i]);
      hOutSyst->SetBinError(1+i,gUncertAbs->GetEYhigh()[i]);
    }
    toFiles(hOutStat);
    toFiles(hOutSyst);
    delete hOutStat;
    delete hOutSyst;
    

    // ****** Print results ***************************************************************

    // Print factors with total uncertainty (stat + syst)
    if( etaBin == 0 ) {
      std::cout << "\\renewcommand{\\arraystretch}{1.5}\n";
      std::cout << "\\begin{tabular}{cr@{ -- }rcc}\n\\toprule\n";
      std::cout << "$|\\eta|$ & \\multicolumn{2}{c}{$\\ptave\\,(\\gevnospace)$} & $\\mean{\\ptave}\\,(\\gevnospace)$ & $\\tailratio(\\tailborder{" << getTailStart(fileNamePrefix) << "})$ \\\\\n\\midrule\n";
    }
    for(int bin = 1; bin <= hScaleNom->GetNbinsX(); ++bin) {
      cout.setf(ios::fixed,ios::floatfield);
      std::cout << std::setprecision(1) << "  $" << util::toTString(binAdm->etaMin(etaBin),1) << "$ -- $" << util::toTString(binAdm->etaMax(etaBin),1) << "$ & $";
      std::cout << std::setprecision(0);
      if( binAdm->ptMin(etaBin,bin-1) < 1000 ) {
	std::cout << " ";
	if( binAdm->ptMin(etaBin,bin-1) < 100 ) {
	  std::cout << " ";
	}
      }
      std::cout << util::toTString(binAdm->ptMin(etaBin,bin-1)) << "$ & $";
      if( binAdm->ptMax(etaBin,bin-1) < 1000 ) {
	std::cout << " ";
	if( binAdm->ptMax(etaBin,bin-1) < 100 ) {
	  std::cout << " ";
	}
      }
      std::cout << util::toTString(binAdm->ptMax(etaBin,bin-1)) << "$ & $";
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
    std::cout << "\\midrule\n";
    if( etaBin == binAdm->nEtaBins()-1 ) {
      std::cout << "\\bottomrule\n\\end{tabular}\n\\renewcommand{\\arraystretch}{1.}\n";
    }
  }

  mcScales->moreLogLabelsX();
  mcScales->noExponentX();
  mcScales->setLogx();
  toFiles(mcScales->canvas(),outNamePrefix);
  delete mcScales;


  ROOT_OUT_FILE->Close();
  delete ROOT_OUT_FILE;

  delete binAdm;
}



// ------------------------------------------------------------------------------------
void plotEvolution() {
  if( SHOW_HEADER ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStyleNoteNoTitle();
  }
  gErrorIgnoreLevel = 1001;

  // Evolution with NVtx
  std::vector<TString> fileNames;
  fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_NVtx00-06_Sig20-Inf_PF.root");
  fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_NVtx07-99_Sig20-Inf_PF.root");
  const TString outNamePrefix = "Tail_163337-180252_NVtxEvolution_Sig20_PF";
  const TString labWin = labelWindow(2.,1000.);

  // Binning
  std::vector<double> xBinEdges;
  for(unsigned int f = 0; f <= fileNames.size(); ++f) {
    xBinEdges.push_back(f);
  }
  std::vector<TString> xBinLabels;
  xBinLabels.push_back("#leq 6");
  xBinLabels.push_back("#geq 7");

  assert( xBinLabels.size()+1 == xBinEdges.size() );

  // ROOT output
  ROOT_OUT_FILE = new TFile(outNamePrefix+".root","RECREATE");
  EPS_OUTPUT = true;

  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");  
  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  
  
  // Read histograms      
  std::vector< std::vector< std::vector<TH1*> > > hFAsymData;
  std::vector< std::vector< std::vector<TH1*> > > hFAsymMC;
  std::vector< std::vector< std::vector<TH1*> > > hFAsymGauss;
  std::vector<double> yMin(binAdm->nPtSoftBins(),1000.);	// per pt3Bin
  std::vector<double> yMax(binAdm->nPtSoftBins(),0.);	// per pt3Bin  
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    
    std::vector< std::vector<TH1*> > hFAsymData1;
    std::vector< std::vector<TH1*> > hFAsymMC1;
    std::vector< std::vector<TH1*> > hFAsymGauss1;
    for(unsigned int pt3Bin = 0; pt3Bin < binAdm->nPtSoftBins(); ++pt3Bin) {
      std::vector<TH1*> hFAsymData2;
      std::vector<TH1*> hFAsymMC2;
      std::vector<TH1*> hFAsymGauss2;
      for(unsigned int f = 0; f < fileNames.size(); ++f) {
	const TString binId = "_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin);
	hFAsymData2.push_back(util::FileOps::readTH1(fileNames.at(f),"hFAsymData"+binId,"hFAsymData"+binId+"_"+util::toTString(f),false));
	hFAsymMC2.push_back(util::FileOps::readTH1(fileNames.at(f),"hFAsymMCSmeared"+binId,"hFAsymMCSmeared"+binId+"_"+util::toTString(f),false));
	hFAsymGauss2.push_back(util::FileOps::readTH1(fileNames.at(f),"hFAsymMCSmearedGauss"+binId,"hFAsymMCSmearedGauss"+binId+"_"+util::toTString(f),false));
	
	double min = 0.;
	double max = 0.;
	util::HistOps::findYRange(hFAsymMC2.back(),min,max);
	if( min < yMin.at(pt3Bin) ) yMin.at(pt3Bin) = min;
	util::HistOps::findYRange(hFAsymData2.back(),min,max);
	if( max > yMax.at(pt3Bin) ) yMax.at(pt3Bin) = max;
      } // end of loop over files
      hFAsymData1.push_back(hFAsymData2);
      hFAsymMC1.push_back(hFAsymMC2);
      hFAsymGauss1.push_back(hFAsymGauss2);
    } // end of loop over pt3 bins
    hFAsymData.push_back(hFAsymData1);
    hFAsymMC.push_back(hFAsymMC1);
    hFAsymGauss.push_back(hFAsymGauss1);
  } // end of loop over eta bins
    
    // Plot evolution
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for(unsigned int pt3Bin = 0; pt3Bin < binAdm->nPtSoftBins(); ++pt3Bin) {
      for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
	const TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+"_Pt3Bin"+util::toTString(pt3Bin);
	TH1* hEvoFAsymData = new TH1D("hEvolutionOfFAsymData"+binId,HEADER+";"+LABEL_NVTX,fileNames.size(),&(xBinEdges.front()));
	for(int bin = 1; bin <= hEvoFAsymData->GetNbinsX(); ++bin) {
	  hEvoFAsymData->GetXaxis()->SetBinLabel(bin,xBinLabels.at(bin-1));
	}
	hEvoFAsymData->SetLabelSize(0.07,"X");
	hEvoFAsymData->GetYaxis()->SetTitle(hFAsymData.at(etaBin).at(pt3Bin).front()->GetYaxis()->GetTitle());
	hEvoFAsymData->SetMarkerStyle(hFAsymData.at(etaBin).at(pt3Bin).front()->GetMarkerStyle());
	hEvoFAsymData->SetMarkerColor(hFAsymData.at(etaBin).at(pt3Bin).front()->GetMarkerColor());
	hEvoFAsymData->SetMarkerSize(hFAsymData.at(etaBin).at(pt3Bin).front()->GetMarkerSize());
	hEvoFAsymData->SetLineWidth(hFAsymData.at(etaBin).at(pt3Bin).front()->GetLineWidth());
	  
	TH1* hEvoFAsymMC = static_cast<TH1*>(hEvoFAsymData->Clone("hEvolutionOfFAsymMCSmeared"+binId));
	hEvoFAsymMC->GetYaxis()->SetTitle(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetYaxis()->GetTitle());
	hEvoFAsymMC->SetMarkerStyle(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetMarkerStyle());
	hEvoFAsymMC->SetMarkerColor(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetMarkerColor());
	hEvoFAsymMC->SetLineStyle(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetLineStyle());
	hEvoFAsymMC->SetLineColor(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetLineColor());
	hEvoFAsymMC->SetFillColor(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetFillColor());
	hEvoFAsymMC->SetFillStyle(hFAsymMC.at(etaBin).at(pt3Bin).front()->GetFillStyle());
	hEvoFAsymMC->SetLineColor(COLOR_FILLED_ASYM_SMEAR);
	
	TH1* hEvoFAsymGauss = static_cast<TH1*>(hEvoFAsymData->Clone("hEvolutionOfFAsymMCSmearedGauss"+binId));
	hEvoFAsymGauss->GetYaxis()->SetTitle(hFAsymGauss.at(etaBin).at(pt3Bin).front()->GetYaxis()->GetTitle());
	hEvoFAsymGauss->SetMarkerStyle(hFAsymGauss.at(etaBin).at(pt3Bin).front()->GetMarkerStyle());
	hEvoFAsymGauss->SetMarkerColor(hFAsymGauss.at(etaBin).at(pt3Bin).front()->GetMarkerColor());
	hEvoFAsymGauss->SetLineStyle(hFAsymGauss.at(etaBin).at(pt3Bin).front()->GetLineStyle());
	hEvoFAsymGauss->SetLineColor(hFAsymGauss.at(etaBin).at(pt3Bin).front()->GetLineColor());

	// Fill plots
	for(unsigned int f = 0; f < fileNames.size(); ++f) {
	  int binEvo = 1+f;
	  int binIn = 1+ptBin;
	  hEvoFAsymData->SetBinContent(binEvo,hFAsymData.at(etaBin).at(pt3Bin).at(f)->GetBinContent(binIn));
	  hEvoFAsymData->SetBinError(binEvo,hFAsymData.at(etaBin).at(pt3Bin).at(f)->GetBinError(binIn));
	  hEvoFAsymMC->SetBinContent(binEvo,hFAsymMC.at(etaBin).at(pt3Bin).at(f)->GetBinContent(binIn));
	  hEvoFAsymMC->SetBinError(binEvo,hFAsymMC.at(etaBin).at(pt3Bin).at(f)->GetBinError(binIn));
	  hEvoFAsymGauss->SetBinContent(binEvo,hFAsymGauss.at(etaBin).at(pt3Bin).at(f)->GetBinContent(binIn));
	  hEvoFAsymGauss->SetBinError(binEvo,hFAsymGauss.at(etaBin).at(pt3Bin).at(f)->GetBinError(binIn));
	}

	// Convert hEvoFAsymMC to graph to be able to plot error band
	TGraphAsymmErrors* gEvoFAsymMC = util::HistOps::getUncertaintyBand(hEvoFAsymMC,COLOR_FILLED_ASYM_SMEAR);      
	gEvoFAsymMC->SetFillStyle(HATCH_STYLE);
	
	// Label
	TPaveText* label = util::LabelFactory::createPaveText(2);
	label->AddText(util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin))+",  "+util::LabelFactory::ptAveCut(binAdm->ptMin(etaBin,ptBin),binAdm->ptMax(etaBin,ptBin)));
	label->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+util::LabelFactory::pt3RelCut(binAdm->ptSoftMax(pt3Bin))+",  "+labWin);
		       
	TLegend* leg = util::LabelFactory::createLegendWithOffset(2,label->GetSize());
	leg->AddEntry(hEvoFAsymData,DATA,"P");
	leg->AddEntry(gEvoFAsymMC,MCSMEAR,"LF");

	double min = 0.6*std::min(hEvoFAsymMC->GetBinContent(hEvoFAsymMC->GetMinimumBin()),
				  hEvoFAsymData->GetBinContent(hEvoFAsymData->GetMinimumBin()));
	double max = 1.10*std::max(hEvoFAsymMC->GetBinContent(hEvoFAsymMC->GetMaximumBin()),
			      hEvoFAsymData->GetBinContent(hEvoFAsymData->GetMaximumBin()));
 	double delta = max - min;
 	hEvoFAsymMC->GetYaxis()->SetRangeUser(min,max+delta);
	TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
	can->cd();
	hEvoFAsymMC->Draw("HIST");
	gEvoFAsymMC->Draw("E2same");
	hEvoFAsymMC->Draw("HISTsame");
	hEvoFAsymData->Draw("PE1same");
	label->Draw("same");
	leg->Draw("same");
	gPad->RedrawAxis();
	toFiles(can,outNamePrefix+binId+"_FAsym");
	
	delete hEvoFAsymData;
	delete hEvoFAsymMC;
	delete hEvoFAsymGauss;
	delete gEvoFAsymMC;
	delete label;
	delete leg;
	delete can;
      } // End of loop over pt bins
    }// End of loop over pt3 bins
  } // End of loop over eta bins
  
  for(unsigned int i = 0; i < hFAsymData.size(); ++i) {
    for(unsigned int j = 0; j < hFAsymData.at(i).size(); ++j) {
      for(unsigned int k = 0; k < hFAsymData.at(i).at(j).size(); ++k) {
	delete hFAsymData.at(i).at(j).at(k);
	delete hFAsymMC.at(i).at(j).at(k);
	delete hFAsymGauss.at(i).at(j).at(k);
      }
    }
  }
  
  ROOT_OUT_FILE->Close();
  delete ROOT_OUT_FILE;
}



// ------------------------------------------------------------------------------------
void plotEvolutionFinalPlots() {
  if( SHOW_HEADER ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStyleNoteNoTitle();
  }
  gErrorIgnoreLevel = 1001;

  std::vector<TString> fileNames;

//   // Evolution with tail region
//   fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_Sig20-Inf_PF_ScaleFactors.root");
//   fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_Sig25-Inf_PF_ScaleFactors.root");
//   fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_Sig30-Inf_PF_ScaleFactors.root");
//   const TString xAxisLabel = "Tail Start "+ASYM_TAIL_START+" ("+COMMON_SIGMA+")";
//   const int nBins = 3;
//   const double xBinEdges[nBins+1] = { 1.75, 2.25, 2.75, 3.25 };
//   const TString binLabels[1] = { "" };
//   const TString outNamePrefix = "Tail_163337-180252_SigEvolution_PF_ScaleFactors";

  // Evolution with NVtx
  fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_NVtx00-06_Sig30-Inf_PF_ScaleFactors.root");
  fileNames.push_back("ScaleFactors_163337-180252_2012-06-12/Tail_163337-180252_NVtx07-99_Sig30-Inf_PF_ScaleFactors.root");
  const TString xAxisLabel = LABEL_NVTX;
  const int nBins = 2;
  const double xBinEdges[nBins+1] = { 0., 1., 2. };
  const TString binLabels[nBins+1] = { "#leq 6", "#geq 7" };
  const TString outNamePrefix = "Tail_163337-180252_NVtxEvolution_Sig30_PF_ScaleFactors";
  const TString labWin = labelWindow(3.,1000.);

  assert( nBins == static_cast<int>(fileNames.size()) );

  // ROOT output
  ROOT_OUT_FILE = new TFile(outNamePrefix+".root","RECREATE");
  EPS_OUTPUT = true;

  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");  
  std::vector<EtaPtBin*> etaPtBins;
  unsigned int nEtaBins = binAdm->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    // Read scale factors 
    std::vector<TH1*> hScaleFactorsStat;
    std::vector<TH1*> hScaleFactorsSyst;
    for(unsigned int f = 0; f < fileNames.size(); ++f) {
      
      hScaleFactorsStat.push_back(util::FileOps::readTH1(fileNames.at(f),"ScaleFactors_StatUncert_Eta"+util::toTString(etaBin),"hScaleFactorsStat"+util::toTString(f)));
      hScaleFactorsSyst.push_back(util::FileOps::readTH1(fileNames.at(f),"ScaleFactors_SystUncert_Eta"+util::toTString(etaBin),"hScaleFactorsSyst"+util::toTString(f)));
    }

    // Plot evolution
    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
      TH1* hEvoFrame = new TH1D("hEvolutionOfScaleFactorsFrame",";"+xAxisLabel+";"+SCALE,nBins,xBinEdges);
      if( binLabels[0] != "" ) {
	for(int bin = 1; bin <= hEvoFrame->GetNbinsX(); ++bin) {
	  hEvoFrame->GetXaxis()->SetBinLabel(bin,binLabels[bin-1]);
	}
	hEvoFrame->SetLabelSize(0.07,"X");
      }
      hEvoFrame->SetLineWidth(LINE_WIDTH);
      hEvoFrame->SetTitle(HEADER);
      hEvoFrame->SetLineStyle(2);
      hEvoFrame->GetYaxis()->SetRangeUser(0.01,3.3);
      TH1* hEvo = static_cast<TH1*>(hEvoFrame->Clone("hEvolutionOfScaleFactors"));
      hEvo->SetLineStyle(1);
      hEvo->SetMarkerStyle(20);
      hEvo->SetMarkerSize(MARKER_SIZE);
      TH1* hEvoSyst = static_cast<TH1*>(hEvo->Clone("hEvoSyst"));
      for(int bin = 1; bin <= hEvo->GetNbinsX(); ++bin) {
	hEvoFrame->SetBinContent(bin,1.);
	TH1* hStat = hScaleFactorsStat.at(bin-1);
	TH1* hSyst = hScaleFactorsSyst.at(bin-1);
	hEvo->SetBinContent(bin,hStat->GetBinContent(ptBin+1));
	hEvo->SetBinError(bin,hStat->GetBinError(ptBin+1));
	hEvoSyst->SetBinContent(bin,hSyst->GetBinContent(ptBin+1));
	hEvoSyst->SetBinError(bin,hSyst->GetBinError(ptBin+1));
      }
      TGraphAsymmErrors* gEvoSyst = util::HistOps::getUncertaintyBand(hEvoSyst,kYellow,1001);      

      // Label
      TPaveText* label = util::LabelFactory::createPaveText(2);
      if( xAxisLabel == LABEL_NVTX ) {
	label->AddText(LUMI+",  "+labWin+",  "+util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
      } else {
	label->AddText(LUMI+",  "+util::LabelFactory::etaCut(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
      }
      label->AddText(util::LabelFactory::ptAveCut(binAdm->ptMin(etaBin,ptBin),binAdm->ptMax(etaBin,ptBin)));

      TCanvas* can = util::HistOps::createTCanvas("can","Scale factor evolution (Eta "+util::toTString(etaBin)+", Pt "+util::toTString(ptBin),500,500);
      can->cd();
      hEvoFrame->Draw("HIST");
      gEvoSyst->Draw("E2same");
      hEvoFrame->Draw("HISTsame");
      hEvo->Draw("PE1same");
      label->Draw("same");
      gPad->RedrawAxis();
      toFiles(can,outNamePrefix+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin));

      delete hEvoFrame;
      delete hEvo;
      delete hEvoSyst;
      delete gEvoSyst;
      delete label;
      delete can;
    } // End of loop over pt bins

    for(unsigned int i = 0; i < hScaleFactorsStat.size(); ++i) {
      delete hScaleFactorsStat.at(i);
      delete hScaleFactorsSyst.at(i);
    }
  } // End of loop over eta bins

  ROOT_OUT_FILE->Close();
  delete ROOT_OUT_FILE;
}
