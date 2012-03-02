// $Id: getTailScalingFactors.C,v 1.14 2012/03/02 14:20:32 mschrode Exp $

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


const bool DEBUG = false;

const bool SHOW_HEADER = false;
const TString LUMI = util::StyleSettings::luminosity(4598);

const double BINLABEL_WIDTH = -0.52;
const double LEG_WIDTH = 0.48;
const TString PT3RELVAR = "p_{T,3} / p^{ave}_{T} threshold";
const double PT3PLOTMAX = 0.23;
const TString FASYM = "f_{asym}";
const TString FASYMMC = "f^{mc}_{asym}";
const TString FASYMDATA = "f^{data}_{asym}";
const TString FASYMGAUSS = "f^{gauss}_{asym}";
const TString FASYMTOY = "f^{toy}_{asym}";
const TString MCSMEAR = "Corrected MC";
const TString SCALE = "#rho_{tail}";
const TString ASYM_TAIL_START = "A_{tail}";
const TString ASYM_TAIL_START_EFF = "#hat{A}_{tail}";

const TString LUMI_LABEL = SHOW_HEADER ? "CMS preliminary, L = "+LUMI+",  #sqrt{s} = 7 TeV" : "#sqrt{s} = 7 TeV,  L = "+LUMI;
const TString HEADER = SHOW_HEADER ? LUMI_LABEL : "";

const int COLOR_GAUSS = 46;
const int COLOR_FILLED_ASYM = 38;
const int COLOR_FILLED_ASYM_SMEAR = 29;
const int COLOR_LINE_ASYM_SMEAR = 30;
const double MARKER_SIZE = SHOW_HEADER ? 1. : 1.4;
const int LINE_WIDTH = 2;



///////////////////////// TYPE DEFINITIONS /////////////////////////////////////////////////


const bool ROOT_OUTPUT = true;
bool EPS_OUTPUT = true;
TFile* ROOT_OUT_FILE = 0;

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

bool completePicName(int padIdx, const sampleTools::BinningAdmin* adm, TString &picName);
void writeLaTeXSlides(const TString &outNameId, const sampleTools::BinningAdmin* adm);
TString labelWindow(double nSigMin, double nSigMax);
double getTailStart(const TString &name);
void toFiles(TNamed* obj, const TString &name = "");
void toFiles(TCanvas* obj, const TString &name = "");



// ------------------------------------------------------------------------------------
class Pt3Bin {
public:
  Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Threshold, double nSigCore, double coreScaleFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel);
  Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, double windowRespMin, double windowRespMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel); 
  Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel);
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

  TH1* readHist(const TString &fileName, const TString &id) const;  
  TH1* readMCTruthResponse(const TString &fileName, const TString &type) const;
  //  void init(double nSigCore, double coreScalingFactor, double tailWindowDataMin, double tailWindowDataMax, double tailWindowMCMin, double tailWindowMCMax);
  void initBinLabel(const sampleTools::BinningAdmin* adm, const TString &jetLabel, bool isMCTruth = false);
  void getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr) const;
};



// ------------------------------------------------------------------------------------
class EtaPtBin {
public:
  EtaPtBin(unsigned int etaBinIdx, unsigned int ptBinIdx, double nSigCore, double coreScaleFactor,const sampleTools::BinningAdmin* adm, const TString &jetLabel);
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
  void addSymMCTruthResponse(const TString &fileName);
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
  const TString jetLabel_;

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
			   bool    variationCoreDown = false,
			   bool    variationExtrapolation = false,
			   bool    variationClosure = false,
			   bool    variationPUDown = false,
			   bool    variationPUUp = false            ) {

  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting up parameters" << std::endl;

  // User specified parameters
  const TString uid                    = "163337-180252";

  const double  minPt3Data   = 0.;
  const bool    fixDataShape = false;
  const double  nSigCore     = 2.;
  const TString jetAlgo      = "PF";
  const bool    archivePlots = true;


  // Sanity checks
  assert( nSigTailStart < nSigTailEnd );
  if( variationCoreUp )
    assert( !variationCoreDown && !variationExtrapolation && !variationClosure && !variationPUDown && !variationPUUp );
  else if( variationCoreDown )
    assert( !variationCoreUp && !variationExtrapolation && !variationClosure && !variationPUDown && !variationPUUp );
  else if( variationExtrapolation )
    assert( !variationCoreUp && !variationCoreDown && !variationClosure && !variationPUDown && !variationPUUp );
  else if( variationClosure )
    assert( !variationCoreUp && !variationCoreDown && !variationExtrapolation && !variationPUDown && !variationPUUp );
  else if( variationPUDown )
    assert( !variationCoreUp && !variationCoreDown && !variationExtrapolation && !variationClosure && !variationPUUp );
  else if( variationPUUp )
    assert( !variationCoreUp && !variationCoreDown && !variationExtrapolation && !variationPUDown && !variationClosure );


  if( variationCoreDown ) 
    std::cout << "\n***** Variation core down ***********************************\n" << std::endl;
  else if( variationCoreUp )
    std::cout << "\n***** Variation core up *************************************\n" << std::endl;
  else if( variationExtrapolation )
    std::cout << "\n***** Variation extrapolation *******************************\n" << std::endl;
  else if( variationClosure )
    std::cout << "\n***** Variation closure *************************************\n" << std::endl;
  else if( variationPUUp )
    std::cout << "\n***** Variation pile-up up **********************************\n" << std::endl;
  else if( variationPUDown )
    std::cout << "\n***** Variation pile-up down ********************************\n" << std::endl;


  if( SHOW_HEADER ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStyleNoteNoTitle();
    //util::StyleSettings::setStylePresentationNoTitle();
  }

  gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
  sampleTools::BinningAdmin* binAdm = new sampleTools::BinningAdmin("BinningAdmin.cfg");

  // Create output directories, file names, and labels
  TString tmpLabelWindow = util::toTString(nSigTailStart);
  if( tmpLabelWindow.Length() == 1 ) tmpLabelWindow += "0";
  tmpLabelWindow += "-";
  if( nSigTailEnd > 10. ) tmpLabelWindow += "Inf";
  else tmpLabelWindow += util::toTString(nSigTailEnd);
  tmpLabelWindow.ReplaceAll(".","");
  TString outLabel = "Tail_"+uid+"_Sig"+tmpLabelWindow+"_"+jetAlgo;
  if( variationCoreDown )
    outLabel += "_VarCoreDown";
  else if( variationCoreUp )
    outLabel += "_VarCoreUp";
  else if( variationExtrapolation )
    outLabel += "_VarExtrapolation";
  else if( variationClosure )
    outLabel += "_VarClosure";
  else if( variationPUUp )
    outLabel += "_VarPUUp";
  else if( variationPUDown )
    outLabel += "_VarPUDown";


  TString jetLabel = "Anti-k_{T} (R=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";




  // +++++ Bins +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating bins" << std::endl;

  std::vector<EtaPtBin*> etaPtBins;
  const unsigned int nEtaBins = binAdm->nEtaBins();
  unsigned int globalBin = 0;
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    double coreScaling = 0.;
    // Core scale factors from runs 163337-167151
    // Version 2011/07/20: closure uncertainty removed
    if( variationCoreDown ) {
      // Uncertainties down
      if( etaBin == 0 )      coreScaling = 0.000;
      else if( etaBin == 1 ) coreScaling = 0.001;
      else if( etaBin == 2 ) coreScaling = 0.032;
      else if( etaBin == 3 ) coreScaling = 0.043;
      else if( etaBin == 4 ) coreScaling = 0.090;
    } else if( variationCoreUp ) {
      // Uncertainties up
      if( etaBin == 0 )      coreScaling = 0.116;
      else if( etaBin == 1 ) coreScaling = 0.115;
      else if( etaBin == 2 ) coreScaling = 0.162;
      else if( etaBin == 3 ) coreScaling = 0.228;
      else if( etaBin == 4 ) coreScaling = 0.489;
    } else {
      // Nominal
      if( etaBin == 0 )      coreScaling = 0.052;
      else if( etaBin == 1 ) coreScaling = 0.057;
      else if( etaBin == 2 ) coreScaling = 0.096;
      else if( etaBin == 3 ) coreScaling = 0.134;
      else if( etaBin == 4 ) coreScaling = 0.288;
    }
    std::cout << "  Eta " << etaBin << ": Using core scale factor " << coreScaling << std::endl;

    for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
      if( DEBUG ) std::cout << "Setting up bin (" << etaBin << ", " << ptBin << ")" << std::endl;

      TString fileNameData = "~/results/ResolutionFit/Run2011_163337-180252_V10/ResTails_PtAveBins_Data2011_PF_L1FastJet_V10_REBINNED.root";
      TString fileNameMC = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";
      if( variationPUDown ) {
     	fileNameMC = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_PUDn_REBINNED.root";
      } else if( variationPUUp ) {
     	fileNameMC = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_PUUp_REBINNED.root";
      }
      TString fileNameMCTruth = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtGenAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";
      
//       TString fileNameData = "~/results/ResolutionFit/Run2011_163337-180252_V10/ResTails_PtAveBins_Data2011_PF_L1FastJet_V10_NVtxLess7_REBINNED.root";
//       TString fileNameMC = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtxLess7_REBINNED.root";
//        if( variationPUDown ) {
//     	fileNameMC = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtxLess7_REBINNED.root";
//           } else if( variationPUUp ) {
//     	fileNameMC = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_NVtxLess7_REBINNED.root";
//           }
//           TString fileNameMCTruth = "~/results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/ResTails_PtGenAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";

      // Create eta-pt bin
      EtaPtBin* bin = new EtaPtBin(etaBin,ptBin,nSigCore,coreScaling,binAdm,jetLabel);

      // Define window
      bin->findWindow(fileNameMC,0,nSigTailStart,nSigTailEnd);
      //bin->setWindow("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_Sig25-Inf_PF.root",nSigTailStart,nSigTailEnd);
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
    TH1* hAsymTailStart = new TH1D("hAsymTailStart_EtaBin"+util::toTString(etaBin),HEADER+";p^{ave}_{T} (GeV);Tail Start "+ASYM_TAIL_START,binAdm->nPtBins(etaBin),&(binAdm->ptBinEdges(etaBin).front()));
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
    TPaveText* label = 0;
    if( SHOW_HEADER ) {
      label = util::LabelFactory::createPaveText(3,BINLABEL_WIDTH);
    } else {
      label = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
      label->AddText(LUMI_LABEL);
    }
    label->AddText(util::LabelFactory::labelJet("ak5pf"));
    label->AddText(util::LabelFactory::labelEta(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
    label->AddText(labelWindow(nSigTailStart,nSigTailEnd));
    
    TLegend* leg = util::LabelFactory::createLegendCol(2,0.7*LEG_WIDTH);
    leg->AddEntry(hAsymTailStart,ASYM_TAIL_START+" ("+util::toTString(nSigTailStart)+" #sigma_{A})","L");
    leg->AddEntry(hAsymTailStartEff,ASYM_TAIL_START_EFF,"L");

    // Plot
    TCanvas* can = new TCanvas("can","",500,500);
    can->cd();
    util::HistOps::setYRangeLinear(hAsymTailStart,label->GetSize(),0.006);
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
      
      TH1* hFAsymData = static_cast<TH1*>(hFAsymMCSmeared->Clone("hFAsymData_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)));
      hFAsymData->SetMarkerStyle(20);
      hFAsymData->SetMarkerSize(MARKER_SIZE);
      
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

      // Ratios to Gauss
      TGraphAsymmErrors* gFAsymDataRelToGauss = util::HistOps::createRatioGraph(gFAsymData,hFAsymMCSmearedGauss);
      TGraphAsymmErrors* gFAsymMCSmearedRelToGauss = util::HistOps::createRatioGraph(gFAsymMCSmeared,hFAsymMCSmearedGauss);
      TH1* hFAsymMCSmearedRelToGauss = util::HistOps::createRatioPlot(hFAsymMCSmeared,hFAsymMCSmearedGauss);
      
      // Labels
      if( SHOW_HEADER ) {
	label = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
      } else {
	label = util::LabelFactory::createPaveText(5,BINLABEL_WIDTH);
	label->AddText(LUMI_LABEL);
      }
      label->AddText(util::LabelFactory::labelJet("ak5pf"));
      label->AddText(util::LabelFactory::labelEta(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
      label->AddText(util::LabelFactory::labelPt3(binAdm->ptSoftMax(pt3Bin)));
      label->AddText(labelWindow(nSigTailStart,nSigTailEnd));
      
      leg = util::LabelFactory::createLegendCol(3,LEG_WIDTH);
      leg->AddEntry(gFAsymData,"Data","P");
      leg->AddEntry(gFAsymMCSmeared,MCSMEAR,"LF");
      leg->AddEntry(hFAsymMCSmearedGauss,"Gauss","L");
      
      // Plot
      can = new TCanvas("can","",500,500);
      can->cd();
      util::HistOps::setYRange(hFAsymMCSmeared,label->GetSize(),0.7*hFAsymMCSmearedGauss->GetBinContent(hFAsymMCSmearedGauss->GetMinimumBin()));
      hFAsymMCSmeared->Draw("HIST");
      hFAsymMCSmearedGauss->Draw("HISTsame");
      gFAsymMCSmeared->Draw("E2same");
      hFAsymMCSmeared->Draw("HISTsame");
      gFAsymData->Draw("PE1same");
      leg->Draw("same");
      label->Draw("same");
      gPad->RedrawAxis();
      can->SetLogx();
      can->SetLogy();	
      toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsym");

      // Plot with ratio relative to Gauss
      delete can;
      can = util::HistOps::createRatioTopCanvas();
      TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
      TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hFAsymMCSmeared);
      TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hFAsymMCSmearedGauss,0.81,2.99);
      bRatioBottomFrame->SetLineStyle(hFAsymMCSmearedGauss->GetLineStyle());
      bRatioBottomFrame->SetLineColor(hFAsymMCSmearedGauss->GetLineColor());
      bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
      bRatioBottomFrame->GetXaxis()->SetNoExponent();
      can->cd();
      util::HistOps::setYRange(bRatioTopFrame,label->GetSize()+2,0.7*hFAsymMCSmearedGauss->GetBinContent(hFAsymMCSmearedGauss->GetMinimumBin()));
      bRatioTopFrame->Draw("HIST");
      hFAsymMCSmearedGauss->Draw("HISTsame");
      gFAsymMCSmeared->Draw("E2same");
      hFAsymMCSmeared->Draw("HISTsame");
      gFAsymData->Draw("PE1same");
      leg->Draw("same");
      label->Draw("same");
      gPad->RedrawAxis();
      can->SetLogx();
      can->SetLogy();	
      bRatioBottomPad->Draw();
      bRatioBottomPad->cd();
      bRatioBottomFrame->Draw("HIST");
      gFAsymMCSmearedRelToGauss->Draw("E2same");
      hFAsymMCSmearedRelToGauss->Draw("HISTsame");
      gFAsymDataRelToGauss->Draw("PE1same");
      bRatioBottomPad->SetLogx();
      bRatioBottomPad->RedrawAxis();
      toFiles(can,outLabel+"_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin)+"_FAsymRelToGaussBottom");

      delete hFAsymData;
      delete hFAsymMCSmeared;
      delete hFAsymMCSmearedGauss;
      delete gFAsymMCSmeared;
      delete gFAsymData;
      delete bRatioBottomPad;
      delete bRatioTopFrame;
      delete bRatioBottomFrame;
      delete gFAsymDataRelToGauss;
      delete gFAsymMCSmearedRelToGauss;
      delete hFAsymMCSmearedRelToGauss;
      delete label;
      delete leg;
      delete can;
    } // End of loop over pt3 bins
    
    // Plots
    can = new TCanvas("can","",500,500);
    
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

  ROOT_OUT_FILE->Close();
  delete ROOT_OUT_FILE;

  printWindowBorders(etaPtBins,binAdm);
  printMCClosure(etaPtBins,binAdm);
  printExtrapolation(etaPtBins,binAdm);

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
  nSigTailStart.push_back(2.5);
  nSigTailStart.push_back(3.0);
  nSigTailStart.push_back(3.5);

  double  nSigTailEnd = 1000.;
  for(std::vector<double>::const_iterator it = nSigTailStart.begin();
      it != nSigTailStart.end(); ++it) {
    EPS_OUTPUT = true ;
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
EtaPtBin::EtaPtBin(unsigned int etaBinIdx, unsigned int ptBinIdx, double nSigCore, double coreScaleFactor, const sampleTools::BinningAdmin* binAdm, const TString &jetLabel)
  : etaBin_(etaBinIdx), ptBin_(ptBinIdx),
    nSigCore_(nSigCore), coreScalingFactor_(coreScaleFactor),
    exMin_(0.), exMax_(0.19),
    binAdm_(binAdm), jetLabel_(jetLabel) {

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

  if( SHOW_HEADER ) {
    binLabel_ = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
  } else {
    binLabel_ = util::LabelFactory::createPaveText(5,BINLABEL_WIDTH);
    binLabel_->AddText(LUMI_LABEL);
  }
  binLabel_->AddText(jetLabel_);
  binLabel_->AddText(util::LabelFactory::labelEta(binAdm_->etaMin(etaBin_),binAdm_->etaMax(etaBin_)));
  binLabel_->AddText(util::LabelFactory::labelPtAve(binAdm_->ptMin(etaBin_,ptBin_),binAdm_->ptMax(etaBin_,ptBin_)));
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
  TLegend* leg = util::LabelFactory::createLegendCol(3,LEG_WIDTH);
  leg->AddEntry(gFTailMC_,"Asymmetry","P");
  leg->AddEntry(fExMC_,"Extrapolation","L");
  leg->AddEntry(gFTailToyMC_,FASYMTOY,"P");

  TCanvas* can = new TCanvas("can","",500,500);
  can->cd();
  TH1* hFrame = new TH1D("hFrame",HEADER+";"+PT3RELVAR,1000,0.,PT3PLOTMAX);
  hFrame->SetNdivisions(505);
  //hFrame->GetYaxis()->SetRangeUser(0.,2.3*gFTailMC_->GetY()[gFTailMC_->GetN()-1]);
  double minY = 0.2*(*std::min_element(gFTailMCGauss_->GetY(),gFTailMCGauss_->GetY()+gFTailMCGauss_->GetN()));
  double maxY = 7.;
  hFrame->GetYaxis()->SetRangeUser(minY,maxY);

  if( hasToyMC() ) {
    hFrame->GetYaxis()->SetTitle(FASYMMC);
    hFrame->Draw();
    fExMC_->Draw("same");
    if( hasSymMCTruth() ) gFTailMCTruth_->Draw("PE1same");
    gFTailToyMC_->Draw("PE1same");
    gFTailMC_->Draw("PE1same");
    binLabel_->DrawClone("same");
    leg->Draw("same");
    gPad->SetLogy(1);
    toFiles(can,outNameId+binId()+"_ExtrapolationMCClosure");
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
  gPad->SetLogy(1);
  toFiles(can,outNameId+binId()+"_Extrapolation");
  delete leg;
  gPad->SetLogy(0);


//   // Extrapolation MC Closure + Data + Shifted Extrapolation
//   if( hasToyMC() ) {
//     hFrame->GetYaxis()->SetTitle("( "+FASYMMC+" - Fit ) / Fit");
//     leg = util::LabelFactory::createLegendCol(5,LEG_WIDTH);
//     leg->AddEntry(gFTailData_,"Asymmetry Data","P");
//     leg->AddEntry(fExData_,"Extrapolation Data","L");
//     leg->AddEntry(gFTailMC_,"Asymmetry MC","P");
//     leg->AddEntry(fExMC_,"Extrapolation MC","L");
//     leg->AddEntry(gFTailToyMC_,"Toy MC","P");
    
//     can->cd();
//     hFrame->Draw();
//     fExMC_->Draw("same");
//     fExData_->Draw("same");
//     gFTailMC_->Draw("PE1same");
//     gFTailData_->Draw("PE1same");
//     if( hasSymMCTruth() ) {
//       gFTailMCTruth_->Draw("PE1same");
//       //gFTailMCTruthNonGauss_->Draw("PE1same");
//     }
//     gFTailToyMC_->Draw("PE1same");
//     binLabel_->DrawClone("same");
//     leg->Draw("same");
//     toFiles(can,outNameId+binId()+"_Extrapolation2");
//     delete leg;
//   }

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
  pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bin,thres,nSigCore_,coreScalingFactor_,tailMinBin_,tailMaxBin_,binAdm_,jetLabel_));
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
  double winMin = 1.+sqrt(2.)*tailWindowMinEff_;
  double winMax = 1.+sqrt(2.)*tailWindowMaxEff_;
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
void EtaPtBin::findWindow(const TString &fileName, unsigned int pt3Bin, double nSigTailWindowMin, double nSigTailWindowMax) {

  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft"+util::toTString(pt3Bin);

  TH1* h = util::FileOps::readTH1(fileName,histName);
  h->GetXaxis()->SetRangeUser(-1.,1);
  double width = 0.;
  double widthErr = 1000.;
  if( util::HistOps::fitCoreWidth(h,nSigCore_,width,widthErr) ) {
    if( width > 2.*h->GetRMS() ) {
      std::cerr << "\n ****** WARNING: Huge fitted width when defining window *********************" << std::endl;
      std::cerr << "  Eta " << etaBin() << ", Pt " << ptBin() << ": ";
      std::cerr << "Fitted Gaussian width = " << width << " +/- " << widthErr << " = " << width/h->GetRMS() << " RMS" << std::endl;
      width = h->GetRMS();
      std::cerr << "  Using RMS = " << h->GetRMS() << " instead" << std::endl;
      std::cerr << " ****************************************************************************\n" << std::endl;
    }
    tailWindowMin_ = nSigTailWindowMin*width;
    tailWindowMax_ = nSigTailWindowMax*width;
    tailMinBin_ = h->FindBin(std::abs(tailWindowMin_));
    tailMaxBin_ = h->FindBin(std::abs(tailWindowMax_));
    tailWindowMinEff_ = h->GetXaxis()->GetBinLowEdge(tailMinBin_);
    tailWindowMaxEff_ = h->GetXaxis()->GetBinUpEdge(tailMaxBin_);
    if( DEBUG ) std::cout << "EtaPtBin::findWindow(): sig(" << etaBin_ << ", " << ptBin_ << "): " << width << " \\pm " << widthErr << std::endl;
  } else {
    std::cerr << "ERROR in EtaPtBin::findWindow(): fit of core region did not converge" << std::endl;
    std::cerr << "  Window is not defined properly" << std::endl;
    tailWindowMinEff_ = 0.;
    tailWindowMaxEff_ = 0.;
    tailWindowMin_ = 0.;
    tailWindowMax_ = 0.;
    tailMinBin_ = 0;
    tailMaxBin_ = 0;
  }
  delete h;
  
  binLabel_->AddText(labelWindow(nSigTailWindowMin,nSigTailWindowMax));
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
  
  binLabel_->AddText(labelWindow(nSigTailWindowMin,nSigTailWindowMax));
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
  extraMCErr_ = fExMC_->GetParError(0)*fExMC_->Eval(0);


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
  extraDataErr_ = fExData_->GetParError(0)*fExData_->Eval(0);
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
Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Threshold, double nSigCore, double coreScaleFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel)
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
    util::HistOps::smearHistogram(hAsymMC_,hAsymMCSmeared_,sig_,coreScalingFactor_);
  } else {
    util::HistOps::smearHistogram(hAsymMC_,hAsymMCSmeared_,0.,0.);
  }

  // Width of smeared asymmetry
  if( !util::HistOps::fitCoreWidth(hAsymMCSmeared_,nSigCore,sigSmeared_,widthErr) ) sigSmeared_ = 0.;


//   hAsymData_->Scale(1./hAsymData_->Integral(tailStartBin_,tailEndBin_,"width"));
//   hAsymMC_->Scale(1./hAsymMC_->Integral(tailStartBin_,tailEndBin_,"width"));
//   hAsymMCSmeared_->Scale(1./hAsymMCSmeared_->Integral(tailStartBin_,tailEndBin_,"width"));


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
Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, double windowRespMin, double windowRespMax, const sampleTools::BinningAdmin* adm, const TString &jetLabel) 
  : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(1000), pt3Thres_(1000.),
    coreScalingFactor_(coreScaleFactor) {

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
    util::HistOps::smearHistogram(h,hSymResp_,width,coreScalingFactor_);
  } else {
    util::HistOps::smearHistogram(h,hSymResp_,0.,0.);
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
Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, const TString &jetLabel) 
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
    util::HistOps::smearHistogram(hAsymMC_,hAsymMCSmeared_,entries,width,coreScalingFactor_);
  } else {
    util::HistOps::smearHistogram(hAsymMC_,hAsymMCSmeared_,entries,0.,0.);
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
  if( SHOW_HEADER ) {
    if( isMCTruth ) binLabel_ = util::LabelFactory::createPaveText(3,BINLABEL_WIDTH);
    else binLabel_ = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
  } else {
    if( isMCTruth ) binLabel_ = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
    else binLabel_ = util::LabelFactory::createPaveText(5,BINLABEL_WIDTH);
    binLabel_->AddText(LUMI_LABEL);
  }
  binLabel_->AddText(jetLabel);
  binLabel_->AddText(util::LabelFactory::labelEta(adm->etaMin(etaBin_),adm->etaMax(etaBin_)));
  binLabel_->AddText(util::LabelFactory::labelPtAve(adm->ptMin(etaBin_,ptBin_),adm->ptMax(etaBin_,ptBin_)));
  //  if( !isMCTruth ) binLabel_->AddText(util::toTString(adm->ptSoftMin(pt3Bin_))+" < p_{T,3} / p^{ave}_{T} < "+util::toTString(adm->ptSoftMax(pt3Bin_)));
  if( !isMCTruth ) binLabel_->AddText(util::LabelFactory::labelPt3(adm->ptSoftMax(pt3Bin_)));
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
  TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft"+util::toTString(pt3Bin_);
  TH1* h = util::FileOps::readTH1(fileName,histName,id+histName);
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
  double nTail = h->Integral(start,end);
  fTail = nTail/nTotal;
  
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
  toFiles(can,outNameId+"PtAsymLog");

  // Linear scale
  util::HistOps::setYRange(hAsymMC_,5);
  hAsymMC_->GetXaxis()->SetRangeUser(-0.4,0.4);
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
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,"Asymmetry","",0.61,1.39);
  bRatioBottomFrame->SetFillColor(0);
  bRatioBottomFrame->SetLineWidth(LINE_WIDTH);
  can->cd();
  bRatioTopFrame->GetXaxis()->SetRangeUser(-0.24,0.24);
  bRatioTopFrame->GetXaxis()->SetNdivisions(505);
  util::HistOps::setYRange(bRatioTopFrame,5);
  bRatioTopFrame->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetRangeUser(-0.24,0.24);
  bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
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

  TLegend* leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hAsymData_,"Data","P");
  leg->AddEntry(hAsymMCSmeared_,MCSMEAR,"F");

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
  gauss->SetParameter(0,1./sqrt(2.*M_PI)/sigmaSmeared());
  gauss->SetParameter(1,0.);
  gauss->SetParameter(2,sigmaSmeared());
  gauss->SetLineWidth(2);
  gauss->SetLineColor(COLOR_GAUSS);
  gauss->SetFillStyle(3004);
  gauss->SetFillColor(gauss->GetLineColor());
  gauss->SetRange(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),gauss->GetX(3E-5,0.,1.));
  
  TLegend* legWin = util::LabelFactory::createLegendCol(4,LEG_WIDTH);
  legWin->AddEntry(hAsymData_,"Data","P");
  legWin->AddEntry(hAsymMCSmeared_,MCSMEAR,"F");
  legWin->AddEntry(gauss,"Asymmetry Gauss","F");
  legWin->AddEntry(arr,labelWindow(nSigTailStart,nSigTailEnd),"L");

  TPaveText* labelFAsym = util::LabelFactory::createPaveTextWithOffset(4,LEG_WIDTH,legWin->GetNRows());
  int nDigits = 4;
  labelFAsym->AddText(FASYMDATA+" = "+util::toTString(fTailData(),nDigits)+" #pm "+util::toTString(fTailDataErr(),nDigits));
  labelFAsym->AddText(FASYMMC+" = "+util::toTString(fTailMCSmeared(),nDigits)+" #pm "+util::toTString(fTailMCSmearedErr(),nDigits));
  labelFAsym->AddText(FASYMGAUSS+" = "+util::toTString(fTailMCSmearedGauss(),nDigits));

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
  hAsymMCSmeared_->GetXaxis()->SetRangeUser(-0.4,0.4);
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
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,"Asymmetry","",0.61,1.39);
  bRatioBottomFrame->SetFillColor(0);
  bRatioBottomFrame->SetLineWidth(LINE_WIDTH);
  can->cd();
  bRatioTopFrame->GetXaxis()->SetRangeUser(-0.24,0.24);
  bRatioTopFrame->GetXaxis()->SetNdivisions(505);
  util::HistOps::setYRange(bRatioTopFrame,5);
  bRatioTopFrame->Draw("HISTE");
  hAsymData_->Draw("PE1same");
  binLabel_->Draw("same");
  leg->Draw("same");
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetRangeUser(-0.24,0.24);
  bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
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
  TCanvas *can = new TCanvas(canName,canName,500,500);

  // Populated x-region
  int xMin = 1;
  int xMax = 1000;
  util::HistOps::findXRange(hPtAveSpecMC_,xMin,xMax);
  xMin = max(1,xMin-10);
  xMax = min(xMax+10,hPtAveSpecMC_->GetNbinsX());

  // Absolute data spectrum
  TLegend* leg = util::LabelFactory::createLegendCol(1,LEG_WIDTH);
  leg->AddEntry(hPtAveSpecData_,"Data N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
  hPtAveSpecData_->GetXaxis()->SetRange(xMin,xMax);
  util::HistOps::setYRange(hPtAveSpecData_,6,3E-1);
  can->cd();
  hPtAveSpecData_->Draw("PE1");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"PtAveSpectrumData");

  // Absolute MC spectrum
  delete leg;
  leg = util::LabelFactory::createLegendCol(1,LEG_WIDTH);
  leg->AddEntry(hPtAveSpecMC_,"MC N = "+util::toTString(hPtAveSpecMC_->GetEntries()),"F");
  util::HistOps::setYRange(hPtAveSpecMC_,6,3E-5);
  can->cd();
  hPtAveSpecMC_->GetXaxis()->SetRange(xMin,xMax);
  hPtAveSpecMC_->Draw("HIST");
  binLabel_->Draw("same");
  leg->Draw("same");
  can->SetLogy(1);
  toFiles(can,outNameId+"PtAveSpectrumMC");

  // Comparison (MC entries scaled to data)
  delete leg;
  leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
  leg->AddEntry(hPtAveSpecData_,"Data, N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
  leg->AddEntry(hPtAveSpecMC_,"MC, N = "+util::toTString(hPtAveSpecMC_->GetEntries()),"F");
  double scaleData = hPtAveSpecData_->Integral("width");
  double scaleMC = hPtAveSpecMC_->Integral("width");
  if( scaleData > 0. && scaleMC > 0. ) {
    hPtAveSpecMC_->Scale(scaleData/scaleMC);
    util::HistOps::setYRange(hPtAveSpecMC_,binLabel_->GetSize(),3E-1);
    can->cd();
    hPtAveSpecMC_->Draw("HISTE");
    hPtAveSpecData_->Draw("PE1same");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    toFiles(can,outNameId+"PtAveSpectra");
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
  leg->AddEntry(hSymResp_,MCSMEAR,"P");
  leg->AddEntry(win,"Tail","F");
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
  TCanvas *can = new TCanvas(canName,canName,500,500);
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
  std::cout << setprecision(1) << "    $" << adm->etaMin(etaBin) << " - " << adm->etaMax(etaBin) << "$ & $";
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
void printWindowBorders(const std::vector<EtaPtBin*> &bins, const sampleTools::BinningAdmin* adm) {
  bool exclRegion = bins.front()->tailWindowMax() < 1.;

  std::cout << "\n\n\n***  Window Borders  ***\n\n";

  std::cout << "\\begin{tabular}[ht]{cccc}\n\\toprule\n";
  std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $A_{0}";
  if( exclRegion ) std::cout << " - A_{1}";
  std::cout << "$ & $A_{0}/\\sigma";
  if( exclRegion ) std::cout << " - A_{1}/\\sigma";
  std::cout << "$ \\\\\n\\midrule" << std::endl;
  for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
    EtaPtBin* bin = *it;
    double min = bin->tailWindowMin();
    double max = bin->tailWindowMax();
    double sigSmear = bin->sigmaSmeared(0);

    printBin(bin->etaBin(),bin->ptBin(),adm);
    std::cout << std::setprecision(3) << " & $" << min;
    if( exclRegion ) {
      std::cout << " - " << max << "$ & $" << min/sigSmear << " - " << max/sigSmear;
    } else {
      std::cout << "$ & $" << min/sigSmear;
    }
    std::cout  << "$ \\\\" << std::endl;
    if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\midrule\n";
  }
  std::cout << "  \\end{tabular}\n\n";


  std::cout << "\n\n\n***  Tail vs Gaussian  ***\n\n";

  std::cout << "\\begin{tabular}[ht]{cccc}\n\\toprule\n";
  std::cout << "$|\\eta|$ & $\\ptave \\,(\\gevnospace)$ & $\\int^{";
  if( exclRegion ) std::cout << "A_{1}";
  else std::cout << "\\infty";
  std::cout << "}_{A_{0}}\\mathcal{G}$ & $\\fasymmc(\\pti{3}/\\ptave = 0.05)$ \\\\ \n\\midrule\n";
  for(EtaPtBinConstIt it = bins.begin(); it != bins.end(); ++it) {
    EtaPtBin* bin = *it;
    double min = bin->tailWindowMin();
    double max = bin->tailWindowMax();
    double sigSmear = bin->sigmaSmeared(0);
    double fGauss = (erf(max/sqrt(2.)/sigSmear) - erf(min/sqrt(2.)/sigSmear))/2.;
    double fAsym = bin->fTailMCSmeared(0);
    double fAsymErr = bin->fTailMCSmearedErr(0);

    printBin(bin->etaBin(),bin->ptBin(),adm);
    std::cout << std::setprecision(3) << " & $";
    std::cout << std::setprecision(4) << fGauss << "$ & $" << fAsym << " \\pm " << fAsymErr << "$ \\\\" << std::endl;
    if( bin->ptBin() == adm->nPtBins(bin->etaBin())-1 ) std::cout << "\\midrule\n";
  }
  std::cout << "  \\end{tabular}\n\n";

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
    label = "Window: "+util::toTString(nSigMin)+" - "+util::toTString(nSigMax);
  else
    label = "Tail: "+ASYM_TAIL_START+" = "+util::toTString(nSigMin);
  label += " #sigma_{A}";

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
//   g->SetFillStyle(0);
//   g->SetFillColor(0);
//   g->SetLineColor(color);

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
    histName = "hScale_Eta"+util::toTString(etaBin);
    TH1* hScaleNom = util::FileOps::readTH1(fileNamePrefix+".root",histName,histName+"_Nominal");

    // Graph of nominal scaling factors with ptAveMean as x
    TGraphAsymmErrors* gNomRatio = nomRatio(hScaleNom,hPtAveMeanData);


    // ****** Variations for systematic uncertainties *************************************

    // Read histograms with variations
    TH1* hScaleVarCoreUp = util::FileOps::readTH1(fileNamePrefix+"_VarCoreUp.root",histName,histName+"_VarCoreUp");
    TH1* hScaleVarCoreDown = util::FileOps::readTH1(fileNamePrefix+"_VarCoreDown.root",histName,histName+"_VarCoreDown");
    TH1* hScaleVarClosure = util::FileOps::readTH1(fileNamePrefix+"_VarClosure.root",histName,histName+"_VarClosure");
    TH1* hScaleVarPUUp = util::FileOps::readTH1(fileNamePrefix+"_VarPUUp.root",histName,histName+"_VarPUUp");
    TH1* hScaleVarPUDown = util::FileOps::readTH1(fileNamePrefix+"_VarPUDown.root",histName,histName+"_VarPUDown");
    TH1* hScaleVarExtra = util::FileOps::readTH1(fileNamePrefix+"_VarExtrapolation.root",histName,histName+"_VarExtrapolation");

    // Get relative uncertainties
    std::vector<TGraphAsymmErrors*> uncerts;
    uncerts.push_back(relUncertainty(hScaleNom,46,1.,hScaleVarCoreUp,hScaleVarCoreDown));
    uncerts.push_back(relUncertainty(hScaleNom,11,1.,hScaleVarPUUp,hScaleVarPUDown));
    uncerts.push_back(relUncertainty(hScaleNom,38,1.,hScaleVarExtra));
    uncerts.push_back(relUncertainty(hScaleNom,8,0.5,hScaleVarClosure));

    // Define labels
    std::vector<TString> uncertLabels;
    uncertLabels.push_back("Resolution");
    uncertLabels.push_back("Pile-up");
    uncertLabels.push_back("Extrapolation");
    uncertLabels.push_back("Non-closure");

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
    TPaveText* label = 0;
    if( SHOW_HEADER ) {
      label = util::LabelFactory::createPaveText(3,BINLABEL_WIDTH);
    } else {
      label = util::LabelFactory::createPaveText(4,BINLABEL_WIDTH);
      label->AddText(LUMI_LABEL);
    }
    label->AddText(util::LabelFactory::labelJet("AK5PF"));
    label->AddText(util::LabelFactory::labelEta(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
    label->AddText(labelWindow(getTailStart(fileNamePrefix),1000.));

    TLegend* legScale = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
    legScale->AddEntry(gNomRatio,"Stat. uncertainty","L");
    legScale->AddEntry(gUncertAbs,"Syst. uncertainty","F");

    TLegend* legUncert = util::LabelFactory::createLegendCol(uncerts.size()+1,LEG_WIDTH);
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      legUncert->AddEntry(uncerts.at(i),uncertLabels.at(i),"F");
    }
    legUncert->AddEntry(gUncertRelTotal,"Total","F");


    // Plot scaling factors and total uncertainty    
    hScaleFrame->SetTitle(HEADER);
    hScaleFrame->GetYaxis()->SetRangeUser(0.01,2.99);
    if( gNomRatio->GetY()[0] > 2. ) hScaleFrame->GetYaxis()->SetRangeUser(0.01,6.99);
    hScaleFrame->GetYaxis()->SetTitle(SCALE);
    hScaleFrame->GetXaxis()->SetMoreLogLabels();
    hScaleFrame->GetXaxis()->SetNoExponent();
    TCanvas* canScale = new TCanvas(outNamePrefix+"_canScale_"+util::toTString(etaBin),"Scale Factors Eta "+util::toTString(etaBin),500,500);
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

    // Plot relative uncertainties
    TText* line = label->GetLine(0);
    line->SetText(line->GetX(),line->GetY(),"#sqrt{s} = 7 TeV");
    TH1* hUncertsFrame = static_cast<TH1D*>(hScaleFrame->Clone("hUncertsFrame"+util::toTString(etaBin)));
    for(int bin = 1; bin <= hUncertsFrame->GetNbinsX(); ++bin) {
      hUncertsFrame->SetBinContent(bin,0.);
    }
    hUncertsFrame->GetYaxis()->SetRangeUser(-0.99,1.49);
    hUncertsFrame->GetXaxis()->SetMoreLogLabels();
    hUncertsFrame->GetXaxis()->SetNoExponent();
    hUncertsFrame->GetYaxis()->SetTitle("Relative Uncertainty");
    TCanvas* canRelUncerts = new TCanvas(outNamePrefix+"_canRelUncerts_"+util::toTString(etaBin),"Relative Uncertainties Eta "+util::toTString(etaBin),500,500);
    canRelUncerts->cd();
    hUncertsFrame->Draw("HIST");
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      uncerts.at(i)->Draw("E2same");
    }
    gUncertRelTotal->Draw("E2same");
    hUncertsFrame->Draw("HISTsame");
    label->Draw("same");
    legUncert->Draw("same");
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
      std::cout << "\\begin{tabular}{ccccc}\n\\hline\n";
      std::cout << "$|\\eta|$ & $\\ptave\\,(\\gevnospace)$ & $\\mean{\\ptave}\\,(\\gevnospace)$ & Scale Factor ($\\pm\\text{stat}{}^{+\\text{syst}}_{-\\text{syst}}$) & Scale Factor ($\\pm\\text{total}$) \\\\\n";
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
    if( etaBin == binAdm->nEtaBins()-1 ) std::cout << "\\end{tabular}\n";
  }

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
  fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_NVtxLess7_Sig25-Inf_PF.root");
  fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_NVtxGreater6_Sig25-Inf_PF.root");
  const TString outNamePrefix = "Tail_163337-180252_NVtxEvolution_Sig25_PF";

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
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for(unsigned int pt3Bin = 0; pt3Bin < binAdm->nPtSoftBins(); ++pt3Bin) {

      // Read histograms      
      std::vector<TH1*> hFAsymData;
      std::vector<TH1*> hFAsymMC;
      std::vector<TH1*> hFAsymGauss;
      for(unsigned int f = 0; f < fileNames.size(); ++f) {
	const TString binId = "_EtaBin"+util::toTString(etaBin)+"_Pt3Bin"+util::toTString(pt3Bin);
	hFAsymData.push_back(util::FileOps::readTH1(fileNames.at(f),"hFAsymData"+binId,"hFAsymData"+binId+"_"+util::toTString(f),false));
	hFAsymMC.push_back(util::FileOps::readTH1(fileNames.at(f),"hFAsymMCSmeared"+binId,"hFAsymMCSmeared"+binId+"_"+util::toTString(f),false));
	hFAsymGauss.push_back(util::FileOps::readTH1(fileNames.at(f),"hFAsymMCSmearedGauss"+binId,"hFAsymMCSmearedGauss"+binId+"_"+util::toTString(f),false));
      }

      // Plot evolution
      for(unsigned int ptBin = 0; ptBin < binAdm->nPtBins(etaBin); ++ptBin) {
	const TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin)+"_Pt3Bin"+util::toTString(pt3Bin);

	TH1* hEvoFAsymData = new TH1D("hEvolutionOfFAsymData"+binId,HEADER+";Number of Vertices",fileNames.size(),&(xBinEdges.front()));
	for(int bin = 1; bin <= hEvoFAsymData->GetNbinsX(); ++bin) {
	  hEvoFAsymData->GetXaxis()->SetBinLabel(bin,xBinLabels.at(bin-1));
	}
	hEvoFAsymData->SetLabelSize(0.07,"X");
	hEvoFAsymData->GetYaxis()->SetTitle(hFAsymData.front()->GetYaxis()->GetTitle());
	hEvoFAsymData->SetMarkerStyle(hFAsymData.front()->GetMarkerStyle());
	hEvoFAsymData->SetMarkerColor(hFAsymData.front()->GetMarkerColor());
	hEvoFAsymData->SetMarkerSize(hFAsymData.front()->GetMarkerSize());
	hEvoFAsymData->SetLineWidth(hFAsymData.front()->GetLineWidth());

	TH1* hEvoFAsymMC = static_cast<TH1*>(hEvoFAsymData->Clone("hEvolutionOfFAsymMCSmeared"+binId));
	hEvoFAsymMC->GetYaxis()->SetTitle(hFAsymMC.front()->GetYaxis()->GetTitle());
	hEvoFAsymMC->SetMarkerStyle(hFAsymMC.front()->GetMarkerStyle());
	hEvoFAsymMC->SetMarkerColor(hFAsymMC.front()->GetMarkerColor());
	hEvoFAsymMC->SetLineStyle(hFAsymMC.front()->GetLineStyle());
	hEvoFAsymMC->SetLineColor(hFAsymMC.front()->GetLineColor());
	hEvoFAsymMC->SetFillColor(hFAsymMC.front()->GetFillColor());
	hEvoFAsymMC->SetFillStyle(hFAsymMC.front()->GetFillStyle());

	TH1* hEvoFAsymGauss = static_cast<TH1*>(hEvoFAsymData->Clone("hEvolutionOfFAsymMCSmearedGauss"+binId));
	hEvoFAsymGauss->GetYaxis()->SetTitle(hFAsymGauss.front()->GetYaxis()->GetTitle());
	hEvoFAsymGauss->SetMarkerStyle(hFAsymGauss.front()->GetMarkerStyle());
	hEvoFAsymGauss->SetMarkerColor(hFAsymGauss.front()->GetMarkerColor());
	hEvoFAsymGauss->SetLineStyle(hFAsymGauss.front()->GetLineStyle());
	hEvoFAsymGauss->SetLineColor(hFAsymGauss.front()->GetLineColor());

	// Fill plots
	for(unsigned int f = 0; f < fileNames.size(); ++f) {
	  int binEvo = 1+f;
	  int binIn = 1+ptBin;
	  hEvoFAsymData->SetBinContent(binEvo,hFAsymData.at(f)->GetBinContent(binIn));
	  hEvoFAsymData->SetBinError(binEvo,hFAsymData.at(f)->GetBinError(binIn));
	  hEvoFAsymMC->SetBinContent(binEvo,hFAsymMC.at(f)->GetBinContent(binIn));
	  hEvoFAsymMC->SetBinError(binEvo,hFAsymMC.at(f)->GetBinError(binIn));
	  hEvoFAsymGauss->SetBinContent(binEvo,hFAsymGauss.at(f)->GetBinContent(binIn));
	  hEvoFAsymGauss->SetBinError(binEvo,hFAsymGauss.at(f)->GetBinError(binIn));
	}

	// Convert hEvoFAsymMC to graph to be able to plot error band
	TGraphAsymmErrors* gEvoFAsymMC = util::HistOps::getUncertaintyBand(hEvoFAsymMC,COLOR_FILLED_ASYM_SMEAR);	

	
	// Label
	TPaveText* label = 0;
	if( SHOW_HEADER ) {
	  label = util::LabelFactory::createPaveText(6,BINLABEL_WIDTH);
	} else {
	  label = util::LabelFactory::createPaveText(7,BINLABEL_WIDTH);
	  label->AddText(LUMI_LABEL);
	}
	label->AddText(util::LabelFactory::labelJet("AK5PF"));
	label->AddText(util::LabelFactory::labelEta(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin)));
	label->AddText(util::LabelFactory::labelPtAve(binAdm->ptMin(etaBin,ptBin),binAdm->ptMax(etaBin,ptBin)));
	label->AddText(util::LabelFactory::labelPt3(binAdm->ptSoftMax(pt3Bin)));
	label->AddText(labelWindow(getTailStart(fileNames.front()),1000.));
		       

	TLegend* leg = util::LabelFactory::createLegendCol(2,LEG_WIDTH);
	leg->AddEntry(hEvoFAsymData,"Data","P");
	leg->AddEntry(gEvoFAsymMC,MCSMEAR,"LF");
	//leg->AddEntry(hEvoFAsymGauss,"Gauss","L");
	
	util::HistOps::setYRange(hEvoFAsymMC,label->GetSize(),0.7*hEvoFAsymGauss->GetBinContent(hEvoFAsymGauss->GetMinimumBin()));
	TCanvas* can = new TCanvas("can","",500,500);
	can->cd();
	hEvoFAsymMC->Draw("HIST");
	//hEvoFAsymGauss->Draw("HISTsame");
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

      for(unsigned int i = 0; i < hFAsymData.size(); ++i) {
	delete hFAsymData.at(i);
	delete hFAsymMC.at(i);
	delete hFAsymGauss.at(i);
      }
    }// End of loop over pt3 bins
  } // End of loop over eta bins
  
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
//   fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_Sig25-Inf_PF_ScaleFactors.root");
//   fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_Sig30-Inf_PF_ScaleFactors.root");
//   fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_Sig35-Inf_PF_ScaleFactors.root");
//   const TString xAxisLabel = "Tail start "+ASYM_TAIL_START+" (#sigma_{A})";
//   const int nBins = 3;
//   const double xBinEdges[nBins+1] = { 2.25, 2.75, 3.25, 3.75 };
//  const TString binLabels[1] = { "" };
//   const TString outNamePrefix = "Tail_163337-180252_SigEvolution_PF_ScaleFactors";

  // Evolution with NVtx
  fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_NVtxLess7_Sig25-Inf_PF_ScaleFactors.root");
  fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_Sig25-Inf_PF_ScaleFactors.root");
  fileNames.push_back("ScaleFactors_163337-180252_2012-02-29/Tail_163337-180252_NVtxGreater6_Sig25-Inf_PF_ScaleFactors.root");
  const TString xAxisLabel = "Number of Vertices";
  const int nBins = 3;
  const double xBinEdges[nBins+1] = { 0., 1., 2., 3. };
  const TString binLabels[nBins+1] = { "#leq 6", "all", "#geq 7" };
  const TString outNamePrefix = "Tail_163337-180252_NVtxEvolution_Sig25_PF_ScaleFactors";

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
      hEvoFrame->GetYaxis()->SetRangeUser(0.01,2.99);
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

      // Label
      TPaveText* label = 0;
      if( SHOW_HEADER ) {
	label = util::LabelFactory::createPaveText(3,-0.9);
      } else {
	label = util::LabelFactory::createPaveText(4,-0.9);
	label->AddText(LUMI_LABEL);
      }
      label->AddText(util::LabelFactory::labelJet("AK5PF"));
      label->AddText(util::LabelFactory::labelEta(binAdm->etaMin(etaBin),binAdm->etaMax(etaBin))+",  "+util::LabelFactory::labelPtAve(binAdm->ptMin(etaBin,ptBin),binAdm->ptMax(etaBin,ptBin)));


      TCanvas* can = new TCanvas("can","Scale factor evolution (Eta "+util::toTString(etaBin)+", Pt "+util::toTString(ptBin),500,500);
      can->cd();
      hEvoFrame->Draw("HIST");
      hEvoSyst->Draw("PEsame");
      hEvo->Draw("PE1same");
      label->Draw("same");
      gPad->RedrawAxis();
      toFiles(can,outNamePrefix+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin));

      delete hEvoFrame;
      delete hEvo;
      delete hEvoSyst;
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
