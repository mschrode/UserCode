// $Id: $

#define UTILS_AS_HEADER_FILE

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/ConfigParser.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"
#include "../util/StyleSettings.h"
#include "../util/LabelFactory.h"



// ==============================================================================

// ------------------------------------------------------------------------------
class SystematicUncertainty {
public:
  SystematicUncertainty(unsigned int nPoints);
  SystematicUncertainty(const TString &label, const TString &fileName, int color = -1, int fillStyle = -1);
  SystematicUncertainty(const TString &label, const TString &fileNameNom, const TString &fileNameVar, int color = -1, int fillStyle = -1);
  SystematicUncertainty(const TString &label, const TString &fileNameNom, const TString &fileNameVarUp, const TString &fileNameVarDown, bool symmetrize, int color = -1, int fillStyle = -1);
  SystematicUncertainty(const TString &label, const TString &fileNameNom, double relUncert, int color = -1, int fillStyle = -1);
  SystematicUncertainty(const TString &label, const std::vector<SystematicUncertainty*> uncerts, int color = -1, int fillStyle = -1);
  SystematicUncertainty(const TString &label, double val, double ptMin, double ptMax, double up, double down, int color = -1, int fillStyle = -1);
  
  TString label() const { return label_; }

  TGraphAsymmErrors* gRelBandSteps(double mean) const;
  TGraphAsymmErrors* gRelBandSmooth(double mean) const;
  TGraphAsymmErrors* gBandSmooth(const TF1* f) const;
  TGraphAsymmErrors* gBandSmooth(const TGraphAsymmErrors* g) const;

  unsigned int nPoints() const { return ptSmooth_.size(); }
  double ptSmooth(unsigned int i) const { return ptSmooth_.at(i); }
  double ptSmoothErrUp(unsigned int i) const { return ptSmoothErrUp_.at(i); }
  double ptSmoothErrDown(unsigned int i) const { return ptSmoothErrDown_.at(i); }
  double ptSteps(unsigned int i) const { return ptSteps_.at(i); }
  double ptStepsErr(unsigned int i) const { return ptStepsErr_.at(i); }
  double relUp(unsigned int i) const { return relErrUp_.at(i); }
  double relDown(unsigned int i) const { return relErrDown_.at(i); }

  
 private:
  const TString label_;

  int color_;
  int fillStyle_;

  std::vector<double> ptSmooth_;
  std::vector<double> ptSmoothErrUp_;
  std::vector<double> ptSmoothErrDown_;
  std::vector<double> ptSteps_;
  std::vector<double> ptStepsErr_;
  std::vector<double> relVal_;
  std::vector<double> relErrUp_;
  std::vector<double> relErrDown_;
};


//  Empty uncertainty (everything initialized to 0)
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(unsigned int nPoints) : label_(""), color_(0), fillStyle_(0) {
  for(unsigned int i = 0; i < nPoints; ++i) {
    ptSmooth_.push_back(0.);
    ptSmoothErrUp_.push_back(0.);
    ptSmoothErrDown_.push_back(0.);
    ptSteps_.push_back(0.);
    ptStepsErr_.push_back(0.);
    relVal_.push_back(0.);
    relErrUp_.push_back(0.);
    relErrDown_.push_back(0.);
  }
}


//  Uncertainty from non-closure (50% of residual)
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(const TString &label, const TString &fileName, int color, int fillStyle) : label_(label) {
  color >= 0 ? color_ = color : color_ = kGray;
  fillStyle >= 0 ? fillStyle_ = fillStyle : fillStyle_ = 1001;  // Filled area

  // Get pt binning, measurement and truth from file
  TH1* hPtMean = util::FileOps::readTH1(fileName,"hPtMean");
  TH1* hRes = util::FileOps::readTH1(fileName,"hResolution");
  TF1* fMCTruth = util::FileOps::readTF1(fileName,"FittedResolution:trueRes");
  for(int bin = 1; bin <= hPtMean->GetNbinsX(); ++bin) {
    double x = hPtMean->GetBinContent(bin);
    ptSmooth_.push_back(x);
    ptSmoothErrDown_.push_back(x-hPtMean->GetXaxis()->GetBinLowEdge(bin));
    ptSmoothErrUp_.push_back(hPtMean->GetXaxis()->GetBinUpEdge(bin)-x);
    ptSteps_.push_back(hPtMean->GetBinCenter(bin));
    ptStepsErr_.push_back(0.5*hPtMean->GetBinWidth(bin));
    relVal_.push_back(0.);
    // Uncertainty is 50 % of difference
    double err = 0.5*std::abs(hRes->GetBinContent(bin)-fMCTruth->Eval(x))/hRes->GetBinContent(bin);
    relErrUp_.push_back(err);
    relErrDown_.push_back(err);
  }
  delete hPtMean;
  delete hRes;
  delete fMCTruth;
}


//  Uncertainty from one-sided variation (becomes symmetrized)
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(const TString &label, const TString &fileNameNom, const TString &fileNameVar, int color, int fillStyle) : label_(label) {

  color >= 0 ? color_ = color : color_ = kGray;
  fillStyle >= 0 ? fillStyle_ = fillStyle : fillStyle_ = 1001;  // Filled area

  // Get pt binning, nominal and varied measurement
  TH1* hPtMean = util::FileOps::readTH1(fileNameNom,"hPtMean");
  TH1* hResNom = util::FileOps::readTH1(fileNameNom,"hResolution","hNominal");
  TH1* hResVar = util::FileOps::readTH1(fileNameVar,"hResolution","hVar");
  for(int bin = 1; bin <= hPtMean->GetNbinsX(); ++bin) {
    double x = hPtMean->GetBinContent(bin);
    ptSmooth_.push_back(x);
    ptSmoothErrDown_.push_back(x-hPtMean->GetXaxis()->GetBinLowEdge(bin));
    ptSmoothErrUp_.push_back(hPtMean->GetXaxis()->GetBinUpEdge(bin)-x);
    ptSteps_.push_back(hPtMean->GetBinCenter(bin));
    ptStepsErr_.push_back(0.5*hPtMean->GetBinWidth(bin));
    relVal_.push_back(0.);
    // Uncertainty is difference
    double err = std::abs(hResNom->GetBinContent(bin)-hResVar->GetBinContent(bin))/hResNom->GetBinContent(bin);
    relErrUp_.push_back(err);
    relErrDown_.push_back(err);
  }
  delete hPtMean;
  delete hResNom;
  delete hResVar;
}  
  

//  Uncertainty from two-sided variation (optionally symmetrized)
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(const TString &label, const TString &fileNameNom, const TString &fileNameVarUp, const TString &fileNameVarDown, bool symmetrize, int color, int fillStyle) : label_(label) {

  color >= 0 ? color_ = color : color_ = kGray;
  fillStyle >= 0 ? fillStyle_ = fillStyle : fillStyle_ = 1001;  // Filled area

  // Get pt binning, nominal and varied measurement
  TH1* hPtMean = util::FileOps::readTH1(fileNameNom,"hPtMean");
  TH1* hResNom = util::FileOps::readTH1(fileNameNom,"hResolution","hNominal");
  TH1* hResVarUp = util::FileOps::readTH1(fileNameVarUp,"hResolution","hVarUp");
  TH1* hResVarDown = util::FileOps::readTH1(fileNameVarDown,"hResolution","hVarDown");
  for(int bin = 1; bin <= hPtMean->GetNbinsX(); ++bin) {
    double x = hPtMean->GetBinContent(bin);
    ptSmooth_.push_back(x);
    ptSmoothErrDown_.push_back(x-hPtMean->GetXaxis()->GetBinLowEdge(bin));
    ptSmoothErrUp_.push_back(hPtMean->GetXaxis()->GetBinUpEdge(bin)-x);
    ptSteps_.push_back(hPtMean->GetBinCenter(bin));
    ptStepsErr_.push_back(0.5*hPtMean->GetBinWidth(bin));
    relVal_.push_back(0.);
    // Uncertainty is difference
    double errUp = std::abs(hResNom->GetBinContent(bin)-hResVarUp->GetBinContent(bin))/hResNom->GetBinContent(bin);
    double errDown = std::abs(hResNom->GetBinContent(bin)-hResVarDown->GetBinContent(bin))/hResNom->GetBinContent(bin);
    if( symmetrize ) {
      errUp = 0.5*(errUp+errDown);
      errDown = errUp;
    }
    relErrUp_.push_back(errUp);
    relErrDown_.push_back(errDown);
  }
  delete hPtMean;
  delete hResNom;
  delete hResVarUp;
  delete hResVarDown;
}



//  Constant relative uncertainty
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(const TString &label, const TString &fileNameNom, double relUncert, int color, int fillStyle) : label_(label) {

  color >= 0 ? color_ = color : color_ = kGray;
  fillStyle >= 0 ? fillStyle_ = fillStyle : fillStyle_ = 1001;  // Filled area

  // Get pt binning, nominal and varied measurement
  TH1* hPtMean = util::FileOps::readTH1(fileNameNom,"hPtMean");
  for(int bin = 1; bin <= hPtMean->GetNbinsX(); ++bin) {
    double x = hPtMean->GetBinContent(bin);
    ptSmooth_.push_back(x);
    ptSmoothErrDown_.push_back(x-hPtMean->GetXaxis()->GetBinLowEdge(bin));
    ptSmoothErrUp_.push_back(hPtMean->GetXaxis()->GetBinUpEdge(bin)-x);
    ptSteps_.push_back(hPtMean->GetBinCenter(bin));
    ptStepsErr_.push_back(0.5*hPtMean->GetBinWidth(bin));
    relVal_.push_back(0.);
    relErrUp_.push_back(relUncert);
    relErrDown_.push_back(relUncert);
  }
  delete hPtMean;
}


//  Total uncertainty (quadratic sum)
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(const TString &label, const std::vector<SystematicUncertainty*> uncerts, int color, int fillStyle) : label_(label) {

  color >= 0 ? color_ = color : color_ = kGray;
  fillStyle >= 0 ? fillStyle_ = fillStyle : fillStyle_ = 1001;  // Filled area

  for(unsigned int i = 0; i < uncerts.at(0)->nPoints(); ++i) {
    ptSmooth_.push_back(uncerts.at(0)->ptSmooth(i));
    ptSmoothErrUp_.push_back(uncerts.at(0)->ptSmoothErrUp(i));
    ptSmoothErrDown_.push_back(uncerts.at(0)->ptSmoothErrDown(i));
    ptSteps_.push_back(uncerts.at(0)->ptSteps(i));
    ptStepsErr_.push_back(uncerts.at(0)->ptStepsErr(i));
    // Quadratic sum of uncertainties
    double relErrUp = 0;
    double relErrDown = 0;
    for(unsigned int k = 0; k < uncerts.size(); ++k) {
      relErrUp += pow(uncerts.at(k)->relUp(i),2.);
      relErrDown += pow(uncerts.at(k)->relDown(i),2.);
    }
    relVal_.push_back(0.);
    relErrUp_.push_back(sqrt(relErrUp));
    relErrDown_.push_back(sqrt(relErrDown));
  }
}


//  Constant relative uncertainty
// ------------------------------------------------------------------------------
SystematicUncertainty::SystematicUncertainty(const TString &label, double ptMin, double ptMax, double val, double up, double down, int color, int fillStyle) : label_(label) {

  color >= 0 ? color_ = color : color_ = kGray;
  fillStyle >= 0 ? fillStyle_ = fillStyle : fillStyle_ = 1001;  // Filled area

  for(int i = 0; i < 2; ++i) {
    if( i == 0 ) {
      ptSmooth_.push_back(ptMin);
      ptSmoothErrUp_.push_back(0.5*(ptMax-ptMin));
      ptSmoothErrDown_.push_back(0.);
    } else {
      ptSmooth_.push_back(ptMax);
      ptSmoothErrDown_.push_back(0.5*(ptMax-ptMin));
      ptSmoothErrUp_.push_back(0.);
    }
    ptSteps_.push_back(0.);
    ptStepsErr_.push_back(0.);
    relVal_.push_back(0.);
    relErrUp_.push_back(std::abs(up-val)/val);
    relErrDown_.push_back(std::abs(val-down)/val);
  }
}





//  Band of relative systematic uncertainty (pt steps)
// ------------------------------------------------------------------------------
TGraphAsymmErrors* SystematicUncertainty::gRelBandSteps(double mean) const {

  std::vector<double> y;
  for(unsigned int i = 0; i < relVal_.size(); ++i) {
    y.push_back(relVal_.at(i)+mean);
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(ptSteps_.size(),
					       &(ptSteps_.front()),&(y.front()),
					       &(ptStepsErr_.front()),&(ptStepsErr_.front()),
					       &(relErrDown_.front()),&(relErrUp_.front()));
  g->SetFillColor(fillStyle_ ? color_ : 0);
  g->SetFillStyle(fillStyle_);
  g->SetLineColor(color_);
  
  return g;  
}


//  Band of relative systematic uncertainty (smooth)
// ------------------------------------------------------------------------------
TGraphAsymmErrors* SystematicUncertainty::gRelBandSmooth(double mean) const {

  std::vector<double> y;
  for(unsigned int i = 0; i < relVal_.size(); ++i) {
    y.push_back(relVal_.at(i)+mean);
  }
  TGraphAsymmErrors* g = new TGraphAsymmErrors(ptSmooth_.size(),
					       &(ptSmooth_.front()),&(y.front()),
					       &(ptSmoothErrDown_.front()),&(ptSmoothErrUp_.front()),
					       &(relErrDown_.front()),&(relErrUp_.front()));
  g->SetFillColor(fillStyle_ ? color_ : 0);
  g->SetFillStyle(fillStyle_);
  g->SetLineColor(color_);
  
  return g;  
}


//  Band of systematic uncertainty around given function
// ------------------------------------------------------------------------------
TGraphAsymmErrors* SystematicUncertainty::gBandSmooth(const TF1* f) const {
  std::vector<double> y;
  std::vector<double> yeu;
  std::vector<double> yed;
  for(unsigned int i = 0; i < ptSmooth_.size(); ++i) {
    double val = f->Eval(ptSmooth_.at(i));
    y.push_back(val);
    yeu.push_back(val*relErrUp_.at(i));
    yed.push_back(val*relErrDown_.at(i));
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(ptSmooth_.size(),
					       &(ptSmooth_.front()),&(y.front()),
					       &(ptSmoothErrDown_.front()),&(ptSmoothErrUp_.front()),
					       &(yed.front()),&(yeu.front()));
  g->SetFillColor(fillStyle_ ? color_ : 0);
  g->SetFillStyle(fillStyle_);
  g->SetLineColor(color_);
  
  return g;  
}


//  Band of systematic uncertainty around TGraph
// ------------------------------------------------------------------------------
TGraphAsymmErrors* SystematicUncertainty::gBandSmooth(const TGraphAsymmErrors* g) const {
  assert( g->GetN() <= static_cast<int>(ptSmooth_.size()) );
  TGraphAsymmErrors* gBand = static_cast<TGraphAsymmErrors*>(g->Clone());
  for(int i = 0; i < gBand->GetN(); ++i) {
    double val = gBand->GetY()[i];
    gBand->SetPointError(i,ptSmoothErrDown_.at(i),ptSmoothErrUp_.at(i),
			 val*relErrDown_.at(i),val*relErrUp_.at(i));
  }
  gBand->SetFillColor(fillStyle_ ? color_ : 0);
  gBand->SetFillStyle(fillStyle_);
  gBand->SetLineColor(color_);
  
  return gBand;  
}







// ==============================================================================


TGraphAsymmErrors* getPhotonJetRatio(const TString &fileName);
TGraphAsymmErrors* getPhotonJetRatioSystematics(const TString &fileName, const TF1* fRatio);
void fitRatio(const TGraphAsymmErrors* gRatio, const TGraphAsymmErrors* gRatioBand, TF1* fit, double &ratio, double &stat, double &systDown, double &systUp, TGraphAsymmErrors* &gRatioStatBand);


// ------------------------------------------------------------------------------
void getResult(const TString &fileName, TGraphAsymmErrors* &gRes, TGraphAsymmErrors* &gExt, TF1* &fMCTruth, TF1* &fPLI, const TString &type) {

  TH1* hPtMean = util::FileOps::readTH1(fileName,"hPtMean");
  TH1* hRes = util::FileOps::readTH1(fileName,"hResolution");
  TH1* hExt = util::FileOps::readTH1(fileName,"hExtrapolationResult");

  std::vector<double> pt;
  std::vector<double> ptErr;
  std::vector<double> res;
  std::vector<double> resErr;
  std::vector<double> ext;
  std::vector<double> extErr;
  for(int bin = 1; bin <= hPtMean->GetNbinsX(); ++bin) {
    double val = hPtMean->GetBinContent(bin);
    pt.push_back((val == val ? val : 0.));
    val = hPtMean->GetBinError(bin);
    ptErr.push_back((val == val ? val : 1000.));
    val = hRes->GetBinContent(bin);
    res.push_back((val == val ? val : 0.));
    val = hRes->GetBinError(bin);
    resErr.push_back((val == val ? val : 1000.));
    val = hExt->GetBinContent(bin);
    ext.push_back((val == val ? val : 0.));
    val = hExt->GetBinError(bin);
    extErr.push_back((val == val ? val : 1000.));
  }
  gRes = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),&(ptErr.front()),&(ptErr.front()),&(resErr.front()),&(resErr.front()));
  gExt = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(ext.front()),&(ptErr.front()),&(ptErr.front()),&(extErr.front()),&(extErr.front()));

  if( type == "mc" ) {
    gRes->SetMarkerStyle(25);
    gRes->SetMarkerColor(kBlue);
    gRes->SetLineColor(kBlue);
    gExt->SetMarkerStyle(26);
    gExt->SetMarkerColor(kBlue);
    gExt->SetLineColor(kBlue);
  } else if( type == "data" ) {
    gRes->SetMarkerStyle(20);
    gRes->SetMarkerColor(kBlack);
    gRes->SetLineColor(kBlack);
    gExt->SetMarkerStyle(23);
    gExt->SetMarkerColor(kBlack);
    gExt->SetLineColor(kBlack);
  }

  fMCTruth = util::FileOps::readTF1(fileName,"FittedResolution:trueRes");
  TString name = fMCTruth->GetName();
  fMCTruth->SetName(name+"_"+type);
  fMCTruth->SetLineWidth(1);
  fMCTruth->SetLineColor(kRed);
  
  fPLI = util::FileOps::readTF1(fileName,"FittedResolution:ptGenAsym");
  name = fPLI->GetName();
  fPLI->SetName(name+"_"+type);
  fPLI->SetLineWidth(1);
  fPLI->SetLineStyle(2);

  delete hPtMean;
  delete hRes;
  delete hExt;
}


// ------------------------------------------------------------------------------
void getResult(const TString &fileName, TGraphAsymmErrors* &gRes, TGraphAsymmErrors* &gExt, const TString &type) {
  TF1* fMCTruth = 0;
  TF1* fPLI = 0;
  getResult(fileName,gRes,gExt,fMCTruth,fPLI,type);
  delete fMCTruth;
  delete fPLI;
}


// ------------------------------------------------------------------------------
TGraphAsymmErrors* getBiasCorrectedMeasurement(const TGraphAsymmErrors* gMeas, const TString &fileNameMC) {
  TH1* hRes = util::FileOps::readTH1(fileNameMC,"hResolution");
  TF1* fMCTruth = util::FileOps::readTF1(fileNameMC,"FittedResolution:trueRes");

  TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(gMeas->Clone());
  int bin = 1;
  for(int i = 0; i < g->GetN(); ++i, ++bin) {
    double x = 0.;
    double y = 0.;
    g->GetPoint(i,x,y);
    if( hRes->GetBinContent(bin) > 0. ) {
      double c = fMCTruth->Eval(x) / hRes->GetBinContent(bin);
      g->SetPoint(i,x,c*y);
    }
  }

  delete hRes;
  delete fMCTruth;

  return g;
}



// ------------------------------------------------------------------------------
void plotResults(unsigned int etaBin = 0, bool showCMSPreliminary = false) {

  //const int etaBin = 0;
  const TString jetAlgo = "PF";
  const TString dir = "results";
  const double xMin = 40.;
  const double xMax = 1100.;
  const double yMin = 1E-3;
  const double yMax = 0.38;
  const double lumi = 33;
  const bool showTitle = showCMSPreliminary;



  // +++++ Parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setting up parameters" << std::endl;
 
  TString nameSuffix = jetAlgo+"_Eta"+util::toTString(etaBin)+"_.root";

  TString nameData = dir+"/Data_"+nameSuffix;
  TString nameMC = dir+"/Closure_"+nameSuffix;
  TString nameExtrapolation = dir+"/VariationExtrapolation_"+nameSuffix;
  TString namePLIUp = dir+"/VariationPLIUp_"+nameSuffix;
  TString namePLIDown = dir+"/VariationPLIDown_"+nameSuffix;
  TString nameJESUp = dir+"/VariationJESUp_"+nameSuffix;
  TString nameJESDown = dir+"/VariationJESDown_"+nameSuffix;
  TString nameSpectrumUp = dir+"/VariationSpectrumUp_"+nameSuffix;
  TString nameSpectrumDown = dir+"/VariationSpectrumDown_"+nameSuffix;

  TString namePhotonJet = "";
  if( jetAlgo == "PF" ) {
    //namePhotonJet = dir+"/PhotonJetRatios_2010-12-17_PF_Eta"+util::toTString(etaBin)+".txt";
  }

  double etaMin = 0.;
  double etaMax = 0.;
  if( etaBin == 0 ) {
    etaMin = 0.;
    etaMax = 1.1;
  } else if( etaBin == 1 ) {
    etaMin = 1.1;
    etaMax = 1.7;
  } else if( etaBin == 2 ) {
    etaMin = 1.7;
    etaMax = 2.3;
  } else if( etaBin == 3 ) {
    etaMin = 2.3;
    etaMax = 5.0;
  }

  TString outNamePrefix = "ResFit_"+jetAlgo+"_Eta"+util::toTString(etaBin)+"_";

  TString jetLabel = "Anti-k_{T} (R=0.5) ";
  if( jetAlgo == "Calo" ) jetLabel += "Calo Jets";
  else if( jetAlgo == "PF" ) jetLabel += "PF Jets";
  else if( jetAlgo == "JPT" ) jetLabel += "JPT Jets";
  jetLabel += ", "+util::toTString(etaMin)+" < |#eta| < "+util::toTString(etaMax);

  TString title = "";
  if( showTitle ) {
    util::StyleSettings::paper();
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleX(0.7);
    gStyle->SetTitleH(0.038);
    title = "CMS preliminary";
  } else {
    util::StyleSettings::paperNoTitle();
  }


  


  // +++++ Closure and Data / MC comparison +++++++++++++++++++++++++++++++++++++

  std::cout << "Creating closure and data / MC comparison plots" << std::endl;

  // Get fitted resolutions
  TGraphAsymmErrors* gMC = 0;
  TGraphAsymmErrors* gMCExt = 0;
  TF1* fMCTruth = 0;
  TF1* fPLI = 0;
  getResult(nameMC,gMC,gMCExt,fMCTruth,fPLI,"mc"); 
  fMCTruth->SetRange(xMin,xMax);
  fPLI->SetRange(xMin,xMax);
  TGraphAsymmErrors* gMCRatio = util::HistOps::createRatioGraph(gMC,fMCTruth);

  TGraphAsymmErrors* gData = 0;
  TGraphAsymmErrors* gDataExt = 0;
  getResult(nameData,gData,gDataExt,"data");
  TGraphAsymmErrors* gDataRatio = util::HistOps::createRatioGraph(gData,fMCTruth);


  // Create labels
  
  double labelStart = 0.75;

  TPaveText* labelMC = util::LabelFactory::createPaveText(2);
  labelMC->AddText("CMS Simulation,  #sqrt{s} = 7 TeV, L = "+util::toTString(lumi)+" pb^{-1}");
  labelMC->AddText(jetLabel);

  TPaveText* labelData = util::LabelFactory::createPaveText(2);
  labelData->AddText("Data,  #sqrt{s} = 7 TeV, L = "+util::toTString(lumi)+" pb^{-1}");
  labelData->AddText(jetLabel);

  TLegend* legMC = util::LabelFactory::createLegendColWithOffset(4,labelStart,2);
  legMC->AddEntry(gMCExt,"Resolution (Extrapolated)","P");
  legMC->AddEntry(fPLI,"Particle Level Imbalance","L");
  legMC->AddEntry(gMC,"Resolution (Corrected)","P");
  legMC->AddEntry(fMCTruth,"Resolution (MC Truth)","L");

  TLegend* legData = util::LabelFactory::createLegendColWithOffset(4,labelStart,2);
  legData->AddEntry(gDataExt,"Resolution (Extrapolated)","P");
  legData->AddEntry(fPLI,"Particle Level Imbalance","L");
  legData->AddEntry(gData,"Resolution (Corrected)","P");
  legData->AddEntry(fMCTruth,"Resolution (MC Truth)","L");



  // Plot closure
  TH1 *tFrame = util::HistOps::createRatioTopFrame(xMin,xMax,yMin,yMax,"#sigma / p^{ref}_{T}");
  tFrame->SetTitle(title);
  tFrame->GetXaxis()->SetMoreLogLabels();
  tFrame->GetXaxis()->SetNoExponent();
  TH1 *bFrame = util::HistOps::createRatioBottomFrame(tFrame,"p^{ref}_{T}","GeV",0.81,1.29);
  bFrame->GetXaxis()->SetMoreLogLabels();
  bFrame->GetXaxis()->SetNoExponent();

  TCanvas *tCanMC = util::HistOps::createRatioTopCanvas();
  TPad *bPadMC = util::HistOps::createRatioBottomPad();
  tCanMC->cd();
  tFrame->Draw();
  fMCTruth->Draw("same");
  fPLI->Draw("same");
  gMCExt->Draw("PE1same");
  gMC->Draw("PE1same");
  labelMC->Draw("same");
  legMC->Draw("same");
  gPad->SetLogx();
  bPadMC->Draw();
  bPadMC->cd();
  bFrame->Draw();
  gMCRatio->Draw("PE1same");
  gPad->SetLogx();
  tCanMC->SaveAs(outNamePrefix+"ClosureMC.eps","eps");

  TCanvas *tCanData = util::HistOps::createRatioTopCanvas();
  TPad *bPadData = util::HistOps::createRatioBottomPad();
  tCanData->Clear();
  tCanData->cd();
  tFrame->Draw();
  fMCTruth->Draw("same");
  fPLI->Draw("same");
  gDataExt->Draw("PE1same");
  gData->Draw("PE1same");
  labelData->Draw("same");
  legData->Draw("same");
  gPad->SetLogx();
  bPadData->Draw();
  bPadData->cd();
  bFrame->Draw();
  gDataRatio->Draw("PE1same");
  gPad->SetLogx();
  //  tCanData->SaveAs(outNamePrefix+"ClosureData.eps","eps");


//   // print bias
//   std::cout << "\n\n";
//   std::cout << " BIAS\n";
//   for(int i = 0; i < gMC->GetN(); ++i) {
//     std::cout << "bias.push_back(" << gMC->GetY()[i] - fMCTruth->Eval(gMC->GetX()[i]) << ");\n";
//   }
//   std::cout << "\n\n";



  // +++++ Systematic uncertainties +++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Systematic uncertaintiy plots" << std::endl;

  TString labelUncertBias = "MC Closure";
  std::vector<SystematicUncertainty*> uncerts;
  uncerts.push_back(new SystematicUncertainty("JEC",nameMC,nameJESUp,nameJESDown,true,14));
  uncerts.push_back(new SystematicUncertainty(labelUncertBias,nameMC,46));
  uncerts.push_back(new SystematicUncertainty("Extrapolation",nameMC,nameExtrapolation,7));
  uncerts.push_back(new SystematicUncertainty("Spectrum",nameMC,nameSpectrumUp,nameSpectrumDown,true,38));
  uncerts.push_back(new SystematicUncertainty("Particle Level Imbalance",nameMC,namePLIUp,namePLIDown,false,8));
  
  SystematicUncertainty* total = new SystematicUncertainty("Total",uncerts,5);


  // Relative systematic uncertainties

  std::vector<TGraphAsymmErrors*> gUncerts;
  TLegend* legRelSyst = util::LabelFactory::createLegendColWithOffset(uncerts.size()+1,-0.7,2);
  gUncerts.push_back(total->gRelBandSteps(0.));
  for(unsigned int i = 0; i < uncerts.size(); ++i) {
    gUncerts.push_back(uncerts.at(i)->gRelBandSteps(0.));
    legRelSyst->AddEntry(gUncerts.back(),uncerts.at(i)->label(),"F");
  }
  legRelSyst->AddEntry(gUncerts.front(),"Total","F");

  TCanvas* canRelSyst = new TCanvas("canRelSyst","Relative Systematic Uncertainties",500,500);
  canRelSyst->cd();
  TH1* hFrameRelSyst = new TH1D("hFrameRelSyst","",1000,xMin,xMax);
  hFrameRelSyst->SetTitle(title);
  hFrameRelSyst->GetXaxis()->SetMoreLogLabels();
  hFrameRelSyst->GetXaxis()->SetNoExponent();
  for(int bin = 1; bin <= hFrameRelSyst->GetNbinsX(); ++bin) {
    hFrameRelSyst->SetBinContent(bin,0.);
  }
  hFrameRelSyst->SetLineStyle(2);
  hFrameRelSyst->GetYaxis()->SetRangeUser(-0.39,0.89);
  hFrameRelSyst->GetXaxis()->SetTitle(bFrame->GetXaxis()->GetTitle());
  hFrameRelSyst->GetYaxis()->SetTitle("Relative Uncertainty");
  hFrameRelSyst->GetYaxis()->SetMoreLogLabels();
  hFrameRelSyst->Draw();
  for(unsigned int i = 0; i < gUncerts.size(); ++i) {
    gUncerts.at(i)->Draw("E2same");
  }
  hFrameRelSyst->Draw("same");
  labelMC->Draw("same");
  legRelSyst->Draw("same");
  gPad->SetLogx();
  canRelSyst->SaveAs(outNamePrefix+"RelativeSystematicUncertainties.eps","eps");




  // +++++ Data / MC ratios +++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Creating data - MC ratio plots" << std::endl;


  // Systematic uncertainty w/o residual bias (cancels in ratio!)
  std::vector<SystematicUncertainty*> uncertsRatio;
  for(unsigned int i = 0; i < uncerts.size(); ++i) {
    if( uncerts.at(i)->label() != labelUncertBias ) {
      std::cout << "Adding uncert " << uncerts.at(i)->label() << std::endl;
      uncertsRatio.push_back(uncerts.at(i));
    }
  }
  SystematicUncertainty* totalRatio = new SystematicUncertainty("Total",uncertsRatio,5);

  TGraphAsymmErrors* gDataMCRatio = util::HistOps::createRatioGraph(gData,gMC);
  TGraphAsymmErrors* gDataMCRatioBand = totalRatio->gBandSmooth(gDataMCRatio);

  TF1* fDataMCRatio = new TF1("fDataMCRatio","pol0",xMin,xMax);
  fDataMCRatio->SetLineWidth(1);
  fDataMCRatio->SetLineColor(28);

  //  gDataMCRatio->Fit(fDataMCRatio,"0");

  // Ratio with combined systematic uncertainties
  double ratio = 0.;
  double ratioStat = 0.;
  double ratioSystDown = 0.;
  double ratioSystUp = 0.;
  TGraphAsymmErrors* gBandRatioStatUncert = 0;
  fitRatio(gDataMCRatio,gDataMCRatioBand,fDataMCRatio,ratio,ratioStat,ratioSystDown,ratioSystUp,gBandRatioStatUncert);

  // Syst band on ratio fit
  SystematicUncertainty* ratioSystUncert = new SystematicUncertainty("Syst. uncertainty",xMin,xMax,ratio,ratio+ratioSystUp,ratio-ratioSystDown,5);
  TGraphAsymmErrors* gBandRatioSystUncert = ratioSystUncert->gBandSmooth(fDataMCRatio);
  

  labelStart = -0.5;

  TPaveText* labelDataMC = util::LabelFactory::createPaveText(2);
  labelDataMC->AddText("#sqrt{s} = 7 TeV, L = "+util::toTString(lumi)+" pb^{-1}");
  labelDataMC->AddText(jetLabel);

  TLegend* legDataMC = util::LabelFactory::createLegendColWithOffset(2,labelStart,2);
  legDataMC->AddEntry(gData,"Resolution (Data)","P");
  legDataMC->AddEntry(gMC,"Resolution (MC)","P");

  TLegend* legDataMCRatio = util::LabelFactory::createLegendColWithOffset(4,labelStart,2);
  legDataMCRatio->AddEntry(gDataMCRatio,"Measurement","P");
  legDataMCRatio->AddEntry(fDataMCRatio,"Fit","L");
  legDataMCRatio->AddEntry(gBandRatioStatUncert,"Stat. Uncertainty","F");
  legDataMCRatio->AddEntry(gBandRatioSystUncert,"Syst. Uncertainty","F");

  TCanvas* tCanDataMC = util::HistOps::createRatioTopCanvas();
  TPad* bPadDataMC = util::HistOps::createRatioBottomPad();
  tCanDataMC->Clear();
  tCanDataMC->cd();
  tFrame->Draw();
  gMC->Draw("PE1same");
  gData->Draw("PE1same");
  labelDataMC->Draw("same");
  legDataMC->Draw("same");
  gPad->SetLogx();
  bPadDataMC->Draw();
  bPadDataMC->cd();
  bFrame->Draw();
  gDataMCRatio->Draw("PE1same");
  gPad->SetLogx();
  tCanDataMC->SaveAs(outNamePrefix+"DataMC.eps","eps");

//   std::cout << "*************************************\n";
//   std::cout << "MC\npt\tres\tstatErr\n";
//   for(int i = 0; i < gMC->GetN(); ++i) {
//     std::cout << gMC->GetX()[i] << " \t " << gMC->GetY()[i] << " +/- " << gMC->GetEYhigh()[i] << std::endl;
//   }
//   std::cout << "\nData\npt\tres\tstatErr\n";
//   for(int i = 0; i < gData->GetN(); ++i) {
//     std::cout << gData->GetX()[i] << " \t " << gData->GetY()[i] << " +/- " << gData->GetEYhigh()[i] << std::endl;
//   }

  TCanvas* canDataMCRatio = new TCanvas("canDataMCRatio","Data MC Ratio",500,500);
  canDataMCRatio->cd();
  TH1* hFrameDataMCRatio = util::HistOps::createRatioFrame(xMin,xMax,0.61,1.89,bFrame->GetXaxis()->GetTitle(),"#sigma(Data) / #sigma(MC)");
  hFrameDataMCRatio->SetTitle(title);
  hFrameDataMCRatio->GetXaxis()->SetMoreLogLabels();
  hFrameDataMCRatio->GetXaxis()->SetNoExponent();
//   hFrameDataMCRatio->Draw();
//   gDataMCRatioBand->Draw("E3same");
//   gDataMCRatio->Draw("PE1same");
//   fDataMCRatio->Draw("same");
//   labelDataMC->Draw("same");
//   legDataMCRatio->Draw("same");
//   gPad->SetLogx();
//   canDataMCRatio->SaveAs(outNamePrefix+"DataMCRatio.eps","eps");

  canDataMCRatio->cd();
  hFrameDataMCRatio->Draw();
  gBandRatioSystUncert->Draw("E3same");
  gBandRatioStatUncert->Draw("E2same");
  fDataMCRatio->Draw("same");
  gDataMCRatio->Draw("PE1same");
  labelDataMC->Draw("same");
  legDataMCRatio->Draw("same");
  gPad->SetLogx();
  gPad->RedrawAxis();
  canDataMCRatio->SaveAs(outNamePrefix+"DataMCRatio.eps","eps");



 

  // +++++ Scaled MC truth ++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Scaled MC truth plots" << std::endl;

  // Multiply formula of MCTruth by fitted ratio data / MC
  TString formulaMCTruth = fMCTruth->GetExpFormula();
  formulaMCTruth = "("+formulaMCTruth+")*["+util::toTString(fMCTruth->GetNpar())+"]";
  TF1* fScaledMCTruth = new TF1("fScaledMCTruth",formulaMCTruth,xMin,xMax);
  for(int i = 0; i < fMCTruth->GetNpar(); ++i) {
    fScaledMCTruth->SetParameter(i,fMCTruth->GetParameter(i));
    fScaledMCTruth->SetParError(i,fMCTruth->GetParError(i));
  }
  fScaledMCTruth->SetParameter(fMCTruth->GetNpar(),fDataMCRatio->GetParameter(0));
  fScaledMCTruth->SetParError(fMCTruth->GetNpar(),fDataMCRatio->GetParError(0));
  fScaledMCTruth->SetLineWidth(1);
  fScaledMCTruth->SetLineColor(kRed);
  fMCTruth->SetLineStyle(2);

  std::vector<SystematicUncertainty*> uncertsScaled;
  for(unsigned int i = 0; i < uncerts.size(); ++i) {
    uncertsScaled.push_back(uncerts.at(i));
  }
  uncertsScaled.push_back(new SystematicUncertainty("RatioStats",nameMC,fDataMCRatio->GetParError(0)));
  SystematicUncertainty* totalScaled = new SystematicUncertainty("Total",uncertsScaled,5);

  TGraphAsymmErrors* gScaledMCTruthBand = totalScaled->gBandSmooth(fScaledMCTruth);
  TGraphAsymmErrors* gScaledMCTruthRelBand = totalScaled->gRelBandSmooth(1.);
  TGraphAsymmErrors* gScaledMCTruthDataRatio = util::HistOps::createRatioGraph(gData,fScaledMCTruth);


  TLegend* legScaledMCTruth = util::LabelFactory::createLegendColWithOffset(4,0.7,2);
  legScaledMCTruth->AddEntry(gData,"Data","P");
  legScaledMCTruth->AddEntry(fMCTruth,"MC Truth","L");
  legScaledMCTruth->AddEntry(fScaledMCTruth,"Scaled MC Truth","L");
  legScaledMCTruth->AddEntry(gScaledMCTruthBand,"Systematic Uncertainty","F");

  TH1 *bFrameSym = util::HistOps::createRatioBottomFrame(tFrame,"p^{ref}_{T}","GeV",0.71,1.29);
  bFrameSym->GetXaxis()->SetMoreLogLabels();
  bFrameSym->GetXaxis()->SetNoExponent();

  TCanvas *tCanScaledMCTruth = util::HistOps::createRatioTopCanvas();
  TPad *bPadScaledMCTruth = util::HistOps::createRatioBottomPad();
  tCanScaledMCTruth->Clear();
  tCanScaledMCTruth->cd();
  tFrame->Draw();
  gScaledMCTruthBand->Draw("E3same");
  //fMCTruth->Draw("same");
  fScaledMCTruth->Draw("same");
  gData->Draw("PE1same");
  labelDataMC->Draw("same");
  legScaledMCTruth->Draw("same");
  gPad->SetLogx();
  bPadScaledMCTruth->Draw();
  bPadScaledMCTruth->cd();
  bFrameSym->Draw();
  gScaledMCTruthRelBand->Draw("E3same");
  gScaledMCTruthDataRatio->Draw("PE1same");
  bFrameSym->Draw("same");
  gPad->SetLogx();
  //tCanScaledMCTruth->SaveAs(outNamePrefix+"ClosureScaledMCTruth.eps","eps");




  // +++++ Scaled MC truth and bias corrected data ++++++++++++++++++++++++++++++

  std::cout << "Scaled MC truth plots (bias corrected measurement)" << std::endl;

  // Systematic uncertainty w/o residual bias
  std::vector<SystematicUncertainty*> uncertsBiasCorr;
//   for(unsigned int i = 0; i < uncerts.size(); ++i) {
//     if( uncerts.at(i)->label() != labelUncertBias ) {
//       std::cout << "Adding uncert " << uncerts.at(i)->label() << std::endl;
//       uncertsBiasCorr.push_back(uncerts.at(i));
//     }
//   }
  uncertsBiasCorr.push_back(new SystematicUncertainty("RatioStats",nameMC,0.5*(ratioSystDown+ratioSystUp)));
  uncertsBiasCorr.push_back(new SystematicUncertainty("RatioStats",nameMC,ratioStat));

  SystematicUncertainty* totalBiasCorr = new SystematicUncertainty("Total",uncertsBiasCorr,5);

  TGraphAsymmErrors* gDataBiasCorr = getBiasCorrectedMeasurement(gData,nameMC);
  TGraphAsymmErrors* gScaledMCTruthBandBiasCorr = totalBiasCorr->gBandSmooth(fScaledMCTruth);
  TGraphAsymmErrors* gScaledMCTruthRelBandBiasCorr = totalBiasCorr->gRelBandSmooth(1.);
  TGraphAsymmErrors* gScaledMCTruthDataBiasCorrRatio = util::HistOps::createRatioGraph(gDataBiasCorr,fScaledMCTruth);


  TLegend* legScaledMCTruthBiasCorr = util::LabelFactory::createLegendColWithOffset(3,0.7,2);
  legScaledMCTruthBiasCorr->AddEntry(gDataBiasCorr,"Bias Corrected Data","P");
  legScaledMCTruthBiasCorr->AddEntry(fScaledMCTruth,"Scaled MC Truth","L");
  legScaledMCTruthBiasCorr->AddEntry(gScaledMCTruthBandBiasCorr,"Systematic Uncertainty","F");

  TCanvas *tCanScaledMCTruthBiasCorr = util::HistOps::createRatioTopCanvas();
  TPad *bPadScaledMCTruthBiasCorr = util::HistOps::createRatioBottomPad();
  tCanScaledMCTruthBiasCorr->Clear();
  tCanScaledMCTruthBiasCorr->cd();
  tFrame->Draw();
  gScaledMCTruthBandBiasCorr->Draw("E3same");
  fScaledMCTruth->Draw("same");
  gDataBiasCorr->Draw("PE1same");
  labelDataMC->Draw("same");
  legScaledMCTruthBiasCorr->Draw("same");
  gPad->SetLogx();
  bPadScaledMCTruthBiasCorr->Draw();
  bPadScaledMCTruthBiasCorr->cd();
  bFrameSym->Draw();
  gScaledMCTruthRelBandBiasCorr->Draw("E3same");
  gScaledMCTruthDataBiasCorrRatio->Draw("PE1same");
  bFrameSym->Draw("same");
  gPad->SetLogx();
  tCanScaledMCTruthBiasCorr->SaveAs(outNamePrefix+"ClosureScaledMCTruthBiasCorr.eps","eps");



  // +++++ Ratio output +++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "\n\nFITTED RATIO\n\n";
  std::cout << std::setprecision(2) << "$" << etaMin << " - " << etaMax << "$ & ";
  std::cout << std::setprecision(4) << " $" << ratio << std::setprecision(3) << " \\pm " << ratioStat << "^{+" << ratioSystUp << "}_{-" << ratioSystDown << "}$ \\\\\n\n\n";

  SystematicUncertainty* uncertRatioOutputJES = new SystematicUncertainty(uncerts.at(0)->nPoints());
  std::vector<SystematicUncertainty*> uncertsRatioOutputOther;
  for(unsigned int i = 0; i < uncerts.size(); ++i) {
    if( uncerts.at(i)->label() != labelUncertBias ) {
      if( uncerts.at(i)->label() == "JEC" ) {
	delete uncertRatioOutputJES;
	uncertRatioOutputJES = uncerts.at(i);
      } else {
	std::cout << "Adding uncert to ratio output: " << uncerts.at(i)->label() << std::endl;
	uncertsRatioOutputOther.push_back(uncerts.at(i));
      }
    }
  }
  SystematicUncertainty* totalRatioOutput = new SystematicUncertainty("Total",uncertsRatioOutputOther,5);
  std::vector<unsigned int> ptBins;
  for(int i = 0; i < gDataMCRatio->GetN(); ++i) {
    if( gDataMCRatio->GetY()[i] == 0. ) {
      std::cerr << "WARNING: problematic ptBin " << i << std::endl;
    } else {
      ptBins.push_back(i);
    }
  }
  std::cout << std::setprecision(5) << std::flush;
  std::cout << "\n\n#Eta " << etaMin << " -- " << etaMax << std::endl;
  std::cout << "MeanPt:" << std::flush;
  for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
    std::cout << "  " << gDataMCRatio->GetX()[*it] << std::flush;
  }
  std::cout << "\nMeanPtError:" << std::flush;
  for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
    std::cout << "  " << gDataMCRatio->GetEXhigh()[*it] << std::flush;
  }
  std::cout << "\nRatio:" << std::flush;
  for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
    std::cout << "  " << gDataMCRatio->GetY()[*it] << std::flush;
  }
  std::cout << "\nRatioError:" << std::flush;
   for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
     std::cout << "  " << gDataMCRatio->GetEYhigh()[*it] << std::flush;
   }
   std::cout << "\nSystJESUp:" << std::flush;
   for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
     std::cout << "  " << gDataMCRatio->GetY()[*it]*uncertRatioOutputJES->relUp(*it) << std::flush;
   }
   std::cout << "\nSystJESDown:" << std::flush;
   for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
     std::cout << "  " << gDataMCRatio->GetY()[*it]*uncertRatioOutputJES->relDown(*it) << std::flush;
   }
   std::cout << "\nSystOtherUp:" << std::flush;
   for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
     std::cout << "  " << gDataMCRatio->GetY()[*it]*totalRatioOutput->relUp(*it) << std::flush;
   }
   std::cout << "\nSystOtherDown:" << std::flush;
   for(std::vector<unsigned int>::const_iterator it = ptBins.begin(); it != ptBins.end(); ++it) {
     std::cout << "  " << gDataMCRatio->GetY()[*it]*totalRatioOutput->relDown(*it) << std::flush;
   }
   std::cout << std::endl;



  // +++++ Data / MC ratio combined Photon-Jet result +++++++++++++++++++++++++++

  if( namePhotonJet != "" ) {
    std::cout << "Creating combined dijet and photon+jet ratio plots" << std::endl;

    double combinedXMin = 20.;

    TGraphAsymmErrors* gPhotonJetRatio = getPhotonJetRatio(namePhotonJet);
    gPhotonJetRatio->SetMarkerStyle(21);
    gPhotonJetRatio->SetMarkerColor(kRed);
    gPhotonJetRatio->SetLineColor(gPhotonJetRatio->GetMarkerColor());
    TGraphAsymmErrors* gCombinedRatio = util::HistOps::combineTGraphs(gPhotonJetRatio,gDataMCRatio);
    gCombinedRatio->SetMarkerStyle(20);
    gCombinedRatio->SetMarkerColor(kBlue);
    gCombinedRatio->SetLineColor(gCombinedRatio->GetMarkerColor());
  
    // Fit ratio to combined data points
    TF1* fCombinedRatio = new TF1("fCombinedRatio","pol0",combinedXMin,xMax);
    fCombinedRatio->SetLineWidth(1);
    fCombinedRatio->SetLineColor(28);
    gCombinedRatio->Fit(fCombinedRatio,"0");

    // Systematic uncertainties in different regions
    TGraphAsymmErrors* gPhotonJetRatioSyst = getPhotonJetRatioSystematics(namePhotonJet,fCombinedRatio);
    gPhotonJetRatioSyst->SetFillColor(46);//gPhotonJetRatio->GetMarkerColor());
    gPhotonJetRatioSyst->SetLineColor(gPhotonJetRatioSyst->GetFillColor());
    //gPhotonJetRatioSyst->SetFillStyle(3005);

    std::vector<SystematicUncertainty*> uncertsCombinedDijet;
    for(unsigned int i = 0; i < uncerts.size(); ++i) {
      if( uncerts.at(i)->label() != labelUncertBias ) {
	std::cout << "Adding uncert " << uncerts.at(i)->label() << std::endl;
	uncertsCombinedDijet.push_back(uncerts.at(i));
      }
    }
    SystematicUncertainty* totalCombinedDijet = new SystematicUncertainty("Total",uncertsCombinedDijet,38);

    TGraphAsymmErrors* gDijetRatioSyst = totalCombinedDijet->gBandSmooth(fCombinedRatio);
    gDijetRatioSyst->SetFillColor(38);//gCombinedRatio->GetMarkerColor());
    gDijetRatioSyst->SetLineColor(gDijetRatioSyst->GetFillColor());
    //gDijetRatioSyst->SetFillStyle(3004);

    TLegend* legCombinedRatio = util::LabelFactory::createLegendColWithOffset(3,-0.6,2);
    legCombinedRatio->AddEntry(gPhotonJetRatio,"#gamma + Jet Measurement","P");
    legCombinedRatio->AddEntry(gCombinedRatio,"Dijet Measurement","P");
    legCombinedRatio->AddEntry(fCombinedRatio,"Fit","L");
  
    TCanvas* canCombinedRatio = new TCanvas("canCombinedRatio","Combined Result",500,500);
    canCombinedRatio->cd();
    TH1* hFrameCombinedRatio = util::HistOps::createRatioFrame(combinedXMin,xMax,0.61,1.89,"p_{T} (GeV)","#sigma(Data) / #sigma(MC)");
    hFrameCombinedRatio->SetTitle(title);
    hFrameCombinedRatio->GetYaxis()->SetMoreLogLabels();
    hFrameCombinedRatio->Draw();
    gDijetRatioSyst->Draw("E3same");
    gPhotonJetRatioSyst->Draw("E3same");
    gCombinedRatio->Draw("PE1same");
    gPhotonJetRatio->Draw("PE1same");
    fCombinedRatio->Draw("same");
    labelDataMC->Draw("same");
    legCombinedRatio->Draw("same");
    gPad->SetLogx();
    canCombinedRatio->SaveAs(outNamePrefix+"CombinedPotonJetDijetRatio.eps","eps");


    // Cite ratio with total uncertainties

    // For now, average (unweighted) dijet uncertainty
    // for all points beyon photon jet
    double statErr = fCombinedRatio->GetParError(0);
    double minPtPhotonJet = 20.;
    double maxPtPhotonJet = gPhotonJetRatio->GetX()[gPhotonJetRatio->GetN()-1];
    double maxPtDijets = gUncerts.front()->GetX()[gUncerts.front()->GetN()-1]+gUncerts.front()->GetEXhigh()[gUncerts.front()->GetN()-1];
    double photonJetErrLow = sqrt( pow(statErr,2.) + pow(gPhotonJetRatioSyst->GetEYlow()[0],2.) );
    double photonJetErrUp = sqrt( pow(statErr,2.) + pow(gPhotonJetRatioSyst->GetEYhigh()[0],2.) );
    double dijetErrUp = 0;
    double dijetErrLow = 0;
    int nPoints = 0;
    for(int i = 0; i < gDijetRatioSyst->GetN(); ++i) {
      if( gDijetRatioSyst->GetX()[i] > maxPtPhotonJet ) {
	nPoints++;
	dijetErrUp += gDijetRatioSyst->GetEYhigh()[i];
	dijetErrLow += gDijetRatioSyst->GetEYlow()[i];
      }
    }
    dijetErrUp /= nPoints;
    dijetErrLow /= nPoints;
    dijetErrUp = sqrt( pow(dijetErrUp,2.) + pow(statErr,2.) );
    dijetErrLow = sqrt( pow(dijetErrLow,2.) + pow(statErr,2.) );

    std::cout << "\n\nCombined measurement (Eta " << etaBin << ")\n";
    std::cout << "  Ratio data / MC                            \t:  " << fCombinedRatio->GetParameter(0) << std::endl;
    std::cout << "  Statistical uncertainty                    \t:  " << statErr << std::endl;
    std::cout << "  Systematic uncertainty (" << minPtPhotonJet << " - " << maxPtPhotonJet << ") \t:  +" << photonJetErrUp << " -" << photonJetErrLow << std::endl;
    std::cout << "  Systematic uncertainty (" << maxPtPhotonJet << " - " << maxPtDijets << ") \t:  +" << dijetErrUp << " -" << dijetErrLow << std::endl;
  } // End of combined photon-jet and dijet section
}



TGraphAsymmErrors* getPhotonJetRatio(const TString &fileName) {

  // Read values from file
  util::ConfigParser parser(fileName.Data());
  std::vector<double> pt = parser.readDoubleVec("MeanPt",":");
  std::vector<double> ptErr = parser.readDoubleVec("MeanError",":");
  std::vector<double> ratio = parser.readDoubleVec("Ratio",":");
  std::vector<double> ratioErr = parser.readDoubleVec("RatioError",":");

  // Create TGraph
  TGraphAsymmErrors* g = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(ratio.front()),
					       &(ptErr.front()),&(ptErr.front()),
					       &(ratioErr.front()),&(ratioErr.front()));

  return g;
}



TGraphAsymmErrors* getPhotonJetRatioSystematics(const TString &fileName, const TF1* fRatio) {

  // Read pt from file
  util::ConfigParser parser(fileName.Data());
  std::vector<double> pt = parser.readDoubleVec("MeanPt",":");
  std::vector<double> ptErr = parser.readDoubleVec("MeanError",":");

  // Store ratio from TF1
  std::vector<double> ratio;
  for(unsigned int i = 0; i < pt.size(); ++i) {
    ratio.push_back(fRatio->Eval(pt.at(i)));
  }  

  // Systematic uncertainties
  std::vector<double> systUp = parser.readDoubleVec("SystUp",":");
  std::vector<double> systDown = parser.readDoubleVec("SystDown",":");
  // In case of constant systematic uncertainty
  // fill up vectors
  if( systUp.size() == 1 && systDown.size() == 1 ) {
    for(unsigned int i = 1; i < pt.size(); ++i) {
      systUp.push_back(systUp.front());
      systDown.push_back(systDown.front());
    }
  }
  // Compute absolute systematic uncertainty
  for(unsigned int i = 0; i < pt.size(); ++i) {
    systUp.at(i) = ratio.at(i)*systUp.at(i);
    systDown.at(i) = ratio.at(i)*systDown.at(i);  
  }

  // Create TGraph
  TGraphAsymmErrors* g = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(ratio.front()),
					       &(ptErr.front()),&(ptErr.front()),
					       &(systDown.front()),&(systUp.front()));

  return g;
}



// Fit ratio and determine statistical and systematic uncertainty
// on ratio
// -------------------------------------------------------------------------------------------------
void fitRatio(const TGraphAsymmErrors* gRatio, const TGraphAsymmErrors* gRatioBand, TF1* fit, double &ratio, double &stat, double &systDown, double &systUp, TGraphAsymmErrors* &gRatioStatBand) {
  assert( gRatio->GetN() == gRatioBand->GetN() );

  // Systematic up and down shift of points (keep stat errors)
  TGraphAsymmErrors* gRatioFit = static_cast<TGraphAsymmErrors*>(gRatio->Clone());
  TGraphAsymmErrors* gUp = static_cast<TGraphAsymmErrors*>(gRatio->Clone());
  TGraphAsymmErrors* gDown = static_cast<TGraphAsymmErrors*>(gRatio->Clone());
  for(int i = 0; i < gRatio->GetN(); ++i) {
    gUp->SetPoint(i,gUp->GetX()[i],gUp->GetY()[i]+gRatioBand->GetEYhigh()[i]);
    gDown->SetPoint(i,gDown->GetX()[i],gDown->GetY()[i]-gRatioBand->GetEYlow()[i]);
  }

  gUp->Fit(fit,"0Q");
  systUp = fit->GetParameter(0);
  gDown->Fit(fit,"0Q");
  systDown = fit->GetParameter(0);
  gRatioFit->Fit(fit,"0Q");
  ratio = fit->GetParameter(0);
  stat =  fit->GetParError(0);
  systUp = std::abs(systUp-ratio);
  systDown = std::abs(systDown-ratio);

  double xMin = fit->GetXmin();
  double xMax = fit->GetXmax();
  double x = 0.5*(xMin+xMax);
  double xe = std::abs(x-xMin);
  gRatioStatBand = new TGraphAsymmErrors(1,&x,&ratio,&xe,&xe,&stat,&stat);
  gRatioStatBand->SetFillStyle(3013);
  gRatioStatBand->SetFillColor(fit->GetLineColor());
  gRatioStatBand->SetLineColor(gRatioStatBand->GetFillColor());

  delete gRatioFit;
  delete gUp;
  delete gDown;
}
