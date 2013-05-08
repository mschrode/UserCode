// $Id: $

#ifndef RESOLUTION_TAILS_ETAPTBIN
#define RESOLUTION_TAILS_ETAPTBIN

#include <vector>
#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"

#include "LittleHelper.h"
#include "Output.h"
#include "Pt3Bin.cc"
#include "Style.h"
#include "../sampleTools/BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"


namespace resolutionTails {
  // ------------------------------------------------------------------------------------
  class EtaPtBin {
  public:
    EtaPtBin(unsigned int etaBinIdx, unsigned int ptBinIdx, double nSigCore, double coreScaleFactor,const sampleTools::BinningAdmin* adm, Output* out, const Style* style);
    ~EtaPtBin();


    unsigned int nPt3Bins() const { return pt3Bins_.size(); }
    unsigned int etaBin() const { return etaBin_; }
    unsigned int ptBin() const { return ptBin_; }
    TString toString() const;
    TString toTexTableCells() const;

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
  
    void plotAsymmetryDistributions(double nSigTailStart, double nSigTailEnd) const;
    void plotSpectra() const;
    void plotExtrapolation() const;
    void plotMCTruth() const;
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
    const Style* style_;

    Output* out_;

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

    TString binId() const { return "EtaBin"+util::toTString(etaBin_)+"_PtBin"+util::toTString(ptBin_); }
    TString binId(unsigned int pt3BinIdx) const { return binId()+"_Pt3Bin"+util::toTString(pt3BinIdx); }
    int nPtAsymBins() const;
  };


  // ------------------------------------------------------------------------------------
  EtaPtBin::EtaPtBin(unsigned int etaBinIdx, unsigned int ptBinIdx, double nSigCore, double coreScaleFactor, const sampleTools::BinningAdmin* binAdm, Output* out, const Style* style)
    : etaBin_(etaBinIdx), ptBin_(ptBinIdx),
      nSigCore_(nSigCore), coreScalingFactor_(coreScaleFactor),
      exMin_(0.), exMax_(0.19),
      binAdm_(binAdm), style_(style),
      out_(out) {

    if( Output::DEBUG ) std::cout << "EtaPtBin::EtaPtBin(): Entering constructor " << toString() << std::endl;


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

    fExMC_ = new TF1("fExMC_"+binId(),"expo",exMin_,exMax_);
    fExMC_->SetLineWidth(style_->lineWidth());


    fExData_ = static_cast<TF1*>(fExMC_->Clone("fExData_"+binId()));
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

    if( Output::DEBUG ) std::cout << "EtaPtBin::EtaPtBin(): Leaving constructor" << std::endl;
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
  void EtaPtBin::plotAsymmetryDistributions(double nSigTailStart, double nSigTailEnd) const {
    unsigned int pt3Bin = 0;
    for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
	it != pt3Bins_.end(); ++it, ++pt3Bin) {
      (*it)->plotAsymmetryDataMC(binId(pt3Bin));
      (*it)->plotAsymmetryDataMCSmeared(binId(pt3Bin),nSigTailStart,nSigTailEnd);
    }
  }


  // ------------------------------------------------------------------------------------
  void EtaPtBin::plotMCTruth() const {
    if( hasToyMC() ) toyMC_->plotToyAsymmetry(binId());
    if( hasSymMCTruth() ) symMCTruth_->plotSymMCTruthResponse(binId());
  }


  // ------------------------------------------------------------------------------------
  void EtaPtBin::plotSpectra() const {
    unsigned int pt3Bin = 0;
    for(std::vector<Pt3Bin*>::const_iterator it = pt3Bins_.begin();
	it != pt3Bins_.end(); ++it, ++pt3Bin) {
      (*it)->plotSpectra(binId(pt3Bin));
    }
  }


  // ------------------------------------------------------------------------------------
  void EtaPtBin::plotExtrapolation() const {
    if( Output::DEBUG ) std::cout << "Entering EtaPtBin::plotExtrapolation()" << std::endl;

    // Extrapolation MC Closure
    TCanvas* can = util::HistOps::createTCanvas("can","",500,500);
    can->cd();
    TH1* hFrame = new TH1D("hFrame",style_->title()+";Threshold "+util::LabelFactory::pt3RelMax(),1000,0.,style_->pt3PlotMax());
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
      leg->AddEntry(gFTailMC_,style_->labelMCSmear(),"PL");
      leg->AddEntry(gFTailMC_,style_->labelFAsymMC()+"("+util::LabelFactory::pt3RelMax()+"#rightarrow0) = "+util::toTString(extraMC(),4)+" #pm "+util::toTString(extraMCErr(),4),"");
      leg->AddEntry(gFTailToyMC_,style_->labelFAsymToy()+" = "+util::toTString(toyMC(),4)+" #pm "+util::toTString(toyMCErr(),4),"P");

      hFrame->GetYaxis()->SetTitle(style_->labelFAsymMC());
      hFrame->GetYaxis()->SetRangeUser(minY,0.7);
      hFrame->Draw();
      fExMC_->Draw("same");
      if( hasSymMCTruth() ) gFTailMCTruth_->Draw("PE1same");
      gFTailToyMC_->Draw("PE1same");
      gFTailMC_->Draw("PE1same");
      binLabel_->DrawClone("same");
      leg->Draw("same");
      gPad->SetLogy(1);
      out_->toFiles(can,binId()+"_ExtrapolationMCClosure");

      delete leg;
    }

    // Extrapolation MC + Data
    hFrame->GetYaxis()->SetRangeUser(minY,maxY);
    hFrame->GetYaxis()->SetTitle(style_->labelFAsym());
    TLegend* leg = util::LabelFactory::createLegendWithOffset(5,binLabel_->GetSize());
    leg->AddEntry(gFTailData_,style_->labelData(),"PL");
    leg->AddEntry(gFTailData_,style_->labelFAsymData()+"("+util::LabelFactory::pt3RelMax()+"#rightarrow0) = "+util::toTString(extraData(),4)+" #pm "+util::toTString(extraDataErr(),4),"");
    leg->AddEntry(gFTailMC_,style_->labelMCSmear(),"PL");
    leg->AddEntry(gFTailMC_,style_->labelFAsymMC()+"("+util::LabelFactory::pt3RelMax()+"#rightarrow0) = "+util::toTString(extraMC(),4)+" #pm "+util::toTString(extraMCErr(),4),"");
    leg->AddEntry(gFTailMCGauss_,style_->labelMCGauss(),"P");

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
    out_->toFiles(can,binId()+"_Extrapolation");
    delete leg;
    gPad->SetLogy(0);

    // Spread of ftail for mc
    double relErr = extraMCErr_/extraMC();
    hFrame->GetYaxis()->SetTitle("( "+style_->labelFAsymMC()+" - Fit ) / Fit");
    hFrame->GetYaxis()->SetRangeUser(std::min(-8.*relErr,0.),15.*relErr);
    for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
      hFrame->SetBinContent(i,0.);
    }
    hFrame->SetLineStyle(2);
    leg = util::LabelFactory::createLegendCol(2,style_->legWidth());
    leg->AddEntry(gFTailSpreadMC_,"Asymmetry MC","P");
    leg->AddEntry(hFrame,"Fit","L");

    can->cd();
    hFrame->Draw();
    gFTailSpreadMC_->Draw("PE1same");
    binLabel_->DrawClone("same");
    leg->Draw("same");
    out_->toFiles(can,binId()+"_SpreadMC");  
    delete leg;

    // Spread of ftail for data
    relErr = extraDataErr_/extraData();
    hFrame->GetYaxis()->SetRangeUser(std::min(-8.*relErr,0.),15.*relErr);
    hFrame->GetYaxis()->SetTitle("( "+style_->labelFAsymData()+" - Fit ) / Fit");
    for(int i = 1; i <= hFrame->GetNbinsX(); ++i) {
      hFrame->SetBinContent(i,0.);
    }
    leg = util::LabelFactory::createLegendCol(2,style_->legWidth());
    leg->AddEntry(gFTailSpreadData_,"Asymmetry Data","P");
    leg->AddEntry(hFrame,"Fit","L");

    can->cd();
    hFrame->Draw();
    gFTailSpreadData_->Draw("PE1same");
    binLabel_->DrawClone("same");
    leg->Draw("same");
    out_->toFiles(can,binId()+"_SpreadData");  
    delete leg;

    delete hFrame;
    delete can;

    if( Output::DEBUG ) std::cout << "Leaving EtaPtBin::plotExtrapolation()" << std::endl;
  }


  // ------------------------------------------------------------------------------------
  void EtaPtBin::plotTailStart() const {
    TH1* hTailStartBins = new TH1I("hTailStartBins_EtaBin"+util::toTString(etaBin())+"_PtBin"+util::toTString(ptBin()),style_->title()+";Pt3Bin",nPt3Bins(),0,nPt3Bins());
    TH1* hTailEndBins = static_cast<TH1*>(hTailStartBins->Clone("hTailEndBins_EtaBin"+util::toTString(etaBin())+"_PtBin"+util::toTString(ptBin())));
    for(unsigned int i = 0; i < nPt3Bins(); ++i) {
      hTailStartBins->SetBinContent(1+i,pt3Bins_.at(i)->tailStartBin());
      hTailEndBins->SetBinContent(1+i,pt3Bins_.at(i)->tailEndBin());
    }
    out_->toFiles(hTailStartBins);
    out_->toFiles(hTailEndBins);
    delete hTailStartBins;
    delete hTailEndBins;
  }


  // ------------------------------------------------------------------------------------
  void EtaPtBin::addPt3Bin(unsigned int pt3Bin, double thres, const TString &fileNameData, const TString &fileNameMC) {
    pt3Bins_.push_back(new Pt3Bin(fileNameData,fileNameMC,etaBin_,ptBin_,pt3Bin,thres,nSigCore_,coreScalingFactor_,tailMinBin_,tailMaxBin_,binAdm_,out_,style_));
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
    toyMC_ = new Pt3Bin(fileName,etaBin_,ptBin_,nSigCore_,coreScalingFactor_,nPtAsymBins(),tailMinBin_,tailMaxBin_,binAdm_,out_,style_);
    hasToyMC_ = true;
  }


  // Find min and max asymmetry
  // The borders are specified in number of sigmas, where sigma is the Gaussian
  // width of the *smeared* asymmetry distribution from fileName
  // ------------------------------------------------------------------------------------
  void EtaPtBin::findWindow(const TString &fileName, unsigned int pt3Bin, double nSigTailWindowMin, double nSigTailWindowMax) {
    if( Output::DEBUG ) std::cout << "EtaPtBin::findWindow(): Entering" << std::endl;

    if( Output::DEBUG ) std::cout << "  Getting asymmetry distributions from file . . . " << std::flush;
    TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin_)+"_Pt"+util::toTString(ptBin_)+"_PtSoft"+util::toTString(pt3Bin);

    TH1* hOrig = util::FileOps::readTH1(fileName,histName);
    hOrig->GetXaxis()->SetRangeUser(-1.,1);
    if( Output::DEBUG ) std::cout << "ok" << std::endl;

    if( Output::DEBUG ) std::cout << "  Smearing asymmetry distributions . . . " << std::flush;

    TH1* h = 0;
    double width = 0.;
    double widthOrig = 0.;
    LittleHelper::correctAsymmetryWidth(hOrig,nSigCore_,coreScalingFactor_,h,widthOrig,width);
    delete hOrig;
    if( Output::DEBUG ) std::cout << "ok" << std::endl;

    if( Output::DEBUG ) std::cout << "  Determining tail borders . . . " << std::flush;
    tailWindowMin_ = nSigTailWindowMin*width;
    tailWindowMax_ = nSigTailWindowMax*width;
    tailMinBin_ = h->FindBin(std::abs(tailWindowMin_));
    tailMaxBin_ = h->FindBin(std::abs(tailWindowMax_));
    tailWindowMinEff_ = h->GetXaxis()->GetBinLowEdge(tailMinBin_);
    tailWindowMaxEff_ = h->GetXaxis()->GetBinUpEdge(tailMaxBin_);
    delete h;
    if( Output::DEBUG ) std::cout << "ok" << std::endl;
  
    binLabel_->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+style_->labelWindow(nSigTailWindowMin,nSigTailWindowMax));

    if( Output::DEBUG ) std::cout << "EtaPtBin::findWindow(): Leaving" << std::endl;
  }


  // Set window borders as specified in previous run and stored
  // in histograms; nSigTailWindow??? are only needed for the
  // labels
  // ------------------------------------------------------------------------------------
  void EtaPtBin::setWindow(const TString &fileName, double nSigTailWindowMin, double nSigTailWindowMax) {
    if( Output::DEBUG ) std::cout << "EtaPtBin::setWindow(): Entering" << std::endl;

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
  
    binLabel_->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+style_->labelWindow(nSigTailWindowMin,nSigTailWindowMax));

    if( Output::DEBUG ) std::cout << "EtaPtBin::setWindow(): Leaving" << std::endl;
  }



  // ------------------------------------------------------------------------------------
  void EtaPtBin::extrapolate(double minPt3Data, bool fixDataShape, bool useExtrapolatedValue, bool mcTruthRef) {
    if( Output::DEBUG ) std::cout << "Entering EtaPtBin::extrapolate()" << std::endl;

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
    style_->applyToMC(gFTailMC_);

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
    gFTailMCGauss_->SetMarkerColor(style_->colorGauss());
    gFTailMCGauss_->SetMarkerSize(style_->markerSize());
    gFTailMCGauss_->SetLineWidth(style_->lineWidth());
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
    fExMC_->SetLineColor(style_->colorLineAsymSmear());

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
    style_->applyToMC(gFTailSpreadMC_);

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
    style_->applyToData(gFTailData_);


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
    style_->applyToData(gFTailSpreadData_);

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

    if( Output::DEBUG ) std::cout << "Leaving EtaPtBin::extrapolate()" << std::endl;
  }


  // ------------------------------------------------------------------------------------
  TString EtaPtBin::toString() const {
    TString id = "(";
    id += etaBin_;
    id += ",";
    id += ptBin_;
    id += ")";

    return id;
  }


  // ------------------------------------------------------------------------------------
  TString EtaPtBin::toTexTableCells() const {
    char num[10];

    TString cell = "    $";
    sprintf(num,"%.1f",binAdm_->etaMin(etaBin_));
    cell += num;
    cell += "$ -- $";
    sprintf(num,"%.1f",binAdm_->etaMax(etaBin_));
    cell += num;
    cell += "$ & $";

    double min = binAdm_->ptMin(etaBin_,ptBin_);
    if( min < 1000 ) cell += " ";
    if( min <  100 ) cell += " ";
    sprintf(num,"%.0f",min);
    cell += num;
    cell += "$ & $";

    double max = binAdm_->ptMax(etaBin_,ptBin_);
    if( max < 1000 ) cell += " ";
    if( max <  100 ) cell += " ";
    sprintf(num,"%.0f",max);
    cell += num;
    cell += "$";

    return cell;
  }
}
#endif
