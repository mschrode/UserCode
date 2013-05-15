// $Id: Pt3Bin.cc,v 1.1 2013/05/08 13:07:29 mschrode Exp $

#ifndef RESOLUTION_TAILS_PT3BIN
#define RESOLUTION_TAILS_PT3BIN


#include <cmath>
#include <iostream>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TString.h"

#include "LittleHelper.h"
#include "Output.h"
#include "Style.h"
#include "../sampleTools/BinningAdmin.h"
#include "../sampleTools/HistNames.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"



namespace resolutionTails {
  // ------------------------------------------------------------------------------------
  class Pt3Bin {
  public:
    Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Threshold, double nSigCore, double coreScaleFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, Output* out, const Style* style);
    Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, Output* out, const Style* style);
    ~Pt3Bin();

    TString binLabel() const;
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

    void plotAsymmetryDataMC(const TString &outNameId) const;
    void plotAsymmetryDataMCSmeared(const TString &outNameId, double nSigTailStart, double nSigTailEnd) const;
    void plotToyAsymmetry(const TString &outNameId) const;
    void plotSpectra(const TString &outNameId) const;


  private:
    const unsigned int etaBin_;
    const unsigned int ptBin_;
    const unsigned int pt3Bin_;
    const double pt3Thres_;
    const Style* style_;

    sampleTools::HistNames hName_;
    Output* out_;

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
    TF1* fGaussMCTruth_;
    TH1* hPtAveSpecData_;
    TH1* hPtAveSpecMC_;

    TPaveText* binLabel_;

    TH1* readHistAsym(const TString &fileName, const TString &id) const;  
    TH1* readMCTruthResponse(const TString &fileName) const;
    void initBinLabel(const sampleTools::BinningAdmin* adm, bool isMCTruth = false);
    void getFTail(const TH1* h, double entries, int start, int end, double &fTail, double &fTailErr) const;
    TString inputDirName() const;
    TString inputHistName(const TString &hName) const;
  };


  // Constructor for asymmetry
  // ------------------------------------------------------------------------------------
  Pt3Bin::Pt3Bin(const TString &fileNameData, const TString &fileNameMC, unsigned int etaBin, unsigned int ptBin, unsigned int pt3Bin, double pt3Threshold, double nSigCore, double coreScaleFactor, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, Output* out, const Style* style)
    : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(pt3Bin), pt3Thres_(pt3Threshold), style_(style), out_(out),
      coreScalingFactor_(coreScaleFactor), 
      tailStartBin_(windowMinBin), tailEndBin_(windowMaxBin) {

    if( Output::DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Entering constructor for asymmetry " << binLabel() << std::endl;


    // Non-read distributions
    hResp_ = 0;
    fGaussMCTruth_ = 0;

    // Get spectra
    if( Output::DEBUG ) std::cout << "  Getting pt spectra from file . . . " << std::flush;
    hPtAveSpecData_ = util::FileOps::readTH1(fileNameData,inputDirName()+inputHistName(hName_.ptAve()),"Data"+inputHistName(hName_.ptAve()));
    util::HistOps::setAxisTitles(hPtAveSpecData_,"p^{ave}_{T}","GeV","events");      
    hPtAveSpecData_->SetTitle(style_->title());
    style_->applyToDataMarker(hPtAveSpecData_);
    ptAveMeanData_ = hPtAveSpecData_->GetMean();
    ptAveMeanDataErr_ = hPtAveSpecData_->GetMeanError();
    hPtAveSpecMC_ = util::FileOps::readTH1(fileNameMC,inputDirName()+inputHistName(hName_.ptAve()),"MC"+inputHistName(hName_.ptAve()));
    util::HistOps::setAxisTitles(hPtAveSpecMC_,"p^{ave}_{T}","GeV","events");      
    hPtAveSpecMC_->SetTitle(style_->title());
    style_->applyToMCFilled(hPtAveSpecMC_);
    if( Output::DEBUG ) std::cout << "ok" << std::endl;

    // Get asymmetry distributions
    if( Output::DEBUG ) std::cout << "  Getting asymmetry distributions from file . . . " << std::flush;
    hAsymData_ = readHistAsym(fileNameData,"Data");
    hAsymMC_ = readHistAsym(fileNameMC,"MC");
    hAsymData_->UseCurrentStyle();
    hAsymData_->SetTitle(style_->title());
    hAsymMC_->UseCurrentStyle();
    hAsymMC_->SetTitle(style_->title());
    if( Output::DEBUG ) std::cout << "ok" << std::endl;
  
    // Smear MC asymmetry
    if( Output::DEBUG ) std::cout << "  Smearing asymmetry distribution . . . " << std::flush;
    sig_ = 0.;
    LittleHelper::correctAsymmetryWidth(hAsymMC_,nSigCore,coreScalingFactor_,hAsymMCSmeared_,sig_,sigSmeared_);
    if( Output::DEBUG ) std::cout << "ok" << std::endl;


    // Total entries in asymmetry histograms
    // (to compute statistical uncertainty on fasym)
    // factor 1/2 bc each entry was filled twice
    double nTotalEntriesMC = hAsymMC_->GetEntries()/2.;
    double nTotalEntriesData = hAsymData_->GetEntries()/2.;

    // Get relative number of entries in tail
    if( Output::DEBUG ) std::cout << "  Computing number of tail events . . . " << std::flush;
    getFTail(hAsymMC_,nTotalEntriesMC,tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
    getFTail(hAsymMCSmeared_,nTotalEntriesMC,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
    double minEff = hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_);
    double maxEff = std::min(1.,hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_));
    // 2*OneSidedGaussianTail/TotalGauss
    fNMCSmearedGauss_ = (erf(maxEff/sqrt(2.)/sigSmeared_) - erf(minEff/sqrt(2.)/sigSmeared_));

    getFTail(hAsymData_,nTotalEntriesData,tailStartBin_,tailEndBin_,fNData_,fNDataErr_);
    if( Output::DEBUG ) std::cout << "ok" << std::endl;


    style_->applyToDataMarker(hAsymData_);
    style_->applyToMCFilled(hAsymMC_);
    style_->applyToMCFilled(hAsymMCSmeared_);
    initBinLabel(adm);

    if( Output::DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Leaving constructor for asymmetry" << std::endl;
  }


  // Constructor for toy asymmetry from MC truth response
  //  fAsymMCSmeared from toy MC (MC truth --> asymmetry --> smearing)
  // ------------------------------------------------------------------------------------
  Pt3Bin::Pt3Bin(const TString &fileName, unsigned int etaBin, unsigned int ptBin, double nSigCore, double coreScaleFactor, int nAsymBins, int windowMinBin, int windowMaxBin, const sampleTools::BinningAdmin* adm, Output* out, const Style* style) 
    : etaBin_(etaBin), ptBin_(ptBin), pt3Bin_(1000), pt3Thres_(1000.), style_(style), out_(out),
      coreScalingFactor_(coreScaleFactor),
      tailStartBin_(windowMinBin), tailEndBin_(windowMaxBin) {

    if( Output::DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Entering constructor for toy asymmetry " << binLabel() << std::endl;

    hAsymData_ = 0;
    hAsymMC_ = 0;
    hAsymMCSmeared_ = 0;
    fGaussMCTruth_ = 0;
    hPtAveSpecData_ = 0;
    hPtAveSpecMC_ = 0;
    sig_ = 0.;
    sigSmeared_ = 0.;
    ptAveMeanData_ = 0.;
    ptAveMeanDataErr_ = 0.;
    fNMCSmearedGauss_ = 0.;


    // Get response from file  
    hResp_ = readMCTruthResponse(fileName);
    double entries = hResp_->GetEntries(); // Use entries of original distribution for statistical precision
    hResp_->UseCurrentStyle();
    hResp_->SetMarkerStyle(20);
    util::HistOps::setAxisTitles(hResp_,"Response","","jets");
    hResp_->SetTitle(style_->title());

    // Generate asymmetry distribution from response
    if( Output::DEBUG ) std::cout << "  Generating toy asymmetry distribution  . . .  " << std::flush;
    hAsymMC_ = new TH1D("ToyAsym_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin),style_->title(),nAsymBins,-1.,1.);
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
    if( Output::DEBUG ) std::cout << "ok" << std::endl;
  
    // Smear toy asymmetry
    if( Output::DEBUG ) std::cout << "  Smearing toy asymmetry distribution  . . .  " << std::flush;
    sig_ = 0.;
    LittleHelper::correctAsymmetryWidth(hAsymMC_,nSigCore,coreScalingFactor_,hAsymMCSmeared_,sig_,sigSmeared_);
    if( Output::DEBUG ) std::cout << "ok" << std::endl;

    if( Output::DEBUG ) std::cout << "  Computing number of tail events in toy asymmetry distribution  . . .  " << std::flush;
    // Get relative number of entries in tail
    getFTail(hAsymMC_,entries,tailStartBin_,tailEndBin_,fNMC_,fNMCErr_);
    getFTail(hAsymMCSmeared_,entries,tailStartBin_,tailEndBin_,fNMCSmeared_,fNMCSmearedErr_);
    if( Output::DEBUG ) std::cout << "ok" << std::endl;

    // Set style
    initBinLabel(adm,true);

    if( Output::DEBUG ) std::cout << "Pt3Bin::Pt3Bin(): Leaving constructor for toy asymmetry" << std::endl;
  }


  // ------------------------------------------------------------------------------------
  void Pt3Bin::initBinLabel(const sampleTools::BinningAdmin* adm, bool isMCTruth) {
    if( isMCTruth ) binLabel_ = util::LabelFactory::createPaveText(1);
    else binLabel_ = util::LabelFactory::createPaveText(2);
    binLabel_->AddText(util::LabelFactory::etaCut(adm->etaMin(etaBin_),adm->etaMax(etaBin_))+",  "+util::LabelFactory::ptAveCut(adm->ptMin(etaBin_,ptBin_),adm->ptMax(etaBin_,ptBin_)));
    if( !isMCTruth ) binLabel_->AddText(util::LabelFactory::deltaPhiCut(2.7)+",  "+util::LabelFactory::pt3RelCut(adm->ptSoftMax(pt3Bin_)));
  }


  // ------------------------------------------------------------------------------------
  Pt3Bin::~Pt3Bin() {
    if( hAsymData_ ) delete hAsymData_;
    if( hAsymMC_ ) delete hAsymMC_;
    if( hAsymMCSmeared_ ) delete hAsymMCSmeared_;
    if( hResp_ ) delete hResp_;
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
    TH1* h = util::FileOps::readTH1(fileName,inputDirName()+inputHistName(hName_.ptAsymAbs()),id+inputHistName(hName_.ptAsymAbs()));
    h->GetXaxis()->SetRangeUser(-1.,1);
    h->Scale(2./h->Integral("width"));
    util::HistOps::setAxisTitles(h,"|Asymmetry|","","events",true);      
    h->SetTitle("");

    return h;
  }


  // Get MC truth response histogram from Kalibri input files
  // ------------------------------------------------------------------------------------
  TH1* Pt3Bin::readMCTruthResponse(const TString &fileName) const {
    TString etaBin = "Eta"+util::toTString(etaBin_);
    TString ptBin = "Pt"+util::toTString(ptBin_);
    TString ptSoftBin = "PtSoft0";
    TString dirName = etaBin+"/"+ptBin+"/"+ptSoftBin+"/";
    TString histName = hName_.respGenBin()+"_"+etaBin+"_"+ptBin+"_"+ptSoftBin;
    TString title = "Response";
    TH1* h = util::FileOps::readTH1(fileName,dirName+histName,histName);
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
    leg->AddEntry(hAsymData_,style_->labelData(),"P");
    leg->AddEntry(hAsymMC_,style_->labelMC(),"F");

    // Log scale
    util::HistOps::setYRange(hAsymMC_,5,3E-5);
    util::HistOps::setYRange(hAsymData_,5,3E-5);
    hAsymMC_->GetXaxis()->SetRangeUser(-1.,1.);

    TString canName = outNameId+"_PtAsym";
    TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);
    can->cd();
    hAsymMC_->GetXaxis()->SetRangeUser(0.,1.);
    hAsymMC_->Draw("HISTE");
    hAsymData_->Draw("PE1same");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    out_->toFiles(can,outNameId+"_PtAsymLog");

    // Linear scale
    util::HistOps::setYRange(hAsymMC_,5);
    hAsymMC_->GetXaxis()->SetRangeUser(0.,0.4);
    can->cd();
    hAsymMC_->Draw("HISTE");
    hAsymData_->Draw("PE1same");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(0);
    out_->toFiles(can,outNameId+"_PtAsym");
  
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
    //   out_->toFiles(can,outNameId+"_PtAsymFit");

    // Ratio data / MC
    TH1 *hRatio = util::HistOps::createRatioPlot(hAsymData_,hAsymMC_);
    //   TH1 *hRatioFrame = util::HistOps::createRatioFrame(hAsymData_,"Data / MC",0.,3.);
    //   can->cd();
    //   hRatioFrame->Draw();
    //   hRatio->Draw("PE1same");
    //   label_->info(etaBin_,ptBin_)->Draw("same");
    //   label_->eta(etaBin_)->Draw("same");
    //   can->SetLogy(0);
    //   out_->toFiles(can,outNameId+"_PtAsymLogRatio");

    //   hRatioFrame->GetXaxis()->SetRangeUser(-0.4,0.4);
    //   can->cd();
    //   hRatioFrame->Draw();
    //   hRatio->Draw("PE1same");
    //   label_->info(etaBin_,ptBin_)->Draw("same");
    //   label_->eta(etaBin_)->Draw("same");
    //   can->SetLogy(0);
    //   out_->toFiles(can,outNameId+"_PtAsymRatio");

    // Bottom ratio plot
    delete can;
    can = util::HistOps::createRatioTopCanvas();
    TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
    TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hAsymMC_);
    TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,static_cast<TString>(hAsymMC_->GetXaxis()->GetTitle()),"","#frac{Data}{Sim.}",0.61,1.39);
    bRatioBottomFrame->SetFillColor(0);
    bRatioBottomFrame->SetLineWidth(style_->lineWidth());
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
    out_->toFiles(can,outNameId+"_PtAsymBottomRatio");

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
    leg->AddEntry(hAsymData_,style_->labelData(),"P");
    leg->AddEntry(hAsymMCSmeared_,style_->labelMCSmear(),"F");

    // Log scale
    util::HistOps::setYRange(hAsymMCSmeared_,5,3E-5);
    util::HistOps::setYRange(hAsymData_,5,3E-5);
    hAsymMCSmeared_->GetXaxis()->SetRangeUser(0.,1.);

    TString canName = outNameId+"_PtAsym";
    TCanvas *can = util::HistOps::createTCanvas(canName,canName,500,500);
    can->cd();
    hAsymMCSmeared_->Draw("HISTE");
    hAsymData_->Draw("PE1same");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    out_->toFiles(can,outNameId+"_PtSmearAsymLog");

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
    gauss->SetLineColor(style_->colorGauss());
    gauss->SetFillStyle(style_->hatchStyle());
    gauss->SetFillColor(gauss->GetLineColor());
    gauss->SetRange(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),gauss->GetX(3E-5,0.,1.));
  
    TLegend* legWin = util::LabelFactory::createLegendWithOffset(4,2);
    legWin->AddEntry(hAsymData_,style_->labelData(),"P");
    legWin->AddEntry(hAsymMCSmeared_,style_->labelMCSmear(),"F");
    legWin->AddEntry(gauss,style_->labelMCGauss(),"F");
    legWin->AddEntry(arr,style_->labelWindow(nSigTailStart,nSigTailEnd),"L");

    TPaveText* labelFAsym = util::LabelFactory::createPaveTextWithOffset(4,0.6,2+legWin->GetNRows());
    int nDigits = 4;
    labelFAsym->AddText(style_->labelFAsymData()+" = "+util::toTString(fTailData(),nDigits)+" #pm "+util::toTString(fTailDataErr(),nDigits));
    labelFAsym->AddText(style_->labelFAsymMC()+" = "+util::toTString(fTailMCSmeared(),nDigits)+" #pm "+util::toTString(fTailMCSmearedErr(),nDigits));
    labelFAsym->AddText(style_->labelFAsymGauss()+" = "+util::toTString(fTailMCSmearedGauss(),nDigits));

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
    out_->toFiles(can,outNameId+"_PtSmearAsymTail");
    delete legWin;
    delete arr;

    // Without data to illustrate window definition
    if( pt3Bin_ == 0 ) {
      legWin = util::LabelFactory::createLegendCol(2,style_->legWidth());
      legWin->AddEntry(hAsymMCSmeared_,style_->labelMCSmear(),"F");
      legWin->AddEntry(win,"Tail","F");

      can->cd();
      hAsymMCSmeared_->Draw("HISTE");
      win->Draw("same");
      binLabel_->Draw("same");
      legWin->Draw("same");
      can->SetLogy(1);
      out_->toFiles(can,outNameId+"_WindowDef");
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
    out_->toFiles(can,outNameId+"_PtSmearAsym");


    // Bottom ratio plot
    TH1 *hRatio = util::HistOps::createRatioPlot(hAsymData_,hAsymMCSmeared_);
    delete can;
    can = util::HistOps::createRatioTopCanvas();
    TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
    TH1 *bRatioTopFrame = util::HistOps::createRatioTopHist(hAsymMCSmeared_);
    TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hAsymMC_,static_cast<TString>(hAsymMC_->GetXaxis()->GetTitle()),"","#frac{Data}{Sim.}",0.61,1.39);
    bRatioBottomFrame->SetFillColor(0);
    bRatioBottomFrame->SetLineWidth(style_->lineWidth());
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
    out_->toFiles(can,outNameId+"_PtSmearAsymBottomRatio");

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
    leg->AddEntry(hPtAveSpecData_,style_->labelData()+",  N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
    hPtAveSpecData_->GetXaxis()->SetRange(xMinBin,xMaxBin);
    util::HistOps::setYRange(hPtAveSpecData_,6,3E-1);
    can->cd();
    hPtAveSpecData_->Draw("PE1");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    out_->toFiles(can,outNameId+"_PtAveSpectrumData");

    // Absolute MC spectrum
    delete leg;
    leg = util::LabelFactory::createLegendWithOffset(1,binLabel_->GetSize());
    leg->AddEntry(hPtAveSpecMC_,style_->labelMC()+",  N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
    util::HistOps::setYRange(hPtAveSpecMC_,6,3E-5);
    can->cd();
    hPtAveSpecMC_->GetXaxis()->SetRange(xMinBin,xMaxBin);
    hPtAveSpecMC_->Draw("HIST");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    out_->toFiles(can,outNameId+"_PtAveSpectrumMC");

    // Comparison (MC entries scaled to data)
    delete leg;
    leg = util::LabelFactory::createLegendWithOffset(2,binLabel_->GetSize());
    leg->AddEntry(hPtAveSpecData_,style_->labelData()+",  N = "+util::toTString(hPtAveSpecData_->GetEntries()),"P");
    leg->AddEntry(hPtAveSpecMC_,style_->labelMC()+",  N = "+util::toTString(hPtAveSpecMC_->GetEntries()),"F");
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
      out_->toFiles(can,outNameId+"_PtAveSpectra");
      // set back MC scale
      hPtAveSpecMC_->Scale(scaleMC/scaleData);
    }  
    delete leg;
    delete can;
  }


  // ------------------------------------------------------------------------------------
  void Pt3Bin::plotToyAsymmetry(const TString &outNameId) const {
    if( Output::DEBUG ) std::cout << "Entering Pt3Bin::plotToyAsymmetry()" << std::endl;

    // Log scale
    TBox* win = new TBox(hAsymMCSmeared_->GetXaxis()->GetBinLowEdge(tailStartBin_),0.,
			 hAsymMCSmeared_->GetXaxis()->GetBinUpEdge(tailEndBin_),10.);
    win->SetLineWidth(1);
    win->SetFillStyle(3444);
    win->SetLineColor(kRed);
    win->SetFillColor(win->GetLineColor());

    TLegend* leg = util::LabelFactory::createLegendCol(2,style_->legWidth());
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
    out_->toFiles(can,outNameId+"_ToyMCTail");

    hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);
    can->cd();
    hAsymMCSmeared_->Draw("HIST");
    binLabel_->Draw("same");
    leg->Draw("same");
    can->SetLogy(1);
    out_->toFiles(can,outNameId+"_ToyMCLog");

    // Linear scale
    util::HistOps::setYRange(hAsymMCSmeared_,4);
    hAsymMCSmeared_->GetXaxis()->SetRangeUser(-0.4,0.4);
    can->cd();
    hAsymMCSmeared_->Draw("HIST");
    binLabel_->Draw("same");
    can->SetLogy(0);
    out_->toFiles(can,outNameId+"_ToyMC");

    hAsymMCSmeared_->GetXaxis()->SetRangeUser(-1.,1.);

    // MC truth response
    delete leg;
    leg = util::LabelFactory::createLegendCol(1,style_->legWidth());
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
    out_->toFiles(can,outNameId+"_MCTruthResponseLog");

    // Linear scale
    util::HistOps::setYRange(hResp_,4);
    hResp_->GetXaxis()->SetRangeUser(0.4,1.6);
    can->cd();
    hResp_->Draw("HIST");
    binLabel_->Draw("same");
    can->SetLogy(0);
    out_->toFiles(can,outNameId+"_MCTruthResponse");

    hResp_->GetXaxis()->SetRangeUser(0.,2.);


    delete win;
    delete leg;
    delete can;

    if( Output::DEBUG ) std::cout << "Leaving Pt3Bin::plotToyAsymmetry()" << std::endl;
  }


  TString Pt3Bin::binLabel() const {
    TString id = "(";
    id += etaBin_;
    id += ",";
    id += ptBin_;
    id += ",";
    id += pt3Bin_;
    id += ")";

    return id;
  }


  TString Pt3Bin::inputDirName() const {
    TString etaBin = "Eta"+util::toTString(etaBin_);
    TString ptBin = "Pt"+util::toTString(ptBin_);
    TString ptSoftBin = "PtSoft"+util::toTString(pt3Bin_);

    return etaBin+"/"+ptBin+"/"+ptSoftBin+"/";
  }


  TString Pt3Bin::inputHistName(const TString &hName) const {
    TString etaBin = "Eta"+util::toTString(etaBin_);
    TString ptBin = "Pt"+util::toTString(ptBin_);
    TString ptSoftBin = "PtSoft"+util::toTString(pt3Bin_);

    return hName+"_"+etaBin+"_"+ptBin+"_"+ptSoftBin;
  }

}
#endif
