#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"



class Measurement {
public:
  Measurement(const TString &name, const TF1 *mcTruth, int color, int marker);

  void readData(const TString &fileName) { read(fileName,"Data"); }
  void readMC(const TString &fileName) { read(fileName,"MC"); }

  TGraphErrors *gData() { return g_; }
  TGraphErrors *gDataRatio() { return gRatio_; }
  TGraphErrors *gMC() { return gMC_; }
  TGraphErrors *gMCRatio() { return gMCRatio_; }

  TString name() const { return name_; }
  int nBins() const { return static_cast<int>(x_.size()); }
  void print() const;


private:
  const TString name_;
  const int color_;
  const int marker_;
  
  std::vector<double> x_;
  std::vector<double> xe_;
  std::vector<double> y_;
  std::vector<double> ye_;
  std::vector<double> xmc_;
  std::vector<double> xmce_;
  std::vector<double> ymc_;
  std::vector<double> ymce_;

  TF1 *mcTruth_;
  TGraphErrors *g_;
  TGraphErrors *gRatio_;
  TGraphErrors *gMC_;
  TGraphErrors *gMCRatio_;

  void createTGraphs();
  void read(const TString &fileName, const TString &sample);
};



Measurement::Measurement(const TString &name, const TF1 *mcTruth, int color, int marker)
  : name_(name), color_(color), marker_(marker), g_(0), gMC_(0) {

  // MC truth resolution
  mcTruth_ = static_cast<TF1*>(mcTruth->Clone(name+"_MCTruth"));

  createTGraphs();
}


void Measurement::createTGraphs() {
  if( g_ ) delete g_;
  if( gRatio_ ) delete gRatio_;
  if( gMC_ ) delete gMC_;
  if( gMCRatio_ ) delete gMCRatio_;

  g_ = new TGraphErrors(x_.size(),&(x_.front()),&(y_.front()),&(xe_.front()),&(ye_.front()));
  g_->SetMarkerStyle(marker_);
  g_->SetLineColor(color_);
  g_->SetMarkerColor(color_);
  gRatio_ = util::HistOps::createRatioGraph(g_,mcTruth_);

  gMC_ = new TGraphErrors(xmc_.size(),&(xmc_.front()),&(ymc_.front()),&(xmce_.front()),&(ymce_.front()));
  gMC_->SetMarkerStyle(marker_+4);
  gMC_->SetLineColor(color_);
  gMC_->SetMarkerColor(color_);
  gMCRatio_ = util::HistOps::createRatioGraph(gMC_,mcTruth_);
}



void Measurement::read(const TString &fileName, const TString &sample) {
  
  std::vector<double> *x = 0;
  std::vector<double> *xe = 0;
  std::vector<double> *y = 0;
  std::vector<double> *ye = 0;

  if( sample == "Data" ) {
    x = &x_;
    xe = &xe_;
    y = &y_;
    ye = &ye_;
  } else if( sample == "MC" ) {
    x = &xmc_;
    xe = &xmce_;
    y = &ymc_;
    ye = &ymce_;
  } else {
    std::cerr << "ERROR: Unknown tag '" << sample << "'" << std::endl;
    exit(-1);
  }

  std::ifstream file;
  file.open(fileName.Data());
  if( file.is_open() ) {
    double tmp = 0.;
    while( !file.eof() ) {
      file >> tmp;
      if( tmp != ( ye->size() ? ye->back() : -1. ) ) {
	x->push_back(tmp);
	file >> tmp;
	xe->push_back(tmp);
	file >> tmp;
	y->push_back(tmp);
	file >> tmp;
	ye->push_back(tmp);
      }
    }
    createTGraphs();
  } else {
    std::cerr << "ERROR opening file '" << fileName << "'\n";
    exit(1);
  }
  file.close();
}


void Measurement::print() const {
  std::cout << "\n\n---------------------------------------------------------------\n";
  std::cout << name() << " (data)" << std::endl;
  for(size_t i = 0; i < x_.size(); ++i) {
    std::cout << x_[i] << " +/- " << xe_[i] << ":  ";
    std::cout << y_[i] << " +/- " << ye_[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << name() << " (MC)" << std::endl;
  for(size_t i = 0; i < xmc_.size(); ++i) {
    std::cout << xmc_[i] << " +/- " << xmce_[i] << ":  ";
    std::cout << ymc_[i] << " +/- " << ymce_[i] << std::endl;
  }
  std::cout << "---------------------------------------------------------------\n\n";
}



void plotCombinedResolution() {
  util::StyleSettings::presentationNoTitle();

  double xMin = 20.;
  double xMax = 250.;


  // MC truth resolution
  TF1 *mcTruth = new TF1("mcTruth","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",xMin,xMax);
  mcTruth->SetLineWidth(2);
  mcTruth->SetLineStyle(2);
  mcTruth->SetParameter(0,3.12704);
  mcTruth->SetParameter(1,1.16413);
  mcTruth->SetParameter(2,0.0344903);
  

  Measurement *mPtBal = new Measurement("#gamma + jet balance",mcTruth,kBlue+1,21);
  mPtBal->readData("resGauss_PtBalance_2010-10-17_Data");
  mPtBal->readMC("resGauss_PtBalance_2010-10-17_MC");
  mPtBal->print();

  Measurement *mMPF = new Measurement("#gamma + jet MPF",mcTruth,kGreen+3,22);
  mMPF->readData("resGauss_MPF_2010-10-17_Data");
  mMPF->readMC("resGauss_MPF_2010-10-17_MC");
  mMPF->print();

  Measurement *mMaxLike = new Measurement("Dijet max likelihood",mcTruth,kRed,20);
  mMaxLike->readData("resGauss_MaxLike_2010-10-17_Data");
  mMaxLike->readMC("resGauss_MaxLike_2010-10-17_MC");
  mMaxLike->print();

//   TPaveText *txt = util::LabelFactory::createPaveText(1,-0.2);
//   txt->AddText("|#eta| < 1.3");

//   TLegend *leg = util::LabelFactory::createLegendCol(7,0.8);
//   leg->AddEntry(g[2],"#gamma + jet balance (data)","P");
//   leg->AddEntry(g[3],"#gamma + jet balance (MC)","P");
//   leg->AddEntry(g[4],"#gamma + jet MPF (data)","P");
//   leg->AddEntry(g[5],"#gamma + jet MPF (MC)","P");
//   leg->AddEntry(g[0],"Dijet max likelihood (data)","P");
//   leg->AddEntry(g[1],"Dijet max likelihood (MC)","P");
//   leg->AddEntry(mc,"MC truth","L");



  TCanvas *can0 = new TCanvas("can0","Resolution",500,500);
  can0->cd();
  TH1 *hFrame = new TH1D("hFrame",";p_{T} (GeV);#sigma / p_{T}",1000,xMin,xMax);
  hFrame->GetYaxis()->SetRangeUser(0.,1.);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->Draw();
  mcTruth->Draw("same");
  mPtBal->gData()->Draw("PE1same");
  mPtBal->gMC()->Draw("PE1same");  
  mMPF->gData()->Draw("PE1same");
  mMPF->gMC()->Draw("PE1same");  
  mMaxLike->gData()->Draw("PE1same");
  mMaxLike->gMC()->Draw("PE1same");  
  can0->SetLogx();

  TCanvas *can1 = new TCanvas("can1","Ratio Meas / MCTruth",500,500);
  can1->cd();
  TH1 *hRatioFrame = util::HistOps::createRatioFrame(hFrame,"#sigma(Meas) / #sigma(MCTruth)",0.6,2.);
  hRatioFrame->GetXaxis()->SetMoreLogLabels();
  hRatioFrame->Draw();
  mPtBal->gDataRatio()->Draw("PE1same");
  mPtBal->gMCRatio()->Draw("PE1same");  
  mMPF->gDataRatio()->Draw("PE1same");
  mMPF->gMCRatio()->Draw("PE1same");  
  mMaxLike->gDataRatio()->Draw("PE1same");
  mMaxLike->gMCRatio()->Draw("PE1same");  
  can1->SetLogx();

  TCanvas *bRatioTopCan = util::HistOps::createRatioTopCanvas();
  TPad *bRatioBottomPad = util::HistOps::createRatioBottomPad();
  TH1 *bRatioTopFrame = util::HistOps::createRatioTopFrame(hFrame);
  TH1 *bRatioBottomFrame = util::HistOps::createRatioBottomFrame(hFrame,"p_{T}","GeV",0.61,1.81);
  
  bRatioTopCan->cd();
  bRatioTopFrame->Draw();
  mcTruth->Draw("same");
  mPtBal->gData()->Draw("PE1same");
  mPtBal->gMC()->Draw("PE1same");  
  mMPF->gData()->Draw("PE1same");
  mMPF->gMC()->Draw("PE1same");  
  mMaxLike->gData()->Draw("PE1same");
  mMaxLike->gMC()->Draw("PE1same");  
  gPad->SetLogx();
  bRatioBottomPad->Draw();
  bRatioBottomPad->cd();
  bRatioBottomFrame->GetXaxis()->SetMoreLogLabels();
  bRatioBottomFrame->Draw();
  mPtBal->gDataRatio()->Draw("PE1same");
  mPtBal->gMCRatio()->Draw("PE1same");  
  mMPF->gDataRatio()->Draw("PE1same");
  mMPF->gMCRatio()->Draw("PE1same");  
  mMaxLike->gDataRatio()->Draw("PE1same");
  mMaxLike->gMCRatio()->Draw("PE1same");  
  gPad->SetLogx();



  // Uncertainties on MaxLike
  TGraphErrors *gMaxLikeRatioData = mMaxLike->gDataRatio();
  
  // Mean scaling factor
  gMaxLikeRatioData->Fit("pol0","0Q");
  double maxLikeScaleMean = gMaxLikeRatioData->GetFunction("pol0")->GetParameter(0);
  double maxLikeScaleUncert = 0.5*(maxLikeScaleMean-1.);
  double maxLikeScaleHigh = gMaxLikeRatioData->GetY()[0];
  double maxLikeScaleLow = gMaxLikeRatioData->GetY()[0];
  for(int i = 1; i < gMaxLikeRatioData->GetN(); ++i) {
    if( gMaxLikeRatioData->GetY()[i] > maxLikeScaleHigh ) maxLikeScaleHigh = gMaxLikeRatioData->GetY()[i];
    else if( gMaxLikeRatioData->GetY()[i] < maxLikeScaleLow ) maxLikeScaleLow = gMaxLikeRatioData->GetY()[i];
  }
  TH1 *hMaxLikeScale = new TH1D("hMaxLikeScale",";p_{T} (GeV);#sigma(Meas) / #sigma(MCTruth)",1000,xMin,xMax);
  hMaxLikeScale->GetYaxis()->SetRangeUser(0.6,2.);
  hMaxLikeScale->GetXaxis()->SetMoreLogLabels();
  hMaxLikeScale->SetLineColor(gMaxLikeRatioData->GetLineColor());
  TH1 *hMaxLikeScaleUncert = static_cast<TH1D*>(hMaxLikeScale->Clone("hMaxLikeScaleUncert"));
  hMaxLikeScaleUncert->SetFillColor(kRed+2);
  TH1 *hMaxLikeScaleHighLow = static_cast<TH1D*>(hMaxLikeScaleUncert->Clone("hMaxLikeScaleHighLow"));

  for(int bin = 1; bin <= hMaxLikeScale->GetNbinsX(); ++bin) {
    hMaxLikeScale->SetBinContent(bin,maxLikeScaleMean);

    hMaxLikeScaleUncert->SetBinContent(bin,maxLikeScaleMean);
    hMaxLikeScaleUncert->SetBinError(bin,maxLikeScaleUncert);

    hMaxLikeScaleHighLow->SetBinContent(bin,maxLikeScaleLow+0.5*(maxLikeScaleHigh-maxLikeScaleLow));
    hMaxLikeScaleHighLow->SetBinError(bin,0.5*(maxLikeScaleHigh-maxLikeScaleLow));
  }

  TCanvas *canMaxLikeExtra1 = new TCanvas("canMaxLikeExtra1","MaxLike: Ratio Meas / MCTruth",500,500);
  canMaxLikeExtra1->cd();
  hMaxLikeScaleHighLow->Draw("E3");
  hMaxLikeScale->Draw("same");
  hRatioFrame->Draw("same");
  mMaxLike->gDataRatio()->Draw("PE1same");
  canMaxLikeExtra1->SetLogx();


  // Extrapolated resolution
  TH1 *hMaxLikeExtra = static_cast<TH1D*>(hFrame->Clone("hMaxLikeExtra"));
  hMaxLikeExtra->SetLineColor(gMaxLikeRatioData->GetLineColor());
  TH1 *hMaxLikeExtraUncert = static_cast<TH1D*>(hMaxLikeExtra->Clone("hMaxLikeExtraUncert"));
  hMaxLikeExtraUncert->SetFillColor(hMaxLikeScaleHighLow->GetFillColor());
  for(int bin = 1; bin <= hMaxLikeExtra->GetNbinsX(); ++bin) {
    double mc = mcTruth->Eval(hMaxLikeExtra->GetBinCenter(bin));
    hMaxLikeExtra->SetBinContent(bin,mc*hMaxLikeScale->GetBinContent(bin));
    hMaxLikeExtraUncert->SetBinContent(bin,mc*hMaxLikeScaleHighLow->GetBinContent(bin));
    hMaxLikeExtraUncert->SetBinError(bin,mc*hMaxLikeScaleHighLow->GetBinError(bin));
  }

  TCanvas *canMaxLikeExtra2 = new TCanvas("canMaxLikeExtra2","MaxLike: Extrapolation",500,500);
  canMaxLikeExtra2->cd();
  hMaxLikeExtraUncert->Draw("E3");
  hMaxLikeExtra->Draw("same");
  mMaxLike->gData()->Draw("PE1same");
  mcTruth->Draw("same");
  canMaxLikeExtra2->SetLogx();
}
