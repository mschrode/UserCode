// $Id: plotCombinedSpectra.C,v 1.8 2011/12/03 13:20:26 mschrode Exp $

// Compare HT and MHT spectra in data with bkg. prediction
//
// The statistical methods used to obtain the approximated
// combined and individual uncertainties are documented in
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/mschrode/notes/AdditionalDocumentation_CMS-SUS-11-004/RA2DataVsCombinedBkgPrediction.pdf?view=log

#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TPad.h"
#include "TString.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


// Event yields and uncertainties
// Baseline 4.65/fb 2011
const double BKG_ZTOINV        = 591.2807;
const double BKG_ZTOINV_STAT   = 16.;
const double BKG_ZTOINV_SYSTUP = 81.;
const double BKG_ZTOINV_SYSTDN = 81.;

const double BKG_LL        = 426.8;
const double BKG_LL_STAT   = 22.4;
const double BKG_LL_SYSTUP = 55.2;
const double BKG_LL_SYSTDN = 53.7;

const double lumiScaleHadTau   = 4650./3943.;
const double BKG_HADTAU        = 427.*lumiScaleHadTau;
const double BKG_HADTAU_STAT   = 10.8*sqrt(lumiScaleHadTau);
const double BKG_HADTAU_SYSTUP = 36.6*sqrt(lumiScaleHadTau);
const double BKG_HADTAU_SYSTDN = 34.9*sqrt(lumiScaleHadTau);

const double BKG_QCD        = 349.02;
const double BKG_QCD_STAT   = 21.38;
const double BKG_QCD_SYSTUP = 188.03;
const double BKG_QCD_SYSTDN = 186.43;

const double BKG_TOTAL = BKG_ZTOINV +
                         BKG_LL  +
                         BKG_HADTAU +
			 BKG_QCD;

double tmp1 = sqrt( pow(BKG_ZTOINV_STAT,2) + pow(0.5*(BKG_ZTOINV_SYSTUP+BKG_ZTOINV_SYSTDN),2) );
double tmp2 = sqrt( pow(BKG_LL_STAT,2) + pow(0.5*(BKG_LL_SYSTUP+BKG_LL_SYSTDN),2) );
double tmp3 = sqrt( pow(BKG_HADTAU_STAT,2) + pow(0.5*(BKG_HADTAU_SYSTUP+BKG_HADTAU_SYSTDN),2) );
double tmp4 = sqrt( pow(BKG_QCD_STAT,2) + pow(0.5*(BKG_QCD_SYSTUP+BKG_QCD_SYSTDN),2) );

const double BKG_TOTAL_UNCERT = sqrt( pow(tmp1,2) + pow(tmp2,2) + pow(tmp3,2) + pow(tmp4,2) );






enum Bkg { ZInv, LostLepton, HadronicTau, QCD };

class HistContainer {
public:
  HistContainer(const TString &id, const TString &dataFileName, const TString &dataHistName, const TString &xTitle, const TString &xUnit, int rebin, double xMax)
    : id_(id), rebin_(rebin), xMax_(xMax) {
    
    std::cout << "Constructing spectra '" << id_ << "'" << std::endl;

    data_ = util::FileOps::readTH1(dataFileName,dataHistName,id_+"data");
    util::HistOps::shiftStartOfXAxisToFirstPopulatedBin(data_);
    data_->UseCurrentStyle();
    data_->SetTitle(util::StyleSettings::title(4600,true)); // Lumi rounded to one digit
    data_->SetMarkerStyle(20);
    data_->SetNdivisions(505);
    if( rebin_ > 1 ) data_->Rebin(rebin_);
    data_->GetXaxis()->SetRange(1,data_->FindBin(xMax_));
    util::HistOps::setAxisTitles(data_,xTitle,xUnit,"events");
    util::HistOps::setYRange(data_,1,3E-1);

    bkgStack_ = new THStack(id_+"BkgStack","");
    totalBkg_ = 0;
  }


  TString bkgLabel(Bkg bkg) const {
    TString label = "";
    if( bkg == ZInv ) label = "Z#rightarrow#nu#bar{#nu}+jets";
    else if( bkg == LostLepton ) label = "W/t#bar{t}#rightarrowe/#mu+X";
    else if( bkg == HadronicTau ) label = "W/t#bar{t}#rightarrow#tau_{h}+X";
    else if( bkg == QCD ) label = "QCD";

    return label;
  }


  int bkgFillColor(Bkg bkg) const {
    int color = kBlack;
    if( bkg == ZInv ) color = kGreen+1;
    else if( bkg == HadronicTau ) color = kYellow;
    else if( bkg == LostLepton ) color = kRed+1;
    else if( bkg == QCD ) color = kRed+3;

    return color;
  }


  TLegend* legend() const {
    //TLegend* leg = util::LabelFactory::createLegendCol(bkgs_.size()+1,0.6);
    TLegend* leg = util::LabelFactory::createLegendColWithOffset(bkgs_.size(),0.6,0.06);
    //leg->AddEntry(data_,"Data ("+util::toTString(data_->Integral(),0)+")","P");
    //    leg->AddEntry(data_,"Data","P");
    //util::LabelFactory::addExtraLegLine(leg,"");
    for(std::vector<Bkg>::const_reverse_iterator rit = bkgs_.rbegin();
	rit != bkgs_.rend(); ++rit) {
      const TH1* h = bkgHists_.find(*rit)->second;
      TString numEvts = util::toTString(h->Integral(1,10000),1);	// Include overflow bins
//       if( *rit == LostLepton ) numEvts = util::toTString(BKG_LL);
//       else if( *rit == HadronicTau ) numEvts = util::toTString(BKG_HADTAU);
//      else if( *rit == ZInv ) numEvts = util::toTString(BKG_ZTOINV);
//      else if( *rit == QCD ) numEvts = util::toTString(BKG_QCD);
      //leg->AddEntry(h,bkgLabel(*rit)+" ("+numEvts+")","F");
      leg->AddEntry(h,bkgLabel(*rit),"F");
    }

    return leg;
  }


  TGraphAsymmErrors* totalBkg() const {
    return util::HistOps::getUncertaintyBand(totalBkg_,kBlue+2,3004);
  }


  TGraphAsymmErrors* bkgForRatio() const {
    TH1* h = static_cast<TH1*>(totalBkg_->Clone("tmp"));
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      double relUncert = h->GetBinError(bin);
      if( h->GetBinContent(bin) > 0. ) relUncert /= h->GetBinContent(bin);
      h->SetBinError(bin,relUncert);
      h->SetBinContent(bin,1.);
    }
    TGraphAsymmErrors* g = util::HistOps::getUncertaintyBand(h,kBlue+2,3004);
    delete h;

    return g;
  }


  TH1* dataForRatio() const {
    TH1* h = static_cast<TH1*>(data_->Clone(id_+"dataInRatio"));
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      double relVal = h->GetBinContent(bin);
      double relUncert = h->GetBinError(bin);
      double pred = totalBkg_->GetBinContent(bin);
      if( pred > 0. ) {
	relVal /= pred;
	relUncert /= pred;
      }
      h->SetBinContent(bin,relVal);
      h->SetBinError(bin,relUncert);
    }

    return h;
  }


  void addBkg(Bkg bkg, const TString &fileName, const TString &histName, double scale = 1.) {
    bkgs_.push_back(bkg);

    TH1* h = util::FileOps::readTH1(fileName,histName,id_+"Bkg"+util::toTString(bkgs_.size()));
    util::HistOps::shiftStartOfXAxisToFirstPopulatedBin(h);
    h->UseCurrentStyle();
    h->SetTitle("");
    h->SetFillColor(bkgFillColor(bkg));
    h->SetLineColor(h->GetFillColor());
    h->SetNdivisions(505);
    if( rebin_ > 1 ) h->Rebin(rebin_);
    h->GetXaxis()->SetRange(1,h->FindBin(xMax_));
    if( scale != 1. ) h->Scale(scale);

    bkgHists_[bkg] = h;
    bkgStack_->Add(h);

    if( totalBkg_ == 0 ) {
      totalBkg_ = static_cast<TH1*>(h->Clone(id_+"TotalBkg"));
      totalBkg_->Sumw2();
      totalBkg_->SetMarkerStyle(0);
      //totalBkg_->SetFillStyle(3004);
      totalBkg_->SetFillStyle(3354);
      totalBkg_->SetFillColor(kBlue+2);
    } else {
      totalBkg_->Add(h);
    }
  }


  void scaleTotalBkgTo(double targetTotalBkg) {
    std::cout << "  Scaling the total background prediction from " << totalBkg_->Integral() << std::flush;
    if( totalBkg_->Integral() ) totalBkg_->Scale(targetTotalBkg/totalBkg_->Integral());
    std::cout << " to " << totalBkg_->Integral() << std::endl;
  }


  void scaleTotalBkgStatUncertaintyTo(const std::vector<double> &regionStartPoints, const std::vector<double> &bkgN, const std::vector<double> &targetStatUncert) {
    std::cout << "  Scaling the total statistical uncertainty of the background" << std::endl;
    for(int bin = 1; bin <= totalBkg_->GetNbinsX(); ++bin) {
      unsigned int i = 0;
      while( i < regionStartPoints.size() && totalBkg_->GetBinCenter(bin) > regionStartPoints.at(i) ) ++i;
      --i;
      double scale = targetStatUncert.at(i)/sqrt(bkgN.at(i));
      totalBkg_->SetBinError(bin,scale*sqrt(totalBkg_->GetBinContent(bin)));
    }
  }


  void setTotalBkgSystUncertainty(const std::vector<double> &bkgN, const std::vector<double> &systUncert) { 
    relBkgSystUncert_ = systUncert;
    for(unsigned int i = 0; i < relBkgSystUncert_.size(); ++i) {
      relBkgSystUncert_.at(i) /= bkgN.at(i);
    }
  }


  void applyTotalBkgUncertainty(const std::vector<double> &regionStartPoints) {
    std::cout << "  Combining statistical and systematic uncertainty on background" << std::endl;
    for(int bin = 1; bin <= totalBkg_->GetNbinsX(); ++bin) {
      unsigned int i = 0;
      while( i < regionStartPoints.size() && totalBkg_->GetBinCenter(bin) > regionStartPoints.at(i) ) ++i;
      --i;
      totalBkg_->SetBinError(bin,sqrt(pow(totalBkg_->GetBinError(bin),2)+pow(totalBkg_->GetBinContent(bin)*relBkgSystUncert_.at(i),2)));
    }
  }



  TH1* data_;
  THStack* bkgStack_;

  
private:
  const TString id_;
  const int rebin_;
  const double xMax_;
  TH1* totalBkg_;
  std::vector<double> relBkgSystUncert_;
  std::vector<Bkg> bkgs_;
  std::map<Bkg,TH1*> bkgHists_;
};


void setBkgUncertainties(const TString &id, std::vector<double> &regionStartPoints, std::vector<double> &bkgN, std::vector<double> &bkgStatUncert, std::vector<double> &bkgSystUncert) {
  if( id == "HT" ) {
    regionStartPoints.push_back(0.);
  } else if( id == "MHT" ) {
    regionStartPoints.push_back(0.);
  } else {
    std::cerr << "ERROR: unknown id '" << id << "'" << std::endl;
    exit(1);
  }

//   // Baseline 1.14/fb 2011
//   const double bkgTotal       = 927.5; 
//   const double bkgTotalUncert = 103.1; 
//   const double bkgStatZToInv  = 12.3;
//   const double bkgStatLL      = 19.8;
//   const double bkgStatHadTau  = 8.0;
//   const double bkgStatQCD     = 35.2;

  // This is the total event yield from the proper combination
  bkgN.push_back(BKG_TOTAL);		
  // Approximated total stat uncertainty as quadratic sum of 
  // individual bkg stat uncertainties
  bkgStatUncert.push_back(sqrt( pow(BKG_ZTOINV_STAT,2)+pow(BKG_LL_STAT,2)+pow(BKG_HADTAU_STAT,2)+pow(BKG_QCD_STAT,2) ));
  // Approximated total systematic uncertainty such that when
  // adding the approximated total stat uncertainty the properly
  // combined total uncertainty is regained
  std::cout << ">>> BKG_TOTAL        = " << BKG_TOTAL << std::endl;
  std::cout << ">>> BKG_TOTAL_UNCERT = " << BKG_TOTAL_UNCERT << std::endl;
  bkgSystUncert.push_back(sqrt( pow(BKG_TOTAL_UNCERT,2) - pow(bkgStatUncert.back(),2) )); 
}
  


void plotSpectrum(const TString &id, const TString &xTitle, const TString &xUnit, int rebin, double xMax) {
  // Data
  HistContainer spectra(id,"RA2Hists/hist_DataBaseline_4p65ifb.root","h_"+id,xTitle,xUnit,rebin,xMax);

  // Add bkgs
  spectra.addBkg(QCD,"RA2Hists/QCD.root","QCD_"+id,1.19); // bias correction scale factor
  spectra.addBkg(LostLepton,"RA2Hists/LostLepton.root","LostLepton_"+id);
  spectra.addBkg(HadronicTau,"RA2Hists/HadTau.root","HadTau_"+id);
  spectra.addBkg(ZInv,"RA2Hists/hist_ZInvPredictedFromGJets_4p65ifb.root","h_RA2_"+id);

  // Scale uncertainties to baseline prediciton
  std::vector<double> regionStartPoints;
  std::vector<double> bkgN;
  std::vector<double> bkgStatUncert;
  std::vector<double> bkgSystUncert;
  setBkgUncertainties(id,regionStartPoints,bkgN,bkgStatUncert,bkgSystUncert);

  std::cout << "Start \tN\t\tstat \tsyst" << std::endl;
  for(unsigned int i = 0; i < regionStartPoints.size(); ++i) {
    std::cout << "  >" << regionStartPoints.at(i) << "\t" << std::flush;
    std::cout << bkgN.at(i) << "\t" << std::flush;
    std::cout << bkgStatUncert.at(i) << "\t" << std::flush;
    std::cout << bkgSystUncert.at(i) << std::endl;
  }
  spectra.scaleTotalBkgTo(bkgN.front());
  spectra.scaleTotalBkgStatUncertaintyTo(regionStartPoints,bkgN,bkgStatUncert);
  spectra.setTotalBkgSystUncertainty(bkgN,bkgSystUncert);
  spectra.applyTotalBkgUncertainty(regionStartPoints);


  // Do plots
//   TCanvas* can = new TCanvas("can"+id,id,500,500);
//   can->cd();
//   spectra.data_->Draw("PE1");
//   spectra.bkgStack_->Draw("HISTsame");
//   spectra.data_->Draw("PE1same");
//   spectra.legend()->Draw("same");
//   can->SetLogy();
//   gPad->RedrawAxis();
//   can->SaveAs("RA2DataVsEstimatedBkg_"+id+".eps","eps");
//   can->SaveAs("RA2DataVsEstimatedBkg_"+id+".root");

  TCanvas* canRatio = util::HistOps::createRatioTopCanvas();
  TPad *bPad = util::HistOps::createRatioBottomPad();
  TH1 *tFrame = util::HistOps::createRatioTopHist(spectra.data_);
  TH1 *bFrame = util::HistOps::createRatioBottomFrame(spectra.data_,spectra.data_->GetXaxis()->GetTitle(),"",0.01,1.99);
  bFrame->GetYaxis()->SetNdivisions(505);
  bFrame->GetYaxis()->SetTitle("Data / Bkg.");

  canRatio->cd();
  util::HistOps::setYRange(tFrame,1,3E-1);
  tFrame->Draw("PE1");
  spectra.bkgStack_->Draw("HISTsame");
  spectra.totalBkg()->Draw("E2same");
  tFrame->Draw("PE1same");

  // Label
  TPaveText* line = util::LabelFactory::createPaveTextWithOffset(1,0.7,0.006);
  line->AddText("Bkg. predicted from data:");
  line->SetTextSize(0.045);
  line->Draw("same");
  spectra.legend()->Draw("same");
  TLegend* leg = util::LabelFactory::createLegendCol(1,-0.3);
  leg->AddEntry(spectra.data_,"Data","P");
  leg->SetTextSize(0.045);
  leg->Draw("same");

  canRatio->SetLogy();
  gPad->RedrawAxis();
  bPad->Draw();
  bPad->cd();
  bFrame->Draw();
  spectra.bkgForRatio()->Draw("E2same");
  spectra.dataForRatio()->Draw("PE1same");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".eps","eps");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".png");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".root");
}



void plotCombinedSpectra() {
  util::StyleSettings::setStylePAS();

  double val = BKG_ZTOINV;
  double unc = sqrt( pow(BKG_ZTOINV_STAT,2) + pow(0.5*(BKG_ZTOINV_SYSTUP+BKG_ZTOINV_SYSTDN),2) );
  std::cout << "ZTOINV\t = " << val << " +/- " << unc << " (" << unc/val << ")" << std::endl;

  val = BKG_LL;
  unc = sqrt( pow(BKG_LL_STAT,2) + pow(0.5*(BKG_LL_SYSTUP+BKG_LL_SYSTDN),2) );
  std::cout << "LL\t = " << val << " +/- " << unc << " (" << unc/val << ")" << std::endl;

  val = BKG_HADTAU;
  unc = sqrt( pow(BKG_HADTAU_STAT,2) + pow(0.5*(BKG_HADTAU_SYSTUP+BKG_HADTAU_SYSTDN),2) );
  std::cout << "HADTAU\t = " << val << " +/- " << unc << " (" << unc/val << ")" << std::endl;

  val = BKG_QCD;
  unc = sqrt( pow(BKG_QCD_STAT,2) + pow(0.5*(BKG_QCD_SYSTUP+BKG_QCD_SYSTDN),2) );
  std::cout << "QCD\t = " << val << " +/- " << unc << " (" << unc/val << ")" << std::endl;

  val = BKG_TOTAL;
  unc = BKG_TOTAL_UNCERT;
  std::cout << "TOTAL\t = " << val << " +/- " << unc << " (" << unc/val << ")" << std::endl;

  plotSpectrum("HT","H_{T}","GeV",2,2300);
  plotSpectrum("MHT","#slash{H}_{T}","GeV",5,950);
}
