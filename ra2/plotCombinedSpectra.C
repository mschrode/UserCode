// $Id: plotCombinedSpectra.C,v 1.2 2011/08/13 21:31:53 mschrode Exp $

// Compare HT and MHT spectra in data with bkg. prediction

#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
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


enum Bkg { ZInv, LostLepton, HadronicTau, QCD };

class HistContainer {
public:
  HistContainer(const TString &id, const TString &dataFileName, const TString &dataHistName, const TString &xTitle, const TString &xUnit, int rebin, double xMax)
    : id_(id), rebin_(rebin), xMax_(xMax) {
    
    std::cout << "Constructing spectra '" << id_ << "'" << std::endl;

    data_ = util::FileOps::readTH1(dataFileName,dataHistName,id_+"data");
    util::HistOps::shiftStartOfXAxisToFirstPopulatedBin(data_);
    data_->UseCurrentStyle();
    data_->SetTitle(util::StyleSettings::title(1140,true));
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
    else if( bkg == HadronicTau ) label = "W/t#bar{t}#rightarrowe/#tau_{h}+X";
    else if( bkg == QCD ) label = "QCD";

    return label;
  }


  int bkgFillColor(Bkg bkg) const {
    int color = kBlack;
    if( bkg == ZInv ) color = kGreen+1;
    else if( bkg == LostLepton ) color = kRed+3;
    else if( bkg == HadronicTau ) color = kRed+1;
    else if( bkg == QCD ) color = kYellow;

    return color;
  }


  TLegend* legend() const {
    TLegend* leg = util::LabelFactory::createLegendCol(bkgs_.size()+1,0.6);
    leg->AddEntry(data_,"Data ("+util::toTString(data_->Integral(),0)+")","P");
    for(std::vector<Bkg>::const_reverse_iterator rit = bkgs_.rbegin();
	rit != bkgs_.rend(); ++rit) {
      const TH1* h = bkgHists_.find(*rit)->second;
      TString numEvts = util::toTString(h->Integral(1,10000),1);	// Include overflow bins
      if( *rit == LostLepton ) numEvts = "243.5";
      else if( *rit == HadronicTau ) numEvts = "263.0";
      leg->AddEntry(h,bkgLabel(*rit)+" ("+numEvts+")","F");
    }

    return leg;
  }


  TH1* totalBkg() const {
    return static_cast<TH1*>(totalBkg_->Clone(id_+"totalBkgInclUncerts"));
  }


  TH1* bkgForRatio() const {
    TH1* h = static_cast<TH1*>(totalBkg_->Clone(id_+"totalBkgInRatioInclUncerts"));
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      double relUncert = h->GetBinError(bin);
      if( h->GetBinContent(bin) > 0. ) relUncert /= h->GetBinContent(bin);
      h->SetBinError(bin,relUncert);
      h->SetBinContent(bin,1.);
    }

    return h;
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
      totalBkg_->SetFillStyle(3004);
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


  void scaleTotalBkgStatUncertaintyTo(double targetStatUncert) {
    double totUncert = sqrt(totalBkg_->Integral());
    double scale = targetStatUncert / totUncert;
    std::cout << "  Scaling the total statistical uncertainty of the background from " << totUncert << " to " << targetStatUncert << std::endl;
    for(int bin = 1; bin <= totalBkg_->GetNbinsX(); ++bin) {
      //      std::cout << "Bin " << bin << ":\tN = " << totalBkg_->GetBinContent(bin) << "\t-->\tsqrt(N) = " << sqrt(totalBkg_->GetBinContent(bin)) << " (" << totalBkg_->GetBinError(bin) << ")" << std::endl;
      totalBkg_->SetBinError(bin,scale*sqrt(totalBkg_->GetBinContent(bin)));
    }
  }


  void setTotalBkgSystUncertainty(double systUncert) { 
    relBkgSystUncert_ = systUncert;
    if( totalBkg_->Integral() > 0. ) relBkgSystUncert_ /= totalBkg_->Integral();
  }


  void applyTotalBkgUncertainty() {
    std::cout << "  Combining statistical and systematic uncertainty on background" << std::endl;
    for(int bin = 1; bin <= totalBkg_->GetNbinsX(); ++bin) {
      totalBkg_->SetBinError(bin,sqrt(pow(totalBkg_->GetBinError(bin),2)+pow(totalBkg_->GetBinContent(bin)*relBkgSystUncert_,2)));
    }
  }




  TH1* data_;
  THStack* bkgStack_;

  
private:
  const TString id_;
  const int rebin_;
  const double xMax_;
  TH1* totalBkg_;
  double relBkgSystUncert_;
  std::vector<Bkg> bkgs_;
  std::map<Bkg,TH1*> bkgHists_;
};


void plotSpectrum(const TString &id, const TString &xTitle, const TString &xUnit, int rebin, double xMax) {
  // Data
  HistContainer spectra(id,"RA2Hists/RA2_BaseLineHists_LPData.root","h_"+id,xTitle,xUnit,rebin,xMax);

  // Add bkgs
  spectra.addBkg(QCD,"RA2Hists/QCD.root","QCD_"+id,1.082*0.997); // Scale factor: lumi, bias correction
  spectra.addBkg(LostLepton,"RA2Hists/LP_All_Final_Data_Residual_App.root","RA2Hists/LP_All_Final_Data_Residual_App.root:///FinalPlotMuCS/"+id+id,1.14/1.0458);	// Scale factor: lumi
  spectra.addBkg(HadronicTau,"RA2Hists/HadronicTau.root","baseline"+id);
  spectra.addBkg(ZInv,"RA2Hists/RA2_ZInvWithPhoton_BaseLineHists_LPData.root","h_RA2_"+id);

  // Scale uncertainties to baseline prediciton
  double bkgN = 927.5;		// Total bkg prediction
  double bkgTotalUncert = 103.1; // Total uncertainty on bkg prediction
  double bkgStatUncert = sqrt( pow(12.3,2)+pow(19.8,2)+pow(8.,2)+pow(35.2,2) );	// Quadratic sum of individual stat uncerts.
  double bkgSystUncert = sqrt( pow(bkgTotalUncert,2) - pow(bkgStatUncert,2) ); // Total uncert \ominus total stat uncert
  spectra.scaleTotalBkgTo(bkgN);
  spectra.scaleTotalBkgStatUncertaintyTo(bkgStatUncert);
  spectra.setTotalBkgSystUncertainty(bkgSystUncert);
  spectra.applyTotalBkgUncertainty();


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

  canRatio->cd();
  util::HistOps::setYRange(tFrame,1,3E-1);
  tFrame->Draw("PE1");
  spectra.bkgStack_->Draw("HISTsame");
  spectra.totalBkg()->Draw("E3same");
  tFrame->Draw("PE1same");
  spectra.legend()->Draw("same");
  canRatio->SetLogy();
  gPad->RedrawAxis();
  bPad->Draw();
  bPad->cd();
  bFrame->Draw();
  spectra.bkgForRatio()->Draw("E3same");
  spectra.dataForRatio()->Draw("PE1same");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".eps","eps");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".root");
}



void plotCombinedSpectra() {
  util::StyleSettings::setStylePAS();

  plotSpectrum("HT","H_{T}","GeV",2,1900);
  plotSpectrum("MHT","#slash{H}_{T}","GeV",5,950);
}
