// $Id: plotCombinedSpectra.C,v 1.1 2011/08/13 21:15:32 mschrode Exp $

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
    data_ = util::FileOps::readTH1(dataFileName,dataHistName,id_+"data");
    data_->UseCurrentStyle();
    data_->SetTitle(util::StyleSettings::title(1140));
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
    else if( bkg == QCD ) color = kYellow+1;

    return color;
  }


  TLegend* legend() const {
    TLegend* leg = util::LabelFactory::createLegendCol(bkgs_.size()+1,0.6);
    leg->AddEntry(data_,"Data ("+util::toTString(data_->Integral(),0)+")","P");
    for(std::vector<Bkg>::const_reverse_iterator rit = bkgs_.rbegin();
	rit != bkgs_.rend(); ++rit) {
      const TH1* h = bkgHists_.find(*rit)->second;
      TString numEvts = util::toTString(h->Integral(),1);
      if( *rit == LostLepton ) numEvts = "243.5";
      leg->AddEntry(h,bkgLabel(*rit)+" ("+numEvts+")","F");
    }

    return leg;
  }


  TH1* dataBkgRatio() const {
    return util::HistOps::createRatioPlot(data_,totalBkg_);
  }


  void addBkg(Bkg bkg, const TString &fileName, const TString &histName, double scale = 1.) {
    bkgs_.push_back(bkg);

    TH1* h = util::FileOps::readTH1(fileName,histName,id_+"Bkg"+util::toTString(bkgs_.size()));
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

    if( totalBkg_ == 0 ) totalBkg_ = static_cast<TH1*>(h->Clone(id_+"TotalBkg"));
    else totalBkg_->Add(h);
  }

  TH1* data_;
  THStack* bkgStack_;

  
private:
  const TString id_;
  const int rebin_;
  const double xMax_;
  TH1* totalBkg_;
  std::vector<Bkg> bkgs_;
  std::map<Bkg,TH1*> bkgHists_;
};


void plotSpectrum(const TString &id, const TString &xTitle, const TString &xUnit, int rebin, double xMax) {
  // Data
  HistContainer spectra(id,"RA2Hists/RA2_BaseLineHists_LPData.root","h_"+id,xTitle,xUnit,rebin,xMax);

  // Add bkgs
  //  spectra.addBkg(QCD,"RA2Hists/.root","");
  spectra.addBkg(LostLepton,"RA2Hists/LP_All_Final_Data_Residual_App.root","RA2Hists/LP_All_Final_Data_Residual_App.root:///FinalPlotMuCS/"+id+id,1.14/1.0458);
  spectra.addBkg(HadronicTau,"RA2Hists/HadronicTau.root","baseline"+id);
  spectra.addBkg(ZInv,"RA2Hists/RA2_ZInvWithPhoton_BaseLineHists_LPData.root","h_RA2_"+id);

TCanvas* can = new TCanvas("can"+id,id,500,500);
  can->cd();
  spectra.data_->Draw("PE1");
  spectra.bkgStack_->Draw("HISTsame");
  spectra.data_->Draw("PE1same");
  spectra.legend()->Draw("same");
  can->SetLogy();
  gPad->RedrawAxis();

  TCanvas* canRatio = util::HistOps::createRatioTopCanvas();
  TPad *bPad = util::HistOps::createRatioBottomPad();
  TH1 *tFrame = util::HistOps::createRatioTopHist(spectra.data_);
  TH1 *bFrame = util::HistOps::createRatioBottomFrame(spectra.data_,spectra.data_->GetXaxis()->GetTitle(),"",0.31,1.69);
  canRatio->cd();
  util::HistOps::setYRange(tFrame,1,3E-1);
  tFrame->Draw("PE1");
  spectra.bkgStack_->Draw("HISTsame");
  tFrame->Draw("PE1same");
  spectra.legend()->Draw("same");
  canRatio->SetLogy();
  gPad->RedrawAxis();
  bPad->Draw();
  bPad->cd();
  bFrame->Draw();
  spectra.dataBkgRatio()->Draw("PE1same");
}



void plotCombinedSpectra() {
  util::StyleSettings::setStylePAS();

  plotSpectrum("HT","H_{T}","GeV",1,1900);
  plotSpectrum("MHT","#slash{H}_{T}","GeV",2,950);
}
