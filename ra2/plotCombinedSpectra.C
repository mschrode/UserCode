// $Id: plotCombinedSpectra.C,v 1.11 2012/03/15 23:08:34 mschrode Exp $

// Compare HT and MHT spectra in data with bkg. prediction

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






// --------------------------------------------------------------
// --------------------------------------------------------------
enum Bkg { ZInv, HadronicTau, LostLepton, QCD }; // Determines plotting order
enum UncertaintyCombinationMethod { BinsStat, InclSystStat };



// --------------------------------------------------------------
class UncertaintyManager {
public:
  // --------------------------------------------------------------
  UncertaintyManager(UncertaintyCombinationMethod combMethod)
    : method_(combMethod) {}

  
  // --------------------------------------------------------------
  UncertaintyCombinationMethod method() const { return method_; }


  // --------------------------------------------------------------
  void scaleUncertainties(std::map<Bkg,TH1*> &bkgs, TH1* &h) const {
    std::cout << "\n  ------ Scaling uncertainties ------------------------" << std::endl;

    for(std::map<Bkg,TH1*>::iterator bit = bkgs.begin(); bit != bkgs.end(); ++bit) {
      std::map<Bkg,double>::const_iterator systIt = systUncerts_.find(bit->first);
      std::map<Bkg,double>::const_iterator statIt = statUncerts_.find(bit->first);
      if( systIt == systUncerts_.end() || statIt == statUncerts_.end() ) {
	std::cerr << "ERROR: missing target uncertainties for bkg '" << bit->first << "'" << std::endl;
	exit(1);
      } else {
	// Total target uncertainty 
	double totS = sqrt( pow(statIt->second,2.) + pow(systIt->second,2.) );
	// Scale factor for bin uncertainties
	double scale = 0.;
	for(int bin = 1; bin <= bit->second->GetNbinsX(); ++bin) {
	  double n = bit->second->GetBinContent(bin);
	  scale += n + n*n;
	}
	scale = totS / sqrt(scale);
	// Scale bin uncertainties
	for(int bin = 1; bin <= bit->second->GetNbinsX(); ++bin) {
	  double n = bit->second->GetBinContent(bin);
	  double e = scale*sqrt(n + n*n);
	  bit->second->SetBinError(bin,e);
	}	  
      }
    }

    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      double totUncert = 0.;
      for(std::map<Bkg,TH1*>::iterator bit = bkgs.begin(); bit != bkgs.end(); ++bit) {
	totUncert += pow(bit->second->GetBinError(bin),2);
      }
      totUncert = sqrt(totUncert);
      h->SetBinError(bin,totUncert);
    }
  }


  // --------------------------------------------------------------
  void scaleUncertainties(TH1* &h) const {
    std::cout << "\n  ------ Scaling uncertainties ------------------------" << std::endl;
    int binMin = 1;
    int binMax = binMin;
    for(size_t intervalIdx = 0; intervalIdx < intervalMax_.size(); ++intervalIdx) {
      double sumErr = 0.;
      for(int bin = binMin; bin <= h->GetNbinsX(); ++bin,++binMax) {
	if( h->GetXaxis()->GetBinUpEdge(bin) > intervalMax_.at(intervalIdx) ) {
	  break;
	} else {
	  sumErr += h->GetBinContent(bin);
	}	
      }
      --binMax;
      sumErr = sqrt(sumErr);
      if( sumErr > 0. ) {
	double scale = targetUncert_.at(intervalIdx)/sumErr;

	std::cout << "     " << h->GetXaxis()->GetBinLowEdge(binMin) << " - " << h->GetXaxis()->GetBinUpEdge(binMax) << " (bins " << binMin << " - " << binMax << "):  target uncertainty = " << targetUncert_.at(intervalIdx) << ",  scale = " << scale << ",  sumErr = " << sumErr << std::endl;

	for(int bin = binMin; bin <= binMax; ++bin) {
	  std::cout << "       > " << bin << ": e = " << scale*sqrt(h->GetBinContent(bin)) << std::endl;
	  h->SetBinError(bin,scale*sqrt(h->GetBinContent(bin)));
	}
      }
      binMin = binMax+1;
      binMax = binMin;
    }
  }


  // --------------------------------------------------------------
  void addTargetUncertainties(Bkg bkg, double syst, double stat) {
    if( statUncerts_.find(bkg) != statUncerts_.end() ) {
      std::cerr << "WARNING: uncertainties already added for bkg '" << bkg << "': Ignoring input" << std::endl;
    } else {
      statUncerts_[bkg] = stat;
      systUncerts_[bkg] = syst;
    }
  }


  // --------------------------------------------------------------
  void addTargetUncertaintyInInterval(double intervalMax, double targetUncert) {
    if( method_ != BinsStat ) {
      std::cerr << "ERROR: incompatible uncertainty combination method" << std::endl;
      exit(1);
    }
    if( intervalMax_.size() > 0 ) {
      if( intervalMax < intervalMax_.back() ) {
	std::cerr << "ERROR in UncertaintyManager::addTargetUncertaintyInInterval(): max=" << intervalMax << " > previous threshold (" << intervalMax_.back() << ")" << std::endl;
      }
    }
    intervalMax_.push_back(intervalMax);
    targetUncert_.push_back(targetUncert);
  }


  // --------------------------------------------------------------
  void printUncertainties(const TH1* h) const {
    std::cout << "\n  ------ Scaled uncertainties -------------------------" << std::endl;
    int binMin = 1;
    int binMax = binMin;
    for(size_t intervalIdx = 0; intervalIdx < intervalMax_.size(); ++intervalIdx) {
      double sumVal = 0.;
      double sumErr = 0.;
      for(int bin = binMin; bin <= h->GetNbinsX(); ++bin,++binMax) {
	if( h->GetXaxis()->GetBinUpEdge(bin) > intervalMax_.at(intervalIdx) ) {
	  break;
	} else {
	  sumVal += h->GetBinContent(bin);
	  sumErr += pow(h->GetBinError(bin),2.);
	}	
      }
      --binMax;
      std::cout << "     " << h->GetXaxis()->GetBinLowEdge(binMin) << " - " << h->GetXaxis()->GetBinUpEdge(binMax) << " (bins " << binMin << " - " << binMax << "):  target uncertainty = " << targetUncert_.at(intervalIdx) << ",  yield = " << sumVal << ",  err = " << sqrt(sumErr) << std::endl;
      binMin = binMax+1;
      binMax = binMin;
    }
  }


  // --------------------------------------------------------------
  void printUncertainties(const std::map<Bkg,TH1*> &bkgs, const TH1* h) const {
    std::cout << "\n  ------ Scaled uncertainties -------------------------" << std::endl;

    for(std::map<Bkg,TH1*>::const_iterator bit = bkgs.begin(); bit != bkgs.end(); ++bit) {
      double totErr = 0.;
      for(int bin = 1; bin <= bit->second->GetNbinsX(); ++bin) {
	totErr += pow(bit->second->GetBinError(bin),2.);
      }
      double syst = systUncerts_.find(bit->first)->second;
      double stat = statUncerts_.find(bit->first)->second;
      std::cout << "     Bkg '" << bit->first << "': target uncertainty = " << sqrt( pow(stat,2.) + pow(syst,2.) ) << " ( = sqrt( " << syst << "^2 + " << stat << "^2 ) ),  err = " << sqrt(totErr) << std::endl;
    }
    double totErr = 0.;
    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      totErr += pow(h->GetBinError(bin),2.);
    }
    std::cout << "     err = " << sqrt(totErr) << std::endl;
  }
  

private:
  const UncertaintyCombinationMethod method_;

  std::map<Bkg,double> statUncerts_;
  std::map<Bkg,double> systUncerts_;  
  std::vector<double> intervalMax_;
  std::vector<double> targetUncert_;
};



// --------------------------------------------------------------
class HistContainer {
public:
  // --------------------------------------------------------------
  HistContainer(const TString &id, const TString &dataFileName, const TString &dataHistName, const TString &xTitle, const TString &xUnit, int rebin, double xMax)
    : id_(id), rebin_(rebin), xMax_(xMax) {
    
    std::cout << "\n\n*****  Constructing " << id << " spectra ********************" << std::endl;

    data_ = util::FileOps::readTH1(dataFileName,dataHistName,id_+"data");
    util::HistOps::shiftStartOfXAxisToFirstPopulatedBin(data_);
    data_->UseCurrentStyle();
    data_->SetTitle(util::StyleSettings::title(4700,true)); // Lumi rounded to one digit
    data_->SetMarkerStyle(20);
    data_->SetNdivisions(505);
    if( rebin_ > 1 ) data_->Rebin(rebin_);
    data_->GetXaxis()->SetRange(1,data_->FindBin(xMax_));
    util::HistOps::setAxisTitles(data_,xTitle,xUnit,"events");
    util::HistOps::setYRange(data_,1,3E-1);

    totalBkg_ = 0;
    signal_ = 0;
  }


  // --------------------------------------------------------------
  TString bkgLabel(Bkg bkg) const {
    TString label = "";
    if( bkg == ZInv ) label = "Z#rightarrow#nu#bar{#nu}+jets";
    else if( bkg == LostLepton ) label = "W/t#bar{t}#rightarrowe/#mu+X";
    else if( bkg == HadronicTau ) label = "W/t#bar{t}#rightarrow#tau_{h}+X";
    else if( bkg == QCD ) label = "QCD";

    return label;
  }


  // --------------------------------------------------------------
  int bkgFillColor(Bkg bkg) const {
    int color = kBlack;
    if( bkg == ZInv ) color = kGreen+1;
    else if( bkg == HadronicTau ) color = kYellow;
    else if( bkg == LostLepton ) color = kRed+1;
    else if( bkg == QCD ) color = kRed+3;

    return color;
  }


  // --------------------------------------------------------------
  TLegend* legend() const {
    TLegend* leg = util::LabelFactory::createLegendColWithOffset(bkgs_.size(),0.6,0.06);
    for(BkgCIt it = bkgs_.begin(); it != bkgs_.end(); ++it) {
      leg->AddEntry(it->second,bkgLabel(it->first),"F");
    }

    return leg;
  }


  // --------------------------------------------------------------
  THStack* bkgStack() const {
    THStack* h = new THStack(id_+"BkgStack"+util::toTString(++COUNT),"");
    for(BkgCRIt rit = bkgs_.rbegin(); rit != bkgs_.rend(); ++rit) {
      TString name = rit->second->GetName();
      h->Add(static_cast<TH1*>(rit->second->Clone(name+"InStack"+util::toTString(COUNT))));
    }

    return h;
  }

  // --------------------------------------------------------------
  TGraphAsymmErrors* bkgUncertaintyBand() const {
    return util::HistOps::getUncertaintyBand(totalBkg_,kBlue+2,3004);
  }

  // --------------------------------------------------------------
  TGraphAsymmErrors* bkgRatio() const {
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


  // --------------------------------------------------------------
  TH1* data() const {
    return static_cast<TH1*>(data_->Clone(id_+"Data"+util::toTString(++COUNT)));
  }

  // --------------------------------------------------------------
  TH1* dataRatio() const {
    TH1* h = static_cast<TH1*>(data_->Clone(id_+"DataRatio"+util::toTString(++COUNT)));
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

  // --------------------------------------------------------------
  TH1* signalExpectation() const {
    return static_cast<TH1*>(signal_->Clone(id_+"SignalExpectation"+util::toTString(++COUNT)));
  }


  // --------------------------------------------------------------
  void print() const {
    std::cout << "\n  ------ Event yields ---------------------------------" << std::endl;
    std::cout << "    Data: " << data_->Integral(1,10000) << std::endl;
    for(BkgCIt it = bkgs_.begin(); it != bkgs_.end(); ++it) {
      std::cout << "     " << bkgLabel(it->first) << ": " << it->second->Integral(1,10000) << std::endl;
    }
  }


  // --------------------------------------------------------------
  void addBkg(Bkg bkg, const TString &fileName, const TString &histName, double scale = 1.) {
    TH1* h = util::FileOps::readTH1(fileName,histName,id_+"Bkg"+util::toTString(bkgs_.size()));
    util::HistOps::shiftStartOfXAxisToFirstPopulatedBin(h);
    h->UseCurrentStyle();
    h->SetTitle("");
    h->SetFillColor(bkgFillColor(bkg));
    h->SetLineColor(h->GetFillColor());
    h->SetNdivisions(505);
    if( rebin_ > 1 && bkg != QCD ) h->Rebin(rebin_);
    h->GetXaxis()->SetRange(1,h->FindBin(xMax_));
    if( scale != 1. ) h->Scale(scale);
    
    if( totalBkg_ == 0 ) {
      totalBkg_ = static_cast<TH1*>(h->Clone(id_+"TotalBkg"));
      //totalBkg_->Sumw2();
      totalBkg_->SetMarkerStyle(0);
      //totalBkg_->SetFillStyle(3004);
      totalBkg_->SetFillStyle(3354);
      totalBkg_->SetFillColor(kBlue+2);
    } else {
      totalBkg_->Add(h);
    }

    bkgs_[bkg] = h;
  }


  // --------------------------------------------------------------
  void scaleBkgUncertainty(const UncertaintyManager* um) {
    if( um->method() == BinsStat ) {
      um->scaleUncertainties(totalBkg_);
      um->printUncertainties(totalBkg_);
    } else if( um->method() == InclSystStat ) {
      um->scaleUncertainties(bkgs_,totalBkg_);
      um->printUncertainties(bkgs_,totalBkg_);
    }
  }


  // --------------------------------------------------------------
  void addSignalExpectation(const TString &fileName, const TString &histName) {
    if( signal_ ) delete signal_;
    signal_ = util::FileOps::readTH1(fileName,histName,id_+"SignalExpectation");
    util::HistOps::shiftStartOfXAxisToFirstPopulatedBin(signal_);
    signal_->UseCurrentStyle();
    signal_->SetTitle("");
    signal_->SetLineColor(600);
    signal_->SetNdivisions(505);
    if( rebin_ > 1 ) signal_->Rebin(rebin_);
    signal_->GetXaxis()->SetRange(1,signal_->FindBin(xMax_));
  }

  
private:
  typedef std::map<Bkg,TH1*> Bkgs;
  typedef std::map<Bkg,TH1*>::const_iterator BkgCIt;
  typedef std::map<Bkg,TH1*>::const_reverse_iterator BkgCRIt;

  mutable unsigned int COUNT;

  const TString id_;
  const int rebin_;
  const double xMax_;

  TH1* data_;
  TH1* totalBkg_;
  Bkgs bkgs_;
  TH1* signal_;
};
  

// --------------------------------------------------------------
void plotSpectrum(const TString &id, const UncertaintyManager* um, const TString &xTitle, const TString &xUnit, int rebin, double xMax) {
  // Data
  HistContainer spectra(id,"RA2Hists/hist_DataBaseline_4p65ifb.root","h_"+id,xTitle,xUnit,rebin,xMax);

  // Add bkgs
  spectra.addBkg(QCD,"RA2Hists/QCD.root","QCD_"+id);
  spectra.addBkg(LostLepton,"RA2Hists/LostLepton.root","LostLepton_"+id,1.04);
  spectra.addBkg(HadronicTau,"RA2Hists/HadTau.root","HadTau_"+id);
  spectra.addBkg(ZInv,"RA2Hists/Paper2011_PreApproval/hist_ZInvPredictedFromGJets_4p65ifb.root","h_RA2_"+id);
  spectra.print();

  // Scale uncertainties
  spectra.scaleBkgUncertainty(um);

  // Add signal point
  spectra.addSignalExpectation("LM5BaselineHisto_4p65_NLO.root","h_"+id);

  // Plots
  TH1* hData = spectra.data();
  TH1* hDataRatio = spectra.dataRatio();
  THStack* sBkg = spectra.bkgStack();
  TGraphAsymmErrors* gBkgUncert = spectra.bkgUncertaintyBand();
  TGraphAsymmErrors* gBkgRatio = spectra.bkgRatio();
  TH1* hSignal = spectra.signalExpectation();

  // Labels
  TPaveText* line = util::LabelFactory::createPaveTextWithOffset(1,0.7,0.006);
  line->AddText("Bkg. predicted from data:");
  line->SetTextSize(0.045);
  TLegend* legBkg = spectra.legend();
  TLegend* legData = util::LabelFactory::createLegendCol(2,-0.3);
  legData->AddEntry(hData,"Data","P");
  legData->AddEntry(hSignal,"LM5","L");
  legData->SetTextSize(0.045);

  TCanvas* canRatio = util::HistOps::createRatioTopCanvas();
  TPad *bPad = util::HistOps::createRatioBottomPad();
  TH1 *tFrame = util::HistOps::createRatioTopHist(hData);
  TH1 *bFrame = util::HistOps::createRatioBottomFrame(hData,0.01,2.49);
  bFrame->GetYaxis()->SetNdivisions(404);
  bFrame->GetYaxis()->SetTitle("Data / Bkg.");

  canRatio->cd();
  util::HistOps::setYRange(tFrame,1,3E-1);
  tFrame->Draw("PE1");
  sBkg->Draw("HISTsame");
  gBkgUncert->Draw("E2same");
  hSignal->Draw("HISTsame");
  tFrame->Draw("PE1same");
  legData->Draw("same");
  line->Draw("same");
  legBkg->Draw("same");
  canRatio->SetLogy();
  gPad->RedrawAxis();
  bPad->Draw();
  bPad->cd();
  bFrame->Draw();
  gBkgRatio->Draw("E2same");
  hDataRatio->Draw("PE1same");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".eps","eps");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".png");
  canRatio->SaveAs("RA2DataVsEstimatedBkg_"+id+".root");
}



void plotCombinedSpectra() {
  util::StyleSettings::setStylePAS();

  // HT spectrum
  UncertaintyManager um(InclSystStat);
   um.addTargetUncertainties(ZInv,81.,16.);
   um.addTargetUncertainties(LostLepton,0.5*(55.2+53.7),22.4);
   um.addTargetUncertainties(HadronicTau,0.5*(38.2+35.6),16.);
   um.addTargetUncertainties(QCD,0.5*(188.03+186.43),21.38);
  plotSpectrum("HT",&um,"H_{T}","GeV",2,2300);

  // MHT spectrum
  plotSpectrum("MHT",&um,"#slash{H}_{T}","GeV",5,950);





//   // HT spectrum
//   UncertaintyManager umHT(BinsStat);
//   umHT.addTargetUncertaintyInInterval(800.,128.98);
//   umHT.addTargetUncertaintyInInterval(1000.,36.29);
//   umHT.addTargetUncertaintyInInterval(1200.,16.88);
//   umHT.addTargetUncertaintyInInterval(1400.,9.36);
//   umHT.addTargetUncertaintyInInterval(100000.,9.6);
//   plotSpectrum("HT",&umHT,"H_{T}","GeV",2,2300);

//   // MHT spectrum
//   UncertaintyManager umMHT;
//   umMHT.addTargetUncertaintyInInterval(350.,131.67);
//   umMHT.addTargetUncertaintyInInterval(500.,31.43);
//   umMHT.addTargetUncertaintyInInterval(600.,8.65);
//   umMHT.addTargetUncertaintyInInterval(100000.,4.22);
//   plotSpectrum("MHT",&umMHT,"#slash{H}_{T}","GeV",5,950);
}
