// $Id: plotCombinedSpectra.C,v 1.13 2012/03/21 20:20:19 mschrode Exp $

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
#include "TVirtualPad.h"

#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"






// --------------------------------------------------------------
// --------------------------------------------------------------
enum Bkg { ZInv, HadronicTau, LostLepton, QCD }; // Determines plotting order
enum UncertaintyCombinationMethod { Stat, PredSystStat };



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
      std::cout << "     Bkg '" << bit->first << "':" << std::endl;
      BkgMap::const_iterator intvIt = intervalMax_.find(bit->first);
      BkgMap::const_iterator predIt = centralPred_.find(bit->first);
      BkgMap::const_iterator systIt = systUncerts_.find(bit->first);
      BkgMap::const_iterator statIt = statUncerts_.find(bit->first);
      if( intvIt == intervalMax_.end() || systIt == systUncerts_.end() || statIt == statUncerts_.end() ) {
	std::cerr << "ERROR: missing target uncertainties for bkg '" << bit->first << "'" << std::endl;
	exit(1);
      } else {
	int binMin = 1;
	int binMax = binMin;
	for(size_t intervalIdx = 0; intervalIdx < intvIt->second.size(); ++intervalIdx) {
	  // Scale bin content
	  double totP = predIt->second.at(intervalIdx);
	  double scalePred = 0.;
	  for(int bin = binMin; bin <= bit->second->GetNbinsX(); ++bin,++binMax) {
	    if( bit->second->GetXaxis()->GetBinUpEdge(bin) > intvIt->second.at(intervalIdx) ) {
	      break;
	    } else {
	      scalePred += bit->second->GetBinContent(bin);
	    }
	  }
	  --binMax;
	  scalePred = totP/scalePred;
	  for(int bin = binMin; bin <= binMax; ++bin) {
	    double n = bit->second->GetBinContent(bin);
	    bit->second->SetBinContent(bin,scalePred*n);
	  }	  

	  // Scale bin uncertainties
	  double totS = sqrt( pow(statIt->second.at(intervalIdx),2.) + pow(systIt->second.at(intervalIdx),2.) );
	  double scaleUncert = 0.;
	  for(int bin = binMin; bin <= binMax; ++bin) {
	    double n = bit->second->GetBinContent(bin);
	    scaleUncert += n + n*n;
	  }	
	  scaleUncert = totS / sqrt(scaleUncert);
	  for(int bin = binMin; bin <= binMax; ++bin) {
	    double n = bit->second->GetBinContent(bin);
	    double e = scaleUncert*sqrt(n + n*n);
	    bit->second->SetBinError(bin,e);
	  }	  
	  binMin = binMax+1;
	  binMax = binMin;

	  std::cout << "       " << bit->second->GetXaxis()->GetBinLowEdge(binMin) << " - " << bit->second->GetXaxis()->GetBinUpEdge(binMax) << " (bins " << binMin << " - " << binMax << ")" << std::endl;
	  std::cout << "         target central value = " << totP << ",  scale = " << scalePred << std::endl;
	  std::cout << "         target uncertainty   = " << totS << ",  scale = " << scaleUncert << std::endl;
	} // End of loop over intervals
      } // End if all bkgs are defined
    } // End of loop over all bkgs

    for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      double totEntries = 0;
      double totUncert = 0.;
      for(std::map<Bkg,TH1*>::iterator bit = bkgs.begin(); bit != bkgs.end(); ++bit) {
	totEntries += bit->second->GetBinContent(bin);
	totUncert += pow(bit->second->GetBinError(bin),2);
      }
      totUncert = sqrt(totUncert);
      h->SetBinContent(bin,totEntries);
      h->SetBinError(bin,totUncert);
    }
  }


  // --------------------------------------------------------------
  void scaleUncertainties(TH1* &h) const {
    std::cout << "\n  ------ Scaling uncertainties ------------------------" << std::endl;
    int binMin = 1;
    int binMax = binMin;
    for(size_t intervalIdx = 0; intervalIdx < intervalMax_.begin()->second.size(); ++intervalIdx) {
      double sumErr = 0.;
      for(int bin = binMin; bin <= h->GetNbinsX(); ++bin,++binMax) {
	if( h->GetXaxis()->GetBinUpEdge(bin) > intervalMax_.begin()->second.at(intervalIdx) ) {
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
	  //	  std::cout << "       > " << bin << ": e = " << scale*sqrt(h->GetBinContent(bin)) << std::endl;
	  h->SetBinError(bin,scale*sqrt(h->GetBinContent(bin)));
	}
      }
      binMin = binMax+1;
      binMax = binMin;
    }
  }


  // --------------------------------------------------------------
  void addTargetUncertainties(double intervalMax, Bkg bkg, double pred, double syst, double stat) {
    if( method_ != PredSystStat ) {
      std::cerr << "ERROR: incompatible uncertainty combination method" << std::endl;
      exit(1);
    }

    if( intervalMax_.find(bkg) == intervalMax_.end() ) {
      intervalMax_[bkg] = std::vector<double>(1,intervalMax);
      centralPred_[bkg] = std::vector<double>(1,pred);
      statUncerts_[bkg] = std::vector<double>(1,stat);
      systUncerts_[bkg] = std::vector<double>(1,syst);
    } else {
      intervalMax_.find(bkg)->second.push_back(intervalMax);
      centralPred_.find(bkg)->second.push_back(pred);
      statUncerts_.find(bkg)->second.push_back(stat);
      systUncerts_.find(bkg)->second.push_back(syst);
    }
  }


  // --------------------------------------------------------------
  void addTargetUncertaintyInInterval(double intervalMax, double targetUncert) {
    if( method_ != Stat ) {
      std::cerr << "ERROR: incompatible uncertainty combination method" << std::endl;
      exit(1);
    }
    if( intervalMax_.size() == 0 ) {
      intervalMax_[QCD] = std::vector<double>(1,intervalMax);
    } else {
      intervalMax_.begin()->second.push_back(intervalMax);
    }
    targetUncert_.push_back(targetUncert);
  }


  // --------------------------------------------------------------
  void printUncertainties(const TH1* h) const {
    std::cout << "\n  ------ Scaled uncertainties -------------------------" << std::endl;
    int binMin = 1;
    int binMax = binMin;
    for(size_t intervalIdx = 0; intervalIdx < intervalMax_.begin()->second.size(); ++intervalIdx) {
      double sumVal = 0.;
      double sumErr = 0.;
      for(int bin = binMin; bin <= h->GetNbinsX(); ++bin,++binMax) {
	if( h->GetXaxis()->GetBinUpEdge(bin) > intervalMax_.begin()->second.at(intervalIdx) ) {
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
      std::cout << "     Bkg '" << bit->first << "':" << std::endl;

      BkgMap::const_iterator intvIt = intervalMax_.find(bit->first);
      BkgMap::const_iterator predIt = centralPred_.find(bit->first);
      BkgMap::const_iterator systIt = systUncerts_.find(bit->first);
      BkgMap::const_iterator statIt = statUncerts_.find(bit->first);

      int binMin = 1;
      int binMax = binMin;
      double totEntriesAllInt = 0.;
      double totErrAllInt = 0.;
      for(size_t intervalIdx = 0; intervalIdx < intvIt->second.size(); ++intervalIdx) {
	double totEntries = 0.;
	double totErr = 0.;
	for(int bin = binMin; bin <= bit->second->GetNbinsX(); ++bin,++binMax) {
	  if( bit->second->GetXaxis()->GetBinUpEdge(bin) > intvIt->second.at(intervalIdx) ) {
	    break;
	  } else {
	    totEntries += bit->second->GetBinContent(bin);
	    totErr += pow(bit->second->GetBinError(bin),2.);
	  }	
	}
	--binMax;
	totEntriesAllInt += totEntries;
	totErrAllInt += totErr;

	double pred = predIt->second.at(intervalIdx);
	double syst = systIt->second.at(intervalIdx);
	double stat = statIt->second.at(intervalIdx);
	std::cout << "       " << h->GetXaxis()->GetBinLowEdge(binMin) << " - " << h->GetXaxis()->GetBinUpEdge(binMax) << " (bins " << binMin << " - " << binMax << ")" << std::endl;
	std::cout << "          target central value = " << pred << ",  val = " << totEntries << std::endl;
	std::cout << "          target uncertainty   = " << sqrt( pow(stat,2.) + pow(syst,2.) ) << " ( = sqrt( " << syst << "^2 + " << stat << "^2 ) ),  err = " << sqrt(totErr) << std::endl;
	binMin = binMax+1;
	binMax = binMin;
      } // End of loop over intervals
      std::cout << "       Total bkg '" << bit->first << "': " << totEntriesAllInt << " +/- " << sqrt(totErrAllInt) << std::endl;
    } // End of loop over bkgs
  }  

private:
  typedef std::map< Bkg, std::vector<double> > BkgMap;

  const UncertaintyCombinationMethod method_;

  BkgMap intervalMax_;  
  BkgMap centralPred_;
  BkgMap statUncerts_;
  BkgMap systUncerts_;  
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
    data_->SetTitle(util::StyleSettings::title(4980,true));
    data_->SetMarkerStyle(20);
    data_->SetNdivisions(505);
    if( rebin_ > 1 ) data_->Rebin(rebin_);
    data_->GetXaxis()->SetRange(1,data_->FindBin(xMax_));
    util::HistOps::setAxisTitles(data_,xTitle,xUnit,"events");
    data_->GetXaxis()->SetTitle(xTitle+" ["+xUnit+"]");
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
    if( um->method() == Stat ) {
      um->scaleUncertainties(totalBkg_);
      um->printUncertainties(totalBkg_);
    } else if( um->method() == PredSystStat ) {
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
void plotSpectrum(const TString &id, const UncertaintyManager* um, const TString &xTitle, const TString &xUnit, int rebin, double xMax, bool recreateOutputFile = false) {
  // Data
  HistContainer spectra(id,"RA2Hists/hist_DataBaseline_4p65ifb.root","h_"+id,xTitle,xUnit,rebin,xMax);

  // Add bkgs
  spectra.addBkg(QCD,"RA2Hists/QCD.root","QCD_"+id);
  spectra.addBkg(LostLepton,"RA2Hists/LostLepton.root","LostLepton_"+id,1.04);
  spectra.addBkg(HadronicTau,"RA2Hists/HadTau.root","HadTau_"+id);
  spectra.addBkg(ZInv,"RA2Hists/hist_ZInvPredictedFromGJets_4p65ifb.root","h_RA2_"+id);
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

  const TString outputName = "RA2DataVsEstimatedBkg.root";
  TFile* f = 0;
  if( recreateOutputFile ) f = new TFile(outputName,"recreate");
  else f = new TFile(outputName,"update");
  f->WriteTObject(hData,id+"Data");
  f->WriteTObject(sBkg,id+"BkgStack");
  f->WriteTObject(hSignal,id+"SignalExpectation");
  f->WriteTObject(hDataRatio,id+"DataRatio");
  f->WriteTObject(gBkgUncert,id+"BkgUncert");
  f->WriteTObject(gBkgRatio,id+"BkgUncertRatio");
  f->WriteTObject(legData,id+"SignalLegend");
  f->WriteTObject(legBkg,id+"BkgLegend");
  f->Close();
  delete f;
}



void plotCombinedSpectra() {
  util::StyleSettings::setStylePAS();

  // HT spectrum
  UncertaintyManager umHT(PredSystStat);
  umHT.addTargetUncertainties(	800	,	ZInv	,	494.1	,	114.7	,	13.3311664905964	);
  umHT.addTargetUncertainties(	1000	,	ZInv	,	74.8	,	31.5	,	5.09607692249636	);
  umHT.addTargetUncertainties(	1200	,	ZInv	,	18.7	,	11	,	2.6134268690744	);
  umHT.addTargetUncertainties(	1400	,	ZInv	,	5.4	,	3.8	,	1.34536240470737	);
  umHT.addTargetUncertainties(	100000	,	ZInv	,	3.2	,	2.2	,	0.7	);
										
  umHT.addTargetUncertainties(	800	,	LostLepton	,	376.8	,	48.5	,	22.8354986807821	);
  umHT.addTargetUncertainties(	1000	,	LostLepton	,	65.4	,	8.55	,	13.6436798555229	);
  umHT.addTargetUncertainties(	1200	,	LostLepton	,	19.2	,	2.7	,	5.37680202350803	);
  umHT.addTargetUncertainties(	1400	,	LostLepton	,	6.2	,	0.8	,	2.30217288664427	);
  umHT.addTargetUncertainties(	100000	,	LostLepton	,	2.6	,	0.4	,	1.5	);
										
  umHT.addTargetUncertainties(	800	,	HadronicTau	,	422.16	,	35.391	,	15.4300226830682	);
  umHT.addTargetUncertainties(	1000	,	HadronicTau	,	65.72	,	6.134	,	5.55847829895917	);
  umHT.addTargetUncertainties(	1200	,	HadronicTau	,	27.11	,	2.8035	,	4.01393821576765	);
  umHT.addTargetUncertainties(	1400	,	HadronicTau	,	6.769	,	0.8705	,	1.70789461033168	);
  umHT.addTargetUncertainties(	100000	,	HadronicTau	,	1.08	,	0.154	,	0.481	);
										
  umHT.addTargetUncertainties(	800	,	QCD	,	120.12	,	77.52	,	12.1070062360602	);
  umHT.addTargetUncertainties(	1000	,	QCD	,	36.24	,	24.3	,	5.41653948568641	);
  umHT.addTargetUncertainties(	1200	,	QCD	,	20.18	,	12.82	,	4.44658295773283	);
  umHT.addTargetUncertainties(	1400	,	QCD	,	11.84	,	7.725	,	3.4410608829255	);
  umHT.addTargetUncertainties(	100000	,	QCD	,	11.9	,	8.25	,	3.8	);

  plotSpectrum("HT",&umHT,"H_{T}","GeV",2,2300,true);
  

  // MHT spectrum
  UncertaintyManager umMHT(PredSystStat);
  umMHT.addTargetUncertainties(	350	,	ZInv	,	424.7	,	108.8	,	11.8570654042221	);
  umMHT.addTargetUncertainties(	500	,	ZInv	,	135.8	,	38.2	,	7.60394634384015	);
  umMHT.addTargetUncertainties(	600	,	ZInv	,	26.9	,	11.7	,	3.38378486313773	);
  umMHT.addTargetUncertainties(	100000	,	ZInv	,	8.8	,	4.5	,	1.72046505340853	);
										
										
  umMHT.addTargetUncertainties(	350	,	LostLepton	,	400.5	,	50.35	,	25.9984614929422	);
  umMHT.addTargetUncertainties(	500	,	LostLepton	,	59.5	,	8.6	,	7.75370878999205	);
  umMHT.addTargetUncertainties(	600	,	LostLepton	,	8.7	,	1.6	,	2.62678510731274	);
  umMHT.addTargetUncertainties(	100000	,	LostLepton	,	1.5	,	0.4	,	1.06301458127347	);
										
										
  umMHT.addTargetUncertainties(	350	,	HadronicTau	,	434.45	,	36.185	,	15.8426153459585	);
  umMHT.addTargetUncertainties(	500	,	HadronicTau	,	73.139	,	6.892	,	5.36030820009447	);
  umMHT.addTargetUncertainties(	600	,	HadronicTau	,	12.25	,	1.6095	,	2.36450438781576	);
  umMHT.addTargetUncertainties(	100000	,	HadronicTau	,	3	,	0.6665	,	1.71172427686237	);
										
										
  umMHT.addTargetUncertainties(	350	,	QCD	,	196.2	,	127.8	,	14.7461859475595	);
  umMHT.addTargetUncertainties(	500	,	QCD	,	3.98	,	2.705	,	2.10309296038002	);
  umMHT.addTargetUncertainties(	600	,	QCD	,	0.09	,	0.095	,	0.29748949561287	);
  umMHT.addTargetUncertainties(	100000	,	QCD	,	0.01	,	0.015	,	0.1	);
  
  plotSpectrum("MHT",&umMHT,"#slash{H}_{T}","GeV",5,950);


}



// Both plots on one canvas, common y axis
void plotCombinedSpectra(const TString &inputPlots) {
  util::StyleSettings::setStylePAS();

  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetTitleSize(0.055,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadTopMargin(0.0);
  gStyle->SetPadBottomMargin(0.0);
  gStyle->SetPadLeftMargin(0.0);
  gStyle->SetPadRightMargin(0.0);

  const double t_ = 0.02;
  const double b_ = 0.15;
  const double y_ = (1.-t_-b_);
  const double l_ = 0.07;
  const double r_ = 0.01;
  const double x_ = (1.-l_-r_)/2.;
  const double fracHeightTopPad = 0.8;
  const double ratioMin = 0.1;
  const double ratioMax = 2.4;

  TCanvas* canComb = new TCanvas("canComb","HT + MHT",500*y_/x_,500);

  TPad* HTPad = new TPad("HTPad","",0.,0.,l_+x_,1.);
  HTPad->SetTopMargin(t_);
  HTPad->SetLeftMargin(2.*l_);
  HTPad->SetBottomMargin(0.2 + 0.8*b_-0.2*t_);
  HTPad->SetRightMargin(0.);
  HTPad->SetFillStyle(0);
  HTPad->SetFrameFillColor(10);
  HTPad->SetFrameBorderMode(0);
  HTPad->SetRightMargin(0.);

  TPad* HTRatioPad = static_cast<TPad*>(HTPad->Clone("HTRatioPad"));
  HTRatioPad->SetTopMargin(t_);
  HTRatioPad->SetBottomMargin(b_);
  HTRatioPad->SetTopMargin(0.8 - 0.8*b_+0.2*t_);

  TH1* HTFrame = new TH1D("HTFrame",";;Events",19,500,2400);
  //TH1* HTFrame = new TH1D("HTFrame",";;Events",1899,501,2399);
  HTFrame->GetYaxis()->SetRangeUser(3E-1,3E4);
  HTFrame->GetXaxis()->SetLabelSize(0);
  HTFrame->SetNdivisions(505);
  HTFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/fracHeightTopPad);
  HTFrame->SetLineWidth(2);

  TH1* HTRatioFrame = static_cast<TH1*>(HTFrame->Clone("HTRatioFrame"));
  HTRatioFrame->GetXaxis()->SetTitle("H_{T} [GeV]");
  HTRatioFrame->GetYaxis()->SetTitle("Data / Pred.");
  HTRatioFrame->GetYaxis()->SetRangeUser(ratioMin,ratioMax);
  HTRatioFrame->GetYaxis()->SetNdivisions(404);
  HTRatioFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/(1.-fracHeightTopPad));
  HTRatioFrame->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
  HTRatioFrame->SetLineStyle(2);
  for(int bin = 1; bin <= HTRatioFrame->GetNbinsX(); ++bin) {
    HTRatioFrame->SetBinContent(bin,1.);
    HTRatioFrame->SetBinError(bin,0.);
  }


  TPad* MHTPad = new TPad("MHTPad","",l_+x_,0.,1,1.);
  MHTPad->SetLeftMargin(0.);
  MHTPad->SetBottomMargin(0.2 + 0.8*b_-0.2*t_);
  MHTPad->SetTopMargin(t_);
  MHTPad->SetRightMargin(2.*r_);
  MHTPad->SetFillStyle(0);
  MHTPad->SetFrameFillColor(10);
  MHTPad->SetFrameBorderMode(0);

  TPad* MHTRatioPad = static_cast<TPad*>(MHTPad->Clone("MHTRatioPad"));
  MHTRatioPad->SetBottomMargin(b_);
  MHTRatioPad->SetTopMargin(0.8 - 0.8*b_+0.2*t_);

  TH1* MHTFrame = new TH1D("MHTFrame","",16,200,1000);
  //TH1* MHTFrame = new TH1D("MHTFrame","",1598,200.5,999.5);
  MHTFrame->GetYaxis()->SetRangeUser(3E-1,3E4);
  MHTFrame->GetXaxis()->SetLabelSize(0);
  MHTFrame->GetYaxis()->SetLabelSize(0);
  MHTFrame->SetNdivisions(505);
  MHTFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/fracHeightTopPad);
  MHTFrame->SetLineWidth(2);

  TH1* MHTRatioFrame = static_cast<TH1*>(MHTFrame->Clone("MHTRatioFrame"));
  MHTRatioFrame->GetXaxis()->SetTitle("#slash{H}_{T} [GeV]");
  MHTRatioFrame->GetYaxis()->SetLabelSize(0);
  MHTRatioFrame->GetYaxis()->SetRangeUser(ratioMin,ratioMax);
  MHTRatioFrame->GetYaxis()->SetNdivisions(404);
  MHTRatioFrame->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/(1.-fracHeightTopPad));
  MHTRatioFrame->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
  MHTRatioFrame->SetLineStyle(2);
  for(int bin = 1; bin <= MHTRatioFrame->GetNbinsX(); ++bin) {
    MHTRatioFrame->SetBinContent(bin,1.);
    MHTRatioFrame->SetBinError(bin,0.);
  }
  double scaleLabels = (l_+x_)/(x_+r_);
  //MHTRatioFrame->GetXaxis()->SetTitleOffset(HTRatioFrame->GetXaxis()->GetTitleOffset()*scaleLabels);
  MHTRatioFrame->GetXaxis()->SetTitleOffset(1.);
  MHTRatioFrame->GetXaxis()->SetTitleSize(HTRatioFrame->GetXaxis()->GetTitleSize()*scaleLabels);
  //MHTRatioFrame->GetXaxis()->SetLabelOffset(HTRatioFrame->GetXaxis()->GetLabelOffset()*scaleLabels);
  MHTRatioFrame->GetXaxis()->SetLabelOffset(0.004);
  MHTRatioFrame->GetXaxis()->SetLabelSize(HTRatioFrame->GetXaxis()->GetLabelSize()*scaleLabels);

  // Get plots from file
  TH1* HTData = util::FileOps::readTH1(inputPlots,"HTData","",false);
  THStack* HTBkg = util::FileOps::readTHStack(inputPlots,"HTBkgStack");
  TH1* HTSignal = util::FileOps::readTH1(inputPlots,"HTSignalExpectation","",false);
  TH1* HTDataRatio = util::FileOps::readTH1(inputPlots,"HTDataRatio","",false);
  TGraphAsymmErrors* HTBkgUncert = util::FileOps::readTGraphAsymmErrors(inputPlots,"HTBkgUncert");
  TGraphAsymmErrors* HTBkgUncertRatio = util::FileOps::readTGraphAsymmErrors(inputPlots,"HTBkgUncertRatio");
  TH1* MHTData = util::FileOps::readTH1(inputPlots,"MHTData","",false);
  THStack* MHTBkg = util::FileOps::readTHStack(inputPlots,"MHTBkgStack");
  TH1* MHTSignal = util::FileOps::readTH1(inputPlots,"MHTSignalExpectation","",false);
  TH1* MHTDataRatio = util::FileOps::readTH1(inputPlots,"MHTDataRatio","",false);
  TGraphAsymmErrors* MHTBkgUncert = util::FileOps::readTGraphAsymmErrors(inputPlots,"MHTBkgUncert");
  TGraphAsymmErrors* MHTBkgUncertRatio = util::FileOps::readTGraphAsymmErrors(inputPlots,"MHTBkgUncertRatio");

  // Labels
  TLegend* legSignal = 0;
  TLegend* legBkg = 0;

  TFile f(inputPlots,"READ");
  f.GetObject("HTSignalLegend",legSignal);
  f.GetObject("HTBkgLegend",legBkg);
  f.Close();

  legSignal->SetTextSize(0.06);
  legSignal->SetX1NDC(0.05);
  legSignal->SetX2NDC(0.3);
  legSignal->SetY1NDC(0.78);
  legSignal->SetY2NDC(0.95);

  legBkg->SetX1NDC(0.50);
  legBkg->SetX2NDC(0.95);
  legBkg->SetY1NDC(0.65);
  legBkg->SetY2NDC(0.93);

  TPaveText* cmsLabel = new TPaveText(0.2,0.88,0.93,0.95,"NDC");
  cmsLabel->SetBorderSize(0);
  cmsLabel->SetFillColor(0);
  cmsLabel->SetTextFont(42);
  cmsLabel->SetTextSize(0.06);
  cmsLabel->AddText("CMS,  L = 4.98 fb^{-1},  #sqrt{s} = 7 TeV");

  TPaveText* cover1 = new TPaveText(0.92,0.09,1.,0.145,"NDC");
  cover1->SetBorderSize(0);
  cover1->SetFillColor(0);
  TPaveText* cover2 = new TPaveText(0.,0.09,0.1,0.145,"NDC");
  cover2->SetBorderSize(0);
  cover2->SetFillColor(0);
  TPaveText* cover3 = new TPaveText(0.1,0.09,0.2,0.145,"NDC");
  cover3->SetBorderSize(0);
  cover3->SetFillColor(0);

  canComb->cd();
  HTPad->Draw();
  HTPad->cd();
  HTFrame->Draw();
  HTBkg->Draw("HISTsame");
  HTBkgUncert->Draw("E2same");
  HTSignal->Draw("HISTsame");
  HTData->Draw("PE1same");
  cmsLabel->Draw("same");
  HTPad->SetLogy();
  HTPad->RedrawAxis();
  canComb->cd();
  HTRatioPad->Draw();
  HTRatioPad->cd();
  HTRatioFrame->Draw();
  HTBkgUncertRatio->Draw("E2same");
  HTDataRatio->Draw("PE1same");
  cover3->Draw("same");

  canComb->cd();
  MHTPad->Draw();
  MHTPad->cd();
  MHTFrame->Draw();
  MHTBkg->Draw("HISTsame");
  MHTBkgUncert->Draw("E2same");
  MHTSignal->Draw("HISTsame");
  MHTData->Draw("PE1same");
  legSignal->Draw("same");
  legBkg->Draw("same");
  MHTPad->SetLogy();
  MHTPad->RedrawAxis();
  canComb->cd();
  MHTRatioPad->Draw();
  MHTRatioPad->cd();
  MHTRatioFrame->Draw();
  MHTBkgUncertRatio->Draw("E2same");
  MHTDataRatio->Draw("PE1same");
  cover1->Draw("same");
  cover2->Draw("same");

  canComb->SaveAs("RA2DataVsEstimatedBkg.eps","eps");
  canComb->SaveAs("RA2DataVsEstimatedBkg.png");
}


