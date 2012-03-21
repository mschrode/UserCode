// $Id: plotCombinedSpectra.C,v 1.12 2012/03/18 20:16:47 mschrode Exp $

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

//   // HT spectrum
//   UncertaintyManager um(InclPredSystStat);
//    um.addTargetUncertainties(ZInv,81.,16.);
//    um.addTargetUncertainties(LostLepton,0.5*(55.2+53.7),22.4);
//    um.addTargetUncertainties(HadronicTau,0.5*(38.2+35.6),16.);
//    um.addTargetUncertainties(QCD,0.5*(188.03+186.43),21.38);
//   plotSpectrum("HT",&um,"H_{T}","GeV",2,2300);

//   // MHT spectrum
//   plotSpectrum("MHT",&um,"#slash{H}_{T}","GeV",5,950);


//   // HT spectrum
//   UncertaintyManager umHT(Stat);
//   umHT.addTargetUncertaintyInInterval(800.,128.98);
//   umHT.addTargetUncertaintyInInterval(1000.,36.29);
//   umHT.addTargetUncertaintyInInterval(1200.,16.88);
//   umHT.addTargetUncertaintyInInterval(1400.,9.36);
//   umHT.addTargetUncertaintyInInterval(100000.,9.6);
//   plotSpectrum("HT",&umHT,"H_{T}","GeV",2,2300);
  
//   // MHT spectrum
//   UncertaintyManager umMHT(Stat);
//   umMHT.addTargetUncertaintyInInterval(350.,131.67);
//   umMHT.addTargetUncertaintyInInterval(500.,31.43);
//   umMHT.addTargetUncertaintyInInterval(600.,8.65);
//   umMHT.addTargetUncertaintyInInterval(100000.,4.22);
//   plotSpectrum("MHT",&umMHT,"#slash{H}_{T}","GeV",5,950);



  // HT spectrum
  UncertaintyManager umHT(PredSystStat);

  umHT.addTargetUncertainties(	800	,       ZInv	,	494.1	,	85.2601313627888	,	13.3311664905964	);
  umHT.addTargetUncertainties(	1000	,	ZInv	,	74.8	,	20.42914584607	,	5.09607692249636	);
  umHT.addTargetUncertainties(	1200	,	ZInv	,	18.7	,	6.63626400921482	,	2.6134268690744	);
  umHT.addTargetUncertainties(	1400	,	ZInv	,	5.4	,	2.70185121722126	,	1.34536240470737	);
  umHT.addTargetUncertainties(	100000	,	ZInv	,	3.2	,	2.2	,	0.7	);
										
  umHT.addTargetUncertainties(	800	,	LostLepton	,      376.8	,	40.8800684930933	,	22.8354986807821	);
  umHT.addTargetUncertainties(	1000	,	LostLepton	,	65.4	,	7.15891053163818	,	13.6436798555229	);
  umHT.addTargetUncertainties(	1200	,	LostLepton	,	19.2	,	1.86279360101972	,	5.37680202350803	);
  umHT.addTargetUncertainties(	1400	,	LostLepton	,	6.2	,	0.58309518948453	,	2.30217288664427	);
  umHT.addTargetUncertainties(	100000	,	LostLepton	,	2.6	,	0.4	,	1.5	);
										
  umHT.addTargetUncertainties(	800	,	HadronicTau	,	422.16	,	27.954195856794	,	15.4300226830682	);
  umHT.addTargetUncertainties(	1000	,	HadronicTau	,	65.72	,	5.00833345535219	,	5.55847829895917	);
  umHT.addTargetUncertainties(	1200	,	HadronicTau	,	27.11	,	2.11494633501656	,	4.01393821576765	);
  umHT.addTargetUncertainties(	1400	,	HadronicTau	,	6.769	,	0.740121611628792	,	1.70789461033168	);
  umHT.addTargetUncertainties(	100000	,	HadronicTau	,	1.08	,	0.159	,	0.481	);
										
  umHT.addTargetUncertainties(	800	,	QCD	,	120.12	,	76.00947572507	,	12.1070062360602	);
  umHT.addTargetUncertainties(	1000	,	QCD	,	36.24	,	24.0133566999701	,	5.41653948568641	);
  umHT.addTargetUncertainties(	1200	,	QCD	,	20.18	,	12.0029204779504	,	4.44658295773283	);
  umHT.addTargetUncertainties(	1400	,	QCD	,	11.84	,	7.60263112349929	,	3.4410608829255	);
  umHT.addTargetUncertainties(	100000	,	QCD	,	11.9	,	8.1	,	3.8	);

  plotSpectrum("HT",&umHT,"H_{T}","GeV",2,2300);
  

  // MHT spectrum
  UncertaintyManager umMHT(PredSystStat);
  umMHT.addTargetUncertainties(	350	,	ZInv	,	424.7	,	83.3394264439107	,	11.8570654042221	);
  umMHT.addTargetUncertainties(	500	,	ZInv	,	135.8	,	27.170940359141	,	7.60394634384015	);
  umMHT.addTargetUncertainties(	600	,	ZInv	,	26.9	,	6.96921803361037	,	3.38378486313773	);
  umMHT.addTargetUncertainties(	100000	,	ZInv	,	8.8	,	3.24499614791759	,	1.72046505340853	);
  
  
  umMHT.addTargetUncertainties(	350	,	LostLepton	,	400.5	,	40.9609570200697	,	25.9984614929422	);
  umMHT.addTargetUncertainties(	500	,	LostLepton	,	59.5	,	6.88912185985993	,	7.75370878999205	);
  umMHT.addTargetUncertainties(	600	,	LostLepton	,	8.7	,	1.0295630140987	,	2.62678510731274	);
  umMHT.addTargetUncertainties(	100000	,	LostLepton	,	1.5	,	0.282842712474619	,	1.06301458127347	);
  
  
  umMHT.addTargetUncertainties(	350	,	HadronicTau	,	434.45	,	27.9305276892507	,	15.8426153459585	);
  umMHT.addTargetUncertainties(	500	,	HadronicTau	,	73.139	,	5.47683932574254	,	5.36030820009447	);
  umMHT.addTargetUncertainties(	600	,	HadronicTau	,	12.25	,	1.09719551584939	,	2.36450438781576	);
  umMHT.addTargetUncertainties(	100000	,	HadronicTau	,	3	,	0.502538555734781	,	1.71172427686237	);
										
  
  umMHT.addTargetUncertainties(	350	,	QCD	,	196.2	,	81.3595108146552	,	14.7461859475595	);
  umMHT.addTargetUncertainties(	500	,	QCD	,	3.98	,	1.47905375155875	,	2.10309296038002	);
  umMHT.addTargetUncertainties(	600	,	QCD	,	0.09	,	0.0616441400296898	,	0.29748949561287	);
  umMHT.addTargetUncertainties(	100000	,	QCD	,	0.01	,	0.02	,	0.1	);
  
  plotSpectrum("MHT",&umMHT,"#slash{H}_{T}","GeV",5,950);


}
