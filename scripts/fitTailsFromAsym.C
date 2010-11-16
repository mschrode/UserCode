// $Id: fitTailsFromAsym.C,v 1.3 2010/11/10 21:04:29 mschrode Exp $

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"





TString gJetAlgo_ = "AK5 Calo-Jets";
TString gDeltaPhiCut_ = "|#Delta#phi| > 2.7";
TString gEtaMin_ = "0";
TString gEtaMax_ = "1.3";
TString gPPCut_ = "p_{||} < 0.1";
TString gUID_ = "uid";



TPaveText* createLabel(double ptAveMin, double ptAveMax, double lumi = -1.) {
  TPaveText *txt = util::LabelFactory::createPaveText(3,-0.68);
  if( lumi > 0. ) txt->AddText("L = "+util::toTString(lumi)+" pb^{-1},  "+gJetAlgo_+",");
  else txt->AddText(gJetAlgo_+",");
  txt->AddText(gDeltaPhiCut_+",  "+gPPCut_+"#upoint#bar{p^{ave}_{T}}");
  txt->AddText(util::toTString(ptAveMin)+" < p^{ave}_{T} < "+util::toTString(ptAveMax)+" GeV");
  return txt;
}

TPaveText* createLabelEta() {
  TPaveText *txt = util::LabelFactory::createPaveText(1,0.32);
  txt->AddText(gEtaMin_+" < |#eta| < "+gEtaMax_);
  return txt;
}

TLegend* createLegend(const TH1* hData, const TH1* hMC, const TF1* tailData = 0, const TF1* tailMC = 0) {
  TLegend *leg = 0;
  if( tailData && tailMC ) {
    leg = util::LabelFactory::createLegendColWithOffset(4,0.3,1);
    leg->AddEntry(hData,"Data","P");
    leg->AddEntry(tailData,"Fit Data","L");
    leg->AddEntry(hMC,"MC","F");
    leg->AddEntry(tailMC,"Fit MC","L");
  } else {
    leg = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
    leg->AddEntry(hData,"Data","P");
    leg->AddEntry(hMC,"MC","F");
  }
  return leg;
}



// --------------------------------------------------
class Sample {
public:
  //  enum SampleType { Data, MC };

  Sample(const TString &fileName, int type, const TString &sampleName, unsigned int startBin, unsigned int endBin, double lumi = -1.);
  ~Sample();

  TString name() const { return name_; }
  unsigned int nBins() const { return hPtAsym_.size(); }
  //  SampleType type() const { return type_; }
  int type() const { return type_; }
  
  TH1 *hPtAve() const { return hist(hPtAve_); }
  TH1 *hPtAsym(unsigned int i) const { return hist(i,hPtAsym_); }
  TH1 *hPtAsymSmeared(unsigned int i) const { return hist(i,hPtAsymSmeared_); }
  TH1 *hTailClean(unsigned int i) const { return hist(i,hTailsClean_); }
  TF1 *fTail(unsigned int i) const;
  
  double widthAsym(unsigned int i) const { return widthAsym_.at(i); }
  double widthAsymErr(unsigned int i) const { return widthAsymErr_.at(i); }  
  double widthAsymSmeared(unsigned int i) const { return widthAsymSmeared_.at(i); }
  double widthAsymErrSmeared(unsigned int i) const { return widthAsymErrSmeared_.at(i); }  

  double nTotal(unsigned int i) const { return nTotal_.at(i); }
  double nTail(unsigned int i) const { return nTail_.at(i); }

  void smearAsymmetry(const Sample* dataSample);


private:
  //const SampleType type_;
  const int type_;
  const TString name_;

  TH1 *hPtAve_;
  util::HistVec hPtJet_;
  util::HistVec hPtAsym_;
  util::HistVec hPtAsymSmeared_;

  util::HistVec hTails_;
  util::HistVec hTailsClean_;
  std::vector<TF1*> fTailFits_;

  std::vector<double> widthAsym_;
  std::vector<double> widthAsymErr_;
  std::vector<double> widthAsymSmeared_;
  std::vector<double> widthAsymErrSmeared_;

  std::vector<double> nTotal_;
  std::vector<double> nTail_;

  mutable unsigned int count_;

  void fitWidth(const util::HistVec &h, std::vector<double> &width, std::vector<double> &widthErr) const;
  double gaussInt(double mean, double sigma, double min, double max) const;
  bool getTails(const util::HistVec &h);
  bool getTail(const TH1* hAbsAsym, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) const;
  TH1* hist(const TH1* h) const;
  TH1* hist(unsigned int i, const util::HistVec &h) const;
  void setStyle(TH1 *h) const;
  void smearAsymmetry(const std::vector<double> &scaling, const util::HistVec &hOrig, util::HistVec &hSmeared, bool abs = false) const;
  void smearAsymmetry(double scaling, const util::HistVec &hOrig, util::HistVec &hSmeared, bool abs = false) const;
  void smearHistogram(const TH1* hOrig, TH1* &hSmeared, double nTotal, double width, double scaling, bool abs = false) const;
};


// --------------------------------------------------
Sample::Sample(const TString &fileName, int type, const TString &sampleName, unsigned int startBin, unsigned int endBin, double lumi)
  : type_(type), name_(sampleName+"_"+util::toTString(type)) {
  count_ = 0;

  // Read histograms from file
  bool ioIsOk = true;
  TFile file(fileName,"READ");
  TString newName;

  // Average pt
  hPtAve_ = 0;
  file.GetObject("hPtAveAbs",hPtAve_);
  if( hPtAve_ ) {
    hPtAve_->SetDirectory(0);
    newName = name()+"_PtAve";
    hPtAve_->SetName(newName);
    util::HistOps::setAxisTitles(hPtAve_,"p^{ave}_{T}","GeV","Events");
    hPtAve_->SetTitle("");
    if( lumi > 0. ) hPtAve_->Scale(lumi);
    setStyle(hPtAve_);
  } else {
    ioIsOk = false;
  }

  // Pt asymmetry
  bool binExists = true;
  unsigned int ptBin = startBin;
  while( ioIsOk && binExists && ptBin < endBin ) {
    TH1 *hPtAsym = 0;
    file.GetObject("hPtAsym_"+util::toTString(ptBin),hPtAsym);
    if( hPtAsym ) {
      hPtAsym->SetDirectory(0);
      newName = name()+"_PtAsym"+util::toTString(ptBin-startBin);
      hPtAsym->SetName(newName);
      hPtAsym->GetXaxis()->SetRangeUser(-0.5,0.5);
      util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","Events",true);      
      hPtAsym->SetTitle("");
      setStyle(hPtAsym);
      hPtAsym_.push_back(hPtAsym);
      nTotal_.push_back(hPtAsym->GetEntries());
      ++ptBin;
    } else {
      binExists = false;
    }    
  }

  file.Close();  

  if( !ioIsOk ) {
    std::cerr << "ERROR reading histograms from file '" << fileName << "'\n";
    exit(1);
  }


  // Fitting width of asymmetry
  fitWidth(hPtAsym_,widthAsym_,widthAsymErr_);

  if( type_ == 0 ) {
    // Extract tails from measured asymmetry
    getTails(hPtAsym_);
  } else {
    // Smearing asymmetry
    smearAsymmetry(0.1,hPtAsym_,hPtAsymSmeared_);
    fitWidth(hPtAsymSmeared_,widthAsymSmeared_,widthAsymErrSmeared_);

    // Extract tails from smeared asymmetry
    getTails(hPtAsymSmeared_);
  }
}



// --------------------------------------------------
Sample::~Sample() {
  //   delete hPtAve_;
  //   for(util::HistIt it = hPtJet_.begin(); it != hPtJet_.end(); ++it) {
  //     delete *it;
  //   }
  //   for(util::HistIt it = hPtAsym_.begin(); it != hPtAsym_.end(); ++it) {
  //     delete *it;
  //   }
  //   for(util::HistIt it = hPtAsymSmeared_.begin(); it != hPtAsymSmeared_.end(); ++it) {
  //     delete *it;
  //   }
  //   for(util::HistIt it = hTails_.begin(); it != hTails_.end(); ++it) {
  //     delete *it;
  //   }
  //   for(util::HistIt it = hTailsClean_.begin(); it != hTailsClean_.end(); ++it) {
  //     delete *it;
  //   }
  //   for(std::vector<TF1*>::iterator it = fTailFits_.begin(); it != fTailFits_.end(); ++it) {
  //     delete *it;
  //   }
}


// --------------------------------------------------
TH1* Sample::hist(const TH1* h) const {
  ++count_;
  TString name = h->GetName();
  name += count_;
  return static_cast<TH1*>(h->Clone(name));
}


// --------------------------------------------------
TH1* Sample::hist(unsigned int i, const util::HistVec &h) const {
  assert( i < nBins() );
  ++count_;
  TString name = h[i]->GetName();
  name += count_;
  return static_cast<TH1*>(h[i]->Clone(name));
}


// --------------------------------------------------
TF1* Sample::fTail(unsigned int i) const {
  assert( i < nBins() );
  ++count_;
  TString name = fTailFits_[i]->GetName();
  name += count_;
  return static_cast<TF1*>(fTailFits_[i]->Clone(name));
}


// --------------------------------------------------
void Sample::fitWidth(const util::HistVec &h, std::vector<double> &width, std::vector<double> &widthErr) const {
  width.clear();
  for(util::HistItConst it = h.begin(); it != h.end(); ++it) {
    double mean = (*it)->GetMean();
    double sig = 1.5*(*it)->GetRMS();
    if( (*it)->Fit("gaus","0QIR","",mean-sig,mean+sig) == 0 ) {
      mean = (*it)->GetFunction("gaus")->GetParameter(1);
      sig = 1.8*(*it)->GetFunction("gaus")->GetParameter(2);
      if( (*it)->Fit("gaus","0QIR","",mean-sig,mean+sig) == 0 ) {
	width.push_back((*it)->GetFunction("gaus")->GetParameter(2));
	widthErr.push_back((*it)->GetFunction("gaus")->GetParError(2));
      }
    } else {
      std::cerr << "WARNING in Sample::fitWidth (" << name() << "): No convergence when fitting width of '" << (*it)->GetName() << "'\n";
      width.push_back(0.);
      widthErr.push_back(100.);
    }
  }
}



// --------------------------------------------------
void Sample::smearAsymmetry(const Sample* dataSample) {
  // Width ratio data / MC
  std::vector<double> scales(dataSample->nBins());
  for(unsigned int i = 0; i < dataSample->nBins(); ++i) {
    scales[i] = ( dataSample->widthAsym(i)/widthAsym(i) - 1. );
  }

  // Smearing asymmetry
  smearAsymmetry(scales,hPtAsym_,hPtAsymSmeared_);
  fitWidth(hPtAsymSmeared_,widthAsymSmeared_,widthAsymErrSmeared_);
  
  // Extract tails from smeared asymmetry
  getTails(hPtAsymSmeared_);
}



// --------------------------------------------------
void Sample::smearAsymmetry(double scaling, const util::HistVec &hOrig, util::HistVec &hSmeared, bool abs) const {
  std::vector<double> scales(hOrig.size(),scaling);
  smearAsymmetry(scales,hOrig,hSmeared,abs);
}
 
 
// --------------------------------------------------
void Sample::smearAsymmetry(const std::vector<double> &scaling, const util::HistVec &hOrig, util::HistVec &hSmeared, bool abs) const {
  for(util::HistIt it = hSmeared.begin(); it != hSmeared.end(); ++it) {
    delete *it;
  }
  hSmeared.clear();
  hSmeared = util::HistVec(hOrig.size());

  for(unsigned int i = 0; i < hOrig.size(); ++i) {
    smearHistogram(hOrig[i],hSmeared[i],nTotal(i),widthAsym(i),scaling[i],abs);
  }
}


// --------------------------------------------------
void Sample::smearHistogram(const TH1* hOrig, TH1* &hSmeared, double nTotal, double width, double scaling, bool abs) const {
  TString name = hOrig->GetName();
  name += "Smeared";
  hSmeared = static_cast<TH1D*>(hOrig->Clone(name));
  hSmeared->Reset();
  setStyle(hSmeared);

  scaling += 1.;
  if( scaling > 1. ) {
    //    std::cout << "\n\nWidth = " << width << std::endl;
    scaling = sqrt( scaling*scaling - 1. )*width;
    //    std::cout << "Scaling: " << scaling << std::endl;
  } else {
    std::cerr << "WARNING in smearHistogram(): scaling = " << scaling << std::endl;
  }
  for(int bin = 1; bin <= hOrig->GetNbinsX(); ++bin) {
    double entries = hOrig->GetBinContent(bin);
    if( entries ) {
      double mean = hOrig->GetBinCenter(bin);
      double norm = 1.;
      if( abs ) norm = 0.5*(1.+erf(mean/sqrt(2.)/scaling));
      for(int i = 1; i <= hSmeared->GetNbinsX(); ++i) {
	double min = hSmeared->GetXaxis()->GetBinLowEdge(i);
	double max = hSmeared->GetXaxis()->GetBinUpEdge(i);
	double weight = gaussInt(mean,scaling,min,max)*entries/norm;
	hSmeared->Fill(hSmeared->GetBinCenter(i),weight);
      }
    }
  }
  for(int bin = 1; bin <= hSmeared->GetNbinsX(); ++bin) {
    hSmeared->SetBinError(bin,sqrt(hSmeared->GetBinContent(bin)/nTotal));
  }
}


// --------------------------------------------------
double Sample::gaussInt(double mean, double sigma, double min, double max) const {
  return 0.5*(erf((max-mean)/sqrt(2.)/sigma) - erf((min-mean)/sqrt(2.)/sigma));
}


// --------------------------------------------------
bool Sample::getTails(const util::HistVec &h) {
  for(util::HistIt it = hTails_.begin(); it != hTails_.end(); ++it) {
    delete *it;
  }
  hTails_.clear();
  hTails_ = util::HistVec(nBins());
  for(util::HistIt it = hTailsClean_.begin(); it != hTailsClean_.end(); ++it) {
    delete *it;
  }
  hTailsClean_.clear();
  hTailsClean_ = util::HistVec(nBins());
  for(std::vector<TF1*>::iterator it = fTailFits_.begin(); it != fTailFits_.end(); ++it) {
    delete *it;
  }
  fTailFits_.clear();
  fTailFits_ = std::vector<TF1*>(nBins());
  nTail_.clear();
  nTail_ = std::vector<double>(nBins());

  bool result = true;
  for(unsigned int i = 0; i < nBins(); ++i) {
    bool tmpResult = getTail(h[i],hTails_[i],hTailsClean_[i],fTailFits_[i]);
    result = result && tmpResult;
    
    nTail_[i] = hTailsClean_[i]->Integral("width");

    //    std::cout << "NTAIL (" << i << "): " << nTail(i) << std::endl; 
  }

  return result;
}


// --------------------------------------------------
bool Sample::getTail(const TH1* hAsym, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) const {
  bool success = false;
  
  TString name = hAsym->GetName();
  name += "_Tails";
  hTail = static_cast<TH1D*>(hAsym->Clone(name));
  hTail->SetLineColor(kBlue);
  hTail->SetMarkerColor(kBlue);
  
  name = hAsym->GetName();
  name += "_TailsClean";
  hTailClean = static_cast<TH1D*>(hAsym->Clone(name));
  hTailClean->Reset();
  
  name = hAsym->GetName();
  name += "_GaussFit";
  gauss = new TF1(name,"gaus",-1.,1.);
  gauss->SetLineWidth(1);
  gauss->SetLineColor(kRed);
  if( hAsym->GetFillColor() > 0 ) gauss->SetLineStyle(2);
  
  double sigma = 2.5*hTail->GetRMS();
  success =  !(hTail->Fit(gauss,"I0QR","",-sigma,sigma));
  if( success ) {
    sigma = 1.8*gauss->GetParameter(2);
    success = !(hTail->Fit(gauss,"I0QR","",-sigma,sigma));
    if( success ) {
      sigma = gauss->GetParameter(2);
      int binMin = 1;
      int binMax = hTail->FindBin(5.*sigma);
      for(int bin = binMin; bin <= binMax; ++bin) {
	double min = hTail->GetXaxis()->GetBinLowEdge(bin);
	double max = hTail->GetXaxis()->GetBinUpEdge(bin);
	double gaussPdf = gauss->Integral(min,max)/hTail->GetBinWidth(1);
	double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
	hTail->SetBinContent(bin,tailPdf);
      }
      for(int bin = 1; bin <= hTailClean->GetNbinsX(); ++bin) {
	if( hTail->GetBinContent(bin) > 0.4*gauss->Eval(hTail->GetBinCenter(bin)) ) {
	  hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	  hTailClean->SetBinError(bin,hTail->GetBinError(bin));
	}
      }
    }
  }
  
  return success;
}



// --------------------------------------------------
void Sample::setStyle(TH1 *h) const {
  //   if( type() == Data ) h->SetMarkerStyle(20);
  //   else if( type() == MC ) h->SetFillColor(5);
  
  h->SetLineWidth(1);
  if( type() == 0 ) {
    h->SetMarkerStyle(20);
  } else if( type() == 1 ) {
    h->SetMarkerStyle(1);
    h->SetFillColor(5);
  }
}





// --------------------------------------------------
class SampleAdmin {
public:
  SampleAdmin(const TString &fileData, const TString &fileMC, std::vector<double> ptBinEdges, double lumi, double ppCut, unsigned int mcStartBin = 0, unsigned int mcEndBin = 0);
  ~SampleAdmin();

  TString name() const { return name_; }
  TString outNamePrefix() const { return outNamePrefix_; }
  unsigned int nPtBins() const { return ptBinEdges_.size()-1; }
  double ptMin() const { return ptBinEdges_.front(); }
  double ptMax() const { return ptBinEdges_.back(); }
  double ptMin(unsigned int i) const { return ptBinEdges_.at(i); }
  double ptMax(unsigned int i) const { return ptBinEdges_.at(i+1); }
  std::vector<double> ptBinEdges() const { return ptBinEdges_; }

  double widthAsymData(unsigned int i) const { return data_->widthAsym(i); }
  double widthAsymErrData(unsigned int i) const { return data_->widthAsymErr(i); }  
  double widthAsymMC(unsigned int i) const { return mc_->widthAsym(i); }
  double widthAsymErrMC(unsigned int i) const { return mc_->widthAsymErr(i); }  
  double widthSmearedAsymMC(unsigned int i) const { return mc_->widthAsymSmeared(i); }
  double widthSmearedAsymErrMC(unsigned int i) const { return mc_->widthAsymErrSmeared(i); }  
  double nTailData(unsigned int i) const { return data_->nTail(i); }
  double nTailMC(unsigned int i) const { return mc_->nTail(i); }
  double nTotalData(unsigned int i) const { return data_->nTotal(i); }
  double nTotalMC(unsigned int i) const { return mc_->nTotal(i); }

  void plotSpectra() const;
  void plotAsymmetry() const;
  void plotSmearedAsymmetry() const;
  void plotTails() const;
  
private:
  const std::vector<double> ptBinEdges_;
  const double lumi_;
  const double ppCut_;

  Sample* data_;
  Sample* mc_;

  TString name_;
  TString outNamePrefix_;

  void plot(const TString &name, TH1 *hData, TH1 *hMC, bool logy = false) const;
};

typedef std::vector<SampleAdmin*>::const_iterator AdminIt;

//sampleName+"_"+util::toTString(ptMin())+"-"+util::toTString(ptMax())+"_"

SampleAdmin::SampleAdmin(const TString &fileData, const TString &fileMC, std::vector<double> ptBinEdges, double lumi, double ppCut, unsigned int mcStartBin, unsigned int mcEndBin)
  : ptBinEdges_(ptBinEdges), lumi_(lumi), ppCut_(ppCut) {

  name_= gUID_+"_Pt"+util::toTString(ptMin())+"-"+util::toTString(ptMax()); 
  outNamePrefix_ = name()+"_";
  std::cout << "Creating samples '" << outNamePrefix() << std::endl;;
  
  // Create samples holding histograms
  TString tmpName = name()+"_"+util::toTString(ptMin())+"-"+util::toTString(ptMax());
  data_ = new Sample(fileData,0,tmpName,0,1000);
  mc_ = new Sample(fileMC,1,tmpName,mcStartBin,mcEndBin,lumi_);
  mc_->smearAsymmetry(data_);

  // Checks
  assert( data_->nBins() == nPtBins() );
  assert( mc_->nBins() == nPtBins() );
  
  //   std::cout << "\n  WIDTH OF PT-ASYMMETRY\n";
  //   for(unsigned int i = 0; i < nPtBins(); ++i) {
  //     std::cout << "    " << i << std::flush;
  //     std::cout << " (" << ptMin(i) << " - " << ptMax(i) << "): " << std::flush;
  //     std::cout << widthAsymData(i) << " +/- " << widthAsymErrData(i) << " (Data), " << widthAsymMC(i) << " +/- " << widthAsymErrMC(i) << " (MC)" << std::endl;
  //   }  
}


SampleAdmin::~SampleAdmin() {
  delete data_;
  delete mc_;
}


void SampleAdmin::plotSpectra() const {
  TH1 *hMC = mc_->hPtAve();
  TH1 *hData = data_->hPtAve();
  util::HistOps::setYRange(hData,2,8E-2);

  TString canName = name()+"PtAve";
  TCanvas *can = new TCanvas(canName,canName,500,500);
  can->cd();
  hData->Draw("PE1");
  hMC->Draw("HISTEsame");
  hData->Draw("PE1same");
  createLabel(ptMin(),ptMax(),lumi_)->Draw("same");
  createLabelEta()->Draw("same");
  createLegend(hData,hMC)->Draw("same");
  can->SetLogy();
  can->SaveAs(outNamePrefix()+"PtAve.eps","eps");
}


void SampleAdmin::plotAsymmetry() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TH1 *hMC = mc_->hPtAsym(i);
    TH1 *hData = data_->hPtAsym(i);

    // Log scale
    util::HistOps::setYRange(hMC,3,3E-5);
    util::HistOps::setYRange(hData,3,3E-5);
    hMC->GetXaxis()->SetRangeUser(-1.,1.);

    TString canName = name()+"PtAsym_"+util::toTString(i);
    TCanvas *can = new TCanvas(canName,canName,500,500);
    can->cd();
    hMC->Draw("HISTE");
    hData->Draw("PE1same");
    createLabel(ptMin(i),ptMax(i),lumi_)->Draw("same");
    createLabelEta()->Draw("same");
    createLegend(hData,hMC)->Draw("same");
    can->SetLogy(1);
    can->SaveAs(outNamePrefix()+"PtAsymLog_"+util::toTString(i)+".eps","eps");

    // Linear scale
    util::HistOps::setYRange(hMC,3);
    hMC->GetXaxis()->SetRangeUser(-0.4,0.4);
    can->cd();
    hMC->Draw("HISTE");
    hData->Draw("PE1same");
    createLabel(ptMin(i),ptMax(i),lumi_)->Draw("same");
    createLabelEta()->Draw("same");
    createLegend(hData,hMC)->Draw("same");
    can->SetLogy(0);
    can->SaveAs(outNamePrefix()+"PtAsym_"+util::toTString(i)+".eps","eps");
  }
}


void SampleAdmin::plotSmearedAsymmetry() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TString canName = name()+"PtAsymSmear_"+util::toTString(i);
    TCanvas *can = new TCanvas(canName,canName,500,500);

    TH1 *hMC = mc_->hPtAsymSmeared(i);
    TH1 *hData = data_->hPtAsym(i);

    // Linear scale    
    util::HistOps::setYRange(hMC,3);
    hMC->GetXaxis()->SetRangeUser(-0.4,0.4);
    can->cd();
    hMC->Draw("HISTE");
    hData->Draw("PE1same");
    createLabel(ptMin(i),ptMax(i),lumi_)->Draw("same");
    createLabelEta()->Draw("same");
    createLegend(hData,hMC)->Draw("same");
    can->SetLogy(0);
    can->SaveAs(outNamePrefix()+"PtSmearAsym_"+util::toTString(i)+".eps","eps");

    // Log scale    
    util::HistOps::setYRange(hMC,3,3E-5);
    util::HistOps::setYRange(hData,3,3E-5);
    hMC->GetXaxis()->SetRangeUser(-1.,1.);
    can->cd();
    hMC->Draw("HISTE");
    hData->Draw("PE1same");
    createLabel(ptMin(i),ptMax(i),lumi_)->Draw("same");
    createLabelEta()->Draw("same");
    createLegend(hData,hMC)->Draw("same");
    can->SetLogy(1);
    can->SaveAs(outNamePrefix()+"PtSmearAsymLog_"+util::toTString(i)+".eps","eps");

    // Superimpose fit
    TF1 *tailMC = mc_->fTail(i);
    TF1 *tailData = data_->fTail(i);
    hMC->GetXaxis()->SetRangeUser(0.,1.);
    can->cd();
    hMC->Draw("HISTE");
    hData->Draw("PE1same");
    createLabel(ptMin(i),ptMax(i),lumi_)->Draw("same");
    createLabelEta()->Draw("same");
    createLegend(hData,hMC,tailData,tailMC)->Draw("same");
    tailMC->Draw("same");
    tailData->Draw("same");
    can->SetLogy(1);
    can->SaveAs(outNamePrefix()+"TailFit_"+util::toTString(i)+".eps","eps");
  }
}


void SampleAdmin::plotTails() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TH1 *hMC = mc_->hTailClean(i);
    TH1 *hData = data_->hTailClean(i);
    util::HistOps::setYRange(hMC,5);
    hMC->GetXaxis()->SetRangeUser(0.,1.);
    TString canName = name()+"Tail_"+util::toTString(i);
    TCanvas *can = new TCanvas(canName,canName,500,500);
    can->cd();
    hMC->Draw("HISTE");
    hData->Draw("PE1same");
    createLabel(ptMin(i),ptMax(i),lumi_)->Draw("same");
    createLabelEta()->Draw("same");
    createLegend(hData,hMC)->Draw("same");
    can->SaveAs(outNamePrefix()+"Tail_"+util::toTString(i)+".eps","eps");
  }
}


// void SampleAdmin::plot(const TString &name, TH1 *hData, TH1 *hMC, bool logy) const {
//   TCanvas *can = new TCanvas(name,name,500,500);
//   can->cd();
//   hMC->Draw("HISTE");
//   hData->Draw("PE1same");
//   if( logy ) can->SetLogy();
// }


void plotDataMCDiff(TH1* hData, TH1* hMC, TPaveText* label1, TPaveText* label2) {
  // Fit width
  double min = hMC->GetXaxis()->GetBinLowEdge(1);
  double max = hMC->GetXaxis()->GetBinUpEdge(hMC->GetNbinsX());
  TF1 *sigmaMC = new TF1("SigmaMC","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",min,max);
  sigmaMC->SetParameter(0,1.);
  sigmaMC->SetParameter(1,1.);
  sigmaMC->SetParameter(2,0.01);
  sigmaMC->SetLineWidth(1);
  sigmaMC->SetLineColor(2);
  sigmaMC->SetLineStyle(2);
  TF1 *sigmaData = static_cast<TF1*>(sigmaMC->Clone("SigmaData"));
  sigmaData->SetLineStyle(1);

  hMC->Fit(sigmaMC,"I0QR");
  hData->Fit(sigmaData,"I0QR");

  TCanvas *canDataMC = new TCanvas("canDataMC","Width Asymmetry",500,500);
  canDataMC->cd();
  hMC->Draw("HISTE");
  sigmaMC->Draw("same");
  hData->Draw("PE1same");
  sigmaData->Draw("same");  
  //   labelAll->Draw("same");
  //   labelEta->Draw("same");
  TLegend* leg = util::LabelFactory::createLegendColWithOffset(4,0.3,1);
  leg->AddEntry(hData,"Data","P");
  leg->AddEntry(sigmaData,"Fit Data","L");
  leg->AddEntry(hMC,"MC","L");
  leg->AddEntry(sigmaMC,"Fit MC","L");
  leg->Draw("same");
  label1->Draw("same");
  label2->Draw("same");
  canDataMC->SetLogx();
  canDataMC->SaveAs(gUID_+"_WidthAsymFit.eps","eps");

  // Parametrise ratio assuming miscalibration
  TH1 *hRatio = util::HistOps::createRatioPlot(hData,hMC);
  TH1 *hRatioFrame = util::HistOps::createRatioFrame(hRatio,"#sigma(Asymmetry) Data / MC",0.8,2.);
  TF1 *fRatio = new TF1("fRatio","sqrt( ([0]/x/x + [1]/x + [2]*[3])/([0]/x/x + [1]/x + [2]) )",min,max);
  fRatio->FixParameter(0,pow(sigmaMC->GetParameter(0),2));
  fRatio->FixParameter(1,pow(sigmaMC->GetParameter(1),2));
  fRatio->FixParameter(2,pow(sigmaMC->GetParameter(2),2));
  fRatio->SetParameter(3,1.);
  fRatio->SetLineWidth(1);
  fRatio->SetLineColor(2);
  hRatio->Fit(fRatio,"I0QRB");

  TCanvas *canRatio = new TCanvas("canRatio","Width Asymmetry Ratio",500,500);
  canRatio->cd();
  hRatioFrame->Draw();
  hRatio->Draw("PE1same");
  fRatio->Draw("same");
  label1->Draw("same");
  label2->Draw("same");
  canRatio->SetLogx();
  canRatio->SaveAs(gUID_+"_WidthAsymRatioFit.eps","eps");

  std::cout << std::endl;
  std::cout << "\\toprule" << std::endl;
  std::cout << "& $b_{0}$ & $b_{1}$ & $b_{2}$ \\\\" << std::endl;
  std::cout << "\\midrule" << std::endl;
  std::cout << "Data";
  for(int i = 0; i < sigmaData->GetNpar(); ++i) {
    std::cout << " & $" << sigmaData->GetParameter(i) << " \\pm " << sigmaData->GetParError(i) << "$";
  }
  std::cout << " \\\\" << std::endl;
  std::cout << "MC";
  for(int i = 0; i < sigmaMC->GetNpar(); ++i) {
    std::cout << " & $" << sigmaMC->GetParameter(i) << " \\pm " << sigmaMC->GetParError(i) << "$";
  }
  std::cout << " \\\\" << std::endl;
  std::cout << "\\bottomrule" << std::endl;
}




void createSlides(const std::vector<SampleAdmin*> &admins) {

  TString name = gUID_+"_Slides.tex";
  std::ofstream oFile(name);

  unsigned int nAdmins = admins.size();

  // N pt bins
  unsigned int nPtBins = 0;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    nPtBins += (*it)->nPtBins();
  }
  // Eta range
  TString etaRange = "$"+gEtaMin_+" < |\\eta| < "+gEtaMax_+"$";


  // PtAve spectra
  oFile << "% ----- PtAve Spectra ---------------------------" << std::endl;
  unsigned int nSlides = nAdmins/3;
  if( nAdmins%3 > 0 ) nSlides++;
  for(unsigned int slide = 0; slide < nSlides; ++slide) {
    oFile << "\n% --------------------------------------------------\n";
    oFile << "\\begin{frame}\n";
    oFile << "  \\frametitle{\\ptave Spectra "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
    oFile << "  \\begin{columns}[T] \n";
    for(int col = 0; col < 3; ++col) {
      oFile << "    \\begin{column}{0.3333\\textwidth} \n";
      oFile << "    \\centering\n";
      unsigned int adminIdx = 3*slide + col;
      if( adminIdx < nAdmins ) {
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "PtAve}\\\\ \n";
      }
      oFile << "  \\end{column} \n";
    }
    oFile << "  \\end{columns} \n";
    oFile << "\\end{frame} \n";
  }


  oFile << "\n\n\n% ----- Pt Asymmetry ---------------------------" << std::endl;
  nSlides = nPtBins/3;
  if( nPtBins%3 > 0 ) nSlides++;
  unsigned int adminIdx = 0;
  unsigned int ptBinIdx = 0;
  for(unsigned int slide = 0; slide < nSlides; ++slide) {
    oFile << "\n% --------------------------------------------------\n";
    oFile << "\\begin{frame}\n";
    oFile << "  \\frametitle{\\pt Asymmetry "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
    oFile << "  \\begin{columns}[T] \n";
    for(int col = 0; col < 3; ++col, ++ptBinIdx) {
      oFile << "    \\begin{column}{0.3333\\textwidth} \n";
      oFile << "    \\centering\n";
      if( !(ptBinIdx < admins[adminIdx]->nPtBins()) ) {
	ptBinIdx = 0;
	adminIdx++;
      }
      if( adminIdx < nAdmins ) {
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "PtAsym_" << ptBinIdx << "}\\\\ \n";
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "PtAsymLog_" << ptBinIdx << "}\\\\ \n";
      }
      oFile << "  \\end{column} \n";
    }
    oFile << "  \\end{columns} \n";
    oFile << "\\end{frame} \n";
  }


  oFile << "\n\n\n% ----- Pt Smeared Asymmetry ---------------------------" << std::endl;
  nSlides = nPtBins/3;
  if( nPtBins%3 > 0 ) nSlides++;
  adminIdx = 0;
  ptBinIdx = 0;
  for(unsigned int slide = 0; slide < nSlides; ++slide) {
    oFile << "\n% --------------------------------------------------\n";
    oFile << "\\begin{frame}\n";
    oFile << "  \\frametitle{Smeared \\pt Asymmetry "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
    oFile << "  \\begin{columns}[T] \n";
    for(int col = 0; col < 3; ++col, ++ptBinIdx) {
      oFile << "    \\begin{column}{0.3333\\textwidth} \n";
      oFile << "    \\centering\n";
      if( !(ptBinIdx < admins[adminIdx]->nPtBins()) ) {
	ptBinIdx = 0;
	adminIdx++;
      }
      if( adminIdx < nAdmins ) {
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "PtSmearAsym_" << ptBinIdx << "}\\\\ \n";
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "PtSmearAsymLog_" << ptBinIdx << "}\\\\ \n";
      }
      oFile << "  \\end{column} \n";
    }
    oFile << "  \\end{columns} \n";
    oFile << "\\end{frame} \n";
  }


  oFile << "\n\n\n% ----- Tails ---------------------------" << std::endl;
  nSlides = nPtBins/3;
  if( nPtBins%3 > 0 ) nSlides++;
  adminIdx = 0;
  ptBinIdx = 0;
  for(unsigned int slide = 0; slide < nSlides; ++slide) {
    oFile << "\n% --------------------------------------------------\n";
    oFile << "\\begin{frame}\n";
    oFile << "  \\frametitle{Tails "+etaRange+" (" << slide+1 << "/" << nSlides << ")}\n";
    oFile << "  \\begin{columns}[T] \n";
    for(int col = 0; col < 3; ++col, ++ptBinIdx) {
      oFile << "    \\begin{column}{0.3333\\textwidth} \n";
      oFile << "    \\centering\n";
      if( !(ptBinIdx < admins[adminIdx]->nPtBins()) ) {
	ptBinIdx = 0;
	adminIdx++;
      }
      if( adminIdx < nAdmins ) {
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "TailFit_" << ptBinIdx << "}\\\\ \n";
	oFile << "      \\includegraphics[width=\\textwidth]{" << admins.at(adminIdx)->outNamePrefix() << "Tail_" << ptBinIdx << "}\\\\ \n";
      }
      oFile << "  \\end{column} \n";
    }
    oFile << "  \\end{columns} \n";
    oFile << "\\end{frame} \n";
  }

  oFile.close();
}


void fitTailsFromAsym() {
  gErrorIgnoreLevel = 1001;
  util::StyleSettings::presentation();

  double ppCut = 0.1;

  gJetAlgo_ = "AK5 Calo-Jets";
  gDeltaPhiCut_ = "|#Delta#phi| > 2.7";
  gEtaMin_ = "0";
  gEtaMax_ = "1.1";
  gPPCut_ = "p_{||} < "+util::toTString(ppCut);
  gUID_ = "Tails_Calo_Eta00-11_Pp10";

  
  std::vector<SampleAdmin*> admins;

  //   // HLT test
  std::vector<double> ptBinEdges;
  
  //90 100 120 150 170 200 250 300 350 400 500 1000
  //0  1   2   3   4   5   6   7   8   9   10  11

  // HLT 70
  ptBinEdges.clear();
  ptBinEdges.push_back(120.);
  ptBinEdges.push_back(150.);
  admins.push_back(new SampleAdmin("input/Tails_Calo_Data_Eta00-11_Pt0120-0150_Pp10/jsResponse.root","input/Tails_Calo_MC_Eta00-11_Pt0090-1000_Pp10/jsResponse.root",ptBinEdges,2.0,ppCut,2,3));

  // HLT 100
  ptBinEdges.clear();
  ptBinEdges.push_back(150.);
  ptBinEdges.push_back(170.);
  ptBinEdges.push_back(200.);
  admins.push_back(new SampleAdmin("input/Tails_Calo_Data_Eta00-11_Pt0150-0200_Pp10/jsResponse.root","input/Tails_Calo_MC_Eta00-11_Pt0090-1000_Pp10/jsResponse.root",ptBinEdges,9.6,ppCut,3,5));

  // HLT 140
  ptBinEdges.clear();
  ptBinEdges.push_back(200.);
  ptBinEdges.push_back(250.);
  ptBinEdges.push_back(300.);
  ptBinEdges.push_back(350.);
  ptBinEdges.push_back(400.);
  ptBinEdges.push_back(500.);
  ptBinEdges.push_back(1000.);
  admins.push_back(new SampleAdmin("input/Tails_Calo_Data_Eta00-11_Pt0200-1000_Pp10/jsResponse.root","input/Tails_Calo_MC_Eta00-11_Pt0090-1000_Pp10/jsResponse.root",ptBinEdges,27.4,ppCut,5,11));



  // Concat pt bin edges
  std::vector<double> binEdgesAll;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    std::vector<double> bins = (*it)->ptBinEdges();
    if( binEdgesAll.size() == 0 ) {
      binEdgesAll = bins;
    } else {
      assert( binEdgesAll.back() == bins.front() );
      binEdgesAll.insert(binEdgesAll.end(),bins.begin()+1,bins.end());
    }
  }
  for(size_t i = 0; i < binEdgesAll.size()-1; ++i) {
    std::cout << "Bin " << i << ": " << binEdgesAll[i] << " - " << binEdgesAll[i+1] << std::endl;
  }

  // Plot spectra and asymmetry distributions
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    (*it)->plotSpectra();
    (*it)->plotAsymmetry();
    (*it)->plotSmearedAsymmetry();
    (*it)->plotTails();
  }

  // Global labels
  TPaveText* labelAll = createLabel(binEdgesAll.front(),binEdgesAll.back());
  TPaveText* labelEta = createLabelEta();

  // Compare width of asymmetry
  TH1 *hWidthAsymData = util::HistOps::createTH1D("hWidthAsymData",binEdgesAll,"p^{ave}_{T}","GeV","#sigma(Asymmetry)");
  hWidthAsymData->GetXaxis()->SetMoreLogLabels();
  TH1 *hWidthAsymMC = static_cast<TH1D*>(hWidthAsymData->Clone("hWidthAsymMC"));
  hWidthAsymData->SetMarkerStyle(20);
  int ptBin = 1;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    for(unsigned int sBin = 0; sBin < (*it)->nPtBins(); ++sBin,++ptBin) {
      hWidthAsymData->SetBinContent(ptBin,(*it)->widthAsymData(sBin));
      hWidthAsymData->SetBinError(ptBin,(*it)->widthAsymErrData(sBin));
      
      hWidthAsymMC->SetBinContent(ptBin,(*it)->widthAsymMC(sBin));
      hWidthAsymMC->SetBinError(ptBin,(*it)->widthAsymErrMC(sBin));
    }
  }
  TCanvas *canWA = new TCanvas("canWA","Width Asymmetry",500,500);
  canWA->cd();
  hWidthAsymMC->GetYaxis()->SetRangeUser(0.,0.4);
  hWidthAsymMC->Draw("HISTE");
  hWidthAsymData->Draw("PE1same");
  labelAll->Draw("same");
  labelEta->Draw("same");
  TLegend* leg = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
  leg->AddEntry(hWidthAsymData,"Data","P");
  leg->AddEntry(hWidthAsymMC,"MC","L");
  leg->Draw("same");
  canWA->SetLogx();
  canWA->SaveAs(gUID_+"_WidthAsym.eps","eps");

  TH1 *hRatioWidthAsym = util::HistOps::createRatioPlot(hWidthAsymData,hWidthAsymMC);
  TH1 *hRatioWidthAsymFrame = util::HistOps::createRatioFrame(hRatioWidthAsym,"#sigma(Asymmetry) Data / MC",0.8,1.5);
  TCanvas *canWAR = new TCanvas("canWAR","Width Asymmetry Ratio",500,500);
  canWAR->cd();
  hRatioWidthAsymFrame->Draw();
  hRatioWidthAsym->Draw("PE1same");
  labelAll->Draw("same");
  labelEta->Draw("same");
  canWAR->SetLogx();
  canWAR->SaveAs(gUID_+"_WidthAsymRatio.eps","eps");

  // Fit asymmetry --> noise term...
  plotDataMCDiff(hWidthAsymData,hWidthAsymMC,labelAll,labelEta);
  


  // Compare width of smeared asymmetry
  TH1 *hWidthSmearedAsymMC = static_cast<TH1D*>(hWidthAsymMC->Clone("hWidthSmearedAsymMC"));
  hWidthAsymMC->SetLineStyle(2);
  ptBin = 1;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    for(unsigned int sBin = 0; sBin < (*it)->nPtBins(); ++sBin,++ptBin) {
      hWidthSmearedAsymMC->SetBinContent(ptBin,(*it)->widthSmearedAsymMC(sBin));
      hWidthSmearedAsymMC->SetBinError(ptBin,(*it)->widthSmearedAsymErrMC(sBin));
    }
  }
  TCanvas *canWSA = new TCanvas("canWSA","Width Smeared Asymmetry",500,500);
  canWSA->cd();
  hWidthSmearedAsymMC->GetYaxis()->SetRangeUser(0.,0.4);
  hWidthSmearedAsymMC->Draw("HISTE");
  hWidthAsymMC->Draw("HISTEsame");
  hWidthAsymData->Draw("PE1same");
  labelAll->Draw("same");
  labelEta->Draw("same");
  TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(3,0.3,1);
  leg2->AddEntry(hWidthAsymData,"Data","P");
  leg2->AddEntry(hWidthAsymMC,"MC","L");
  leg2->AddEntry(hWidthSmearedAsymMC,"MC (Smeared)","L");
  leg2->Draw("same");
  canWSA->SetLogx();
  canWSA->SaveAs(gUID_+"_WidthSmearedAsym.eps","eps");

  TH1 *hRatioWidthSmearAsym = util::HistOps::createRatioPlot(hWidthAsymData,hWidthSmearedAsymMC);
  hRatioWidthSmearAsym->SetName("hRatioWidthSmearAsym");
  TCanvas *canWSAR = new TCanvas("canWSAR","Width Smearerd Asymmetry Ratio",500,500);
  canWSAR->cd();
  hRatioWidthAsymFrame->Draw();
  hRatioWidthAsym->SetMarkerStyle(24);
  hRatioWidthAsym->Draw("PE1same");
  hRatioWidthSmearAsym->Draw("PE1same");
  labelAll->Draw("same");
  labelEta->Draw("same");
  TLegend* leg4 = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
  leg4->AddEntry(hRatioWidthAsym,"Regular MC","P");
  leg4->AddEntry(hRatioWidthSmearAsym,"Smeared MC","P");
  leg4->Draw("same");
  canWSAR->SetLogx();
  canWSAR->SaveAs(gUID_+"_WidthSmearedAsymRatio.eps","eps");



  ///////////////////////////////////////////////////////////////
  //   Number of tail events and scaling factors
  ///////////////////////////////////////////////////////////////

  // Different error modes
  unsigned int nErrorModes = 2;
  util::HistVec hNTailData(nErrorModes);
  util::HistVec hNTailMC(nErrorModes);

  for(unsigned int i = 0; i < nErrorModes; ++i) {
    // Declare histograms of absolute number of tail events
    hNTailData[i] = util::HistOps::createTH1D("hNTailData"+util::toTString(i),binEdgesAll,
					      "p^{ave}_{T}","GeV","Number Tail Events");
    hNTailData[i]->SetMarkerStyle(20+i);
    util::HistOps::setStyleColor(hNTailData[i],i);
    hNTailMC[i] = static_cast<TH1D*>(hNTailData[i]->Clone("hNTailMC"+util::toTString(i)));
    hNTailMC[i]->SetFillColor(5);
  }

  // Fill histograms with absolute number (from data x-sec)
  std::cout << "\n\nNUMBER OF TAILS EVENTS\n";
  ptBin = 1;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    for(unsigned int sBin = 0; sBin < (*it)->nPtBins(); ++sBin,++ptBin) {
      const SampleAdmin *admin = *it;
      double nTotal = admin->nTotalData(sBin);

      // Uncertainty mode 0
      hNTailData[0]->SetBinContent(ptBin,nTotal*admin->nTailData(sBin));
      hNTailData[0]->SetBinError(ptBin,sqrt(hNTailData[0]->GetBinContent(ptBin)));
      hNTailMC[0]->SetBinContent(ptBin,nTotal*admin->nTailMC(sBin));
      hNTailMC[0]->SetBinError(ptBin,sqrt(hNTailMC[0]->GetBinContent(ptBin)));


      // Uncertainty mode 1
      hNTailData[1]->SetBinContent(ptBin,nTotal*admin->nTailData(sBin));
      hNTailData[1]->SetBinError(ptBin,0.);
      hNTailMC[1]->SetBinContent(ptBin,nTotal*admin->nTailMC(sBin));
      double sigMCStat = hNTailMC[1]->GetBinContent(ptBin)/sqrt(admin->nTotalMC(sBin)*admin->nTailMC(sBin)); 
      double sigPred = sqrt(hNTailMC[1]->GetBinContent(ptBin));
      hNTailMC[1]->SetBinError(ptBin,sqrt(sigMCStat*sigMCStat + sigPred*sigPred));

      std::cout << ptBin << " (" << admin->ptMin(sBin) << " - " << admin->ptMax(sBin) << "): \t\t" << hNTailData[0]->GetBinContent(ptBin) << " +/- " << hNTailData[0]->GetBinError(ptBin) << " (" << hNTailData[1]->GetBinError(ptBin) << "),  \t" << hNTailMC[0]->GetBinContent(ptBin) << " +/- " << hNTailMC[0]->GetBinError(ptBin) << " (" << hNTailMC[1]->GetBinError(ptBin) << ")" << std::endl;
    }
  }
  
  // Plot numbers
  TCanvas *canNTail = new TCanvas("canNTail","Number of tail events",500,500);
  canNTail->cd();
  hNTailMC[0]->GetYaxis()->SetRangeUser(0,500.);
  hNTailMC[0]->GetXaxis()->SetMoreLogLabels();
  hNTailMC[0]->Draw("HIST");
  hNTailData[0]->Draw("PE1same");
//   for(unsigned int i = 0; i < nErrorModes; ++i) {
//     hNTailMC[i]->Draw("HISTsame");
//     hNTailData[i]->Draw("PE1same");
//   }
  labelAll->Draw("same");
  labelEta->Draw("same");
  TLegend* leg3 = util::LabelFactory::createLegendColWithOffset(2,0.3,1);
  leg3->AddEntry(hNTailData[0],"Data","P");
  leg3->AddEntry(hNTailMC[0],"MC","F");
  leg3->Draw("same");
  canNTail->SetLogx();
  canNTail->SaveAs(gUID_+"_NumTailEvts.eps","eps");

  
  util::HistVec hScalFac(nErrorModes);
  for(unsigned int i = 0; i < nErrorModes; ++i) {
    hScalFac[i] = util::HistOps::createRatioPlot(hNTailData[i],hNTailMC[i],"");
  }
  TH1 *hScalFacFrame = util::HistOps::createRatioFrame(hScalFac[0],"Tail Scaling Factor",2);
  TCanvas *canScalFac = new TCanvas("canScalFac","Scaling factors",500,500);
  canScalFac->cd();
  hScalFacFrame->GetYaxis()->SetRangeUser(0.,10.);
  hScalFacFrame->GetXaxis()->SetMoreLogLabels();
  hScalFacFrame->Draw();
  for(unsigned int i = 0; i < 1; ++i) {
    hScalFac[i]->Draw("PE1same");
  }
  labelAll->Draw("same");
  labelEta->Draw("same");
  canScalFac->SetLogx();
  canScalFac->SaveAs(gUID_+"_ScalingFactors.eps","eps");


  // Write selected histograms to file
  TFile outFile((gUID_+".root"),"RECREATE");
  outFile.WriteTObject(hRatioWidthAsymFrame);
  outFile.WriteTObject(hRatioWidthAsym);
  outFile.WriteTObject(hRatioWidthSmearAsym);
  outFile.WriteTObject(hScalFacFrame);
  outFile.WriteTObject(hScalFac[0]);

  outFile.Close();


  createSlides(admins);
}

