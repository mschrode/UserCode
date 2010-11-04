// $Id: $

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TRandom3.h"

#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"




// --------------------------------------------------
class Sample {
public:
  //  enum SampleType { Data, MC };

  //Sample(const TString &fileName, SampleType type, const TString &name);
  Sample(const TString &fileName, int type, const TString &sampleName);
  ~Sample();

  TString name() const { return name_; }
  unsigned int nBins() const { return hPtAsym_.size(); }
  //  SampleType type() const { return type_; }
  int type() const { return type_; }
  
  TH1 *hPtAve() const { return hist(hPtAve_); }
  TH1 *hPtAsym(unsigned int i) const { return hist(i,hPtAsym_); }
  TH1 *hPtAsymSmeared(unsigned int i) const { return hist(i,hPtAsymSmeared_); }
  TH1 *hPtAbsAsym(unsigned int i) const { return hist(i,hPtAbsAsym_); }
  TH1 *hPtAbsAsymSmeared(unsigned int i) const { return hist(i,hPtAbsAsymSmeared_); }
  TH1 *hTailClean(unsigned int i) const { return hist(i,hTailsClean_); }
  TF1 *fTail(unsigned int i) const;
  
  double widthAsym(unsigned int i) const { return widthAsym_.at(i); }
  double widthAsymErr(unsigned int i) const { return widthAsymErr_.at(i); }  
  double widthAsymSmeared(unsigned int i) const { return widthAsymSmeared_.at(i); }
  double widthAsymErrSmeared(unsigned int i) const { return widthAsymErrSmeared_.at(i); }  

  double nTail(unsigned int i) const { return nTail_.at(i); }

private:
  //const SampleType type_;
  const int type_;
  const TString name_;

  TRandom3 *rand_;

  TH1 *hPtAve_;
  util::HistVec hPtJet_;
  util::HistVec hPtAsym_;
  util::HistVec hPtAsymSmeared_;
  util::HistVec hPtAbsAsym_;
  util::HistVec hPtAbsAsymSmeared_;

  util::HistVec hTails_;
  util::HistVec hTailsClean_;
  std::vector<TF1*> fTailFits_;

  std::vector<double> widthAsym_;
  std::vector<double> widthAsymErr_;
  std::vector<double> widthAsymSmeared_;
  std::vector<double> widthAsymErrSmeared_;

  std::vector<double> nTail_;

  mutable unsigned int count_;

  void fitWidth(const util::HistVec &h, std::vector<double> &width, std::vector<double> &widthErr) const;
  double gaussInt(double mean, double sigma, double min, double max) const;
  bool getTails(const util::HistVec &h);
  bool getTail(const TH1* hAbsAsym, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) const;
  TH1* hist(const TH1* h) const;
  TH1* hist(unsigned int i, const util::HistVec &h) const;
  void setStyle(TH1 *h) const;
  void smearAsymmetry(double scaling, const util::HistVec &hOrig, util::HistVec &hSmeared, bool abs = false) const;
  void smearHistogram(const TH1* hOrig, TH1* &hSmeared, double width, double scaling, bool abs = false) const;
};


// --------------------------------------------------
Sample::Sample(const TString &fileName, int type, const TString &sampleName) : type_(type), name_(sampleName+"_"+util::toTString(type)) {
  rand_ = new TRandom3(0);
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
    setStyle(hPtAve_);
  } else {
    ioIsOk = false;
  }

//   // Average pt
//   for(int i = 0; ioIsOk && (i < 4); ++i) {
//     TH1 *hPtJet = 0;
//     file.GetObject("hPtJet"+util::toTString(i),hPtJet);
//     if( hPtJet ) {
//       hPtJet->SetDirectory(0);
//       newName = name()+"_PtJet"+util::toTString(i);
//       hPtJet->SetName(newName);
//       setStyle(hPtJet);
//       hPtJet_.push_back(hPtJet);
//     } else {
//       ioIsOk = false;
//     }
//   }

  // Pt asymmetry
  bool binExists = true;
  int ptBin = 0;
  int entries = 0;
  while( ioIsOk && binExists ) {
    TH1 *hPtAsym = 0;
    file.GetObject("hPtAsym_"+util::toTString(ptBin),hPtAsym);
    if( hPtAsym ) {
      hPtAsym->SetDirectory(0);
      newName = name()+"_PtAsym"+util::toTString(ptBin);
      hPtAsym->SetName(newName);
      hPtAsym->GetXaxis()->SetRangeUser(-0.5,0.5);
      util::HistOps::setAxisTitles(hPtAsym,"Asymmetry","","Events");
      hPtAsym->Scale(hPtAsym->GetEntries()*hPtAsym->GetBinWidth(1));
      entries += hPtAsym->Integral();
      setStyle(hPtAsym);
      hPtAsym_.push_back(hPtAsym);
      ++ptBin;
    } else {
      binExists = false;
    }    
  }

  //  std::cout << ">>>> " << entries << " : " << hPtAve_->Integral() << std::endl;

  // Absolute pt asymmetry
  binExists = true;
  ptBin = 0;
  while( ioIsOk && binExists ) {
    TH1 *hPtAbsAsym = 0;
    file.GetObject("hPtAsymBiased_"+util::toTString(ptBin),hPtAbsAsym);
    if( hPtAbsAsym ) {
      ++ptBin;
      hPtAbsAsym->SetDirectory(0);
      newName = name()+"_PtAbsAsym"+util::toTString(ptBin);
      hPtAbsAsym->SetName(newName);
      hPtAbsAsym->Scale(hPtAbsAsym->GetEntries()*hPtAbsAsym->GetBinWidth(1));
      setStyle(hPtAbsAsym);
      hPtAbsAsym_.push_back(hPtAbsAsym);
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
    smearAsymmetry(0.1,hPtAbsAsym_,hPtAbsAsymSmeared_,true);
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
//   for(util::HistIt it = hPtAbsAsym_.begin(); it != hPtAbsAsym_.end(); ++it) {
//     delete *it;
//   }
//   for(util::HistIt it = hPtAbsAsymSmeared_.begin(); it != hPtAbsAsymSmeared_.end(); ++it) {
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
  delete rand_;
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
    if( (*it)->Fit("gaus","0QIR","",(*it)->GetMean()-2.*(*it)->GetRMS(),(*it)->GetMean()+2.*(*it)->GetRMS()) == 0 ) {
      width.push_back((*it)->GetFunction("gaus")->GetParameter(2));
      widthErr.push_back((*it)->GetFunction("gaus")->GetParError(2));
    } else {
      std::cerr << "WARNING in Sample::fitWidth (" << name() << "): No convergence when fitting width of '" << (*it)->GetName() << "'\n";
      width.push_back(0.);
      widthErr.push_back(100.);
    }
  }
}


// --------------------------------------------------
void Sample::smearAsymmetry(double scaling, const util::HistVec &hOrig, util::HistVec &hSmeared, bool abs) const {
  for(util::HistIt it = hSmeared.begin(); it != hSmeared.end(); ++it) {
    delete *it;
  }
  hSmeared.clear();
  hSmeared = util::HistVec(nBins());

  for(unsigned int i = 0; i < nBins(); ++i) {
    smearHistogram(hOrig[i],hSmeared[i],widthAsym(i),scaling,abs);
  }
}


// --------------------------------------------------
void Sample::smearHistogram(const TH1* hOrig, TH1* &hSmeared, double width, double scaling, bool abs) const {
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
    hSmeared->SetBinError(bin,sqrt(hSmeared->GetBinContent(bin)));
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
    nTail_[i] = hTailsClean_[i]->Integral();
    std::cout << "NTAIL (" << i << "): " << nTail_[i] << std::endl; 
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
  
  double sigma = hTail->GetRMS();
  success =  !(hTail->Fit(gauss,"I0Q","",0.,2.5*sigma));
  if( success ) {
    success = !(hTail->Fit(gauss,"I0Q","",0.,2.*sigma));
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
  SampleAdmin(const TString &sampleName, const TString &fileData, const TString &fileMC, std::vector<double> ptBinEdges, double lumi, unsigned int runStart, unsigned int runEnd, double ppCut);
  ~SampleAdmin();

  TString name() const { return name_; }
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

  void plotSpectra() const;
  void plotAsymmetry() const;
  void plotSmearedAsymmetry() const;
  void plotAbsAsymmetry() const;
  void plotSmearedAbsAsymmetry() const;
  void plotTails() const;
  
private:
  const TString name_;
  const std::vector<double> ptBinEdges_;
  const double lumi_;
  const unsigned int runStart_;
  const unsigned int runEnd_;
  const double ppCut_;

  Sample* data_;
  Sample* mc_;

  void plot(const TString &name, TH1 *hData, TH1 *hMC, bool logy = false) const;
};

typedef std::vector<SampleAdmin*>::const_iterator AdminIt;

SampleAdmin::SampleAdmin(const TString &sampleName, const TString &fileData, const TString &fileMC, std::vector<double> ptBinEdges, double lumi, unsigned int runStart, unsigned int runEnd, double ppCut)
  : name_(sampleName), ptBinEdges_(ptBinEdges), lumi_(lumi), runStart_(runStart), runEnd_(runEnd), ppCut_(ppCut) {
  std::cout << "Creating samples for range " << runStart << " - " << runEnd << std::endl;
  
  // Create samples holding histograms
  TString tmpName = name()+"_"+util::toTString(ptMin())+"-"+util::toTString(ptMax());
  data_ = new Sample(fileData,0,tmpName);
  mc_ = new Sample(fileMC,1,tmpName);

  std::cout << "  Data '" << fileData << "': " << data_->nBins() << " ptBins" << std::endl;
  std::cout << "  MC '" << fileData << "': " << mc_->nBins() << " ptBins" << std::endl;

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
  // PtAve spectrum
  plot(name()+"PtAve",data_->hPtAve(),mc_->hPtAve());
}


void SampleAdmin::plotAsymmetry() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TString histName = name()+"PtAsym"+util::toTString(i);
    plot(histName,data_->hPtAsym(i),mc_->hPtAsym(i));
  }
}


void SampleAdmin::plotSmearedAsymmetry() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TString canName = name()+"PtAsymSmear"+util::toTString(i);
    TCanvas *can = new TCanvas(canName,canName,500,500);
    can->cd();
    mc_->hPtAsymSmeared(i)->Draw("HISTE");
    data_->hPtAsym(i)->Draw("PE1same");

    canName = name()+"PtAsymSmear"+util::toTString(i)+"Log";
    can = new TCanvas(canName,canName,500,500);
    can->cd();
    TH1 *hMC = mc_->hPtAsymSmeared(i);
    hMC->GetXaxis()->SetRangeUser(-1.,1.);
    hMC->GetYaxis()->SetRangeUser(8E-2,10.*hMC->GetMaximum());
    hMC->Draw("HISTE");
    data_->hPtAsym(i)->Draw("PE1same");
    can->SetLogy();

    canName += "Fit";
    can = new TCanvas(canName,canName,500,500);
    can->cd();
    hMC->Draw("HISTE");
    data_->hPtAsym(i)->Draw("PE1same");
    mc_->fTail(i)->Draw("same");
    data_->fTail(i)->Draw("same");
    can->SetLogy();
  }
}


void SampleAdmin::plotAbsAsymmetry() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TString histName = name()+"PtAbsAsym"+util::toTString(i);
    plot(histName,data_->hPtAbsAsym(i),mc_->hPtAbsAsym(i),true);
  }
}


void SampleAdmin::plotSmearedAbsAsymmetry() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TString histName = name()+"PtAbsAsym"+util::toTString(i);
    plot(histName,data_->hPtAbsAsym(i),mc_->hPtAbsAsymSmeared(i),true);
  }
}


void SampleAdmin::plotTails() const {
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    TString canName = name()+"Tails"+util::toTString(i);
    plot(canName,data_->hTailClean(i),mc_->hTailClean(i));
  }
}


void SampleAdmin::plot(const TString &name, TH1 *hData, TH1 *hMC, bool logy) const {
  TCanvas *can = new TCanvas(name,name,500,500);
  can->cd();
  hMC->Draw("HISTE");
  hData->Draw("PE1same");
  if( logy ) can->SetLogy();
}





void fitTailsFromAsym() {
  util::StyleSettings::presentation();

  std::vector<SampleAdmin*> admins;

  // HLT test
  std::vector<double> ptBinEdges;
  ptBinEdges.push_back(50.);
  ptBinEdges.push_back(80.);
  ptBinEdges.push_back(120.);
  ptBinEdges.push_back(250.);
  ptBinEdges.push_back(500.);

  admins.push_back(new SampleAdmin("Test1","jsResponse.root","jsResponse.root",ptBinEdges,10.,123456,123456,0.1));

//   ptBinEdges.clear();
//   ptBinEdges.push_back(500.);
//   ptBinEdges.push_back(600.);
//   ptBinEdges.push_back(700.);
//   ptBinEdges.push_back(800.);
//   ptBinEdges.push_back(1000.);

//   admins.push_back(new SampleAdmin("Test2","jsResponse.root","jsResponse.root",ptBinEdges,20.,123456,123456,0.1));



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
//     (*it)->plotAbsAsymmetry();
//     (*it)->plotSmearedAbsAsymmetry();
    (*it)->plotTails();
  }


  // Compare width of asymmetry
  TH1 *hWidthAsymData = util::HistOps::createTH1D("hWidthAsymData",binEdgesAll,"p^{ave}_{T}","GeV","#sigma(Asymmetry)");
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
  hWidthAsymMC->GetYaxis()->SetRangeUser(0.,1.09);
  hWidthAsymMC->Draw("HISTE");
  hWidthAsymData->Draw("PE1same");


  // Compare width of smeared asymmetry
  TH1 *hWidthSmearedAsymMC = static_cast<TH1D*>(hWidthAsymMC->Clone("hWidthSmearedAsymMC"));
  hWidthSmearedAsymMC->SetYTitle("#sigma(Smeared Asymmetry)");
  ptBin = 1;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    for(unsigned int sBin = 0; sBin < (*it)->nPtBins(); ++sBin,++ptBin) {
      hWidthSmearedAsymMC->SetBinContent(ptBin,(*it)->widthSmearedAsymMC(sBin));
      hWidthSmearedAsymMC->SetBinError(ptBin,(*it)->widthSmearedAsymErrMC(sBin));
    }
  }
  TCanvas *canWSA = new TCanvas("canWSA","Width Smeared Asymmetry",500,500);
  canWSA->cd();
  hWidthSmearedAsymMC->GetYaxis()->SetRangeUser(0.,1.09);
  hWidthSmearedAsymMC->Draw("HISTE");
  hWidthAsymData->Draw("PE1same");


  // Number of tail events and scaling factors
  TH1 *hNumTailEvtsData = util::HistOps::createTH1D("hNumTailEvtsData",binEdgesAll,"p^{ave}_{T}","GeV","Number Tail Events");
  TH1 *hNumTailEvtsMC = static_cast<TH1D*>(hNumTailEvtsData->Clone("hNumTailEvtsMC"));
  hNumTailEvtsData->SetMarkerStyle(20);
  hNumTailEvtsMC->SetFillColor(5);
  ptBin = 1;
  for(AdminIt it = admins.begin(); it != admins.end(); ++it) {
    for(unsigned int sBin = 0; sBin < (*it)->nPtBins(); ++sBin,++ptBin) {
      std::cout << sBin << " --> " << (*it)->nTailData(sBin) << std::endl;
      hNumTailEvtsData->SetBinContent(ptBin,(*it)->nTailData(sBin));
      hNumTailEvtsMC->SetBinContent(ptBin,(*it)->nTailMC(sBin));
    }
  }
  
  TCanvas *canNumTailEvt = new TCanvas("canNumTailEvt","Number of tail events",500,500);
  canNumTailEvt->cd();
  hNumTailEvtsMC->Draw("HISTE");
  hNumTailEvtsData->Draw("PE1same");
  
  TH1 *hScalFac = util::HistOps::createRatioPlot(hNumTailEvtsData,hNumTailEvtsMC,"");
  TH1 *hScalFacFrame = util::HistOps::createRatioFrame(hScalFac,"Tail Scaling Factor",2);
  TCanvas *canScalFac = new TCanvas("canScalFac","Scaling factors",500,500);
  canScalFac->cd();
  hScalFacFrame->Draw();
  hScalFac->Draw("PE1same");
}



// START OLD SCRIPT

// // --------------------------------------------------
// bool getScaleFactorIntegral(const TH1* hTailData, const TH1* hTailMC, double &scale, double &uncert) {
//   bool success = false;
//   double nData = hTailData->Integral();
//   double nMC = hTailMC->Integral();
// //   std::cout << "N(data) = " << nData << std::endl;
// //   std::cout << "N(MC) = " << nMC << std::endl;
//   if( nData > 0. && nMC > 0. ) {
//     success = true;
//     scale = nData / nMC;
//     double errData = sqrt(nData);
//     double errMC = sqrt(nMC);
//     uncert = sqrt( errData*errData/nMC/nMC + errMC*errMC*nData*nData/nMC/nMC/nMC/nMC );
//   } else {
//     scale = 1.;
//     uncert = 1000.;
//   }

//   return success;
// }



// // --------------------------------------------------
// bool getScaleFactorRatioFit(const TH1* hTailData, const TH1* hTailMC, TH1* &hRatio, TF1 *&fit, double &scale, double &uncert) {
//   hRatio = util::HistOps::createRatioPlot(hTailData,hTailMC,"Data / MC (Asymmetry Tail)",1);
//   TString name = hRatio->GetName();
//   name += "Fit";
//   fit = new TF1(name,"pol0",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   bool success = !(hRatio->Fit(fit,"0QR"));
//   if( success ) {
//     scale = fit->GetParameter(0);
//     uncert = fit->GetParError(0);
//   } else {
//     scale = 1.;
//     uncert = 1000.;
//   }

//   return success;
// }



// // --------------------------------------------------
// bool getTail(TH1* hAbsAsym, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) {
//   bool success = false;

//   TString name = hAbsAsym->GetName();
//   name += "_Tails";
//   hTail = static_cast<TH1D*>(hAbsAsym->Clone(name));
//   hTail->SetLineColor(kBlue);
//   hTail->SetMarkerColor(kBlue);

//   name = hAbsAsym->GetName();
//   name += "_TailsClean";
//   hTailClean = static_cast<TH1D*>(hAbsAsym->Clone(name));
//   hTailClean->Reset();
//   hTailClean->SetLineColor(kBlue+3);
//   hTailClean->SetMarkerColor(kBlue+3);

//   name = hAbsAsym->GetName();
//   name += "_GaussFit";
//   gauss = new TF1(name,"gaus",0.,2.);
//   gauss->SetLineWidth(1);
//   gauss->SetLineColor(kRed);

//   double sigma = hTail->GetRMS();
//   success =  !(hTail->Fit(gauss,"I0Q","",0.,2.5*sigma));
//   if( success ) {
//     success = !(hTail->Fit(gauss,"I0Q","",0.,2.*sigma));
//     if( success ) {
//       sigma = gauss->GetParameter(2);
//       int binMin = 1;
//       int binMax = hTail->FindBin(5.*sigma);
//       for(int bin = binMin; bin <= binMax; ++bin) {
// 	double min = hTail->GetXaxis()->GetBinLowEdge(bin);
// 	double max = hTail->GetXaxis()->GetBinUpEdge(bin);
// 	double gaussPdf = gauss->Integral(min,max)/hTail->GetBinWidth(1);
// 	double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
// 	hTail->SetBinContent(bin,tailPdf);
//       }
//       for(int bin = 1; bin <= hTailClean->GetNbinsX(); ++bin) {
// 	if( hTail->GetBinContent(bin) > 0.4*gauss->Eval(hTail->GetBinCenter(bin)) ) {
// 	  hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
// 	  hTailClean->SetBinError(bin,hTail->GetBinError(bin));
// 	}
//       }
//     }
//   }
  
//   return success;
// }



// // --------------------------------------------------
// void readTailHistograms(const TString &name, util::HistVec &hData, util::HistVec &hMC) {
//   std::cout << "Getting histograms from file '" << name << "'" << std::endl;

//   TFile file(name,"READ");

//   TString histNameBase = "hPtBiasAsymCorr";
//   bool histExists = true;
//   int ptBin = 3;
//   while( histExists ) {
//     TH1 *hDataTmp = 0;
//     TH1 *hMCTmp = 0;
//     file.GetObject(histNameBase+"0"+util::toTString(ptBin),hDataTmp);
//     file.GetObject(histNameBase+"1"+util::toTString(ptBin)+"Smeared",hMCTmp);
//     //file.GetObject(histNameBase+"1"+util::toTString(ptBin),hMCTmp);
//     if( hDataTmp && hMCTmp ) {
//       ptBin++;
//       hDataTmp->SetDirectory(0);
//       hMCTmp->SetDirectory(0);
//       hData.push_back(hDataTmp);
//       hMC.push_back(hMCTmp);
//     } else {
//       histExists = false;
//     }    
//   }
//   file.Close();  
// }




// // --------------------------------------------------
// void fitTailsFromAsym(const TString &name = "Test_Calo_Eta13-30_Pt3Rel10_.root") {
//   util::StyleSettings::presentationNoTitle();
//   // Do not print ROOT message if eps file has been created
//   gErrorIgnoreLevel = 1001;

// //   TString outPrefix = "TailFit_WithSmearing_"+name;
// //   for(int i = 0; i < 5; ++i) outPrefix.Chop();

//   TString outPrefix = "Bla_";

  
//   // Create histograms of tails
//   util::HistVec hAbsAsymData;
//   util::HistVec hAbsAsymMC;
//   readTailHistograms(name,hAbsAsymData,hAbsAsymMC);

//   // Scale factors
//   std::vector<double> ptBinEdges;
// //   ptBinEdges.push_back(80.);
// //   ptBinEdges.push_back(100.);
// //   ptBinEdges.push_back(120.);
//   ptBinEdges.push_back(140.);
//   ptBinEdges.push_back(170.);
//   ptBinEdges.push_back(200.);
//   ptBinEdges.push_back(250.);
//   ptBinEdges.push_back(300.);
//   ptBinEdges.push_back(400.);
//   ptBinEdges.push_back(1000.);

//   if( ptBinEdges.size() != (hAbsAsymData.size()+1) ) {
//     std::cerr << "ERROR: check number of pt bins:" << std::endl;
//     std::cerr << "  hAbsAsymData.size() = " << hAbsAsymData.size() << std::endl;
//     std::cerr << "  ptBinEdges.size() = " << ptBinEdges.size() << std::endl;
//     exit(-1);
//   }
  
//   TH1 *hScalesInt = new TH1D("hScalesInt",";p^{ave}_{T} (GeV);Scaling Factor",ptBinEdges.size()-1,&(ptBinEdges.front()));
//   hScalesInt->SetMarkerStyle(20);

//   for(size_t ptBin = 0; ptBin < hAbsAsymData.size(); ++ptBin) {
//     std::cout << "BIN " << ptBin << ": " << ptBinEdges[ptBin] << " - " << ptBinEdges[ptBin+1] << std::endl;

//     // Get tails
//     TH1 *hTailData = 0;
//     TH1 *hTailCleanData = 0;
//     TF1 *fGaussData = 0;
//     if( !getTail(hAbsAsymData[ptBin],hTailData,hTailCleanData,fGaussData) ) {
//       std::cerr << "ERROR getting tail from '" << hAbsAsymData[ptBin]->GetName() << "'" << std::endl;
//       //exit(-1);
//     }
//     TH1 *hTailMC = 0;
//     TH1 *hTailCleanMC = 0;
//     TF1 *fGaussMC = 0;
//     if( !getTail(hAbsAsymMC[ptBin],hTailMC,hTailCleanMC,fGaussMC) ) {
//       std::cerr << "ERROR getting tail from '" << hAbsAsymMC[ptBin]->GetName() << "'" << std::endl;
//       //exit(-1);
//     }
//     fGaussMC->SetLineStyle(2);

//     std::cout << "  N(Data) = " << hAbsAsymData[ptBin]->Integral() << ", N(MC) = " << hAbsAsymMC[ptBin]->Integral() << std::endl;

//     // Plot asymmetry and tails
//     TCanvas *can1 = new TCanvas("can1","AbsAsym "+util::toTString(ptBin),500,500);
//     can1->cd();
//     util::HistOps::setYRange(hAbsAsymMC[ptBin],1,3E-1);
//     hAbsAsymMC[ptBin]->Draw("HISTE");
//     hAbsAsymData[ptBin]->Draw("PE1same");
//     fGaussMC->Draw("same");
//     fGaussData->Draw("same");
//     can1->SetLogy(1);
//     can1->SaveAs(outPrefix+"AbsAsym_PtBin"+util::toTString(ptBin)+".eps","eps");
//     delete can1;

//     TCanvas *can2 = new TCanvas("can2","AbsAsym "+util::toTString(ptBin),500,500);
//     can2->cd();
//     util::HistOps::setYRange(hTailCleanMC);
//     hTailCleanMC->Draw("HISTE");
//     hTailCleanData->Draw("PE1same");
//     can2->SetLogy(0);
//     can2->SaveAs(outPrefix+"Tails_PtBin"+util::toTString(ptBin)+".eps","eps");
//     delete can2;

//     // Scale factor from integral
//     double scaleInt = 0.;
//     double uncertInt = 0.;
//     if( getScaleFactorIntegral(hTailCleanData,hTailCleanMC,scaleInt,uncertInt) ) {
//       std::cout << "  Scale(Int) = " << scaleInt << " +/- " << uncertInt << std::endl;
//       hScalesInt->SetBinContent(1+ptBin,scaleInt);
//       hScalesInt->SetBinError(1+ptBin,uncertInt);
//     }

//     // Fit scale factor
//     TH1* hRatio = 0;
//     TF1 *fit = 0;
//     double scaleFit = 0.;
//     double uncertFit = 0.;
//     if( getScaleFactorRatioFit(hTailCleanData,hTailCleanMC,hRatio,fit,scaleFit,uncertFit) ) {
//       std::cout << "  Scale(Fit) = " << scaleFit << " +/- " << uncertFit << std::endl;
//       TCanvas *can = new TCanvas("can","Scales "+util::toTString(ptBin),500,500);
//       can->cd();
//       util::HistOps::createRatioFrame(hRatio,"Scale factor",1)->Draw();
//       hRatio->Draw("PE1same");
//       fit->Draw("same");
//       can->SetLogy(0);
//       can->SaveAs(outPrefix+"Scale_PtBin"+util::toTString(ptBin)+".eps","eps");
//       delete can;
//     }
//   }

//   TCanvas *can = new TCanvas("can","Scaling factors",500,500);
//   can->cd();
//   util::HistOps::createRatioFrame(hScalesInt,1)->Draw();
//   hScalesInt->Draw("PE1same");
//   can->SaveAs(outPrefix+"ScalingFactors.eps","eps");
//   //  delete can;

//   std::cout << std::endl;
//   for(int bin = 1; bin <= hScalesInt->GetNbinsX(); ++bin) {
//     std::cout << bin << " (" << ptBinEdges[bin-1] << " - " << ptBinEdges[bin] << "):" << hScalesInt->GetBinContent(bin) << std::endl;
//   }
// }


