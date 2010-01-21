// $Id: ControlPlotsProfile.cc,v 1.2 2010/01/14 13:13:17 mschrode Exp $

#include "ControlPlotsProfile.h"

#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"

#include "CalibData.h"
#include "ControlPlotsFunction.h"



//! \param config Configuration file
//! \param function The functions which return the x, y, and binning
//!                 variables from an event (see \p ControlPlotsFunction)
// ----------------------------------------------------------------   
ControlPlotsProfile::ControlPlotsProfile(const ControlPlotsConfig *config, const ControlPlotsFunction *function)
  : config_(config), function_(function) {
  if( !function_->isInit() ) {
    std::cerr << "ERROR: Initialization of ControlPlotsProfile with non-initialized Function.\n";
  }

  // Set up Bin objects
  bins_ = std::vector<Bin*>(config_->nBins());
  for(int i = 0; i < config_->nBins(); i++) {
    bins_.at(i) = new Bin(i,config_->binEdges()->at(i),
			  config_->binEdges()->at(i+1),
			  config_);
  }
  for(size_t i = 0; i < bins_.size(); i++) {
    assert( bins_.size() > 0 && bins_.at(i)->min() < bins_.at(i)->max() );
  }
  
  // Create histogram for x spectrum
  std::string name = config_->name();
  name += "_" + config_->xVariable() + "Spectrum";
  hXSpectrum_ = new TH1D(name.c_str(),"",
			 config_->nXBins(),
			 config_->xMin(),
			 config_->xMax());
  hXSpectrum_->SetXTitle((config_->xTitle()).c_str());
  hXSpectrum_->SetYTitle("Number of events");
  hXSpectrum_->GetXaxis()->SetNdivisions(505);
}



// ----------------------------------------------------------------   
ControlPlotsProfile::~ControlPlotsProfile() {
  for(std::vector<Bin*>::iterator it = bins_.begin(); 
      it != bins_.end(); it++) {
    delete *it;
  }
  bins_.clear();
  delete hXSpectrum_;
}



//! Draws the 2D histograms, the (zoomed) profile histograms, and if
//! specified the y distributions and writes them to .eps file. The
//! histograms are also stored in a .root file.
// ----------------------------------------------------------------   
void ControlPlotsProfile::draw() {
  TCanvas *c1 = new TCanvas("c1","",500,500);
  c1->cd();
  std::string fileName;

  // Draw 2D histograms
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    ControlPlotsConfig::CorrectionTypeIt corrTypeIt = config_->correctionTypesBegin();
    for(; corrTypeIt != config_->correctionTypesEnd(); corrTypeIt++) {
      c1->Clear();
      TH2D *h = (*binIt)->hYvsX(*corrTypeIt);
      h->Draw("COLZ");
      if( config_->logX() ) c1->SetLogx(1);
      c1->RedrawAxis();

      fileName = config_->outDirName() + "/";
      fileName += (*binIt)->hist2DFileName(*corrTypeIt) + "." + config_->outFileType();
      c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }
  }

  // Draw profile histograms
  TLegend *leg = bins_.front()->createLegend();
  TLine *hLine = bins_.front()->createHorizontalLine();
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    ControlPlotsConfig::ProfileTypeIt profTypeIt = config_->profileTypesBegin();
    for(; profTypeIt != config_->profileTypesEnd(); profTypeIt++) {
      c1->Clear();
      bool firstHist = true;
      ControlPlotsConfig::CorrectionTypeIt corrTypeIt = config_->correctionTypesBegin();
      for(; corrTypeIt != config_->correctionTypesEnd(); corrTypeIt++) {
	TH1D *h = (*binIt)->hXProfile(*corrTypeIt,*profTypeIt);
	if( firstHist ) {
	  h->Draw("PE1");
	  firstHist = false;
	} else {
	  h->Draw("PE1same");
	}
	config_->toRootFile(h);
      }
      hLine->Draw("same");
      leg->Draw("same");

      if( config_->logX() ) c1->SetLogx(1);
      c1->RedrawAxis();

      fileName = config_->outDirName() + "/";
      fileName += (*binIt)->profileFileName(*profTypeIt) + "." + config_->outFileType();
      c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }
  }  

  // Draw profile histograms (zoom on y axis)
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    ControlPlotsConfig::ProfileTypeIt profTypeIt = config_->profileTypesBegin();
    for(; profTypeIt != config_->profileTypesEnd(); profTypeIt++) {
      c1->Clear();
      bool firstHist = true;
      ControlPlotsConfig::CorrectionTypeIt corrTypeIt = config_->correctionTypesBegin();
      for(; corrTypeIt != config_->correctionTypesEnd(); corrTypeIt++) {
	TH1D *h = (*binIt)->hXProfile(*corrTypeIt,*profTypeIt);
	h->GetYaxis()->SetRangeUser(config_->yMinZoom(),config_->yMaxZoom());
	if( firstHist ) {
	  h->Draw("PE1");
	  firstHist = false;
	} else {
	  h->Draw("PE1same");
	}
      }
      hLine->Draw("same");
      leg->Draw("same");

      if( config_->logX() ) c1->SetLogx(1);
      c1->RedrawAxis();

      fileName = config_->outDirName() + "/";
      fileName += (*binIt)->profileFileName(*profTypeIt) + "_zoom." + config_->outFileType();
      c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }
  }  

  // Draw distributions
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    ControlPlotsConfig::CorrectionTypeIt corrTypeIt = config_->distributionCorrectionTypesBegin();
    for(; corrTypeIt != config_->distributionCorrectionTypesEnd(); corrTypeIt++) {
      for(int n = 0; n < (*binIt)->nDistributions(); n++) {
	c1->Clear();
	TH1D *h = (*binIt)->hYDistribution(n,*corrTypeIt);
	h->Draw("HIST");
	c1->SetLogx(0);
	c1->RedrawAxis();

	fileName = config_->outDirName() + "/";
	fileName += (*binIt)->distributionFileName(n,*corrTypeIt) + "." + config_->outFileType();
	c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
      }
    }
  }

  // Draw x spectrum
  c1->Clear();
  c1->cd();
  hXSpectrum_->Draw();
  c1->SetLogx(0);
  c1->SetLogy(1);

  fileName = config_->outDirName() + "/";
  fileName += hXSpectrum_->GetName();
  fileName += "." + config_->outFileType();
  c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
  config_->toRootFile(hXSpectrum_);

  // Clean up
  delete c1;
  delete leg;
  delete hLine;
}



//! This method finds the bin into which the event \p evt falls
//! from \p ControlPlotsFunction::binValue(evt). The 2D histograms
//! in this bin are filled for the different correction types as
//! y vs x with the weight w. The values of x and y are given by
//! \p ControlPlotsFunction::xValue(evt) and
//! \p ControlPlotsFunction::yValue(evt).
//! \sa Bin::fill(), findBin()
// ----------------------------------------------------------------   
void ControlPlotsProfile::fill(const Event * evt) {
  double x = function_->xValue(evt);
  hXSpectrum_->Fill(x,evt->GetWeight());

  int bin = findBin(evt);
  if( bin >= 0 ) {
    ControlPlotsConfig::CorrectionTypeIt corrTypeIt = config_->correctionTypesBegin();
    for(; corrTypeIt != config_->correctionTypesEnd(); corrTypeIt++) {
      double y = function_->yValue(evt,*corrTypeIt);
      if( bins_.at(bin)->fill(x,y,evt->GetWeight(),*corrTypeIt) ) {
	std::cerr << "ERROR when filling YvsX histograms of CorrectionType '" << *corrTypeIt << "'\n";
      }
    }
  }
}



//! see also \p ControlPlotsProfile::Bin::fitProfiles()
// ----------------------------------------------------------------   
void ControlPlotsProfile::fitProfiles() {
  for(std::vector<Bin*>::iterator it = bins_.begin(); it != bins_.end(); it++) {
    (*it)->fitProfiles();
  }
}



//! The value of the binning variable for this event \p evt
//! is given by \p ControlPlotsFunction::binValue(evt).
// ----------------------------------------------------------------   
int ControlPlotsProfile::findBin(const Event * evt) const {
  int bin = -1;
  double binValue = function_->binValue(evt);
  if( binValue >= config_->min() && binValue <= config_->max() ) {
    for(int i = 0; i < static_cast<int>(bins_.size()); i++) {
      bin = i;
      if( binValue <= bins_.at(i)->max() ) break;
    }
  }
  return bin;
}



//! \param binIdx Index of this bin
//! \param min Minimum value of the binning variable in this bin
//! \param max Maximum value of the binning variable in this bin
//! \param config Configuration parameters
// ----------------------------------------------------------------   
ControlPlotsProfile::Bin::Bin(int binIdx, double min, double max, const ControlPlotsConfig *config)
  : idx_(binIdx), min_(min), max_(max), config_(config), nCallsFitProfiles_(0) {

  // Create one 2D histogram per correction type
  ControlPlotsConfig::CorrectionTypeIt corrIt;
  for(corrIt = config_->correctionTypesBegin();
      corrIt != config_->correctionTypesEnd(); corrIt++) {
    std::string name = config_->name();
    name += "_"+config_->yVariable()+"Vs"+config_->xVariable();
    name += "_"+config_->correctionTypeName(*corrIt);
    name += "_"+config_->binName(idx_);

    TH2D * h = new TH2D(name.c_str(),"",
			config_->nXBins(),&(config_->xBinEdges()->front()),
			config_->nYBins(),config_->yMin(),config_->yMax());
    h->SetTitle((config_->binTitle(min,max)).c_str());
    h->SetXTitle((config_->xTitle()).c_str());
    h->SetYTitle((config_->yTitle()).c_str());
    h->SetMarkerColor(config_->color(*corrIt));
    hYvxX_[*corrIt] = h;
  }
}



// ----------------------------------------------------------------   
ControlPlotsProfile::Bin::~Bin() {
  for(std::map<ControlPlotsConfig::CorrectionType,TH2D*>::iterator it = hYvxX_.begin();
      it != hYvxX_.end(); it++) {
    delete it->second;
  }
  std::map< ControlPlotsConfig::CorrectionType, std::map< ControlPlotsConfig::ProfileType, TH1D*> >::iterator corrIt = hXProfile_.begin();
  for(; corrIt != hXProfile_.end(); corrIt++) {
    std::map< ControlPlotsConfig::ProfileType, TH1D*>::iterator profIt = corrIt->second.begin();
    for(; profIt != corrIt->second.end(); profIt++) {
      delete profIt->second;
    }
    std::vector<TH1D*>::iterator distIt = hYDistributions_[corrIt->first].begin();
    for(; distIt != hYDistributions_[corrIt->first].end(); distIt++) {
      delete *distIt;
    }
  }
}



// ----------------------------------------------------------------   
TH2D *ControlPlotsProfile::Bin::hYvsX(ControlPlotsConfig::CorrectionType corrType) {
  TH2D * h = 0;
  std::map<ControlPlotsConfig::CorrectionType,TH2D*>::iterator it = hYvxX_.find(corrType);
  if( it != hYvxX_.end() ) {
    h = it->second;
  }

  return h;
}



// ----------------------------------------------------------------   
TH1D *ControlPlotsProfile::Bin::hXProfile(ControlPlotsConfig::CorrectionType corrType, ControlPlotsConfig::ProfileType profType) {
  TH1D * h = 0;
  std::map< ControlPlotsConfig::CorrectionType, std::map< ControlPlotsConfig::ProfileType, TH1D* > >::iterator corrIt = hXProfile_.find(corrType);
  if( corrIt != hXProfile_.end() ) {
    std::map< ControlPlotsConfig::ProfileType, TH1D* >::iterator profIt = corrIt->second.find(profType);
    if( profIt != corrIt->second.end() ) {
      h = profIt->second;
    }
  }

  return h;
}



// ----------------------------------------------------------------   
TH1D *ControlPlotsProfile::Bin::hYDistribution(int n, ControlPlotsConfig::CorrectionType corrType) {
  TH1D * h = 0;
  std::map< ControlPlotsConfig::CorrectionType, std::vector< TH1D* > >::iterator corrIt = hYDistributions_.find(corrType);
  if( corrIt != hYDistributions_.end() ) {
    if( n >=0 && n < nDistributions() ) {
      h = corrIt->second.at(n);
    }
  }

  return h;
}



// ----------------------------------------------------------------   
int ControlPlotsProfile::Bin::fill(double x, double y, double w, ControlPlotsConfig::CorrectionType type) {
  int status = -1;
  std::map<ControlPlotsConfig::CorrectionType,TH2D*>::iterator it = hYvxX_.find(type);
  if( it != hYvxX_.end() ) {
    it->second->Fill(x,y,w);
    status = 0;
  }
    
  return status;
}



//!  Creates the y distributions for different x bins from the
//!  2D histograms. Different quantities are calculated from 
//!  these distributions and stored in the profiles vs x:
//!  - Mean
//!  - StandardDeviation
//!  - GaussFitMean
//!  - GaussFitWidth
//!  - Median
//!  - Chi2
//!  - Probability
//!  - Quantiles
//---------------------------------------------------------------
int ControlPlotsProfile::Bin::fitProfiles() {
  nCallsFitProfiles_++;
  if( nCallsFitProfiles_ > 1 ) {
    std::cerr << "WARNING: FitProfiles has been called already\n";
    // Deleting profiles
    std::map< ControlPlotsConfig::CorrectionType, std::map< ControlPlotsConfig::ProfileType, TH1D*> >::iterator corrIt = hXProfile_.begin();
    for(; corrIt != hXProfile_.end(); corrIt++) {
      std::map< ControlPlotsConfig::ProfileType, TH1D*>::iterator profIt = corrIt->second.begin();
      for(; profIt != corrIt->second.end(); profIt++) {
	delete profIt->second;
      }
      std::vector<TH1D*>::iterator distIt = hYDistributions_[corrIt->first].begin();
      for(; distIt != hYDistributions_[corrIt->first].end(); distIt++) {
	delete *distIt;
      }
      corrIt->second.clear();
    }
    hXProfile_.clear();
    hYDistributions_.clear();
  }
  int status = 0; // No error handling implement yet

  // Create profile histograms
  // Loop over CorrectionTypes
  for(std::map<ControlPlotsConfig::CorrectionType,TH2D*>::iterator corrIt = hYvxX_.begin();
      corrIt != hYvxX_.end(); corrIt++) {

    hYDistributions_[corrIt->first] = std::vector<TH1D*>(corrIt->second->GetNbinsX());
    std::map< ControlPlotsConfig::ProfileType, TH1D*> profMap;
    // Loop over *all possible* profile quantities
    for(int i = 0; i < ControlPlotsConfig::nProfileTypes; i++) {
      ControlPlotsConfig::ProfileType profType = static_cast<ControlPlotsConfig::ProfileType>(i);
      std::string name = corrIt->second->GetName();
      name += "_"+config_->profileTypeName(profType);
      TH1D *h = new TH1D(name.c_str(),"",
			 corrIt->second->GetNbinsX(),
			 corrIt->second->GetXaxis()->GetXbins()->GetArray());
      h->GetYaxis()->SetRangeUser(config_->yMin(),config_->yMax());
      h->SetTitle((config_->binTitle(min(),max())).c_str());
      h->SetXTitle((config_->xTitle()).c_str());
      h->SetYTitle((config_->yProfileTitle(profType)).c_str());
      h->SetMarkerStyle(config_->markerStyle(corrIt->first));
      h->SetMarkerColor(config_->color(corrIt->first));
      h->SetLineColor(config_->color(corrIt->first));	    

      // Add histogram to map
      profMap[profType] = h;
    } // End of loop over profile quantities

    hXProfile_[corrIt->first] = profMap;
  } // End of loop over CorrectionTypes

    // Fill profile histograms
    // Loop over CorrectionTypes
  for(std::map<ControlPlotsConfig::CorrectionType,TH2D*>::iterator corrIt = hYvxX_.begin();
      corrIt != hYvxX_.end(); corrIt++) {
    // Temp histogram to store slice of 2D hists
    TH1D *htemp = new TH1D("htemp","",
			   corrIt->second->GetNbinsY(),
			   corrIt->second->GetYaxis()->GetXmin(),
			   corrIt->second->GetYaxis()->GetXmax());
    htemp->Sumw2();

    const int nq = 2;
    double yq[2],xq[2];
    xq[0] = 0.5;
    xq[1] = 0.90;
    // Loop over x bins
    for(int xBin = 1; xBin <= corrIt->second->GetNbinsX(); xBin++) {
      htemp->Reset();
      // Copy y bin content in this x bins to htemp
      for(int yBin = 1; yBin <= corrIt->second->GetNbinsY(); yBin++) {
	htemp->SetBinContent(yBin,corrIt->second->GetBinContent(corrIt->second->GetBin(xBin,yBin)));
	htemp->SetBinError(yBin,corrIt->second->GetBinError(xBin,yBin));
      }
      // Store distribution
      char name[100];
      sprintf(name,"%s_Distribution_%i",corrIt->second->GetName(),xBin);
      TH1D *h = static_cast<TH1D*>(htemp->Clone(name));
      h->SetTitle((config_->xBinTitle(xBin-1,min(),max())).c_str());
      h->SetXTitle((config_->yTitle()).c_str());
      h->SetYTitle("N");
      h->SetLineColor(config_->color(corrIt->first));	   
      hYDistributions_[corrIt->first].at(xBin-1) = h;
        
      double mean = htemp->GetMean(); 
      double meanerror = htemp->GetMeanError();
      double width = htemp->GetRMS();
      if(width < 0.1) width = 0.1;
      if(htemp->GetSumOfWeights() <= 0) {
	continue; 
      } else {
	htemp->Fit("gaus","QNO","", mean - 3 * width,mean + 3 * width);
	TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
	mean = f->GetParameter(1);
	meanerror = f->GetParError(1);
	width = f->GetParameter(2);
	if(width < 0.05) width = 0.05;
	if( (htemp->Fit(f,"LLQNO","goff",mean - 1.5 * width, mean + 1.5 * width) == 0) ) {
	  mean = f->GetParameter(1);
	  meanerror = f->GetParError(1);
	  width = f->GetParameter(2);
	  
	  hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->SetBinContent(xBin,mean);
	  hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->SetBinError(xBin,meanerror);
	  hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetBinContent(xBin,width/mean);
	  hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetBinError(xBin,f->GetParError(2)/mean);
	}

	hXProfile_[corrIt->first][ControlPlotsConfig::Chi2]->SetBinContent(xBin, f->GetChisquare() / f->GetNumberFreeParameters());
	hXProfile_[corrIt->first][ControlPlotsConfig::Chi2]->SetBinError(xBin, 0.01);
	hXProfile_[corrIt->first][ControlPlotsConfig::Probability]->SetBinContent(xBin, f->GetProb());
	hXProfile_[corrIt->first][ControlPlotsConfig::Probability]->SetBinError(xBin, 0.01);
	mean = htemp->GetMean();
	meanerror = htemp->GetMeanError();
	width = htemp->GetRMS();
	hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->SetBinContent(xBin,mean);
	hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->SetBinError(xBin,meanerror);
	hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetBinContent(xBin,width/mean); 
	hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetBinError(xBin,htemp->GetRMSError()/mean);
	htemp->GetQuantiles(nq,yq,xq);
	hXProfile_[corrIt->first][ControlPlotsConfig::Median]->SetBinContent(xBin,yq[0]);
	hXProfile_[corrIt->first][ControlPlotsConfig::Median]->SetBinError(xBin,0.0001);
	hXProfile_[corrIt->first][ControlPlotsConfig::Quantiles]->SetBinContent(xBin,yq[1]/yq[0]-1);
	hXProfile_[corrIt->first][ControlPlotsConfig::Quantiles]->SetBinError(xBin,0.0001);
	delete f;
      }
    } // End of loop over x bins
    delete htemp;
  } // End of loop over CorrectionTypes

  return status;
}



// ----------------------------------------------------------------   
std::string ControlPlotsProfile::Bin::hist2DFileName(ControlPlotsConfig::CorrectionType type) const {
  std::string name = config_->name();
  name += "_2D_"+config_->binName(idx_);
  name += "_"+config_->correctionTypeName(type);

  return name;
}



// ----------------------------------------------------------------   
std::string ControlPlotsProfile::Bin::profileFileName(ControlPlotsConfig::ProfileType type) const {
  std::string name = config_->name();
  name += "_"+config_->profileTypeName(type);
  name += "_"+config_->binName(idx_);

  return name;
}



// ----------------------------------------------------------------   
std::string ControlPlotsProfile::Bin::distributionFileName(int xBin, ControlPlotsConfig::CorrectionType type) const {
  std::string name = config_->name();
  name += "_Dist_"+config_->binName(idx_);
  name += "_"+config_->correctionTypeName(type);
  name += "_"+config_->xBinName(xBin);

  return name;
}



// ----------------------------------------------------------------   
TLine *ControlPlotsProfile::Bin::createHorizontalLine() const {
  double y = 0.;
  if( (config_->yVariable()).find("Response") != std::string::npos )
    y = 1.;

  TLine *line = new TLine(config_->xMin(),y,config_->xMax(),y);
  line->SetLineStyle(2);
  line->SetLineColor(1);

  return line;
}



// ----------------------------------------------------------------   
TLegend *ControlPlotsProfile::Bin::createLegend() {
  size_t nEntries = hXProfile_.size();
  TLegend * leg = 0;
  if( nEntries ) {
    leg = new TLegend(0.5,0.8-nEntries*0.07,0.8,0.8);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);

    ControlPlotsConfig::ProfileTypeIt profTypeIt = config_->profileTypesBegin();
    ControlPlotsConfig::CorrectionTypeIt corrTypeIt = config_->correctionTypesBegin();
    for(; corrTypeIt != config_->correctionTypesEnd(); corrTypeIt++) {
      leg->AddEntry(hXProfile(*corrTypeIt,*profTypeIt),
		    (config_->legendLabel(*corrTypeIt)).c_str(),
		    "PL");
    }
  }

  return leg;
}
