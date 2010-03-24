// $Id: ControlPlotsJetSmearing.cc,v 1.16 2010/02/25 15:28:18 mschrode Exp $

#include "ControlPlotsJetSmearing.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TROOT.h"


#include "SmearData.h"
#include "SmearDiJet.h"
#include "SmearPhotonJet.h"
#include "Jet.h"



//!  \brief Constructor
//!
//!  By default, all controlplots are written to the
//!  directory \p ./controlplots/
//!
//!  \param configfile Name of the configuration file
//!  \param data The data
//!  \param param The parametrization
// --------------------------------------------------
ControlPlotsJetSmearing::ControlPlotsJetSmearing(const std::string& configfile, const std::vector<Event*> * data, TParameters * param, const std::string &outDir)
  : data_(data),
    config_(new ConfigFile(configfile.c_str())),
    param_(param),
    respNBins_(60),
    respMin_(0.),
    respMax_(3.),
    dir_(outDir)
{
  // Override possible existing root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"RECREATE");
  rootfile.Close();
  setGStyle();
  rand_ = new TRandom3(0);
}

ControlPlotsJetSmearing::~ControlPlotsJetSmearing() {
  delete rand_;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::makePlots() const {
  if( config_->read<bool>("create dijet plots",false) )
    plotDijets();
  if( config_->read<bool>("create logP plots",false) )
    plotLogP();
  if( config_->read<bool>("create mean response and resolution plots",false) )
    plotMeanResponseAndResolution();
  if( config_->read<bool>("create parameter error plots",false) )
    plotParameters();
  if( config_->read<bool>("create parameter scan plots",false) )
    plotParameterScan();
  if( config_->read<bool>("create response plots",false) )
    plotResponse();
}



//!  \brief Draw response control plots for events
//!         of type \p SmearDiJet
// --------------------------------------------------
void ControlPlotsJetSmearing::plotResponse() const
{
  std::cout << "Creating response control plots\n";

  // --- Create histograms of response and spectrum ---------------------
  std::string param = config_->read<std::string>("Parametrization Class","");
  double tMin = config_->read<double>("DiJet integration min",0.);
  double tMax = config_->read<double>("DiJet integration max",1.);
  double rMin = config_->read<double>("Response pdf min",0.);
  double rMax = config_->read<double>("Response pdf max",2.);

  // Create pt binning for response function evaluation
  std::string binning = config_->read<std::string>("plots pt binning","");

  std::string binningVar = "ptGen";
  if( binning.find("ptDijet") != std::string::npos ) binningVar = "ptDijet";

  std::vector<double> ptBinEdges(2,0.);
  if( binning.find("binning") != std::string::npos ) {
    ptBinEdges.clear();
    ptBinEdges = bag_of<double>(config_->read<std::string>("plots pt bin edges","100 500"));
  } else if( binning.find("cuts") != std::string::npos ) {
    if( binningVar == "ptGen" ) {
      ptBinEdges.clear();
      ptBinEdges.push_back(config_->read<double>("Et genJet min",100));
      ptBinEdges.push_back(config_->read<double>("Et genJet max",500));
    } else if( binningVar == "ptDijet" ) {
      ptBinEdges.clear();
      ptBinEdges.push_back(config_->read<double>("Et min cut on dijet",100));
      ptBinEdges.push_back(config_->read<double>("Et max cut on dijet",500));
    }
  }
  int nPtBins = static_cast<int>(ptBinEdges.size()-1);
  std::vector<double> ptBinCenters(nPtBins);
  for(int i = 0; i < nPtBins; i++) {
    ptBinCenters[i] = 0.5 * ( ptBinEdges[i] + ptBinEdges[i+1] );
  }

  std::vector<TH1F*> hRespMeasAbs(nPtBins);   // The measured response ptJet / ptGen absolute entries
  std::vector<TH1F*> hRespMeas(nPtBins);      // The measured response ptJet / ptGen
  std::vector<TH1F*> hRespMCPtHat(nPtBins);      // The measured response ptJet / ptGen
  std::vector<TH1F*> hRespFitStart(nPtBins);  // The response pdf with start values
  std::vector<TH1F*> hRespFit(nPtBins);       // The fitted response pdf
  std::vector<TH1F*> hRespFitErrStat(nPtBins);// The fitted response pdf with fitted errors
  std::vector<TH1F*> hRespFitStep(nPtBins);   // Step function part of the response pdf
  std::vector<TH1F*> hRespFitGaus(nPtBins);   // Gauss part of the response pdf
  std::vector<TH1F*> hRespFitSum(nPtBins);    // Sum of step and Gauss part
  std::vector<TH1F*> hRatio(nPtBins);
  std::vector<TH1F*> hAbsResp(nPtBins);
  std::vector<TH1F*> hAbsRespFit(nPtBins);
  std::vector<TH1F*> hPtGenAbsBins(nPtBins);
  std::vector<TH1F*> hPtGenAsym(nPtBins);
  std::vector<TH1F*> hPtAsym(nPtBins);
  std::vector<TH1F*> hFitPtAsym(nPtBins);

  TH1F * hTruthPDF = 0;      // Truth pdf
  TH1F * hTruthPDFErrStat = 0;      // Truth pdf
  TH1F * hPtGenAbs = 0;         // PtGen spectrum
  TH1F * hPtGen = 0;         // PtGen spectrum
  TH1F * hPtGenJet1 = 0;         // PtGen spectrum
  TH1F * hPtHat = 0;         // PtHat spectrum
  TH1F * hPtDijet = 0;       // Dijet spectrum


  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    std::string name = "hRespMeasAbs_" + toString(ptBin);
    hRespMeasAbs[ptBin] = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{gen}_{T};dN / dR",
				   respNBins_,respMin_,respMax_);
    hRespMeasAbs[ptBin]->SetLineWidth(2);

    name = "hRespMeas_" + toString(ptBin);
    hRespMeas[ptBin] = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{gen}_{T};1 / N  dN / dR",
				respNBins_,respMin_,respMax_);
    hRespMeas[ptBin]->Sumw2();
    hRespMeas[ptBin]->SetLineWidth(2);

    name = "hRespMCPtHat_" + toString(ptBin);
    hRespMCPtHat[ptBin] = static_cast<TH1F*>(hRespMeas[ptBin]->Clone(name.c_str()));
    hRespMCPtHat[ptBin]->SetTitle(";R = p^{jet}_{T} / #hat{p}_{T};1 / N  dN / dR");

    name = "hRespFit_" + toString(ptBin);
    hRespFit[ptBin] = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1 / N  dN / dR",
			       5*respNBins_,respMin_,respMax_);
    hRespFit[ptBin]->SetLineColor(2);
    hRespFit[ptBin]->SetLineWidth(2);

    name = "hRespFitErrStat_" + toString(ptBin);
    hRespFitErrStat[ptBin] = static_cast<TH1F*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitErrStat[ptBin]->SetFillColor(45);

    name = "hRespFitStart_" + toString(ptBin);
    hRespFitStart[ptBin] = static_cast<TH1F*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitStart[ptBin]->SetLineStyle(2);

    name = "hRespFitStep_" + toString(ptBin);
    hRespFitStep[ptBin] = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1 / N  dN / dR",
				  config_->read<int>("Response pdf nsteps",10),
				  config_->read<double>("Response pdf min",0.),
				  config_->read<double>("Response pdf max",1.8));
    hRespFitStep[ptBin]->Sumw2();
    hRespFitStep[ptBin]->SetLineColor(9);
    hRespFitStep[ptBin]->SetLineWidth(2);

    name = "hRespFitGaus_" + toString(ptBin);
    hRespFitGaus[ptBin] = static_cast<TH1F*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitGaus[ptBin]->SetLineColor(8);

    name = "hRespFitSum_" + toString(ptBin);
    hRespFitSum[ptBin] = static_cast<TH1F*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitSum[ptBin]->SetLineColor(1);

    name = "hRatio_" + toString(ptBin);
    hRatio[ptBin] = new TH1F(name.c_str(),";Prediction / truth",
				  respNBins_,respMin_,respMax_);
    hRatio[ptBin]->SetLineWidth(2);

    name = "hPtGenAbs_" + toString(ptBin);
    hPtGenAbsBins[ptBin] = new TH1F(name.c_str(),";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)"
				    ,25,0.8*ptBinEdges[ptBin],1.1*ptBinEdges[ptBin+1]);
    hPtGenAbsBins[ptBin]->GetXaxis()->SetNdivisions(505);
    hPtGenAbsBins[ptBin]->SetLineWidth(2);

    name = "hAbsResp_" + toString(ptBin);
    double min = ptBinCenters[ptBin] - 5*1.25*sqrt(ptBinCenters[ptBin]);
    double max = ptBinCenters[ptBin] + 5*1.25*sqrt(ptBinCenters[ptBin]);
    hAbsResp[ptBin] = new TH1F(name.c_str(),";p_{T} (GeV);1 / N  dN / p_{T}  1 / (GeV)"
			       ,respNBins_,min,max);
    hAbsResp[ptBin]->Sumw2();
    hAbsResp[ptBin]->SetLineWidth(2);

    name = "hAbsRespFit_" + toString(ptBin);
    hAbsRespFit[ptBin] = new TH1F(name.c_str(),
				  ";p_{T} (GeV);1 / N  dN / p_{T}  1 / (GeV)"
				  ,5*respNBins_,min,max);
    hAbsRespFit[ptBin]->SetLineColor(2);
    hAbsRespFit[ptBin]->SetLineWidth(2);

    name = "hPtGenAsym_" + toString(ptBin);
    hPtGenAsym[ptBin] = new TH1F(name.c_str(),
				 ";p^{gen}_{T} asymmetry;",
				 30,-0.3,0.3);
    hPtGenAsym[ptBin]->Sumw2();
    hPtGenAsym[ptBin]->SetLineWidth(2);

    name = "hPtAsym_" + toString(ptBin);
    hPtAsym[ptBin] = new TH1F(name.c_str(),
			      ";p_{T} asymmetry;",
			      30,-0.4,0.4);
    hPtAsym[ptBin]->Sumw2();
    hPtAsym[ptBin]->SetLineWidth(2);

    name = "hFitPtAsym_" + toString(ptBin);
    hFitPtAsym[ptBin] = new TH1F(name.c_str(),
				 ";p_{T} asymmetry;",
				 5*respNBins_,-0.4,0.4);
    hFitPtAsym[ptBin]->Sumw2();
    hFitPtAsym[ptBin]->SetLineWidth(2);
    hFitPtAsym[ptBin]->SetLineColor(2);
  }

  hPtGenAbs = new TH1F("hPtGenAbs",";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
		       40,ptBinEdges.front(),ptBinEdges.back());
  hPtGenAbs->GetXaxis()->SetNdivisions(505);
  hPtGenAbs->SetLineWidth(2);

  hPtGen = static_cast<TH1F*>(hPtGenAbs->Clone("hPtGen"));
  hPtGen->SetTitle(";p^{gen}_{T} (GeV);1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
  hPtGen->Sumw2();

  hPtGenJet1 = static_cast<TH1F*>(hPtGen->Clone("hPtGenJet1"));
  hPtGenJet1->SetTitle(";Jet 1 p^{gen}_{T} (GeV);1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
  
  hPtHat = static_cast<TH1F*>(hPtGen->Clone("hPtHat"));
  hPtHat->SetTitle(";#hat{p}_{T} (GeV);1 / N  dN / d#hat{p}_{T}  1 / (GeV)");

  hPtDijet = static_cast<TH1F*>(hPtGen->Clone("hPtDijet"));
  hPtDijet->SetTitle(";p^{dijet}_{T} (GeV);1 / N  dN / dp^{dijet}_{T}  1 / (GeV)");

  hTruthPDF = new TH1F("hTruthPDF",";p^{true}_{T} (GeV);1 / N  dN / dp^{true}_{T}  1 /  (GeV)",
		       5*hPtGen->GetNbinsX(),tMin,tMax);
  hTruthPDF->SetLineColor(2);
  hTruthPDF->SetLineWidth(2);

  hTruthPDFErrStat = static_cast<TH1F*>(hTruthPDF->Clone("hTruthPDFErrStat"));
  hTruthPDFErrStat->SetFillColor(45);


  // --- Fill histograms of measured response --------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      hPtHat->Fill( dijet->ptHat(), dijet->GetWeight() );

      const Jet *j1 = dijet->jet1();
      const Jet *j2 = dijet->jet2();
      if( rand_->Uniform() > 0.5 ) {
	j1 = dijet->jet2();
	j2 = dijet->jet1();
      }
      double ptGenAsym = j1->genPt() - j2->genPt();
      ptGenAsym /= j1->genPt() + j2->genPt();
      double ptAsym = j1->pt() - j2->pt();
      ptAsym /= j1->pt() + j2->pt();

      for(int i = 0; i < 2; i++) {        // Loop over both jets
	const Jet * jet = dijet->jet1();
	if( i == 1 ) jet = dijet->jet2();

	hPtGenAbs->Fill( jet->genPt(), dijet->GetWeight() );
	hPtGen->Fill( jet->genPt(), dijet->GetWeight() );
	if( i == 0 ) hPtGenJet1->Fill( jet->genPt(), dijet->GetWeight() );

	for(int bin = 0; bin < nPtBins; bin++) {
	  double var = 0.;
	  if( binningVar == "ptGen" ) var = jet->genPt();
	  else if( binningVar == "ptDijet" ) var = dijet->dijetPt();
	  if( ptBinEdges[bin] <= var && var < ptBinEdges[bin+1] ) {
	    hRespMeasAbs[bin]->Fill( jet->pt() / jet->genPt(), dijet->GetWeight() );
	    hRespMeas[bin]->Fill( jet->pt() / jet->genPt(), dijet->GetWeight() );
	    hAbsResp[bin]->Fill( jet->pt(), dijet->GetWeight() );
	    hPtGenAbsBins[bin]->Fill( jet->genPt(), dijet->GetWeight() );
	    hPtGenAsym[bin]->Fill( ptGenAsym, dijet->GetWeight() );
	    hPtAsym[bin]->Fill( ptAsym, dijet->GetWeight() );
	    hRespMCPtHat[bin]->Fill( jet->pt()/dijet->ptHat(), dijet->GetWeight() );
	    continue;
	  }
	}
      }
      hPtDijet->Fill( dijet->dijetPt(), dijet->GetWeight() );
      //     } else if( (*datait)->GetType() == TypeSmearPhotonJet )  {  // Select GammaJet events
      //       SmearPhotonJet * gammaJet = static_cast<SmearPhotonJet*>(*datait);  
      
      //       hPtHat->Fill( gammaJet->ptHat(), gammaJet->GetWeight() );
      //       const Jet * jet = gammaJet->jet();
      
      //       hPtGenAbs->Fill( jet->genPt(), gammaJet->GetWeight() );
      //       //hPtGen->Fill( jet->genPt(), gammaJet->GetWeight() );
      
      //       for(int i = 0; i < nPtBins; i++) {
      // 	double var = 0.;
      // 	if( binningVar == "ptGen" ) var = jet->genPt();
      // 	if( ptBinEdges[i] <= var && var < ptBinEdges[i+1] ) {
      // 	  hRespMeasAbs[i]->Fill( jet->pt() / jet->genPt(), gammaJet->GetWeight() );
      // 	  hRespMeas[i]->Fill( jet->pt() / jet->genPt(), gammaJet->GetWeight() );
      // 	  hPtGenAbsBins[i]->Fill( jet->genPt(), gammaJet->GetWeight() );
      // 	  continue;
      // 	}
      //       }
    }
  } // End of loop over data
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    normHist(hRespMeas[ptBin],"width");
    normHist(hRespMCPtHat[ptBin],"width");
    normHist(hAbsResp[ptBin],"width");
    normHist(hPtGenAsym[ptBin],"width");
    normHist(hPtAsym[ptBin],"width");
  }
  normHist(hPtGen,tMin,tMax,"width");
  normHist(hPtGenJet1,tMin,tMax,"width");
  normHist(hPtHat,tMin,tMax,"width");
  normHist(hPtDijet,tMin,tMax,"width");


  // --- Fill histograms of fitted response ----------------
  // Get parameters
  std::vector<double> fittedPar(param_->GetNumberOfParameters());
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    fittedPar.at(i) = param_->GetPars()[i];
  }

  std::vector<double> scale = bag_of<double>(config_->read<string>("jet parameter scales",""));
  std::vector<double> startParJet = bag_of<double>(config_->read<string>("jet start values",""));
  std::vector<double> auxPar = bag_of<double>(config_->read<string>("mean response parameters","1 0"));
  SmearData * smearData = dynamic_cast<SmearData*>(data_->front());
  if( smearData ) {
    // Loop over ptBins
    for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
      // Interpolated response function
      for(int bin = 1; bin <= hRespFit[ptBin]->GetNbinsX(); bin++) {
	double r = hRespFit[ptBin]->GetBinCenter(bin);
	double pt = hPtGenAbsBins[ptBin]->GetMean();
	double val = smearData->pdfResp(r,pt);
	hRespFit[ptBin]->SetBinContent(bin,val);
	hRespFitErrStat[ptBin]->SetBinContent(bin,val);
	hRespFitErrStat[ptBin]->SetBinError(bin,smearData->pdfRespError(r,pt));

//   	if( ptBin == 0 ) {
//   	  std::cout << r << ":  " << smearData->respPDF(r,ptBinCenters[ptBin]) << " +/- " << smearData->respPDFError(r,ptBinCenters[ptBin]) << std::endl;
//   	}
      }
      for(int bin = 1; bin <= hFitPtAsym[ptBin]->GetNbinsX(); bin++) {
	double a = hFitPtAsym[ptBin]->GetBinCenter(bin);
	double pt = hPtGenAbsBins[ptBin]->GetMean();
	hFitPtAsym[ptBin]->SetBinContent(bin,smearData->pdfDijetAsym(a,pt));
      }
      // Interpolated pdf of measured jet pt
      for(int bin = 1; bin <= hAbsRespFit[ptBin]->GetNbinsX(); bin++) {
	double ptMeas = hAbsRespFit[ptBin]->GetBinCenter(bin);
	double ptTrue = hPtGenAbsBins[ptBin]->GetMean();
	hAbsRespFit[ptBin]->SetBinContent(bin,smearData->pdfPtMeas(ptMeas,ptTrue));
      }

      // Ratio plot
      for(int bin = 1; bin <= hRatio[ptBin]->GetNbinsX(); bin++) {
	double rMin = hRatio[ptBin]->GetBinLowEdge(bin);
	double rMax = rMin + hRatio[ptBin]->GetBinWidth(bin);
	int min = hRespFit[ptBin]->FindBin(rMin);
	int max = hRespFit[ptBin]->FindBin(rMax);
	double pred = (hRespFit[ptBin]->Integral(min,max))/(1+max-min);
	double truth = hRespMeas[ptBin]->GetBinContent(bin);
	double ratio = 0.;
	double error = 0.;
	if( truth ) {
	  ratio = pred / truth;
	  error = pred / truth / truth * hRespMeas[ptBin]->GetBinError(bin);
	}
	hRatio[ptBin]->SetBinContent(bin,ratio);
	hRatio[ptBin]->SetBinError(bin,error);
      }

      // Interpolated fit function with start values
      // Copy start values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = startParJet.at(i);
      }
      // Plot response function
      for(int bin = 1; bin <= hRespFitStart[ptBin]->GetNbinsX(); bin++) {
	double pt = hPtGenAbsBins[ptBin]->GetMean();
	double r = hRespFitStart[ptBin]->GetBinCenter(bin);
	hRespFitStart[ptBin]->SetBinContent(bin,smearData->pdfResp(r,pt));
      }
      // Copy back fitted values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = fittedPar.at(i);
      }

      // In case of interpolated step + gauss parametrization
      //      if( param == "SmearParametrizationStepGaussInter" ) {
// 	// Step part of fit function
// 	for(int bin = 1; bin <= hRespFitStep[ptBin]->GetNbinsX(); bin++) {
// 	  double val  = scale.at(bin+1)*(smearData->respPar(bin+1));
// 	  hRespFitStep[ptBin]->SetBinContent(bin,val);
// 	}
// 	normHist(hRespFitStep[ptBin],"width");
// 	hRespFitStep[ptBin]->Scale(1. - scale.at(0)*(smearData->respPar(0)));
	
// 	// Gauss part of fit function
// 	for(int bin = 1; bin <= hRespFitGaus[ptBin]->GetNbinsX(); bin++) {
// 	  // Mean
// 	  double mu = auxPar.at(0);
// 	  // Width
// // 	  double a1 = scale.at(1)*(smearData->respPar()[1]);
// // 	  double a2 = scale.at(2)*(smearData->respPar()[2]);
// // 	  double a3 = scale.at(3)*(smearData->respPar()[3]);
// // 	  double sigma = sqrt( a1*a1/ptBinCenters[ptBin]/ptBinCenters[ptBin]
// // 			       + a2*a2/ptBinCenters[ptBin] + a3*a3 );
// 	  double sigma = scale[1]*(smearData->respPar(1));
// 	  // pdf
// 	  double c     = scale.at(0)*(smearData->respPar(0));
// 	  double r     = hRespFitGaus[ptBin]->GetBinCenter(bin);
// 	  double val   = c * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
// 	  hRespFitGaus[ptBin]->SetBinContent(bin,val);
// 	}
      
// 	// Sum
// 	for(int binGaus = 1; binGaus <= hRespFitGaus[ptBin]->GetNbinsX(); binGaus++) {
// 	  int binStep = hRespFitStep[ptBin]->FindBin(hRespFitGaus[ptBin]->GetBinCenter(binGaus));
// 	  double val = hRespFitStep[ptBin]->GetBinContent(binStep)
// 	    + hRespFitGaus[ptBin]->GetBinContent(binGaus);
// 	  hRespFitSum[ptBin]->SetBinContent(binGaus,val);
// 	}
//       } else if( param == "SmearParametrizationCrystalBall"
// 		 ||  param == "SmearParametrizationCrystalBallPt"
// 		 ||  param == "SmearParametrizationGauss") {
//       } else {
// 	std::cout << "WARNING: No controlplots implemented for parametrization '" << param << "'\n";
//       }
    } // End of loop over ptBins
  } // End if( smearData )


  // --- Fill histograms of fitted truth spectrum -----------

  // Fill histogram of assumed dijet truth pdf
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    if( (*datait)->GetType() == TypeSmearDiJet ) {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      
      for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
	double t = hTruthPDF->GetBinCenter(bin);
	hTruthPDF->SetBinContent(bin,dijet->pdfPtTrue(t));
      }
      for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
	hTruthPDFErrStat->SetBinContent(bin,hTruthPDF->GetBinContent(bin));
	double t = hTruthPDFErrStat->GetBinCenter(bin);
	hTruthPDFErrStat->SetBinError(bin,dijet->pdfPtTrueError(t));
      }
      break;
    }
  }


  // --- Set x-axis ranges ----------------------------------
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    hRespMeasAbs[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);   
    hRespMeas[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);      
    hRespMCPtHat[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);  
    hRespFitStart[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);  
    hRespFit[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);       
    hRespFitErrStat[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);
    hRespFitStep[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);   
    hRespFitGaus[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);   
    hRespFitSum[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);    
    hRatio[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);
  }
  

  // --- Set y-axis ranges ----------------------------------
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    double minFac = 0.;
    double maxFac = 1.6;
    
    double min = 0.;
    double max = 0.;
    findYRange(hRespMeas[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hRespMeas[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hRespMCPtHat[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hRespMCPtHat[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hAbsResp[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hAbsResp[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hPtAsym[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hPtAsym[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hPtGenAsym[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hPtGenAsym[ptBin]->GetYaxis()->SetRangeUser(min,max);

    setYRange(hRespMeasAbs[ptBin],0.5,50.);
  }
  setYRange(hPtDijet, 0.5, 100.);
  setYRange(hPtGen, 0.5, 100.);
  setYRange(hPtGenJet1, 0.5, 100.);
  setYRange(hPtHat, 0.5, 100.);


  // --- Plot histograms -----------------------------------
  // Label bins
  std::vector<TLegend*> legPtRange(nPtBins);
  std::vector<TLegend*> legPtRangeAndCenters(nPtBins);
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    legPtRange[ptBin] = new TLegend(0.23,0.72,0.78,0.8);
    legPtRange[ptBin]->SetBorderSize(0);
    legPtRange[ptBin]->SetFillColor(0);
    legPtRange[ptBin]->SetTextFont(42);

    legPtRangeAndCenters[ptBin] = new TLegend(0.23,0.65,0.8,0.8);
    legPtRangeAndCenters[ptBin]->SetBorderSize(0);
    legPtRangeAndCenters[ptBin]->SetFillColor(0);
    legPtRangeAndCenters[ptBin]->SetTextFont(42);

    std::string binVar;
    if( binningVar == "ptDijet" ) binVar = "p^{dijet}_{T}";
    else if( binningVar == "ptGen" ) binVar = "p^{gen}_{T}";

    std::string label = toString(ptBinEdges[ptBin])
      + " < " + binVar + " < "
      + toString(ptBinEdges[ptBin+1])
      + " GeV";
    legPtRange[ptBin]->AddEntry(hRespMeas[ptBin],label.c_str(),"L");
    legPtRangeAndCenters[ptBin]->AddEntry(hRespMeas[ptBin],label.c_str(),"L");
    label = "p_{T} = " + toString(hPtGenAbsBins[ptBin]->GetMean()) + " GeV";
    legPtRangeAndCenters[ptBin]->AddEntry(hRespFit[ptBin],label.c_str(),"L");
  }

  // Write histos to ps file
  TPostScript * const ps = new TPostScript((dir_+"/jsResponse.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Jet Response",0,0,600,600);

  TLegend *legFitStart = new TLegend(0.23,0.5,0.5,0.65);
  legFitStart->SetBorderSize(0);
  legFitStart->SetFillColor(0);
  legFitStart->SetTextFont(42);
  legFitStart->AddEntry(hRespFitStart[0],"At start","L");
  legFitStart->AddEntry(hRespFit[0],"After fit","L");
  
  int logy = 0;
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    // Measured and fitted response
//     ps->NewPage();
//     c1->cd();
//     hRespMeasAbs[ptBin]->Draw();
//     legPtRange[ptBin]->Draw("same");
//     c1->SetLogy();
//     c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hRespMeas[ptBin]->Draw();
//     legPtRange[ptBin]->Draw("same");
//     c1->SetLogy();
//     c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->Draw();
    hRespFitErrStat[ptBin]->Draw("E3 same");
    hRespMeas[ptBin]->Draw("same");
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(logy);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->Draw();
    hRespFit[ptBin]->Draw("Lsame");
    hRespFitStart[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    legPtRangeAndCenters[ptBin]->Draw("same");
    legFitStart->Draw("same");
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->Draw();
    if( param == "SmearParametrizationStepGaussInter" ) {
      hRespFitStep[ptBin]->Draw("same");
      hRespFitGaus[ptBin]->Draw("same");
      //hRespFitSum[ptBin]->Draw("same");
    }
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(logy);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMCPtHat[ptBin]->Draw();
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(logy);
    c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hAbsResp[ptBin]->Draw();
//     //    hAbsRespFit[ptBin]->Draw("Lsame");
//     gPad->RedrawAxis();
//     legPtRangeAndCenters[ptBin]->Draw("same");
//     c1->SetLogy();
//     c1->Draw();

    ps->NewPage();
    c1->cd();
    hPtGenAsym[ptBin]->Draw();
    legPtRange[ptBin]->Draw("same");
    c1->SetLogy(logy);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hPtAsym[ptBin]->Draw();
    hFitPtAsym[ptBin]->Draw("Lsame");
    legPtRangeAndCenters[ptBin]->Draw("same");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hRatio[ptBin]->Draw();
//     c1->SetLogy(0);
//     c1->Draw();
  }

  // Truth spectrum
  delete legPtRange[0];
  legPtRange[0] = new TLegend(0.23,0.72,0.78,0.8);
  legPtRange[0]->SetBorderSize(0);
  legPtRange[0]->SetFillColor(0);
  legPtRange[0]->SetTextFont(42);
  std::string binVar = config_->read<std::string>("plots pt binning","");
  if( binVar.find("ptDijet") != std::string::npos ) binVar = "p^{dijet}_{T}";
  else if( binVar.find("ptGen") != std::string::npos ) binVar = "p^{gen}_{T}";
  else binVar = "p^{gen}_{T}";
  std::string label = toString(ptBinEdges.front())
    + " < " + binVar + " < "
    + toString(ptBinEdges.back())
    + " GeV";
  legPtRange[0]->AddEntry(hPtGen,label.c_str(),"L");

  ps->NewPage();
  c1->cd();
  hPtGenAbs->Draw();
  legPtRange[0]->Draw("same");
  c1->SetLogy();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  hPtGen->Draw();
  //hTruthPDFErrStat->Draw("E3same");
  //hPtGen->Draw("same");
  hTruthPDF->Draw("Lsame");
  gPad->RedrawAxis();
  legPtRange[0]->Draw("same");
  c1->SetLogy();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  double tmpMin = 0.;
  double tmpMax = 0.;
  findYRange(hPtGen,tmpMin,tmpMax);
  tmpMax *= 1.4;
  hPtGen->GetYaxis()->SetRangeUser(0.,tmpMax);
  hPtGen->Draw();
  hTruthPDF->Draw("same");
  gPad->RedrawAxis();
  legPtRange[0]->Draw("same");
  c1->SetLogy(0);
  c1->SetLogx(0);
  c1->Draw();
  c1->SetLogx(0);

  std::vector<TObject*> objs;
  objs.clear();
  objs.push_back(hPtHat);
  objs.push_back(hTruthPDF);
  drawPSPage(ps,c1,objs,"",true);

  objs.clear();
  objs.push_back(hPtDijet);
  objs.push_back(hTruthPDF);
  drawPSPage(ps,c1,objs,"",true);


  // Write histos to root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    rootfile.WriteTObject(hRespMeasAbs[ptBin]);
    rootfile.WriteTObject(hRespMeas[ptBin]);
    rootfile.WriteTObject(hRespMCPtHat[ptBin]);
    rootfile.WriteTObject(hRespFit[ptBin]);
    rootfile.WriteTObject(hAbsResp[ptBin]);
    rootfile.WriteTObject(hRespFitErrStat[ptBin]);
    rootfile.WriteTObject(hRespFitStart[ptBin]);
    rootfile.WriteTObject(hRespFitStep[ptBin]);
    rootfile.WriteTObject(hRespFitGaus[ptBin]);
    rootfile.WriteTObject(hRespFitSum[ptBin]);
    rootfile.WriteTObject(hAbsRespFit[ptBin]);
    rootfile.WriteTObject(hRatio[ptBin]);
    rootfile.WriteTObject(hPtGenAsym[ptBin]);
    rootfile.WriteTObject(hPtAsym[ptBin]);
    rootfile.WriteTObject(hFitPtAsym[ptBin]);
  }
  rootfile.WriteTObject(hPtGenAbs);
  rootfile.WriteTObject(hPtGen);
  rootfile.WriteTObject(hPtGenJet1);
  rootfile.WriteTObject(hPtHat);
  rootfile.WriteTObject(hPtDijet);
  rootfile.WriteTObject(hTruthPDF);
  rootfile.WriteTObject(hTruthPDFErrStat);

  rootfile.Close();


  // --- Clean up ------------------------------------------
  for(int ptBin = 0; ptBin < nPtBins; ptBin++) {
    delete hRespMeasAbs[ptBin];
    delete hRespMeas[ptBin];
    delete hRespMCPtHat[ptBin];
    delete hAbsResp[ptBin];
    delete hRespFit[ptBin];
    delete hRespFitErrStat[ptBin];
    delete hRespFitStart[ptBin];
    delete hRespFitStep[ptBin];
    delete hRespFitGaus[ptBin];
    delete hRespFitSum[ptBin];
    delete legPtRangeAndCenters[ptBin];
    delete hAbsRespFit[ptBin];
    delete hRatio[ptBin];
    delete hPtGenAbsBins[ptBin];
    delete hPtGenAsym[ptBin];
    delete hPtAsym[ptBin];
    delete hFitPtAsym[ptBin];
  }
  delete legFitStart;
  delete hPtGenAbs;
  delete hPtGen;
  delete hPtGenJet1;
  delete hPtHat;
  delete hPtDijet;
  delete hTruthPDF;
  delete hTruthPDFErrStat;
  delete c1;
  delete ps;
}



//!  \brief Draw control plots for events
//!         of type \p SmearDiJet
// --------------------------------------------------
void ControlPlotsJetSmearing::plotDijets() const
{
  std::cout << "Creating dijet control plots\n";

  // --- Create histograms --------------------------

  // Find pt ranges
  double minPtHat    = 10000.;
  double maxPtHat    = 0.;
  double minGenJetPt = 10000.;
  double maxGenJetPt = 0.;
  double minCalJetPt = 10000.;
  double maxCalJetPt = 0.;
  double minDijetPt  = 10000.;
  double maxDijetPt  = 0.;
  double min3rdJetPt = 10000.;
  double max3rdJetPt = 0.;
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      // Loop over both jets
      for(int i = 0; i < 2; i++) {        
	const Jet * jet = dijet->jet1();
	if( i == 1 ) jet = dijet->jet2();

	if( jet->genPt() < minGenJetPt ) minGenJetPt = jet->genPt();
	if( jet->genPt() > maxGenJetPt ) maxGenJetPt = jet->genPt();

	if( jet->pt() < minCalJetPt ) minCalJetPt = jet->pt();
	if( jet->pt() > maxCalJetPt ) maxCalJetPt = jet->pt();
      }

      if( dijet->ptHat() < minPtHat ) minPtHat = dijet->ptHat();
      if( dijet->ptHat() > maxPtHat ) maxPtHat = dijet->ptHat();

      if( dijet->dijetPt() < minDijetPt ) minDijetPt = dijet->dijetPt();
      if( dijet->dijetPt() > maxDijetPt ) maxDijetPt = dijet->dijetPt();

      const Jet * jet3 = dijet->jet3();
      if( jet3->pt() < min3rdJetPt ) min3rdJetPt = jet3->pt();
      if( jet3->pt() > max3rdJetPt ) max3rdJetPt = jet3->pt();
    }
  }

  // Pt distributions
  TH1F * hPtHat = new TH1F("hPtHat",";#hat{p}_{T} (GeV);dN / d#hat{p}_{T}  1 / (GeV)",
			   50,0.9*minPtHat, 1.1*maxPtHat);
  hPtHat->SetLineWidth(2);
  hPtHat->SetNdivisions(505);
  hPtHat->Sumw2();

  std::vector<TH1F*> hGenJetPt;
  std::vector<TH1F*> hCalJetPt;
  std::string genJetNames[3] = { "hGenJetPtBothJets",
				 "hGenJetPtJet1", 
				 "hGenJetPtJet2" };
  std::string calJetNames[3] = { "hCalJetPtBothJets",
				 "hCalJetPtJet1", 
				 "hCalJetPtJet2" };
  int color[3] = { 1, 2, 4 };
  for(int i = 0; i < 3; i++) {
    TH1F * h = 0;

    h = new TH1F(genJetNames[i].c_str(),";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
		 50,0.9*minGenJetPt,1.1*maxGenJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    h->SetNdivisions(505);
    hGenJetPt.push_back(h);

    h = new TH1F(calJetNames[i].c_str(),";p^{jet}_{T} (GeV);dN / dp^{jet}_{T}  1 / (GeV)",
		 50,0.9*minCalJetPt,1.1*maxCalJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    h->SetNdivisions(505);
    hCalJetPt.push_back(h);
  }

  TH2F * hCalJet2vsCalJet1Pt = new TH2F("hCalJet2vsCalJet1Pt",
					";p^{jet1}_{T} (GeV);p^{jet2}_{T} (GeV)",
					50,0.9*minCalJetPt,1.1*maxCalJetPt,
					50,0.9*minCalJetPt,1.1*maxCalJetPt);

  TH1F * hDijetPt = new TH1F("hDijetPt",";p^{dijet}_{T} (GeV);1 / N  dN / dp^{dijet}_{T}  1 / (GeV)",
			     50,0.9*minDijetPt, 1.1*maxDijetPt);
  hDijetPt->SetLineWidth(2);
  hDijetPt->Sumw2();

  TH1F * h3rdJetPt = new TH1F("h3rdJetPt",";p^{jet3}_{T} (GeV);1 / N  dN / dp^{jet3}_{T}  1 / (GeV)",
			     50,0.9*min3rdJetPt, 1.1*max3rdJetPt);
  h3rdJetPt->SetLineWidth(2);
  h3rdJetPt->Sumw2();

  TH2F * h3rdJetvsDijetPt = new TH2F("h3rdJetvsDijetPt",
				     ";p^{dijet}_{T} (GeV);p^{jet3}_{T} (GeV)",
				     50,0.9*minDijetPt, 1.1*maxDijetPt,
				     50,0.9*min3rdJetPt, 1.1*max3rdJetPt);

  TH1F * hRel3rdJetPt = new TH1F("hRel3rdJetPt",
				 ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};1 / N  dN / dp^{jet3}_{T,rel}",
			     50,0,1.4);
  hRel3rdJetPt->SetLineWidth(2);
  hRel3rdJetPt->Sumw2();

  TH1F * hDeltaPhi = new TH1F("hDeltaPhi",";#Delta#phi;1 / N  dN / d#Delta#phi",
			      25,1.5,M_PI);
  hDeltaPhi->SetLineWidth(2);
  hDeltaPhi->Sumw2();

  TH2F * hDeltaPhivsRel3rdJetPt = new TH2F("hDeltaPhivsRel3rdJetPt",
					   ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};#Delta#phi",
					   25,0,1,25,1.5,M_PI);

  // Response correlations
  TH2F * hRvsDeltaPhi = new TH2F("hRvsDeltaPhi",";#Delta#phi;p^{jet}_{T} / p^{gen}_{T}",
				 50,1.5,M_PI,50,0,2);
  hRvsDeltaPhi->SetMarkerStyle(7);
  TH2F * hRvsRel3rdJetPt = new TH2F("hRvsRel3rdJetPt",
				    ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};p^{jet}_{T} / p^{gen}_{T}",
				    50,0,1.4,50,0,2);
  hRvsRel3rdJetPt->SetMarkerStyle(7);
  TH2F * hRvs3rdJetPt = new TH2F("hRvs3rdJetPt",
				 ";p^{jet3} (GeV);p^{jet}_{T} / p^{gen}_{T}",
				 50,0.9*min3rdJetPt, 1.1*max3rdJetPt,50,0,2);
  hRvs3rdJetPt->SetMarkerStyle(7);
  TH2F * hRvsEMF = new TH2F("hRvsEMF",";EMF;p^{jet}_{T} / p^{gen}_{T}",
			    50,0,1,50,0,2);
  hRvsEMF->SetMarkerStyle(7);
  TH2F * hRvsDeltaR = new TH2F("hRvsDeltaR",";#Delta R(jet,genJet);p^{jet}_{T} / p^{gen}_{T}",
			       25,0,0.4,50,0,2);
  hRvsDeltaR->SetMarkerStyle(7);
  hRvsDeltaR->GetXaxis()->SetNdivisions(505);
    



  // --- Fill histograms ----------------------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    if( (*datait)->GetType() == TypeSmearDiJet )  { // Select DiJet events
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      const Jet * jet1 = dijet->jet1();
      const Jet * jet2 = dijet->jet2();
      const Jet * jet3 = dijet->jet3();

      double weight = dijet->GetWeight();

      double dPhi = std::abs(TVector2::Phi_mpi_pi( jet1->phi() - jet2->phi() ));
      hDeltaPhi->Fill( dPhi, weight );

      hRvsDeltaPhi->Fill( dPhi, jet1->pt() / jet1->genPt(), weight );
      hRvsDeltaPhi->Fill( dPhi, jet2->pt() / jet2->genPt(), weight );

      hRvsEMF->Fill( jet1->EmEt() / jet1->pt(), jet1->pt() / jet1->genPt(), weight );
      hRvsEMF->Fill( jet2->EmEt() / jet2->pt(), jet2->pt() / jet2->genPt(), weight );

      hRvsRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), jet1->pt() / jet1->genPt(), weight );
      hRvsRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), jet2->pt() / jet2->genPt(), weight );

      hRvs3rdJetPt->Fill( jet3->pt(), jet1->pt() / jet1->genPt(), weight );
      hRvs3rdJetPt->Fill( jet3->pt(), jet2->pt() / jet2->genPt(), weight );

      hRvsDeltaR->Fill( jet1->dR(), jet1->pt() / jet1->genPt(), weight );
      hRvsDeltaR->Fill( jet2->dR(), jet2->pt() / jet2->genPt(), weight );

      hPtHat->Fill( dijet->ptHat(), weight );

      hGenJetPt.at(0)->Fill( jet1->genPt(), weight );
      hGenJetPt.at(0)->Fill( jet2->genPt(), weight );
      hGenJetPt.at(1)->Fill( jet1->genPt(), weight );
      hGenJetPt.at(2)->Fill( jet2->genPt(), weight );

      hCalJetPt.at(0)->Fill( jet1->pt(), weight );
      hCalJetPt.at(0)->Fill( jet2->pt(), weight );
      hCalJetPt.at(1)->Fill( jet1->pt(), weight );
      hCalJetPt.at(2)->Fill( jet2->pt(), weight );

      hCalJet2vsCalJet1Pt->Fill( jet1->pt(), jet2->pt(), weight );
      hDijetPt->Fill( dijet->dijetPt(), weight );

      h3rdJetPt->Fill( jet3->pt(), weight );
      h3rdJetvsDijetPt->Fill( dijet->dijetPt(), jet3->pt(), weight );
      hRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), weight );
      hDeltaPhivsRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), dPhi, weight );
    }
  }


  // Normalizing histograms
//   for(size_t i = 0; i < hGenJetPt.size(); i++) {
//     normHist( hGenJetPt.at(i) );
//     normHist( hCalJetPt.at(i) );
//   }
//  normHist( hPtHat );
//  normHist( hDijetPt );
//  normHist( h3rdJetPt );
//  normHist( hRel3rdJetPt );
//  normHist( hDeltaPhi );


  // --- Plot histograms ----------------------------
  TPostScript * const ps = new TPostScript((dir_+"/jsDijets.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","Dijets",0,0,600,600);

  // Applied cuts
  TPaveText * appCuts = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  appCuts->SetFillColor(0);
  appCuts->SetTextFont(42);
  //  appCuts->SetTextSize(0.1);
  appCuts->SetTextAlign(12);
  appCuts->SetBorderSize(0);
  appCuts->AddText(("E^{jet}_{T} > "
		    + toString(config_->read<double>("Et cut on jet",0.)) ).c_str());
  appCuts->AddText(("|#eta| < "
		    + toString(config_->read<double>("Eta cut on jet",0.)) ).c_str());
  appCuts->AddText((toString(config_->read<double>("Et min cut on dijet",0.))
		    + " < E^{dijet}_{T} < "
		    + toString(config_->read<double>("Et max cut on dijet",0.))).c_str());
  appCuts->AddText(("E^{jet3}_{T} < "
		    + toString(config_->read<double>("Et cut on n+1 Jet",0.))
		    + " or E^{jet3}_{T} / E^{dijet}_{T} < "
		    + toString(config_->read<double>("Relative n+1 Jet Et Cut",0.)) ).c_str());
  appCuts->AddText(("#Delta#phi > "
		    + toString(config_->read<double>("Min Delta Phi",0.)) ).c_str());
  appCuts->AddText((toString(config_->read<double>("Min had fraction",0.))
		    + " < E^{had}_{T} / E^{jet}_{T} < "
		    + toString(config_->read<double>("Max had fraction")) ).c_str());
  appCuts->AddText((toString(config_->read<double>("Et genJet min",0.))
		    + " < E^{gen}_{T} < "
		    + toString(config_->read<double>("Et genJet max",0.)) ).c_str());
  appCuts->AddText(("#Delta R(jet,genJet) < "
		    + toString(config_->read<double>("DeltaR cut on jet matching")) ).c_str());
  drawPSPage(ps,c1,appCuts,"",false);

  drawPSPage(ps,c1,hPtHat,"",true);
  
  TLegend * leg = new TLegend(0.7,0.75,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hGenJetPt.at(1),"Jet 1","L");
  leg->AddEntry(hGenJetPt.at(2),"Jet 2","L");

  std::vector<TObject*> objs;
  for(size_t i = 1; i < hGenJetPt.size(); i++) {
    objs.push_back(hGenJetPt.at(i));
  }
  objs.push_back(leg);
  drawPSPage(ps,c1,hGenJetPt.at(0),"",true);
  drawPSPage(ps,c1,objs,"",true);

  objs.clear();
  for(size_t i = 1; i < hCalJetPt.size(); i++) {
    objs.push_back(hCalJetPt.at(i));
  }
  objs.push_back(leg);
  drawPSPage(ps,c1,hCalJetPt.at(0),"",true);
  drawPSPage(ps,c1,objs,"",true);

  drawPSPage(ps,c1,hCalJet2vsCalJet1Pt);
  drawPSPage(ps,c1,hDijetPt,"",true);
  drawPSPage(ps,c1,h3rdJetPt,"",true);
  drawPSPage(ps,c1,h3rdJetvsDijetPt);
  drawPSPage(ps,c1,hRel3rdJetPt,"",true);
  drawPSPage(ps,c1,hDeltaPhi,"",true);
  drawPSPage(ps,c1,hDeltaPhivsRel3rdJetPt);
  drawPSPage(ps,c1,hRvsRel3rdJetPt,"COLZ",true);
  drawPSPage(ps,c1,hRvs3rdJetPt,"COLZ",true);
  drawPSPage(ps,c1,hRvsDeltaPhi,"COLZ",true);
  drawPSPage(ps,c1,hRvsEMF,"COLZ",true);
  drawPSPage(ps,c1,hRvsDeltaR,"COLZ",true);

  ps->Close();

  // Clean up
  for(size_t i = 0; i < hGenJetPt.size(); i++) {
    delete hGenJetPt.at(i);
    delete hCalJetPt.at(i);
  }
  delete appCuts;
  delete leg;
  delete hCalJet2vsCalJet1Pt;
  delete hPtHat;
  delete hDijetPt;
  delete h3rdJetPt;
  delete h3rdJetvsDijetPt;
  delete hRel3rdJetPt;
  delete hDeltaPhi;
  delete hDeltaPhivsRel3rdJetPt;
  delete hRvsRel3rdJetPt;
  delete hRvs3rdJetPt;
  delete hRvsDeltaPhi;
  delete hRvsEMF;
  delete hRvsDeltaR;
  delete c1;
  delete ps;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plotParameters() const {
  std::cout << "Creating parameter control plots\n";

  // ----- Quantities defining the parametrization -----
  std::vector<double> scale = bag_of<double>(config_->read<string>("jet parameter scales",""));
  int nPar = param_->GetNumberOfParameters();


  // ----- Create histograms -----
  // Fitted parameter values with errors
  TH1D *hPars = new TH1D("hParameters",";Parameter index;Parameter value",nPar,-0.5,nPar-0.5);
  hPars->SetMarkerStyle(20);
  hPars->SetNdivisions(nPar);
  // Fitted absolute parameter values (par*scale) with errors
  TH1D *hAbsPars = static_cast<TH1D*>(hPars->Clone("hAbsoluteParameters"));
  hAbsPars->SetYTitle("Absolute parameter value");
  // Relative parameter errors
  TH1D *hRelParErrors = static_cast<TH1D*>(hPars->Clone("hRelativeParameterErrors"));
  hRelParErrors->SetYTitle("Relative parameter error");
  // Global parameter correlations
  TH1D *hGlobalCorr = static_cast<TH1D*>(hPars->Clone("hGlobalParameterCorrelations"));
  hGlobalCorr->SetYTitle("Global parameter correlation");
  // Parameter correlations
  TH2D *hParCorr = new TH2D("hParameterCorrelations",";Parameter index;Parameter index",
			nPar,-0.5,nPar-0.5,nPar,-0.5,nPar-0.5);
  hParCorr->SetNdivisions(nPar,"XY");


  // ----- Fill histograms -----
  for(int i = 0; i < nPar; i++) {
    int bin = 1+i;

    hPars->SetBinContent(bin,param_->GetPars()[i]);
    hPars->SetBinError(bin,param_->GetErrors()[i]);

    hAbsPars->SetBinContent(bin,scale[i]*param_->GetPars()[i]);
    hAbsPars->SetBinError(bin,scale[i]*param_->GetErrors()[i]);
  
    hRelParErrors->SetBinContent(bin,param_->GetErrors()[i]/param_->GetPars()[i]);

    hGlobalCorr->SetBinContent(bin,param_->GetGlobalCorrCoeff()[i]);
  }
  for(int i = 0; i < nPar; i++) {
    for(int j = 0; j < i+1; j++) {
      int idx = (i*i + i)/2 + j;
      double corr = 0.;
      if( param_->GetErrors()[i] && param_->GetErrors()[j] ) {
	corr = param_->GetCovCoeff()[idx] / param_->GetErrors()[i] / param_->GetErrors()[j];
      }
      hParCorr->Fill(i,j,corr);
      if( i != j ) {
	hParCorr->Fill(j,i,corr);
      }
    }
  }


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsParameters.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","Parameters",0,0,600,600);

  drawPSPage(ps,c1,hPars,"PE1");
  drawPSPage(ps,c1,hAbsPars,"PE1");
  drawPSPage(ps,c1,hRelParErrors,"P");
  drawPSPage(ps,c1,hGlobalCorr,"P");
  drawPSPage(ps,c1,hParCorr,"COLZ");

  ps->Close();

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hPars);
  rootfile.WriteTObject(hAbsPars);
  rootfile.WriteTObject(hRelParErrors);
  rootfile.WriteTObject(hGlobalCorr);
  rootfile.WriteTObject(hParCorr);


  // ----- Clean up -----
  delete hPars;
  delete hAbsPars;
  delete hRelParErrors;
  delete hGlobalCorr;
  delete hParCorr;
  delete c1;
  delete ps;
}



//! For each pair (i,j) of free parameters in the fit the
//! likelihood is plotted in the (i,j) plane. The parameters
//! are varied in 3 steps of 2*sigma (error from the fit)
//! below and above the fitted parameter value.
//!
//! The plots are written to the files "jsResponse.root"
//! and "jsParameterScans.ps".
// --------------------------------------------------
void ControlPlotsJetSmearing::plotParameterScan() const {
  std::cout << "Creating parameter scan control plots\n";

  // ----- Set up quantities -----
  int n = 4;
  int nSteps = 2*n+1;
  double nSigma = 2;
  int nPar = param_->GetNumberOfParameters();
  std::vector<TH2D*> hParScans2D;
  std::vector<TLine*> lines;

  // Store likelihood for original parameter values
  double offset = 0.;
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    if( (*dataIt)->GetType() == TypeSmearDiJet )  { // Select DiJet events
      offset += (*dataIt)->chi2();
    }
  }     

  // ----- Vary parameters and calculate likelihood -----
  // Store original parameter values
  std::vector<double> origPars(nPar);
  for(int i = 0; i < nPar; i++) {
    origPars[i] = param_->GetPars()[i];
  }
  // Outer loop over parameters
  for(int i = 0; i < nPar; i++) {
    if( param_->isFixedPar(i) ) continue;
    double idVal = nSigma*param_->GetErrors()[i];
    if( idVal == 0 ) idVal = 0.1;
			
    // Inner loop over parameters
    for(int j = 0; j < i; j++) {
      if( param_->isFixedPar(j) ) continue;
      double jdVal = nSigma*param_->GetErrors()[j];
      if( jdVal == 0 ) jdVal = 0.046;
      // Create histogram of likelihood from i vs j
      TString name = "hParScan2D";
      name += hParScans2D.size();
      TString title = "- #Deltaln(L);Parameter ";
      title += i;
      if( param_->parName(i) != "" ) title += " (" + param_->parName(i) + ")";
      title += ";Parameter ";
      title += j;
      if( param_->parName(j) != "" ) title += " (" + param_->parName(j) + ")";
      hParScans2D.push_back(new TH2D(name,title,
				     nSteps,origPars[i]-(n+0.5)*idVal,origPars[i]+(n+0.5)*idVal,
				     nSteps,origPars[j]-(n+0.5)*jdVal,origPars[j]+(n+0.5)*jdVal));
			    
      // Vary parameters i and j
      for(int is = 0; is < nSteps; is++) {
	double iPar = origPars[i] + (is-n)*idVal; 
	param_->GetPars()[i] = iPar;
	for(int js = 0; js < nSteps; js++) {
	  double jPar = origPars[j] + (js-n)*jdVal; 
	  param_->GetPars()[j] = jPar;
	  double deltaLkh = 0.;
	  // Calculate likelihood for varied parameters
	  if( is != n || js != n ) {
	    for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
	      if( (*dataIt)->GetType() == TypeSmearDiJet )  { // Select DiJet events
		deltaLkh += (*dataIt)->chi2();
	      }
	    }    
	    deltaLkh -= offset;
	  }
	  hParScans2D.back()->Fill(iPar,jPar,deltaLkh);
	}
      }
      // Reset parameters to original values
      param_->GetPars()[i] = origPars[i];
      param_->GetPars()[j] = origPars[j];
    } // End of inner loop over parameters
  } // End of outer loop over parameters

  // Project out 1D parameter scans
  std::vector<TH1D*> hParScansTmp(nPar);
  int parScan2Didx = -1;
  int parScan1Didx = -1;

  // Create histograms of likelihood from i
  for(int i = 0; i < nPar; i++) {
    hParScansTmp[i] = 0;
    if( param_->isFixedPar(i) ) continue;

    // Inner loop over parameters
    for(int j = 0; j < i; j++) {
      if( param_->isFixedPar(j) ) continue;
      parScan2Didx++;
      const TH2D *h2 = hParScans2D[parScan2Didx];

      // Does 1D hist for parameter j exist?
      if( hParScansTmp[j] == 0 ) {
	TString name = "hParScan";
	name += ++parScan1Didx;
	TString title = ";";
	title += h2->GetYaxis()->GetTitle();
	title += ";- #Deltaln(L)";
	TH1D *h = new TH1D(name,title,
			   h2->GetNbinsY(),
			   h2->GetYaxis()->GetXmin(),
			   h2->GetYaxis()->GetXmax());
	// Copy y bin content in the central x bin
	for(int yBin = 1; yBin <= h2->GetNbinsY(); yBin++) {
	  h->SetBinContent(yBin,h2->GetBinContent(h2->GetBin(n+1,yBin)));
	}
	hParScansTmp[j] = h;
      }
      // Does 1D hist for parameter i exist?
      if( hParScansTmp[i] == 0 ) {
	TString name = "hParScan";
	name += ++parScan1Didx;
	TString title = ";";
	title += h2->GetXaxis()->GetTitle();
	title += ";- #Deltaln(L)";
	TH1D *h = new TH1D(name,title,
			   h2->GetNbinsX(),
			   h2->GetXaxis()->GetXmin(),
			   h2->GetXaxis()->GetXmax());
	// Copy x bin content in the central y bin
	for(int xBin = 1; xBin <= h2->GetNbinsX(); xBin++) {
	  h->SetBinContent(xBin,h2->GetBinContent(h2->GetBin(xBin,n+1)));
	}
	hParScansTmp[i] = h;
      }
    }
  }
  std::vector<TH1D*> hParScans;
  for(int i = 0; i < nPar; i++) {
    if( param_->isFixedPar(i) ) continue;
    hParScans.push_back(hParScansTmp[i]);
  }
  hParScansTmp.clear();


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsParameterScans.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","ParameterScans",0,0,600,600);

  for(size_t i = 0; i < hParScans.size(); i++) {
    hParScans[i]->SetMarkerStyle(20);
    drawPSPage(ps,c1,hParScans[i],"P",true);
  }
  for(size_t i = 0; i < hParScans2D.size(); i++) {
    drawPSPage(ps,c1,hParScans2D[i],"COLZ",true);
  }

  ps->Close();

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(size_t i = 0; i < hParScans.size(); i++) {
    rootfile.WriteTObject(hParScans[i]);
  }
  for(size_t i = 0; i < hParScans2D.size(); i++) {
    rootfile.WriteTObject(hParScans2D[i]);
  }


  // ----- Clean up -----
  for(size_t i = 0; i < hParScans2D.size(); i++) {
    delete hParScans2D[i];
  }
  for(size_t i = 0; i < hParScans.size(); i++) {
    delete hParScans[i];
  }
  for(size_t i = 0; i < lines.size(); i++) {
    delete lines[i];
  }
}  



//! These are the distributions of the negative logarithm
//! of the probability density of each event multiplied
//! by the event weight. That are the summands each
//! event adds to the negative log-likelihood of the fit.
//! (For a \f$\chi^{2}\f$-fit these would be the pull
//! distributions.)
//!
//! The plots are written to the files "jsResponse.root"
//! and "jsLogP.ps".
// --------------------------------------------------
void ControlPlotsJetSmearing::plotLogP() const {
  std::cout << "Creating -log(P) control plots\n";

  // ----- Fill vectors of -log(P) -----
  std::vector<double> logPstart;
  std::vector<double> logPend;
  std::vector<double> logPWstart;
  std::vector<double> logPWend;

  // Loop over data and fill -log(P) with fitted
  // parameter values
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    // Select DiJet events
    if( (*dataIt)->GetType() == TypeSmearDiJet )  {
      logPWend.push_back((*dataIt)->chi2());
      logPend.push_back((*dataIt)->chi2()/(*dataIt)->GetWeight());
    }
  }

  // Store fitted parameters
  std::vector<double> fittedPar(param_->GetNumberOfParameters());
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    fittedPar[i] = param_->GetPars()[i];
  }
  // Copy start values into parameter array
  std::vector<double> startParJet = bag_of<double>(config_->read<string>("jet start values",""));
  std::vector<double> startParGlobal = bag_of<double>(config_->read<string>("global jet start values",""));
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    if( i < param_->GetNumberOfJetParameters() )
      param_->GetPars()[i] = startParJet.at(i);
    else
      param_->GetPars()[i] = startParGlobal[i-param_->GetNumberOfJetParameters()];
  }
  // Loop over data and fill -log(P) with start
  // parameter values
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    // Select DiJet events
    if( (*dataIt)->GetType() == TypeSmearDiJet )  {
      logPWstart.push_back((*dataIt)->chi2());
      logPstart.push_back((*dataIt)->chi2()/(*dataIt)->GetWeight());
    }
  }
  // Copy back fitted values into parameter array
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    param_->GetPars()[i] = fittedPar.at(i);
  }


  // ----- Fill histograms of -log(P) -----
  std::sort(logPstart.begin(),logPstart.end());
  std::sort(logPend.begin(),logPend.end());
  double max = logPstart.back() > logPend.back() ? logPstart.back() : logPend.back();
  TH1F *hLogPstart = new TH1F("hLogPstart",";- ln(P)",100,0.,1.1*max);
  hLogPstart->SetLineWidth(2);
  hLogPstart->SetLineColor(4);
  TH1F *hLogPend = static_cast<TH1F*>(hLogPstart->Clone("hLogPend"));
  hLogPend->SetLineColor(2);

  std::sort(logPWstart.begin(),logPWstart.end());
  std::sort(logPWend.begin(),logPWend.end());
  max = logPWstart.back() > logPWend.back() ? logPWstart.back() : logPWend.back();
  TH1F *hLogPWstart = new TH1F("hLogPWstart",";- w #upoint ln(P)",100,0.,1.1*max);
  hLogPWstart->SetLineWidth(2);
  hLogPWstart->SetLineColor(4);
  TH1F *hLogPWend = static_cast<TH1F*>(hLogPWstart->Clone("hLogPWend"));
  hLogPWend->SetLineColor(2);

  for(size_t i = 0; i < logPend.size(); i++) {
    hLogPstart->Fill(logPstart[i]);
    hLogPend->Fill(logPend[i]);
    hLogPWstart->Fill(logPWstart[i]);
    hLogPWend->Fill(logPWend[i]);
  }
  logPstart.clear();
  logPend.clear();
  logPWstart.clear();
  logPWend.clear();

  max = hLogPend->GetMaximum() > hLogPstart->GetMaximum() ?
    hLogPend->GetMaximum() : hLogPstart->GetMaximum();
  hLogPstart->GetYaxis()->SetRangeUser(5E-2,5.*max);
  hLogPend->GetYaxis()->SetRangeUser(5E-2,5.*max);

  max = hLogPWend->GetMaximum() > hLogPWstart->GetMaximum() ?
    hLogPWend->GetMaximum() : hLogPWstart->GetMaximum();
  hLogPWstart->SetMaximum(5.*max);
  hLogPWend->SetMaximum(5.*max);


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsLogP.ps").c_str(),111);
  TCanvas * c1 = new TCanvas("c1","LogP",0,0,600,600);

  TLegend * leg = new TLegend(0.23,0.64,0.78,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hLogPstart,"Before fit","L");
  leg->AddEntry(hLogPend,"After fit","L");
  
  ps->NewPage();
  c1->cd()->SetLogy();
  hLogPstart->Draw();
  hLogPend->Draw("same");
  gPad->RedrawAxis();
  leg->Draw("same");
  c1->Draw();

  ps->NewPage();
  c1->cd()->SetLogy();
  hLogPWstart->Draw();
  hLogPWend->Draw("same");
  gPad->RedrawAxis();
  leg->Draw("same");
  c1->Draw();
  ps->Close();

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hLogPstart);
  rootfile.WriteTObject(hLogPend);
  rootfile.WriteTObject(hLogPWstart);
  rootfile.WriteTObject(hLogPWend);

  // ----- Clean up -----
  delete hLogPstart;
  delete hLogPend;
  delete hLogPWstart;
  delete hLogPWend;
  delete leg;
  delete c1;
  delete ps;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plotMeanResponseAndResolution() const {
  std::cout << "Creating response and resolution fits\n";

  // ----- Create histograms -----

  // Response vs log ptgen
  std::vector<double> ptBinEdges = bag_of<double>(config_->read<std::string>("plots pt bin edges","100 500"));
  TH2D* hRespVsLogPtGen = new TH2D("hRespVsLogPtGen",
				   ";p^{gen}_{T} (GeV);p^{jet}_{T} / p^{gen}_{T}",
				   (ptBinEdges.size()-1),&(ptBinEdges.front()),51,0.,2.);
  hRespVsLogPtGen->SetNdivisions(505);
  hRespVsLogPtGen->Sumw2();

  TH1D* hLogPtGen = new TH1D("hLogPtGen",";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
			     (ptBinEdges.size()-1),&(ptBinEdges.front()));
  hLogPtGen->SetNdivisions(505);


  // Response vs ptgen
  std::cout << ">>> " << ptBinEdges.front() << " - " << ptBinEdges.back() << std::endl;
  TH2D* hRespVsPtGen = new TH2D("hRespVsPtGen",
				";p^{gen}_{T} (GeV);p^{jet}_{T} / p^{gen}_{T}",
				30,ptBinEdges.front(),ptBinEdges.back(),51,0.,2.);
  hRespVsPtGen->SetNdivisions(505);
  hRespVsPtGen->Sumw2();

  TH1D* hPtGen = new TH1D("hPtGen",";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
			  500,ptBinEdges.front(),ptBinEdges.back());
  hPtGen->SetNdivisions(505);


  // ----- Fill histograms -----
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait); 

      hLogPtGen->Fill(dijet->jet1()->genPt(),dijet->GetWeight());
      hPtGen->Fill(dijet->jet1()->genPt(),dijet->GetWeight());

      hLogPtGen->Fill(dijet->jet2()->genPt(),dijet->GetWeight());
      hPtGen->Fill(dijet->jet2()->genPt(),dijet->GetWeight());

      hRespVsLogPtGen->Fill(dijet->jet1()->genPt(),dijet->jet1()->pt()/dijet->jet1()->genPt(),dijet->GetWeight());
      hRespVsLogPtGen->Fill(dijet->jet2()->genPt(),dijet->jet2()->pt()/dijet->jet2()->genPt(),dijet->GetWeight());
      hRespVsPtGen->Fill(dijet->jet1()->genPt(),dijet->jet1()->pt()/dijet->jet1()->genPt(),dijet->GetWeight());
      hRespVsPtGen->Fill(dijet->jet2()->genPt(),dijet->jet2()->pt()/dijet->jet2()->genPt(),dijet->GetWeight());
    }
  }


  // ----- Fit profiles -----
  std::vector<TH1D*> hResp(2);
  std::vector<TH1D*> hReso(2);
  for(int i = 0; i < 2; i++) {
    TH2D *h2 = hRespVsPtGen;
    if( i == 1 ) h2 = hRespVsLogPtGen;

    if( h2->GetXaxis()->GetXbins()->GetSize() == h2->GetNbinsX() +1) {
      TString name = "hResp";
      name += i;
      hResp[i] = new TH1D(name,";p^{gen}_{T} (GeV);< p^{jet}_{T} / p^{gen}_{T} >",
			  h2->GetNbinsX(),h2->GetXaxis()->GetXbins()->GetArray());

      name = "hReso";
      name += i;
      hReso[i] = new TH1D(name,
			  ";p^{gen}_{T} (GeV);#sigma(p^{jet}_{T} / p^{gen}_{T}) / < p^{jet}_{T} / p^{gen}_{T} >",
			  h2->GetNbinsX(),h2->GetXaxis()->GetXbins()->GetArray());

    } else {
      TString name = "hResp";
      name += i;
      hResp[i] = new TH1D(name,";p^{gen}_{T} (GeV);< p^{jet}_{T} / p^{gen}_{T} >",
			  h2->GetNbinsX(),h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax());

      name = "hReso";
      name += i;
      hReso[i] = new TH1D(name,
			  ";p^{gen}_{T} (GeV);#sigma(p^{jet}_{T} / p^{gen}_{T}) / < p^{jet}_{T} / p^{gen}_{T} >",
			  h2->GetNbinsX(),h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax());
    }
    hResp[i]->SetMarkerStyle(20);
    hResp[i]->SetLineWidth(2);
    hResp[i]->Sumw2();

    hReso[i]->SetMarkerStyle(20);
    hReso[i]->SetLineWidth(2);
    hReso[i]->Sumw2();


    // Get 1D slices and get mean, sigma
    TH1D* hSlice = new TH1D("hSlice","",h2->GetNbinsY(),h2->GetYaxis()->GetXmin(),
			    h2->GetYaxis()->GetXmax());
    hSlice->Sumw2();
    for(int xBin = 1; xBin <= h2->GetNbinsX(); xBin++) {
      hSlice->Reset();
      for(int yBin = 1; yBin <= h2->GetNbinsY(); yBin++) {
	hSlice->SetBinContent(yBin,h2->GetBinContent(h2->GetBin(xBin,yBin)));
	hSlice->SetBinError(yBin,h2->GetBinError(h2->GetBin(xBin,yBin)));
      }  

      double mean       = hSlice->GetMean();
      double width      = hSlice->GetRMS();
      if( width < 0.1 ) width = 0.1;
      hSlice->Fit("gaus","QN","",mean-1.5*width,mean+1.5*width);
      TF1 *f = static_cast<TF1*>(gROOT->GetFunction("gaus"));
      mean = f->GetParameter(1);
      double meanError = f->GetParError(1);
      width = f->GetParameter(2);
      double widthError = f->GetParError(2);
      
      hResp[i]->SetBinContent(xBin,mean);
      hResp[i]->SetBinError(xBin,meanError);
      if( mean ) {
	hReso[i]->SetBinContent(xBin,width/mean);
	hReso[i]->SetBinError(xBin,widthError/mean);
      }
    }
    delete hSlice;
    
    hResp[i]->GetYaxis()->SetRangeUser(0.9,1.2);
    hReso[i]->GetYaxis()->SetRangeUser(0.,0.4);
  }


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsMeanRespAndReso.ps").c_str(),111);
  TCanvas * c1 = new TCanvas("c1","MeanRespAndReso",0,0,600,600);
  drawPSPage(ps,c1,hPtGen);
  drawPSPage(ps,c1,hLogPtGen,"",true);
  for(size_t i = 0; i < hResp.size(); i++) {
    if( i == 1 ) gPad->SetLogx();
    drawPSPage(ps,c1,hResp[i],"P");
    if( i == 1 ) gPad->SetLogx();
    drawPSPage(ps,c1,hReso[i],"P");
  }

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(size_t i = 0; i < hResp.size(); i++) {
    rootfile.WriteTObject(hResp[i]);
    rootfile.WriteTObject(hReso[i]);
  }
  rootfile.WriteTObject(hLogPtGen);
  rootfile.WriteTObject(hPtGen);


  // ----- Clean up -----
  for(size_t i = 0; i < hResp.size(); i++) {
    delete hResp[i];
    delete hReso[i];
  }
  delete hLogPtGen;
  delete hPtGen;
  delete c1;
  delete ps;
}



//!  \brief Draw TObject on one page of a ps file
//!
//!  Opens a new page in the PostScript file \p ps
//!  and draws the TObject \p obj on it.
//!
//!  \param ps     A new page is added to this file
//!                containing the TObject \p obj
//!  \param can    The \p obj is drawn on this canvas
//!  \param obj    The TObject to be drawn
//!  \param option Draw options
//!  \param log    Sets log-scale on last axis if true
// --------------------------------------------------
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, std::string option, bool log) const {
  std::vector<TObject*> objs;
  objs.push_back(obj);
  drawPSPage(ps,can,objs,option,log);
}



//!  \brief Draw TObjects on one page of a ps file
//!
//!  Opens a new page in the PostScript file \p ps
//!  and draws all TObjects in \p obs on it. If \p objs
//!  contains more than one TObject, the draw option
//!  "same" is automatically used.
//!
//!  \param ps     A new page is added to this file
//!                containing all TObjects in \p objs
//!  \param can    The \p objs are drawn on this canvas
//!  \param objs   Contains the objects that are drawn
//!  \param option Draw options (except for "same",
//!                see above)
//!  \param log    Sets log-scale on last axis if true
// --------------------------------------------------
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, std::string option, bool log) const {
  ps->NewPage();
  can->cd();
  for( size_t i = 0; i < objs.size(); i++ ) {
    if( i == 0 ) objs.at(i)->Draw(option.c_str());
    else         objs.at(i)->Draw((option.append("same")).c_str());
  }
  gPad->RedrawAxis();
  std::string hist = objs.at(0)->ClassName();
  if( hist.compare("TH1F") == 0 ) {
    if( log ) can->SetLogy(1);
    else      can->SetLogy(0);
  } else if( hist.compare("TH2F") == 0 ) {
    can->SetLogy(0);
    if( log ) can->SetLogz(1);
    else      can->SetLogz(0);
  }

  can->Draw();
}



//!  \brief Find y-axis range
//!
//!  Sets \p min and \p max to the minimum (non-zero) and
//!  maximum bin content of \p h, respectively.
// --------------------------------------------------
void ControlPlotsJetSmearing::findYRange(const TH1F * h, double& min, double& max) const {
  min = 10000.;
  max = 0.;
  for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
    double val = h->GetBinContent(bin);
    if( val > 0. && val < min ) min = val;
    if( val > 0. && val > max ) max = val;
  }
  if( min > max ) {
    min = 1E-3;
    max = 1;
  }
}



//!  \brief Set default \p gStyle options
// --------------------------------------------------
void ControlPlotsJetSmearing::setGStyle() const
{
  gStyle->SetPalette(1);

  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);

  // For the statistics box:
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("neMR");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.86);              
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.22);

  // For the legend
  gStyle->SetLegendBorderSize(1);

  //  Margins
  // -------------------------------------------
  gStyle->SetPadTopMargin(0.16);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(0.16);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.12);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.515);
  gStyle->SetTitleH(0.06);
  gStyle->SetTitleXOffset(0);
  gStyle->SetTitleYOffset(0);
  gStyle->SetTitleBorderSize(0);

  // For the axis labels:
  //  For the axis labels and titles
  // -------------------------------------------
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  // For the axis labels:
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  
  // For the axis titles:
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}
 
 

//!  \brief Adjust y-axis range
//!
//!  Sets the y-axis range of \p h from
//!  <tt> c1 * min</tt> to <tt> c2 * max</tt>,
//!  where \p min and \p max are the minimal and
//!  the maximal bin non-zero content, respectively.
//!  If <tt>min < minLimit</tt>, \p minLimit is used
//!  instead as minimum.
// --------------------------------------------------
void ControlPlotsJetSmearing::setYRange(TH1F * h, double c1, double c2, double minLimit) const {
  double min = 0.;
  double max = 0.;
  findYRange(h,min,max);
  min *= c1;
  max *= c2;
  if( min < minLimit ) min = minLimit;
  h->GetYaxis()->SetRangeUser( min, max );
}

void ControlPlotsJetSmearing::normHist(TH1F *h, std::string option) const { 
  if( h->Integral(option.c_str()) ) h->Scale(1./h->Integral(option.c_str())); 
}

void ControlPlotsJetSmearing::normHist(TH1F *h, double min, double max, std::string option) const { 
  double norm = h->Integral(h->FindBin(min),h->FindBin(max),option.c_str());
  if( norm ) h->Scale(1./norm);
}

 //! \brief Convert to std::string
 //!
 //! Converts a given object to an std::string
 //! using std::stringstream.
 //!
 //! \param t Object to be converted to a string
 //! \return String representation of t
 // --------------------------------------------------
 template <class T> std::string ControlPlotsJetSmearing::toString(const T& t) const {
   std::stringstream ss;
   ss << t;
   return ss.str();
 }
