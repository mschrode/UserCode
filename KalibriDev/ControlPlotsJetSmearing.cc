// $Id: ControlPlotsJetSmearing.cc,v 1.10 2010/01/26 17:49:22 mschrode Exp $

#include "ControlPlotsJetSmearing.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"


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
ControlPlotsJetSmearing::ControlPlotsJetSmearing(const std::string& configfile, const std::vector<Event*> * data, TParameters * param)
  : data_(data),
    config_(new ConfigFile(configfile.c_str())),
    param_(param),
    respNBins_(40),
    respMin_(0.),
    respMax_(3.),
    dir_("./controlPlots")
{
  // Override possible existing root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"RECREATE");
  rootfile.Close();
  setGStyle();
}



//!  \brief Draw response control plots for events
//!         of type \p SmearDiJet
// --------------------------------------------------
void ControlPlotsJetSmearing::plotResponse() const
{
  std::cout << "Creating response control plots\n";

  // --- Create histograms of response and spectrum ---------------------
  double tMin = config_->read<double>("DiJet integration min",0.);
  double tMax = config_->read<double>("DiJet integration max",1.);

  // The vectors contain the distribution for
  // 0 - Resp function evaluated at mean of ptGen distribution
  // Larger entries: Measured reponse in a certain ptGen bin; resp
  //                 function evaluated at mean of ptGen spectrum 
  //                 of that bin
  int nPlotBins = 4;
  std::vector<TH1F*> hRespMeasAbs(nPlotBins);   // The measured response ptJet / ptGen absolute entries
  std::vector<TH1F*> hRespMeas(nPlotBins);      // The measured response ptJet / ptGen
  std::vector<TH1F*> hRespFitStart(nPlotBins);  // The response pdf with start values
  std::vector<TH1F*> hRespFit(nPlotBins);       // The fitted response pdf
  std::vector<TH1F*> hRespFitErrStat(nPlotBins);// The fitted response pdf with fitted errors
  std::vector<TH1F*> hRespFitStep(nPlotBins);   // Step function part of the response pdf
  std::vector<TH1F*> hRespFitGaus(nPlotBins);   // Gauss part of the response pdf
  std::vector<TH1F*> hRespFitSum(nPlotBins);    // Sum of step and Gauss part
  std::vector<TH1F*> hRatio(nPlotBins);
  TH1F * hTruthPDF = 0;      // Truth pdf
  TH1F * hTruthPDFErrStat = 0;      // Truth pdf
  TH1F * hPtGen = 0;         // PtGen spectrum
  TH1F * hPtHat = 0;         // PtHat spectrum
  TH1F * hPtDijet = 0;       // Dijet spectrum

  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    std::string name = "hRespMeasAbs_" + toString(plotBin);
    hRespMeasAbs.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{gen}_{T};dN / dR",
			    respNBins_,respMin_,respMax_);
    hRespMeasAbs.at(plotBin)->SetLineWidth(2);

    name = "hRespMeas_" + toString(plotBin);
    hRespMeas.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{gen}_{T};1/(Nw)  dN / dR",
			       respNBins_,respMin_,respMax_);
    hRespMeas.at(plotBin)->Sumw2();
    hRespMeas.at(plotBin)->SetLineWidth(2);

    name = "hRespFit_" + toString(plotBin);
    hRespFit.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",
			      5*respNBins_,respMin_,respMax_);
    hRespFit.at(plotBin)->SetLineColor(2);
    hRespFit.at(plotBin)->SetLineWidth(2);

    name = "hRespFitErrStat_" + toString(plotBin);
    hRespFitErrStat.at(plotBin) = static_cast<TH1F*>(hRespFit[plotBin]->Clone(name.c_str()));
    hRespFitErrStat.at(plotBin)->SetFillColor(45);

    name = "hRespFitStart_" + toString(plotBin);
    hRespFitStart.at(plotBin) = static_cast<TH1F*>(hRespFit[plotBin]->Clone(name.c_str()));
    hRespFitStart.at(plotBin)->SetLineStyle(2);

    name = "hRespFitStep_" + toString(plotBin);
    hRespFitStep.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",
				  config_->read<int>("Response pdf nsteps",10),
				  config_->read<double>("Response pdf min",0.),
				  config_->read<double>("Response pdf max",1.8));
    hRespFitStep.at(plotBin)->Sumw2();
    hRespFitStep.at(plotBin)->SetLineColor(9);
    hRespFitStep.at(plotBin)->SetLineWidth(2);

    name = "hRespFitGaus_" + toString(plotBin);
    hRespFitGaus.at(plotBin) = static_cast<TH1F*>(hRespFit[plotBin]->Clone(name.c_str()));
    hRespFitGaus.at(plotBin)->SetLineColor(8);

    name = "hRespFitSum_" + toString(plotBin);
    hRespFitSum.at(plotBin) = static_cast<TH1F*>(hRespFit[plotBin]->Clone(name.c_str()));
    hRespFitSum.at(plotBin)->SetLineColor(1);

    name = "hRatio_" + toString(plotBin);
    hRatio.at(plotBin) = new TH1F(name.c_str(),";Prediction / truth",
				  respNBins_,respMin_,respMax_);
    hRatio.at(plotBin)->SetLineWidth(2);
  }
  
  hPtGen = new TH1F("hPtGen",";p^{gen}_{T} (GeV);1/(Nw)  dN / dp^{gen}_{T}  1 / (GeV)"
		    ,25,0.8*tMin,1.1*tMax);
  hPtGen->GetXaxis()->SetNdivisions(505);
  hPtGen->Sumw2();
  hPtGen->SetLineWidth(2);
  
  hPtHat = new TH1F("hPtHat",";#hat{p}_{T} (GeV);1/(Nw)  dN / d#hat{p}_{T}  1 / (GeV)",
		    25,0.8*tMin,1.1*tMax);
  hPtHat->GetXaxis()->SetNdivisions(505);
  hPtHat->Sumw2();
  hPtHat->SetLineWidth(2);

  hPtDijet = new TH1F("hPtDijet",";p^{dijet}_{T} (GeV);1/(Nw)  dN / dp^{dijet}_{T}  1 / (GeV)",
		      25,0.8*tMin,1.1*tMax);
  hPtDijet->GetXaxis()->SetNdivisions(505);
  hPtDijet->Sumw2();
  hPtDijet->SetLineWidth(2);

  hTruthPDF = new TH1F("hTruthPDF",";p^{true}_{T} (GeV);1 / (Nw)  dN / dp^{true}_{T}  1 /  (GeV)",
		       5*respNBins_,0.8*tMin,1.1*tMax);
  hTruthPDF->SetLineColor(2);
  hTruthPDF->SetLineWidth(2);

  hTruthPDFErrStat = static_cast<TH1F*>(hTruthPDF->Clone("hTruthPDFErrStat"));
  hTruthPDFErrStat->SetFillColor(45);


  // --- Define ptDijet bins -------------
  int nPtGenBins = nPlotBins - 1;
  double deltaPtGen = (tMax - tMin) / nPtGenBins;
  std::vector<double> ptGenBinEdges(nPtGenBins+1);
  for(size_t i = 0; i < ptGenBinEdges.size(); i++) {
    ptGenBinEdges.at(i) = tMin + i*deltaPtGen;
  }


  // --- Fill histograms of measured response --------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      hPtHat->Fill( dijet->ptHat(), dijet->GetWeight() );
      for(int i = 0; i < 2; i++) {        // Loop over both jets
	const Jet * jet = dijet->jet1();
	if( i == 1 ) jet = dijet->jet2();

	hPtGen->Fill( jet->genPt(), dijet->GetWeight() );

	hRespMeasAbs.at(0)->Fill( jet->pt() / jet->genPt(), dijet->GetWeight() );
	hRespMeas.at(0)->Fill( jet->pt() / jet->genPt(), dijet->GetWeight() );

	for(int i = 0; i < nPtGenBins; i++) {
	  if( ptGenBinEdges.at(i) <= jet->genPt() && jet->genPt() < ptGenBinEdges.at(i+1) ) {
	    hRespMeasAbs.at(i+1)->Fill( jet->pt() / jet->genPt(), dijet->GetWeight() );
	    hRespMeas.at(i+1)->Fill( jet->pt() / jet->genPt(), dijet->GetWeight() );
	    continue;
	  }
	}
      }
      hPtDijet->Fill( dijet->dijetPt(), dijet->GetWeight() );
    }
  } // End of loop over data
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    normHist(hRespMeas.at(plotBin),"width");
  }
  normHist(hPtGen,"width");
  normHist(hPtHat,"width");
  normHist(hPtDijet,"width");


  // --- Fill histograms of fitted response ----------------
  std::string param = config_->read<std::string>("Parametrization Class","");

  // Get parameters
  std::vector<double> fittedPar(param_->GetNumberOfParameters());
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    fittedPar.at(i) = param_->GetPars()[i];
  }

  std::vector<double> scale = bag_of<double>(config_->read<string>("jet parameter scales",""));
  std::vector<double> startParJet = bag_of<double>(config_->read<string>("jet start values",""));
  std::vector<double> startParGlobal = bag_of<double>(config_->read<string>("global jet start values",""));
  std::vector<double> auxPar = bag_of<double>(config_->read<string>("mean response parameters","1 0"));

  // Find bin centers for response function evaluation
  std::vector<double> ptGenBinCenters(nPlotBins);
  ptGenBinCenters.at(0) = hPtGen->GetMean();
  for(int i = 1; i < nPlotBins; i++) {
    int j = i - 1;
    ptGenBinCenters.at(i) = 0.5 * ( ptGenBinEdges.at(j) + ptGenBinEdges.at(j+1) );
  }


  SmearData * smearData = dynamic_cast<SmearData*>(data_->front());
  if( smearData ) {
    // Loop over plotBins
    for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
      // Interpolated response function
      for(int bin = 1; bin <= hRespFit[plotBin]->GetNbinsX(); bin++) {
	double r = hRespFit[plotBin]->GetBinCenter(bin);
	hRespFit[plotBin]->SetBinContent(bin,smearData->respPDF(r,ptGenBinCenters[plotBin]));
	hRespFitErrStat[plotBin]->SetBinContent(bin,smearData->respPDF(r,ptGenBinCenters[plotBin]));
	hRespFitErrStat[plotBin]->SetBinError(bin,smearData->respPDFSigma(r,ptGenBinCenters[plotBin]));

//  	if( plotBin == 0 ) {
//  	  std::cout << r << ":  " << smearData->respPDF(r,ptGenBinCenters.at(plotBin)) << " +/- " << smearData->respPDFSigma(r,ptGenBinCenters.at(plotBin)) << std::endl;
//  	}
      }
      // Ratio plot
      for(int bin = 1; bin <= hRatio[plotBin]->GetNbinsX(); bin++) {
	double rMin = hRatio[plotBin]->GetBinLowEdge(bin);
	double rMax = rMin + hRatio[plotBin]->GetBinWidth(bin);
	int min = hRespFit[plotBin]->FindBin(rMin);
	int max = hRespFit[plotBin]->FindBin(rMax);
	double pred = (hRespFit[plotBin]->Integral(min,max))/(1+max-min);
	double truth = hRespMeas[plotBin]->GetBinContent(bin);
	double ratio = 0.;
	double error = 0.;
	if( truth ) {
	  ratio = pred / truth;
	  error = pred / truth / truth * hRespMeas[plotBin]->GetBinError(bin);
	}
	hRatio[plotBin]->SetBinContent(bin,ratio);
	hRatio[plotBin]->SetBinError(bin,error);
      }

      // Interpolated fit function with start values
      // Copy start values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	if( i < param_->GetNumberOfJetParameters() )
	  param_->GetPars()[i] = startParJet.at(i);
	else
	  param_->GetPars()[i] = startParGlobal.at(i-param_->GetNumberOfJetParameters());
      }
      // Plot response function
      for(int bin = 1; bin <= hRespFitStart.at(plotBin)->GetNbinsX(); bin++) {
	double r = hRespFitStart.at(plotBin)->GetBinCenter(bin);
	hRespFitStart.at(plotBin)->SetBinContent(bin,smearData->respPDF(r,ptGenBinCenters.at(plotBin)));
      }
      // Copy back fitted values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = fittedPar.at(i);
      }

      // In case of interpolated step + gauss parametrization
      if( param == "SmearParametrizationStepGaussInter" ) {
	// Step part of fit function
	for(int bin = 1; bin <= hRespFitStep.at(plotBin)->GetNbinsX(); bin++) {
	  double val  = scale.at(bin+1)*(smearData->respPar()[bin+1]);
	  hRespFitStep.at(plotBin)->SetBinContent(bin,val);
	}
	normHist(hRespFitStep.at(plotBin),"width");
	hRespFitStep.at(plotBin)->Scale(1. - scale.at(0)*(smearData->respPar()[0]));
	
	// Gauss part of fit function
	for(int bin = 1; bin <= hRespFitGaus.at(plotBin)->GetNbinsX(); bin++) {
	  // Mean
	  double mu = auxPar.at(0);
	  // Width
// 	  double a1 = scale.at(1)*(smearData->respPar()[1]);
// 	  double a2 = scale.at(2)*(smearData->respPar()[2]);
// 	  double a3 = scale.at(3)*(smearData->respPar()[3]);
// 	  double sigma = sqrt( a1*a1/ptGenBinCenters.at(plotBin)/ptGenBinCenters.at(plotBin)
// 			       + a2*a2/ptGenBinCenters.at(plotBin) + a3*a3 );
	  double sigma = scale[1]*(smearData->respPar()[1]);
	  // pdf
	  double c     = scale.at(0)*(smearData->respPar()[0]);
	  double r     = hRespFitGaus.at(plotBin)->GetBinCenter(bin);
	  double val   = c * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	  hRespFitGaus.at(plotBin)->SetBinContent(bin,val);
	}
      
	// Sum
	for(int binGaus = 1; binGaus <= hRespFitGaus.at(plotBin)->GetNbinsX(); binGaus++) {
	  int binStep = hRespFitStep.at(plotBin)->FindBin(hRespFitGaus.at(plotBin)->GetBinCenter(binGaus));
	  double val = hRespFitStep.at(plotBin)->GetBinContent(binStep)
	    + hRespFitGaus.at(plotBin)->GetBinContent(binGaus);
	  hRespFitSum.at(plotBin)->SetBinContent(binGaus,val);
	}
      } else if( param == "SmearParametrizationCrystalBall" ) {
      } else {
	std::cout << "WARNING: No controlplots implemented for parametrization '" << param << "'\n";
      }
    } // End of loop over plotBins
  } // End if( smearData )


  // --- Fill histograms of fitted truth spectrum -----------

  // Fill histogram of assumed truth pdf
  // and fit with 1/x^n function
  DataIt datait = data_->begin();
  while( (*datait)->GetType() != TypeSmearDiJet  &&  datait != data_->end() ) datait++;
  if( datait != data_->end() ) {
    SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

    for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
      double t = hTruthPDF->GetBinCenter(bin);
      hTruthPDF->SetBinContent(bin,dijet->truthPDF(t));
    }
    normHist(hTruthPDF,"width");
    for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
      hTruthPDFErrStat->SetBinContent(bin,hTruthPDF->GetBinContent(bin));
      double t = hTruthPDFErrStat->GetBinCenter(bin);
      hTruthPDFErrStat->SetBinError(bin,dijet->truthPDFSigma(t));
    }
  }


  // --- Find populated x-axis ranges -----------------------
  int maxBin = 0;
  int minBin = 1;
  for(int bin = 1; bin <= hRespMeas.at(0)->GetNbinsX(); bin++) {
    if( hRespMeas.at(0)->GetBinContent(bin) > 0 ) maxBin = bin;
    if( minBin > maxBin ) minBin = bin;
  }
  if( maxBin < hRespMeas.at(0)->GetNbinsX() ) maxBin++;
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    hRespMeas.at(plotBin)->GetXaxis()->SetRange(minBin,maxBin);
    hRespMeasAbs.at(plotBin)->GetXaxis()->SetRange(minBin,maxBin);
    hRatio.at(plotBin)->GetXaxis()->SetRange(minBin,maxBin);
  }

  hPtGen->GetXaxis()->SetRangeUser(tMin,tMax);
  hPtHat->GetXaxis()->SetRangeUser(tMin,tMax);
  

  // --- Set y-axis ranges ----------------------------------
  double yMin = 10000.;
  double yMax = 0.;
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    double min = 0.;
    double max = 0.;
    findYRange(hRespMeas.at(plotBin),min,max);
    min *= 0.5;
    max *= 80.;
    if( min < yMin ) yMin = min;
    if( max > yMax ) yMax = max;
    if( yMin < 8E-5 ) yMin = 8E-5;
  }
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    hRespMeas.at(plotBin)->GetYaxis()->SetRangeUser(yMin,yMax);
    setYRange(hRespMeasAbs.at(plotBin),0.5,50.);
  }
  setYRange(hPtDijet, 0.5, 100.);
  setYRange(hPtGen, 0.5, 100.);
  setYRange(hPtHat, 0.5, 100.);


  // --- Plot histograms -----------------------------------
  // Label bins
  double ptDijetMin = config_->read<double>("Et min cut on dijet",-1.);
  double ptDijetMax = config_->read<double>("Et max cut on dijet",-1.);

  std::vector<TLegend*> legPtRange(nPlotBins);
  std::vector<TLegend*> legPtRangeAndCenters(nPlotBins);
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    legPtRange.at(plotBin) = new TLegend(0.23,0.72,0.78,0.8);
    legPtRange.at(plotBin)->SetBorderSize(0);
    legPtRange.at(plotBin)->SetFillColor(0);
    legPtRange.at(plotBin)->SetTextFont(42);

    legPtRangeAndCenters.at(plotBin) = new TLegend(0.23,0.65,0.8,0.8);
    legPtRangeAndCenters.at(plotBin)->SetBorderSize(0);
    legPtRangeAndCenters.at(plotBin)->SetFillColor(0);
    legPtRangeAndCenters.at(plotBin)->SetTextFont(42);

    std::string label;
    if( plotBin == 0 ) {
      //      label = toString(ptDijetMin) + " < p^{dijet}_{T} < " + toString(ptDijetMax) + " GeV";
      label = toString(ptGenBinEdges.front()) + " < p^{gen}_{T} < " + toString(ptGenBinEdges.back()) + " GeV";
      legPtRange.at(plotBin)->AddEntry(hRespMeas.at(plotBin),label.c_str(),"L");
      legPtRangeAndCenters.at(plotBin)->AddEntry(hRespMeas.at(plotBin),label.c_str(),"L");
    } else {
      int i = plotBin - 1;
      label = toString(ptGenBinEdges.at(i)) + " < p^{gen}_{T} < " + toString(ptGenBinEdges.at(1+i)) + " GeV";
      legPtRange.at(plotBin)->AddEntry(hRespMeas.at(plotBin),label.c_str(),"L");
      legPtRangeAndCenters.at(plotBin)->AddEntry(hRespMeas.at(plotBin),label.c_str(),"L");
    }
    label = "p_{T} = " + toString(ptGenBinCenters.at(plotBin)) + " GeV";
    legPtRangeAndCenters.at(plotBin)->AddEntry(hRespFit.at(plotBin),label.c_str(),"L");
  }

  // Write histos to ps file
  TPostScript * const ps = new TPostScript((dir_+"/jsResponse.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Jet Response",0,0,600,600);

  TLegend *legFitStart = new TLegend(0.23,0.5,0.5,0.65);
  legFitStart->SetBorderSize(0);
  legFitStart->SetFillColor(0);
  legFitStart->SetTextFont(42);
  legFitStart->AddEntry(hRespFitStart.at(0),"At start","L");
  legFitStart->AddEntry(hRespFit.at(0),"After fit","L");

  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    // Measured and fitted response
    ps->NewPage();
    c1->cd();
    hRespMeasAbs.at(plotBin)->Draw();
    legPtRange.at(plotBin)->Draw("same");
    c1->SetLogy();
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas.at(plotBin)->Draw();
    legPtRange.at(plotBin)->Draw("same");
    c1->SetLogy();
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas.at(plotBin)->Draw();
    hRespFitErrStat.at(plotBin)->Draw("E3 same");
    hRespMeas.at(plotBin)->Draw("same");
    hRespFit.at(plotBin)->Draw("Lsame");
    legPtRangeAndCenters.at(plotBin)->Draw("same");
    c1->SetLogy();
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas.at(plotBin)->Draw();
    hRespFit.at(plotBin)->Draw("Lsame");
    hRespFitStart.at(plotBin)->Draw("Lsame");
    c1->SetLogy();
    legPtRangeAndCenters.at(plotBin)->Draw("same");
    legFitStart->Draw("same");
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas.at(plotBin)->Draw();
    if( param == "SmearParametrizationStepGaussInter" ) {
      hRespFitStep.at(plotBin)->Draw("same");
      hRespFitGaus.at(plotBin)->Draw("same");
      //hRespFitSum.at(plotBin)->Draw("same");
    }
    hRespFit.at(plotBin)->Draw("Lsame");
    legPtRangeAndCenters.at(plotBin)->Draw("same");
    c1->SetLogy();
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRatio[plotBin]->Draw();
    c1->SetLogy(0);
    c1->Draw();
  }

  // Truth spectrum
  double n = param_->GetGlobalJetParRef()[0];
  TLegend *legPtGen = new TLegend(0.4,0.67,0.8,0.8);
  legPtGen->SetBorderSize(0);
  legPtGen->SetFillColor(0);
  legPtGen->SetTextFont(42);
  char entry[50];
  sprintf(entry,"#propto 1 / (p^{gen}_{T})^{%.1f}",n);
  legPtGen->AddEntry(hTruthPDF,entry,"L");

  ps->NewPage();
  c1->cd();
  hPtGen->Draw();
  hTruthPDFErrStat->Draw("E3same");
  hPtGen->Draw("same");
  hTruthPDF->Draw("Lsame");
  legPtRange.at(0)->Draw("same");
  c1->SetLogy();
  c1->Draw();

  std::vector<TObject*> objs;
  objs.clear();
  objs.push_back(hPtHat);
  objs.push_back(hTruthPDF);
  objs.push_back(legPtRange.at(0));
  drawPSPage(ps,c1,objs,"",true);

  objs.clear();
  objs.push_back(hPtDijet);
  objs.push_back(hTruthPDF);
  objs.push_back(legPtRange.at(0));
  drawPSPage(ps,c1,objs,"",true);


  // Write histos to root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    rootfile.WriteTObject(hRespMeasAbs.at(plotBin));
    rootfile.WriteTObject(hRespMeas.at(plotBin));
    rootfile.WriteTObject(hRespFit.at(plotBin));
    rootfile.WriteTObject(hRespFitErrStat.at(plotBin));
    rootfile.WriteTObject(hRespFitStart.at(plotBin));
    rootfile.WriteTObject(hRespFitStep.at(plotBin));
    rootfile.WriteTObject(hRespFitGaus.at(plotBin));
    rootfile.WriteTObject(hRespFitSum.at(plotBin));
    rootfile.WriteTObject(hRatio.at(plotBin));
  }
  rootfile.WriteTObject(hPtGen);
  rootfile.WriteTObject(hPtHat);
  rootfile.WriteTObject(hPtDijet);
  rootfile.WriteTObject(hTruthPDF);
  rootfile.WriteTObject(hTruthPDFErrStat);

  rootfile.Close();


  // --- Clean up ------------------------------------------
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    delete hRespMeasAbs.at(plotBin);
    delete hRespMeas.at(plotBin);
    delete hRespFit.at(plotBin);
    delete hRespFitErrStat.at(plotBin);
    delete hRespFitStart.at(plotBin);
    delete hRespFitStep.at(plotBin);
    delete hRespFitGaus.at(plotBin);
    delete hRespFitSum.at(plotBin);
    delete legPtRangeAndCenters.at(plotBin);
    delete hRatio[plotBin];
  }
  delete legFitStart;
  delete hPtGen;
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
  TH1F * hPtHat = new TH1F("hPtHat",";#hat{p}_{T} (GeV);1 / N  dN / d#hat{p}_{T}  1 / (GeV)",
			   50,0.9*minPtHat, 1.1*maxPtHat);
  hPtHat->SetLineWidth(2);
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

    h = new TH1F(genJetNames[i].c_str(),";p^{gen}_{T} (GeV);1 / N  dN / dp^{gen}_{T}  1 / (GeV)",
		 50,0.9*minGenJetPt,1.1*maxGenJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    hGenJetPt.push_back(h);

    h = new TH1F(calJetNames[i].c_str(),";p^{jet}_{T} (GeV);1 / N  dN / dp^{jet}_{T}  1 / (GeV)",
		 50,0.9*minCalJetPt,1.1*maxCalJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
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
  for(size_t i = 0; i < hGenJetPt.size(); i++) {
    normHist( hGenJetPt.at(i) );
    normHist( hCalJetPt.at(i) );
  }
  normHist( hPtHat );
  normHist( hDijetPt );
  normHist( h3rdJetPt );
  normHist( hRel3rdJetPt );
  normHist( hDeltaPhi );



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

void ControlPlotsJetSmearing::normHist(TH1F * h, std::string option) const
{ 
  if( h->Integral(option.c_str()) ) h->Scale(1./h->Integral(option.c_str())); 
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
