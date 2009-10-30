// $Id: ControlPlotsJetSmearing.cc,v 1.3 2009/10/30 15:47:14 mschrode Exp $

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

#include "SmearData.h"
#include "SmearDiJet.h"
#include "SmearPhotonJet.h"


//!  \brief Constructor
//!
//!  By default, all controlplots are written to the
//!  directory \p ./controlplots/
//!
//!  \param configfile Name of the configuration file
//!  \param data The data
//!  \param param The parametrization
// --------------------------------------------------
ControlPlotsJetSmearing::ControlPlotsJetSmearing(const std::string& configfile, const std::vector<TData*> * data, TParameters * param)
  : data_(data),
    config_(new ConfigFile(configfile.c_str())),
    param_(param),
    respNBins_(80),
    respMin_(0.),
    respMax_(3.),
    dir_("./controlplots")
{
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
  // 0 - Resp function evaluated at mean of intervall (tMin,tMax)
  // 1 - Resp function evaluated at mean of ptGen distribution
  // Larger entries: Measured reponse in a certain ptGen bin; resp
  //                 function evaluated at mean of ptGen spectrum 
  //                 of that bin
  int nPlotBins = 6;
  std::vector<TH1F*> hRespMeasAbs(nPlotBins);   // The measured response ptJet / ptGen absolute entries
  std::vector<TH1F*> hRespMeas(nPlotBins);      // The measured response ptJet / ptGen
  std::vector<TH1F*> hRespFitStart(nPlotBins);  // The response pdf with start values
  std::vector<TH1F*> hRespFit(nPlotBins);       // The fitted response pdf
  std::vector<TH1F*> hRespFitStep(nPlotBins);   // Step function part of the response pdf
  std::vector<TH1F*> hRespFitGaus(nPlotBins);   // Gauss part of the response pdf
  std::vector<TH1F*> hRespFitSum(nPlotBins);    // Sum of step and Gauss part
  TH1F * hTruthPDF = 0;      // Truth pdf
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
    hRespFit.at(plotBin)->Sumw2();

    name = "hRespFitStart_" + toString(plotBin);
    hRespFitStart.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",
				   5*respNBins_,respMin_,respMax_);
    hRespFitStart.at(plotBin)->SetLineColor(2);
    hRespFitStart.at(plotBin)->SetLineWidth(2);
    hRespFitStart.at(plotBin)->SetLineStyle(2);
    hRespFitStart.at(plotBin)->Sumw2();

    name = "hRespFitStep_" + toString(plotBin);
    hRespFitStep.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",
				  config_->read<int>("Response pdf nsteps",10),
				  config_->read<double>("Response pdf min",0.),
				  config_->read<double>("Response pdf max",1.8));
    hRespFitStep.at(plotBin)->Sumw2();
    hRespFitStep.at(plotBin)->SetLineColor(9);
    hRespFitStep.at(plotBin)->SetLineWidth(2);

    name = "hRespFitGaus_" + toString(plotBin);
    hRespFitGaus.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",
					5*respNBins_,respMin_,respMax_);
    hRespFitGaus.at(plotBin)->Sumw2();
    hRespFitGaus.at(plotBin)->SetLineColor(8);
    hRespFitGaus.at(plotBin)->SetLineWidth(2);

    name = "hRespFitSum_" + toString(plotBin);
    hRespFitSum.at(plotBin) = new TH1F(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",
				 5*respNBins_,respMin_,respMax_);
    hRespFitSum.at(plotBin)->Sumw2();
    hRespFitSum.at(plotBin)->SetLineColor(1);
    hRespFitSum.at(plotBin)->SetLineWidth(2);
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


  // --- Define ptDijet bins -------------
  //! \todo Should be bins of ptGen
  double ptDijetMin = config_->read<double>("Et min cut on dijet",-1.);
  double ptDijetMax = config_->read<double>("Et max cut on dijet",-1.);
  int nPtDijetBins = nPlotBins - 2;
  double deltaPtDijet = (ptDijetMax - ptDijetMin) / nPtDijetBins;
  std::vector<double> ptDijetBinEdges(nPtDijetBins+1);
  for(size_t i = 0; i < ptDijetBinEdges.size(); i++) {
    ptDijetBinEdges.at(i) = ptDijetMin + i*deltaPtDijet;
  }


  // --- Fill histograms of measured response --------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      for(int i = 0; i < 2; i++) {        // Loop over both jets
	TJet * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 ) jet = static_cast<TJet*>(dijet->GetSecondMess());

	hPtGen->Fill( jet->genPt, dijet->GetWeight() );
	hPtHat->Fill( jet->ptHat, dijet->GetWeight() );

	hRespMeasAbs.at(0)->Fill( jet->pt / jet->genPt, dijet->GetWeight() );
	hRespMeasAbs.at(1)->Fill( jet->pt / jet->genPt, dijet->GetWeight() );

	hRespMeas.at(0)->Fill( jet->pt / jet->genPt, dijet->GetWeight() );
	hRespMeas.at(1)->Fill( jet->pt / jet->genPt, dijet->GetWeight() );

	for(int i = 0; i < nPtDijetBins; i++) {
	  if( ptDijetBinEdges.at(i) <= dijet->dijetPt() && dijet->dijetPt() < ptDijetBinEdges.at(i+1) ) {
	    hRespMeasAbs.at(i)->Fill( jet->pt / jet->genPt, dijet->GetWeight() );
	    hRespMeas.at(i)->Fill( jet->pt / jet->genPt, dijet->GetWeight() );
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
  std::vector<double> meanRespPar = bag_of<double>(config_->read<string>("mean response parameters","1 0"));

  // Find bin centers for response function evaluation
  std::vector<double> ptGenBinCenters(nPlotBins);
  ptGenBinCenters.at(0) = 0.5 * ( tMin + tMax );
  ptGenBinCenters.at(1) = hPtGen->GetMean();
  for(int i = 2; i < nPlotBins; i++) {
    int j = i - 2;
    ptGenBinCenters.at(i) = 0.5 * ( ptDijetBinEdges.at(j) + ptDijetBinEdges.at(j+1) );
  }


  SmearData * smearData = dynamic_cast<SmearData*>(data_->front());
  if( smearData ) {
    // Loop over plotBins
    for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
      // Interpolated response function
      for(int bin = 1; bin <= hRespFit.at(plotBin)->GetNbinsX(); bin++) {
	double r = hRespFit.at(plotBin)->GetBinCenter(bin);
	hRespFit.at(plotBin)->SetBinContent(bin,smearData->respPDF(r,ptGenBinCenters.at(plotBin)));
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
	  double val  = scale.at(3+bin)*(smearData->GetRespPar()[3+bin]);
	  hRespFitStep.at(plotBin)->SetBinContent(bin,val);
	}
	normHist(hRespFitStep.at(plotBin),"width");
	hRespFitStep.at(plotBin)->Scale(1. - scale.at(0)*(smearData->GetRespPar()[0]));
	
	// Gauss part of fit function
	for(int bin = 1; bin <= hRespFitGaus.at(plotBin)->GetNbinsX(); bin++) {
	  // Mean
	  double a1 = meanRespPar.at(0);
	  double a2 = meanRespPar.at(1);
	  double mu = a1 + a2*ptGenBinCenters.at(plotBin);
	  // Width
	  a1 = scale.at(1)*smearData->GetRespPar()[1];
	  a2 = scale.at(2)*smearData->GetRespPar()[2];
	  double a3 = scale.at(3)*smearData->GetRespPar()[3];
	  double sigma = sqrt( a1*a1/ptGenBinCenters.at(plotBin)/ptGenBinCenters.at(plotBin)
			       + a2*a2/ptGenBinCenters.at(plotBin) + a3*a3 );
	  // pdf
	  double c     = scale.at(0)*(smearData->GetRespPar()[0]);
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
      } else {
	std::cout << "WARNING: No controlplots implemented for parametrization '" << param << "'\n";
      }
    } // End of loop over plotBins
  } // End if( smearData )


  // --- Fill histograms of fitted truth spectrum -----------

  // Fill histogram of assumed truth pdf
  // and fit with 1/x^n function
  std::cout << "Filling truth pdf\n";
  DataIt datait = data_->begin();
  while( (*datait)->GetType() != TypeSmearDiJet  &&  datait != data_->end() ) datait++;
  if( datait != data_->end() ) {
    SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

    for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
      double t = hTruthPDF->GetBinCenter(bin);
      hTruthPDF->SetBinContent(bin,dijet->truthPDF(t));
    }
    normHist(hTruthPDF,"width");
  }


  // --- Find populated x-axis ranges -----------------------
  int maxBin = 0;
  int minBin = 1;
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    int maxBinTmp = 0;
    int minBinTmp = 1;
    for(int bin = 1; bin <= hRespMeas.at(plotBin)->GetNbinsX(); bin++) {
      if( hRespMeas.at(plotBin)->GetBinContent(bin) > 0 ) maxBinTmp = bin;
      if( minBinTmp > maxBinTmp ) minBinTmp = bin;
    }
    if( maxBinTmp < hRespMeas.at(plotBin)->GetNbinsX() ) maxBinTmp++;
    if( minBinTmp < minBin ) minBin = minBinTmp;
    if( maxBinTmp > maxBin ) maxBin = maxBinTmp;
  }
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    hRespMeas.at(plotBin)->GetXaxis()->SetRange(minBin,maxBin);
    hRespMeasAbs.at(plotBin)->GetXaxis()->SetRange(minBin,maxBin);
  }
  
  maxBin = 0;
  minBin = 1;
  for(int bin = 1; bin <= hPtGen->GetNbinsX(); bin++) {
    if( hPtGen->GetBinContent(bin) > 0 ) maxBin = bin;
    if( minBin > maxBin ) minBin = bin;
  }
  if( maxBin < hPtGen->GetNbinsX() ) maxBin++;
  hPtGen->GetXaxis()->SetRange(minBin,maxBin);
  hPtHat->GetXaxis()->SetRange(minBin,maxBin);
  

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

    std::string label = toString(ptDijetMin) + " < p^{dijet}_{T} < " + toString(ptDijetMax) + " GeV";
    if( plotBin < 2 ) {
      legPtRange.at(plotBin)->AddEntry(hRespMeas.at(plotBin),label.c_str(),"L");
      legPtRangeAndCenters.at(plotBin)->AddEntry(hRespMeas.at(plotBin),label.c_str(),"L");
    } else {
      int i = plotBin - 2;
      label = toString(ptDijetBinEdges.at(i)) + " < p^{dijet}_{T} < " + toString(ptDijetBinEdges.at(1+i)) + " GeV";
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
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"RECREATE");
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    rootfile.WriteTObject(hRespMeasAbs.at(plotBin));
    rootfile.WriteTObject(hRespMeas.at(plotBin));
    rootfile.WriteTObject(hRespFit.at(plotBin));
    rootfile.WriteTObject(hRespFitStart.at(plotBin));
    rootfile.WriteTObject(hRespFitStep.at(plotBin));
    rootfile.WriteTObject(hRespFitGaus.at(plotBin));
    rootfile.WriteTObject(hRespFitSum.at(plotBin));
  }
  rootfile.WriteTObject(hPtGen);
  rootfile.WriteTObject(hPtHat);
  rootfile.WriteTObject(hPtDijet);
  rootfile.WriteTObject(hTruthPDF);

  rootfile.Close();


  // --- Clean up ------------------------------------------
  for(int plotBin = 0; plotBin < nPlotBins; plotBin++) {
    delete hRespMeasAbs.at(plotBin);
    delete hRespMeas.at(plotBin);
    delete hRespFit.at(plotBin);
    delete hRespFitStart.at(plotBin);
    delete hRespFitStep.at(plotBin);
    delete hRespFitGaus.at(plotBin);
    delete hRespFitSum.at(plotBin);
    delete legPtRangeAndCenters.at(plotBin);
  }
  delete legFitStart;
  delete hPtGen;
  delete hPtHat;
  delete hPtDijet;
  delete hTruthPDF;
  delete c1;
  delete ps;
}



//!  \brief Plot the negative log-likelihood for different
//!         parameter values around the fitted value
//!  \param pars Indices of the parameters for which a
//!              scan plot is done
// --------------------------------------------------
void ControlPlotsJetSmearing::plotParameterScan(const vector<unsigned int>& pars) const {
  std::cout << "Creating parameter scan control plots\n";

  // Create one histogram of likelihood per parameter
  vector<TH1F*> hLikelihood;
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    char name[50];
    sprintf(name,"hNLogL%i",i);
    char title[50];
    sprintf(title,";Parameter %i;-ln(L)",i);
    double min = param_->GetPars()[i] - 0.5;
    if( min < 0. ) min = 0.;
    double max = param_->GetPars()[i] + 0.5;
    TH1F * h = new TH1F(name,title,10,min,max);
    h->SetMarkerStyle(24);
    hLikelihood.push_back(h);
  }

  // Do parameter scan
  for(unsigned int i = 0; i < pars.size(); i++) {
    int idx = pars.at(i);
    if( idx < param_->GetNumberOfParameters() ) {
      std::cout << "  Par " << idx << "... " << std::flush;
      double oldpar = param_->GetPars()[idx];
      TH1F * h = hLikelihood.at(idx);
      for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
	param_->GetPars()[idx] = h->GetBinCenter(bin);
	double fsum = 0.;
	for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
	  fsum += (*datait)->chi2();
	}
	h->SetBinContent(bin,fsum);
      }
      param_->GetPars()[idx] = oldpar;
      std::cout << "ok" << std::endl;
    }
  }

  // Plot histograms
  TPostScript * const ps = new TPostScript((dir_+"/jsParameterScan.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Parameter Scan",0,0,600,600);
  ps->NewPage();
  c1->cd();
  for(unsigned int i = 0; i < pars.size(); i++) {
    int idx = pars.at(i);
    if( idx < param_->GetNumberOfParameters() ) {
      TH1F * h = hLikelihood.at(idx);
      h->Draw("P");
      TLine * line = new TLine(param_->GetPars()[idx],h->GetMinimum(),
			       param_->GetPars()[idx],h->GetMaximum());
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw("same");
      c1->Draw();
      ps->NewPage();
      delete line;
    }
  }
  ps->Close();

  // Clean up
  for(vector<TH1F*>::iterator it = hLikelihood.begin(); it != hLikelihood.end(); it++) {
    delete *it;
  }
  hLikelihood.clear();
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
	TJet         * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 )   jet = static_cast<TJet*>(dijet->GetSecondMess());

	if( jet->genPt < minGenJetPt ) minGenJetPt = jet->genPt;
	if( jet->genPt > maxGenJetPt ) maxGenJetPt = jet->genPt;

	if( jet->pt < minCalJetPt ) minCalJetPt = jet->pt;
	if( jet->pt > maxCalJetPt ) maxCalJetPt = jet->pt;
      }

      if( dijet->ptHat() < minPtHat ) minPtHat = dijet->ptHat();
      if( dijet->ptHat() > maxPtHat ) maxPtHat = dijet->ptHat();

      if( dijet->dijetPt() < minDijetPt ) minDijetPt = dijet->dijetPt();
      if( dijet->dijetPt() > maxDijetPt ) maxDijetPt = dijet->dijetPt();

      TJet * jet3 = static_cast<TJet*>(dijet->GetThirdMess());
      if( jet3->pt < min3rdJetPt ) min3rdJetPt = jet3->pt;
      if( jet3->pt > max3rdJetPt ) max3rdJetPt = jet3->pt;
    }
  }

  // Pt distributions
  TH1F * hPtHat = new TH1F("hPtHat",";#hat{p}_{T} (GeV);1 / N  dN / d#hat{p}_{T}  1 / (GeV)",
			   50,0.9*minPtHat, 1.1*maxPtHat);
  hPtHat->SetLineWidth(2);
  hPtHat->Sumw2();

  vector<TH1F*> hGenJetPt;
  vector<TH1F*> hCalJetPt;
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

      TJet * jet1 = static_cast<TJet*>(dijet->GetMess());
      TJet * jet2 = static_cast<TJet*>(dijet->GetSecondMess());
      TJet * jet3 = static_cast<TJet*>(dijet->GetThirdMess());

      double weight = dijet->GetWeight();

      double dPhi = std::abs(TVector2::Phi_mpi_pi( jet1->phi - jet2->phi ));
      hDeltaPhi->Fill( dPhi, weight );

      hRvsDeltaPhi->Fill( dPhi, jet1->pt / jet1->genPt, weight );
      hRvsDeltaPhi->Fill( dPhi, jet2->pt / jet2->genPt, weight );

      hRvsEMF->Fill( jet1->EMF / jet1->pt, jet1->pt / jet1->genPt, weight );
      hRvsEMF->Fill( jet2->EMF / jet2->pt, jet2->pt / jet2->genPt, weight );

      hRvsRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt(), jet1->pt / jet1->genPt, weight );
      hRvsRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt(), jet2->pt / jet2->genPt, weight );

      hRvs3rdJetPt->Fill( jet3->pt, jet1->pt / jet1->genPt, weight );
      hRvs3rdJetPt->Fill( jet3->pt, jet2->pt / jet2->genPt, weight );

      hRvsDeltaR->Fill( jet1->dR, jet1->pt / jet1->genPt, weight );
      hRvsDeltaR->Fill( jet2->dR, jet2->pt / jet2->genPt, weight );

      hPtHat->Fill( dijet->ptHat(), weight );

      hGenJetPt.at(0)->Fill( jet1->genPt, weight );
      hGenJetPt.at(0)->Fill( jet2->genPt, weight );
      hGenJetPt.at(1)->Fill( jet1->genPt, weight );
      hGenJetPt.at(2)->Fill( jet2->genPt, weight );

      hCalJetPt.at(0)->Fill( jet1->pt, weight );
      hCalJetPt.at(0)->Fill( jet2->pt, weight );
      hCalJetPt.at(1)->Fill( jet1->pt, weight );
      hCalJetPt.at(2)->Fill( jet2->pt, weight );

      hCalJet2vsCalJet1Pt->Fill( jet1->pt, jet2->pt, weight );
      hDijetPt->Fill( dijet->dijetPt(), weight );

      h3rdJetPt->Fill( jet3->pt, weight );
      h3rdJetvsDijetPt->Fill( dijet->dijetPt(), jet3->pt, weight );
      hRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt(), weight );
      hDeltaPhivsRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt(), dPhi, weight );
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

  vector<TObject*> objs;
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



//!  \brief Plot the mean response and resolution vs ptGen
//!
//!  The mean is determined as mean of the distribution
//!  and as mean of a Gauss fit. The relative and
//!  absolute resolution are determined as RMS and width
//!  of a Gauss fit.
//!
//!  The response is fitted with a linear function of ptGen,
//!  the resolution with the usual calorimeter resolution
//!  parametrization.
// --------------------------------------------------
void ControlPlotsJetSmearing::plotMeanResponseAndResolution() const {
  std::cout << "Creating mean response and resolution control plots\n";

  // Find pt ranges
  double minCalJetPt = 10000.;
  double maxCalJetPt = 0.;
  double minGenJetPt = 10000.;
  double maxGenJetPt = 0.;
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      // Loop over both jets
      for(int i = 0; i < 2; i++) {        
	TJet         * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 )   jet = static_cast<TJet*>(dijet->GetSecondMess());

	if( jet->genPt < minGenJetPt ) minGenJetPt = jet->genPt;
	if( jet->genPt > maxGenJetPt ) maxGenJetPt = jet->genPt;

	if( jet->pt < minCalJetPt ) minCalJetPt = jet->pt;
	if( jet->pt > maxCalJetPt ) maxCalJetPt = jet->pt;
      }
    }
  }


  // Create histograms
  TH2F * hRespVsPtGen = new TH2F("hRespVsPtGen",
				 ";p^{gen}_{T} (GeV);p^{jet}_{T} / p^{gen}_{T}",
				 25,0.9*minGenJetPt,1.1*maxGenJetPt,51,0,2);
  hRespVsPtGen->SetNdivisions(505);
  hRespVsPtGen->Sumw2();

  TH2F * hPtJetVsPtGenJet = new TH2F("hPtJetVsPtGenJet",
			       ";p^{gen}_{T} (GeV);p^{jet}_{T}",
			       25,0.9*minGenJetPt,1.1*maxGenJetPt,
			       25,0.9*minCalJetPt,1.1*maxCalJetPt);
  hPtJetVsPtGenJet->SetNdivisions(505);
  hPtJetVsPtGenJet->Sumw2();

  // Fill histograms
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      // Loop over both jets
      for(int i = 0; i < 2; i++) {        
	TJet         * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 )   jet = static_cast<TJet*>(dijet->GetSecondMess());

	hRespVsPtGen->Fill( jet->genPt, jet->pt / jet->genPt, dijet->GetWeight() );
	hPtJetVsPtGenJet->Fill( jet->genPt, jet->pt, dijet->GetWeight() );
      }
    }
  }


  // Get mean response and resolution vs ptgen
  std::vector<TH1F*> hRes;
  fitSlices(hRespVsPtGen,hRes);

  std::vector<TH1F*> hPtJet;
  fitSlices(hPtJetVsPtGenJet,hPtJet);


  // Fit mean response
  std::vector<TF1*> fResp(4);
  for(int i = 0; i < 4; i++) {
    TF1 * f = 0;
    if( i%2 == 0 ) {
      f = new TF1(("fResp"+toString(i)).c_str(),"pol1",1.1*minGenJetPt,0.9*maxGenJetPt);
      f->SetParameter(0,1);
      f->SetParameter(1,0);
    } else {
      f = new TF1(("fResp"+toString(i)).c_str(),"sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		  1.1*minGenJetPt,0.9*maxGenJetPt);
      f->SetParameter(0,4.4);
      f->SetParameter(1,1.1);
      f->SetParameter(2,0.0);
    }
    f->SetLineWidth(2);
    f->SetLineColor(2);
    hRes.at(i)->Fit(f,"0QLR");
    fResp.at(i) = f;
  }


  // Fit resolution
  std::vector<TF1*> fReso(2);
  for(int i = 0; i < 2; i++) {
    TF1 * f = new TF1(("fReso"+toString(i)).c_str(),"sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		1.1*minGenJetPt,0.9*maxGenJetPt);
    f->SetParameter(0,4.4);
    f->SetParameter(1,1.1);
    f->SetParameter(2,0.0);
    f->SetLineWidth(2);
    f->SetLineColor(2);
    hPtJet.at(2*i+1)->Fit(f,"0QLR");
    fReso.at(i) = f;
  }


  // Draw mean response and resolution vs ptgen
  TPostScript * const ps = new TPostScript((dir_+"/jsMeanResponse.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","Mean response",0,0,600,600);

  drawPSPage(ps,c1,hRespVsPtGen,"COLZ");
  drawPSPage(ps,c1,hPtJetVsPtGenJet,"COLZ");

  std::vector<TObject*> objs;
  for(size_t i = 0; i < hRes.size(); i++) {
    TLine * line = 0;
    TLegend * fitstat = new TLegend(0.21,0.7,0.81,0.8);
    fitstat->SetBorderSize(0);
    fitstat->SetFillColor(0);
    fitstat->SetTextFont(42);
    fitstat->SetTextAlign(12);

    hRes.at(i)->SetNdivisions(505);

    objs.clear();
    objs.push_back(hRes.at(i));
    objs.push_back(fResp.at(i));

    if( i%2 == 0 ) {
      line = new TLine(hRes.at(i)->GetXaxis()->GetXmin(),1.,
		       hRes.at(i)->GetXaxis()->GetXmax(),1.);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->SetLineColor(4);
      objs.push_back(line);

      char label[100];
      sprintf(label,"< p^{jet}_{T} / p^{gen}_{T} > = %.3f %.0f #upoint10^{-6} p^{gen}_{T}",
	      fResp.at(i)->GetParameter(0),
	      1E6*(fResp.at(i)->GetParameter(1)));
      fitstat->AddEntry(fResp.at(i),label,"L");
    } else {
      char label[150];
      sprintf(label,"#frac{#sigma(p^{jet}_{T}/p^{gen}_{T})}{<p^{jet}_{T}/p^{gen}_{T}>} = #frac{%.2f}{p^{gen}_{T}} #oplus #frac{%.2f}{#sqrt{p^{gen}_{T}}} #oplus %.2f",
	      fResp.at(i)->GetParameter(0),
	      fResp.at(i)->GetParameter(1),
	      fResp.at(i)->GetParameter(2));
      fitstat->AddEntry(fResp.at(i),label,"L");
    }
    objs.push_back(fitstat);
    drawPSPage(ps,c1,objs,"PE");

    if( line ) delete line;
    delete fitstat;
  }

  int j = 0;
  for(size_t i = 1; i < hPtJet.size(); i += 2, j++) {
    hPtJet.at(i)->SetNdivisions(505);

    TLegend * fitstat = new TLegend(0.21,0.7,0.81,0.8);
    fitstat->SetBorderSize(0);
    fitstat->SetFillColor(0);
    fitstat->SetTextFont(42);
    fitstat->SetTextAlign(12);
    char label[100];
    sprintf(label,"#frac{#sigma(p^{jet}_{T})}{<p^{jet}_{T}>} = #frac{%.1f}{p^{gen}_{T}} #oplus #frac{%.1f}{#sqrt{p^{gen}_{T}}} #oplus %.1f",
	    fReso.at(j)->GetParameter(0),
	    fReso.at(j)->GetParameter(1),
	    fReso.at(j)->GetParameter(2));
    fitstat->AddEntry(fReso.at(j),label,"L");

    objs.clear();
    objs.push_back(hPtJet.at(i));
    objs.push_back(fReso.at(j));
    objs.push_back(fitstat);

    drawPSPage(ps,c1,objs,"PE");
    delete fitstat;
  }


  // Clean up
  for(size_t i = 0; i < hRes.size(); i++) {
    delete hRes.at(i);
    delete hPtJet.at(i);
    delete fResp.at(i);
    if( i < fReso.size() ) delete fReso.at(i);
  }
  delete hRespVsPtGen;
  delete hPtJetVsPtGenJet;
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
  vector<TObject*> objs;
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
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, vector<TObject*> objs, std::string option, bool log) const {
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




//!  \brief Get different x-profiles a 2D histogram
//!
//!  Calculates the 1D projections of \p h2 along the
//!  x-axis for different x bins. Different quantities
//!  of these projections are calculated and shown
//!  versus x in the 1D histograms \p hFit
//!   - 0: Mean
//!   - 1: RMS
//!   - 2: Mean of central (\f$ \pm3\sigma \f$) Gauss fit
//!   - 3: Width of central Gauss fit
//!
//!  The binning and the x-axis title of the histograms
//!  in \p hFit is the same as the x binning of \p h2, the
//!  y-axis title and the histogram title are adjusted
//!  to the displayed quantity.
// --------------------------------------------------
void ControlPlotsJetSmearing::fitSlices(const TH2F * h2, std::vector<TH1F*>& hFit) const {
  // Reset result vector
  const int nResHist = 4;
  hFit.clear();
  hFit = std::vector<TH1F*>(nResHist);

  // Create result histograms
  std::string name   = h2->GetName();
  std::string xTitle = h2->GetXaxis()->GetTitle();
  std::string yTitle = h2->GetYaxis()->GetTitle();
  std::string quant[nResHist] = { "Mean", "RMS", "GaussMean", "GaussSigma" };
  std::string title[nResHist] = { "Mean", "RMS", "Mean of Gauss fit", "#sigma of Gauss fit" };
  
  for(int i = 0; i < nResHist; i++) {
    TH1F * h = 0;
    if( i % 2 == 0 ) {
      h = new TH1F((name+quant[i]).c_str(),
		   (title[i]+";"+xTitle+";< "+yTitle+" >").c_str(),
		   h2->GetNbinsX(),h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax());
      h->GetYaxis()->SetRangeUser(0.95,1.05);
    } else {
      h = new TH1F((name+quant[i]).c_str(),
		   (title[i]+";"+xTitle+";#sigma("+yTitle+") / < "+yTitle+" >").c_str(),
		   h2->GetNbinsX(),h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax());
      h->GetYaxis()->SetRangeUser(0.,0.15);
    }
    h->SetMarkerStyle(20);
    h->SetLineWidth(2);
    h->SetMarkerSize(1.2);
    h->Sumw2();
    hFit.at(i) = h;
  }

  // Get 1D slices and get mean, sigma etc
  TH1F* hSlice = new TH1F("hSlice","",h2->GetNbinsY(),h2->GetYaxis()->GetXmin(),
  			   h2->GetYaxis()->GetXmax());
  hSlice->Sumw2();
  for(int xBin = 1; xBin <= h2->GetNbinsX(); xBin++) {
    hSlice->Reset();
    for(int yBin = 1; yBin <= h2->GetNbinsY(); yBin++) {
      hSlice->SetBinContent(yBin,h2->GetBinContent(h2->GetBin(xBin,yBin)));
      hSlice->SetBinError(yBin,h2->GetBinError(h2->GetBin(xBin,yBin)));
    }  

    double mean       = hSlice->GetMean();
    double meanError  = hSlice->GetMeanError();
    double width      = hSlice->GetRMS();
    double widthError = hSlice->GetRMSError();

    hFit.at(0)->SetBinContent(xBin,mean);
    hFit.at(0)->SetBinError(xBin,meanError);
    if( mean ) {
      hFit.at(1)->SetBinContent(xBin,width/mean);
      hFit.at(1)->SetBinError(xBin,widthError/mean);
    }

    if( width < 0.1 ) width = 0.1;
    hSlice->Fit("gaus","QN","",mean-3*width,mean+3*width);
    TF1 * f = static_cast<TF1*>(gROOT->GetFunction("gaus"));
    mean       = f->GetParameter(1);
    meanError  = f->GetParError(1);
    width      = f->GetParameter(2);
    widthError = f->GetParError(2);

    hFit.at(2)->SetBinContent(xBin,mean);
    hFit.at(2)->SetBinError(xBin,meanError);
    if( mean ) {
      hFit.at(3)->SetBinContent(xBin,width/mean);
      hFit.at(3)->SetBinError(xBin,widthError/mean);
    }
  }
  delete hSlice;
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
