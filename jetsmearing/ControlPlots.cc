// $Id: ControlPlots.cc,v 1.4 2009/05/06 12:17:32 mschrode Exp $

#include "ControlPlots.h"

#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TPostScript.h"
#include "TStyle.h"

#include "Event.h"
#include "NJetEvent.h"
#include "PhotonJetEvent.h"
#include "HistOps.h"



namespace js
{
  // --------------------------------------------------
  ControlPlots::ControlPlots(const Data& data)
    : mDijetNBins(100), mDijetMin(50), mDijetMax(1000),
      mPhotonJetNBins(100), mPhotonJetMin(50), mPhotonJetMax(1000),
      mRespNBins(30), mRespMin(0.), mRespMax(6.),
      mDir("./controlplots"), mRootFileName("Plots.root"),
      mFileNameSuffix("")
  {
    // Copy pointers to events over to local
    // Data object
    for(DataIt it = data.begin(); it != data.end(); it++)
      {
	mData.push_back(*it);
      }

    SetGStyle();
  }



  // --------------------------------------------------
  ControlPlots::~ControlPlots() {;}



  // --------------------------------------------------
  void ControlPlots::PlotDijets() const
  {
    std::cout << "Creating dijet control plots... " << std::flush;

    // Create histograms
    TH1F * hDijetPtTrue       = new TH1F("hDijetPtTrue",";p^{true}_{T} (GeV)",mDijetNBins,mDijetMin,mDijetMax);
    TH1F * hDijetPtMeas       = new TH1F("hDijetPtMeas",";p^{jet}_{T} (GeV)",mDijetNBins,mDijetMin,mDijetMax);
    TH1F * hDijetDeltaPtMeas  = new TH1F("hDijetDeltaPtMeas",
					 ";(p^{jet}_{T,0} - p^{jet}_{T,1}) / (p^{jet}_{T,0} + p^{jet}_{T,1})",
					 mDijetNBins,-0.5,0.5);
    TH1F * hDijetDeltaPtTrue  = new TH1F("hDijetDeltaPtTrue",
					 ";(p^{true}_{T,0} - p^{true}_{T,1}) / (p^{true}_{T,0} + p^{true}_{T,1})",
					 mDijetNBins,-0.5,0.5);
    TH1F * hDijetRelPt        = new TH1F("hDijetRelPt",";p^{jet}_{T} / p^{true}_{T}",mDijetNBins,0,2);

    TH1F * hDijetPhiTrue      = new TH1F("hDijetPhiTrue",";#phi^{true}",50,-M_PI,M_PI);
    TH1F * hDijetPhiMeas      = new TH1F("hDijetPhiMeas",";#phi^{jet}",50,-M_PI,M_PI);
    TH1F * hDijetDeltaPhiMeas = new TH1F("hDijetDeltaPhiMeas",";#Delta#phi^{meas}",50,M_PI-1,M_PI+1);
    TH1F * hDijetDeltaPhiTrue = new TH1F("hDijetDeltaPhiTrue",";#Delta#phi^{true}",50,M_PI-1,M_PI+1);
    TH1F * hDijetRelPhi       = new TH1F("hDijetRelPhi",";#phi^{jet} / #phi^{true}",50,0.5,1.5);
    

    // Fill histograms
    for(DataIt datait = mData.begin(); datait != mData.end(); datait++) // Loop over Data
      {
	// Select DiJet events
	if( (*datait)->Type() == "DiJetEvent" )
	  {
	    const DiJetEvent * njevt = static_cast<const DiJetEvent*>(*datait);

	    // Loop over both jets
	    JetIt jetit = njevt->Begin();
	    for(; jetit != njevt->End(); jetit++)
	      {
		double t = (*jetit)->PtTrue();
		double m = (*jetit)->PtMeas();

		hDijetPtTrue->Fill(t);
		hDijetPtMeas->Fill(m);
		hDijetRelPt->Fill(m/t);		    

		m = (*jetit)->PhiMeas();
		t = (*jetit)->PhiTrue();
		    
		hDijetPhiTrue->Fill(t);
		hDijetPhiMeas->Fill(m);
		hDijetRelPhi->Fill(fabs(m/t));
	      }

	    double dif = njevt->PtTrue(0) - njevt->PtTrue(1);
	    double sum = njevt->PtTrue(0) + njevt->PtTrue(1);
	    hDijetDeltaPtTrue->Fill( dif/sum );

	    dif = njevt->PtMeas(0) - njevt->PtMeas(1);
	    sum = njevt->PtMeas(0) + njevt->PtMeas(1);
	    hDijetDeltaPtMeas->Fill( dif/sum );

	    dif = njevt->PhiTrue(0) - njevt->PhiTrue(1);
	    hDijetDeltaPhiTrue->Fill( dif );

	    dif = njevt->PhiMeas(0) - njevt->PhiMeas(1);
	    hDijetDeltaPhiMeas->Fill( dif );
	  }
      } // End loop over Data


    // Write histos to ps file
    hDijetPhiTrue->GetYaxis()->SetRangeUser(0,1.2*hDijetPhiTrue->GetMaximum());
    hDijetPhiMeas->GetYaxis()->SetRangeUser(0,1.2*hDijetPhiMeas->GetMaximum());

    std::vector<TH1F*> hists;
    hists.push_back(hDijetPtTrue);
    hists.push_back(hDijetPtMeas);
    hists.push_back(hDijetDeltaPtMeas);
    hists.push_back(hDijetDeltaPtTrue);
    hists.push_back(hDijetRelPt);
    hists.push_back(hDijetPhiTrue);
    hists.push_back(hDijetPhiMeas);
    hists.push_back(hDijetDeltaPhiMeas);
    hists.push_back(hDijetDeltaPhiTrue);
    hists.push_back(hDijetRelPhi);

    TPostScript * const ps = new TPostScript((mDir+"/js_"+mFileNameSuffix+"_Dijets.ps").c_str(),111);
    TCanvas *c1 = new TCanvas("c1","Dijets",0,0,600,600);
    c1->cd();

    std::vector<TH1F*>::iterator it = hists.begin();
    for( ; it != hists.end(); it++)
      {
	(*it)->Draw();
	c1->Draw();
	ps->NewPage();
      }
    ps->Close();

    // Write histos to root file
    util::HistOps::WriteToRootFile(hists,mDir+"/js_"+mFileNameSuffix+"_"+mRootFileName);

    // Clean up
    for( ; it != hists.end(); it++)
      {
	delete *it;
      }
    hists.clear();
    delete ps;
    delete c1;


    std::cout << "ok" << std::endl;
  }



  // --------------------------------------------------
  void ControlPlots::PlotPhotonJets() const
  {
    std::cout << "Creating photon-jet control plots... " << std::flush;

    // Create histograms
    TH1F * hPhotonJetPtTrue       = new TH1F("hPhotonJetPtTrue",";p^{#gamma}_{T} (GeV)",mPhotonJetNBins,mPhotonJetMin,mPhotonJetMax);
    TH1F * hPhotonJetPtMeas       = new TH1F("hPhotonJetPtMeas",";p^{jet}_{T} (GeV)",mPhotonJetNBins,mPhotonJetMin,mPhotonJetMax);
    TH1F * hPhotonJetRelPt        = new TH1F("hPhotonJetRelPt",";p^{jet}_{T} / p^{#gamma}_{T}",mPhotonJetNBins,0,2);

    TH1F * hPhotonJetPhiTrue      = new TH1F("hPhotonJetPhiTrue",";#phi^{true}",50,-M_PI,M_PI);
    TH1F * hPhotonJetPhiMeas      = new TH1F("hPhotonJetPhiMeas",";#phi^{jet}",50,-M_PI,M_PI);
    TH1F * hPhotonJetRelPhi       = new TH1F("hPhotonJetRelPhi",";#phi^{jet} / #phi^{true}",50,0.5,1.5);
    

    // Fill histograms
    for(DataIt datait = mData.begin(); datait != mData.end(); datait++) // Loop over Data
      {
	// Select DiJet events
	if( (*datait)->Type() == "PhotonJetEvent" )
	  {
	    const PhotonJetEvent * evt = static_cast<const PhotonJetEvent*>(*datait);

	    double t = evt->PtPhoton();
	    double m = evt->PtMeas();

	    hPhotonJetPtTrue->Fill(t);
	    hPhotonJetPtMeas->Fill(m);
	    hPhotonJetRelPt->Fill(m/t);		    

	    m = evt->PhiMeas();
	    t = evt->PhiTrue();
		    
	    hPhotonJetPhiTrue->Fill(t);
	    hPhotonJetPhiMeas->Fill(m);
	    hPhotonJetRelPhi->Fill(fabs(m/t));
	  }
      } // End loop over Data


    // Write histos to ps file
    hPhotonJetPhiTrue->GetYaxis()->SetRangeUser(0,1.2*hPhotonJetPhiTrue->GetMaximum());
    hPhotonJetPhiMeas->GetYaxis()->SetRangeUser(0,1.2*hPhotonJetPhiMeas->GetMaximum());

    std::vector<TH1F*> hists;
    hists.push_back(hPhotonJetPtTrue);
    hists.push_back(hPhotonJetPtMeas);
    hists.push_back(hPhotonJetRelPt);
    hists.push_back(hPhotonJetPhiTrue);
    hists.push_back(hPhotonJetPhiMeas);
    hists.push_back(hPhotonJetRelPhi);

    TPostScript * const ps = new TPostScript((mDir+"/js_"+mFileNameSuffix+"_PhotonJets.ps").c_str(),111);
    TCanvas *c1 = new TCanvas("c1","PhotonJets",0,0,600,600);
    c1->cd();

    std::vector<TH1F*>::iterator it = hists.begin();
    for( ; it != hists.end(); it++)
      {
	(*it)->Draw();
	c1->Draw();
	ps->NewPage();
      }
    ps->Close();

    // Write histos to root file
    util::HistOps::WriteToRootFile(hists,mDir+"/js_"+mFileNameSuffix+"_"+mRootFileName);

    // Clean up
    for( ; it != hists.end(); it++)
      {
	delete *it;
      }
    hists.clear();
    delete ps;
    delete c1;


    std::cout << "ok" << std::endl;
  }



  // --------------------------------------------------
  void ControlPlots::PlotResponse(TObject * pdf) const
  {
    std::cout << "Creating response control plots... " << std::flush;

    // Create histograms
    TH1F* hRespMeas = new TH1F("hRespMeas",";p^{jet}_{T} / p^{true}_{T}",mRespNBins,mRespMin,mRespMax);
    hRespMeas->Sumw2();

    //    TH1F* hRelDiff = new TH1F("hRelDiff",";( p^{jet}_{T} - p^{true}_{T} ) / p^{true}_{T}",mDiffNBins,mDiffMin,mDiffMax);
    //hRespMeas->Sumw2();


    // Fill histograms
    for(DataIt datait = mData.begin(); datait != mData.end(); datait++)
      {
	// NJet events
	if( (*datait)->Type() == "NJetEvent"  ||  (*datait)->Type() == "DiJetEvent" )
	  {
	    const NJetEvent *njevt = static_cast<const NJetEvent*>(*datait);

	    // Loop over jets in this event
	    JetIt jetit = njevt->Begin();
	    for(; jetit != njevt->End(); jetit++)
	      {
		double m = (*jetit)->PtMeas();
		double t = (*jetit)->PtTrue();
		hRespMeas->Fill(m/t);
	      }
	  }

	// PhotonJet events
	if( (*datait)->Type() == "PhotonJetEvent" )
	  {
	    const PhotonJetEvent *evt = static_cast<const PhotonJetEvent*>(*datait);

	    double m = evt->PtMeas();
	    double t = evt->PtTrue();
	    hRespMeas->Fill(m/t);
	  }
      }

    // Normalize histogram
    double norm = hRespMeas->Integral("width");
    hRespMeas->Scale(1./norm);

    // Find populated x-axis range
    int maxBin = 1;
    for(int bin = 1; bin <= hRespMeas->GetNbinsX(); bin++)
      {
	if( hRespMeas->GetBinContent(bin) > 0 ) maxBin = bin;
      }
    if( maxBin < hRespMeas->GetNbinsX() ) maxBin++;
    hRespMeas->GetXaxis()->SetRange(1,maxBin);

    // In case of histogramed pdf, fill a histogram
    TH1F * hRespFit = new TH1F("hRespFit","",10*mRespNBins,mRespMin,mRespMax);
    std::string type = pdf->ClassName();
    if( type == "TH1D" )
      {
	TH1D * hHistPDF = dynamic_cast<TH1D*>(pdf);
	for(int bin = 1; bin <= hRespFit->GetNbinsX(); bin++)
	  {
	    double r = hRespFit->GetBinCenter(bin);
	    hRespFit->SetBinContent(bin,hHistPDF->Interpolate(r));
	  }
	norm = hRespFit->Integral("width");
	hRespFit->Scale(1./norm);
      }

    // Write histos to eps file
    TCanvas *c1 = new TCanvas("c1","Jet Response",0,0,600,600);
    c1->cd();
    hRespMeas->Draw();
    if( type == "TH1D" ) hRespFit->Draw("Lsame");
    else pdf->Draw("same");
    //    c1->SetGrid();
    c1->SetLogy();
    c1->SaveAs((mDir+"/js_"+mFileNameSuffix+"_JetResponse.eps").c_str());
    delete c1;

    // Write histos to root file
    std::vector<TH1F*> hists;
    hists.push_back(hRespMeas);
    util::HistOps::WriteToRootFile(hists,mDir+"/js_"+mFileNameSuffix+"_"+mRootFileName);

    // Clean up
    delete hRespMeas;
    delete hRespFit;
    delete pdf;

    std::cout << "ok" << std::endl;
  }



  // --------------------------------------------------
  void ControlPlots::SetGStyle() const
  {
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

    // Margins:
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadRightMargin(0.04);

    // For the Global title:
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFont(42,"");
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.58);
    gStyle->SetTitleH(0.05);
    gStyle->SetTitleXOffset(0);
    gStyle->SetTitleYOffset(0);
    gStyle->SetTitleBorderSize(0);

    // For the axis titles:
    gStyle->SetTitleColor(1,"XYZ");
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetTitleSize(0.04,"XYZ");
    gStyle->SetTitleXOffset(1.5);
    gStyle->SetTitleYOffset(2.0);

    // For the axis labels:
    gStyle->SetLabelColor(1,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetLabelOffset(0.007,"XYZ");
    gStyle->SetLabelSize(0.04,"XYZ");

    // For the axis:
    gStyle->SetAxisColor(1,"XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03,"XYZ");
    gStyle->SetNdivisions(510,"XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
  }
}
