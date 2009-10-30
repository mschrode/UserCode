#include "ControlPlotsComparison.h"

#include <iostream>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TList.h>
#include <TPostScript.h>
#include <TString.h>
#include <TStyle.h>



//---------------------------------------------------------------
//  Constructor: Comparison of existing controlplots
//  using 'Compare(...)', which takes the root-files storing
//  the control plots to be compared as parameter. The output
//  is in .ps-format only.
//---------------------------------------------------------------
TControlPlotsComparison::TControlPlotsComparison()
{
  SetGStyle();
}


//---------------------------------------------------------------
//  Compar controlplots stored in root file. The output is in
//  .ps-format. The method takes a vector 'fileName' with the
//  names of the root-files to be compared. Optionally, a vector
//  'descr' with comments/descriptions to each root file can be
//  specified. The descriptions will be displayed in the legend
//  of the comparison plots. Default behaviour is to use the
//  file name as description of each file.
//---------------------------------------------------------------
void TControlPlotsComparison::CompareControlPlots(const std::vector<std::string> &fileName, std::vector<std::string> descr) const
{
  std::vector<TFile*> file;
  bool compPlotsGammaJet = true;
  bool compPlotsDiJet = true;
  bool descrExists = ( descr.size() == fileName.size() );
  if( !descrExists ) descr.clear();
  for(std::vector<std::string>::const_iterator it = fileName.begin(); it < fileName.end(); it++)
    {
      std::cout << "Opening file " << *it << "... " << std::flush;
      TFile *f = new TFile( (*it).c_str(),"READ");
      file.push_back( f );


      // Check which directories exist
      TString dirName = f->GetName();
      dirName += ":/GammaJet";
      compPlotsGammaJet = compPlotsGammaJet && f->GetDirectory( dirName );

      dirName = f->GetName();
      dirName += ":/DiJet";
      compPlotsDiJet = compPlotsDiJet && f->GetDirectory( dirName );


      // If no description exist, use file names
      if( !descrExists ) descr.push_back(*it);

      std::cout << "ok" << std::endl;
    }


  // Compare existing controlplots
  if( compPlotsGammaJet )
    {
      std::cout << "Comparing gamma jet control plots... " << std::flush;
      CompareControlPlotsGammaJet(file, descr);
      std::cout << "ok" << std::endl;
    }
  if( compPlotsDiJet )
    {
      std::cout << "Comparing dijet control plots... " << std::flush;
      CompareControlPlotsDiJet(file, descr);
      std::cout << "ok" << std::endl;
    }


  // Clean up
  for(std::vector<TFile*>::const_iterator it = file.begin(); it < file.end(); it++)
    {
      (*it)->Close();
    }
  file.clear();
}


//---------------------------------------------------------------
//   Compare DiJet Control Plots
//---------------------------------------------------------------
void TControlPlotsComparison::CompareControlPlotsDiJet(const std::vector<TFile*> &file, const std::vector<std::string> &descr) const
{
  int nFILES = static_cast<int>(file.size());

  // Declare histo pointers

  // Control quantities of B vs eta, pt
  TH1F **hists_eta[8][8];
  TH1F **hists_pt[2][8];
  for(int i = 0; i < 8; i++)
    {
      for(int j = 0; j < 8; j++)
	{
	  hists_eta[i][j] = new TH1F*[nFILES];
	  if ( i < 2 ) hists_pt[i][j] = new TH1F*[nFILES];
	}
    }

  // Get histos from file
  bool histosExist = true;
  for(int f = 0; f < nFILES; f++) // Loop over files
    {
      TString dirName = file.at(f)->GetName();
      dirName += ":/DiJet/";
      for(int i = 0; i < 8; i++)
	{
	  for(int j = 0; j < 8; j++)
	    {
	      TString name = dirName;
	      name += "hBeta";
	      name += i;
	      name += "_result";
	      name += j;
	      hists_eta[i][j][f] = static_cast<TH1F*>(file.at(f)->Get(name));
	      if( hists_eta[i][j][f] == 0 )
		{
		  histosExist = false;
		  std::cerr << "ERROR: " << name << " does not exist!" << std::endl;
		}

	      if( i < 2 )
		{
		  name = dirName;
		  name += "hBpt";
		  name += i;
		  name += "_result";
		  name += j;
		  hists_pt[i][j][f] = static_cast<TH1F*>(file.at(f)->Get(name));
		  if( hists_pt[i][j][f] == 0 )
		    {
		      histosExist = false;
		      std::cerr << "ERROR: " << name << " does not exist!" << std::endl;
		    }
		}
	    }
	}
    } // End of loop over files


  if( histosExist ) // Do comparison plots, if all histos exist
    {
      // Compare plots from different files
      TCanvas * const c1 = new TCanvas("c1","",600,600);
      TPostScript * const ps = new TPostScript("controlplotsDiJet_Comp.ps",111);
      TLegend * const leg = new TLegend(0.65,(0.9-0.05*nFILES),0.96,0.9);
      leg->SetFillColor(0);


      // Control quantities of B vs eta
      for(int i = 0; i < 8; i++) // Loop over B and Escale bins
	{
	  TString title = ",  after fit";
	  if( i%2 == 1 )  title = ",  before fit";

	  for(int j = 0; j < 4; j++) // Loop over some control quantities
	    {
	      leg->Clear();
	      for(int f = 0; f < nFILES; f++) // Loop over files
		{
		  c1->cd();
		  hists_eta[i][j][f]->SetStats(0);
		  hists_eta[i][j][f]->SetTitle(hists_eta[i][j][f]->GetTitle() + title);
		  if( f == 0 )
		    {
		      hists_eta[i][j][f]->Draw("P");
		    }
		  else
		    {
		      hists_eta[i][j][f]->SetMarkerStyle(hists_eta[i][j][f]->GetMarkerStyle()+4);
		      hists_eta[i][j][f]->SetMarkerColor(hists_eta[i][j][f]->GetMarkerColor()+1);
		      hists_eta[i][j][f]->SetLineColor(hists_eta[i][j][f]->GetLineColor()+1);
		      hists_eta[i][j][f]->Draw("P SAME");
		    }
		  leg->AddEntry(hists_eta[i][j][f],descr.at(f).c_str(),"P");
		}
	      leg->Draw("same");
	      c1->Draw();
	      ps->NewPage();
	    }
	}  // End of loop over Pt ratios and Escale bins


      // Control quantities of B vs Pt
      for(int i = 0; i < 2; i++) // Loop over B
	{
	  TString title = ",  after fit";
	  if( i%2 == 1 )  title = ",  before fit";

	  for(int j = 0; j < 4; j++) // Loop over some control quantities
	    {
	      leg->Clear();
	      for(int f = 0; f < nFILES; f++) // Loop over files
		{
		  c1->cd();
		  hists_pt[i][j][f]->SetStats(0);
		  hists_pt[i][j][f]->SetTitle(hists_pt[i][j][f]->GetTitle() + title);
		  if( f == 0 )
		    {
		      hists_pt[i][j][f]->Draw("P");
		    }
		  else
		    {
		      hists_pt[i][j][f]->SetMarkerStyle(hists_pt[i][j][f]->GetMarkerStyle()+4);
		      hists_pt[i][j][f]->SetMarkerColor(hists_pt[i][j][f]->GetMarkerColor()+1);
		      hists_pt[i][j][f]->SetLineColor(hists_pt[i][j][f]->GetLineColor()+1);
		      hists_pt[i][j][f]->Draw("P SAME");
		    }
		  leg->AddEntry(hists_pt[i][j][f],descr.at(f).c_str(),"P");
		}
	      leg->Draw("same");
	      c1->Draw();
	      ps->NewPage();
	    }
	}  // End of loop over B


      ps->Close();

      delete ps;
      delete leg;
      delete c1;
    }
}




//---------------------------------------------------------------
//   Compare Gamma-Jet Control Histograms
//---------------------------------------------------------------
void TControlPlotsComparison::CompareControlPlotsGammaJet(const std::vector<TFile*> &file, const std::vector<std::string> &descr) const
{
  int nFILES = static_cast<int>(file.size());

  // Declare histo pointers

  // Control quantities of Pt ratios vs eta, pt
  TH1F **hists_eta[12][8];
  TH1F **hists_pt[12][8];
  for(int i = 0; i < 12; i++)
    {
      for(int j = 0; j < 8; j++)
	{
	  hists_eta[i][j] = new TH1F*[nFILES];
	  hists_pt[i][j] = new TH1F*[nFILES];
	}
    }

  // Get histos from file
  bool histosExist = true;
  for(int f = 0; f < nFILES; f++) // Loop over files
    {
      TString dirName = file.at(f)->GetName();
      dirName += ":/GammaJet/";
      for(int i = 0; i < 12; i++)
	{
	  for(int j = 0; j < 8; j++)
	    {
	      TString name = dirName;
	      name += "heta";
	      name += i;
	      name += "_result";
	      name += j;
	      hists_eta[i][j][f] = static_cast<TH1F*>(file.at(f)->Get(name));
	      if( hists_eta[i][j][f] == 0 )
		{
		  histosExist = false;
		  std::cerr << "ERROR: " << name << " does not exist!" << std::endl;
		}

	      name = dirName;
	      name += "hpt_uncorr";
	      name += i;
	      name += "_result";
	      name += j;
	      hists_pt[i][j][f] = static_cast<TH1F*>(file.at(f)->Get(name));
	      if( hists_pt[i][j][f] == 0 )
		{
		  histosExist = false;
		  std::cerr << "ERROR: " << name << " does not exist!" << std::endl;
		}
	    }
	}
    } // End of loop over files


  if( histosExist ) // Do comparison plots, if all histos exist
    {
      // Compare plots from different files
      TCanvas * const c1 = new TCanvas("c1","",600,600);
      TPostScript * const ps = new TPostScript("controlplotsGammaJet_Comp.ps",111);
      TLegend * const leg = new TLegend(0.65,(0.9-0.05*nFILES),0.96,0.9);
      leg->SetFillColor(0);


      // Control quantities of Pt ratios vs eta
      for(int i = 0; i < 12; i++) // Loop over Pt ratios and Egamma bins
	{
	  if( i == 2  ||  i == 5  ||  i == 8  ||  i == 11 ) continue;

	  TString ptratio = ",  p^{jet}_{T}/E^{#gamma}_{T}";
	  if( i%3 == 1 )  ptratio = ",  p_{T}^{cor. jet}/E_{T}^{#gamma}";

	  for(int j = 0; j < 4; j++) // Loop over some control quantities
	    {
	      leg->Clear();
	      for(int f = 0; f < nFILES; f++) // Loop over files
		{
		  c1->cd();
		  hists_eta[i][j][f]->SetStats(0);
		  hists_eta[i][j][f]->SetTitle(hists_eta[i][j][f]->GetTitle() + ptratio);
		  if( f == 0 )
		    {
		      hists_eta[i][j][f]->Draw("P");
		    }
		  else
		    {
		      hists_eta[i][j][f]->SetMarkerStyle(hists_eta[i][j][f]->GetMarkerStyle()+4);
		      hists_eta[i][j][f]->SetMarkerColor(hists_eta[i][j][f]->GetMarkerColor()+1);
		      hists_eta[i][j][f]->SetLineColor(hists_eta[i][j][f]->GetLineColor()+1);
		      hists_eta[i][j][f]->Draw("P SAME");
		    }
		  leg->AddEntry(hists_eta[i][j][f],descr.at(f).c_str(),"P");
		}
	      leg->Draw("same");
	      c1->Draw();
	      ps->NewPage();
	    }
	}  // End of loop over Pt ratios and Egamma bins


      // Control quantities of Pt ratios vs jet Pt
      for(int i = 0; i < 12; i++) // Loop over Pt ratios and eta bins
	{
	  if( i == 2  ||  i == 5  ||  i == 8  ||  i == 11 ) continue;

	  TString ptratio = ",  p^{jet}_{T}/E^{#gamma}_{T}";
	  if( i%3 == 1 )  ptratio = ",  p_{T}^{cor. jet}/E_{T}^{#gamma}";

	  for(int j = 0; j < 4; j++) // Loop over some control quantities
	    {
	      leg->Clear();
	      for(int f = 0; f < nFILES; f++) // Loop over files
		{
		  c1->cd();
		  hists_pt[i][j][f]->SetStats(0);
		  hists_pt[i][j][f]->SetTitle(hists_pt[i][j][f]->GetTitle() + ptratio);
		  if( f == 0 )
		    {
		      hists_pt[i][j][f]->Draw("P");
		    }
		  else
		    {
		      hists_pt[i][j][f]->SetMarkerStyle(hists_pt[i][j][f]->GetMarkerStyle()+4);
		      hists_pt[i][j][f]->SetMarkerColor(hists_pt[i][j][f]->GetMarkerColor()+1);
		      hists_pt[i][j][f]->SetLineColor(hists_pt[i][j][f]->GetLineColor()+1);
		      hists_pt[i][j][f]->Draw("P SAME");
		    }
		  leg->AddEntry(hists_pt[i][j][f],descr.at(f).c_str(),"P");
		}
	      leg->Draw("same");
	      c1->Draw();
	      ps->NewPage();
	    }
	}  // End of loop over Pt ratios and Egamma bins


      ps->Close();

      delete ps;
      delete leg;
      delete c1;
    }
}




//---------------------------------------------------------------
// Set style option for ps output.
//---------------------------------------------------------------
void TControlPlotsComparison::SetGStyle() const
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
  gStyle->SetOptStat("neMR");
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
