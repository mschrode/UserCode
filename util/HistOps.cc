// $Id: $

#include "HistOps.h"

#include <iostream>


//!  \brief Write histograms a .root file
//!
//!  If a file with that name already exists, it is open
//!  and the objects are added to the file's content.
//!
//!  \param hists Pointers to the histograms to be written to file
//!  \param fileName Name of file to which the histograms are written
//!  \return Returns 0 if operation fails, some positive int else.
//---------------------------------------------------------------
int util::HistOps::WriteToRootFile(const std::vector<TH1F*>& hists, TString fileName)
{
  std::vector<TObject*> obj;
  for(std::vector<TH1F*>::const_iterator it = hists.begin(); it < hists.end(); it++)
    {
      obj.push_back(*it);
    }

  return WriteToRootFileGeneric(obj,fileName);
}



//---------------------------------------------------------------
void util::HistOps::DrawRatioPlotGeneric(TCanvas& can, const TH1F& h1, const TObject& h2, TString drawOption, bool drawRatioError)
{
  char name[100];
  sprintf(name,"hTop1_%s",h1.GetName());
  TH1F *hTop1 = static_cast<TH1F*>(h1.Clone(name));
  hTop1->GetXaxis()->SetLabelSize(0);
  hTop1->GetXaxis()->SetTitle("");

  sprintf(name,"hBottom_%s",h1.GetName());
  TH1F *hBottom = static_cast<TH1F*>(h1.Clone(name));
  hBottom->SetYTitle("");
  hBottom->GetYaxis()->SetNdivisions(505);
  hBottom->GetYaxis()->SetTickLength(0.08);

  sprintf(name,"hTop2_%s",h2.GetName());
  TObject *hTop2 = h2.Clone(name);

  TString className = h2.ClassName();
  if(  className.CompareTo("TH1F") == 0  )
    {
      TH1F *h = static_cast<TH1F*>(hTop2);
      h->GetXaxis()->SetLabelSize(0);
      h->GetXaxis()->SetTitle("");
      hBottom->Divide(h);
      if( drawRatioError )
	{
	  for(int bin = 1; bin <= hBottom->GetNbinsX(); bin++)
	    {
	      float sigma2 = hTop1->GetBinError(bin)*hTop1->GetBinError(bin) /
		( h->GetBinContent(bin)*h->GetBinContent(bin) )
		+ h->GetBinError(bin)*h->GetBinError(bin) *
		hTop1->GetBinContent(bin)*hTop1->GetBinContent(bin)  /
		( h->GetBinContent(bin)*h->GetBinContent(bin)*
		  h->GetBinContent(bin)*h->GetBinContent(bin) );
	      hBottom->SetBinError(bin,sqrt(sigma2)); 
	    }
	}
    }
  else if(  className.CompareTo("TF1") == 0  )
    {
      TF1 *f = static_cast<TF1*>(hTop2);
      f->GetXaxis()->SetLabelSize(0);
      f->GetXaxis()->SetTitle("");
      for(int bin = 1; bin <= hTop1->GetNbinsX(); bin++)
	{
	  hBottom->SetBinContent( bin, hBottom->GetBinContent(bin)/f->Eval(hBottom->GetBinCenter(bin)) );
	  if( drawRatioError )
	    {
	      float sigma2 = hTop1->GetBinError(bin)*hTop1->GetBinError(bin) /
		( f->Eval(hTop1->GetBinCenter(bin))*f->Eval(hTop1->GetBinCenter(bin)) );
	      hBottom->SetBinError(bin,sqrt(sigma2)); 
	    }
	}
    }

  can.cd();
  can.Clear();
  gPad->SetBottomMargin(0.4);
  hTop1->Draw(drawOption);
  if(  className.CompareTo("TH1F") == 0  ) hTop2->Draw(drawOption+"SAME");
  if(  className.CompareTo("TF1") == 0  ) hTop2->Draw("SAME");
  TPad *bottomPad = new TPad("bottomPad","",0,0,1,1);
  bottomPad->SetFillStyle(0);
  bottomPad->SetFrameFillColor(10);
  bottomPad->SetFrameBorderMode(0);
  bottomPad->SetTopMargin(0.6);
  bottomPad->SetGridy();
  bottomPad->Draw();
  bottomPad->cd();
  if( drawRatioError ) hBottom->Draw("PE");
  else hBottom->Draw("P");

  can.Update();
}



//!  \brief Fit several quantities from profiles of of a TH2F histogram.
//!
//!  Creates histograms with theses quantities (in brackets
//!  the basename of that histogram):
//!   - 0 - Mean value (Mean)
//!   - 1 - RMS (Standard deviation) / Mean value (Res)
//!   - 2 - Mean of Gauss fit (MeanGaus)
//!   - 3 - Sigma of Gauss fit / Mean of Gauss fit (ResGaus)
//!
//!  \param h2 2D histogram from which the profiles are generated
//!  \param namePrefix Prefix added to the histogram's basenames
//!  \param nameSuffix Suffix added to the histogram's basenames
//!  \return Vector of pointers to the created histograms
// ------------------------------------------------------
std::vector<TH1F*> util::HistOps::FitMean(const TH2F *h2, const std::string namePrefix, const std::string nameSuffix)
{
  // Set up 1D histos for mean values
  std::vector<TH1F*> hists;
  std::string name = namePrefix;
  name.append("Mean");
  name.append(nameSuffix);
  TH1F *h = new TH1F(name.c_str(),"",h2->GetNbinsX(),
		     h2->GetXaxis()->GetXmin(),
		     h2->GetXaxis()->GetXmax());
  hists.push_back(h);
  hists.back()->SetTitle("Mean");
  hists.back()->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());

  name = namePrefix;
  name.append("Res"); 
  name.append(nameSuffix);
  hists.push_back(static_cast<TH1F*>(h->Clone(name.c_str())));
  hists.back()->SetTitle("Normalized standard deviation");

  name = namePrefix;
  name.append("MeanGaus"); 
  name.append(nameSuffix);
  hists.push_back(static_cast<TH1F*>(h->Clone(name.c_str())));
  hists.back()->SetTitle("Mean of Gauss fit");

  name = namePrefix;
  name.append("ResGaus"); 
  name.append(nameSuffix);
  hists.push_back(static_cast<TH1F*>(h->Clone(name.c_str())));
  hists.back()->SetTitle("Normalized sigma of Gauss fit");

  // Loop over x-bins of h2, fit mean values
  // and fill them into hists
  TH1F *hProfile = new TH1F("hProfile","",h2->GetNbinsY(),
			    h2->GetYaxis()->GetXmin(),
			    h2->GetYaxis()->GetXmax());
  for(int binX = 1; binX <= h2->GetNbinsX(); binX++) // Loop over x-bins
    {
      hProfile->Reset();
      for(int binY = 1; binY <= h2->GetNbinsY(); binY++) // Loop over y-bins
	{
	  hProfile->SetBinContent( binY,h2->GetBinContent( h2->GetBin(binX,binY) ) );
	  hProfile->SetBinError(   binY,h2->GetBinError(   h2->GetBin(binX,binY) ) );
	}

      // Mean value
      hists.at(0)->SetBinContent( binX, hProfile->GetMean()      );
      hists.at(0)->SetBinError(   binX, hProfile->GetMeanError() );

      // RMS / Mean value
      if(hProfile->GetMean())
	{
	  hists.at(1)->SetBinContent( binX, (hProfile->GetRMS())/(hProfile->GetMean()) );
	  double err2 = pow(hProfile->GetRMSError(),2) / pow(hProfile->GetMean(),2);
	  err2       += pow(hProfile->GetRMS(),4)      / pow(hProfile->GetMean(),4);
	  hists.at(1)->SetBinError( binX, sqrt(err2) );

// 	  hists.at(1)->SetBinContent( binX, hProfile->GetRMS() );
// 	  hists.at(1)->SetBinError(   binX, hProfile->GetRMSError() );
	}
      
      // Mean and sigma of gausfit
      // 1) Fit gaussian over whole range
      if( hProfile->GetEntries() )
	{
	  int fitOk = hProfile->Fit("gaus","IL0Q");
	  if( fitOk == 0 )
	    {
	      TF1 *fit = hProfile->GetFunction("gaus");
	      
	      // 2) Fit gaussian again in range +/-sigma
	      // around mean of first fit
	      fitOk = hProfile->Fit("gaus","IL0Q","",
				    fit->GetParameter(1) - 1.5*fit->GetParameter(2),
				    fit->GetParameter(1) + 1.5*fit->GetParameter(2) );
	      if( fitOk == 0 )
		{
		  fit = hProfile->GetFunction("gaus");
		  hists.at(2)->SetBinContent( binX, fit->GetParameter(1) );
		  hists.at(2)->SetBinError(   binX, fit->GetParError(1)  );
		  if( fit->GetParameter(1) )
		    {
		      hists.at(3)->SetBinContent( binX, (fit->GetParameter(2))/(fit->GetParameter(1)) );
		      double err2 = pow(fit->GetParError(2),2)  / pow(fit->GetParameter(1),2);
		      err2       += pow(fit->GetParameter(2),4) / pow(fit->GetParameter(1),4);
		      hists.at(3)->SetBinError( binX, sqrt(err2) );

// 		      hists.at(3)->SetBinContent( binX, fit->GetParameter(2) );
// 		      hists.at(3)->SetBinError(   binX, fit->GetParError(2) );
		    }
		}
	    }
	}
      else
	{
	  hists.at(2)->SetBinContent( binX, 0 );
	  hists.at(3)->SetBinContent( binX, 0 );
	}
    }
  // Clean up
  delete hProfile;

  return hists;
}



TH1F* util::HistOps::GetRelDiff(const TH1F *hBase, const TH1F* hist)
{
  TString name = hist->GetName();
  name        += "Diff";
  TH1F *hDiff  = static_cast<TH1F*>(hist->Clone(name));
  hDiff->Reset();
  hDiff->SetTitle("Relative difference");
  hDiff->SetYTitle("");
  hDiff->SetMarkerStyle(hist->GetMarkerStyle());
  hDiff->SetMarkerColor(hist->GetMarkerColor());
  hDiff->SetLineColor(hist->GetLineColor());
  for(int bin = 1; bin <= hDiff->GetNbinsX(); bin++)
    {
      float x1     = hBase->GetBinContent(bin);
      float x2     = hist->GetBinContent(bin);
      float sigma1 = hBase->GetBinError(bin);
      float sigma2 = hist->GetBinError(bin);
      if( x1 )
	{
	  hDiff->SetBinContent(bin,(x2 - x1)/x1);
	  hDiff->SetBinError(bin,sqrt(pow(sigma2/x1,2) + pow(x2*sigma1/x1/x1,2)));
	}
    }

  return hDiff;
}



//!  \brief Get histogram with relative difference of two or more other histograms
//!
//!  \note Histograms must have the same binning.
//!
//!  \param hBase Histogram to which ratio is calculated
//!  \param hists Histograms from which ratio is calculated
//!  \return Histograms with ratios
// ------------------------------------------------------
std::vector<TH1F*> util::HistOps::GetRelDiff(const TH1F *hBase, const std::vector<TH1F*> hists)
{
  std::vector<TH1F*> result;
  for(std::vector<TH1F*>::const_iterator it = hists.begin(); it != hists.end(); it++)
    {
      result.push_back(util::HistOps::GetRelDiff(hBase,*it));
    }

  return result;
}



// ------------------------------------------------------
int util::HistOps::WriteToRootFileGeneric(const std::vector<TObject*>& obj, TString fileName)
{
  int ok = 1;
  TFile outFile(fileName,"UPDATE");
  for(std::vector<TObject*>::const_iterator it = obj.begin(); it < obj.end(); it++)
    {
      ok *= outFile.WriteTObject( *it );
      if( !ok ) std::cerr << "Error writing '" << (*it)->GetName() << "' to file " << fileName << "." << std::endl;
    }
  outFile.Close();

  return ok;
}
