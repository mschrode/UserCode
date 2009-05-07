// $Id: run.cc,v 1.4 2009/05/05 13:58:37 mschrode Exp $

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"

#include "ControlPlots.h"
#include "Event.h"
#include "EventGenerator.h"
#include "Fitter.h"


void WriteParameters(const std::vector<double>& par, const std::string& head)
{
  std::string name = "parameters.txt";
  ofstream outfile(name.c_str(),std::ios::app);
  if(outfile.is_open())
    {
      outfile << head << " : ";
      for(unsigned int i = 0; i < par.size(); i++)
	{
	  outfile << par.at(i) << "  ";
	}
      outfile << "\n";
      outfile.close();
    }
}


int main()
{
  std::string model = "ThreeGauss";

  // Parameter vector
  std::vector<double> par;

  // Set start parameters
  if( model == "FermiTail" )
    {
      par.push_back(0.9);
      par.push_back(0.1);
      par.push_back(0.05);
    }
  else if( model == "TwoGauss" )
    {
      par.push_back(0.24);
      par.push_back(0.07);
      par.push_back(2.06);
      par.push_back(0.90);
    }
  else if( model == "ThreeGauss" )
    {
      par.push_back(0.2);  // Sigma of main Gaussian
      par.push_back(0.002);  // Norm of 1. Gaussian tail
      par.push_back(0.6);   // Mean of 1. Gaussian tail
      par.push_back(0.2);   // Sigma of 1. Gaussian tail
      par.push_back(0.05);  // Norm of 2. Gaussian tail
      par.push_back(2.3);   // Mean of 2. Gaussian tail
      par.push_back(1.0);   // Sigma of 2. Gaussian tail
    }
  else if( model == "ExpTail" )
    {
      par.push_back(0.95);
      par.push_back(0.2);
      par.push_back(1.2);
      par.push_back(1.5);
    }


  std::vector<std::string> ptbins;
  //ptbins.push_back("15_20");
  //ptbins.push_back("20_30");
  ptbins.push_back("30_50");
   ptbins.push_back("50_80");
   ptbins.push_back("80_120");
  ptbins.push_back("120_170");
  ptbins.push_back("170_230");
  ptbins.push_back("230_300");
  ptbins.push_back("300_380");
  ptbins.push_back("380_470");
  ptbins.push_back("470_600");
  ptbins.push_back("600_800");

  for(unsigned int i = 0; i < ptbins.size(); i++)
    {
      std::cout << std::endl << "PT BIN " << ptbins.at(i) << std::endl;

      // Generate events
      double min = 100.;
      double max = 1000;
      std::vector<double> parTruth;
      std::vector<double> parResp;
      parResp.push_back(0.05);

      TH1F * h = 0;
      TFile * f = new TFile("input/RespDiJetHisto.root","READ");
      std::string name = "hRespCorr_"+ptbins.at(i);
      f->GetObject(name.c_str(),h);
      if( h == 0 )
	{
	  f->Close();
	  delete f;
	  continue;
	}
      h->SetDirectory(0);
      f->Close();
      delete f;

      js::EventGenerator gen(min,max,"Uniform",parTruth,"Histogram",parResp);
      gen.SetRespHist(h);

      js::Data data;

      js::Data dataDiJet = gen.GenerateDijetEvents(1000);
      data.insert(data.end(),dataDiJet.begin(),dataDiJet.end());

      js::Data dataPhotonJet = gen.GeneratePhotonJetEvents(1000);
      data.insert(data.end(),dataPhotonJet.begin(),dataPhotonJet.end());

      // Fit
      js::Fitter fitter(data,min,max,model,par);
      for(int iter = 0; iter < 2; iter++)
	{
	  std::cout << ">>>>>>>> Iteration " << iter << std::endl;
	  fitter.Fit();
	  
	  // Copy parameters
	  for(int p = 0; p < fitter.GetNPar(); p++)
	    {
	      par.at(p) = fitter.GetPar(p);
	    }
	}

      // Write parameters to file
      WriteParameters(par,"PT BIN "+ptbins.at(i));

      // Make control plots
      if( true )
	{
	  js::ControlPlots plots(data);
	  plots.SetFileNameSuffix(model+"_"+ptbins.at(i));
	  plots.SetRespBinning(40,0,6);
	  //	  plots.PlotDijets();
	  //	  plots.PlotPhotonJets();
	  plots.PlotResponse(fitter.GetTF1());
	}

      // Clean up
      std::vector<js::Event*>::iterator it = data.begin();
      for(; it != data.end(); it++)
	{
	  delete *it;
	}
      data.clear();
    }

  return 0;
}
