#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"

#include "ControlPlots.h"
#include "EventGenerator.h"
#include "Fitter.h"


int main()
{
  // Generate events
  double min = 100.;
  double max = 1000;
  std::vector<double> parTruth;
  std::vector<double> parResp;
//   parResp.push_back(1.2);
//   parResp.push_back(0.1);
  parResp.push_back(0.05);
  TH1F * h = 0;
  TFile * f = new TFile("input/RespHisto.root","READ");
  f->GetObject("hRespCorr",h);
  if( h == 0 ) return 1;
  js::EventGenerator gen(min,max,"Uniform",parTruth,"Histogram",parResp);
  gen.SetRespHist(h);
  js::Data data = gen.GenerateDijetEvents(1000);

  // Fit
  std::vector<double> par;
//   par.push_back(0.9);
//   par.push_back(0.1);
//   par.push_back(0.05);

  par.push_back(0.21);
  par.push_back(0.04);
  par.push_back(1.5);
  par.push_back(0.4);
  js::Fitter fitter(data,min,max,"TwoGauss",par);
  fitter.Fit();

  // Make control plots
  js::ControlPlots plots(data);
  plots.PlotDijets();
  plots.PlotResponse(fitter.GetTF1());

  return 0;
}
