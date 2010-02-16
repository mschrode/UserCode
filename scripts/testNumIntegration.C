#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"



double func(const double *x, const double *par) {
  //  return par[0]*pow(x[0]+par[1],2) + par[2];
  //return exp(-par[0]*x[0]);
  return par[0]*sin(par[1]*(x[0]-par[2]));
}

double func(double x, const std::vector<double> &par) {
  return func(&x,&(par.front()));
}


double integral(double a, double b, const std::vector<double> &par) {
  //  return par[0]*(pow(par[1]+b,3)-pow(par[1]+a,3))/3. + par[2]*(b-a);
  //return (exp(-par[0]*a)-exp(-par[0]*b))/par[0];
  return (cos(par[1]*(a-par[2]))-cos(par[1]*(b-par[2])))*par[0]/par[1];
}


double numIntegralStripes(double a, double b, const std::vector<double> &par, int n) {
  double sum = 0.;
  double h = (b-a)/n;
  for(int i = 0; i < n; i++) {
    double x = a + (i+0.5)*h;
    sum += func(x,par);
  }
  sum *= h;
  
  return sum;
}


double numIntegralTrapezium(double a, double b, const std::vector<double> &par, int n) {
  double sum = 0.5*(func(a,par)+func(b,par));
  double h = (b-a)/n;
  for(int i = 1; i < n; i++) {
    sum += func(a+i*h,par);
  }
  sum *= h;
  
  return sum;
}


double numIntegralSimpsons(double a, double b, const std::vector<double> &par, int n) {
  double h = b - a;
  double iVal = 0.;
  double iValOld = 1.;
  int nIter = 0;
  int maxNIter = log(n)/log(3.) - 1.;
  if( maxNIter < 1 ) maxNIter = 1;
  std::vector<double> vals;    
  std::vector<double> valsOld;

  while( nIter < maxNIter ) {
    iValOld = iVal;
    iVal = 0.;
    valsOld = vals;
    vals.clear();
    h /= 3.;    // In each iteration, split h into 3 new intervals
    
    // Loop over nodes xi i.e. interval borders
    for(int i = 0; i <= pow(3.0,nIter+1); ++i){
      // Calculate value only at new nodes
      if(nIter == 0 || i % 3 != 0) {
	vals.push_back(func(a+i*h,par));
      } else {
	vals.push_back(valsOld.at(i/3));
      }
    }
    
    // Sum up weighted function values
    for(size_t i = 0; i < vals.size(); i++)	{
      double w = 1.;                       // Weight w from Simpson's rule
      if( i > 0 && i < (vals.size() - 1) ) { // w = 1 for x0 and last node
	if( i % 3 == 0 ) {                   // w = 2 for x3, x6, ...
	  w = 2.;
	} else {
	  w = 3.;
	}
      }
      iVal += w * (vals[i]);  // Sum up weighted function values
    }
    iVal *= (3. * h / 8.);    // Apply overall normalization
    nIter++;
  }

  return iVal;
}


void compareMethods(double a, double b, const std::vector<double> &par) {
  TF1 *f = new TF1("f",func,a,b,par.size());
  for(size_t i = 0; i < par.size(); i++) {
    f->SetParameter(i,par[i]);
  }
  TCanvas *can1 = new TCanvas("can1","Function",500,500);
  can1->cd();
  f->Draw();

  std::cout << "\nINTEGRAL " << a << " to " << b << std::endl;
  std::cout << "Method        \tValue"<< std::endl;
  std::cout << "Analytic      \t" << integral(a,b,par) << std::endl;
  std::cout << "ROOT          \t" << f->Integral(a,b) << std::endl;
  std::cout << "Stripes (5)   \t" << numIntegralStripes(a,b,par,5) << std::endl;
  std::cout << "Trapezium (5) \t" << numIntegralTrapezium(a,b,par,5) << std::endl;
  std::cout << "Simpson's (5) \t" << numIntegralSimpsons(a,b,par,5) << std::endl;
  std::cout << "Stripes (50)  \t" << numIntegralStripes(a,b,par,50) << std::endl;
  std::cout << "Trapezium (50)\t" << numIntegralTrapezium(a,b,par,50) << std::endl;
  std::cout << "Simpson's (50)\t" << numIntegralSimpsons(a,b,par,50) << std::endl;

  std::vector<double> nVal;
  nVal.push_back(50.);
  nVal.push_back(100.);
  nVal.push_back(200.);
  nVal.push_back(500.);
  nVal.push_back(1000.);
  double trueInt = integral(a,b,par);
  std::vector<double> dIntStripes(nVal.size());
  std::vector<double> dIntTrapezium(nVal.size());
  std::vector<double> dIntSimpsons(nVal.size());
  for(size_t i = 0; i < nVal.size(); i++) {
    int n = static_cast<int>(nVal[i]);
    dIntStripes[i] = (numIntegralStripes(a,b,par,n) - trueInt)/trueInt; 
    dIntTrapezium[i] = (numIntegralTrapezium(a,b,par,n) - trueInt)/trueInt;
    dIntSimpsons[i] = (numIntegralSimpsons(a,b,par,n) - trueInt)/trueInt;
  }

  TGraph *gStripes = new TGraph(dIntStripes.size(),&(nVal.front()),&(dIntStripes.front()));
  gStripes->SetLineColor(2);
  gStripes->SetMarkerColor(2);
  TGraph *gTrapezium = new TGraph(dIntTrapezium.size(),&(nVal.front()),&(dIntTrapezium.front()));
  gTrapezium->SetLineColor(4);
  gTrapezium->SetMarkerColor(4);
  TGraph *gSimpsons = new TGraph(dIntSimpsons.size(),&(nVal.front()),&(dIntSimpsons.front()));
 
  TCanvas *can2 = new TCanvas("can2","Relative difference",500,500);
  can2->cd();
  TH1F *hFrame = new TH1F("hFrame","Number of intervals",1000,0.5*nVal.front(),1.05*nVal.back());
  hFrame->GetYaxis()->SetRangeUser(0.9*dIntStripes.front(),1.1*dIntTrapezium.front());
  hFrame->Draw();
  gStripes->Draw("L");
  gTrapezium->Draw("Lsame");
  gSimpsons->Draw("Lsame");
}


void testNumIntegration() {
  std::vector<double> par(3);
  par[0] = 1.2;
  par[1] = 1.1;
  par[2] = -2.;
  compareMethods(5.3,17.2,par);
}
