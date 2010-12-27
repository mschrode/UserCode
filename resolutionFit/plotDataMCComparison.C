#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../util/HistOps.h"
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


const bool pf_ = true;
const unsigned int etaBin = 0;
const unsigned int interpolationMode = 1;



void fill(std::vector<double> &ptData, std::vector<double> &ptDataErr,
	  std::vector<double> &resData, std::vector<double> &resDataErr,
	  std::vector<double> &ptMC, std::vector<double> &ptMCErr,
	  std::vector<double> &resMC, std::vector<double> &resMCErr) {


  if( pf_ ) {
    if( etaBin == 0 ) {
      ptMC.push_back(58.497);
      ptMCErr.push_back(0.101635);
      resMC.push_back(0.118233);
      resMCErr.push_back(0.000145685);
      ptMC.push_back(73.4567);
      ptMCErr.push_back(0.0918879);
      resMC.push_back(0.114351);
      resMCErr.push_back(0.000210336);
      ptMC.push_back(88.4721);
      ptMCErr.push_back(0.0936784);
      resMC.push_back(0.107541);
      resMCErr.push_back(0.000310972);
      ptMC.push_back(103.442);
      ptMCErr.push_back(0.0904839);
      resMC.push_back(0.100316);
      resMCErr.push_back(0.000426455);
      ptMC.push_back(118.551);
      ptMCErr.push_back(0.0904728);
      resMC.push_back(0.0949578);
      resMCErr.push_back(0.000545242);
      ptMC.push_back(135.698);
      ptMCErr.push_back(0.100805);
      resMC.push_back(0.0866737);
      resMCErr.push_back(0.000595124);
      ptMC.push_back(155.579);
      ptMCErr.push_back(0.0956096);
      resMC.push_back(0.0821818);
      resMCErr.push_back(0.000743109);
      ptMC.push_back(175.479);
      ptMCErr.push_back(0.098313);
      resMC.push_back(0.0818675);
      resMCErr.push_back(0.00100935);
      ptMC.push_back(199.495);
      ptMCErr.push_back(0.121463);
      resMC.push_back(0.0752986);
      resMCErr.push_back(0.00101327);
      ptMC.push_back(229.391);
      ptMCErr.push_back(0.119654);
      resMC.push_back(0.0731143);
      resMCErr.push_back(0.0013306);
      ptMC.push_back(266.655);
      ptMCErr.push_back(0.153905);
      resMC.push_back(0.0691148);
      resMCErr.push_back(0.00136133);
      ptMC.push_back(317.139);
      ptMCErr.push_back(0.156735);
      resMC.push_back(0.0637833);
      resMCErr.push_back(0.00188231);
      ptMC.push_back(367.087);
      ptMCErr.push_back(0.155769);
      resMC.push_back(0.0603893);
      resMCErr.push_back(0.00256069);
      ptMC.push_back(433.519);
      ptMCErr.push_back(0.21931);
      resMC.push_back(0.0591722);
      resMCErr.push_back(0.00260624);
      ptMC.push_back(584.954);
      ptMCErr.push_back(0.451223);
      resMC.push_back(0.0543258);
      resMCErr.push_back(0.00276959);

      ptData.push_back(57.1173);
      ptDataErr.push_back(0.65706);
      resData.push_back(0.152034);
      resDataErr.push_back(0.0145124);
      ptData.push_back(73.217);
      ptDataErr.push_back(0.281124);
      resData.push_back(0.120034);
      resDataErr.push_back(0.00530784);
      ptData.push_back(88.49);
      ptDataErr.push_back(0.409836);
      resData.push_back(0.103954);
      resDataErr.push_back(0.00691724);
      ptData.push_back(103.695);
      ptDataErr.push_back(0.157742);
      resData.push_back(0.0983386);
      resDataErr.push_back(0.00257517);
      ptData.push_back(118.324);
      ptDataErr.push_back(0.208103);
      resData.push_back(0.0957513);
      resDataErr.push_back(0.00329241);
      ptData.push_back(135.089);
      ptDataErr.push_back(0.201415);
      resData.push_back(0.094938);
      resDataErr.push_back(0.00255272);
      ptData.push_back(154.868);
      ptDataErr.push_back(0.270037);
      resData.push_back(0.0929255);
      resDataErr.push_back(0.00344508);
      ptData.push_back(175.061);
      ptDataErr.push_back(0.265093);
      resData.push_back(0.0878645);
      resDataErr.push_back(0.00308302);
      ptData.push_back(199.136);
      ptDataErr.push_back(0.428186);
      resData.push_back(0.0808423);
      resDataErr.push_back(0.00307998);
      ptData.push_back(229.425);
      ptDataErr.push_back(0.398127);
      resData.push_back(0.0715905);
      resDataErr.push_back(0.00263551);
      ptData.push_back(266.054);
      ptDataErr.push_back(0.699017);
      resData.push_back(0.0766838);
      resDataErr.push_back(0.00287992);
      ptData.push_back(317.718);
      ptDataErr.push_back(0.96152);
      resData.push_back(0.0583359);
      resDataErr.push_back(0.00328466);
      ptData.push_back(365.809);
      ptDataErr.push_back(1.20588);
      resData.push_back(0.0679423);
      resDataErr.push_back(0.00678646);
      ptData.push_back(431.464);
      ptDataErr.push_back(2.77005);
      resData.push_back(0.0731373);
      resDataErr.push_back(0.0060503);
      ptData.push_back(584.347);
      ptDataErr.push_back(9.42048);
      resData.push_back(0.0556621);
      resDataErr.push_back(0.00460951);

    } else if( etaBin == 1 ) {
      ptMC.push_back(58.1691);
ptMCErr.push_back(0.260036);
resMC.push_back(0.135855);
resMCErr.push_back(0.000414643);
ptMC.push_back(73.3287);
ptMCErr.push_back(0.242896);
resMC.push_back(0.122355);
resMCErr.push_back(0.000561018);
ptMC.push_back(88.8965);
ptMCErr.push_back(0.232443);
resMC.push_back(0.0968478);
resMCErr.push_back(0.000757249);
ptMC.push_back(103.87);
ptMCErr.push_back(0.237017);
resMC.push_back(0.0969345);
resMCErr.push_back(0.00107816);
ptMC.push_back(118.703);
ptMCErr.push_back(0.229104);
resMC.push_back(0.0995849);
resMCErr.push_back(0.00153906);
ptMC.push_back(142.094);
ptMCErr.push_back(0.367687);
resMC.push_back(0.0949359);
resMCErr.push_back(0.00118197);
ptMC.push_back(186.135);
ptMCErr.push_back(0.429616);
resMC.push_back(0.0834621);
resMCErr.push_back(0.00191719);
ptMC.push_back(244.085);
ptMCErr.push_back(0.540785);
resMC.push_back(0.0787844);
resMCErr.push_back(0.00279943);
ptMC.push_back(345.223);
ptMCErr.push_back(0.863075);
resMC.push_back(0.0662931);
resMCErr.push_back(0.00391661);

      ptData.push_back(53.3834);
ptDataErr.push_back(0.91443);
resData.push_back(0.177974);
resDataErr.push_back(0.028147);
ptData.push_back(70.7698);
ptDataErr.push_back(0.635486);
resData.push_back(0.130876);
resDataErr.push_back(0.01449);
ptData.push_back(89.6339);
ptDataErr.push_back(0.897656);
resData.push_back(0.0723452);
resDataErr.push_back(0.0121724);
ptData.push_back(103.217);
ptDataErr.push_back(0.357416);
resData.push_back(0.119546);
resDataErr.push_back(0.00853183);
ptData.push_back(116.28);
ptDataErr.push_back(0.500219);
resData.push_back(0.120272);
resDataErr.push_back(0.0108239);
ptData.push_back(141.087);
ptDataErr.push_back(0.795346);
resData.push_back(0.110689);
resDataErr.push_back(0.00531411);
ptData.push_back(184.837);
ptDataErr.push_back(1.02443);
resData.push_back(0.107283);
resDataErr.push_back(0.00607363);
ptData.push_back(245.144);
ptDataErr.push_back(1.76513);
resData.push_back(0.063704);
resDataErr.push_back(0.00623297);
ptData.push_back(344.725);
ptDataErr.push_back(6.05185);
resData.push_back(0.0666585);
resDataErr.push_back(0.00753492);

    } else if( etaBin == 2 ) {
      
ptMC.push_back(58.4463);
ptMCErr.push_back(0.300917);
resMC.push_back(0.114213);
resMCErr.push_back(0.000453183);
ptMC.push_back(74.4145);
ptMCErr.push_back(0.276893);
resMC.push_back(0.103272);
resMCErr.push_back(0.000630599);
ptMC.push_back(89.3857);
ptMCErr.push_back(0.275997);
resMC.push_back(0.0918158);
resMCErr.push_back(0.000891293);
ptMC.push_back(104.563);
ptMCErr.push_back(0.30251);
resMC.push_back(0.0769392);
resMCErr.push_back(0.00128737);
ptMC.push_back(119.355);
ptMCErr.push_back(0.287478);
resMC.push_back(0.0829106);
resMCErr.push_back(0.00172885);
ptMC.push_back(143.184);
ptMCErr.push_back(0.494183);
resMC.push_back(0.0670889);
resMCErr.push_back(0.00128691);
ptMC.push_back(186.372);
ptMCErr.push_back(0.621386);
resMC.push_back(0.0588537);
resMCErr.push_back(0.00217789);
ptMC.push_back(249.33);
ptMCErr.push_back(0.843428);
resMC.push_back(0.0603056);
resMCErr.push_back(0.00227022);

ptData.push_back(57.2962);
ptDataErr.push_back(1.36532);
resData.push_back(0.164893);
resDataErr.push_back(0.170505);
ptData.push_back(70.4808);
ptDataErr.push_back(0.96599);
resData.push_back(0.119739);
resDataErr.push_back(0.0174033);
ptData.push_back(90.3149);
ptDataErr.push_back(2.08433);
resData.push_back(nan);
resDataErr.push_back(0.0275772);
ptData.push_back(105.302);
ptDataErr.push_back(0.60548);
resData.push_back(0.0589376);
resDataErr.push_back(0.00628139);
ptData.push_back(117.679);
ptDataErr.push_back(0.770668);
resData.push_back(0.116786);
resDataErr.push_back(0.0165647);
ptData.push_back(143.089);
ptDataErr.push_back(1.0602);
resData.push_back(0.0644072);
resDataErr.push_back(0.00528558);
ptData.push_back(185.283);
ptDataErr.push_back(1.9398);
resData.push_back(0.0708145);
resDataErr.push_back(0.00924586);
ptData.push_back(247.366);
ptDataErr.push_back(3.04332);
resData.push_back(0.0751586);
resDataErr.push_back(0.00532112);

      
    } else if( etaBin == 3 ) {

ptMC.push_back(56.8718);
ptMCErr.push_back(0.205635);
resMC.push_back(0.120576);
resMCErr.push_back(0.000316705);
ptMC.push_back(72.6559);
ptMCErr.push_back(0.235905);
resMC.push_back(0.0951205);
resMCErr.push_back(0.000459364);
ptMC.push_back(86.2279);
ptMCErr.push_back(0.271733);
resMC.push_back(0.105181);
resMCErr.push_back(0.000895434);
ptMC.push_back(103.675);
ptMCErr.push_back(0.403937);
resMC.push_back(0.0647532);
resMCErr.push_back(0.00104957);
ptMC.push_back(117.708);
ptMCErr.push_back(0.471907);
resMC.push_back(0.0776072);
resMCErr.push_back(0.00212727);
ptMC.push_back(136.283);
ptMCErr.push_back(0.852915);
resMC.push_back(0.0811158);
resMCErr.push_back(0.00244864);
ptMC.push_back(177.901);
ptMCErr.push_back(1.82514);
resMC.push_back(0.0405981);
resMCErr.push_back(0.0048831);
ptMC.push_back(229.653);
ptMCErr.push_back(3.34474);
resMC.push_back(0.0744346);
resMCErr.push_back(0.0138573);

ptData.push_back(57.4192);
ptDataErr.push_back(1.22295);
resData.push_back(0.120944);
resDataErr.push_back(0.0201999);
ptData.push_back(70.0778);
ptDataErr.push_back(0.718827);
resData.push_back(0.0985142);
resDataErr.push_back(0.0128748);
ptData.push_back(87.6139);
ptDataErr.push_back(2.39541);
resData.push_back(0);
resDataErr.push_back(100000);
ptData.push_back(97.7769);
ptDataErr.push_back(0.633409);
resData.push_back(0.146733);
resDataErr.push_back(0.0190621);
ptData.push_back(116.993);
ptDataErr.push_back(1.45501);
resData.push_back(0.063973);
resDataErr.push_back(0.0153374);
ptData.push_back(138.801);
ptDataErr.push_back(1.19177);
resData.push_back(0.0493269);
resDataErr.push_back(0.00796352);
ptData.push_back(183.095);
ptDataErr.push_back(2.93815);
resData.push_back(0.0490037);
resDataErr.push_back(0.011322);
ptData.push_back(224.652);
ptDataErr.push_back(0);
resData.push_back(0.0971119);
resDataErr.push_back(0.0463379);

    }
  }
}



TGraphAsymmErrors* ratioGammaJet() {
  std::vector<double> pt;
  std::vector<double> ptErr;
  std::vector<double> ratio;
  std::vector<double> ratioErr;

  if( etaBin == 0 ) {
    pt.push_back(60.);
    ptErr.push_back(0.);
    ratio.push_back(1.102);
    ratioErr.push_back(0.0523501);
  } else if( etaBin == 1 ) {
    
    
  } else if( etaBin == 2 ) {
    
    
  } else if( etaBin == 3 ) {
    
    
  }

  return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(ratio.front()),&(ptErr.front()),&(ptErr.front()),&(ratioErr.front()),&(ratioErr.front()));			      
}



TF1* fitResolution(TGraphAsymmErrors* g, const TString &name, double min, double max) {
  TF1* fit = 0;
  if( interpolationMode == 0 ) {
    fit = new TF1(name,"sqrt(sq([0]/x) + sq([1])/x + sq([2]))",min,max);
    fit->FixParameter(0,0.);
    fit->SetParameter(1,1.);
    fit->SetParameter(2,0.01);
  } else if( interpolationMode == 1 ) {
    fit = new TF1(name,"sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",min,max);
    fit->SetParameter(0,-0.34);
    fit->SetParameter(1,0.3);
    fit->FixParameter(2,0.);
    fit->SetParameter(3,0.47);
  }
  fit->SetLineWidth(1);
  g->Fit(fit,"BNR");

  return fit;
}


TF1* fitAddConst(TGraphAsymmErrors* g, const TString &name, const TF1* mcFit) {
  TF1* fit = 0;
  if( interpolationMode == 0 ) {
    fit = new TF1(name,"sqrt(sq([0]/x) + sq([1])/x + sq([2]) + sq([3]))",mcFit->GetXmin(),mcFit->GetXmax());
    fit->FixParameter(0,mcFit->GetParameter(0));
    fit->FixParameter(1,mcFit->GetParameter(1));
    fit->FixParameter(2,mcFit->GetParameter(2));
    fit->SetParameter(3,mcFit->GetParameter(3));
  } else if( interpolationMode == 1 ) {
    fit = new TF1(name,"sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",mcFit->GetXmin(),mcFit->GetXmax());
    fit->FixParameter(0,mcFit->GetParameter(0));
    fit->FixParameter(1,mcFit->GetParameter(1));
    fit->SetParameter(2,0.);
    fit->FixParameter(3,mcFit->GetParameter(3));
  }
  fit->SetLineWidth(1);
  g->Fit(fit,"BNR");

  return fit;
}


void plotDataMCComparison() {

  util::StyleSettings::presentationNoTitle();

  // Input
  double ptMin = 55.;
  double ptMax = 1000.;

  TString outNamePrefix = "ResFit_DataMCComparison_";
  if( pf_ ) outNamePrefix += "PF_";
  else outNamePrefix += "Calo_";
  outNamePrefix += "Eta0_";

  std::vector<double> ptData;
  std::vector<double> ptDataErr;
  std::vector<double> resData;
  std::vector<double> resDataErr;
  std::vector<double> ptMC;
  std::vector<double> ptMCErr;
  std::vector<double> resMC;
  std::vector<double> resMCErr;

  fill(ptData,ptDataErr,resData,resDataErr,ptMC,ptMCErr,resMC,resMCErr);



  // Resolution
  TGraphAsymmErrors* gData = new TGraphAsymmErrors(ptData.size(),&(ptData.front()),&(resData.front()),
						   &(ptDataErr.front()),&(ptDataErr.front()),
						   &(resDataErr.front()),&(resDataErr.front()));
  gData->SetMarkerStyle(25);
  gData->SetMarkerColor(2);
  gData->SetLineColor(2);

  TGraphAsymmErrors* gMC = new TGraphAsymmErrors(ptMC.size(),&(ptMC.front()),&(resMC.front()),
						 &(ptMCErr.front()),&(ptMCErr.front()),
						 &(resMCErr.front()),&(resMCErr.front()));
  gMC->SetMarkerStyle(21);

  // Fit resolution
  TF1* fitMC = fitResolution(gMC,"fitData",ptMin,ptMax);
  fitMC->SetLineStyle(2);
  fitMC->SetLineColor(kBlack);

  TF1* fitData1 = fitResolution(gData,"fitData1",ptMin,ptMax);
  fitData1->SetLineColor(kRed);
  TF1* fitData2 = fitAddConst(gData,"fitData2",fitMC);
  fitData2->SetLineColor(kBlue);


  // Ratios
  TGraphAsymmErrors* gRatio = util::HistOps::createRatioGraph(gData,gMC);
  gRatio->SetMarkerStyle(21);
  gRatio->SetMarkerColor(1);
  gRatio->SetLineColor(1);

  TGraphAsymmErrors* gRatioGJ = ratioGammaJet();
  gRatioGJ->SetMarkerStyle(20);
  gRatioGJ->SetMarkerColor(kRed);
  gRatioGJ->SetLineColor(gRatioGJ->GetMarkerColor());

  TF1* fitRatio1 = 0;
  if( interpolationMode == 0 ) {
    fitRatio1 = new TF1("fitRatio1","sqrt(sq([0]/x) + sq([1])/x + sq([2])) / sqrt(sq([3]/x) + sq([4])/x + sq([5]))",ptMin,ptMax);
  } else if( interpolationMode == 1 ) {
    fitRatio1 = new TF1("fitRatio1","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2])) / sqrt(((TMath::Sign(1,[4])*sq([4]/x))+(sq([5])*(x^([7]-1))))+sq([6]))",ptMin,ptMax);
  }
  int nPar = fitData1->GetNpar();
  for(int i = 0; i < nPar; ++i) {
    fitRatio1->SetParameter(i,fitData1->GetParameter(i));
    fitRatio1->SetParError(i,fitData1->GetParError(i));
    fitRatio1->SetParameter(nPar+i,fitMC->GetParameter(i));
    fitRatio1->SetParError(nPar+i,fitMC->GetParError(i));
  }
  fitRatio1->SetLineWidth(1);
  fitRatio1->SetLineColor(fitData1->GetLineColor());


  TF1* fitRatio2 = 0;
  if( interpolationMode == 0 ) {
    fitRatio2 = new TF1("fitRatio2","sqrt(sq([0]/x) + sq([1])/x + sq([2]) + sq([3])) / sqrt(sq([4]/x) + sq([5])/x + sq([6]))",ptMin,ptMax);
    for(int i = 0; i < 3; ++i) {
      fitRatio2->SetParameter(i,fitMC->GetParameter(i));
      fitRatio2->SetParError(i,fitMC->GetParError(i));
      fitRatio2->SetParameter(4+i,fitMC->GetParameter(i));
      fitRatio2->SetParError(4+i,fitMC->GetParError(i));
    }
    fitRatio2->SetParameter(3,fitData2->GetParameter(3));
    fitRatio2->SetParError(3,fitData2->GetParError(3));
  } else if( interpolationMode == 1 ) {
    fitRatio2 = new TF1("fitRatio2","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2])) / sqrt(((TMath::Sign(1,[4])*sq([4]/x))+(sq([5])*(x^([7]-1))))+sq([6]))",ptMin,ptMax);
    for(int i = 0; i < fitData2->GetNpar(); ++i) {
      fitRatio2->SetParameter(i,fitData2->GetParameter(i));
      fitRatio2->SetParError(i,fitData2->GetParError(i));
      fitRatio2->SetParameter(fitData2->GetNpar()+i,fitMC->GetParameter(i));
      fitRatio2->SetParError(fitData2->GetNpar()+i,fitMC->GetParError(i));
    }
  } 
  fitRatio2->SetLineWidth(1);
  fitRatio2->SetLineColor(fitData2->GetLineColor());


  // fit ratio directly (include gamma jet data)
  TF1* fitRatio3 = 0;
  if( interpolationMode == 0 ) {
    fitRatio3 = new TF1("fitRatio3","sqrt(sq([0]/x) + sq([1])/x + sq([2]+[3])) / sqrt(sq([0]/x) + sq([1])/x + sq([2]))",ptMin,ptMax);
    fitRatio3->SetParameter(0,1.);
    fitRatio3->SetParameter(1,0.8);
    fitRatio3->SetParameter(2,0.03);
  } else if( interpolationMode == 1 ) {
    fitRatio3 = new TF1("fitRatio3","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2])) / sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1)))))",ptMin,ptMax);
    fitRatio3->SetParameter(0,1.);
    fitRatio3->SetParameter(1,0.8);
    fitRatio3->SetParameter(2,0.02);
    fitRatio3->SetParameter(3,0.);
  } 
  fitRatio3->SetLineWidth(1);
  fitRatio3->SetLineColor(kGreen);


  // constant ratio
  TF1* fitRatio4 = new TF1("fitRatio4","pol0",ptMin,ptMax);
  fitRatio4->SetLineWidth(1);
  fitRatio4->SetLineColor(kBlue);


  // Combined dijet and gammajet ratio graph
  TGraphAsymmErrors* gCombined = util::HistOps::combineTGraphs(gRatioGJ,gRatio);
  gCombined->SetMarkerStyle(22);
  gCombined->Fit(fitRatio3,"0R");
  gCombined->Fit(fitRatio4,"0R");



  std::cout << endl;
  std::cout << endl;
  for(int i = 0; i < fitRatio3->GetNpar(); ++i) {
    std::cout << i << ": " << fitRatio3->GetParameter(i) << std::endl;
  }
  
  



  //   TF1* fitRatio2 = new TF1("fitRatio2","[0] - exp(x/[1])",ptMin,ptMax);
  //   fitRatio2->SetParameter(0,1.2);
  //   fitRatio2->SetParameter(1,-100.);
  //   fitRatio2->SetLineWidth(1);
  //   fitRatio2->SetLineColor(2);
  //   fitRatio2->SetLineStyle(1);  

  //   std::cout << "\nfitRatio1: " << std::endl;
  //   gRatio->Fit(fitRatio1,"R0");
  //   std::cout << "\n\nfitRatio2: " << std::endl;
  //   gRatio->Fit(fitRatio2,"R0");

  //   std::cout << "\nPar0 = " << fitRatio2->GetParameter(0) << " (" << fitRatio2->GetParameter(0)-fitRatio2->GetParError(0) << ", " << fitRatio2->GetParameter(0)+fitRatio2->GetParError(0) << ")" << std::endl;
  //   std::cout << "Par1 = " << fitRatio2->GetParameter(1) << " (" << fitRatio2->GetParameter(1)-fitRatio2->GetParError(1) << ", " << fitRatio2->GetParameter(1)+fitRatio2->GetParError(1) << ")" << std::endl;

  

  // Resolution plots
  TPaveText* label = util::LabelFactory::createPaveText(1);
  if( pf_ ) label->AddText("AK5 PF-Jets,  |#eta| < 1.1");
  else label->AddText("AK5 Calo-Jets,  |#eta| < 1.1");

  TLegend* legRes = util::LabelFactory::createLegendColWithOffset(5,-0.8,1);
  legRes->AddEntry(gMC,"Measurement (MC)","P");  
  legRes->AddEntry(fitMC,"Interpolation (MC)","L");  
  legRes->AddEntry(gData,"Measurement (Data)","P");
  legRes->AddEntry(fitData1,"Interpolation (Data)","L");  
  legRes->AddEntry(fitData2,"MC + C^{2} (Data)","L");  

  TLegend* legRatio = util::LabelFactory::createLegendColWithOffset(2,-0.8,1);
  legRatio->AddEntry(fitRatio1,"Ratio of interpolations","L");
  legRatio->AddEntry(fitRatio2,"MC + C^{2}","L");



  TCanvas* canRes = new TCanvas("canRes","Resolution",500,500);
  canRes->cd();
  TH1* hFrame = util::HistOps::createTH1D("hFrame",1000,ptMin,ptMax,"p^{ref}_{T}","GeV","#sigma / p^{ref}_{T}");
  hFrame->GetYaxis()->SetRangeUser(1E-3,0.29);
  hFrame->GetXaxis()->SetMoreLogLabels();
  //   hFrame->Draw();
  //   gMC->Draw("PE1same");
  //   gData->Draw("PE1same");
  //   label->Draw("same");
  //   legRes->Draw("same");
  //   canRes->SetLogx();
  //   canRes->SaveAs(outNamePrefix+"Resolution.eps","eps");

  canRes->cd();
  hFrame->Draw();
  gMC->Draw("PE1same");
  gData->Draw("PE1same");
  fitMC->Draw("same");
  fitData1->Draw("same");
  fitData2->Draw("same");
  label->Draw("same");
  legRes->Draw("same");
  canRes->SetLogx();
  canRes->SaveAs(outNamePrefix+"ResolutionFit.eps","eps");


  TCanvas* canRatio = new TCanvas("canRatio","Ratio",500,500);
  canRatio->cd();
  TH1* hRatioFrame = util::HistOps::createRatioFrame(hFrame,"#sigma(Data)  /  #sigma(MC)",0.81,1.69);
  hRatioFrame->GetXaxis()->SetMoreLogLabels();
  hRatioFrame->Draw();
  fitRatio1->Draw("same");
  fitRatio2->Draw("same");
  fitRatio3->Draw("same");
  fitRatio4->Draw("same");
  //  gCombined->Draw("PE1same");
  gRatio->Draw("PE1same");
  gRatioGJ->Draw("PE1same");
  label->Draw("same");
  legRatio->Draw("same");
  canRatio->SetLogx();
  canRatio->SaveAs(outNamePrefix+"Ratio.eps","eps");
}



