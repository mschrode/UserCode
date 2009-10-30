//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventProcessor.cc,v 1.2 2009/07/23 13:46:20 mschrode Exp $
//   
#include "EventProcessor.h"

#include "ConfigFile.h"
#include "Parameters.h"
#include "CalibData.h"

#include "TH1F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPostScript.h"

#include <iostream>
#include <algorithm>

using std::cout;
using std::endl;



EventProcessor::EventProcessor(const std::string& configfile, TParameters* param)
  : par_(param),etCutOnGamma_(0),etCutOnJet_(0),flattenSpectra_(0)
{  
  ConfigFile config(configfile.c_str());

  flattenSpectra_ = config.read<int>("Flatten Spectra",0); 
  
  //last minute kinematic cuts
  etCutOnJet_   = config.read<double>("Et cut on jet",0.0); 
  etCutOnGamma_ = config.read<double>("Et cut on gamma",0.0); 
  //relative sample weight 
  relWeight_[2]        = config.read<double>("Gamma-Jet weight",1.0);
  relWeight_[1]        = config.read<double>("Track-Tower weight",1.0);
  relWeight_[3]        = config.read<double>("Track-Cluster weight",1.0);
  relWeight_[4]        = config.read<double>("Di-Jet weight",1.0);
  relWeight_[5]        = config.read<double>("Multi-Jet weight",1.0);
  //relWeight_[2]        = config.read<double>("Z-Jet weight",1.0);
}
 
EventProcessor::~EventProcessor()
{
}
  

int EventProcessor::process(std::vector<TData*>& data)
{
  if(flattenSpectra_) {
    std::cout << "Flatten spectra for ";
    std::cout << flattenSpectra(data) << " events\n";
  }
  //BalanceSpectra(data);
  return data.size();
}
 
bool EventProcessor::NotBalancedRejection::operator()(TData *d) {
  bool result = false;
  TAbstractData *ad = dynamic_cast<TAbstractData*>(d);
  if(! ad) return true;
  
  if(ad->GetType()!=GammaJet ||
     ad->GetMess()->pt<_min ||
     ad->GetMess()->pt>_max ||
     ad->GetMess()->pt==0.0 ) result = true;
  else
    result = (1.0-ad->GetTruth()/ad->GetMess()->pt) >
      _cut[(int)(ad->GetMess()->pt-_min)];
  return result;
}

int EventProcessor::getSpectraBin(double m1, double m2=0., double m3=0.)
{
  //pt
  int bin1, bins1 = 7000; 
  double min1    = 0.;
  double max1    = 7000.;
  if      (m1<min1) bin1=0;
  else if (m1>max1) bin1=bins1+1;
  else              bin1=(int)(((m1-min1)/max1)*(double)bins1);
  //eta
  int bin2=0, bins2 = par_->GetEtaGranularityJet()*par_->GetPhiGranularityJet();
  double min2    = 0.;
  double max2    = 82.;
  if      (m2<min2) bin2=0;
  else if (m2>max2) bin2=bins2;
  else              bin2=(int)(((m2-min2)/max2)*(double)bins2);
  //EMF
  int bin3=0, bins3 = 100;
  double min3    = 0.;
  double max3    = 1.;
  if      (m3<min3) bin3=0;
  else if (m3>max3) bin3=bins3;
  else              bin3=(int)(((m3-min3)/max3)*(double)bins3);

  return bin1 *bins2*bins3 + bin2 *bins3 + bin3;
}

double EventProcessor::gaussStep(double *x, double *par)
{
  return par[2]/(par[1]*2.5)*exp(-(x[0]-par[0])*(x[0]-par[0])/(2.0*par[1]*par[1]))
    *
    (1.0-1.0/(1.0+exp((par[3]-x[0])/par[4]) ) );  	
}

void EventProcessor::fitWithoutBottom(TH1 * hist, TF1 * func, double bottom)
{
  TH1F * result=(TH1F*)hist->Clone();
  double maximum = hist->GetMaximum();
  int min=0, max=0;
  for (int i=0; i<hist->GetNbinsX(); ++i)
    if (hist->GetBinContent(i)>bottom*maximum){
      result->SetBinContent(i,hist->GetBinContent(i));
      max=i;
      if (min==0.) min=i;
    }      
  func->SetRange(hist->GetXaxis()->GetXmin()+(double)min/(double)hist->GetNbinsX()*(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin()),
		 hist->GetXaxis()->GetXmin()+(double)max/(double)hist->GetNbinsX()*(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin()));
  result->Fit("gauss_step","LLQNO","");
}

int EventProcessor::flattenSpectra(std::vector<TData*>& data)
{
  int numProcEvts = 0; // Number of processed events

  double allweights=0;
  //@@ Replace "7" with some more meaningful variable name!!!
  for(int i=0;i<7;++i)
    allweights += relWeight_[i];

  map<int,double> weights[7];
  double tot[7];
  double alltotal=0;
  for(int i=0;i<7;++i)
    tot[i]=0.;
 
  for (DataConstIter it = data.begin(); it!=data.end(); ++it) {
    TAbstractData* dt = dynamic_cast<TAbstractData*>(*it);
    if(! dt) continue;
    if( dt->GetType() == ParLimit ) continue;
    numProcEvts++;
    alltotal+=dt->GetWeight();
    for (int type=0; type<7; ++type){
      if (dt->GetType()!=type) continue;
      if (dt->GetType()==InvMass) continue;
      double em = 0.;
      double had = 0.;
      int index=0;
      double min_tower_dr = 10.;
      //double ptmax=0;
       
      TLorentzVector Ljet(0.,0.,0.,0.);
      Ljet.SetPtEtaPhiE(dt->GetMess()->pt,dt->GetMess()->eta,dt->GetMess()->phi,dt->GetMess()->E);
      for(std::vector<TAbstractData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
	em  += (*t)->GetMess()->EMF;
	had += (*t)->GetMess()->HadF;
	had += (*t)->GetMess()->OutF;
	TLorentzVector Ltower(0.,0.,0.,0.);
	Ltower.SetPtEtaPhiE((*t)->GetMess()->pt,(*t)->GetMess()->eta,(*t)->GetMess()->phi,(*t)->GetMess()->E);
	double dr = Ltower.DeltaR(Ljet);
	if (dr<min_tower_dr) {
	  index = (*t)->GetIndex();
	  min_tower_dr = dr;
	}
      }
      //int bin = getSpectraBin( dt->GetScale(), index, em/(em+had)  );
      //int bin = getSpectraBin( dt->GetScale(), index );
      int bin = getSpectraBin( dt->GetScale() );
      double error = 1.;//dt->GetParametrizedErr( &dt->GetMess()->pt );
      weights[type][bin]+=dt->GetWeight()/error;
      tot[type]+=dt->GetWeight()/error;
    }
  }
  for (int type=0; type<7; ++type){
    if (tot[type]!=0.)
      for (DataIter it = data.begin(); it!=data.end(); ++it) {
	TAbstractData* dt = dynamic_cast<TAbstractData*>(*it);
	if(! dt) continue;
	if (dt->GetType()!=type) continue;
	if (dt->GetType()==InvMass) continue;
	 
	 
	double em = 0;
	double had = 0;
	int index=0;
	double min_tower_dr = 10.;
	TLorentzVector Ljet(0,0,0,0);
	 
	Ljet.SetPtEtaPhiE(dt->GetMess()->pt,dt->GetMess()->eta,dt->GetMess()->phi,dt->GetMess()->E);
	for(std::vector<TAbstractData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
	  em  += (*t)->GetMess()->EMF;
	  had += (*t)->GetMess()->HadF;
	  had += (*t)->GetMess()->OutF;
	  TLorentzVector Ltower(0.,0.,0.,0.);
	  Ltower.SetPtEtaPhiE((*t)->GetMess()->pt,(*t)->GetMess()->eta,(*t)->GetMess()->phi,(*t)->GetMess()->E);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    index = (*t)->GetIndex();
	    min_tower_dr = dr;
	  }
	}
	//int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had) );
	//int bin = GetSpectraBin( dt->GetScale(), index );
	int bin = getSpectraBin( dt->GetScale() );
	//dt->SetWeight(1);
	 
	 
	//dt->SetWeight(dt->GetWeight()/weights[type][bin] * (double(tot[type]) / double(weights[type].size())));
	 
	dt->setWeight(dt->GetWeight()/weights[type][bin] * (alltotal / double(weights[type].size()) )  * relWeight_[type] / allweights );
	 
	 
	//dt->SetWeight((1./weights[bin]) * (double(tot) / weights.size()));
      }
     
    /*
    // Old version using leading tower bin as jet bin
    for (DataConstIter it = data.begin(); it!=data.end(); ++it) {
    if (dt->GetType()!=type) continue;
    double em = 0;
    double had = 0;
    int index=0;
    double ptmax=0;
    for(std::vector<TData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
    em  += (*t)->GetMess()[1];
    had += (*t)->GetMess()[2];
    had += (*t)->GetMess()[3];
    if ((*t)->GetMess()[1]>ptmax) {
    ptmax=(*t)->GetMess()[1];
    index=(*t)->GetIndex();
    }
    }
    //int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had)  );
    int bin = GetSpectraBin( dt->GetScale(), index );
    //int bin = GetSpectraBin( dt->GetScale() );
    weights[bin]+=dt->GetWeight();
    tot+=dt->GetWeight();
    }
	
    if (tot!=0.)
    for (DataIter it = data.begin(); it!=data.end(); ++it) {
    if (dt->GetType()!=type) continue;
    
    double em = 0;
    double had = 0;
    int index=0;
    double ptmax=0;
    for(std::vector<TData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
    em  += (*t)->GetMess()[1];
    had += (*t)->GetMess()[2];
    had += (*t)->GetMess()[3];
    if ((*t)->GetMess()[1]>ptmax) {
    ptmax=(*t)->GetMess()[1];
    index=(*t)->GetIndex();
    }
    }

    //int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had) );
    int bin = GetSpectraBin( dt->GetScale(), index );
    //int bin = GetSpectraBin( dt->GetScale() );
    //dt->SetWeight(1);
    //dt->SetWeight(1./weights[bin]);
    dt->SetWeight((1./weights[bin]) * (double(tot) / weights.size()));
    }
    */

  } 

  return numProcEvts;
}
  
//further weighting.............................................................
void EventProcessor::balanceSpectra(std::vector<TData*>& data)
{
  cout<<"...further weighting"<<endl;
  double min = etCutOnGamma_;
  double max = 100.; //GeV
  int nbins = (int)(max-min);//one bin per GeV
  if (nbins<2) return;
  double EMF[nbins];
  double TOT[nbins];
  
  TCanvas * c1 = new TCanvas("controlplots","",600,600);
  TPostScript ps("balance_spectra.ps",111);
  TH1F * gauss_forpt[nbins];
  TH1F * gauss_forpt_truth[nbins];
  gauss_forpt[0] = new TH1F("hgauss","pT bin[20..21GeV];#frac{pT jet - pT truth}{pT jet}",600,-3,3);
  gauss_forpt_truth[0] = new TH1F("hgauss_truth","pT bin[20..21GeV];pT truth",400,0,200);
  char * name = new char[100];
  for(int i = 1 ; i < nbins; ++i) {
    gauss_forpt[i] = (TH1F*)gauss_forpt[0]->Clone();
    sprintf(name,"pT bin[%d..%dGeV]",(int)min+i,(int)min+i+1);
    gauss_forpt[i]->SetTitle(name);
    gauss_forpt_truth[i] = (TH1F*)gauss_forpt_truth[0]->Clone();
    sprintf(name,"pT bin[%d..%dGeV]",(int)min+i,(int)min+i+1);
    gauss_forpt_truth[i]->SetTitle(name);
  }

  cout<<"...fill truth histograms for each jet-pT bin"<<endl;
  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data.begin(); 
        i != data.end() ; ++i )  {
    TAbstractData* jg = dynamic_cast<TAbstractData*>(*i);
    if(! jg) continue;
    if (jg->GetType()!=GammaJet) continue;
    
    //double etjetcor = jg->GetParametrizedMess();
    if(jg->GetMess()->pt>min && jg->GetMess()->pt<max) {
      gauss_forpt[(int)(jg->GetMess()->pt-min)]->Fill( (jg->GetMess()->pt-jg->GetTruth())/jg->GetMess()->pt,jg->GetWeight() );
      gauss_forpt_truth[(int)(jg->GetMess()->pt-min)]->Fill( jg->GetTruth(),jg->GetWeight() );
      EMF[(int)(jg->GetMess()->pt-min)] += jg->GetWeight()*jg->GetMess()->EMF;
      TOT[(int)(jg->GetMess()->pt-min)] += jg->GetWeight()*jg->GetMess()->pt;
    }      
  }

  cout<<"...fit the truth distributions"<<endl;
  double edge;
  TF1 * f = new TF1("gauss_step",gaussStep,-3,3,5);
  double * cuts = new double[nbins];
  TText * text = new TText();
  text->SetTextSize(0.03);
  text->SetTextColor(2);

  for(int i = 0; i < nbins; ++i) {
    if ( (i+1)%(nbins/10)==0) cout << (100*(i+1)/nbins)<<"% events weighted" << endl;  
    //TF1 *f=0;
    //gauss_forpt[i]->Fit("gaus","LLQNO","");
    //f = (TF1*)gROOT->GetFunction("gaus")->Clone();
    edge = 1.0-etCutOnGamma_/(((double)i)+min+0.5);
    f->SetParameters(-1.,2.0,3.0, edge, 0.0001);
    f->FixParameter(3, edge);
    f->FixParameter(4, 0.0001);
    fitWithoutBottom(gauss_forpt[i], f);
    //bla->Fit("gauss_step","LLQNO","");
    //delete bla;

    gauss_forpt[i]->Draw("h");
    f->SetLineColor(2);
    f->Draw("same");
    sprintf(name,"mean %f",f->GetParameter(0));
    text->DrawText(1.4,0.7*gauss_forpt[i]->GetMaximum(),name);

    sprintf(name,"average truth %f", (double)i+0.5+min-((double)i+0.5+min)*f->GetParameter(0));
    text->DrawText(0.50,0.65*gauss_forpt[i]->GetMaximum(),name);
    c1->Draw();

    if (TOT[i]!=0.0) EMF[i] = EMF[i]/TOT[i];
    sprintf(name,"average/truth %+f",(double)i+0.5+min-((double)i+0.5+min)*f->GetParameter(0) / 
	    (1.3*((double)i+0.5+min)   )); //-EMF[i])  );
    text->DrawText(0.50,0.60*gauss_forpt[i]->GetMaximum(),name);
    c1->Draw();


    gauss_forpt_truth[i]->Draw("h");
    c1->Draw();
    
    cout<<"bin "<<i
	<<": mean="<<f->GetParameter(0)
	<<", sigma="<<f->GetParameter(1)
	<<", height="<<f->GetParameter(2)
	<<", edge("<<edge<<")="<<f->GetParameter(3)
	<<", width-edge="<<f->GetParameter(4)
	<<endl;
    cuts[i] = f->GetParameter(0)-fabs(f->GetParameter(0)-f->GetParameter(3))+0.2;
  }
  delete f;
  ps.Close();
 
  DataIter beg = partition(data.begin(), data.end(), 
                           NotBalancedRejection(cuts, min, max));
  cout<<"...remove " << int(data.end()-beg) << " events which are not 'balanced'";
  for(DataIter i = beg ; i != data.end() ; ++i) {
    cout<<".";
    delete *i;
  }
  cout<<endl;
  data.erase(beg,data.end());

  cout<<"...cleaning up"<<endl;
  for(int i = 0; i < nbins; ++i){
    delete gauss_forpt[i];
  }  
  delete [] cuts;
  delete name;
}
