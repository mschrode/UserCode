#include "Plot.h"
#include "Event.h"

#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>

//Function to fill a given histogram 'h'
void FillHist(TH1*h, const std::vector<Event*>& evts, float(*var)(const Event*) )
{
  //loop over all events 'evts'
  for (std::vector<Event*>::const_iterator it=evts.begin(); it!=evts.end(); ++it) {
    //Apply a selection on HT:
    if ( (*it)->HT()<350.) continue;
    
    //Fill the current event '*it' with weight '(*it)->EvtWgt' into the histogram 'h'
    h->Fill( var(*it), (*it)->EvtWgt );
  }  
}

//Global functions used by the function 'FillHist'
float HT(const Event * evt){ return evt->HT(); }
float MHT(const Event * evt){ return evt->MHT(); }
float Jet1Pt(const Event * evt){ return evt->recoJetPt[0]; }
float Jet2Pt(const Event * evt){ return evt->recoJetPt[1]; }
float Jet3Pt(const Event * evt){ return evt->recoJetPt[2]; }
float JetMultiplicity(const Event * evt){ int N=0; for (int i=0;i<evt->NrecoJet;++i) if (evt->recoJetPt[i]>50.) ++N;return N; }


int main(int argc, char* argv[])
{
  //Make the standard root histograms less ugly
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  TCanvas c1;
  c1.SetLogy(1); //Set log y-axis

  //STL container for all events
  std::vector<Event*> eventsA, eventsB;

  //Read the events from the root n-tuple into the event containers
  ReadEvents( "data/QCDcontrol_data.root", eventsA);
  ReadEvents( "data/QCDcontrol_rs.root", eventsB );

  //Fill histogram HT, plot both samples into one plot
  TLegend * leg = new  TLegend(0.50, 0.50, 0.9, 0.9);
  { TH1F* ht_A = new TH1F("ht_A",";HT (GeV);events",100, 0,2000);
    TH1F* ht_B = (TH1F*)ht_A->Clone();  
    FillHist(ht_A, eventsA, HT);
    FillHist(ht_B, eventsB, HT);
    ht_B->SetMarkerStyle(8);
    ht_B->Scale( ht_A->Integral()/ht_B->Integral() ); //normalize histogram B to hist A
    if (ht_B->GetMaximum()>ht_A->GetMaximum()) ht_A->SetMaximum( ht_B->GetMaximum() );
    ht_A->Draw("h");        //draw has histogram
    ht_B->Draw("pe, same"); //draw as points with errors into same canvas
    leg->SetBorderSize(1);leg->SetFillColor(0);
    leg->SetHeader("Comparison Plots");
    leg->AddEntry(ht_A, "Sample A", "l");
    leg->AddEntry(ht_B, "Sample B", "pe");
    leg->Draw();           //draw legend
    c1.SaveAs("HT.pdf");   //write plot to PDF file
  }  

  return 0;
}
