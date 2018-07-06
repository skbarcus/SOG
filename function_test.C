#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <math.h>

void function_test() 
{

  Double_t xval[100];
  Double_t yval[100];

  for(Int_t i=0;i<100;i++)
    {
      xval[i] = i;
      yval[i] = 2.*i+11.;
    }

  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  graph = new TGraph(100,xval,yval);
  graph->SetMarkerColor(1);
  graph->SetMarkerSize(1.);
  graph->SetMarkerStyle(20);
  graph->Draw("");

  Double_t fn1(Double_t *x, Double_t *par)
  {
    Double_t val = 0.; 
    
    val = x[0]*par[0]+par[1];
    
    return val;
  }

  TF1 *f1 = new TF1("f1",fn1, 0, 100,2);
  //graph->Fit(f1,"0");

  Double_t fn2(Double_t *y, Double_t *par)
  {
    Double_t val = 0.; 
    
    val = f1->Eval(y[0]);
    
    return val;
  }

  TF1 *f2 = new TF1("f2",fn2, 0, 100,2);
  graph->Fit(f2,"");

}
