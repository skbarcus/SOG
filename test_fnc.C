#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

Double_t num = 5.;
Double_t array[3] = {2.,4.,8.};

void test_fnc()
{

  //Double_t n = 5.;
  TF1 *rand = new TF1("rand","x",0,10);

  TH1 *fnc = new TH1D("fcn", "fcn", 10., 0., 5.);

  for(Int_t i=0;i<11;i++)
    {
      fnc->SetBinContent(i+1.,rand->GetRandom());
    }

  Double_t fitf(Double_t *x,Double_t *par)//,Double_t*R)
  {
    //Double_t n = 5.;
    Double_t fitval = num*par[0]*x[0];
    Double_t y = num;
    cout<<"y = "<<y<<"   par[0] = "<<par[0]<<"   num = "<<num<<endl;
    cout<<array[0]<<" "<<array[1]<<" "<<array[2]<<endl;
    return fitval;
  }

  TF1 *f1 = new TF1("fit",fitf,0.,1.,1.0);

  fit->SetParameter(0,1.);
  //fit->SetParameter(1,2.);

  fcn->Fit("fit","");
}
