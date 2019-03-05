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
#include "TMinuit.h"

#include <TMath.h>
#include "Math/IFunction.h"
#include <cmath>

#include "Math/Functor.h"
#include "Math/RichardsonDerivator.h"

Double_t pi = 3.141592654;

Double_t r = 2.;
Double_t xmin = 0.;
Double_t xmax = 10.;
Double_t xmin_sph = 0.;
Double_t xmax_sph = 3.;

Double_t hard_sph(Double_t *x,Double_t *par)
{
    Bool_t reject;
    reject = kTRUE;
    //reject = kFALSE;
    if (reject && x[0] > xmin_sph && x[0] < xmax_sph) 
      {
	TF1::RejectPoint();
	return 0;
      }
    Double_t fitval = 1/pow(2*pi,3.)*(3./4.)*pi*pow(r,3.)*;
    return fitval;
}

void FF_Demo() 
{
  TF1 *fhard_sph = new TF1("fhard_sph",hard_sph,xmin,xmax,0);

  TCanvas* c1=new TCanvas("c1");
  c1->Divide(1,3);
  
  c1->cd(1);
  //A function to sample
  TF1 *fsin = new TF1("fsin", "sin(x)+sin(2*x)+sin(0.5*x)+1", 0, 4*TMath::Pi());
  fsin->Draw();
  
  Int_t n=25;
  TH1D *hsin = new TH1D("hsin", "hsin", n+1, 0, 4*TMath::Pi());
  Double_t x;
  
  //Fill the histogram with function values
  for (Int_t i=0; i<=n; i++){
    x = (Double_t(i)/n)*(4*TMath::Pi());
    hsin->SetBinContent(i+1, fsin->Eval(x));
  }
  hsin->Draw("same");
  fsin->GetXaxis()->SetLabelSize(0.05);
  fsin->GetYaxis()->SetLabelSize(0.05);
  
  c1->cd(2);
  //Compute the transform and look at the magnitude of the output
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = hsin->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");
  hm->Draw();
  //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
  //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
  
  hm->SetStats(kFALSE);
  hm->GetXaxis()->SetLabelSize(0.05);
  hm->GetYaxis()->SetLabelSize(0.05);

  Double_t re, im;
  //That's the way to get the current transform object:
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  //Use the following method to get just one point of the output
  fft->GetPointComplex(0, re, im);

  //Use the following method to get the full output:
  Double_t *re_full = new Double_t[n];
  Double_t *im_full = new Double_t[n];
  fft->GetPointsComplex(re_full,im_full);
  
  c1->cd(3);
  //Now let's make a backward transform:
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  //Let's look at the output
  hb = TH1::TransformHisto(fft_back,hb,"Re");
  hb->SetTitle("The backward transform result");
  hb->Draw();
  //NOTE: here you get at the x-axes number of bins and not real values
  //(in this case 25 bins has to be rescaled to a range between 0 and 4*Pi;
  //also here the y-axes has to be rescaled (factor 1/bins)
  hb->SetStats(kFALSE);
  hb->GetXaxis()->SetLabelSize(0.05);
  hb->GetYaxis()->SetLabelSize(0.05);
  delete fft_back;
  fft_back=0;
}
