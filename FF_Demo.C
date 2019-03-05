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

#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"

Double_t pi = 3.141592654;

Double_t range = 1000.;
Int_t n=10000;
Double_t r = 2.5;
Double_t xmin = 0.;
Double_t xmax = 10.;
Double_t xmin_sph = 0.;
Double_t xmax_sph = 5.;

Double_t hard_sph(Double_t *x,Double_t *par)
{
    Bool_t reject;
    //reject = kTRUE;
    //reject = kFALSE;
    //if (reject && x[0] > xmin_sph && x[0] < xmax_sph) 
    if (x[0] < xmin_sph || x[0] > xmax_sph) 
      {
	TF1::RejectPoint();
	return 0;
      }
    Double_t fitval = 1/pow(2*pi,3.)*(3./4.)*pi*pow(r,3.);
    return fitval;
}

void FF_Demo() 
{
  //TF1 *fblank = new TF1("fblank","0",xmin,xmax,0);
  TF1 *fhard_sph = new TF1("fhard_sph",hard_sph,xmin,xmax,0);

  TCanvas* c1=new TCanvas("c1");
  c1->Divide(1,2);
  //fblank->GetHistogram()->SetTitle("Hard Sphere Charge Distribution");
  c1->cd(1);
  //fblank->Draw();

  //A function to sample
  TF1 *fsin = new TF1("fsin", hard_sph, 0, range,0);
  //fsin->GetHistogram()->SetTitle("Hard Sphere Charge Distribution");
  fsin->GetHistogram()->GetYaxis()->SetLabelOffset(10);
  fsin->GetHistogram()->GetXaxis()->SetLabelOffset(10);
  //fsin->Draw();
  
  TH1D *hsin = new TH1D("hsin", "Hard Sphere Charge Distribution", n+1, 0, range);
  Double_t x;
  
  //Fill the histogram with function values
  for (Int_t i=0; i<=n; i++)
    {
      x = (Double_t(i)/n)*(range);
      hsin->SetBinContent(i+1, fsin->Eval(x));
    }
  hsin->Draw();
  //gStyle->SetOptStat(11);
  //gStyle->SetTitleX(0.5);
  //gStyle->SetTitleAlign(23);
  gStyle->SetOptStat(0);
  //hsin->SetTitle("Hard Sphere Charge Distribution");
  hsin->GetXaxis()->SetTitle("r (arbitrary units)");
  hsin->GetXaxis()->CenterTitle(true);
  hsin->GetXaxis()->SetLabelSize(0.05);
  hsin->GetXaxis()->SetTitleSize(0.06);
  hsin->GetXaxis()->SetTitleOffset(0.75);
  hsin->GetYaxis()->SetTitle("Charge (arbitrary units)");
  hsin->GetYaxis()->CenterTitle(true);
  hsin->GetYaxis()->SetLabelSize(0.05);
  hsin->GetYaxis()->SetTitleSize(0.06);
  hsin->GetYaxis()->SetTitleOffset(0.75);
  //hsin->GetXaxis()->SetLabelOffset(10);
  //hsin->GetYaxis()->SetLabelOffset(10);
  hsin->GetXaxis()->SetNdivisions(0);
  hsin->GetYaxis()->SetNdivisions(0);
  //hsin->GetXaxis()->SetLabelSize(0);
  //hsin->GetYaxis()->SetLabelSize(0);
  
  c1->cd(2);
  //Compute the transform and look at the magnitude of the output
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = hsin->FFT(hm, "MAG");
  hm->SetTitle("Form Factor of Hard Sphere Charge Distribution");
  hm->GetXaxis()->SetTitle("q (arbitrary units)");
  hm->GetXaxis()->CenterTitle(true);
  hm->GetXaxis()->SetLabelSize(0.05);
  hm->GetXaxis()->SetTitleSize(0.06);
  hm->GetXaxis()->SetTitleOffset(0.75);
  hm->GetYaxis()->SetTitle("|F_{ch}(q^{2})|");
  hm->GetYaxis()->CenterTitle(true);
  hm->GetYaxis()->SetLabelSize(0.05);
  hm->GetYaxis()->SetTitleSize(0.06);
  hm->GetYaxis()->SetTitleOffset(0.75);
  //hm->GetXaxis->SetTitle("q");
  hm->Draw();
  //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
  //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
  
  hm->SetStats(kFALSE);
  hm->GetXaxis()->SetNdivisions(0);
  hm->GetYaxis()->SetNdivisions(0);
  //hm->GetXaxis()->SetLabelSize(0.05);
  //hm->GetYaxis()->SetLabelSize(0.05);
  //hm->GetXaxis()->SetLabelOffset(10);
  //hm->GetYaxis()->SetLabelOffset(10);

  /*
  //Now rescale the x and y axes back to their correct values after the FT.
  TH1F *hFT_rescaled = new TH1F("hFT_rescaled","hFT_rescaled",1000,0.,n/range);
  //cout<<hm->GetMaximumBin()<<endl;
  cout<<hm->GetXaxis()->GetBinCenter(hm->GetMaximumBin())<<endl;
  for(Int_t i=0;i<1000;i++)
    {
      hFT_rescaled->SetBinContent(i+1, hm->GetBinContent(i+1)*1/pow(n,0.5));
    }

  c1->cd(3);
  hFT_rescaled->Draw();
  */
  /*
  c1->cd(3);
  ScaleXaxis(hm,ScaleX);
  ScaleYaxis(hm,ScaleY);
  hm->Draw();
  */
  /*
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
  */
}
