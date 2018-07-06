#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

//gRoot->Reset();
Double_t ga2 = pow(0.8*pow(2./3., 0.5), 2.);
Double_t m = 3.06;
//Double_t R[12] = {0.1*m,.5*m,.9*m,1.3*m,1.6*m,2.*m,2.4*m,2.9*m,3.4*m,4.*m,4.6*m,5.2*m};
//Double_t R[12] = {0.1,.5,.9,1.3,1.6,2.,2.4,2.9,3.4,4.,4.6,5.2};
//Double_t R[12] = {0.};
//Double_t R[12] = {3.,6.,9.,12.,15.,18.,21.,24.,27.,30.,33.,36.};
//Double_t R[12] = {3.05991, 6.11991, 9.17988, 12.2398, 15.2998, 18.3597, 21.4197, 24.4796, 27.5395, 30.5995, 33.6594, 36.7194};    //Roughly even increments of 3.06.
Double_t R[12] = {1.*m, 2.*m, 3.*m, 4.*m, 5.*m, 6.*m, 7.*m, 8.*m, 9.*m, 10.*m, 11.*m, 12.*m};
Int_t userand = 0;
Int_t npar = 12;
Double_t fitmin = 0.;
Double_t fitmax = 1.;
Double_t fitrange = abs(fitmax-fitmin);

void fft_sphere()
{
//This tutorial illustrates the Fast Fourier Transforms interface in ROOT.
//FFT transform types provided in ROOT:
// - "C2CFORWARD" - a complex input/output discrete Fourier transform (DFT)
//                  in one or more dimensions, -1 in the exponent
// - "C2CBACKWARD"- a complex input/output discrete Fourier transform (DFT)
//                  in one or more dimensions, +1 in the exponent
// - "R2C"        - a real-input/complex-output discrete Fourier transform (DFT)
//                  in one or more dimensions,
// - "C2R"        - inverse transforms to "R2C", taking complex input
//                  (storing the non-redundant half of a logically Hermitian array)
//                  to real output
// - "R2HC"       - a real-input DFT with output in Â¡ÃˆhalfcomplexÂ¡Ã‰ format,
//                  i.e. real and imaginary parts for a transform of size n stored as
//                  r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
// - "HC2R"       - computes the reverse of FFTW_R2HC, above
// - "DHT"        - computes a discrete Hartley transform
// Sine/cosine transforms:
//  DCT-I  (REDFT00 in FFTW3 notation)
//  DCT-II (REDFT10 in FFTW3 notation)
//  DCT-III(REDFT01 in FFTW3 notation)
//  DCT-IV (REDFT11 in FFTW3 notation)
//  DST-I  (RODFT00 in FFTW3 notation)
//  DST-II (RODFT10 in FFTW3 notation)
//  DST-III(RODFT01 in FFTW3 notation)
//  DST-IV (RODFT11 in FFTW3 notation)
//First part of the tutorial shows how to transform the histograms
//Second part shows how to transform the data arrays directly
//Authors: Anna Kreshuk and Jens Hoffmann


//********* Histograms ********//


   //prepare the canvas for drawing
   TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 800, 600);
   myc->SetFillColor(45);
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.67,0.49,0.99);
   TPad *c1_2 = new TPad("c1_2", "c1_2",0.51,0.67,0.99,0.99);
   TPad *c1_3 = new TPad("c1_3", "c1_3",0.01,0.34,0.49,0.65);
   TPad *c1_4 = new TPad("c1_4", "c1_4",0.51,0.34,0.99,0.65);
   TPad *c1_5 = new TPad("c1_5", "c1_5",0.01,0.01,0.49,0.32);
   TPad *c1_6 = new TPad("c1_6", "c1_6",0.51,0.01,0.99,0.32);
   c1_1->Draw();
   c1_2->Draw();
   c1_3->Draw();
   c1_4->Draw();
   c1_5->Draw();
   c1_6->Draw();
   c1_1->SetFillColor(30);
   c1_1->SetFrameFillColor(42);
   c1_2->SetFillColor(30);
   c1_2->SetFrameFillColor(42);
   c1_3->SetFillColor(30);
   c1_3->SetFrameFillColor(42);
   c1_4->SetFillColor(30);
   c1_4->SetFrameFillColor(42);
   c1_5->SetFillColor(30);
   c1_5->SetFrameFillColor(42);
   c1_6->SetFillColor(30);
   c1_6->SetFrameFillColor(42);

   c1_1->cd();
   TH1::AddDirectory(kFALSE);

   Double_t pi = 3.141592653;
   Double_t ymax = 100.;
   Double_t ymin = 0.;
   Double_t range = ymax-ymin;

   //A function to sample
   TF1 *fsin = new TF1("fsin", "(x<2.&&x>0)*6 + \ (x>2.)*0", ymin, ymax);

   //TF1 *fsin = new TF1("fsin", "cos(2*pi*2*x)", ymin, ymax);
   //TF1 *fsin = new TF1("fsin", "(x<2)*6 + \ (x>2)*0", 0, 4);
   fsin->Draw();

   Int_t n=1000;
   TH1 *hsin = new TH1D("hsin", "hsin", n+1, ymin, ymax);
   Double_t x;

   //Fill the histogram with function values
   for (Int_t i=0; i<=n; i++)
     {
       x = ymin + (Double_t(i)/n)*range;
       hsin->SetBinContent(i+1, fsin->Eval(x));
       //cout<<"hsin["<<i<<"] = "<<hsin[i]<<endl;
       /*
	 Int_t Entries = hsin->GetEntries();
	 cout<<"Entries = "<<Entries<<endl;
	 Double_t Value = hsin->GetBinContent(i+1);
	 cout<<"Value = "<<Value<<endl;
       */
     }
   hsin->Draw("same");
   fsin->GetXaxis()->SetLabelSize(0.05);
   fsin->GetYaxis()->SetLabelSize(0.05);

   c1_2->cd();
   //Compute the transform and look at the magnitude of the output
   TH1 *hm =0;
   //TH1 *hm = new TH1D("hm", "hm", n+1, 0, 4);
   TVirtualFFT::SetTransform(0);
   hm = hsin->FFT(hm, "mag ex");
   
   
   //Double_t *bins = new Double_t[25];
   for (Int_t i=0;i<n+1;i++)
     {

       //Int_t Entries = hm->GetEntries();
       //cout<<"Entries = "<<Entries<<endl;
       Double_t Value = hm->GetBinContent(i+1);
       //cout<<"Bin # = "<<(i+1)<<"   Value = "<<Value<<endl;

       //bins[i] = hm[i];
       //cout<<"bins["<<i<<"] = "<<bins[i]<<endl;
       //cout<<"hm["<<i<<"] = "<<endl;//hm[i]<<endl;
     }

   //Double_t *bins1 = hm->GetArray();

   hm->SetTitle("Magnitude of the 1st transform");
   hm->Draw();
   //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
   //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!

   hm->SetStats(kFALSE);
   hm->GetXaxis()->SetLabelSize(0.05);
   hm->GetYaxis()->SetLabelSize(0.05);
   //hm->GetXaxis()->SetRange(0,4);
   c1_3->cd();
   //Look at the phase of the output
   TH1 *hp = 0;
   hp = hsin->FFT(hp, "PH");
   hp->SetTitle("Phase of the 1st transform");
   hp->Draw();
   hp->SetStats(kFALSE);
   hp->GetXaxis()->SetLabelSize(0.05);
   hp->GetYaxis()->SetLabelSize(0.05);

   //Look at the DC component and the Nyquist harmonic:
   Double_t re, im;
   //That's the way to get the current transform object:
   TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
   c1_4->cd();
   //Use the following method to get just one point of the output
   fft->GetPointComplex(0, re, im);
   printf("1st transform: DC component: %f\n", re);
   fft->GetPointComplex(n/2+1, re, im);
   printf("1st transform: Nyquist harmonic: %f\n", re);

   //Use the following method to get the full output:
   Double_t *re_full = new Double_t[n];
   Double_t *im_full = new Double_t[n];
   fft->GetPointsComplex(re_full,im_full);

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

//********* Data array - same transform ********//

   //Allocate an array big enough to hold the transform output
   //Transform output in 1d contains, for a transform of size N,
   //N/2+1 complex numbers, i.e. 2*(N/2+1) real numbers
   //our transform is of size n+1, because the histogram has n+1 bins

   Double_t *in = new Double_t[2*((n+1)/2+1)];
   Double_t re_2,im_2;
   for (Int_t i=0; i<=n; i++)
     {
       x = (Double_t(i)/n)*(4*TMath::Pi());
       in[i] =  fsin->Eval(x);
     }

   //Make our own TVirtualFFT object (using option "K")
   //Third parameter (option) consists of 3 parts:
   //-transform type:
   // real input/complex output in our case
   //-transform flag:
   // the amount of time spent in planning
   // the transform (see TVirtualFFT class description)
   //-to create a new TVirtualFFT object (option "K") or use the global (default)
   Int_t n_size = n+1;
   TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
   if (!fft_own) return;
   fft_own->SetPoints(in);
   fft_own->Transform();

   //Copy all the output points:
   fft_own->GetPoints(in);
   //Draw the real part of the output
   c1_5->cd();
   TH1 *hr = 0;
   hr = TH1::TransformHisto(fft_own, hr, "RE");
   hr->SetTitle("Real part of the 3rd (array) tranfsorm");
   hr->Draw();
   hr->SetStats(kFALSE);
   hr->GetXaxis()->SetLabelSize(0.05);
   hr->GetYaxis()->SetLabelSize(0.05);
   c1_6->cd();
   TH1 *him = 0;
   him = TH1::TransformHisto(fft_own, him, "IM");
   him->SetTitle("Im. part of the 3rd (array) transform");
   him->Draw();
   him->SetStats(kFALSE);
   him->GetXaxis()->SetLabelSize(0.05);
   him->GetYaxis()->SetLabelSize(0.05);

   myc->cd();
   //Now let's make another transform of the same size
   //The same transform object can be used, as the size and the type of the transform
   //haven't changed
   TF1 *fcos = new TF1("fcos", "cos(x)+cos(0.5*x)+cos(2*x)+1", 0, 4*TMath::Pi());
   for (Int_t i=0; i<=n; i++)
     {
       x = (Double_t(i)/n)*(4*TMath::Pi());
       in[i] =  fcos->Eval(x);
     }
   fft_own->SetPoints(in);
   fft_own->Transform();
   fft_own->GetPointComplex(0, re_2, im_2);
   printf("2nd transform: DC component: %f\n", re_2);
   fft_own->GetPointComplex(n/2+1, re_2, im_2);
   printf("2nd transform: Nyquist harmonic: %f\n", re_2);
   delete fft_own;
   delete [] in;
   delete [] re_full;
   delete [] im_full;



  //Make a new canvas to plot data.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  //gROOT->Reset();    //THIS RESETS GLOBAL VARIABLES!!!

 //TH1 *FFT = new TH1D("FFT", "FFT", n+1, -n/(1./(3./n)), n/(1./(3./n)));
  //TH1 *FFT = new TH1D("FFT", "FFT", (n)/2., -(n+1)/(range*2.), (n+1)/(range*2.));
  TH1 *FFT = new TH1D("FFT", "FFT", (n)/2., 0., (n+1)/(ymax*2.));

  Int_t inflection = (n)/2.;
  cout<<"Inflection = "<<inflection<<endl;

  /*
  //Negative Fourier frequencies.
  for(Int_t i=0;i<n+2;i++)
    {
      if(i>inflection+1)
	{
	  FFT->SetBinContent(i-1-inflection,hm->GetBinContent(i)/(1./(range/(n+1.))));
	}
    }
  */

  //Positive Fourier frequencies.
  for(Int_t i=0;i<n/2.+1..;i++)
    {
      FFT->SetBinContent(i,hm->GetBinContent(i+1)/(1./(range/(n+1.))));
    }
  
  FFT->Draw();
  
  for(Int_t i=0;i<n/2.+1;i++)
    {
      //cout<<"FFT["<<i<<"]   Value = "<<FFT->GetBinContent(i)<<endl;
    }
  

  //Generate random R values from 0-10.
  if(userand == 1)
    {
      /*
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand = new TF1("rand","x",0.,10.);
      for(Int_t i=0;i<12;i++)
	{
	  R[i]=rand->GetRandom();
	  cout<<"R["<<i<<"] ="<<R[i]<<endl;
	}
      */

      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand = new TF1("rand","x",0.,1.);
      R[0] = rand->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand1 = new TF1("rand1","x",R[0],R[0]+1.);
      R[1] = rand1->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand2 = new TF1("rand2","x",R[1],R[1]+1.);
      R[2] = rand2->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand3 = new TF1("rand3","x",R[2],R[2]+1.);
      R[3] = rand3->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand4 = new TF1("rand4","x",R[3],R[3]+1.);
      R[4] = rand4->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand5 = new TF1("rand5","x",R[4],R[4]+1.);
      R[5] = rand5->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand6 = new TF1("rand6","x",R[5],R[5]+1.);
      R[6] = rand6->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand7 = new TF1("rand7","x",R[6],R[6]+1.);
      R[7] = rand7->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand8 = new TF1("rand8","x",R[7],R[7]+1.);
      R[8] = rand8->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand9 = new TF1("rand9","x",R[8],R[8]+1.);
      R[9] = rand9->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand10 = new TF1("rand10","x",R[9],R[9]+1.);
      R[10] = rand10->GetRandom();
      gRandom->SetSeed(0);                    //Sets new random seed.
      TF1 *rand11 = new TF1("rand11","x",R[10],R[10]+1.);
      R[11] = rand11->GetRandom();
    }


  //Fit the FFT histogram ( later this will be real cross section data).
  //Define function fitf for the SOG fit with a variable number of parameters npar.
  Double_t fitf(Double_t *x,Double_t *par)
  {

    /*
    Double_t fitval = exp(-0.25*pow(x[0],2.)*ga2) * ( par[0]/(1.+2*pow(par[12],2.)/ga2) * (cos(x[0]*par[12]) + (2*pow(par[12],2.)/ga2) * (sin(x[0]*par[12])/(x[0]*par[12]))) + par[1]/(1.+2*pow(par[13],2.)/ga2) * (cos(x[0]*par[13]) + (2*pow(par[13],2.)/ga2) * (sin(x[0]*par[13])/(x[0]*par[13]))) + par[2]/(1.+2*pow(par[14],2.)/ga2) * (cos(x[0]*par[14]) + (2*pow(par[14],2.)/ga2) * (sin(x[0]*par[14])/(x[0]*par[14]))) + par[3]/(1.+2*pow(par[15],2.)/ga2) * (cos(x[0]*par[15]) + (2*pow(par[15],2.)/ga2) * (sin(x[0]*par[15])/(x[0]*par[15]))) + par[4]/(1.+2*pow(par[16],2.)/ga2) * (cos(x[0]*par[16]) + (2*pow(par[16],2.)/ga2) * (sin(x[0]*par[16])/(x[0]*par[16]))) + par[5]/(1.+2*pow(par[17],2.)/ga2) * (cos(x[0]*par[17]) + (2*pow(par[17],2.)/ga2) * (sin(x[0]*par[17])/(x[0]*par[17]))) + par[6]/(1.+2*pow(par[18],2.)/ga2) * (cos(x[0]*par[18]) + (2*pow(par[18],2.)/ga2) * (sin(x[0]*par[18])/(x[0]*par[18]))) + par[7]/(1.+2*pow(par[19],2.)/ga2) * (cos(x[0]*par[19]) + (2*pow(par[19],2.)/ga2) * (sin(x[0]*par[19])/(x[0]*par[19]))) + par[8]/(1.+2*pow(par[20],2.)/ga2) * (cos(x[0]*par[20]) + (2*pow(par[20],2.)/ga2) * (sin(x[0]*par[20])/(x[0]*par[20]))) + par[9]/(1.+2*pow(par[21],2.)/ga2) * (cos(x[0]*par[21]) + (2*pow(par[21],2.)/ga2) * (sin(x[0]*par[21])/(x[0]*par[21]))) + par[10]/(1.+2*pow(par[22],2.)/ga2) * (cos(x[0]*par[22]) + (2*pow(par[22],2.)/ga2) * (sin(x[0]*par[22])/(x[0]*par[22]))) + par[11]/(1.+2*pow(par[23],2.)/ga2) * (cos(x[0]*par[23]) + (2*pow(par[23],2.)/ga2) * (sin(x[0]*par[23])/(x[0]*par[23]))) );
    */

    
    Double_t fitval = 0.;
    Double_t tempfitval = 0.;
    for(Int_t i=0;i<npar;i++)
      {
	tempfitval = ( par[i]/(1.+2*pow(R[i],2.)/ga2) * (cos(x[0]*R[i]) + (2*pow(R[i],2.)/ga2) * (sin(x[0]*R[i])/(x[0]*R[i]))) );
	fitval = fitval + tempfitval;
      }

    //Calibrate evenly spaced step sizes. //No fit for some reason.
    /*for(Int_t i=0;i<npar;i++)
      {
	tempfitval = ( par[i]/(1.+2*pow(R[i]+(i+1)*par[npar],2.)/ga2) * (cos(x[0]*(R[i]+(i+1)*par[npar])) + (2*pow(R[i]+(i+1)*par[npar],2.)/ga2) * (sin(x[0]*(R[i]+(i+1)*par[npar]))/(x[0]*(R[i]+(i+1)*par[npar])))) );
	fitval = fitval + tempfitval;
	}*/

    //Calibrate evenly spaced step sizes with Rmax being calibrated.
    /*for(Int_t i=0;i<npar;i++)
      {
	tempfitval = ( par[i]/(1.+2*pow(par[npar]/(npar-i),2.)/ga2) * (cos(x[0]*par[npar]/(npar-i)) + (2*pow(par[npar]/(npar-i),2.)/ga2) * (sin(x[0]*par[npar]/(npar-i))/(x[0]*par[npar]/(npar-i)))) );
	fitval = fitval + tempfitval;
	}*/

    //Calibrate using m as a parameter with Amroun's R[i]s.
    /*for(Int_t i=0;i<npar;i++)
      {
	tempfitval = ( par[i]/(1.+2*pow(R[i]*par[npar],2.)/ga2) * (cos(x[0]*R[i]*par[npar]) + (2*pow(R[i]*par[npar],2.)/ga2) * (sin(x[0]*R[i]*par[npar])/(x[0]*R[i]*par[npar]))) );
	fitval = fitval + tempfitval;
	}
    */
    fitval = exp(-0.25*pow(x[0],2.)*ga2) * fitval;


    //Double_t fitval = exp(-0.25*pow(x[0],2.)*ga2) * ( par[0]/(1.+2*pow(R[0],2.)/ga2) * (cos(x[0]*R[0]) + (2*pow(R[0],2.)/ga2) * (sin(x[0]*R[0])/(x[0]*R[0]))) + par[1]/(1.+2*pow(R[1],2.)/ga2) * (cos(x[0]*R[1]) + (2*pow(R[1],2.)/ga2) * (sin(x[0]*R[1])/(x[0]*R[1]))) + par[2]/(1.+2*pow(R[2],2.)/ga2) * (cos(x[0]*R[2]) + (2*pow(R[2],2.)/ga2) * (sin(x[0]*R[2])/(x[0]*R[2]))) + par[3]/(1.+2*pow(R[3],2.)/ga2) * (cos(x[0]*R[3]) + (2*pow(R[3],2.)/ga2) * (sin(x[0]*R[3])/(x[0]*R[3]))) + par[4]/(1.+2*pow(R[4],2.)/ga2) * (cos(x[0]*R[4]) + (2*pow(R[4],2.)/ga2) * (sin(x[0]*R[4])/(x[0]*R[4]))) + par[5]/(1.+2*pow(R[5],2.)/ga2) * (cos(x[0]*R[5]) + (2*pow(R[5],2.)/ga2) * (sin(x[0]*R[5])/(x[0]*R[5]))) + par[6]/(1.+2*pow(R[6],2.)/ga2) * (cos(x[0]*R[6]) + (2*pow(R[6],2.)/ga2) * (sin(x[0]*R[6])/(x[0]*R[6]))) + par[7]/(1.+2*pow(R[7],2.)/ga2) * (cos(x[0]*R[7]) + (2*pow(R[7],2.)/ga2) * (sin(x[0]*R[7])/(x[0]*R[7]))) + par[8]/(1.+2*pow(R[8],2.)/ga2) * (cos(x[0]*R[8]) + (2*pow(R[8],2.)/ga2) * (sin(x[0]*R[8])/(x[0]*R[8]))) + par[9]/(1.+2*pow(R[9],2.)/ga2) * (cos(x[0]*R[9]) + (2*pow(R[9],2.)/ga2) * (sin(x[0]*R[9])/(x[0]*R[9]))) + par[10]/(1.+2*pow(R[10],2.)/ga2) * (cos(x[0]*R[10]) + (2*pow(R[10],2.)/ga2) * (sin(x[0]*R[10])/(x[0]*R[10]))) + par[11]/(1.+2*pow(R[11],2.)/ga2) * (cos(x[0]*R[11]) + (2*pow(R[11],2.)/ga2) * (sin(x[0]*R[11])/(x[0]*R[11]))) );
     

    return fitval;
  }
				 
  TF1 *f1 = new TF1("fit",fitf,fitmin,fitmax,npar);
  //fit->SetParameter(0,1.);
  //fit->SetParameters(1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.);    //12+ set parameters -> Error. Not sure why.

  /*
  Double_t m = 6.5;
  Double_t R[12] = {0.1*m,.5*m,.9*m,1.3*m,1.6*m,2.*m,2.4*m,2.9*m,3.4*m,4.*m,4.6*m,5.2*m};
  fit->SetParameter(12,R[0]);
  fit->SetParameter(13,R[1]);
  fit->SetParameter(14,R[2]);
  fit->SetParameter(15,R[3]);
  fit->SetParameter(16,R[4]);
  fit->SetParameter(17,R[5]);
  fit->SetParameter(18,R[6]);
  fit->SetParameter(19,R[7]);
  fit->SetParameter(20,R[8]);
  fit->SetParameter(21,R[9]);
  fit->SetParameter(22,R[10]);
  fit->SetParameter(23,R[11]);
  */
  //Set initial spacing between R[i] values.
  //fit->SetParameter(12,15.);

  //Make sure all parameter's are positive as per I. Sick pg 513.
  /*
  fit->SetParLimits(0,0,10000);
  fit->SetParLimits(1,0,10000);
  fit->SetParLimits(2,0,10000);
  fit->SetParLimits(3,0,10000);
  fit->SetParLimits(4,0,10000);
  fit->SetParLimits(5,0,10000);
  fit->SetParLimits(6,0,10000);
  fit->SetParLimits(7,0,10000);
  fit->SetParLimits(8,0,10000);
  fit->SetParLimits(9,0,10000);
  fit->SetParLimits(10,0,10000);
  fit->SetParLimits(11,0,10000);
  */
  FFT->Fit("fit","R");
  //f1->Draw("same");

  Double_t chi2 = fit->GetChisquare();
  cout<<"chi^2 = "<<chi2<<endl;

  //Plot the individual gaussians making up the total fit. 
  Double_t Q0 = fit->GetParameter(0);
  Double_t Q1 = fit->GetParameter(1);
  Double_t Q2 = fit->GetParameter(2);
  Double_t Q3 = fit->GetParameter(3);
  Double_t Q4 = fit->GetParameter(4);
  Double_t Q5 = fit->GetParameter(5);
  Double_t Q6 = fit->GetParameter(6);
  Double_t Q7 = fit->GetParameter(7);
  Double_t Q8 = fit->GetParameter(8);
  Double_t Q9 = fit->GetParameter(9);
  Double_t Q10 = fit->GetParameter(10);
  Double_t Q11 = fit->GetParameter(11);

  //Function for individual Gaussians.
  Double_t fitg(Double_t *x,Double_t *par)
  {
    Double_t ga2 = pow(0.8*pow(2./3., 0.5), 2.);
    Double_t fitval =exp(-0.25*pow(x[0],2.)*ga2) * ( par[0]/(1.+2*pow(par[1],2.)/ga2) * (cos(x[0]*par[1]) + (2*pow(par[1],2.)/ga2) * (sin(x[0]*par[1])/(x[0]*par[1]))) );
    return fitval;
  }

  //Plot individual Gaussians with their fit parameters. 
  TF1 *g0 = new TF1("g0", fitg, ymin, ymax,2.);
  g0->SetParameters(Q0,R[0]);
  g0->SetLineColor(1);
  g0->Draw("same");
  TF1 *g1 = new TF1("g1", fitg, ymin, ymax,2.);
  g1->SetParameters(Q1,R[1]);
  g1->SetLineColor(2);
  g1->Draw("Same");
  TF1 *g2 = new TF1("g2", fitg, ymin, ymax,2.);
  g2->SetParameters(Q2,R[2]);
  g2->SetLineColor(3);
  g2->Draw("Same");
  TF1 *g3 = new TF1("g3", fitg, ymin, ymax,2.);
  g3->SetParameters(Q3,R[3]);
  g3->SetLineColor(4);
  g3->Draw("Same");
  TF1 *g4 = new TF1("g4", fitg, ymin, ymax,2.);
  g4->SetParameters(Q4,R[4]);
  g4->SetLineColor(5);
  g4->Draw("Same");
  TF1 *g5 = new TF1("g5", fitg, ymin, ymax,2.);
  g5->SetParameters(Q5,R[5]);
  g5->SetLineColor(6);
  g5->Draw("Same");
  TF1 *g6 = new TF1("g6", fitg, ymin, ymax,2.);
  g6->SetParameters(Q6,R[6]);
  g6->SetLineColor(7);
  g6->Draw("Same");
  TF1 *g7 = new TF1("g7", fitg, ymin, ymax,2.);
  g7->SetParameters(Q7,R[7]);
  g7->SetLineColor(8);
  g7->Draw("Same");
  TF1 *g8 = new TF1("g8", fitg, ymin, ymax,2.);
  g8->SetParameters(Q8,R[8]);
  g8->SetLineColor(9);
  g8->Draw("Same");
  TF1 *g9 = new TF1("g9", fitg, ymin, ymax,2.);
  g9->SetParameters(Q9,R[9]);
  g9->SetLineColor(10);
  g9->Draw("Same");
  TF1 *g10 = new TF1("g10", fitg, ymin, ymax,2.);
  g10->SetParameters(Q10,R[10]);
  g10->SetLineColor(11);
  g10->Draw("Same");
  TF1 *g11 = new TF1("g11", fitg, ymin, ymax,2.);
  g11->SetParameters(Q11,R[11]);
  g11->SetLineColor(12);
  g11->Draw("Same");

  for(Int_t i=0;i<12;i++)
    {
      cout<<"R["<<i<<"] = "<<R[i]<<endl;
    }


  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();

  //Fill a new histogram with data using the fit function. This will then be inverse Fourier transformed back to be suce the original function is recoverable. 
  TH1 *hfit = new TH1D("hfit", "hfit", n*1.+1, ymin, ymax);
  //Fill the histogram with function values
  for (Int_t i=0; i<=n*1.; i++)
    {
      x = ymin + (Double_t(i)/(n*1.))*fitrange;//range;  //Only fill bins up to max of fitted range.
      //cout<<x<<endl;
      if(x==0.)
	{
	  x = 0.0000001; //Inf/NaN at zero for some resaon.
	}
      hfit->SetBinContent(i+1, fit->Eval(x));
      //hfit->SetBinContent(i+1, cos(500.*x)); //Test function.
      hfit->GetEntries();
      //cout<<hfit->GetBinContent(i+1)<<endl;
      //cout<<fit->Eval(x)<<endl;

    }
  //hfit->SetFillColor(17);
  hfit->Draw("");

  //cout<<"fit(1.5) = "<<fit->Eval(1.5)<<endl;
  //Do an inverse FFT on hfit to get back to the charge density.





  //Make a new canvas to plot data.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();

  /*
  //Look at the DC component and the Nyquist harmonic:
  Double_t re, im;
  //That's the way to get the current transform object:
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  c1_4->cd();
  //Use the following method to get just one point of the output
  fft->GetPointComplex(0, re, im);
  printf("1st transform: DC component: %f\n", re);
  fft->GetPointComplex(n/2+1, re, im);
  printf("1st transform: Nyquist harmonic: %f\n", re);
  */


  
  //Create arrays for real and complex imnputs for inverse FFT. 
  Double_t *re_full = new Double_t[n*1.];
  Double_t *im_full = new Double_t[n*1.];
  
  //Fill the real and complex arrays. The complex array is all zeros since we have real data. The real data is from the histo fit. 
  for(Int_t i=0;i<n*1.;i++)
    {
      re_full[i] = hfit->GetBinContent(i+1);
      im_full[i] = 0;
    }
  
  TVirtualFFT *iFFT = TVirtualFFT::FFT(1, &n, "C2R M K");
  iFFT->SetPointsComplex(re_full,im_full);
  iFFT->Transform();
  TH1 *hcharge = 0;
  //Let's look at the output
  hcharge = TH1::TransformHisto(iFFT,hcharge,"re");     //Since there are no complex numbers here re or mag are equivalent.
  hcharge->SetTitle("The backward transform result");
  hcharge->Draw();
  //NOTE: here you get at the x-axes number of bins and not real values
  //(in this case 25 bins has to be rescaled to a range between 0 and 4*Pi;
  //also here the y-axes has to be rescaled (factor 1/bins)
  hcharge->SetStats(kFALSE);
  hcharge->GetXaxis()->SetLabelSize(0.05);
  hcharge->GetYaxis()->SetLabelSize(0.05);
  delete fft_back;
  fft_back=0;

  //Make a new canvas to plot data.
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();

  TH1 *FFTBack = new TH1D("FFTBack", "FFTBack", n*1.+1, -(n*1.+1)/(ymax*2.), (n*1.+1)/(ymax*2.));   

  //Int_t inflection = (n)/2.;
  //cout<<"Inflection = "<<inflection<<endl;

  
  //Negative Fourier frequencies.
  for(Int_t i=0;i<n*1.+2;i++)
    {
      if(i>inflection+1)
	{
	  FFTBack->SetBinContent(i-1-inflection,hcharge->GetBinContent(i)/(1./(fitrange/(n+1.))));//range/(n+1.))));
	}
    }
  

  //Positive Fourier frequencies.
  for(Int_t i=0;i<n*1.;i++)
    {
      if(i<=inflection)
	{
	  FFTBack->SetBinContent(i+inflection+1,hcharge->GetBinContent(i+1)/(1./(fitrange/(n+1.))));
	}
    }
  
  FFTBack->Draw();
  
  for(Int_t i=0;i<n/2.+1;i++)
    {
      //cout<<"FFT["<<i<<"]   Value = "<<FFT->GetBinContent(i)<<endl;
    }

  


  /*
  TH1 *hcharge =0;
  TVirtualFFT::SetTransform(0);
  hcharge = hfit->FFT(hcharge, "re ex");
  hcharge->Draw("same");
  */
}
