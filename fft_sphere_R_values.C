#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

//gRoot->Reset();
Double_t ga2 = pow(0.8*pow(2./3., 0.5), 2.);
Double_t m = 6.;
Double_t R[12] = {0.1*m,.5*m,.9*m,1.3*m,1.6*m,2.*m,2.4*m,2.9*m,3.4*m,4.*m,4.6*m,5.2*m};
Int_t loop = 0;
Int_t loopmax = 10;
Double_t steploopmax = 1.;
Int_t npar = 12;
Double_t  avgchi2 = 0.;
Double_t  avgchi2temp = 0.;

void fft_sphere_R_values()
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

/*
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
*/

   Double_t pi = 3.141592653;
   Double_t ymax = 100.;
   Double_t ymin = 0.;
   Double_t range = ymax-ymin;

   //A function to sample
   TF1 *fsin = new TF1("fsin", "(x<2.0&&x>0)*6 + \ (x>2.)*0", ymin, ymax);

   //TF1 *fsin = new TF1("fsin", "cos(2*pi*2*x)", ymin, ymax);
   //TF1 *fsin = new TF1("fsin", "(x<2)*6 + \ (x>2)*0", 0, 4);
   //fsin->Draw();

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
   //hsin->Draw("same");
   //fsin->GetXaxis()->SetLabelSize(0.05);
   //fsin->GetYaxis()->SetLabelSize(0.05);

   //c1_2->cd();
   //Compute the transform and look at the magnitude of the output
   TH1 *hm =0;
   //TH1 *hm = new TH1D("hm", "hm", n+1, 0, 4);
   TVirtualFFT::SetTransform(0);
   hm = hsin->FFT(hm, "mag ex");
   
   /*
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
   */
   /*
   //Double_t *bins1 = hm->GetArray();

   //hm->SetTitle("Magnitude of the 1st transform");
   //hm->Draw();
   //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
   //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!

   //hm->SetStats(kFALSE);
   //hm->GetXaxis()->SetLabelSize(0.05);
   //hm->GetYaxis()->SetLabelSize(0.05);
   //hm->GetXaxis()->SetRange(0,4);
   //c1_3->cd();
   //Look at the phase of the output
   TH1 *hp = 0;
   hp = hsin->FFT(hp, "PH");
   //hp->SetTitle("Phase of the 1st transform");
   //hp->Draw();
   //hp->SetStats(kFALSE);
   //hp->GetXaxis()->SetLabelSize(0.05);
   //hp->GetYaxis()->SetLabelSize(0.05);

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

//********* Data array - same transform ********

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
*/


  //Make a new canvas to plot data.
  //TCanvas* c1=new TCanvas("c1");
  //c1->SetGrid();
  //gROOT->Reset();    //THIS RESETS GLOBAL VARIABLES!!!

  //File to write output of good R values to.
  std::ofstream output ("R_Values.txt", std::ofstream::out);
  output<<"rChi^2    R[1]     R[2]    R[3]    R[4]    R[5]    R[6]    R[7]    R[8]    R[9]    R[10]    R[11]   Q0   Q1   Q2   Q3   Q4   Q5   Q6   Q7   Q8   Q9   Q10   Q11   "<<endl;

 //TH1 *FFT = new TH1D("FFT", "FFT", n+1, -n/(1./(3./n)), n/(1./(3./n)));
  //TH1 *FFT = new TH1D("FFT", "FFT", (n)/2., -(n+1)/(range*2.), (n+1)/(range*2.));
  TH1 *FFT = new TH1D("FFT", "FFT", (n)/2., 0., (n+1)/(ymax*2.));

  Int_t inflection = (n)/2.;
  //cout<<"Inflection = "<<inflection<<endl;

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
  
  //FFT->Draw();
  
  /*
  for(Int_t i=0;i<n/2.+1;i++)
    {
      //cout<<"FFT["<<i<<"]   Value = "<<FFT->GetBinContent(i)<<endl;
    }
  */

  for(Double_t j=0.;j<steploopmax;j++)
    //while(steploop<steploopmax)
    {
      loop = 0;         //Reset inner while loop's variable.
      avgchi2 = 0.;
      avgchi2temp = 0.;
      
      //Loop through various R[i] values while calculating and outputting the 'good' chi^2 values.
      while(loop<loopmax)
	{
	  //Generate random R values from 0-10.
	  /*gRandom->SetSeed(0);                    //Sets new random seed.
	    TF1 *rand = new TF1("rand","x",0.,10.);
	    for(Int_t i=0;i<12;i++)
	    {
	    R[i]=rand->GetRandom();
	    cout<<"R["<<i<<"] ="<<R[i]<<endl;
	    }*/
	  
	  Double_t step = 3.06+0.22*(j/steploopmax); //Max distance of R[i+1] from R[i].
	  Double_t d = 3.06+0.22*(j/steploopmax)-0.1;    //Min distance of R[i+1] from R[i].
	  //Generate random but increasing R[i] values. 
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",d,step);
	  R[0] = rand->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",R[0]+d,R[0]+step);
	  R[1] = rand1->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",R[1]+d,R[1]+step);
	  R[2] = rand2->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",R[2]+d,R[2]+step);
	  R[3] = rand3->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",R[3]+d,R[3]+step);
	  R[4] = rand4->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",R[4]+d,R[4]+step);
	  R[5] = rand5->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",R[5]+d,R[5]+step);
	  R[6] = rand6->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",R[6]+d,R[6]+step);
	  R[7] = rand7->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",R[7]+d,R[7]+step);
	  R[8] = rand8->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand9 = new TF1("rand9","x",R[8]+d,R[8]+step);
	  R[9] = rand9->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand10 = new TF1("rand10","x",R[9]+d,R[9]+step);
	  R[10] = rand10->GetRandom();
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand11 = new TF1("rand11","x",R[10]+d,R[10]+step);
	  R[11] = rand11->GetRandom();
	  
	  //Fit the FFT histogram (later this will be real cross section data).
	  //Define function fitf for the SOG fit with a variable number of parameters npar.
	  Double_t fitf(Double_t *x,Double_t *par)
	  {   
	    Double_t fitval = 0.;
	    Double_t tempfitval = 0.;
	    for(Int_t i=0;i<npar;i++)
	      {
		tempfitval = ( par[i]/(1.+2*pow(R[i],2.)/ga2) * (cos(x[0]*R[i]) + (2*pow(R[i],2.)/ga2) * (sin(x[0]*R[i])/(x[0]*R[i]))) );
		fitval = fitval + tempfitval;
	      }

	    fitval = exp(-0.25*pow(x[0],2.)*ga2) * fitval;
	    
	    return fitval;
	  }
	  
	  TF1 *f1 = new TF1("fit",fitf,0.,1.,12.0);
	  //fit->SetParameter(0,1.);
	  fit->SetParameters(1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.);    //12+ set parameters -> Error. Not sure why.
	  
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
	  
	  FFT->Fit("fit","R0q");   //R fits a certain range only. 0 supresses drawing the fit on the histogram.
	  //f1->Draw("same");
	  
	  Double_t chi2 = fit->GetChisquare();
	  cout<<"chi^2 = "<<chi2<<endl;
	  avgchi2temp = avgchi2temp + chi2;
	  
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
	  
	  if(fit->GetChisquare() < 0.9)
	    {
	      //output<<chi2<<" "<<Q0<<" "<<Q1<<" "<<Q2<<" "<<Q3<<" "<<Q4<<" "<<Q5<<" "<<Q6<<" "<<Q7<<" "<<Q8<<" "<<Q9<<" "<<Q10<<" "<<Q11<<endl;
	      output<<chi2<<" "<<R[0]<<" "<<R[1]<<" "<<R[2]<<" "<<R[3]<<" "<<R[4]<<" "<<R[5]<<" "<<R[6]<<" "<<R[7]<<" "<<R[8]<<" "<<R[9]<<" "<<R[10]<<" "<<R[11]<<" "<<Q0<<" "<<Q1<<" "<<Q2<<" "<<Q3<<" "<<Q4<<" "<<Q5<<" "<<Q6<<" "<<Q7<<" "<<Q8<<" "<<Q9<<" "<<Q10<<" "<<Q11<<endl;
	    }
	  
	  //output.close();
	  
	  
	  
	  
	  loop++;
	}//End while loop that generates and checks fits for random R[i] values.
      
      avgchi2 = avgchi2temp/loopmax;
      cout<<"R[i] step size = "<<step<<"   Average Reduced Chi ^2 = "<<avgchi2<<endl;
      output<<"R[i] step size = "<<step<<"   Average Reduced Chi ^2 = "<<avgchi2<<endl;

    }//End outer for loop to increment step sizes of R[i] values. 
  
  output.close();
}
