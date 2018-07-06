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

Double_t pi = 3.141592654;
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1.0/0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 

Double_t gamma = 0.8*pow(2.0/3.0,0.5);                  //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
Double_t ymin = 0.0;
Double_t ymax = 1000.;
Double_t range = fabs(ymax - ymin);
Int_t n = 10000;
Int_t ndim = n+1;
Double_t truncate = 25.;                 //Truncate the histogram before inverse FFT. [fm^-2]

Double_t R[12] = {0.1, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9, 3.4, 4.0, 4.6, 5.2};  //Radii [fm].
Double_t QH3ch[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121, 0.0, 0.0};

void fft_charge_ff_H3() 
{
  //Make a new canvas to plot data.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  
  //Define fit function for H3 charge FF as a function of Q^2.
  Double_t fitch(Double_t *Q2, Double_t *par)
  {
    Double_t fitval = 0.;
    Double_t sumH3chtemp = 0.;
    
    for(Int_t j=0; j<12; j++)
      {
	sumH3chtemp = (QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[j])/(pow(Q2[0],0.5)*R[j])) );
	fitval = fitval + sumH3chtemp;
      }
    
    fitval = fabs( fitval ) * exp(-0.25*Q2[0]*pow(gamma,2.0));
    return fitval;
  }
  
  TF1 *H3chFF = new TF1("H3chFF",fitch, ymin, ymax,1.);
  //H3chFF->SetParLimits(0,0.00001,100.);
  c1->SetLogy();
  H3chFF->SetTitle("H3 Charge Form Factor");
  H3chFF->GetXaxis()->SetTitle("Q^{2} [fm^{-2}]");
  H3chFF->GetYaxis()->SetTitle("Fch(Q)");
  H3chFF->Draw();

//Fill a new histogram with data using the fit function. This will then be inverse Fourier transformed to view the charge distribution.
  TH1 *hfit = new TH1D("hH3chFF", "hH3chFF", n+1, ymin, ymax);
  //Fill the histogram with function values
  for (Int_t i=0; i<=n; i++)
    {
      Double_t x = 0.;
      x = ymin + (Double_t(i)/(n))*(range+0.);//range;  //Only fill bins up to max of fitted range.
      //cout<<x<<endl;
      if(x==0.)
	{
	  x = 0.0000001; //Inf/NaN at zero.
	}
      
      if(x<=ymax && x<truncate)
	{
	  hH3chFF->SetBinContent(i+1, H3chFF->Eval(x));
	}
      
      if(x>ymax)
	{
	  //hfit->SetBinContent(i+1, (0.8*cos(12.*x-3.1)+0.5)*exp(-x*.1));
	  //hH3chFF->SetBinContent(i+1, 0.);
	}
      
      hH3chFF->GetEntries();
    }
  //hfit->SetFillColor(17);
  hH3chFF->Draw("same");
  
  //Inverse Fourier transform the H3 charge FF to get the charge distribution of H3.
  //Create arrays for real and complex imnputs for inverse FFT. 
  Double_t *re_full = new Double_t[n+1];
  Double_t *im_full = new Double_t[n+1];
  
  //Fill the real and complex arrays. The complex array is all zeros since we have real data. The real data is from the histo fit. 
  for(Int_t i=0;i<n+1;i++)
    {
      re_full[i] = hH3chFF->GetBinContent(i+1);
      im_full[i] = 0;
      //cout<<"re_full["<<i<<"] = "<<re_full[i]<<endl;
    }

  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();

  TVirtualFFT *iFFT = TVirtualFFT::FFT(1, &ndim, "C2R M K");   //Re and Mag look the same for C2R which makes sense. C2CBACKWARD mag definitely different from C2R and looks wrong. Same for C2CBackward re although first min x position looks ok. Stick to C2R mag it appears correct. 
  iFFT->SetPointsComplex(re_full,im_full);
  iFFT->Transform();
  TH1 *hcharge = 0;
  //TH1 *hcharge = new TH1D("hcharge", "hcharge", nback*10.+1, ymin, ymax);
  //Let's look at the output
  hcharge = TH1::TransformHisto(iFFT,hcharge,"mag");     //Not totally sure if we want re or mag.
  hcharge->SetTitle("H3 Charge Distribution");
  //hcharge->Draw();
  //NOTE: here you get at the x-axes number of bins and not real values
  //(in this case 25 bins has to be rescaled to a range between 0 and 4*Pi;
  //also here the y-axes has to be rescaled (factor 1/bins)
  hcharge->SetStats(kFALSE);
  hcharge->GetXaxis()->SetLabelSize(0.05);
  hcharge->GetYaxis()->SetLabelSize(0.05);
  delete iFFT;
  iFFT=0;

  //Rebin the inverse FT result to compensate for ROOT's weird output. 
  TH1 *H3ch_dist = new TH1D("H3ch_dist", "H3ch_dist", n+1, -(n+1)/(ymax*2.), (n+1)/(ymax*2.));   

  Int_t inflection = (n)/2.;
  //cout<<"inflection = "<<inflectionback<<endl;
  
  //Negative Fourier frequencies.
  for(Int_t i=0;i<n+2;i++)
    {
      if(i>inflection+1)
	{
	  H3ch_dist->SetBinContent(i-1-inflection,hcharge->GetBinContent(i)/(1./(range/(n+1.))));//range/(n+1.))));
	  //cout<<"i - inflectionback - 1 = "<<i<<" - "<<inflectionback<<"-1 = "<<i-inflectionback-1<<endl;
	}
    }
  
  //Positive Fourier frequencies.
  for(Int_t i=0;i<n;i++)
    {
      if(i<=inflection)
	{
	  H3ch_dist->SetBinContent(i+inflection+1,hcharge->GetBinContent(i+1)/(1./(range/(n+1.))));
	  //cout<<i+inflection+1<<endl;
	}
    }
  c2->SetLogy();
  H3ch_dist->SetTitle("H3 Charge Distribution");
  H3ch_dist->GetXaxis()->SetTitle("r");
  H3ch_dist->GetYaxis()->SetTitle("#rho(r)");
  H3ch_dist->Draw("");


  //Now (forward) Fourier transform the H3 charge distribution to recover the H3 charge FF. 

  //Make a new canvas to plot data.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  
  //Compute the transform and look at the magnitude of the output
  TH1 *FFback =0;
  //TH1 *hm = new TH1D("hm", "hm", n+1, 0, 4);
  TVirtualFFT::SetTransform(0);
  FFback = H3ch_dist->FFT(FFback, "mag ex");
  //FFback->Draw("");

  //Rebin the FT result to compensate for ROOT's weird output. 
  TH1 *hH3chFFback = new TH1D("hH3chFFback", "hH3chFFback", n+1,-range/2.,range/2.);//-(n+1)/(ymax*2.), (n+1)/(ymax*2.)); 

  //Negative Fourier frequencies.
  for(Int_t i=0;i<n+2;i++)
    {
      if(i>inflection+1)
	{
	  hH3chFFback->SetBinContent(i-1-inflection,FFback->GetBinContent(i)*(1./(range)) );//range/(n+1.))));
	  //cout<<"i - inflectionback - 1 = "<<i<<" - "<<inflectionback<<"-1 = "<<i-inflectionback-1<<endl;
	}
    }
  
  //Positive Fourier frequencies.
  for(Int_t i=0;i<n;i++)
    {
      if(i<=inflection)
	{
	  hH3chFFback->SetBinContent(i+inflection+1,FFback->GetBinContent(i+1)*(1./(range)) );
	  //cout<<i+inflection+1<<endl;
	}
    }
  c3->SetLogy();
  hH3chFFback->SetTitle("H3 Charge Form Factor Transformed back (FF->iFFT(FF)->FFT(iFFT(FF))=FF)");
  hH3chFFback->GetXaxis()->SetTitle("Q^{2} [fm^{-2}]");
  hH3chFFback->GetYaxis()->SetTitle("Fch(Q)");
  hH3chFFback->Draw("");
}
