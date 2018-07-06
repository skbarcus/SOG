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
Double_t alpha = 1.0/137.0;              //Fine structure constant.

Int_t npar = 24;                         //Number of parameters in fit.
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
Int_t z = 12;
Double_t MtC = 12.0*0.9315;              //Mass of C12 in GeV.
Double_t gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
Double_t E0 = 0.3745;                    //Initial e- energy GeV.
Double_t Ef = 0.;                        //Final e- energy GeV.
Double_t ymin = 30.;
Double_t ymax = 100.;
Double_t range = fabs(ymax - ymin);
Int_t n = 10000;
Int_t ndim = n+1;
Double_t truncate = 25.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 2;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char* str[100];                          //Variable to read lines of the data file.
Float_t thetatemp,qefftemp,sigexptemp;
Double_t theta[100];                     //Angle in degrees.
Double_t qeff[100];                      //q effective in fm^-1.
Double_t sigexp[100];                    //Sigma experimental (cross section). Not sure on units yet.

Double_t m = 2.5;
Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
Double_t QH3ch[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121, 0.0, 0.0};

void C12_STANF_374_SOG() 
{
  //Make a new canvas to plot data.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  //Read in data from text file.
  //Open file.
  FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/STANF_374.txt","r");

  //Read in data.
  while (1) {
    //Skips the first 5 lines of the file. 
    if (nlines < skip)
      {
	fgets(str,100,fp);
	nlines++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(fp,"%f %f %f",&thetatemp , &qefftemp, &sigexptemp);
	if (ncols < 0) break;   
	cout<<thetatemp<<"   "<<qefftemp<<"   "<<sigexptemp<<endl;
	theta[nlines-skip] = thetatemp;
	qeff[nlines-skip] = qefftemp;
	sigexp[nlines-skip] = sigexptemp;
	//Fill histograms with x and y data.
	//h1->Fill(x);
	//h2->Fill(x,y);
	//h2->SetMarkerSize(5);
	//Fill ntuple with x and y data.
	//ntuple->Fill(x,y);
	//Count the number of lines in the file. 
	nlines++;
      }
  }
  cout<<"Number of lines = "<<nlines<<endl;
  fclose(fp);

 //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph = new TGraph(nlines-skip,theta,sigexp);
  //Draw the new TGraph called graph on the canvas. 
  graph->Draw("");
  c1->SetLogy();
  //Set X axis
  //graph->GetXaxis()->SetLimits(-12,12);
  //Set Y axis Min and Max (not sure why different from X).
  //graph->SetMinimum(0);
  //graph->SetMaximum(120);
  graph->SetLineWidth(1);
  graph->SetLineColor(4);
  graph->SetFillColor(0);
  graph->SetMarkerColor(1);
  graph->SetMarkerSize(1.);
  graph->SetMarkerStyle(20);
  graph->SetTitle("Main Title; X Axis Title; Y Axis Title");
  //graph_expected.SetFillColor(kYellow);
  //graph_expected.DrawClone("E3AL"); // E3 draws the band

  // Draw the Legend
  TLegend leg(0.9,.7,.56,.9,"Legend Title");
  leg.SetFillColor(0);
  leg.AddEntry(graph,"Curve Name");
  leg.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.



  
  //Define fit function for H3 charge FF as a function of Q^2.
  Double_t fitch(Double_t *angle, Double_t *par)
  {
    Double_t fitval = 0.;
    Double_t sumH3chtemp = 0.;
    Ef = E0/(1.0+2.0*E0*pow(sin(angle[0]*deg2rad/2.0),2.0)/MtC);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle[0]*deg2rad/2.0),2.0) * GeV2fm;
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*6.*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(12.0,1.0/3.0))) ,2.0);   //A=6 Z=12
    
    for(Int_t i=0; i<ngaus; i++)
      {
	//sumH3chtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Fit R[i] values and Qi.
	sumH3chtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );

	fitval = fitval + sumH3chtemp;
      }
    
    fitval = fabs( fitval ) * exp(-0.25*Q2eff*pow(gamma,2.0));
    return fitval;
  }
  
  
  TF1 *H3chFF = new TF1("H3chFF",fitch, ymin, ymax,npar);
  //H3chFF->SetParLimits(0,0.00001,100.);
  //c1->SetLogy();
  
  /*
  H3chFF->SetParameter(0,1.);
  H3chFF->SetParameter(1,1.);
  H3chFF->SetParameter(2,1.);
  H3chFF->SetParameter(3,1.);
  H3chFF->SetParameter(4,1.);
  H3chFF->SetParameter(5,1.);
  H3chFF->SetParameter(6,1.);
  H3chFF->SetParameter(7,1.);
  H3chFF->SetParameter(8,1.);
  H3chFF->SetParameter(9,1.);
  H3chFF->SetParameter(10,1.);
  H3chFF->SetParameter(11,1.);
*/
  
  H3chFF->SetParameter(12,R[0]);
  H3chFF->SetParameter(13,R[1]);
  H3chFF->SetParameter(14,R[2]);
  H3chFF->SetParameter(15,R[3]);
  H3chFF->SetParameter(16,R[4]);
  H3chFF->SetParameter(17,R[5]);
  H3chFF->SetParameter(18,R[6]);
  H3chFF->SetParameter(19,R[7]);
  H3chFF->SetParameter(20,R[8]);
  H3chFF->SetParameter(21,R[9]);
  H3chFF->SetParameter(22,R[10]);
  H3chFF->SetParameter(23,R[11]);
  
  /*
  H3chFF->SetParLimits(0,0.000001,10000);
  H3chFF->SetParLimits(1,0.000001,10000);
  H3chFF->SetParLimits(2,0.000001,10000);
  H3chFF->SetParLimits(3,0.000001,10000);
  H3chFF->SetParLimits(4,0.000001,10000);
  H3chFF->SetParLimits(5,0.000001,10000);
  H3chFF->SetParLimits(6,0.000001,10000);
  H3chFF->SetParLimits(7,0.000001,10000);
  H3chFF->SetParLimits(8,0.000001,10000);
  H3chFF->SetParLimits(9,0.000001,10000);
  H3chFF->SetParLimits(10,0.000001,10000);
  H3chFF->SetParLimits(11,0.000001,10000);

  H3chFF->SetParLimits(12,0.000001,10000);
  H3chFF->SetParLimits(13,0.000001,10000);
  H3chFF->SetParLimits(14,0.000001,10000);
  H3chFF->SetParLimits(15,0.000001,10000);
  H3chFF->SetParLimits(16,0.000001,10000);
  H3chFF->SetParLimits(17,0.000001,10000);
  H3chFF->SetParLimits(18,0.000001,10000);
  H3chFF->SetParLimits(19,0.000001,10000);
  H3chFF->SetParLimits(20,0.000001,10000);
  H3chFF->SetParLimits(21,0.000001,10000);
  H3chFF->SetParLimits(22,0.000001,10000);
  H3chFF->SetParLimits(23,0.000001,10000);
*/

  //graph->Fit(H3chFF,"0");
  graph->Draw("same");
  H3chFF->SetLineColor(3);
  H3chFF->Draw("same");
  cout<<H3chFF->Eval(42.)<<endl;

  
 //Define fit function for H3 magnetic FF as a function of Q^2.
  Double_t fitm(Double_t *Q2, Double_t *par)
  {
    Double_t fitval = 0.;
    Double_t sumH3mtemp = 0.;
    
    for(Int_t i=0; i<ngaus; i++)
      {
	//sumH3chtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Fit R[i] values and Qi.
	sumH3mtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*par[ngaus+i])/(pow(Q2[0],0.5)*par[ngaus+i])) );

	fitval = fitval + sumH3mtemp;
      }
    
    fitval = fabs( fitval ) * exp(-0.25*Q2[0]*pow(gamma,2.0));
    return fitval;
  }
  
  
  TF1 *H3mFF = new TF1("H3mFF",fitm, ymin, ymax,npar);
  //H3chFF->SetParLimits(0,0.00001,100.);
  //c1->SetLogy();
  
  /*
  H3mFF->SetParameter(0,1.);
  H3mFF->SetParameter(1,1.);
  H3mFF->SetParameter(2,1.);
  H3mFF->SetParameter(3,1.);
  H3mFF->SetParameter(4,1.);
  H3mFF->SetParameter(5,1.);
  H3mFF->SetParameter(6,1.);
  H3mFF->SetParameter(7,1.);
  H3mFF->SetParameter(8,1.);
  H3mFF->SetParameter(9,1.);
  H3mFF->SetParameter(10,1.);
  H3mFF->SetParameter(11,1.);
*/
  
  H3mFF->SetParameter(12,R[0]);
  H3mFF->SetParameter(13,R[1]);
  H3mFF->SetParameter(14,R[2]);
  H3mFF->SetParameter(15,R[3]);
  H3mFF->SetParameter(16,R[4]);
  H3mFF->SetParameter(17,R[5]);
  H3mFF->SetParameter(18,R[6]);
  H3mFF->SetParameter(19,R[7]);
  H3mFF->SetParameter(20,R[8]);
  H3mFF->SetParameter(21,R[9]);
  H3mFF->SetParameter(22,R[10]);
  H3mFF->SetParameter(23,R[11]);
  
  /*
  H3mFF->SetParLimits(0,0.000001,10000);
  H3mFF->SetParLimits(1,0.000001,10000);
  H3mFF->SetParLimits(2,0.000001,10000);
  H3mFF->SetParLimits(3,0.000001,10000);
  H3mFF->SetParLimits(4,0.000001,10000);
  H3mFF->SetParLimits(5,0.000001,10000);
  H3mFF->SetParLimits(6,0.000001,10000);
  H3mFF->SetParLimits(7,0.000001,10000);
  H3mFF->SetParLimits(8,0.000001,10000);
  H3mFF->SetParLimits(9,0.000001,10000);
  H3mFF->SetParLimits(10,0.000001,10000);
  H3mFF->SetParLimits(11,0.000001,10000);

  H3mFF->SetParLimits(12,0.000001,10000);
  H3mFF->SetParLimits(13,0.000001,10000);
  H3mFF->SetParLimits(14,0.000001,10000);
  H3mFF->SetParLimits(15,0.000001,10000);
  H3mFF->SetParLimits(16,0.000001,10000);
  H3mFF->SetParLimits(17,0.000001,10000);
  H3mFF->SetParLimits(18,0.000001,10000);
  H3mFF->SetParLimits(19,0.000001,10000);
  H3mFF->SetParLimits(20,0.000001,10000);
  H3mFF->SetParLimits(21,0.000001,10000);
  H3mFF->SetParLimits(22,0.000001,10000);
  H3mFF->SetParLimits(23,0.000001,10000);
*/

  graph->Fit(H3mFF,"0");
  //graph->Draw("same");
  H3mFF->SetLineColor(4);
  H3mFF->Draw("same");
  cout<<H3mFF->Eval(42.)<<endl;





  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  
Double_t mottxs(Double_t *angle2, Double_t *par)
  {
    Double_t val = 0.; 
    Double_t angle = 0.;
    Ef = E0/(1.0+2.0*E0*pow(sin(angle2[0]*deg2rad/2.0),2.0)/MtC);

    val = (  (pow(6.,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle2[0]*deg2rad/2.0),4.0)))*pow(cos(angle2[0]*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

    return val;
  }

  TF1 *fmottxs = new TF1("fmottxs",mottxs, ymin, ymax,1);
  cout<<"!!! = "<<fmottxs->Eval(40.)<<endl;
  //fmottxs->Draw();

 Double_t xs(Double_t *angle3, Double_t *par)
  {
    Double_t val = 0.;
    Ef = E0/(1.0+2.0*E0*pow(sin(angle3[0]*deg2rad/2.0),2.0)/MtC);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle3[0]*deg2rad/2.0),2.0) * GeV2fm;
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*6.*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(12.0,1.0/3.0))) ,2.0);   //A=6 Z=12
                
    Double_t W = E0 - Ef;
    //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
    Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    Double_t eta = 1.0 + Q2eff/(4.0*pow(MtC,2.0)*GeV2fm);        //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.

    //val = 1.1*fmottxs->Eval(angle3[0]);
    val = fmottxs->Eval(angle3[0]) * (1./eta) * ( (Q2/q2_3)*H3chFF->Eval(angle3[0]) ); //magnetic moment for C12 is 0 -> no mag part of XS.

    return val;
  }
 
 c2->SetLogy(); 
 graph->Draw();
 TF1 *fxs = new TF1("fxs",xs, ymin, ymax,24);

  
 fxs->SetParameter(0,1.);
 fxs->SetParameter(1,1.);
 fxs->SetParameter(2,1.);
 fxs->SetParameter(3,1.);
 fxs->SetParameter(4,1.);
 fxs->SetParameter(5,1.);
 fxs->SetParameter(6,1.);
 fxs->SetParameter(7,1.);
 fxs->SetParameter(8,1.);
 fxs->SetParameter(9,1.);
 fxs->SetParameter(10,1.);
 fxs->SetParameter(11,1.);
  
 
 fxs->SetParameter(12,R[0]);
 fxs->SetParameter(13,R[1]);
 fxs->SetParameter(14,R[2]);
 fxs->SetParameter(15,R[3]);
 fxs->SetParameter(16,R[4]);
 fxs->SetParameter(17,R[5]);
 fxs->SetParameter(18,R[6]);
 fxs->SetParameter(19,R[7]);
 fxs->SetParameter(20,R[8]);
 fxs->SetParameter(21,R[9]);
 fxs->SetParameter(22,R[10]);
 fxs->SetParameter(23,R[11]);
  
 graph->Fit(fxs);
 fxs->Draw("same");


  /*
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
  */
}
