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
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t e = 1.60217662E-19;             //Electron charge C.
Double_t alpha = 0.0072973525664;//1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.
Double_t mu3H = 2.9788*(3.0/1.0);

Int_t loops = 1;
const Int_t datapts = 259;//248
Int_t userand = 3;                       //0 = use predetermined Ri from Amroun. 1 = use random Ri in generated in a range around Amroun's. 2 = use random Ri generated in increments of 0.1 with larger possible spacing at greater radii. 3 = use predetermined Ri for the purposes of trying to tune the fit by hand.
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 0;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t fft = 0;                           //0 = don't use FFT to try to get a charge radii. 1 = do use FFT to extract a charge radii.
Int_t Amroun_Qi = 0;                     //1 = Override fitted Qi and use Amroun's values.
Int_t showplots = 1;
Int_t useFB = 1;                         //Turn on Fourier Bessel fit.
Int_t useFB_GM = 1;                      //0 = Turn on Fourier Bessel fit just for GE. 1 = Turn on Fourier Bessel fit attempting GE and GM.
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 8;                        //Number of Gaussians used to fit data.
Int_t ngaus_Amroun = 10;                        //Number of Gaussians used to fit data from Amroun.
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 1.;                         //Atomic number H3.
Double_t A = 3.;                        //Mass number H3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
Double_t Mt3H = 3.0160492*0.9315;         //Mass of H3 in GeV.
Double_t gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
//Double_t E0 = 0.5084;                    //Initial e- energy GeV.
Double_t Ef = 0.;                        //Final e- energy GeV.
Double_t ymin = 30.;//30
Double_t ymax = 100.;//100
Double_t yminFF = 0.;//30
Double_t ymaxFF = 6.;
Double_t range = fabs(ymaxFF - yminFF);
Int_t n = 10000;
Int_t ndim = n+1;
Int_t npdraw = 10001;                     //Number of points to be used when drawing a function.
Double_t truncate = 100.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 2.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char* str[1000];                          //Variable to read lines of the data file.
Float_t thetatemp,qefftemp,sigexptemp,uncertaintytemp,E0temp;
Float_t theta[1000];                     //Angle in degrees.
Float_t qeff[1000];                      //q effective in fm^-1.
Float_t sigexp[1000];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[1000];
Float_t E0[1000];
Float_t Q2[datapts];

Int_t Amroun_pts = 57;
Int_t Collard_pts = 118;
Int_t Szlata_pts = 22;
Int_t Dunn_pts = 27;
Int_t Camsonne_pts = 16;
Int_t Nakagawa_pts = 5;
Int_t my_pts = 1;

Double_t m = 2.;
//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].//0.8,1.4,1.7
//Double_t R[12] = {0.2,0.7,1.3,1.5,2.1,2.8,3.6,4.2,5.2,0.,0.,0.};//{0.3,0.7,1.3,2.,2.7,3.6,4.3,5.5,0.,0.,0.,0.};//{0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit

Double_t R[12] = {0.3,0.8,1.4,1.9,2.5,3.3,4.1,4.8}; //My Fit.
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit

Double_t Qich[12] = {0.151488,0.348372,0.29635,0.0978631,0.121983,0.0242654,0.049329,4.40751e-11}; //My Fit
Double_t Qim[12] = {0.190646,0.301416,0.318972,0.159433,0.173933,0.106361,0.0665564,0.0148866};//My Fit

//Double_t Qich[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121}; //Amroun 3H
//Double_t Qim[12] = {0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077}; //Amroun 3H

Double_t Qich_Amroun[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121};//Amroun 3H
Double_t Qim_Amroun[12] = {0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077};//Amroun 3H
Double_t av[24] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t averr[24] = {};
Double_t Qicherr[12]={}; 
Double_t Qimerr[12]={};
Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};
Double_t Chi2_FB[datapts]={};
Double_t residual_FB[datapts]={};
Double_t FBfit[datapts]={};

//Adding FF world data for plots.
const Int_t size1 = 100;
Int_t nlines1,nlines2,nlines3,nlines4,nlines5,nlines6,nlines7,nlines8,nlines9,nlines10,nlines11,nlines12,nlines13,nlines14,nlines15,nlines16;
Int_t skip1 = 0,skip2 = 0,skip3 = 0,skip4 = 0,skip5 = 0,skip6 = 0,skip7 = 0,skip8 = 0,skip9 = 0,skip10 = 0,skip11 = 0,skip12 = 0,skip13 = 0,skip14 = 0,skip15 = 0,skip16 = 0;
char* str1[1000],str2[1000],str3[1000],str4[1000],str5[1000],str6[1000],str7[1000],str8[1000],str9[1000],str10[1000],str11[1000],str12[1000],str13[1000],str14[1000],str15[1000],str16[1000];
Float_t q2_col65[size1],fch_col65[size1],dq2_col65[size1],dfch_col65[size1],fm_col65[size1],dfm_col65[size1];
Float_t q2_col65_temp,fch_col65_temp,dfch_col65_temp,fm_col65_temp,dfm_col65_temp;
Float_t q2_fch_bec84[size1],fch_bec84[size1],dq2_fch_bec84[size1],dfch_bec84[size1];
Float_t q2_fch_bec84_temp,fch_bec84_temp,dfch_bec84_temp;
Float_t q2_fm_bec84[size1],fm_bec84[size1],dq2_fm_bec84[size1],dfm_bec84[size1];
Float_t q2_fm_bec84_temp,fm_bec84_temp,dfm_bec84_temp;

void Plot_FFs_3H() 
{
  //Open and read in files for Collard 1965.
  FILE *fp1;
  fp1 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Collard1965_3H.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines1 < skip1)
	{
	  fgets(str1,1000,fp1);
	  nlines1++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp1,"%f %f %f %f %f", &q2_col65_temp, &fch_col65_temp, &dfch_col65_temp, &fm_col65_temp, &dfm_col65_temp);
	  if (ncols < 0) break;    
  
	  q2_col65[nlines1-skip1] = q2_col65_temp;
	  fch_col65[nlines1-skip1] = fch_col65_temp;
	  dfch_col65[nlines1-skip1] = dfch_col65_temp;
	  fm_col65[nlines1-skip1] = fm_col65_temp;
	  dfm_col65[nlines1-skip1] = dfm_col65_temp;
	  dq2_col65[nlines1-skip1] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines1++;
	}
    }
  fclose(fp1);
  cout<<"nlines1 = "<<nlines1<<endl;

  //Open and read in files for Beck 1984 Fch.
  FILE *fp2;
  fp2 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Beck1984_Fch.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines2 < skip2)
	{
	  fgets(str2,1000,fp2);
	  nlines2++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp2,"%f %f %f", &q2_fch_bec84_temp, &fch_bec84_temp, &dfch_bec84_temp);
	  if (ncols < 0) break;    
  
	  q2_fch_bec84[nlines2-skip2] = q2_fch_bec84_temp;
	  fch_bec84[nlines2-skip2] = fch_bec84_temp;
	  dfch_bec84[nlines2-skip2] = dfch_bec84_temp;
	  dq2_fch_bec84[nlines2-skip2] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines2++;
	}
    }
  fclose(fp2);
  cout<<"nlines2 = "<<nlines2<<endl;

  //Open and read in files for Beck 1984 Fm.
  FILE *fp3;
  fp3 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Beck1984_Fm.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines3 < skip3)
	{
	  fgets(str3,1000,fp3);
	  nlines3++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp3,"%f %f %f", &q2_fm_bec84_temp, &fm_bec84_temp, &dfm_bec84_temp);
	  if (ncols < 0) break;    
  
	  q2_fm_bec84[nlines3-skip3] = q2_fm_bec84_temp;
	  fm_bec84[nlines3-skip3] = fm_bec84_temp;
	  dfm_bec84[nlines3-skip3] = dfm_bec84_temp;
	  dq2_fm_bec84[nlines3-skip3] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines3++;
	}
    }
  fclose(fp3);
  cout<<"nlines3 = "<<nlines3<<endl;

  TCanvas* cFch=new TCanvas("cFch");
  cFch->SetGrid();
  cFch->SetLogy();

  //Plot Charge FF Fch(Q^2) fm^-2.
  Double_t ChFF_Q2(Double_t *Q2, Double_t *par)
  {
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;

    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 	
	//Use SOG fit for C12 Qi coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Convert to fm. Not sure I need to do this.
	//sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0]*GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0]*GeV2fm,0.5)*R[i])/(pow(Q2[0]*GeV2fm,0.5)*R[i])) );
	sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
	fitch = fitch + sumchtemp;
      }
    //Convert to fm. Not sure I need to do this.
    //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
    fitch = fitch * exp(-0.25*Q2[0]*pow(gamma,2.0));
    fitch = fabs(fitch);
    return fitch;
  }

  TF1 *fChFF = new TF1("fChFF",ChFF_Q2,yminFF,ymaxFF+54,1);
  fChFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
  fChFF->Draw("L");
  cFch->SetTitle("Charge Form Factor");
  //fChFF->SetTitle("C12 Charge Form Factor","#Q^2 (#fm^-2)","#F_{Ch}(q)");
  fChFF->SetTitle("^{3}H Charge Form Factor");
  fChFF->GetHistogram()->GetYaxis()->SetTitle("|F_{ch}(q^{2})|");
  fChFF->GetHistogram()->GetYaxis()->CenterTitle(true);
  fChFF->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  fChFF->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fChFF->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fChFF->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
  fChFF->GetHistogram()->GetXaxis()->CenterTitle(true);
  fChFF->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  fChFF->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  fChFF->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);


  //Plot Amroun's charge FF. No idea whay I can't just redefine Qi from ChFF_Q2.
  Double_t ChFF_Q2_Amroun(Double_t *Q2, Double_t *par)
  {
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;

    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus_Amroun; i++)
      { 	
	//Use SOG fit for C12 Qi coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Convert to fm. Not sure I need to do this.
	//sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0]*GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0]*GeV2fm,0.5)*R[i])/(pow(Q2[0]*GeV2fm,0.5)*R[i])) );
	sumchtemp = (Qich_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
	
	fitch = fitch + sumchtemp;
      }
    //Convert to fm. Not sure I need to do this.
    //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
    fitch = fitch * exp(-0.25*Q2[0]*pow(gamma,2.0));
    fitch = fabs(fitch);
    return fitch;
  }

  TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,yminFF,ymaxFF+54,1);
  //cout<<fChFF_Amroun->Eval(35)<<endl;
  fChFF_Amroun->SetNpx(npdraw);
  fChFF_Amroun->SetLineColor(4);
  fChFF_Amroun->Draw("L same");

  TGraphErrors *gr1 = new TGraphErrors (nlines1, q2_col65, fch_col65, dq2_col65, dfch_col65); 
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);
  gr1->Draw("same p");

  TGraphErrors *gr2 = new TGraphErrors (nlines1, q2_fch_bec84, fch_bec84, dq2_fch_bec84, dfch_bec84); 
  gr2->SetMarkerColor(1);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1);
  gr2->Draw("same p");

  auto ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
  ChFF_leg->AddEntry("fChFF","New ^{3}H |F_{ch}(q^{2})| Fit","l");
  ChFF_leg->AddEntry("fChFF_Amroun","^{3}H |F_{ch}(q^{2})| Fit from Amroun et al.","l");
  ChFF_leg->AddEntry(gr1,"Collard 1965","p");
  ChFF_leg->AddEntry(gr2,"Beck 1984","p");
  ChFF_leg->Draw();


  TCanvas* cFm=new TCanvas("cFm");
  cFm->SetGrid();
  cFm->SetLogy();
  
  //Plot magnetic FF(Q^2) fm^-2.
  Double_t MFF_Q2(Double_t *Q2, Double_t *par)
  {
    Double_t fitm = 0.;
    Double_t summtemp = 0.;
       
    //Define SOG for magnetic FF.
    for(Int_t i=0; i<ngaus; i++)
      { 	
	//Use SOG fit for C12 Qi coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
	summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
	fitm = fitm + summtemp;
      }
       
    fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
    fitm = fabs(fitm);
    return fitm;
  }

  TF1 *fMFF = new TF1("fMFF",MFF_Q2,0.,60.,1);
  fMFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
  fMFF->Draw("L");
  fMFF->SetTitle("^{3}H Magnetic Form Factor");
  fMFF->GetHistogram()->GetYaxis()->SetTitle("|F_{m}(q^{2})|");
  fMFF->GetHistogram()->GetYaxis()->CenterTitle(true);
  fMFF->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  fMFF->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fMFF->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fMFF->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
  fMFF->GetHistogram()->GetXaxis()->CenterTitle(true);
  fMFF->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  fMFF->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  fMFF->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

  //Define Amroun's FF. Not sure why I can't just change Qi in MFF_Q2.
  Double_t MFF_Q2_Amroun(Double_t *Q2, Double_t *par)
  {
    Double_t fitm = 0.;
    Double_t summtemp = 0.;
       
    //Define SOG for magnetic FF.
    for(Int_t i=0; i<ngaus_Amroun; i++)
      { 	
	//Use SOG fit for C12 Qi coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
	summtemp = (Qim_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
	   
	fitm = fitm + summtemp;
      }
       
    fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
    fitm = fabs(fitm);
    return fitm;
  }

  TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,ymaxFF+54,1);
  //cout<<fMFF_Amroun->Eval(30)<<endl;
  fMFF_Amroun->SetNpx(npdraw);
  fMFF_Amroun->SetLineColor(4);
  fMFF_Amroun->Draw("L same");

 TGraphErrors *gr3 = new TGraphErrors (nlines1, q2_col65, fm_col65, dq2_col65, dfm_col65); 
  gr3->SetMarkerColor(4);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(1);
  gr3->Draw("same p");

  TGraphErrors *gr4 = new TGraphErrors (nlines3, q2_fm_bec84, fm_bec84, dq2_fm_bec84, dfm_bec84); 
  gr4->SetMarkerColor(1);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(1);
  gr4->Draw("same p");

  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  MFF_leg->AddEntry("fMFF","New ^{3}H |F_{m}(q^{2})| Fit","l");
  MFF_leg->AddEntry("fMFF_Amroun","^{3}H |F_{m}(q^{2})| Fit from Amroun et al.","l");
  MFF_leg->AddEntry(gr3,"Collard 1965","p");
  MFF_leg->AddEntry(gr4,"Beck 1984","p");
  MFF_leg->Draw();
}
