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
Double_t muHe3 = -2.1275*(3.0/2.0); //3He -2.1275*(3.0/2.0). Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive. //3H 2.9788*(3.0/1.0)
Int_t target = 1; //3He = 0. 3H = 1. 

const Int_t nfunc = 3000;
Double_t maxchi2 = 603;//3H 611.70 n=7 100//3H 603 n=8 100//3H 604 n=9 100//3H 603 n=10 100//3H 602 n=11 100//3He 765 n=8 100 //3He 521 n=9 100 //3He 519 n=10 100 //3He 503 n=11 100 //3He 501 n=12 100//3He 500 n=13 100//My old point for combined 3He 505, 3H 603   //Max chi2 value above which fits are removed from the analysis.
Double_t Qim_range = 50.; //Determines the amount above or below 1 the sum of the magnetic Qi may have and be accepted. (Note Qich is consistently close to 1 so it is not cut on.
Int_t loops = 1;
Int_t current_loop = 0;
const Int_t datapts = 259;//248
const Int_t size = 3000;//248
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
Int_t ngaus_Amroun = 12;                        //Number of Gaussians used to fit data from Amroun.
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //3.0160293*0.9315 Mass of He3 in GeV. 3.0160492*0.9315 mass 3H.
Double_t gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
//Double_t E0 = 0.5084;                    //Initial e- energy GeV.
Double_t Ef = 0.;                        //Final e- energy GeV.
Double_t ymin = 30.;//30
Double_t ymax = 100.;//100
Double_t yminFF = 0.0001;//30
Double_t ymaxFF = 6.;
Double_t range = fabs(ymaxFF - yminFF);
Int_t n = 10000;
Int_t ndim = n+1;
Int_t npdraw = 1001;                     //Number of points to be used when drawing a function.
Double_t transparency = 0.2;//0.05              //Sets the transparency level of the multiplot lines.
Double_t linewidth = 2.;
Double_t truncate = 100.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 1.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char* str[1000];                          //Variable to read lines of the data file.
Float_t Chi2[size],rChi2[size],BIC[size],AIC[size],Qichtot[size],Qimtot[size],R0[size],R1[size],R2[size],R3[size],R4[size],R5[size],R6[size],R7[size],R8[size],R9[size],R10[size],R11[size],R12[size],Q0ch[size],Q1ch[size],Q2ch[size],Q3ch[size],Q4ch[size],Q5ch[size],Q6ch[size],Q7ch[size],Q8ch[size],Q9ch[size],Q10ch[size],Q11ch[size],Q12ch[size],Q0m[size],Q1m[size],Q2m[size],Q3m[size],Q4m[size],Q5m[size],Q6m[size],Q7m[size],Q8m[size],Q9m[size],Q10m[size],Q11m[size],Q12m[size];
Float_t Chi2temp,rChi2temp,BICtemp,AICtemp,Qichtottemp,Qimtottemp,R0temp,R1temp,R2temp,R3temp,R4temp,R5temp,R6temp,R7temp,R8temp,R9temp,R10temp,R11temp,R12temp,Q0chtemp,Q1chtemp,Q2chtemp,Q3chtemp,Q4chtemp,Q5chtemp,Q6chtemp,Q7chtemp,Q8chtemp,Q9chtemp,Q10chtemp,Q11chtemp,Q12chtemp,Q0mtemp,Q1mtemp,Q2mtemp,Q3mtemp,Q4mtemp,Q5mtemp,Q6mtemp,Q7mtemp,Q8mtemp,Q9mtemp,Q10mtemp,Q11mtemp,Q12mtemp;
Float_t theta[datapts];                     //Angle in degrees.
Float_t qeff[datapts];                      //q effective in fm^-1.
Float_t sigexp[datapts];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[datapts];
Float_t E0[datapts];
Float_t Q2[datapts];
Float_t Rmulti[size][15];
Float_t Qichmulti[size][15];
Float_t Qimmulti[size][15];
Float_t rms_deriv[size];                          //Charge radius defined as -6*derivative of Ch FF at Q^2=0.
Float_t rms_int[size];                          //Charge radius defined in the integral manner.

Int_t Amroun_pts = 57;
Int_t Collard_pts = 118;
Int_t Szlata_pts = 22;
Int_t Dunn_pts = 27;
Int_t Camsonne_pts = 16;
Int_t Nakagawa_pts = 5;
Int_t my_pts = 1;

Double_t R[15] = {0.1,0.7,1.3,2.,2.7,3.6,4.4,5.6,0.,0.,0.,0.};//7
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t Qich[15] = {0.0784469,0.247165,0.406019,0.120177,0.137968,4.57535E-11,0.0200847,2.63439E-9,0.,0.,0.,0.};//7
Double_t Qim[15] = {0.0770148,0.298502,0.282963,0.175066,0.0769078,0.0381075,0.0899692,0.0675,0.,0.,0.,0.};

//Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};//3He
//Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};//3He

Double_t Qich_Amroun[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121};//Amroun 3H
Double_t Qim_Amroun[12] = {0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077};//Amroun 3H

Double_t av[24] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t averr[24] = {};
Double_t Qicherr[15]={}; 
Double_t Qimerr[15]={};
//Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};
Double_t Chi2_FB[datapts]={};
Double_t residual_FB[datapts]={};
Double_t FBfit[datapts]={};
Int_t first = 0;                     //Counter to check if this is the first curve to be plotted since it needs draw("L") not draw("L SAME").
Int_t total_funcs = 0;
Double_t Rimax = 6.4;
Int_t Ri_divisions = 64;
Double_t BIC_tot = 0.;               //Total value for Bayesian information criterion.
Double_t AIC_tot = 0.;               //Total value for Akaike information criterion.
Double_t Qich_tot = 0.;              //Total value for sum of Qi values.
Double_t Qim_tot = 0.;
Double_t Chi2_tot = 0.;              //Total value for sum of Chi^2 and rChi^2.
Double_t rChi2_tot = 0;
Double_t min_chi2 = 10000.;              //Variable to store the minimum chi2 fit result.
Int_t min_chi2_fit = 0;                  //Fit number of the min chi2 fit.
Int_t rep_chi2_fit = 556;
Double_t rms_deriv_tot = 0;                //Total value of rms charge radius defined as -6*derivative of Ch FF at Q^2=0.
Double_t rms_deriv_min = 1.85;//1.85
Double_t rms_deriv_max = 2.1;//2.1
Int_t rms_deriv_divisions = 100;
Double_t rms_int_tot = 0;                //Total value of rms charge radius integral method.
Double_t rms_int_min = rms_deriv_min;
Double_t rms_int_max = rms_deriv_max;
Int_t rms_int_divisions = rms_deriv_divisions;
Double_t gaus_min = 0.;              //Function range for gaussian fits of individual Ri histograms.
Double_t gaus_max = 7.;
Double_t Qi_ch_min = -0.1;             //Create ranges for Qi distribution plots.
Double_t Qi_ch_max = 0.5;
Double_t Qi_m_min = -0.1;
Double_t Qi_m_max = 0.5;
Double_t Qi_ch_divisions = 60.;
Double_t Qi_m_divisions = 60.;
Double_t gaus_Qi_min = -0.1;              //Function range for gaussian fits of individual Qi histograms.
Double_t gaus_Qi_max = 0.6;

Int_t z = 0; //Counter for main loop.

//Additional data file (error bands) variables.
const Int_t size1 = 200;
Int_t skip1 = 0;                         //Gives number of lines to skip at top of data file. 
Int_t nlines1 = 0;                        //Counts number of lines in the data file. 
Int_t ncols1;                             //Set how many columns of data we have in the data file.
char* str1[1000];                         //Variable to read lines of the data file.
Float_t x_3H_Fch_up[size1],y_3H_Fch_up[size1];    //Arrays for x,y of the 3H Fch upper error band.
Float_t x_3H_Fch_up_temp,y_3H_Fch_up_temp;


//Plot Charge FF Fch(Q^2) fm^-2.
Double_t ChFF_Q2(Double_t *Q2, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
  //cout<<"Current Loop = "<<current_loop<<endl;
  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      sumchtemp = (par[i]/(1.0+2.0*pow(par[i+ngaus],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*par[i+ngaus]) + (2.0*pow(par[i+ngaus],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*par[i+ngaus])/(pow(Q2[0],0.5)*par[i+ngaus])) );
	
      fitch = fitch + sumchtemp;
    }

  fitch = fitch * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}

//Plot Amroun's charge FF. No idea whay I can't just redefine Qi from ChFF_Q2.
Double_t ChFF_Q2_Amroun(Double_t *Q2, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus_Amroun; i++)
    { 	
      sumchtemp = (Qich_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
	
      fitch = fitch + sumchtemp;
    }

  fitch = fitch * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}

Double_t MFF_Q2(Double_t *Q2, Double_t *par)
{
  Double_t fitm = 0.;
  Double_t summtemp = 0.;
  //cout<<"Current Loop = "<<current_loop<<endl;
  //Define SOG for magnetic FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      summtemp = (par[i]/(1.0+2.0*pow(par[i+ngaus],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*par[i+ngaus]) + (2.0*pow(par[i+ngaus],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*par[i+ngaus])/(pow(Q2[0],0.5)*par[i+ngaus])) );
	
      fitm = fitm + summtemp;
      //cout<<"Loop "<<i<<" fitm = "<<fitm<<endl;
    }

  fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitm = fabs(fitm);
  //cout<<"fitm = "<<fitm<<endl;
  return fitm;
}

//Define Amroun's FF. Not sure why I can't just change Qi in MFF_Q2.
Double_t MFF_Q2_Amroun(Double_t *Q2, Double_t *par)
{
  Double_t fitm = 0.;
  Double_t summtemp = 0.;
       
  //Define SOG for magnetic FF.
  for(Int_t i=0; i<ngaus_Amroun; i++)
    { 	 
      summtemp = (Qim_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
	   
      fitm = fitm + summtemp;
    }
       
  fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitm = fabs(fitm);
  return fitm;
}

//Define the charge density from I. Sick. 
Double_t rho_ch(Double_t *r, Double_t *par)
{
  Double_t rho = 0;
  Double_t rho_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_temp = par[i]/( 1+2*pow(par[i+ngaus],2.)/pow(gamma,2.) ) * (  exp( -pow((r[0]-par[i+ngaus]),2.)/pow(gamma,2.) ) + exp( -pow((r[0]+par[i+ngaus]),2.)/pow(gamma,2.) )  );
      rho = rho + rho_temp;
    }

  rho = Z/(2*pow(pi,1.5)*pow(gamma,3.)) * rho; //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho;
}

Double_t ChFF_Deriv(Double_t Q2) 
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
   
  //Define SOG for charge FF.

  fitch = (Q0ch[z]/(1.0+2.0*pow(R0[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R0[z]) + (2.0*pow(R0[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R0[z])/(pow(Q2,0.5)*R0[z])) ) 
    + (Q1ch[z]/(1.0+2.0*pow(R1[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R1[z]) + (2.0*pow(R1[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R1[z])/(pow(Q2,0.5)*R1[z])) )
    + (Q2ch[z]/(1.0+2.0*pow(R2[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R2[z]) + (2.0*pow(R2[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R2[z])/(pow(Q2,0.5)*R2[z])) )
    + (Q3ch[z]/(1.0+2.0*pow(R3[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R3[z]) + (2.0*pow(R3[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R3[z])/(pow(Q2,0.5)*R3[z])) )
    + (Q4ch[z]/(1.0+2.0*pow(R4[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R4[z]) + (2.0*pow(R4[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R4[z])/(pow(Q2,0.5)*R4[z])) )
    + (Q5ch[z]/(1.0+2.0*pow(R5[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R5[z]) + (2.0*pow(R5[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R5[z])/(pow(Q2,0.5)*R5[z])) )
    + (Q6ch[z]/(1.0+2.0*pow(R6[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R6[z]) + (2.0*pow(R6[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R6[z])/(pow(Q2,0.5)*R6[z])) )
    + (Q7ch[z]/(1.0+2.0*pow(R7[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R7[z]) + (2.0*pow(R7[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R7[z])/(pow(Q2,0.5)*R7[z])) );
    //+ (Q8ch[z]/(1.0+2.0*pow(R8[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R8[z]) + (2.0*pow(R8[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R8[z])/(pow(Q2,0.5)*R8[z])) )
    //+ (Q9ch[z]/(1.0+2.0*pow(R9[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R9[z]) + (2.0*pow(R9[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R9[z])/(pow(Q2,0.5)*R9[z])) )
    //+ (Q10ch[z]/(1.0+2.0*pow(R10[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R10[z]) + (2.0*pow(R10[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R10[z])/(pow(Q2,0.5)*R10[z])) )
    //+ (Q11ch[z]/(1.0+2.0*pow(R11[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R11[z]) + (2.0*pow(R11[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R11[z])/(pow(Q2,0.5)*R11[z])) );
  // + (Q12ch[z]/(1.0+2.0*pow(R12[z],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2,0.5)*R12[z]) + (2.0*pow(R12[z],2.0)/pow(gamma,2.0)) * (sin(pow(Q2,0.5)*R12[z])/(pow(Q2,0.5)*R12[z])) );//Need to make this smart badly. Add loop and set the pars to the Ri and Qi.
 
  fitch = fitch * exp(-0.25*Q2*pow(gamma,2.0));
  //fitch = fabs(fitch);
  return fitch;
}

//Create a function that can be integrated to check that the normilaization to Ze is correct.
Double_t rho_ch_int(Double_t *r, Double_t *par)
{
  Double_t rho_int = 0;
  Double_t rho_int_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_int_temp = par[i]/( 1+2*pow(par[i+ngaus],2.)/pow(gamma,2.) ) * (  exp( -pow((r[0]-par[i+ngaus]),2.)/pow(gamma,2.) ) + exp( -pow((r[0]+par[i+ngaus]),2.)/pow(gamma,2.) )  );
      rho_int = rho_int + rho_int_temp;
    }

  rho_int = Z/(2*pow(pi,1.5)*pow(gamma,3.)) * rho_int * 4*pi*pow(r[0],2.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho_int;
}

//Create a function to calculate rms radius.
Double_t rho_rms(Double_t *r, Double_t *par)
{
  Double_t rho_rms = 0;
  Double_t rho_rms_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_rms_temp = par[i]/( 1+2*pow(par[ngaus+i],2.)/pow(gamma,2.) ) * (  exp( -pow((r[0]-par[ngaus+i]),2.)/pow(gamma,2.) ) + exp( -pow((r[0]+par[ngaus+i]),2.)/pow(gamma,2.) )  );
      rho_rms = rho_rms + rho_rms_temp;
    }

  rho_rms = Z/(2*pow(pi,1.5)*pow(gamma,3.)) * rho_rms * 4*pi*pow(r[0],4.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho_rms;
}

//Create Gaussian fit for the elastic peak.
Double_t fit_gaus(Double_t *x,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((x[0]-par[1])/par[2]),2));
  return fitval;
}

//Create a Poisson fit.
Double_t poisson(Double_t*x,Double_t*par)                                         
{                                                                              
  return par[0]*TMath::Poisson(x[0],par[1]);
}  

void Multifit_FF_Plots() 
{

  FILE *fp;
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Chi2.txt","r");
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_BS_300_Ri_Chi2.txt","r");
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_Ri_Fits_180_9_25_2018.txt","r");
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Combined_Ri_Fits.txt","r");//My old point.
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_Closer_Ri_Fits_n=10_525_10_1_2018.txt","r");
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_Ri_Fits_n=9_300_9_28_2018.txt","r");
    
  if(target == 0)
    {
      //3He
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=13_100_12_21_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=12_100_12_17_2018.txt","r");
      fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=12_1352_12_22_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=11_100_12_11_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=10_100_12_11_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=9_100_12_11_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=8_100_12_12_2018.txt","r");
    }
  
  if(target == 1)
    {
      //3H
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_3H_Ri_Fits_n=10_100_10_15_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=7_100_12_19_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=8_100_12_12_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=8_Wider_Ri_100_12_20_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=8_2600_12_22_2018.txt","r");
      fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=8_Short_12_22_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=9_100_12_12_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=10_100_10_15_2018.txt","r");
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=11_100_12_12_2018.txt","r");
    }

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines < skip)
	{
	  fgets(str,1000,fp);
	  nlines++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  
	  //Reading in the normal datafile with max ngaus=12.
	  ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp, &R3temp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &Q0chtemp, &Q1chtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp, &Q11chtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp, &Q8mtemp, &Q9mtemp);
	  ncols = ncols + fscanf(fp,"%f %f", &Q10mtemp, &Q11mtemp);
	  
	  //If using n=13 need to turn on reading the n=13 parameter columns.
	  /*
	    ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp, &R3temp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &R12temp, &Q0chtemp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q1chtemp, &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q11chtemp, &Q12chtemp, &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f", &Q8mtemp, &Q9mtemp, &Q10mtemp, &Q11mtemp, &Q12mtemp);
	  */
	  
	  //cout<<"ncols = "<<ncols<<endl;
	  if (ncols < 0) break;    
	  
	  Chi2[nlines-skip] = Chi2temp;
	  rChi2[nlines-skip] = rChi2temp;
	  BIC[nlines-skip] = BICtemp;
	  AIC[nlines-skip] = AICtemp;
	  Qichtot[nlines-skip] = Qichtottemp;
	  Qimtot[nlines-skip] = Qimtottemp;
	  R0[nlines-skip] = R0temp;
	  R1[nlines-skip] = R1temp;
	  R2[nlines-skip] = R2temp;
	  R3[nlines-skip] = R3temp;
	  R4[nlines-skip] = R4temp;
	  R5[nlines-skip] = R5temp;
	  R6[nlines-skip] = R6temp;
	  R7[nlines-skip] = R7temp;
	  R8[nlines-skip] = R8temp;
	  R9[nlines-skip] = R9temp;
	  R10[nlines-skip] = R10temp;
	  R11[nlines-skip] = R11temp;
	  //R12[nlines-skip] = R12temp;
	  Q0ch[nlines-skip] = Q0chtemp;
	  Q1ch[nlines-skip] = Q1chtemp;
	  Q2ch[nlines-skip] = Q2chtemp;
	  Q3ch[nlines-skip] = Q3chtemp;
	  Q4ch[nlines-skip] = Q4chtemp;
	  Q5ch[nlines-skip] = Q5chtemp;
	  Q6ch[nlines-skip] = Q6chtemp;
	  Q7ch[nlines-skip] = Q7chtemp;
	  Q8ch[nlines-skip] = Q8chtemp;
	  Q9ch[nlines-skip] = Q9chtemp;
	  Q10ch[nlines-skip] = Q10chtemp;
	  Q11ch[nlines-skip] = Q11chtemp;
	  //Q12ch[nlines-skip] = Q12chtemp;
	  Q0m[nlines-skip] = Q0mtemp;
	  Q1m[nlines-skip] = Q1mtemp;
	  Q2m[nlines-skip] = Q2mtemp;
	  Q3m[nlines-skip] = Q3mtemp;
	  Q4m[nlines-skip] = Q4mtemp;
	  Q5m[nlines-skip] = Q5mtemp;
	  Q6m[nlines-skip] = Q6mtemp;
	  Q7m[nlines-skip] = Q7mtemp;
	  Q8m[nlines-skip] = Q8mtemp;
	  Q9m[nlines-skip] = Q9mtemp;
	  Q10m[nlines-skip] = Q10mtemp;
	  Q11m[nlines-skip] = Q11mtemp;
	  //Q12m[nlines-skip] = Q12mtemp;
	  Rmulti[nlines-skip][0] = R0temp;
	  Rmulti[nlines-skip][1] = R1temp;
	  Rmulti[nlines-skip][2] = R2temp;
	  Rmulti[nlines-skip][3] = R3temp;
	  Rmulti[nlines-skip][4] = R4temp;
	  Rmulti[nlines-skip][5] = R5temp;
	  Rmulti[nlines-skip][6] = R6temp;
	  Rmulti[nlines-skip][7] = R7temp;
	  Rmulti[nlines-skip][8] = R8temp;
	  Rmulti[nlines-skip][9] = R9temp;
	  Rmulti[nlines-skip][10] = R10temp;
	  Rmulti[nlines-skip][11] = R11temp;
	  //Rmulti[nlines-skip][12] = R12temp;
	  Qichmulti[nlines-skip][0] = Q0chtemp;
	  Qichmulti[nlines-skip][1] = Q1chtemp;
	  Qichmulti[nlines-skip][2] = Q2chtemp;
	  Qichmulti[nlines-skip][3] = Q3chtemp;
	  Qichmulti[nlines-skip][4] = Q4chtemp;
	  Qichmulti[nlines-skip][5] = Q5chtemp;
	  Qichmulti[nlines-skip][6] = Q6chtemp;
	  Qichmulti[nlines-skip][7] = Q7chtemp;
	  Qichmulti[nlines-skip][8] = Q8chtemp;
	  Qichmulti[nlines-skip][9] = Q9chtemp;
	  Qichmulti[nlines-skip][10] = Q10chtemp;
	  Qichmulti[nlines-skip][11] = Q11chtemp;
	  //Qichmulti[nlines-skip][12] = Q12chtemp;
	  Qimmulti[nlines-skip][0] = Q0mtemp;
	  Qimmulti[nlines-skip][1] = Q1mtemp;
	  Qimmulti[nlines-skip][2] = Q2mtemp;
	  Qimmulti[nlines-skip][3] = Q3mtemp;
	  Qimmulti[nlines-skip][4] = Q4mtemp;
	  Qimmulti[nlines-skip][5] = Q5mtemp;
	  Qimmulti[nlines-skip][6] = Q6mtemp;
	  Qimmulti[nlines-skip][7] = Q7mtemp;
	  Qimmulti[nlines-skip][8] = Q8mtemp;
	  Qimmulti[nlines-skip][9] = Q9mtemp;
	  Qimmulti[nlines-skip][10] = Q10mtemp;
	  Qimmulti[nlines-skip][11] = Q11mtemp;
	  //Qimmulti[nlines-skip][12] = Q12mtemp;
	  
	  //cout<<"!!! Chi2["<<nlines-skip<<"]"<<Chi2[nlines-skip]<<endl;

	  nlines++;
	}
    }
  
  //Add in summing up individual Qich,m to account for the mistake where higher n fits leave values for lower n fits.
  for(Int_t i=0;i<(nlines-skip);i++)
    {
      
      Qichtot[i] = 0.;
      Qimtot[i] = 0.;
      for(Int_t j=0;j<ngaus;j++)
	{
	  Qichtot[i] = Qichtot[i] + Qichmulti[i][j];
	  Qimtot[i] = Qimtot[i] + Qimmulti[i][j];
	  //cout<<"!!!!!!!!!!!!!!!"<<Qichtot[i]<<endl;
	}
    }
  //Print the data read from the file.
  if(showplots == 1)
    { 
      for(int i=0; i<(nlines-skip); i++)
	{
	  cout<<"Fit# = "<<i<<"  Chi2 = "<<Chi2[i]<<"  rChi2 = "<<rChi2[i]<<"  BIC = "<<BIC[i]<<"  AIC = "<<AIC[i]<<"  Qichtot = "<<Qichtot[i]<<"  Qimtot = "<<Qimtot[i]<<endl;
	  cout<<endl;
	  cout<<"R[0] = "<<R0[i]<<"  R[1] = "<<R1[i]<<"  R[2] = "<<R2[i]<<"  R[3] = "<<R3[i]<<"  R[4] = "<<R4[i]<<"  R[5] = "<<R5[i]<<"  R[6] = "<<R6[i]<<"  R[7] = "<<R7[i]<<"  R[8] = "<<R8[i]<<"  R[9] = "<<R9[i]<<"  R[10] = "<<R10[i]<<"  R[11] = "<<R11[i]<<endl;
	  cout<<endl;
	  cout<<"Qch[0] = "<<Q0ch[i]<<"  Qch[1] = "<<Q1ch[i]<<"  Qch[2] = "<<Q2ch[i]<<"  Qch[3] = "<<Q3ch[i]<<"  Qch[4] = "<<Q4ch[i]<<"  Qch[5] = "<<Q5ch[i]<<"  Qch[6] = "<<Q6ch[i]<<"  Qch[7] = "<<Q7ch[i]<<"  Qch[8] = "<<Q8ch[i]<<"  Qch[9] = "<<Q9ch[i]<<"  Qch[10] = "<<Q10ch[i]<<"  Qch[11] = "<<Q11ch[i]<<endl;
	  cout<<endl;
	  cout<<"Qm[0] = "<<Q0m[i]<<"  Qm[1] = "<<Q1m[i]<<"  Qm[2] = "<<Q2m[i]<<"  Qm[3] = "<<Q3m[i]<<"  Qm[4] = "<<Q4m[i]<<"  Qm[5] = "<<Q5m[i]<<"  Qm[6] = "<<Q6m[i]<<"  Qm[7] = "<<Q7m[i]<<"  Qm[8] = "<<Q8m[i]<<"  Qm[9] = "<<Q9m[i]<<"  Qm[10] = "<<Q10m[i]<<"  Qm[11] = "<<Q11m[i]<<endl;
	  cout<<"-----------------------------------------------------"<<endl;
	}  
      /*
	for(int i=0; i<(nlines-skip); i++)
	{
	for(int j=0; j<(ngaus); j++)
	{
	cout<<Rmulti[i][j]<<endl;
	}
	cout<<"-----------------------------------------------------"<<endl;
	}
	for(int i=0; i<(nlines-skip); i++)
	{
	  for(int j=0; j<(ngaus); j++)
	    {
	      cout<<Qichmulti[i][j]<<endl;
	    }
	  cout<<"-----------------------------------------------------"<<endl;
	}
      for(int i=0; i<(nlines-skip); i++)
	{
	  for(int j=0; j<(ngaus); j++)
	    {
	      cout<<Qimmulti[i][j]<<endl;
	    }
	  cout<<"-----------------------------------------------------"<<endl;
	}
      */
      
      cout<<"Number of lines = "<<nlines<<endl;
    }
  fclose(fp);
  
  //Open and read in files for Amroun's Error Bands.
  //3H upper Fch error band.
  FILE *fp1;
  fp1 = fopen("/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_Amroun_Error_Band_Upper.txt","r");

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
	  ncols1 = fscanf(fp1,"%f %f", &x_3H_Fch_up_temp, &y_3H_Fch_up_temp);
	  
	  if (ncols1 < 0) break;    
	  
	  x_3H_Fch_up[nlines1-skip1] = x_3H_Fch_up_temp;
	  y_3H_Fch_up[nlines1-skip1] = y_3H_Fch_up_temp;
	  //cout<<"!!! x_3H_Fch_up["<<nlines-skip<<"] = "<<x_3H_Fch_up[nlines-skip]<<"   x_3H_Fch_up_temp["<<nlines-skip<<"] = "<<x_3H_Fch_up_temp<<endl;
	  nlines1++;
	}
    }
  fclose(fp1);

  for(Int_t i=0;i<nlines1;i++)
    {
      cout<<"x_3H_Fch_up["<<i<<"] = "<<x_3H_Fch_up[i]<<"   y_3H_Fch_up["<<i<<"] = "<<y_3H_Fch_up[i]<<endl;
    }

  //Now plot all of the curves on one canvas to form an error band.
  TCanvas* cFch=new TCanvas("cFch");
  cFch->SetGrid();
  cFch->SetLogy();
  cFch->SetTitle("Charge Form Factor");

  TCanvas* cCh_den=new TCanvas("cCh_den");
  cCh_den->SetGrid();
  cCh_den->SetTitle("Charge Density");

  cFch->cd();

  //Define a few analysis histograms to store the Ri values and look for patterns.
  TH1F *hRi = new TH1F("hRi", "Ri Distribution", Ri_divisions, 0., Rimax);

  TH1F **hRi_ind = new TH1F*[ngaus];
  for(Int_t i=0;i<ngaus;i++)
    {
      hRi_ind[i] = new TH1F(Form("hR%d_ind",i), Form("R%d Distributions",i), Ri_divisions, 0., Rimax);
    }

  TH1F **hRi_sep = new TH1F*[ngaus-1];
  for(Int_t i=0;i<(ngaus-1);i++)
    {
      hRi_sep[i] = new TH1F(Form("hR%d_sep",i), Form("R%d to R%d Separation",i,i+1), 15, 0., 1.5);
    }

  TH1F **Qi_ch = new TH1F*[ngaus];
  for(Int_t i=0;i<ngaus;i++)
    {
      Qi_ch[i] = new TH1F(Form("Qi_ch_%d",i), Form("Qi_ch_%d Distributions",i), Qi_ch_divisions, Qi_ch_min, Qi_ch_max);
    }

  TH1F **Qi_m = new TH1F*[ngaus];
  for(Int_t i=0;i<ngaus;i++)
    {
      Qi_m[i] = new TH1F(Form("Qi_m_%d",i), Form("Qi_m_%d Distributions",i), Qi_m_divisions, Qi_m_min, Qi_m_max);
    }
  
  //Define array of TF1 to plot the various fits with.
  TF1 **fChFF = new TF1*[nfunc];
  //TF1 **fChFF_Deriv = new TF1*[nfunc];
  TF1 **frho_ch_int = new TF1*[nfunc];
  TF1 **frho_rms = new TF1*[nfunc];
  TF1 **frho_ch = new TF1*[nfunc];

  //for(current_loop=0; current_loop<1; current_loop++)//nlines-skip;current_loop++)
  for(z=0; z<nlines-skip; z++)
    {
      //cout<<"loop = "<<current_loop<<endl;

      if(current_loop==0)
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);
	  fChFF[current_loop] = new TF1(fname,ChFF_Q2,yminFF,ymaxFF+54,2*ngaus);//was 20.
	  //char fname1[20];
	  //sprintf(fname1,"f%d",z);
	  //fChFF_Deriv[current_loop] = new TF1(fname1,ChFF_Deriv,yminFF,ymaxFF+54,20);
	  char fname1[20];
	  sprintf(fname1,"f%d",z);
	  frho_ch_int[current_loop] = new TF1(fname1,rho_ch_int,0.,10.,2*ngaus+1);
	  char fname2[20];
	  sprintf(fname2,"f%d",z);
	  frho_rms[current_loop] = new TF1(fname2,rho_rms,0.,10.,2*ngaus+1);
	  char fname3[20];
	  sprintf(fname3,"f%d",z);
	  frho_ch[current_loop] = new TF1(fname3,rho_ch,0.,10.,2*ngaus+1);

	  //Set parameters for the various functions.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      fChFF[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      fChFF[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      //fChFF_Deriv[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      //fChFF_Deriv[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      frho_ch_int[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      frho_ch_int[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      frho_rms[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      frho_rms[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      frho_ch[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      frho_ch[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	    }

	  fChFF[current_loop]->SetLineColorAlpha(2,transparency);
	  fChFF[current_loop]->SetLineWidth(linewidth);
	  frho_ch[current_loop]->SetLineColorAlpha(2,transparency);
	  frho_ch[current_loop]->SetLineWidth(linewidth);

	  if(Chi2[current_loop]<maxchi2 && Qimtot[current_loop]>(1-Qim_range) && Qimtot[current_loop]<(1+Qim_range))
	    {
	      cFch->cd();
	      fChFF[current_loop]->Draw("L");
	      cCh_den->cd();
	      frho_ch[current_loop]->Draw("L");

	      if(Chi2[current_loop]<min_chi2)
		{
		  cout<<"*****Hello there!*****"<<endl;
		  min_chi2 = Chi2[current_loop];
		  min_chi2_fit = current_loop;
		}

	      for(Int_t i=0;i<ngaus;i++)
		{
		  hRi->Fill(Rmulti[current_loop][i]);
		  hRi_ind[i]->Fill(Rmulti[current_loop][i]);
		  Qi_ch[i]->Fill(Qichmulti[current_loop][i]);
		  Qi_m[i]->Fill(Qimmulti[current_loop][i]);
		}
	      for(Int_t i=0;i<(ngaus-1);i++)
		{
		  hRi_sep[i]->Fill(Rmulti[current_loop][i+1]-Rmulti[current_loop][i]);
		}

	      //Calculate the charge radii.
	      double x0 = 0.0015;
	      ROOT::Math::Functor1D f1D(&ChFF_Deriv); //Fail fChFF[current_loop]. Fail ChFF_Q2.
	      ROOT::Math::RichardsonDerivator rd;
	      rd.SetFunction(f1D);
	      cout<<"First Derivative:   "<<rd.Derivative1(x0)<<"   -6*dFC(0)/dQ^2 = "<<-6*rd.Derivative1(x0)<<"   rms radius = "<<pow(-6*rd.Derivative1(x0),0.5)<<endl;

	      cout<<"Sum Qich = "<<Qichtot[current_loop]<<endl;
	      first = 1; //No longer first plot.
	      BIC_tot = BIC_tot + BIC[current_loop];         //Sum of BIC.
	      AIC_tot = AIC_tot + AIC[current_loop];         //Sum of AIC.
	      Qich_tot = Qich_tot + Qichtot[current_loop];   //Sum of Qich total values.
	      Qim_tot = Qim_tot + Qimtot[current_loop];      //Sum of Qim total values.
	      Chi2_tot = Chi2_tot + Chi2[current_loop];      //Sum of Chi^2 total values.
	      rChi2_tot = rChi2_tot + rChi2[current_loop];   //Sum of rChi^2 total values.

	      cout<<"Total Charge = "<<frho_ch_int[current_loop]->Integral(0.,10.)<<"e."<<endl;
	      cout<<"rms radius (integral method) = "<<pow( frho_rms[current_loop]->Integral(0.0,10.)/frho_ch_int[current_loop]->Integral(0.0,10.) ,0.5)<<endl;
	      rms_deriv[z] = pow(-6*rd.Derivative1(x0),0.5);
	      rms_deriv_tot = rms_deriv[z] + rms_deriv_tot;
	      rms_int[z] = pow( frho_rms[current_loop]->Integral(0.0,10.)/frho_ch_int[current_loop]->Integral(0.0,10.) ,0.5);
	      rms_int_tot = rms_int[z] + rms_int_tot;
       
	      total_funcs++;
	    }
	  fChFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  if(target == 0)
	    {
	      fChFF[current_loop]->SetTitle("^{3}He Charge Form Factor");
	    }
	  if(target == 1)
	    {
	      fChFF[current_loop]->SetTitle("^{3}H Charge Form Factor");
	    }
	  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetTitle("|F_{ch}(q^{2})|");
	  fChFF[current_loop]->GetHistogram()->GetYaxis()->CenterTitle(true);
	  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
	  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
	  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
	  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
	  fChFF[current_loop]->GetHistogram()->GetXaxis()->CenterTitle(true);
	  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
	  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
	  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

	  frho_ch[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  if(target == 0)
	    {
	      frho_ch[current_loop]->SetTitle("^{3}He Charge Density");
	    }
	  if(target == 1)
	    {
	      frho_ch[current_loop]->SetTitle("^{3}H Charge Density");
	    }
	  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetTitle("e/fm^{3}");
	  frho_ch[current_loop]->GetHistogram()->GetYaxis()->CenterTitle(true);
	  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
	  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
	  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
	  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetTitle("r (fm)");
	  frho_ch[current_loop]->GetHistogram()->GetXaxis()->CenterTitle(true);
	  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
	  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
	  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

	}
      else
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);
	  fChFF[current_loop] = new TF1(fname,ChFF_Q2,yminFF,ymaxFF+54,2*ngaus);
	  //char fname1[20];
	  //sprintf(fname1,"f%d",z);
	  //fChFF_Deriv[current_loop] = new TF1(fname1,ChFF_Deriv,yminFF,ymaxFF+54,20);
	  char fname1[20];
	  sprintf(fname1,"f%d",z);
	  frho_ch_int[current_loop] = new TF1(fname1,rho_ch_int,0.,10.,2*ngaus+1);
	  char fname2[20];
	  sprintf(fname2,"f%d",z);
	  frho_rms[current_loop] = new TF1(fname2,rho_rms,0.,10.,2*ngaus+1);
	  char fname3[20];
	  sprintf(fname3,"f%d",z);
	  frho_ch[current_loop] = new TF1(fname3,rho_ch,0.,10.,2*ngaus+1);

	  //Set parameters for the various functions.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      fChFF[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      fChFF[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      //fChFF_Deriv[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      //fChFF_Deriv[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      frho_ch_int[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      frho_ch_int[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      frho_rms[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      frho_rms[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	      frho_ch[current_loop]->SetParameter(i,Qichmulti[current_loop][i]);
	      frho_ch[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	    }

	  fChFF[current_loop]->SetLineColorAlpha(2,transparency);
	  fChFF[current_loop]->SetLineWidth(linewidth);
	  frho_ch[current_loop]->SetLineColorAlpha(2,transparency);
	  frho_ch[current_loop]->SetLineWidth(linewidth);

	  fChFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  frho_ch[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.

	  if(Chi2[current_loop]<maxchi2 && Qimtot[current_loop]>(1-Qim_range) && Qimtot[current_loop]<(1+Qim_range))
	    {
	      if(Chi2[current_loop]<min_chi2)
		{
		  cout<<"*****Hello there!*****"<<endl;
		  min_chi2 = Chi2[current_loop];
		  min_chi2_fit = current_loop;
		}
	      if(first==0)
		{
		  cFch->cd();
		  fChFF[current_loop]->Draw("L");
		  if(target == 0)
		    {
		      fChFF[current_loop]->SetTitle("^{3}He Charge Form Factor");
		    }
		  if(target == 1)
		    {
		      fChFF[current_loop]->SetTitle("^{3}H Charge Form Factor");
		    }
		  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetTitle("|F_{ch}(q^{2})|");
		  fChFF[current_loop]->GetHistogram()->GetYaxis()->CenterTitle(true);
		  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
		  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
		  fChFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
		  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
		  fChFF[current_loop]->GetHistogram()->GetXaxis()->CenterTitle(true);
		  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
		  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
		  fChFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);
		  
		  cCh_den->cd();
		  frho_ch[current_loop]->Draw("L");
		  if(target == 0)
		    {
		      frho_ch[current_loop]->SetTitle("^{3}He Charge Density");
		    }
		  if(target == 1)
		    {
		      frho_ch[current_loop]->SetTitle("^{3}H Charge Density");
		    }
		  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetTitle("e/fm^{3}");
		  frho_ch[current_loop]->GetHistogram()->GetYaxis()->CenterTitle(true);
		  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
		  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
		  frho_ch[current_loop]->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
		  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetTitle("r (fm)");
		  frho_ch[current_loop]->GetHistogram()->GetXaxis()->CenterTitle(true);
		  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
		  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
		  frho_ch[current_loop]->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

		  first = 1; //No longer first plot.
		}
	      else
		{
		  cFch->cd();
		  fChFF[current_loop]->Draw("L SAME");

		  cCh_den->cd();
		  frho_ch[current_loop]->Draw("L SAME");
		  //frho_ch_int[current_loop]->Draw("L SAME");
		}
	      
	      for(Int_t i=0;i<ngaus;i++)
		{
		  hRi->Fill(Rmulti[current_loop][i]);
		  hRi_ind[i]->Fill(Rmulti[current_loop][i]);
		  Qi_ch[i]->Fill(Qichmulti[current_loop][i]);
		  Qi_m[i]->Fill(Qimmulti[current_loop][i]);
		}
	      for(Int_t i=0;i<(ngaus-1);i++)
		{
		  hRi_sep[i]->Fill(Rmulti[current_loop][i+1]-Rmulti[current_loop][i]);
		}
	      
	      //Calculate the charge radii.
	      double x0 = 0.0015;
	      ROOT::Math::Functor1D f1D(&ChFF_Deriv); //Fail fChFF[current_loop]. Fail ChFF_Q2.
	      ROOT::Math::RichardsonDerivator rd;
	      rd.SetFunction(f1D);
	      cout<<"First Derivative:   "<<rd.Derivative1(x0)<<"   -6*dFC(0)/dQ^2 = "<<-6*rd.Derivative1(x0)<<"   rms radius = "<<pow(-6*rd.Derivative1(x0),0.5)<<endl;

	      cout<<"Sum Qich = "<<Qichtot[current_loop]<<endl;
	      Qich_tot = Qich_tot + Qichtot[current_loop];   //Sum of Qich total values.
	      Qim_tot = Qim_tot + Qimtot[current_loop];      //Sum of Qim total values.
	      BIC_tot = BIC_tot + BIC[current_loop];         //Sum of BIC.
	      AIC_tot = AIC_tot + AIC[current_loop];         //Sum of AIC.
	      Chi2_tot = Chi2_tot + Chi2[current_loop];      //Sum of Chi^2 total values.
	      rChi2_tot = rChi2_tot + rChi2[current_loop];   //Sum of rChi^2 total values.

	      cout<<"Total Charge = "<<frho_ch_int[current_loop]->Integral(0.,10.)<<"e."<<endl;
	      cout<<"rms radius (integral method) = "<<pow( frho_rms[current_loop]->Integral(0.0,10.)/frho_ch_int[current_loop]->Integral(0.0,10.) ,0.5)<<endl; 

	      rms_deriv[z] = pow(-6*rd.Derivative1(x0),0.5);
	      rms_deriv_tot = rms_deriv[z] + rms_deriv_tot;
	      rms_int[z] = pow( frho_rms[current_loop]->Integral(0.0,10.)/frho_ch_int[current_loop]->Integral(0.0,10.) ,0.5);
	      rms_int_tot = rms_int[z] + rms_int_tot;

	      total_funcs++;
	      cout<<"Loop = "<<z<<endl;
	    }
	}
//cout<<"fChFF->Eval(0) = "<<fChFF[current_loop]->Eval(10.)<<endl;
//cout<<"fChFF->Eval(0) = "<<fChFF[current_loop]->Eval(0.0001)<<endl;
      //cout<<"loop before ++ = "<<current_loop<<endl;
      /*
      //Calculate the charge radii.
      double x0 = 0.0015;
      ROOT::Math::Functor1D f1D(&ChFF_Deriv); //Fail fChFF[current_loop]. Fail ChFF_Q2.
      ROOT::Math::RichardsonDerivator rd;
      rd.SetFunction(f1D);
      cout<<"First Derivative:   "<<rd.Derivative1(x0)<<"   -6*dFC(0)/dQ^2 = "<<-6*rd.Derivative1(x0)<<"   rms radius = "<<pow(-6*rd.Derivative1(x0),0.5)<<endl;
      rms[z] = pow(-6*rd.Derivative1(x0),0.5);
      rms_tot = rms[z] + rms_tot;
      */
      
      if(current_loop<(nlines-skip-1))
	{
	  current_loop++;
	}
      //cout<<"loop after ++ = "<<current_loop<<endl;
    }

  if(target == 0)
    {
      TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,yminFF,35.,1);
    }
  if(target == 1)
    {
      TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,yminFF,25.,1);
    }
  fChFF_Amroun->SetNpx(npdraw);
  fChFF_Amroun->SetLineColor(4);
  fChFF_Amroun->SetLineWidth(linewidth);
  cFch->cd();
  fChFF_Amroun->Draw("L SAME");
  auto ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
  if(target == 0)
    {
      ChFF_leg->AddEntry(fChFF[0],"New ^{3}He |F_{ch}(q^{2})| Fits","l");
      ChFF_leg->AddEntry("fChFF_Amroun","^{3}He |F_{ch}(q^{2})| Fit from Amroun et al.","l");
    }
  if(target == 1)
    {
      ChFF_leg->AddEntry(fChFF[0],"New ^{3}H |F_{ch}(q^{2})| Fits","l");
      ChFF_leg->AddEntry("fChFF_Amroun","^{3}H |F_{ch}(q^{2})| Fit from Amroun et al.","l");
    }
  ChFF_leg->Draw();

  //Now draw Amroun's error bands.
  //3H Fch upper error band.
  TGraph *gr1 = new TGraph (61, x_3H_Fch_up, y_3H_Fch_up);
  gr1->Draw("SAME C*");


  TCanvas* cFm=new TCanvas("cFm");
  cFm->SetGrid();
  cFm->SetLogy();
  cFm->SetTitle("Magnetic Form Factor");

  //for(current_loop=0; current_loop<1; current_loop++)//nlines-skip;current_loop++)
  //Reset current loop and first counter.
  current_loop = 0;
  first = 0;

  //Define array of TF1 to plot the various fits with.
  TF1 **fMFF = new TF1*[nfunc];

  for(Int_t z=0; z<nlines-skip; z++)
    {
      //cout<<"loop = "<<current_loop<<endl;

      if(current_loop==0)
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);
	  fMFF[current_loop] = new TF1(fname,MFF_Q2,yminFF,ymaxFF+54,2*ngaus);

	  //Set parameters for the various functions.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      fMFF[current_loop]->SetParameter(i,Qimmulti[current_loop][i]);
	      fMFF[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	    }

	  fMFF[current_loop]->SetLineColorAlpha(2,transparency);
	  fMFF[current_loop]->SetLineWidth(linewidth);

	  if(Chi2[current_loop]<maxchi2 && Qimtot[current_loop]>(1-Qim_range) && Qimtot[current_loop]<(1+Qim_range))
	    {
	      fMFF[current_loop]->Draw("L");
	      first = 1; //No longer first plot.
	      cout<<"Sum Qim = "<<Qimtot[current_loop]<<endl;
	    }

	  //fMFF[current_loop]->SetLineColor(3);
	  fMFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  if(target == 0)
	    {
	      fMFF[current_loop]->SetTitle("^{3}He Magnetic Form Factor");
	    }
	  if(target == 1)
	    {
	      fMFF[current_loop]->SetTitle("^{3}H Magnetic Form Factor");
	    }
	  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetTitle("|F_{m}(q^{2})|");
	  fMFF[current_loop]->GetHistogram()->GetYaxis()->CenterTitle(true);
	  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
	  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
	  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
	  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
	  fMFF[current_loop]->GetHistogram()->GetXaxis()->CenterTitle(true);
	  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
	  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
	  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);
	}
      else
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);

	  fMFF[current_loop] = new TF1(fname,MFF_Q2,yminFF,ymaxFF+54,2*ngaus);

	  //Set parameters for the various functions.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      fMFF[current_loop]->SetParameter(i,Qimmulti[current_loop][i]);
	      fMFF[current_loop]->SetParameter(i+ngaus,Rmulti[current_loop][i]);
	    }

	  fMFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  fMFF[current_loop]->SetLineColorAlpha(2,transparency);
	  fMFF[current_loop]->SetLineWidth(linewidth);

	  if(Chi2[current_loop]<maxchi2 && Qimtot[current_loop]>(1-Qim_range) && Qimtot[current_loop]<(1+Qim_range))
	    {
	      if(first==0)
		{
		  fMFF[current_loop]->Draw("L");
		  if(target == 0)
		    {
		      fMFF[current_loop]->SetTitle("^{3}He Magnetic Form Factor");
		    }
		  if(target == 1)
		    {
		      fMFF[current_loop]->SetTitle("^{3}H Magnetic Form Factor");
		    }
		  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetTitle("|F_{m}(q^{2})|");
		  fMFF[current_loop]->GetHistogram()->GetYaxis()->CenterTitle(true);
		  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
		  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
		  fMFF[current_loop]->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
		  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
		  fMFF[current_loop]->GetHistogram()->GetXaxis()->CenterTitle(true);
		  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
		  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
		  fMFF[current_loop]->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);
		  first = 1; //No longer first plot.
		}
	      else
		{
		  fMFF[current_loop]->Draw("L SAME");
		}
	      cout<<"Sum Qim = "<<Qimtot[current_loop]<<endl;
	    }
	}
      //fMFF[2]->Draw("L");
      //cout<<"fMFF->Eval(0) = "<<fMFF[current_loop]->Eval(0.0001)<<endl;
      //cout<<"fMFF->Eval(0) = "<<fMFF[current_loop]->Eval(0.01)<<endl;
      //cout<<"loop before ++ = "<<current_loop<<endl;
      if(current_loop<(nlines-skip-1))
	{
	  current_loop++;
	}
      //cout<<"loop after ++ = "<<current_loop<<endl;
    }
  /*
  fMFF[0]->Draw();
  for(Int_t i=1;i<10;i++)
    {
      fMFF[i]->Draw("lsame");
    }
  */



  /*
    TF1 *fMFF = new TF1("fMFF",MFF_Q2,0.,70.,1);
    fMFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
    fMFF->Draw("L");
    fMFF->SetTitle("^{3}He Magnetic Form Factor");
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
  */

  if(target == 0)
    {
      TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,35.,1);
    }
  if(target == 1)
    {
      TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,31.,1);
    }
  //cout<<fMFF_Amroun->Eval(30)<<endl;
  fMFF_Amroun->SetNpx(npdraw);
  fMFF_Amroun->SetLineColor(4);
  fMFF_Amroun->SetLineWidth(linewidth);
  //fMFF_Amroun->SetLineColorAlpha(4,0.1);
  fMFF_Amroun->Draw("L SAME");
  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  if(target == 0)
    {
      MFF_leg->AddEntry(fMFF[0],"New ^{3}He |F_{m}(q^{2})| Fits","l");
      MFF_leg->AddEntry("fMFF_Amroun","^{3}He |F_{m}(q^{2})| Fit from Amroun et al.","l");
    }
  if(target == 1)
    {
      MFF_leg->AddEntry(fMFF[0],"New ^{3}H |F_{m}(q^{2})| Fits","l");
      MFF_leg->AddEntry("fMFF_Amroun","^{3}H |F_{m}(q^{2})| Fit from Amroun et al.","l");
    }
  MFF_leg->Draw();
  
  //Plot distribution of charge radii derivative method.
  TCanvas* crms_deriv=new TCanvas("crms_deriv");
  crms_deriv->SetGrid();
  crms_deriv->SetTitle("Charge Radii (Derivative Method)");

  hrms_deriv = new TH1F("hrms_deriv", "Charge Radii Distribution (Derivative Method)", rms_deriv_divisions, rms_deriv_min, rms_deriv_max);
  for(Int_t i=0;i<(nlines-1);i++)
    {
      if(rms_deriv[i]!=0)
	{
	  hrms_deriv->Fill(rms_deriv[i]);
	}
      cout<<"rms_deriv["<<i<<"] = "<<rms_deriv[i]<<endl;
    }
  hrms_deriv->GetYaxis()->CenterTitle(true);
  hrms_deriv->GetYaxis()->SetLabelSize(0.05);
  hrms_deriv->GetYaxis()->SetTitleSize(0.06);
  hrms_deriv->GetYaxis()->SetTitleOffset(0.75);
  hrms_deriv->GetXaxis()->CenterTitle(true);
  hrms_deriv->GetXaxis()->SetLabelSize(0.05);
  hrms_deriv->GetXaxis()->SetTitleSize(0.06);
  hrms_deriv->GetXaxis()->SetTitleOffset(0.75);
  hrms_deriv->SetLineWidth(linewidth);
  hrms_deriv->GetXaxis()->SetTitle("r (fm)");
  hrms_deriv->Draw();

  TF1 *fgaus1 = new TF1("fgaus1",fit_gaus,rms_deriv_min,rms_deriv_max,3);
  fgaus1->SetParameter(0,200.);//20
  fgaus1->SetParameter(1,hrms_deriv->GetMean());
  fgaus1->SetParameter(2,hrms_deriv->GetStdDev());
  //fgaus1->SetParLimits(0,100.,200.);//3H 39.098
  //fgaus1->SetParLimits(1,2.,2.03);//3H 2.03474
  //fgaus1->SetParLimits(2,0.009,0.015);//3H -0.0212239
  hrms_deriv->Fit("fgaus1","R same M");
  cout<<"Gaussian Fit of Charge Radii (Derivative Method): Chi^2 = "<<fgaus1->GetChisquare()<<"   nDOF = "<<fgaus1->GetNDF()<<"   Fit Probablility = "<<fgaus1->GetProb()<<endl;
  cout<<"Gaussian Height = "<<fgaus1->GetParameter(0)<<"   Gaussian Mean = "<<fgaus1->GetParameter(1)<<"   Gaussian Standard Deviation = "<<fgaus1->GetParameter(2)<<endl;

  //Plot distribution of charge radii integral method.
  TCanvas* crms_int=new TCanvas("crms_int");
  crms_int->SetGrid();
  crms_int->SetTitle("Charge Radii (Integral Method)");

  hrms_int = new TH1F("hrms_int", "Charge Radii Distribution (Integral Method)", rms_deriv_divisions, rms_deriv_min, rms_deriv_max);
  for(Int_t i=0;i<(nlines-1);i++)
    {
      if(rms_deriv[i]!=0)
	{
	  hrms_int->Fill(rms_int[i]);
	}
    }
  hrms_int->GetYaxis()->CenterTitle(true);
  hrms_int->GetYaxis()->SetLabelSize(0.05);
  hrms_int->GetYaxis()->SetTitleSize(0.06);
  hrms_int->GetYaxis()->SetTitleOffset(0.75);
  hrms_int->GetXaxis()->CenterTitle(true);
  hrms_int->GetXaxis()->SetLabelSize(0.05);
  hrms_int->GetXaxis()->SetTitleSize(0.06);
  hrms_int->GetXaxis()->SetTitleOffset(0.75);
  hrms_int->SetLineWidth(linewidth);
  hrms_int->SetLineWidth(linewidth);
  hrms_int->GetXaxis()->SetTitle("r (fm)");
  hrms_int->Draw();

  TF1 *fgaus2 = new TF1("fgaus2",fit_gaus,rms_int_min,rms_int_max,3);
  fgaus2->SetParameter(0,1.);//(0,12.);
  fgaus2->SetParameter(1,hrms_int->GetMean());//(1,1.9);
  fgaus2->SetParameter(2,hrms_int->GetStdDev());     //Messes up fit.
  hrms_int->Fit("fgaus2","R same M");
  cout<<"Gaussian Fit of Charge Radii (Integral Method): Chi^2 = "<<fgaus2->GetChisquare()<<"   nDOF = "<<fgaus2->GetNDF()<<"   Fit Probablility = "<<fgaus2->GetProb()<<endl;
  cout<<"Gaussian Height = "<<fgaus2->GetParameter(0)<<"   Gaussian Mean = "<<fgaus2->GetParameter(1)<<"   Gaussian Standard Deviation = "<<fgaus2->GetParameter(2)<<endl;

  //Print total number of fits suviving the chi2 cut.
  cout<<"# of fits below Chi^2 cutoff = "<<total_funcs<<endl;
  cout<<"Average sum of Chi^2 = "<<Chi2_tot/total_funcs<<".   Average sum of rChi^2 = "<<rChi2_tot/total_funcs<<"."<<endl;
  cout<<"Average sum of Qich = "<<Qich_tot/total_funcs<<".   Average sum of Qim = "<<Qim_tot/total_funcs<<"."<<endl;
  cout<<"Average fit BIC = "<<BIC_tot/total_funcs<<".   Average fit AIC = "<<AIC_tot/total_funcs<<"."<<endl;
  cout<<"Average charge rms (derivative method) = "<<rms_deriv_tot/total_funcs<<endl;
  cout<<"Average charge rms (integral method) = "<<rms_int_tot/total_funcs<<endl;

  //Plot Ri on one histogram.
  TCanvas* cRi_all=new TCanvas("cR_all");
  cRi_all->SetGrid();
  cRi_all->SetTitle("Ri Distribution");
  hRi->Draw();
  hRi->SetLineWidth(2);

  //Plot individual Ri histograms.
  TCanvas* cRi_ind=new TCanvas("cRi_ind");
  cRi_ind->SetGrid();
  cRi_ind->SetTitle("Ri Distributions");
  cRi_ind->Divide(ngaus/2,2);

  TF1 **fgaus_ri = new TF1*[ngaus];

  for(Int_t i=0;i<ngaus;i++)
    {
      cRi_ind->cd(i+1);
      hRi_ind[i]->Draw();
      hRi_ind[i]->SetLineWidth(2);

      char fname5[20];
      sprintf(fname5,"fgaus_r%d",i);
      fgaus_ri[i] = new TF1(fname5,fit_gaus,gaus_min,gaus_max,3);

      fgaus_ri[i]->SetParameter(0,500.);
      fgaus_ri[i]->SetParameter(1,hRi_ind[i]->GetMean());
      fgaus_ri[i]->SetParameter(2,hRi_ind[i]->GetStdDev());
      hRi_ind[i]->Fit(Form("fgaus_r%d",i),"R same M Q 0");
    }

  //Plot Ri separation histograms.
  TCanvas* cRi_sep=new TCanvas("cRi_sep");
  cRi_sep->SetGrid();
  cRi_sep->SetTitle("Ri Distributions");
  cRi_sep->Divide(ngaus/2,2);

  for(Int_t i=0;i<(ngaus-1);i++)
    {
      cRi_sep->cd(i+1);
      hRi_sep[i]->Draw();
      hRi_sep[i]->SetLineWidth(2);
    }

  //Plot individual Qi_ch histograms.
  TCanvas* cQi_ch=new TCanvas("cQi_ch");
  cQi_ch->SetGrid();
  cQi_ch->SetTitle("Qi_ch Distributions");
  cQi_ch->Divide(ngaus/2,2);

  TF1 **fgaus_Qich = new TF1*[ngaus];

  for(Int_t i=0;i<ngaus;i++)
    {
      cQi_ch->cd(i+1);
      Qi_ch[i]->Draw();
      Qi_ch[i]->SetLineWidth(2);

      char fname6[20];
      sprintf(fname6,"fgaus_Q%d_ch",i);
      fgaus_Qich[i] = new TF1(fname6,fit_gaus,gaus_Qi_min,gaus_Qi_max,3);

      fgaus_Qich[i]->SetParameter(0,200.);
      fgaus_Qich[i]->SetParameter(1,Qi_ch[i]->GetMean());
      fgaus_Qich[i]->SetParameter(2,Qi_ch[i]->GetStdDev());
      Qi_ch[i]->Fit(Form("fgaus_Q%d_ch",i),"R same M Q 0");
    }

  //Plot individual Qi_m histograms.
  TCanvas* cQi_m=new TCanvas("cQi_m");
  cQi_m->SetGrid();
  cQi_m->SetTitle("Qi_m Distributions");
  cQi_m->Divide(ngaus/2,2);

  TF1 **fgaus_Qim = new TF1*[ngaus];

  for(Int_t i=0;i<ngaus;i++)
    {
      cQi_m->cd(i+1);
      Qi_m[i]->Draw();
      Qi_m[i]->SetLineWidth(2);

      char fname7[20];
      sprintf(fname7,"fgaus_Q%d_m",i);
      fgaus_Qim[i] = new TF1(fname7,fit_gaus,gaus_Qi_min,gaus_Qi_max,3);

      fgaus_Qim[i]->SetParameter(0,200.);
      fgaus_Qim[i]->SetParameter(1,Qi_m[i]->GetMean());
      //fgaus_Qim[i]->SetParLimits(1,0.,0.4);
      fgaus_Qim[i]->SetParameter(2,Qi_m[i]->GetStdDev());
      Qi_m[i]->Fit(Form("fgaus_Q%d_m",i),"R same M Q 0");
    }

  /*
  //Average Ri and Qi values. Doesn't seem to be very useful.
  for(Int_t i=0;i<ngaus;i++)
    {
      cout<<"R"<<i<<" mean = "<<fgaus_ri[i]->GetParameter(1)<<endl;
    }

  for(Int_t i=0;i<ngaus;i++)
    {
      cout<<"Q"<<i<<"ch mean = "<<fgaus_Qich[i]->GetParameter(1)<<endl;
    }

  for(Int_t i=0;i<ngaus;i++)
    {
      cout<<"Q"<<i<<"m mean = "<<fgaus_Qim[i]->GetParameter(1)<<endl;
    }
  */

  //Create an output file to store a single fit result in.
  std::ofstream output ("Representative_Fit.txt", std::ofstream::out);

  output<<"Fit#   Chi2   rChi2   BIC   AIC    Qichtot    Qimtot  R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch    Q0m    Q1m    Q2m    Q3m    Q4m    Q5m    Q6m    Q7m    Q8m    Q9m    Q10m    Q11m"<<endl;

  //Fill text file with fit#, fcn (chi2), Qichtot, Qimtot, Ri, Qich, Qim.
  output<<rep_chi2_fit<<" "<<Chi2[rep_chi2_fit]<<" "<<rChi2[rep_chi2_fit]<<" "<<BIC[rep_chi2_fit]<<" "<<AIC[rep_chi2_fit]<<" "<<Qichtot[rep_chi2_fit]<<" "<<Qimtot[rep_chi2_fit]<<" "<<Rmulti[rep_chi2_fit][0]<<" "<<Rmulti[rep_chi2_fit][1]<<" "<<Rmulti[rep_chi2_fit][2]<<" "<<Rmulti[rep_chi2_fit][3]<<" "<<Rmulti[rep_chi2_fit][4]<<" "<<Rmulti[rep_chi2_fit][5]<<" "<<Rmulti[rep_chi2_fit][6]<<" "<<Rmulti[rep_chi2_fit][7]<<" "<<Rmulti[rep_chi2_fit][8]<<" "<<Rmulti[rep_chi2_fit][9]<<" "<<Rmulti[rep_chi2_fit][10]<<" "<<Rmulti[rep_chi2_fit][11]<<" "<<Qichmulti[rep_chi2_fit][0]<<" "<<Qichmulti[rep_chi2_fit][1]<<" "<<Qichmulti[rep_chi2_fit][2]<<" "<<Qichmulti[rep_chi2_fit][3]<<" "<<Qichmulti[rep_chi2_fit][4]<<" "<<Qichmulti[rep_chi2_fit][5]<<" "<<Qichmulti[rep_chi2_fit][6]<<" "<<Qichmulti[rep_chi2_fit][7]<<" "<<Qichmulti[rep_chi2_fit][8]<<" "<<Qichmulti[rep_chi2_fit][9]<<" "<<Qichmulti[rep_chi2_fit][10]<<" "<<Qichmulti[rep_chi2_fit][11]<<" "<<Qimmulti[rep_chi2_fit][0]<<" "<<Qimmulti[rep_chi2_fit][1]<<" "<<Qimmulti[rep_chi2_fit][2]<<" "<<Qimmulti[rep_chi2_fit][3]<<" "<<Qimmulti[rep_chi2_fit][4]<<" "<<Qimmulti[rep_chi2_fit][5]<<" "<<Qimmulti[rep_chi2_fit][6]<<" "<<Qimmulti[rep_chi2_fit][7]<<" "<<Qimmulti[rep_chi2_fit][8]<<" "<<Qimmulti[rep_chi2_fit][9]<<" "<<Qimmulti[rep_chi2_fit][10]<<" "<<Qimmulti[rep_chi2_fit][11]<<endl;

  for(Int_t i=0;i<ngaus;i++)
    {
      cout<<"R"<<i<<" min Chi^2 = "<<Rmulti[min_chi2_fit][i]<<endl;
    }

  for(Int_t i=0;i<ngaus;i++)
    {
      cout<<"Q"<<i<<"ch min Chi^2 = "<<Qichmulti[min_chi2_fit][i]<<endl;
    }

  for(Int_t i=0;i<ngaus;i++)
    {
      cout<<"Q"<<i<<"m min Chi^2  = "<<Qimmulti[min_chi2_fit][i]<<endl;
    }

  cout<<"Min Chi^2 fit number = "<<min_chi2_fit<<". Min Chi^2 value = "<<min_chi2<<"."<<endl;
  /*
  //Test Poisson.
  TCanvas* cpoisson=new TCanvas("cpoisson");
  Qi_m[1]->Draw();
  Qi_m[1]->SetLineWidth(2);
  fpoisson = new TF1("fpoisson",poisson,-1.,1.,2);
  fpoisson->SetParameter(0,10.);
  fpoisson->SetParameter(1,0.1);
  fpoisson->Draw("same");
  //Qi_m[1]->Fit(Form("fpoisson",i),"R same M");
  */
}
