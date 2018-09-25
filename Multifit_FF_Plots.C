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

const Int_t nfunc = 1000;
Int_t loops = 1;
Int_t current_loop = 0;
const Int_t datapts = 259;//248
const Int_t size = 1000;//248
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
Int_t ngaus = 10;                        //Number of Gaussians used to fit data.
Int_t ngaus_Amroun = 12;                        //Number of Gaussians used to fit data from Amroun.
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
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
Double_t truncate = 100.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 1.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char* str[1000];                          //Variable to read lines of the data file.
Float_t Chi2[size],rChi2[size],BIC[size],AIC[size],Qichtot[size],Qimtot[size],R0[size],R1[size],R2[size],R3[size],R4[size],R5[size],R6[size],R7[size],R8[size],R9[size],R10[size],R11[size],Q0ch[size],Q1ch[size],Q2ch[size],Q3ch[size],Q4ch[size],Q5ch[size],Q6ch[size],Q7ch[size],Q8ch[size],Q9ch[size],Q10ch[size],Q11ch[size],Q0m[size],Q1m[size],Q2m[size],Q3m[size],Q4m[size],Q5m[size],Q6m[size],Q7m[size],Q8m[size],Q9m[size],Q10m[size],Q11m[size];
Float_t Chi2temp,rChi2temp,BICtemp,AICtemp,Qichtottemp,Qimtottemp,R0temp,R1temp,R2temp,R3temp,R4temp,R5temp,R6temp,R7temp,R8temp,R9temp,R10temp,R11temp,Q0chtemp,Q1chtemp,Q2chtemp,Q3chtemp,Q4chtemp,Q5chtemp,Q6chtemp,Q7chtemp,Q8chtemp,Q9chtemp,Q10chtemp,Q11chtemp,Q0mtemp,Q1mtemp,Q2mtemp,Q3mtemp,Q4mtemp,Q5mtemp,Q6mtemp,Q7mtemp,Q8mtemp,Q9mtemp,Q10mtemp,Q11mtemp;
Float_t theta[datapts];                     //Angle in degrees.
Float_t qeff[datapts];                      //q effective in fm^-1.
Float_t sigexp[datapts];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[datapts];
Float_t E0[datapts];
Float_t Q2[datapts];
Float_t Rmulti[size][12];
Float_t Qichmulti[size][12];
Float_t Qimmulti[size][12];

Int_t Amroun_pts = 57;
Int_t Collard_pts = 118;
Int_t Szlata_pts = 22;
Int_t Dunn_pts = 27;
Int_t Camsonne_pts = 16;
Int_t Nakagawa_pts = 5;
Int_t my_pts = 1;

Double_t R[12] = {0.1,0.7,1.3,2.,2.7,3.6,4.4,5.6,0.,0.,0.,0.};//7
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t Qich[12] = {0.0784469,0.247165,0.406019,0.120177,0.137968,4.57535E-11,0.0200847,2.63439E-9,0.,0.,0.,0.};//7
Double_t Qim[12] = {0.0770148,0.298502,0.282963,0.175066,0.0769078,0.0381075,0.0899692,0.0675,0.,0.,0.,0.};
Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
Double_t av[24] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t averr[24] = {};
Double_t Qicherr[12]={}; 
Double_t Qimerr[12]={};
//Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};
Double_t Chi2_FB[datapts]={};
Double_t residual_FB[datapts]={};
Double_t FBfit[datapts]={};

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
/*
//Plot Charge FF Fch(Q^2) fm^-2.
Double_t ChFF_Q2(Double_t *Q2, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
  //cout<<"Current Loop = "<<current_loop<<endl;
  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      sumchtemp = (Qichmulti[current_loop][i]/(1.0+2.0*pow(Rmulti[current_loop][i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*Rmulti[current_loop][i]) + (2.0*pow(Rmulti[current_loop][i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*Rmulti[current_loop][i])/(pow(Q2[0],0.5)*Rmulti[current_loop][i])) );
	
      fitch = fitch + sumchtemp;
    }

  fitch = fitch * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}
*/
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
/*
//Plot magnetic FF(Q^2) fm^-2.
Double_t MFF_Q2(Double_t *Q2, Double_t *par)
{
  Double_t fitm = 0.;
  Double_t summtemp = 0.;
       
  //Define SOG for magnetic FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
      fitm = fitm + summtemp;
    }
       
  fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitm = fabs(fitm);
  return fitm;
}
*/

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
    }

  fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
  fitm = fabs(fitm);
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

void Multifit_FF_Plots() 
{
  FILE *fp;
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Chi2.txt","r");
  //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_BS_300_Ri_Chi2.txt","r");
  fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Save_Ri_Fits_180_9_25_2018.txt","r");
    
  //Read in data.
  while (1) {
    //Skips the first 5 lines of the file. 
    if (nlines < skip)
      {
	fgets(str,1000,fp);
	nlines++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp, &R3temp);
	ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &Q0chtemp, &Q1chtemp);
	ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp, &Q11chtemp);
	ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp, &Q8mtemp, &Q9mtemp);
	ncols = ncols + fscanf(fp,"%f %f", &Q10mtemp, &Q11mtemp);
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

	nlines++;
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

  //Now plot all of the curves on one canvas to form an error band.

  TCanvas* cFch=new TCanvas("cFch");
  cFch->SetGrid();
  cFch->SetLogy();
  cFch->SetTitle("Charge Form Factor");

  //Define array of TF1 to plot the various fits with.
  TF1 **fChFF = new TF1*[nfunc];

  //for(current_loop=0; current_loop<1; current_loop++)//nlines-skip;current_loop++)
  for(Int_t z=0; z<nlines-skip; z++)
    {
      //cout<<"loop = "<<current_loop<<endl;

      if(current_loop==0)
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);
	  fChFF[current_loop] = new TF1(fname,ChFF_Q2,yminFF,ymaxFF+54,20);

	  fChFF[current_loop]->SetParameter(0,Qichmulti[current_loop][0]);
	  fChFF[current_loop]->SetParameter(1,Qichmulti[current_loop][1]);
	  fChFF[current_loop]->SetParameter(2,Qichmulti[current_loop][2]);
	  fChFF[current_loop]->SetParameter(3,Qichmulti[current_loop][3]);
	  fChFF[current_loop]->SetParameter(4,Qichmulti[current_loop][4]);
	  fChFF[current_loop]->SetParameter(5,Qichmulti[current_loop][5]);
	  fChFF[current_loop]->SetParameter(6,Qichmulti[current_loop][6]);
	  fChFF[current_loop]->SetParameter(7,Qichmulti[current_loop][7]);
	  fChFF[current_loop]->SetParameter(8,Qichmulti[current_loop][8]);
	  fChFF[current_loop]->SetParameter(9,Qichmulti[current_loop][9]);
	  fChFF[current_loop]->SetParameter(10,Rmulti[current_loop][0]);
	  fChFF[current_loop]->SetParameter(11,Rmulti[current_loop][1]);
	  fChFF[current_loop]->SetParameter(12,Rmulti[current_loop][2]);
	  fChFF[current_loop]->SetParameter(13,Rmulti[current_loop][3]);
	  fChFF[current_loop]->SetParameter(14,Rmulti[current_loop][4]);
	  fChFF[current_loop]->SetParameter(15,Rmulti[current_loop][5]);
	  fChFF[current_loop]->SetParameter(16,Rmulti[current_loop][6]);
	  fChFF[current_loop]->SetParameter(17,Rmulti[current_loop][7]);
	  fChFF[current_loop]->SetParameter(18,Rmulti[current_loop][8]);
	  fChFF[current_loop]->SetParameter(19,Rmulti[current_loop][9]);

	  fChFF[current_loop]->Draw("L");
	  fChFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  fChFF[current_loop]->SetTitle("^{3}He Charge Form Factor");
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
	}
      else
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);
	  fChFF[current_loop] = new TF1(fname,ChFF_Q2,yminFF,ymaxFF+54,20);

	  fChFF[current_loop]->SetParameter(0,Qichmulti[current_loop][0]);
	  fChFF[current_loop]->SetParameter(1,Qichmulti[current_loop][1]);
	  fChFF[current_loop]->SetParameter(2,Qichmulti[current_loop][2]);
	  fChFF[current_loop]->SetParameter(3,Qichmulti[current_loop][3]);
	  fChFF[current_loop]->SetParameter(4,Qichmulti[current_loop][4]);
	  fChFF[current_loop]->SetParameter(5,Qichmulti[current_loop][5]);
	  fChFF[current_loop]->SetParameter(6,Qichmulti[current_loop][6]);
	  fChFF[current_loop]->SetParameter(7,Qichmulti[current_loop][7]);
	  fChFF[current_loop]->SetParameter(8,Qichmulti[current_loop][8]);
	  fChFF[current_loop]->SetParameter(9,Qichmulti[current_loop][9]);
	  fChFF[current_loop]->SetParameter(10,Rmulti[current_loop][0]);
	  fChFF[current_loop]->SetParameter(11,Rmulti[current_loop][1]);
	  fChFF[current_loop]->SetParameter(12,Rmulti[current_loop][2]);
	  fChFF[current_loop]->SetParameter(13,Rmulti[current_loop][3]);
	  fChFF[current_loop]->SetParameter(14,Rmulti[current_loop][4]);
	  fChFF[current_loop]->SetParameter(15,Rmulti[current_loop][5]);
	  fChFF[current_loop]->SetParameter(16,Rmulti[current_loop][6]);
	  fChFF[current_loop]->SetParameter(17,Rmulti[current_loop][7]);
	  fChFF[current_loop]->SetParameter(18,Rmulti[current_loop][8]);
	  fChFF[current_loop]->SetParameter(19,Rmulti[current_loop][9]);

	  fChFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  fChFF[current_loop]->Draw("L SAME");
	}
      cout<<"fChFF->Eval(0) = "<<fChFF[current_loop]->Eval(0.0001)<<endl;
      //cout<<"loop before ++ = "<<current_loop<<endl;
      if(current_loop<(nlines-skip-1))
	{
	  current_loop++;
	}
      //cout<<"loop after ++ = "<<current_loop<<endl;
    }

  TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,yminFF,ymaxFF+54,1);
  fChFF_Amroun->SetNpx(npdraw);
  fChFF_Amroun->SetLineColor(4);
  fChFF_Amroun->Draw("L SAME");
  auto ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
  ChFF_leg->AddEntry(fChFF[0],"New ^{3}He |F_{ch}(q^{2})| Fit","l");
  ChFF_leg->AddEntry("fChFF_Amroun","^{3}He |F_{ch}(q^{2})| Fit from Amroun et al. [4]","l");
  ChFF_leg->Draw();

  TCanvas* cFm=new TCanvas("cFm");
  cFm->SetGrid();
  cFm->SetLogy();
  cFm->SetTitle("Magnetic Form Factor");

  //for(current_loop=0; current_loop<1; current_loop++)//nlines-skip;current_loop++)
  //Reset current loop.
  current_loop = 0;

  //Define array of TF1 to plot the various fits with.
  TF1 **fMFF = new TF1*[nfunc];

  for(Int_t z=0; z<nlines-skip; z++)
    {
      //cout<<"loop = "<<current_loop<<endl;

      if(current_loop==0)
	{
	  char fname[20];
	  sprintf(fname,"f%d",z);
	  fMFF[current_loop] = new TF1(fname,MFF_Q2,yminFF,ymaxFF+54,20);

	  fMFF[current_loop]->SetParameter(0,Qimmulti[current_loop][0]);
	  fMFF[current_loop]->SetParameter(1,Qimmulti[current_loop][1]);
	  fMFF[current_loop]->SetParameter(2,Qimmulti[current_loop][2]);
	  fMFF[current_loop]->SetParameter(3,Qimmulti[current_loop][3]);
	  fMFF[current_loop]->SetParameter(4,Qimmulti[current_loop][4]);
	  fMFF[current_loop]->SetParameter(5,Qimmulti[current_loop][5]);
	  fMFF[current_loop]->SetParameter(6,Qimmulti[current_loop][6]);
	  fMFF[current_loop]->SetParameter(7,Qimmulti[current_loop][7]);
	  fMFF[current_loop]->SetParameter(8,Qimmulti[current_loop][8]);
	  fMFF[current_loop]->SetParameter(9,Qimmulti[current_loop][9]);
	  fMFF[current_loop]->SetParameter(10,Rmulti[current_loop][0]);
	  fMFF[current_loop]->SetParameter(11,Rmulti[current_loop][1]);
	  fMFF[current_loop]->SetParameter(12,Rmulti[current_loop][2]);
	  fMFF[current_loop]->SetParameter(13,Rmulti[current_loop][3]);
	  fMFF[current_loop]->SetParameter(14,Rmulti[current_loop][4]);
	  fMFF[current_loop]->SetParameter(15,Rmulti[current_loop][5]);
	  fMFF[current_loop]->SetParameter(16,Rmulti[current_loop][6]);
	  fMFF[current_loop]->SetParameter(17,Rmulti[current_loop][7]);
	  fMFF[current_loop]->SetParameter(18,Rmulti[current_loop][8]);
	  fMFF[current_loop]->SetParameter(19,Rmulti[current_loop][9]);

	  fMFF[current_loop]->Draw("L");

	  //fMFF[current_loop]->SetLineColor(3);
	  fMFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  fMFF[current_loop]->SetTitle("^{3}He Magnetic Form Factor");
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

	  fMFF[current_loop] = new TF1(fname,MFF_Q2,yminFF,ymaxFF+54,20);
	  fMFF[current_loop]->SetParameter(0,Qimmulti[current_loop][0]);
	  fMFF[current_loop]->SetParameter(1,Qimmulti[current_loop][1]);
	  fMFF[current_loop]->SetParameter(2,Qimmulti[current_loop][2]);
	  fMFF[current_loop]->SetParameter(3,Qimmulti[current_loop][3]);
	  fMFF[current_loop]->SetParameter(4,Qimmulti[current_loop][4]);
	  fMFF[current_loop]->SetParameter(5,Qimmulti[current_loop][5]);
	  fMFF[current_loop]->SetParameter(6,Qimmulti[current_loop][6]);
	  fMFF[current_loop]->SetParameter(7,Qimmulti[current_loop][7]);
	  fMFF[current_loop]->SetParameter(8,Qimmulti[current_loop][8]);
	  fMFF[current_loop]->SetParameter(9,Qimmulti[current_loop][9]);
	  fMFF[current_loop]->SetParameter(10,Rmulti[current_loop][0]);
	  fMFF[current_loop]->SetParameter(11,Rmulti[current_loop][1]);
	  fMFF[current_loop]->SetParameter(12,Rmulti[current_loop][2]);
	  fMFF[current_loop]->SetParameter(13,Rmulti[current_loop][3]);
	  fMFF[current_loop]->SetParameter(14,Rmulti[current_loop][4]);
	  fMFF[current_loop]->SetParameter(15,Rmulti[current_loop][5]);
	  fMFF[current_loop]->SetParameter(16,Rmulti[current_loop][6]);
	  fMFF[current_loop]->SetParameter(17,Rmulti[current_loop][7]);
	  fMFF[current_loop]->SetParameter(18,Rmulti[current_loop][8]);
	  fMFF[current_loop]->SetParameter(19,Rmulti[current_loop][9]);
	  fMFF[current_loop]->SetNpx(npdraw);   //Sets number of points to use when drawing the function.
	  fMFF[current_loop]->Draw("L SAME");
	}
      //fMFF[2]->Draw("L");
      cout<<"fMFF->Eval(0) = "<<fMFF[current_loop]->Eval(0.0001)<<endl;
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

  TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,ymaxFF+54,1);
  //cout<<fMFF_Amroun->Eval(30)<<endl;
  fMFF_Amroun->SetNpx(npdraw);
  fMFF_Amroun->SetLineColor(4);
  fMFF_Amroun->Draw("L same");
  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  MFF_leg->AddEntry(fMFF[0],"New ^{3}He |F_{m}(q^{2})| Fit","l");
  MFF_leg->AddEntry("fMFF_Amroun","^{3}He |F_{m}(q^{2})| Fit from Amroun et al. [4]","l");
  MFF_leg->Draw();

}
