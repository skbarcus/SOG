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

#include <vector>  //Needed to use C++ vectors.
//#include <algorithm>
//#include <iostream>
//#include <utility>

Double_t pi = 3.141592654;
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t e = 1.60217662E-19;             //Electron charge C.
Double_t alpha = 0.0072973525664;//1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.

Int_t Rep_Fit = 30;//122    //Fit# from Multifit requires -1. From Plot_FFs or Chi2_Sorted Fit# is just the number.
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
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
Int_t ngaus_Amroun = 8;                        //Number of Gaussians used to fit data from Amroun.
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z_3He = 2.;                         //Atomic number 3He.
Double_t Z_3H = 1.;                         //Atomic number 3H.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
Double_t gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
Double_t theta = 0.;//21.04;
Double_t theta_cor = 0.;                //Theta that corrects for the Q^2eff adjustment. Basically when we plot the XS and FFs the Q2[0] is really Q^2eff if we don't use this theta_cor. This variable is for the slightly smaller theta representing the real scattering angle.
Double_t bin_min_Q2 = 30.5589;
Double_t bin_max_Q2 = 41.247;
Double_t E0 = 3.356;                    //Initial e- energy GeV.
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
//Float_t theta[1000];                     //Angle in degrees.
Float_t qeff[1000];                      //q effective in fm^-1.
Float_t sigexp[1000];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[1000];
//Float_t E0[1000];
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

//Double_t R[12] = {0.1,0.6,0.7,1.3,1.4,2.,2.7,3.4,4.3,5.2,0.,0.};//51
//Double_t R[12] = {0.0001,0.6,0.7,1.2,1.3,1.9,2.6,3.3,4.2,4.4,0.,0.};//48
//Double_t R[12] = {0.3,0.6,1.,1.4,1.8,2.5,3.2,4.1,4.2,0.,0.,0.};//42
//Double_t R[12] = {0.3,0.6,1.,1.3,1.7,2.4,3.1,3.8,4.4,0.,0.,0.};//40
//Double_t R[12] = {0.3,0.7,1.1,1.2,1.8,2.5,3.2,4.,4.9,0.,0.,0.};//38
//Double_t R[12] = {0.2,0.5,1.,1.5,2.,2.6,3.5,4.,5.,0.,0.,0.};//34
//Double_t R[12] = {0.2,0.7,1.3,1.5,2.1,2.8,3.6,4.2,5.2,0.,0.,0.};//32
//Double_t R[12] = {0.2,0.6,1.1,1.4,1.9,2.6,3.5,4.1,5.1,0.,0.,0.};//31
//Double_t R[12] = {0.2,0.7,1.3,1.5,2.1,2.8,3.6,4.1,5.1,0.,0.,0.};//28
//Double_t R[12] = {0.3,0.8,1.1,1.5,2.2,2.8,3.4,4.4,4.8,0.,0.,0.};//27
//Double_t R[12] = {0.3,0.8,1.3,1.6,2.1,2.7,3.3,4.1,4.9,0.,0.,0.};//25
//Double_t R[12] = {0.2,0.7,1.3,1.7,2.2,2.8,3.5,4.5,4.4,0.,0.,0.};//24
//Double_t R[12] = {0.3,0.8,1.1,1.4,1.9,2.5,3.3,4.1,5.1,0.,0.,0.};//23
//Double_t R[12] = {0.3,0.8,1.3,1.7,2.3,3.,3.9,4.9,0.,0.,0.,0.};//18
//Double_t R[12] = {0.2,0.7,1.3,1.9,2.6,3.4,4.3,5.1,0.,0.,0.,0.};//11
//Double_t R[12] = {0.1,0.7,1.3,2.,2.7,3.6,4.4,5.6,0.,0.,0.,0.};//7
//Double_t R[12] = {0.1,0.7,1.3,1.9,2.7,3.2,4.5,5.5,0.,0.,0.,0.};//5
//Double_t R[12] = {0.1,0.6,0.8,1.3,1.6,2.2,3.,3.8,4.5,5.6,0.,0.};//4.
//Double_t R[12] = {0.1,0.6,0.9,1.4,1.8,2.2,3,3.8,4.6,4.5,0.,0.};//9/12/18 3.
//Double_t R[12] = {0.1,0.6,1.,1.5,2.1,2.4,3.,3.7,4.4,4.7,0.,0.};//9/12/18 pretty good 2.
//Double_t R[12] = {0.1,0.6,1.,1.5,2.1,2.4,3.1,3.9,4.5,4.8,0.,0.};//9/12/18 pretty good 1.


//Double_t R[12] = {0.190672,0.693568,1.06127,1.5556,1.97331,2.42858,3.13643,3.87028,4.60607,5.52267};//Avg. Ri 357 fits x^2<505.
//Double_t R[12] = {0.3,0.7,0.9,1.4,1.7,2.2,2.9,3.6,4.3,4.9};//Min Chi^2 of 357 (#556 combined). 

Double_t R[12] = {0.3, 0.7, 0.9, 1.1, 1.5, 1.6, 2.2, 2.7, 3.3, 4.2, 4.3, 4.8};  //My fit 3He
Double_t R_Amroun[12] = {0.3,0.8,1.4,1.9,2.5,3.3,4.1,4.8}; //My 3H Fit


Double_t Qich[12] = {0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.98962e-12};//My 3He fit
Double_t Qim[12] = {0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.77602e-12,0.0240927,8.94934e-12};//My 3He fit

Double_t Qich_Amroun[12] = {0.151488,0.348372,0.29635,0.0978631,0.121983,0.0242654,0.049329,4.40751e-11};//My 3H fit
Double_t Qim_Amroun[12] = {0.190646,0.301416,0.318972,0.159433,0.173933,0.106361,0.0665564,0.0148866};//My 3H fit
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

const Int_t size = 1000;
Float_t Fit[size], Chi2[size],rChi2[size],BIC[size],AIC[size],Qichtot[size],Qimtot[size],R0[size],R1[size],R2[size],R3[size],R4[size],R5[size],R6[size],R7[size],R8[size],R9[size],R10[size],R11[size],Q0ch[size],Q1ch[size],Q2ch[size],Q3ch[size],Q4ch[size],Q5ch[size],Q6ch[size],Q7ch[size],Q8ch[size],Q9ch[size],Q10ch[size],Q11ch[size],Q0m[size],Q1m[size],Q2m[size],Q3m[size],Q4m[size],Q5m[size],Q6m[size],Q7m[size],Q8m[size],Q9m[size],Q10m[size],Q11m[size];
Float_t Fittemp,Chi2temp,rChi2temp,BICtemp,AICtemp,Qichtottemp,Qimtottemp,R0temp,R1temp,R2temp,R3temp,R4temp,R5temp,R6temp,R7temp,R8temp,R9temp,R10temp,R11temp,Q0chtemp,Q1chtemp,Q2chtemp,Q3chtemp,Q4chtemp,Q5chtemp,Q6chtemp,Q7chtemp,Q8chtemp,Q9chtemp,Q10chtemp,Q11chtemp,Q0mtemp,Q1mtemp,Q2mtemp,Q3mtemp,Q4mtemp,Q5mtemp,Q6mtemp,Q7mtemp,Q8mtemp,Q9mtemp,Q10mtemp,Q11mtemp;
Float_t Rmulti[size][12];
Float_t Qichmulti[size][12];
Float_t Qimmulti[size][12];

//Adding FF world data for plots.
const Int_t size1 = 100;
Int_t nlines1,nlines2,nlines3,nlines4,nlines5,nlines6,nlines7,nlines8,nlines9,nlines10,nlines11,nlines12,nlines13,nlines14,nlines15,nlines16;
Int_t skip1 = 0,skip2 = 0,skip3 = 0,skip4 = 0,skip5 = 0,skip6 = 0,skip7 = 0,skip8 = 0,skip9 = 0,skip10 = 0,skip11 = 0,skip12 = 0,skip13 = 0,skip14 = 0,skip15 = 0,skip16 = 0;
char* str1[1000],str2[1000],str3[1000],str4[1000],str5[1000],str6[1000],str7[1000],str8[1000],str9[1000],str10[1000],str11[1000],str12[1000],str13[1000],str14[1000],str15[1000],str16[1000];
Float_t q2_fch_cam16[size1],fch_cam16[size1],dq2_fch_cam16[size1],dfch_cam16[size1],q2_fm_cam16[size1],fm_cam16[size1],dq2_fm_cam16[size1],dfm_cam16[size1];
Float_t q2_fch_cam16_temp,fch_cam16_temp,dfch_cam16_temp,q2_fm_cam16_temp,fm_cam16_temp,dfm_cam16_temp;
Float_t q2_col65[size1],fch_col65[size1],dq2_col65[size1],dfch_col65[size1],fm_col65[size1],dfm_col65[size1];
Float_t q2_col65_temp,fch_col65_temp,dfch_col65_temp,fm_col65_temp,dfm_col65_temp;
Float_t q2_sza77[size1],fch_sza77[size1],dq2_sza77[size1],dfch_sza77[size1];
Float_t q2_sza77_temp,fch_sza77_temp,dfch_sza77_temp;
Float_t q2_arn78[size1],fch_arn78[size1],dq2_arn78[size1],dfch_arn78[size1];
Float_t q2_arn78_temp,fch_arn78_temp,dfch_arn78_temp;
Float_t q2_dun83[size1],fch_dun83[size1],dq2_dun83[size1],dfch_dun83[size1];
Float_t q2_dun83_temp,fch_dun83_temp,dfch_dun83_temp;
Float_t q2_cav82[size1],fm_cav82[size1],dq2_cav82[size1],dfm_cav82[size1];
Float_t q2_cav82_temp,fm_cav82_temp,dfm_cav82_temp;
Float_t q2_nak01[size1],fm_nak01[size1],dq2_nak01[size1],dfm_nak01[size1];
Float_t q2_nak01_temp,fm_nak01_temp,dfm_nak01_temp;

void Compare_FFs() 
{
  //Read in parameters for the representative fit.
  FILE *fp;
  fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Combined_Ri_Fits.txt","r");

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
	ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp, &R3temp);
	ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &Q0chtemp, &Q1chtemp);
	ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp, &Q11chtemp);
	ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp, &Q8mtemp, &Q9mtemp);
	ncols = ncols + fscanf(fp,"%f %f", &Q10mtemp, &Q11mtemp);

	//Only if using representative fit text file.
	/*
	  ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Fittemp, &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R3temp, &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &Q0chtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q1chtemp, &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q11chtemp, &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp, &Q8mtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f", &Q9mtemp, &Q10mtemp, &Q11mtemp);
	*/
	  //cout<<"ncols = "<<ncols<<endl;
	  if (ncols < 0) break;    
	  
	  //Fit[nlines-skip] = Fittemp;               //Only if using representative fit text file. 
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


  //Open and read in files for Camsonne 2016 Fch.
  FILE *fp1;
  fp1 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Camsonne2016_Fch.txt","r");

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
	  ncols = fscanf(fp1,"%f %f %f", &q2_fch_cam16_temp, &fch_cam16_temp, &dfch_cam16_temp);
	  if (ncols < 0) break;    
  
	  q2_fch_cam16[nlines1-skip1] = q2_fch_cam16_temp;
	  fch_cam16[nlines1-skip1] = fch_cam16_temp;
	  dfch_cam16[nlines1-skip1] = dfch_cam16_temp;
	  dq2_fch_cam16[nlines1-skip1] = 0.;
	  //cout<<"q2_fch_cam16[["<<nlines1-skip1<<"] = "<<q2_fch_cam16[nlines1-skip1]<<"   fch_cam16["<<nlines1-skip1<<"] = "<<fch_cam16[nlines1-skip1]<<"   dfch_cam16["<<nlines1-skip1<<"] = "<<dfch_cam16[nlines1-skip1]<<endl;
	  nlines1++;
	}
    }
  fclose(fp1);
  cout<<"nlines1 = "<<nlines1<<endl;

  //Open and read in files for Camsonne 2016 Fm.
  FILE *fp2;
  fp2 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Camsonne2016_Fm.txt","r");

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
	  ncols = fscanf(fp2,"%f %f %f", &q2_fm_cam16_temp, &fm_cam16_temp, &dfm_cam16_temp);
	  if (ncols < 0) break;    
  
	  q2_fm_cam16[nlines2-skip2] = q2_fm_cam16_temp;
	  fm_cam16[nlines2-skip2] = fm_cam16_temp;
	  dfm_cam16[nlines2-skip2] = dfm_cam16_temp;
	  dq2_fm_cam16[nlines2-skip2] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines2++;
	}
    }
  fclose(fp2);
  cout<<"nlines2 = "<<nlines2<<endl;

  //Open and read in files for Collard 1965.
  FILE *fp3;
  fp3 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Collard1965_3He.txt","r");

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
	  ncols = fscanf(fp3,"%f %f %f %f %f", &q2_col65_temp, &fch_col65_temp, &dfch_col65_temp, &fm_col65_temp, &dfm_col65_temp);
	  if (ncols < 0) break;    
  
	  q2_col65[nlines3-skip3] = q2_col65_temp;
	  fch_col65[nlines3-skip3] = fch_col65_temp;
	  dfch_col65[nlines3-skip3] = dfch_col65_temp;
	  fm_col65[nlines3-skip3] = fm_col65_temp;
	  dfm_col65[nlines3-skip3] = dfm_col65_temp;
	  dq2_col65[nlines3-skip3] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines3++;
	}
    }
  fclose(fp3);
  cout<<"nlines3 = "<<nlines3<<endl;

  //Open and read in files for Szalata Fch 1977.
  FILE *fp4;
  fp4 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Szalata1977_Fch.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines4 < skip4)
	{
	  fgets(str4,1000,fp4);
	  nlines4++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp4,"%f %f %f", &q2_sza77_temp, &fch_sza77_temp, &dfch_sza77_temp);
	  if (ncols < 0) break;    
  
	  q2_sza77[nlines4-skip4] = q2_sza77_temp;
	  fch_sza77[nlines4-skip4] = fch_sza77_temp;
	  dfch_sza77[nlines4-skip4] = dfch_sza77_temp;
	  dq2_sza77[nlines4-skip4] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines4++;
	}
    }
  fclose(fp4);
  cout<<"nlines4 = "<<nlines4<<endl;

  //Open and read in files for Arnold Fch 1978.
  FILE *fp5;
  fp5 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Arnold1978_Fch.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines5 < skip5)
	{
	  fgets(str5,1000,fp5);
	  nlines5++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp5,"%f %f %f", &q2_arn78_temp, &fch_arn78_temp, &dfch_arn78_temp);
	  if (ncols < 0) break;    
  
	  q2_arn78[nlines5-skip5] = q2_arn78_temp;
	  fch_arn78[nlines5-skip5] = fch_arn78_temp;
	  dfch_arn78[nlines5-skip5] = dfch_arn78_temp;
	  dq2_arn78[nlines5-skip5] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines5++;
	}
    }
  fclose(fp5);
  cout<<"nlines5 = "<<nlines5<<endl;

  //Open and read in files for Dunn Fch 1983.
  FILE *fp6;
  fp6 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Dunn1983_Fch.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines6 < skip6)
	{
	  fgets(str6,1000,fp6);
	  nlines6++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp6,"%f %f %f", &q2_dun83_temp, &fch_dun83_temp, &dfch_dun83_temp);
	  if (ncols < 0) break;    
  
	  q2_dun83[nlines6-skip6] = q2_dun83_temp;
	  fch_dun83[nlines6-skip6] = fch_dun83_temp;
	  dfch_dun83[nlines6-skip6] = dfch_dun83_temp;
	  dq2_dun83[nlines6-skip6] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines6++;
	}
    }
  fclose(fp6);
  cout<<"nlines6 = "<<nlines6<<endl;

  //Open and read in files for Cavedon Fm 1982.
  FILE *fp7;
  fp7 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Cavedon1982_Fm.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines7 < skip7)
	{
	  fgets(str7,1000,fp7);
	  nlines7++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp7,"%f %f %f", &q2_cav82_temp, &fm_cav82_temp, &dfm_cav82_temp);
	  if (ncols < 0) break;    
  
	  q2_cav82[nlines7-skip7] = q2_cav82_temp;
	  fm_cav82[nlines7-skip7] = fm_cav82_temp;
	  dfm_cav82[nlines7-skip7] = dfm_cav82_temp;
	  dq2_cav82[nlines7-skip7] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines7++;
	}
    }
  fclose(fp7);
  cout<<"nlines7 = "<<nlines7<<endl;

  //Open and read in files for Nakagawa Fm 2001.
  FILE *fp8;
  fp8 = fopen("/home/skbarcus/Tritium/Analysis/SOG/Nakagawa2001_Fm.txt","r");

  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines8 < skip8)
	{
	  fgets(str8,1000,fp8);
	  nlines8++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp8,"%f %f %f", &q2_nak01_temp, &fm_nak01_temp, &dfm_nak01_temp);
	  if (ncols < 0) break;    
  
	  q2_nak01[nlines8-skip8] = q2_nak01_temp;
	  fm_nak01[nlines8-skip8] = fm_nak01_temp;
	  dfm_nak01[nlines8-skip8] = dfm_nak01_temp;
	  dq2_nak01[nlines8-skip8] = 0.;
	  //cout<<"q2_fm_cam16[["<<nlines2-skip2<<"] = "<<q2_fm_cam16[nlines2-skip2]<<"   fm_cam16["<<nlines2-skip2<<"] = "<<fm_cam16[nlines2-skip2]<<"   dfm_cam16["<<nlines2-skip2<<"] = "<<dfm_cam16[nlines2-skip2]<<endl;
	  nlines8++;
	}
    }
  fclose(fp8);
  cout<<"nlines8 = "<<nlines8<<endl;

  //Set SOG parameters to the representative fit's parameters.
  for(Int_t i=0;i<ngaus;i++)
    {
      /*
      R[i] = Rmulti[Rep_Fit][i];
      Qich[i] = Qichmulti[Rep_Fit][i];
      Qim[i] = Qimmulti[Rep_Fit][i];
      */
    }

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
	//cout<<"Qich["<<i<<"] = "<<Qich[i]<<endl;
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

  TF1 *fChFF = new TF1("fChFF",ChFF_Q2,0.,35.,1);
  fChFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
  fChFF->Draw("L");
  cFch->SetTitle("Charge Form Factor");
  //fChFF->SetTitle("C12 Charge Form Factor","#Q^2 (#fm^-2)","#F_{Ch}(q)");
  fChFF->SetTitle("^{3}He Charge Form Factor");
  fChFF->GetHistogram()->GetYaxis()->SetTitle("|F_{ch}(Q^{2})|");
  fChFF->GetHistogram()->GetYaxis()->CenterTitle(true);
  fChFF->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  fChFF->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fChFF->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fChFF->GetHistogram()->GetXaxis()->SetTitle("Q^{2} (fm^{-2})");
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

  TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,0.,35.,1);
  //cout<<fChFF_Amroun->Eval(35)<<endl;
  fChFF_Amroun->SetNpx(npdraw);
  fChFF_Amroun->SetLineColor(4);
  fChFF_Amroun->Draw("L same");

  /*
  TGraphErrors *gr1 = new TGraphErrors (nlines1, q2_fch_cam16, fch_cam16, dq2_fch_cam16, dfch_cam16); 
  gr1->SetMarkerColor(1);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);
  gr1->Draw("same p");

  TGraphErrors *gr3 = new TGraphErrors (nlines3, q2_col65, fch_col65, dq2_col65, dfch_col65); 
  gr3->SetMarkerColor(4);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(1);
  gr3->Draw("same p");

  TGraphErrors *gr5 = new TGraphErrors (nlines4, q2_sza77, fch_sza77, dq2_sza77, dfch_sza77); 
  gr5->SetMarkerColor(3);
  gr5->SetMarkerStyle(20);
  gr5->SetMarkerSize(1);
  gr5->Draw("same p");

  TGraphErrors *gr6 = new TGraphErrors (nlines5, q2_arn78, fch_arn78, dq2_arn78, dfch_arn78); 
  gr6->SetMarkerColor(kGreen+2);
  gr6->SetMarkerStyle(20);
  gr6->SetMarkerSize(1);
  gr6->Draw("same p");

  TGraphErrors *gr7 = new TGraphErrors (nlines6, q2_dun83, fch_dun83, dq2_dun83, dfch_dun83); 
  gr7->SetMarkerColor(6);
  gr7->SetMarkerStyle(20);
  gr7->SetMarkerSize(1);
  gr7->Draw("same p");
*/

  auto ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
  ChFF_leg->AddEntry("fChFF","New ^{3}He |F_{ch}(q^{2})| Fit","l");
  ChFF_leg->AddEntry("fChFF_Amroun","New ^{3}H |F_{ch}(q^{2})| Fit","l");
  //ChFF_leg->AddEntry(gr3,"Collard 1965","p");
  //ChFF_leg->AddEntry(gr5,"Szalata 1977","p");
  //ChFF_leg->AddEntry(gr6,"Arnold 1978","p");
  //ChFF_leg->AddEntry(gr7,"Dunn 1983","p");
  //ChFF_leg->AddEntry(gr1,"Camsonne 2016","p");
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

  TF1 *fMFF = new TF1("fMFF",MFF_Q2,0.,35.,1);
  fMFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
  fMFF->Draw("L");
  fMFF->SetTitle("^{3}He Magnetic Form Factor");
  fMFF->GetHistogram()->GetYaxis()->SetTitle("|F_{m}(Q^{2})|");
  fMFF->GetHistogram()->GetYaxis()->CenterTitle(true);
  fMFF->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  fMFF->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fMFF->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fMFF->GetHistogram()->GetXaxis()->SetTitle("Q^{2} (fm^{-2})");
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

  TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,0.,35.,1);
  //cout<<fMFF_Amroun->Eval(30)<<endl;
  fMFF_Amroun->SetNpx(npdraw);
  fMFF_Amroun->SetLineColor(4);
  fMFF_Amroun->Draw("L same");

  /*
  TGraphErrors *gr2 = new TGraphErrors (nlines2, q2_fm_cam16, fm_cam16, dq2_fm_cam16, dfm_cam16); 
  gr2->SetMarkerColor(1);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1);
  gr2->Draw("same p");

  TGraphErrors *gr4 = new TGraphErrors (nlines3, q2_col65, fm_col65, dq2_col65, dfm_col65); 
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(1);
  gr4->Draw("same p");

  TGraphErrors *gr8 = new TGraphErrors (nlines7, q2_cav82, fm_cav82, dq2_cav82, dfm_cav82); 
  gr8->SetMarkerColor(2);
  gr8->SetMarkerStyle(20);
  gr8->SetMarkerSize(1);
  gr8->Draw("same p");

  TGraphErrors *gr9 = new TGraphErrors (nlines8, q2_nak01, fm_nak01, dq2_nak01, dfm_nak01); 
  gr9->SetMarkerColor(7);
  gr9->SetMarkerStyle(20);
  gr9->SetMarkerSize(1);
  gr9->Draw("same p");
  */

  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  MFF_leg->AddEntry("fMFF","New ^{3}He |F_{m}(q^{2})| Fit","l");
  MFF_leg->AddEntry("fMFF_Amroun","New ^{3}H |F_{m}(q^{2})| Fit","l");
  //MFF_leg->AddEntry(gr4,"Collard 1965","p");
  //MFF_leg->AddEntry(gr8,"Cavedon 1982 (Amroun 1994)","p");
  //MFF_leg->AddEntry(gr9,"Nakagawa 2001","p");
  //MFF_leg->AddEntry(gr2,"Camsonne 2016","p");
  MFF_leg->Draw();

  //Define the charge density from I. Sick. 
  Double_t rho_3he(Double_t *r, Double_t *par)
  {
    Double_t rho = 0;
    Double_t rho_temp = 0;
    
    for(Int_t i=0;i<ngaus;i++)
      {
	rho_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(gamma,2.) )  );
	rho = rho + rho_temp;
      }
    
    rho = Z_3He/(2*pow(pi,1.5)*pow(gamma,3.)) * rho; //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.
    
    return rho;
  }

  TCanvas* crho=new TCanvas("crho");
  crho->SetGrid();

  TF1 *f3He_Rho = new TF1("f3He_Rho",rho_3he,0.,5,1);
  f3He_Rho->SetNpx(npdraw);
  f3He_Rho->SetLineColor(2);
  f3He_Rho->Draw("L");
  f3He_Rho->SetTitle("Charge Densities");
  f3He_Rho->GetHistogram()->GetYaxis()->SetTitle("#rho e/fm^{3}");
  f3He_Rho->GetHistogram()->GetYaxis()->CenterTitle(true);
  f3He_Rho->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  f3He_Rho->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  f3He_Rho->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  f3He_Rho->GetHistogram()->GetXaxis()->SetTitle("r (fm)");
  f3He_Rho->GetHistogram()->GetXaxis()->CenterTitle(true);
  f3He_Rho->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  f3He_Rho->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  f3He_Rho->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

  //Define the charge density from I. Sick. 
  Double_t rho_3h(Double_t *r, Double_t *par)
  {
    Double_t rho = 0;
    Double_t rho_temp = 0;
    
    for(Int_t i=0;i<ngaus;i++)
      {
	rho_temp = Qich_Amroun[i]/( 1+2*pow(R_Amroun[i],2.)/pow(gamma,2.) ) * (  exp( -pow((r[0]-R_Amroun[i]),2.)/pow(gamma,2.) ) + exp( -pow((r[0]+R_Amroun[i]),2.)/pow(gamma,2.) )  );
	rho = rho + rho_temp;
      }
    
    rho = Z_3H/(2*pow(pi,1.5)*pow(gamma,3.)) * rho; //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.
    
    return rho;
  }

  TF1 *f3H_Rho = new TF1("f3H_Rho",rho_3h,0.,5,1);
  f3H_Rho->SetNpx(npdraw);
  f3H_Rho->SetLineColor(4);
  f3H_Rho->Draw("L SAME");

  auto Rho_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  Rho_leg->AddEntry("f3He_Rho","New ^{3}He Charge Density Fit","l");
  Rho_leg->AddEntry("f3H_Rho","New ^{3}H Charge Density Fit","l");
  Rho_leg->Draw();

  TCanvas* cxs=new TCanvas("cxs");
  cxs->SetGrid();
  cxs->SetLogy();

  Double_t xs(Double_t *Q2, Double_t *par)
  {
    Double_t val = 0.;
    Double_t mottxs = 0.;
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;
    Double_t fitm = 0.;
    Double_t summtemp = 0.;
    
    //All Q[0] were Q2eff previously from SOG fit code.

    //theta = 2*TMath::ASin(  pow( (1/(4*pow(E0,2.)/Q2[0]-2*E0/MtHe3)) , 0.5 )  );
    theta = 2*TMath::ASin(  pow( (1/(4*pow(E0,2.)*GeV2fm/Q2[0]-2*E0/MtHe3)) , 0.5 )  );
    theta_cor = 2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/   pow(  pow(Q2[0],0.5)/(1+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.,1./3.)))  ,2.)   )-(2.*E0/MtHe3))) , 0.5 )  );

    Ef = E0/(1.0+2.0*E0*pow(sin(theta/2.0),2.0)/MtHe3);

    //theta = 2*TMath::ASin( pow( Q2[0]/(4*E0*Ef*GeV2fm), 0.5 ) );

    //Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
    //Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
    Double_t W = E0 - Ef;
    //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
    Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2[0]  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    Double_t eta = 1.0 + Q2[0]/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2. 
    Double_t Qtot = 1.0;
    Double_t Qtemp = 0.;

    //Calculate Mott XS.
    mottxs = (  (pow(Z_3He,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta/2.0),4.0)))*pow(cos(theta/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

    //Calculate XS from FFs.
    val = mottxs * (1./eta) * ( (Q2[0]/q2_3)*pow(ChFF_Q2(Q2,par),2.) + (pow(muHe3,2.0)*Q2[0]/(2*pow(MtHe3,2)*GeV2fm))*(0.5*Q2[0]/q2_3 + pow(tan(theta/2),2))*pow(MFF_Q2(Q2,par),2.) ); 
    
    return val;
    
    //return MFF_Q2_Amroun(Q2,par) + ChFF_Q2_Amroun(Q2,par);
  }

  TF1 *fxs = new TF1("fxs",xs,yminFF,ymaxFF+54,1);
  //fMFF_Amroun->Draw("L");
  fxs->SetNpx(npdraw);
  fxs->Draw("L");
  fxs->SetTitle("^{3}He Cross Section at E_{0} = 3.356 Gev");
  fxs->GetHistogram()->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (fm^{2}/sr)");
  fxs->GetHistogram()->GetYaxis()->CenterTitle(true);
  fxs->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  fxs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fxs->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fxs->GetHistogram()->GetXaxis()->SetTitle("Q^{2} (fm^{-2})");
  fxs->GetHistogram()->GetXaxis()->CenterTitle(true);
  fxs->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  fxs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  fxs->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

  //Plot Amroun's XS as well.
  Double_t xs_Amroun(Double_t *Q2, Double_t *par)
  {
    Double_t val = 0.;
    Double_t mottxs = 0.;
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;
    Double_t fitm = 0.;
    Double_t summtemp = 0.;
    
    //All Q[0] were Q2eff previously from SOG fit code.

    theta = 2*TMath::ASin(  pow( (1/(4*pow(E0,2.)*GeV2fm/Q2[0]-2*E0/MtHe3)) , 0.5 )  );
    theta_cor = 2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/   pow(  pow(Q2[0],0.5)/(1+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.,1./3.)))  ,2.)   )-(2.*E0/MtHe3))) , 0.5 )  );
    //theta = 2*TMath::ASin( pow( Q2[0]/(4*E0*Ef*GeV2fm), 0.5 ) );

    Ef = E0/(1.0+2.0*E0*pow(sin(theta/2.0),2.0)/MtHe3);

    //theta = 2*TMath::ASin( pow( Q2[0]/(4*E0*Ef*GeV2fm), 0.5 ) );

    //Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
    //Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
    Double_t W = E0 - Ef;
    //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
    Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2[0]  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    Double_t eta = 1.0 + Q2[0]/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2. 
    Double_t Qtot = 1.0;
    Double_t Qtemp = 0.;

    //Calculate Mott XS.
    mottxs = (  (pow(Z_3H,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta/2.0),4.0)))*pow(cos(theta/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

    //Calculate XS from FFs.
    val = mottxs * (1./eta) * ( (Q2[0]/q2_3)*pow(ChFF_Q2_Amroun(Q2,par),2.) + (pow(muHe3,2.0)*Q2[0]/(2*pow(MtHe3,2.)*GeV2fm))*(0.5*Q2[0]/q2_3 + pow(tan(theta/2.),2.))*pow(MFF_Q2_Amroun(Q2,par),2.) ); 
    
    return val;
    
    //return MFF_Q2_Amroun(Q2,par) + ChFF_Q2_Amroun(Q2,par);
  }

  TF1 *fxs_Amroun = new TF1("fxs_Amroun",xs_Amroun,yminFF,ymaxFF+54,1);
  //fMFF_Amroun->Draw("L");
  fxs_Amroun->SetNpx(npdraw);
  fxs_Amroun->SetLineColor(4);
  fxs_Amroun->Draw("L SAME");

  Double_t my_q2[1], my_dq2[1], my_xs[1], my_dxs[1];
  my_q2[0] = 34.19;
  my_dq2[0] = 0;
  my_xs[0] = 1.335e-10;
  my_dxs[0] = 0.086e-10;

  TGraphErrors *gr10 = new TGraphErrors (1, my_q2, my_xs, my_dq2, my_dxs); 
  gr10->SetMarkerColor(1);
  gr10->SetMarkerStyle(20);
  gr10->SetMarkerSize(1);
  //gr10->GetHistogram()->SetLineWidth(6);
  gr10->Draw("same p");
  
  //Plot my data point.
  m1 = new TMarker(34.0981, 1.33459E-10, 20);
  m1->SetMarkerColor(kBlack);//(kOrange+7);
  m1->SetMarkerSize(1);
  //m1->Draw();

  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  MFF_leg->AddEntry("fxs","New ^{3}He Cross Section","l");
  MFF_leg->AddEntry("fxs_Amroun","^{3}He Cross Section from Amroun 1994","l");
  //MFF_leg->AddEntry(m1,"^{3}He Cross Section from Experiment E08-014","p");
  MFF_leg->AddEntry(gr10,"^{3}He Cross Section from Experiment E08-014","p");
  MFF_leg->Draw();

  //Show two lines representing the edges of our bin in Q^2.
  TLine *line = new TLine(bin_min_Q2,0.,bin_min_Q2,fxs->Eval(bin_min_Q2));
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  //line->Draw();

  TLine *line1 = new TLine(bin_max_Q2,0.,bin_max_Q2,fxs->Eval(bin_max_Q2));
  line1->SetLineColor(kBlack);
  line1->SetLineWidth(2);
  //line1->Draw();

  using namespace std;

  //Sort Chi2 array by ascending Chi^2.
  std::vector<float> Chi2_Sorted;
  //std::pair<float,float> Test;
  vector< pair <double,int> > Test;
  pair<double,int> pair;
  //std::vector< pair <float,float> > Test;
  //cout<<Chi2_Sorted.size()<<endl;

  for(Int_t i=0;i<nlines-skip;i++)
    {
      pair.first = Chi2[i];
      pair.second = i;
      Test.push_back(pair);

      //Test.first = Chi2[i];
      //Test.second = i;
      //Test.push_back( make_pair(Chi2[i],i) ); 
      //Test.push_back( pair<double,double>(Chi2[i],i) ); 
      //Test[i].first = Chi2[i];
      //Test[i].second = i;
      Chi2_Sorted.push_back(Chi2[i]);
    }  

  //Sort array in ascending order.
  sort(Chi2_Sorted.begin(),Chi2_Sorted.end());           
  for(Int_t i=0;i<nlines-skip;i++)
    {
      //cout<<"Chi2_Sorted["<<i<<"] = "<<Chi2_Sorted[i]<<endl;
    } 

  //sort(Test.begin(),Test.end());              //Will only work if compile first with root -l Plot_FFs.C+ and it's too irritating to bother making this compile nicely.
  /*
  std::sort(Test.begin(), Test.end(),
	    [](const pair& x, const pair& y) {
	      // compare second value
	      if (x.first < y.first)
		return x.first < y.first;
	    });
  */
   for(Int_t i=0;i<nlines-skip;i++)
    {
      //cout<<Test[i].second<<"    "<<Test[i].first<<endl;
      //cout<<"Fit Number = "<<Test[i].second<<"   Chi2_Sorted["<<i<<"] = "<<Test[i].first<<endl;
    }

  Double_t test_point = 35.8152;//35.8152;//35.7582;
  cout<<"My XS at "<<test_point<<" = "<<fxs->Eval(test_point)<<" fm^2/sr = "<<fxs->Eval(test_point)*1E4<<" ub/sr"<<endl;
  cout<<"Amroun's XS at test point = "<<fxs_Amroun->Eval(test_point)<<" fm^2/sr = "<<fxs_Amroun->Eval(test_point)*1E4<<" ub/sr"<<endl;
  cout<<"E0 = "<<E0<<"   Ef = "<<E0/(1.0+2.0*E0*pow(sin((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),2.0)/MtHe3)<<endl;
  cout<<"Amroun's theta = "<<2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ) * 180./pi<<endl;
  cout<<"Amroun's theta Q^2eff Correction = "<<2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/   pow(  pow(test_point,0.5)/(1+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.,1./3.)))  ,2.)   )-(2.*E0/MtHe3))) , 0.5 )  ) * 180./pi<<endl;
  cout<<"Amroun Mott = "<<(  (pow(Z_3He,2.)*((E0/(1.0+2.0*E0*pow(sin((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),2.0)/MtHe3))/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),4.0)))*pow(cos((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),2.0)  ) * 1.0/25.7<<endl;
  cout<<"fChFF(test_point) = "<<fChFF_Amroun->Eval(test_point)<<endl;
  cout<<"fMFF(test_point) = "<<fMFF_Amroun->Eval(test_point)<<endl;
  cout<<"******************************************************************"<<endl;

  Double_t test_Q2 = 35.;           //Close starting point for correct Q^2.
  Double_t model_xs = 0.;           //Model value of xs at test_Q2.
  Double_t exp_xs = 1.33459E-10;    //fm^2/sr. My experimental XS result.
  Double_t window = 0.5E-15;          //Window of acceptable Q^2 model agreement with exp XS.
  //Find where model XS equals my experimental XS and that's the Q^2 to place it at.

  //Set initial model xs guess.
  model_xs = fxs->Eval(test_Q2);
  cout<<"Initial model_xs = "<<model_xs<<" at Q^2 = "<<test_Q2<<" fm^-2."<<endl;
  
  /*
  //Increment guess up or down until close enough to the exp. xs.
  while(model_xs>(exp_xs+window) || (model_xs<exp_xs-window))
    {
      //If model xs too high.
      if(model_xs>(exp_xs+window))
	{
	  test_Q2 = test_Q2 + 0.00001;
	}
      //If model xs too low.
      if(model_xs<(exp_xs-window))
	{
	  test_Q2 = test_Q2 - 0.00001;
	}
      model_xs = fxs->Eval(test_Q2);
      //cout<<"model_xs = "<<model_xs<<" at Q^2 = "<<test_Q2<<" fm^-2 = "<<test_Q2*0.0389<<" GeV^2."<<endl;
    }
  cout<<"model_xs = "<<model_xs<<" fm^2/sr at Q^2 = "<<test_Q2<<" fm^-2 = "<<test_Q2*0.0389<<" GeV^2."<<endl;

  cout<<"Theta = "<<2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_Q2)-(2.*E0/MtHe3))) , 0.5 )  ) * 180./pi<<endl;
  */
  //cout<<"Theta = "<<2*TMath::ASin(  pow( (1/(4*pow(E0,2.)/test_Q2-2*E0/MtHe3)) , 0.5 )  )<<endl;

  //Get weighted average of Q^2 (weighted by XS value) in the bin containing my XS.
  
  Double_t dphi_cut_min = -0.03;        //Min and Max dphi cuts. Deviation from spectrometer angle.
  Double_t dphi_cut_max = 0.03;
  Double_t dphi_min = TMath::ATan(dphi_cut_min);     //Convert cut to actual radian value of the deviation from spectrometer angle.
  Double_t dphi_max = TMath::ATan(dphi_cut_max);
  Double_t spect_ang = 21.04;                        //Spectrometer angle in degrees.
  Double_t ang_min = spect_ang*pi/180.+dphi_min;     //Min and max bin angles in radians.
  Double_t ang_max = spect_ang*pi/180.+dphi_max;
  Double_t Q2_min = 4*pow(E0,2.)*pow(sin(ang_min/2.),2.) / (1+2*E0*pow(sin(ang_min/2.),2.)/MtHe3) *GeV2fm;  //Min and max Q^2 for bin (fm^-2).
  Double_t Q2_max = 4*pow(E0,2.)*pow(sin(ang_max/2.),2.) / (1+2*E0*pow(sin(ang_max/2.),2.)/MtHe3) *GeV2fm;
  Int_t nbin = 30000;                                 //Number of divisions to apply to the bin when calculating weighted average.
  Double_t bin_size = (Q2_min+Q2_max)/nbin;          //Size of each bin division.
  Double_t Q2_weight = 0.;                           //Individual Q^2 values weighted by XS.
  Double_t Q2_weight_tot = 0.;                       //Sum of individual Q^2 values weighted by XS.
  Double_t weight_tot = 0.;                          //Sum of all weighting (XS) values.
  Double_t Q2_cor = 0.;                              //Bin centering corrected Q^2 fm^-2.

  //Print values for the bin.
  cout<<"phi cut phi analyzer variable: "<<dphi_cut_min<<" to "<<dphi_cut_max<<endl;
  cout<<"phi cut radians: "<<dphi_min<<" to "<<dphi_max<<endl;
  cout<<"Spectrometer angle = "<<spect_ang<<" degrees"<<endl;
  cout<<"Min bin angle = "<<ang_min<<" radians = "<<ang_min*180/pi<<" degrees.   Max bin angle = "<<ang_max<<" radians = "<<ang_max*180/pi<<"."<<endl;
  cout<<"Min bin Q^2 = "<<Q2_min/GeV2fm<<" GeV^2 = "<<Q2_min<<" fm^-2.   Max bin Q^2 = "<<Q2_max/GeV2fm<<" GeV^2 = "<<Q2_max<<"."<<endl;

  //Find weighted average.
  for(Int_t i=0;i<nbin;i++)
    {
      Q2_weight = (Q2_min+i*bin_size)*fxs->Eval(Q2_min+i*bin_size);
      Q2_weight_tot = Q2_weight_tot + Q2_weight;
      weight_tot = weight_tot + fxs->Eval(Q2_min+i*bin_size);
    }
  //Calculate weighted average by multiplying each Q^2 value in the bin by the XS value (weight) and summing these weighted Q^2. Then divide by the sum of the weights.
  Q2_cor = Q2_weight_tot/weight_tot;
  cout<<"Average Q^2 for the bin = "<<((Q2_min+Q2_max)/2)/GeV2fm<<" GeV^2 = "<<(Q2_min+Q2_max)/2<<" fm^-2.   Corrected bin center = "<<Q2_cor/GeV2fm<<" GeV^2 = "<<Q2_cor<<" fm^-2."<<endl;
}
