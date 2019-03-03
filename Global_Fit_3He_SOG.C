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
#include "TCanvas.h" //Needed to compile canvas objects.
#include "TStopwatch.h" //Needed to compile stopwatch objects.
#include "TRandom.h" //Needed to compile random objects.
#include "TMarker.h" //Needed to compile marker objects.
#include "TLegend.h" //Needed to compile legend objects.
#include "TStyle.h" //Needed to compile gstyle objects.
#include "TSystem.h" //Needed to compile gsystem objects.

#include <TMath.h>
#include "Math/IFunction.h"
#include <cmath>
#include "Math/SpecFunc.h" //Access special functions like spherical Bessel.

#include "Math/Functor.h"
#include "Math/RichardsonDerivator.h"

Double_t xyz = 5.5;
//Double_t truncate = 5.5;

Double_t pi = 3.141592654;
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t e = 1.60217662E-19;             //Electron charge [C].
Double_t e2_nuclear = 1.4399643929E-3;             //Electron charge squared in nuclear units [GeV * fm].
Double_t alpha = 0.0072973525664;//1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.
//Double_t muHe3 = 2.9788*(3.0/1.0); //3H
//Double_t mu3H = 2.9788*(3.0/1.0); //Magnetic moment of trinucleon (H3 or He3). NIST: http://physics.nist.gov/cgi-bin/cuu/Results?search_for=magnet+moment   //MCEEP Code for H3 and He3 eleastic FFs has magnetic moments multiplied by 3.0/Z. I don't know why but it works. Maybe it's a factor of A/Z?

Int_t loops = 1;
Int_t userand = 3;                       //0 = use predetermined Ri from Amroun. 1 = use random Ri in generated in a range around Amroun's. 2 = use random Ri, ngaus=12, generated in increments of 0.1 with larger possible spacing at greater radii. 3 = use predetermined Ri for the purposes of trying to tune the fit by hand. 4 = ngaus=8. 5 = ngaus=9. 6 = ngaus=10. 7 = ngaus=11.
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 1;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t fft = 0;                           //0 = don't use FFT to try to get a charge radii. 1 = do use FFT to extract a charge radii.
Int_t Amroun_Qi = 0;                     //1 = Override fitted Qi and use Amroun's values.
Int_t showplots = 1;
Int_t useFB = 0;                         //Turn on Fourier Bessel fit.
Int_t useFB_FM = 1;                      //0 = Turn on Fourier Bessel fit just for FC. 1 = Turn on Fourier Bessel fit attempting FC and FM.
Int_t improve = 0;                       //1 = use mnimpr() to check for other minima around the one MIGRAD finds.
Int_t MINOS = 0;                         //1 = use MINOS to calculate parameter errors. With ERRordef=30, npar=24, 10000 calls took about 1.5 hours and gave results only slightly different from intial parameter errors given. Several pars were hitting limits. 
Int_t optimize_Ri = 0;                   //1 = Have code loop over each Ri value shifting it 0.1 higher and 0.1 lower until chi2 stops improving.
Int_t bootstrap = 0;                     //0 = No bootstrapping. 1 = Using a fixed Ri set randomly select points in the dataset a number of times equal to the number of points in the dataset and then use those points for a fit.
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
Int_t ngaus_Amroun = 12;
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
//Double_t Z = 1.;                         //Atomic number H3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
//Double_t MtHe3 = 3.0160492*0.9315;         //Mass of H3 in GeV.
//Double_t Mt3H = 3.0160492*0.9315;         //Mass of H3 in GeV.
Double_t Gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
//Double_t E0 = 0.5084;                    //Initial e- energy GeV.
Double_t Ef = 0.;                        //Final e- energy GeV.
Double_t ymin = 30.;//30
Double_t ymax = 100.;//100
Double_t yminFF = 0.0001;//30
Double_t ymaxFF = 6.;
Double_t range = fabs(ymaxFF - yminFF);
Int_t n = 10000;
Int_t ndim = n+1;
Int_t npdraw = 10001;                     //Number of points to be used when drawing a function.
Double_t Truncate = 100.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 2.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char str[1000];                          //Variable to read lines of the data file..
Float_t thetatemp,qefftemp,sigexptemp,uncertaintytemp,E0temp;
Float_t theta[1000];                     //Angle in degrees.
Float_t qeff[1000];                      //q effective in fm^-1.
Float_t sigexp[1000];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[1000];
Float_t E0[1000];

Int_t Amroun_pts = 57;                 //Dropped two points with no energy value give.
Int_t Collard_pts = 118;
Int_t Szlata_pts = 22;
Int_t Dunn_pts = 27;
Int_t Camsonne_pts = 18;               //Dropped two points with crazy Chi^2 values. Should reevaluate eventually.
Int_t Nakagawa_pts = 5;
Int_t my_pts = 1;
Int_t Arnold_pts = 11;                 //These XSs had to be calculated from A^1/2 function.
const Int_t datapts = 259;//Amroun_pts+Collard_pts+Szlata+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts;//248,257

Double_t m = 2.;

Double_t R[12] = {0.3,0.7,0.9,1.1,1.5,1.9,2.2,2.7,3.3,4.2,4.3,4.8};//Final 3He representative fit.

//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
//Double_t R[12] = {0.1,0.6,0.7,1.3,1.4,2.,2.7,3.4,4.3,5.2,0.,0.};//51
//Double_t R[12] = {0.2,0.7,1.3,1.5,2.1,2.8,3.6,4.2,5.2,0.,0.,0.};
//Double_t R[12] = {0.1,0.6,1.,1.5,2.1,2.4,3.,3.7,4.4,4.7,0.,0.};//9/12/18 pretty good 2.
//Double_t R[12] = {0.1,0.6,1.,1.5,2.0,2.4,3.1,3.9,4.5,4.8,0.,0.}; //9/12/18 pretty good 1. Also these Ri work well with Qi 2. Very smooth Fm.
//Double_t R[12] = {0.1,0.5,0.9,1.3,1.5,1.9,2.3,2.8,3.1,3.8,4.3,5.};//Best n=12 from Amroun starting.
//Double_t R[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_init[12] = {};
Double_t R_best[12] = {};
Double_t R_best_chi2 = 0;
//Double_t Qich[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};//Amroun
//Double_t Qim[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
//Double_t Qich[12] = {0.0289116,0.176012,0.227652,0.18408,0.186425,0.093576,0.0329847,0.052855,0.016552,0.00285043,0.00615456,1.58614E-11};
//Double_t Qim[12] = {0.0585454,0.160715,0.222426,0.156211,0.191486,0.125172,0.00162158,0.0602476,0.0228372,6.94389E-13,0.0192538,1.03119E-11};

//3He
Double_t Qich[12] = {0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.98962e-12};//3He final #30.
Double_t Qim[12] = {0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.77602e-12,0.0240927,8.94934e-12};//3He final #30.

//Double_t Qich[12] = {0.0440183,0.116665,0.202577,0.26934,0.0690628,0.179559,0.0854789,0.0318623,0.00963141,9.4369e-16,0.,0.};//9/15/18 51.
//Double_t Qim[12] = {0.0725889,0.0926902,0.209459,0.231037,0.079899,0.188591,0.0836548,0.0440054,0.0235449,2.42131e-10,0.,0.};
//Double_t Qich[12] = {0.0877489,0.157569,0.251965,0.153987,0.18569,0.114409,0.042258,0.0141493,5.25852e-12,0.,0.,0.};//9/14/18 42
//Double_t Qim[12] = {0.151938,0.0687805,0.318107,0.0627664,0.234607,0.100287,0.0562909,0.00777628,0.018992,0.,0.,0.};
//Double_t Qich[12] = {0.087498,0.157869,0.245378,0.101457,0.223892,0.123962,0.0483213,0.014929,0.00494579,0.,0.,0.};//9/14/18 40
//Double_t Qim[12] = {0.148784,0.0789473,0.304254,0.026821,0.248647,0.120223,0.05251,0.0205607,0.0136972,0.,0.,0.};
//Double_t Qich[12] = {0.0639697,0.315117,0.225156,0.13069,0.163794,0.0767875,0.0269792,0.00477479,1.95333e-11,0.,0.,0.};//9/14/18 modified 32
//Double_t Qim[12] = {0.056945,0.243831,0.158645,0.0456532,0.046148,0.0984145,0.0466099,0.0416049,0.0220442,0.,0.,0.};
//Double_t Qich[12] = {0.0509747,0.238208,0.239747,0.142914,0.189107,0.104061,0.0383111,0.0033425,0.000442121,0.,0.,0.};//9/13/18 31
//Double_t Qim[12] = {0.0993896,0.17783,0.271219,0.116363,0.185706,0.131835,0.0496486,0.0399545,0.0223591,0.,0.,0.};
//Double_t Qich[12] = {0.0639787,0.315092,0.225231,0.131057,0.163557,0.077119,0.027274,0.00349874,2.54989e-09,0.,0.,0.};//9/13/18 28
//Double_t Qim[12] = {0.107046,0.252678,0.310642,0.0425304,0.184205,0.0965293,0.0376706,0.0422947,0.0229587,0.,0.,0.};
//Double_t Qich[12] = {0.108444,0.321339,0.0334387,0.3054,0.13961,0.0594638,0.0319809,0.00944856,9.35541e-12,0.,0.,0.};//9/13/18 27
//Double_t Qim[12] = {0.160582,0.187699,0.164107,0.20848,0.164083,0.023726,0.0274688,1.31014e-08,1.78631e-09,0.,0.,0.};
//Double_t Qich[12] = {0.110049,0.316031,0.156167,0.186552,0.11344,0.0849853,0.0285813,0.0108707,2.74628e-11,0.,0.,0.};//9/13/18 25
//Double_t Qim[12] = {0.153001,0.248979,0.219662,0.100067,0.168079,0.0703413,0.0666113,0.0560969,0.0270054,0.,0.,0.};
//Double_t Qich[12] = {0.0644352,0.311729,0.274891,0.128099,0.126056,0.0656117,0.0319948,0.00588861,4.54212e-10,0.,0.,0.};//9/13/18 24
//Double_t Qim[12] = {0.10636,0.252116,0.309826,0.0695926,0.146614,0.0447945,0.0110169,1.55139e-11,1.19383e-10,0.,0.,0.};
//Double_t Qich[12] = {0.0649212,0.308528,0.299338,0.187702,0.102192,0.0379261,0.00712276,2.74347e-12,0.,0.,0.,0.};//11
//Double_t Qim[12] = {0.0997782,0.274768,0.285966,0.191182,0.109338,0.0645759,0.0435615,0.0126547,0.,0.,0.,0.};
//Double_t Qich[12] = {0.0896211,0.234994,0.367953,0.193344,0.0460444,0.0763543,1.05355e-12,1.20741e-10,0.,0.,0.,0.};//5
//Double_t Qim[12] = {0.0738571,0.287688,0.219872,0.125082,0.000139601,7.65055e-13,0.0676226,0.0808108,0.,0.,0.,0.};
//Double_t Qich[12] = {0.0440183,0.116665,0.202577,0.26934,0.0690628,0.179559,0.0854789,0.0318623,0.00963141,9.4369e-16,0.,0.};//9/15/18 51.
//Double_t Qim[12] = {0.0725889,0.0926902,0.209459,0.231037,0.079899,0.188591,0.0836548,0.0440054,0.0235449,2.42131e-10,0.,0.};
//Double_t Qich[12] = {0.0419664,0.220685,0.160807,0.293756,0.0525773,0.147678,0.0700733,0.0169877,0.00368665,6.48132E-12,0.,0.};//9/12/18 3.
//Double_t Qim[12] = {0.0701842,0.207855,0.149425,0.271277,0.064722,0.153486,0.083882,0.0352331,0.0250242,1.03941E-8,0.,0.};
//Double_t Qich[12] = {0.0411639,0.236983,0.204183,0.276745,0.0977119,0.0724304,0.0541865,0.0207177,1.26197e-07,0.00411213,0.,0.};//9/12/18 pretty good 2.
//Double_t Qim[12] = {0.07112,0.213554,0.205701,0.238441,0.162279,0.017472,0.0869983,0.0239365,0.0295892,5.58285e-10,0.,0.};
//Double_t Qich[12] = {0.0411535,0.237047,0.203676,0.277975,0.092553,0.0821181,0.0562533,0.0157997,2.61438e-10,0.00176867,0.,0.};//9/12/18 pretty good 1.
//Double_t Qim[12] = {0.0692335,0.221858,0.189445,0.254013,0.127296,0.0508983,0.0726007,0.0160226,0.0224739,4.38656e-10,0.,0.}

Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};//3He
Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};//3He

//3H
//Double_t Qich_Amroun[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121, 0.0, 0.0};//3H
//Double_t Qim_Amroun[12] = {0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077, 0.0, 0.0};//3H
Double_t Qich_best[12] = {};
Double_t Qim_best[12] = {};
Double_t av[24] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t averr[24] = {};
Double_t Qicherr[12]={}; 
Double_t Qimerr[12]={};
Float_t Q2[datapts];
Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};
Double_t Chi2_FB[datapts]={};
Double_t residual_FB[datapts]={};
Double_t FBfit[datapts]={};
Double_t E0_bs[datapts] = {};
Double_t theta_bs[datapts] = {};
Double_t sigexp_bs[datapts] = {};
Double_t uncertainty_bs[datapts] = {};
Double_t Q2_bs[datapts] = {};

Double_t  Qichtot = 0.;
Double_t  Qimtot = 0.;
Double_t amin = 0.;
Double_t maxQ2 = 0.;
TMarker *m1,*m2,*m3,*m4,*m5,*m6,*m7,*m8;

//Define SOG FFs and XS.
Double_t XS(float E0, float theta, Double_t *par)
{
  //Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
  //Double_t value = par[0] * x*x + par[1];
  Double_t val = 0.;
  Double_t mottxs = 0.;
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
  Double_t fitm = 0.;
  Double_t summtemp = 0.;
  
  Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
  Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
  Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
  
  Double_t W = E0 - Ef;
  //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
  Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
  Double_t eta = 1.0 + Q2eff/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.
  
  Double_t Qtot = 1.0;
  Double_t Qtemp = 0.;

  /*
    for(Int_t i=0;i++,ngaus)
    {
    Qtemp = Qtemp + par[i];
    }*/
  
  //Calculate Mott XS.
  mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  
  //if(par[0]+par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]+par[8]+par[9]+par[10]+par[11] == 1.)
  //{
  //cout<<"***************************************************************************"<<endl;
  /*
    for(Int_t i=0;i<2*ngaus;i++)
    {
    cout<<"par["<<i<<"] = "<<par[i]<<endl;
    }*/

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 
      //Fit just the Qi values using predetermined R[i] values.
      sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
	
      fitch =  fitch + sumchtemp;
      //cout<<"fitch["<<i<<"] = "<<fitch<<endl;
    }

  //}
  //fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
  fitch =  fitch * exp(-0.25*Q2eff*pow(Gamma,2.0));
  
  //if(par[ngaus+0]+par[ngaus+1]+par[ngaus+2]+par[ngaus+3]+par[ngaus+4]+par[ngaus+5]+par[ngaus+6]+par[ngaus+7]+par[ngaus+8]+par[ngaus+9]+par[ngaus+10]+par[ngaus+11] == 1.)
  //{
  //Define SOG for magnetic FF.
  for(Int_t i=0; i<ngaus; i++)
    {
      //Fit just the Qi values using predetermined R[i] values.
      summtemp = (par[ngaus+i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );	
      
      fitm = fitm + summtemp;
      //cout<<"fitm["<<i<<"] = "<<fitm<<endl;
    }
  //}
  fitm = fitm * exp(-0.25*Q2eff*pow(Gamma,2.0));   //For some reason had fabs(fitm).

  /*
    cout<<"E0 = "<<E0<<"   theta = "<<theta<<endl;
    cout<<"Ef = "<<Ef<<"   Q2 = "<<Q2<<"   Q2eff = "<<Q2eff<<"   W = "<<W<<"   q2_3 = "<<q2_3<<"   eta = "<<eta<<endl;
    cout<<"Mott XS = "<<mottxs<<endl;
    cout<<"fitch = "<<fitch<<"   fitm = "<<fitm<<endl;
  */
  /*  
      for(Int_t i=0;i<2*ngaus;i++)
      {
      cout<<"par["<<i<<"] = "<<par[i]<<endl;
      }
  */
  /*
    for(Int_t i=0;i<ngaus;i++)
    {
    cout<<"R["<<i<<"] = "<<R[i]<<endl;
    }
  */

  val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + (pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + pow(tan(theta*deg2rad/2),2))*pow(fitm,2.) ); 
  //cout<<"XS = "<<val<<endl;
  return val;
}
/*
Double_t FB(float E0, float theta, Double_t *par)
{
  Double_t val = 0.;

  Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
  Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
  Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
  Double_t FB_FC_sum = 0.;
  Double_t FB_FC_temp = 0.;
  Double_t FB_FM_sum = 0.;
  Double_t FB_FM_temp = 0.;
  Double_t R_FB = 5.;  //fm
  Double_t mottxs = 0.;
  Double_t tau = 0;
  Double_t rho = E0/MtHe3;

  //Calculate Mott XS.
  mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/GeV2fm;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  //cout<<"MottXS = "<<mottxs<<endl;
  //Calculate tau.
  tau = Q2eff/(4*pow(MtHe3,2.)*GeV2fm);

  //Calculate FC.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_FC_temp = ( -4 * par[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_FC_sum = FB_FC_sum + FB_FC_temp;
    }

  //Calculate FM.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_FM_temp = ( -4 * par[i-1+nFB] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_FM_sum = FB_FM_sum + FB_FM_temp;
    }

  //val = FB_sum;
  if(useFB_FM == 0)
    {
      //val = pow(Z,2.) * mottxs * pow(FB_FC_sum,2.)/(1+tau); //Just FC (FM=0). This was for the recoil proton being measure.
      val = mottxs * pow(FB_FC_sum,2.)/(1+tau);
    }
  if(useFB_FM == 1)
    {
      //val = pow(Z,2.) * mottxs * (   (  pow(FB_FC_sum,2.)+tau*pow(FB_FM_sum,2.)  )/(1+tau) + 2*tau*pow(FB_FM_sum,2.)*pow(1/(1+rho),2.)*pow(1/tan(theta*deg2rad),2.)  );  //Try to fit FC and FM. This was for the recoil proton being measure.
      val = mottxs * (   (  pow(FB_FC_sum,2.)+tau*pow(FB_FM_sum,2.)  )/(1+tau) + 2*tau*pow(FB_FM_sum,2.)*pow(tan(theta*deg2rad/2.),2.)  );
    }
  return val;
}
*/

Double_t FB(float E0, float theta, Double_t *par)
{
  Double_t val = 0.;

  Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
  Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
  Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
  Double_t FB_FC_sum = 0.;
  Double_t FB_FC_temp = 0.;
  Double_t FB_FM_sum = 0.;
  Double_t FB_FM_temp = 0.;
  Double_t R_FB = 5.;  //fm
  Double_t mottxs = 0.;
  Double_t tau = 0;
  Double_t rho = E0/MtHe3;

  Double_t eta = 1 + (2*E0/MtHe3) * pow(sin(theta*deg2rad/2.0),2.);//Beck 1984.

  //Double_t W = E0 - Ef;
  //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
  //Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
  //Double_t eta = 1.0 + Q2eff/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.

  //Calculate Mott XS.
  //mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/GeV2fm;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7. SOG formula.
  //mottxs = ( pow(Z,2.)/(4*pow(E0,2.)) ) * pow(cos(theta*deg2rad/2.0),2.)/pow(sin(theta*deg2rad/2.0),4.) * 1.0/GeV2fm;//Beck 1984.
  mottxs = ( pow(Z,2.)*pow(e2_nuclear,2.)/(4*pow(E0,2.)) ) * pow(cos(theta*deg2rad/2.0),2.)/pow(sin(theta*deg2rad/2.0),4.);//Beck 1984.
  //cout<<"MottXS = "<<mottxs<<endl;
  //Calculate tau.
  tau = Q2eff/(4*pow(MtHe3,2.)*GeV2fm);

  //Calculate FC.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_FC_temp = ( -4 * par[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_FC_sum = FB_FC_sum + FB_FC_temp;
    }

  //Calculate FM.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_FM_temp = ( -4 * par[i-1+nFB] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_FM_sum = FB_FM_sum + FB_FM_temp;
    }

  //val = FB_sum;
  if(useFB_FM == 0)
    {
      //val = pow(Z,2.) * mottxs * pow(FB_FC_sum,2.)/(1+tau); //Just FC (FM=0). This was for the recoil proton being measure.
      val = mottxs * pow(FB_FC_sum,2.)/(1+tau);
    }
  if(useFB_FM == 1)
    {
      //val = pow(Z,2.) * mottxs * (   (  pow(FB_FC_sum,2.)+tau*pow(FB_FM_sum,2.)  )/(1+tau) + 2*tau*pow(FB_FM_sum,2.)*pow(1/(1+rho),2.)*pow(1/tan(theta*deg2rad),2.)  );  //Try to fit FC and FM. This was for the recoil proton being measure.
      //val = mottxs * (   (  pow(FB_FC_sum,2.)+tau*pow(FB_FM_sum,2.)  )/(1+tau) + 2*tau*pow(FB_FM_sum,2.)*pow(tan(theta*deg2rad/2.),2.)  );
      //val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(FB_FC_sum,2.) + (pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + pow(tan(theta*deg2rad/2),2))*pow(FB_FM_sum,2.) ); //SOG XS formula.
      val = mottxs * (1/eta) * (  pow(FB_FC_sum,2.)/(1+tau) + tau*pow(muHe3,2)*pow(FB_FM_sum,2.)*( 1/(1+tau) + 2*pow(tan(theta*deg2rad/2),2.) )  ); //Beck 1984 XS formula.
    }
  return val;
}

//Create a Chi^2 function to minimize. 
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //const Int_t nbins = datapts;//177
  //Int_t i;
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  Double_t res;
  if(bootstrap == 0)
    {
      for(Int_t i=0;i<datapts;i++) 
	{
	  delta  = (sigexp[i]-XS(E0[i],theta[i],par))/uncertainty[i];
	  chisq += delta*delta;
	  Chi2[i] = delta*delta;
	  //residual[i] = (sigexp[i] - XS(E0[i],theta[i],par))/sigexp[i]; 
	  //residual[i] = fabs(sigexp[i] - XS(E0[i],theta[i],par))/XS(E0[i],theta[i],par);
	  residual[i] = (sigexp[i] - XS(E0[i],theta[i],par))/XS(E0[i],theta[i],par); 
	  xsfit[i] = XS(E0[i],theta[i],par);
	  //cout<<"xsfit["<<i<<"] = "<<xsfit[i]<<endl;
	}
    }
  if(bootstrap == 1)
    {
      for(Int_t i=0;i<datapts;i++)  
	{
	  delta  = (sigexp_bs[i]-XS(E0_bs[i],theta_bs[i],par))/uncertainty_bs[i];
	  chisq += delta*delta;
	  Chi2[i] = delta*delta;
	  //residual[i] = (sigexp_bs[i] - XS(E0_bs[i],theta_bs[i],par))/sigexp_bs[i];
	  //residual[i] = fabs(sigexp_bs[i] - XS(E0_bs[i],theta_bs[i],par))/XS(E0_bs[i],theta_bs[i],par);
	  residual[i] = (sigexp_bs[i] - XS(E0_bs[i],theta_bs[i],par))/XS(E0_bs[i],theta_bs[i],par);
	  xsfit[i] = XS(E0_bs[i],theta_bs[i],par);
	  //cout<<"xsfit["<<i<<"] = "<<xsfit[i]<<endl;
	}
    }
  f = chisq;
}

//Create a Chi^2 function to minimize with the Fourrier-Bessel series. 
void fcn_FB(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //const Int_t nbins = datapts;//177
  //Int_t i;
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  Double_t res;
  for(Int_t i=0;i<datapts;i++) 
    //for(Int_t i=0;i<226;i++) 
    {
      delta  = (sigexp[i]-FB(E0[i],theta[i],par))/uncertainty[i];
      chisq += delta*delta;
      Chi2_FB[i] = delta*delta;
      //residual_FB[i] = (sigexp[i] - FB(E0[i],theta[i],par))/sigexp[i];
      //residual_FB[i] = fabs(sigexp[i] - FB(E0[i],theta[i],par))/FB(E0[i],theta[i],par);
      residual_FB[i] = (sigexp[i] - FB(E0[i],theta[i],par))/FB(E0[i],theta[i],par);
      FBfit[i] = FB(E0[i],theta[i],par);
      //cout<<"FBfit["<<i<<"] = "<<FBfit[i]<<endl;
    }
  f = chisq;
}

//Define FFs for plotting purposes.
//Plot Charge FF Fch(Q) fm^-1.
Double_t ChFF(Double_t *Q, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      //Use SOG fit for C12 Qi coefficients and R[i] values. 
      //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

      //Convert to fm. Not sure I need to do this.
      //sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*pow(GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(Q[0]*pow(GeV2fm,0.5)*R[i])/(Q[0]*pow(GeV2fm,0.5)*R[i])) );
      sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
	
      fitch = fitch + sumchtemp;
    }
  //Convert to fm. Not sure I need to do this.
  //fitch = fitch * exp(-0.25*pow(Q[0]*pow(GeV2fm,0.5),2.)*pow(gamma,2.0));
  fitch = fitch * exp(-0.25*pow(Q[0],2.)*pow(Gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}

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
      sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
      fitch = fitch + sumchtemp;
    }
  //Convert to fm. Not sure I need to do this.
  //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
  fitch = fitch * exp(-0.25*Q2[0]*pow(Gamma,2.0));
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
      //Use SOG fit for C12 Qi coefficients and R[i] values. 
      //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

      //Convert to fm. Not sure I need to do this.
      //sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0]*GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0]*GeV2fm,0.5)*R[i])/(pow(Q2[0]*GeV2fm,0.5)*R[i])) );
      sumchtemp = (Qich_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
	
      fitch = fitch + sumchtemp;
    }
  //Convert to fm. Not sure I need to do this.
  //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
  fitch = fitch * exp(-0.25*Q2[0]*pow(Gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}

//Plot magnetic FF(Q) fm^-1.
Double_t MFF(Double_t *Q, Double_t *par)
{
  Double_t fitm = 0.;
  Double_t summtemp = 0.;
  
  //Define SOG for magnetic FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      //Use SOG fit for C12 Qi coefficients and R[i] values. 
      //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
      
      summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
      
      fitm = fitm + summtemp;
    }
  
  fitm = fitm * exp(-0.25*pow(Q[0],2.)*pow(Gamma,2.0));
  fitm = fabs(fitm);
  return fitm;
}

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
      
      summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
      
      fitm = fitm + summtemp;
    }
  
  fitm = fitm * exp(-0.25*Q2[0]*pow(Gamma,2.0));
  fitm = fabs(fitm);
  return fitm;
}

//Reset Ri to equal Amroun's values.
//for(Int_t i=0;i<12;i++)
// {
//	 R[i] = R_Amroun[i];
// }

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
      
      summtemp = (Qim_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
      
      fitm = fitm + summtemp;
    }
  
  fitm = fitm * exp(-0.25*Q2[0]*pow(Gamma,2.0));
  fitm = fabs(fitm);
  return fitm;
}

void print_fit()//Never finished since I'm not sure this will be useful.
{
  TCanvas* c_FF=new TCanvas("c_FF");
  c_FF->Divide(1,2);
  c_FF->cd(1)->SetLogy();
  c_FF->cd(1)->SetGrid();
  c_FF->cd(1);
  
  TF1 *fChFF1 = new TF1("fChFF1",ChFF_Q2,yminFF,ymaxFF+54,1);
  fChFF1->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
  fChFF1->Draw("L");
  c_FF->SetTitle("Current Fit Form Factors");
  fChFF1->SetTitle("^{3}He Charge Form Factor");
  fChFF1->GetHistogram()->GetYaxis()->SetTitle("|F_{ch}(Q^{2})|");
  fChFF1->GetHistogram()->GetYaxis()->CenterTitle(true);
  fChFF1->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  fChFF1->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fChFF1->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fChFF1->GetHistogram()->GetXaxis()->SetTitle("Q^{2} (fm^{-2})");
  fChFF1->GetHistogram()->GetXaxis()->CenterTitle(true);
  fChFF1->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  fChFF1->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  fChFF1->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);

  c_FF->cd(2)->SetLogy();
  c_FF->cd(2)->SetGrid();
  c_FF->cd(2);
  /*
    TF1 *fMFF = new TF1("fMFF",MFF_Q2,yminFF,ymaxFF+54,1);
    fMFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
    fMFF->Draw("L");
    //c4->SetTitle("He3 Magnetic Form Factor");
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

    c_FF->SaveAs("Test_Print.png");
    //c1->SaveAs(Form("./%i/%s/pictures_png/%i_%s_wire%i_wire%i.png",run,plane,run,plane,wire1,wire2));
    */
}

Double_t fitg(Double_t *Q, Double_t *par)
{
  Double_t val = 0.;
  
  //Show Gaussian Part of FFs.
  val = (par[0]/(1.0+2.0*pow(par[1],2.0)/pow(Gamma,2.0))) * ( cos(Q[0]*par[1]) + (2.0*pow(par[1],2.0)/pow(Gamma,2.0)) * (sin(Q[0]*par[1])/(Q[0]*par[1])) );
  
  val = val * exp(-0.25*pow(Q[0],2.)*pow(Gamma,2.0));

  return val;
}

Double_t fitg_rho(Double_t *r, Double_t *par)
{
  Double_t val = 0.;

  //Show Gaussian part of rho.
  val = par[0]/( 1+2*pow(par[1],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-par[1]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+par[1]),2.)/pow(Gamma,2.) )  );

  val = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * val;

  return val;
}

//Define the charge density from I. Sick. 
Double_t rho_ch(Double_t *r, Double_t *par)
{
  Double_t rho = 0;
  Double_t rho_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(Gamma,2.) )  );
      rho = rho + rho_temp;
    }

  rho = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * rho; //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho;
}

//Create a function that can be integrated to check that the normilaization to Ze is correct.
Double_t rho_ch_int(Double_t *r, Double_t *par)
{
  Double_t rho_int = 0;
  Double_t rho_int_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_int_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(Gamma,2.) )  );
      rho_int = rho_int + rho_int_temp;
    }

  rho_int = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * rho_int * 4*pi*pow(r[0],2.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho_int;
}

//Create a function to calculate rms radius.
Double_t rho_rms(Double_t *r, Double_t *par)
{
  Double_t rho_rms = 0;
  Double_t rho_rms_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_rms_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(Gamma,2.) )  );
      rho_rms = rho_rms + rho_rms_temp;
    }

  rho_rms = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * rho_rms * 4*pi*pow(r[0],4.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho_rms;
}

Double_t ChFF_Deriv(Double_t Q2) 
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
   
  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    {
      sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2,0.5)*R[i])/(pow(Q2,0.5)*R[i])) );
      fitch = fitch + sumchtemp;
    }
  fitch = fitch * exp(-0.25*Q2*pow(Gamma,2.0));
  //fitch = fabs(fitch);
  return fitch;
}

//Plot FC and FM for FB fits.
//Calculate FC.
Double_t FB_FC(Double_t *Q2, Double_t *par)
{
  Double_t FB_FC_temp = 0;
  Double_t FB_FC_sum = 0;
  Double_t R_FB = 5.;
  
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_FC_temp = ( -4 * av[i-1] * sin( pow(Q2[0],0.5) * R_FB ) ) / ( pow(Q2[0],0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2[0] - pow(i*pi/R_FB,2.))  );
      FB_FC_sum = FB_FC_sum + FB_FC_temp;
      //cout<<"av["<<i-1<<"] = "<<av[i-1]<<endl;
    }
  return FB_FC_sum;
}

//Plot FC and FM for FB fits.
//Calculate FM.
Double_t FB_FM(Double_t *Q2, Double_t *par)
{
  Double_t FB_FM_temp = 0;
  Double_t FB_FM_sum = 0;
  Double_t R_FB = 5.;
  
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_FM_temp = ( -4 * av[i-1+nFB] * sin( pow(Q2[0],0.5) * R_FB ) ) / ( pow(Q2[0],0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2[0] - pow(i*pi/R_FB,2.))  );
      FB_FM_sum = FB_FM_sum + FB_FM_temp;
    }
  return FB_FM_sum;
}

//Define the charge density from I. Sick. 
Double_t rho_ch_FB(Double_t *r, Double_t *par)
{
  Double_t rho = 0;
  Double_t rho_temp = 0;
  Double_t R_FB = 5.;
   
  for(Int_t i=0;i<nFB;i++)
    {
      rho_temp = av[i] * ROOT::Math::sph_bessel(0,(i*pi/R_FB)*r[0]);
      rho = rho + rho_temp;
    }

  return rho;
}

//Create a function to calculate rms radius.
Double_t FB_rho_rms(Double_t *r, Double_t *par)
{
  Double_t rho = 0;
  Double_t rho_temp = 0;
  Double_t R_FB = 5.;
   
  for(Int_t i=0;i<nFB;i++)
    {
      rho_temp = av[i] * ROOT::Math::sph_bessel(0,(i*pi/R_FB)*r[0]);
      rho = rho + rho_temp;
    }

  rho = rho * 4*pi*pow(r[0],4.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho;
}

//Create a function that can be integrated to check that the normilaization to Ze is correct.
Double_t FB_rho_ch_int(Double_t *r, Double_t *par)
{
  Double_t rho = 0;
  Double_t rho_temp = 0;
  Double_t R_FB = 5.;
   
  for(Int_t i=0;i<nFB;i++)
    {
      rho_temp = av[i] * ROOT::Math::sph_bessel(0,(i*pi/R_FB)*r[0]);
      rho = rho + rho_temp;
    }

  rho = rho * 4*pi*pow(r[0],2.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho;
}

Double_t test(Double_t X)
{
  cout<<pi<<endl;
  Double_t val = X;
  return val;
}

void Global_Fit_3He_SOG()
{

  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Make a new canvas to plot data.
  //TCanvas* c1=new TCanvas("c1");
  //c1->SetGrid();

  //Read in data from text file.
  //Open file.
  FILE *fp;
  if(usedifmin == 0)
    {
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_640.txt","r");
      fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data.txt","r");
    }

  if(usedifmin == 1)
    {//File *fp = new File;
      //File *fp = new File(“a.txt”);
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_640.txt","r");
      fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data.txt","r");
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data_No_New.txt","r");
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data_Removed_Bad_Chi2.txt","r");
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data_Low_Chi2_Only.txt","r");
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data_Short.txt","r");
    }

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
	ncols = fscanf(fp,"%f %f %f %f",&E0temp, &thetatemp, &sigexptemp, &uncertaintytemp);
	if (ncols < 0) break;   
	//cout<<"ncols = "<<ncols<<endl;
	//cout<<thetatemp<<"   "<<qefftemp<<"   "<<sigexptemp<<"   "<<uncertaintytemp<<endl;
	E0[nlines-skip] = E0temp;
	theta[nlines-skip] = thetatemp;
	sigexp[nlines-skip] = sigexptemp;
        uncertainty[nlines-skip] = uncertaintytemp; 

	Q2[nlines-skip] = 4 * E0[nlines-skip] * (E0[nlines-skip]/(1.0+2.0*E0[nlines-skip]*pow(sin(theta[nlines-skip]*deg2rad/2.0),2.0)/MtHe3)) * pow(sin(theta[nlines-skip]*deg2rad/2.0),2.) * GeV2fm;

	nlines++;
      }
  }

  //Print the data read from the file.
  if(showplots == 1)
    { 
      for(int i=0; i<(nlines-skip); i++)
	{
	  cout<<"E0["<<i<<"] = "<<E0[i]<<"   theta["<<i<<"] = "<<theta[i]<<"   Q^2["<<i<<"] = "<<Q2[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   uncertainty["<<i<<"] = "<<uncertainty[i]<<endl;
	}
      
      cout<<"Number of lines = "<<nlines<<endl;
    }
  fclose(fp);

  //Create an output file to store fit results.
  std::ofstream output ("Ri_Chi2.txt", std::ofstream::out);
  output<<"Chi2   rChi2   BIC   AIC    Qichtot    Qimtot  R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch    Q0m    Q1m    Q2m    Q3m    Q4m    Q5m    Q6m    Q7m    Q8m    Q9m    Q10m    Q11m"<<endl;

  //Begin loop over fit with different Ri values each time.
  for(Int_t q=0;q<loops;q++)
    {
      //If using a bootstrap create a separate set of the main variable arrays and fill them randomly from the data.
      if(bootstrap == 1)
	{
	  for(Int_t i=0;i<datapts;i++)
	    {
	      //Select random datapoints from the dataset. Allow selecting one point multiple times. 
	      gRandom->SetSeed(0);                    //Sets new random seed.
	      TF1 *rand = new TF1("rand","x",0.,datapts-1);
	      Int_t j=0;
	      j = TMath::Nint(rand->GetRandom());
	     
	      E0_bs[i] = E0[j];
	      theta_bs[i] = theta[j];
	      sigexp_bs[i] = sigexp[j];
	      uncertainty_bs[i] = uncertainty[j];
	      Q2_bs[i] = Q2[j];
	      cout<<"j = "<<j<<"   E0_bs["<<i<<"] = "<<E0_bs[i]<<"   theta_bs["<<i<<"] = "<<theta_bs[i]<<"   Q^2_bs["<<i<<"] = "<<Q2_bs[i]<<"   sigexp_bs["<<i<<"] = "<<sigexp_bs[i]<<"   uncertainty_bs["<<i<<"] = "<<uncertainty_bs[i]<<endl;
	    }
	}
     
      if(userand == 1)
	{
	  //Generate random R[i] values. 
	  Double_t d = 0.25;//0.49
	  Double_t step = 0.5;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",0.,.01);
	  R[0] = 0.1;//rand->GetRandom();
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
	}

      if(userand == 2)
	{
	  //Generate random R[i] values. 
	  Double_t d = 0.49;
	  Double_t step = 0.5;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",0.,.01);
	  R[0] = 0.1;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",3.,4.);
	  R[1] = TMath::Nint(rand1->GetRandom())/10.+R[0];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",3.,4.);
	  R[2] = TMath::Nint(rand2->GetRandom())/10.+R[1];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",3.,4.);
	  R[3] = TMath::Nint(rand3->GetRandom())/10.+R[2];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",3.,4.);
	  R[4] = TMath::Nint(rand4->GetRandom())/10.+R[3];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",3.,4.);
	  R[5] = TMath::Nint(rand5->GetRandom())/10.+R[4];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",3.,4.);
	  R[6] = TMath::Nint(rand6->GetRandom())/10.+R[5];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",5.,6.);
	  R[7] = TMath::Nint(rand7->GetRandom())/10.+R[6];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",5.,6.);
	  R[8] = TMath::Nint(rand8->GetRandom())/10.+R[7];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand9 = new TF1("rand9","x",5.,6.);
	  R[9] = TMath::Nint(rand9->GetRandom())/10.+R[8];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand10 = new TF1("rand10","x",5.,6.);
	  R[10] = TMath::Nint(rand10->GetRandom())/10.+R[9];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand11 = new TF1("rand11","x",5.,6.);
	  R[11] = TMath::Nint(rand11->GetRandom())/10.+R[10];
	}
     
      if(userand == 3)
	{
	  //Generate random R[i] values. 
	  Double_t d = 0.49;
	  Double_t step = 0.5;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",0.,.01);
	  R[0] = 0.3;//0.1;//0.1;//0.1;//0.1;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",3.,4.);
	  R[1] = 0.7;//0.5;//0.5;//0.5;//0.5;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",3.,4.);
	  R[2] = 0.9;//0.9;//0.9;//0.9;//0.9;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",3.,4.);
	  R[3] = 1.1;//1.3;//1.3;//1.3;//1.3;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",3.,4.);
	  R[4] = 1.5;//1.7;//1.5;//1.6;//1.5;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",3.,4.);
	  R[5] = 1.6;//2.3;//1.9;//2.1;//1.9;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",3.,4.);
	  R[6] = 2.2;//2.4;//2.3;//2.4;//2.3;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",5.,6.);
	  R[7] = 2.7;//3.;//2.8;//2.9;//2.8;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",5.,6.);
	  R[8] = 3.3;//3.6;//3.1;//3.4;//3.1;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand9 = new TF1("rand9","x",5.,6.);
	  R[9] = 4.2;//3.9;//3.8;//3.8;//3.8;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand10 = new TF1("rand10","x",5.,6.);
	  R[10] = 4.3;//4.3;//4.3;//4.6;//4.3;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand11 = new TF1("rand11","x",5.,6.);
	  R[11] = 4.8;//5.2;//5.;//5.;//5.;
	}
     
      if(userand == 4) //ngaus = 8
	{
	  //Generate random R[i] values. 
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",2.,3.);
	  R[0] = TMath::Nint(rand->GetRandom())/10.;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",5.,6.);
	  R[1] = TMath::Nint(rand1->GetRandom())/10.+R[0];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",5.,6.);
	  R[2] = TMath::Nint(rand2->GetRandom())/10.+R[1];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",5.,6.);
	  R[3] = TMath::Nint(rand3->GetRandom())/10.+R[2];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",5.,6.);
	  R[4] = TMath::Nint(rand4->GetRandom())/10.+R[3];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",8.,9.);
	  R[5] = TMath::Nint(rand5->GetRandom())/10.+R[4];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",8.,9.);
	  R[6] = TMath::Nint(rand6->GetRandom())/10.+R[5];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",8.,9.);
	  R[7] = TMath::Nint(rand7->GetRandom())/10.+R[6];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	}

      if(userand == 5) //ngaus = 9
	{
	  //Generate random R[i] values. 
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",2.,3.);
	  R[0] = TMath::Nint(rand->GetRandom())/10.;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",3.,4.);
	  R[1] = TMath::Nint(rand1->GetRandom())/10.+R[0];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",3.,4.);
	  R[2] = TMath::Nint(rand2->GetRandom())/10.+R[1];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",3.,4.);
	  R[3] = TMath::Nint(rand3->GetRandom())/10.+R[2];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",3.,4.);
	  R[4] = TMath::Nint(rand4->GetRandom())/10.+R[3];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",7.,8.);
	  R[5] = TMath::Nint(rand5->GetRandom())/10.+R[4];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",7.,8.);
	  R[6] = TMath::Nint(rand6->GetRandom())/10.+R[5];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",7.,8.);
	  R[7] = TMath::Nint(rand7->GetRandom())/10.+R[6];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",7.,8.);
	  R[8] = TMath::Nint(rand8->GetRandom())/10.+R[7];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	}

      if(userand == 6) //ngaus = 10
	{
	  //Generate random R[i] values. 
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",1.,2.);
	  R[0] = TMath::Nint(rand->GetRandom())/10.;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",4.,5.);
	  R[1] = TMath::Nint(rand1->GetRandom())/10.+R[0];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",4.,5.);
	  R[2] = TMath::Nint(rand2->GetRandom())/10.+R[1];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",4.,5.);
	  R[3] = TMath::Nint(rand3->GetRandom())/10.+R[2];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",4.,5.);
	  R[4] = TMath::Nint(rand4->GetRandom())/10.+R[3];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",4.,5.);
	  R[5] = TMath::Nint(rand5->GetRandom())/10.+R[4];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",7.,8.);
	  R[6] = TMath::Nint(rand6->GetRandom())/10.+R[5];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",7.,8.);
	  R[7] = TMath::Nint(rand7->GetRandom())/10.+R[6];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",7.,8.);
	  R[8] = TMath::Nint(rand8->GetRandom())/10.+R[7];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand9 = new TF1("rand9","x",7.,8.);
	  R[9] = TMath::Nint(rand9->GetRandom())/10.+R[8];
	}

      if(userand == 7) //ngaus = 11
	{
	  //Generate random R[i] values. 
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",1.,2.);
	  R[0] = TMath::Nint(rand->GetRandom())/10.;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",3.,4.);
	  R[1] = TMath::Nint(rand1->GetRandom())/10.+R[0];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",3.,4.);
	  R[2] = TMath::Nint(rand2->GetRandom())/10.+R[1];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",3.,4.);
	  R[3] = TMath::Nint(rand3->GetRandom())/10.+R[2];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",3.,4.);
	  R[4] = TMath::Nint(rand4->GetRandom())/10.+R[3];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",3.,4.);
	  R[5] = TMath::Nint(rand5->GetRandom())/10.+R[4];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",6.,7.);
	  R[6] = TMath::Nint(rand6->GetRandom())/10.+R[5];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",6.,7.);
	  R[7] = TMath::Nint(rand7->GetRandom())/10.+R[6];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",6.,7.);
	  R[8] = TMath::Nint(rand8->GetRandom())/10.+R[7];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand9 = new TF1("rand9","x",6.,7.);
	  R[9] = TMath::Nint(rand9->GetRandom())/10.+R[8];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand10 = new TF1("rand10","x",6.,7.);
	  R[10] = TMath::Nint(rand10->GetRandom())/10.+R[9];
	}

      //Add a constant to each value of Ri.
      /* 
	 for(Int_t i=0;i<ngaus;i++)
	 {
	 R[i] = R[i] + 0.01;
	 }
      */

      //Print the Ri used for the minimization. 
      if(showplots == 1)
	{
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      cout<<"R["<<i<<"] = "<<R[i]<<endl;
	    }
	}
      //Initiate Minuit for minimization.
      TMinuit *gMinuit = new TMinuit(24);  //initialize TMinuit with a maximum of 24 params
      gMinuit->SetFCN(fcn);

      Double_t arglist[10];
      Int_t ierflg = 0;

      arglist[0] = 30.; //1 is for simple chi^2. For multiparameter errors this needs to be increased. 30 adds ~ 1 min.
      gMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //Set the ERRordef or UP. 
  
      //Set step sizes.
      static Double_t stepsize[4] = {0.001 , 0.1 , 0.01 , 0.001};
  
      //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
      for(Int_t i=0;i<ngaus;i++)
	{
	  gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], 0.,1.,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], 0.,0.,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.001,Qich[i]+0.001,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.000000000001,Qich[i]+0.000000000001,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.05,Qich[i]+0.05,ierflg);
	}
      for(Int_t i=0;i<ngaus;i++)
	{
	  gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], 0.,1.,ierflg);
	  //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], 0.,0.,ierflg);
	  //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], Qim[i]-0.001,Qim[i]+0.001,ierflg);
	  //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], Qim[i]-0.000000000001,Qim[i]+0.000000000001,ierflg);
	  //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], Qim[i]-0.05,Qim[i]+0.05,ierflg);
	}
  
      // Now ready for minimization step
      arglist[0] = 10000.;//Max calls. 50000.
      arglist[1] = 0.1;//Tolerance for convergance. 1 seems to give the same results as 0.1.
      //cout<<"Sup1"<<endl;
      gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
      if(improve == 1)
	{
	  gMinuit->mnimpr(); //Check for other minima to see if we're trapped in a local minima.
	}

      //Implement MINOS to calculate multiparameter errors. arglist[0] is giving # of calls again.
      if(MINOS == 1)
	{
	  gMinuit->mnexcm("MINOS", arglist, 1, ierflg);
	}

      //cout<<"Sup2"<<endl;
      // Print results
      Double_t edm,errdef;//Moved amin to global.
      Int_t nvpar,nparx,icstat;
      gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
      //gMinuit->mnprin(3,amin);

      //Attempt to optimize the Ri by varying them systematically.
      if(optimize_Ri == 1)
	{
	  Int_t chi2_init = amin;
	  R_best_chi2 = amin;
	  cout<<"Initial Chi^2 = "<<amin<<endl;
	  R_init[0] = R[0];
	  R_best[0] = R[0];

	  //Set initial best Ri and Qi values to the initial Ri values.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      gMinuit->GetParameter(i,Qich[i],Qicherr[i]);
	      gMinuit->GetParameter(ngaus+i,Qim[i],Qimerr[i]);
	    }
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      R_best[i] = R[i];
	      Qich_best[i] = Qich[i];
	      Qim_best[i] = Qim[i];
	    }

	  for(Int_t i=0;i<ngaus;i++)
	    {
	      Int_t chi2_better = 1;  //Test if the change improved chi2.
	      R_init[i] = R[i];
	      cout<<"*********************************************************"<<endl;
	      cout<<"Optimizing Initial R["<<i<<"] = "<<R[i]<<endl;
	  
	      //Check for better Ri values below the initial Ri value.
	      while(chi2_better == 1)
		{
		  //amin = 0.;
		  //Protect against R[i]=0. while still basically testing R[i]=0.
		  if(R[i]==0.1)
		    {
		      R[i] = 0.0001; 
		    }
		  else
		    {
		      R[i] = R[i] - 0.1;
		    }
		  cout<<"R["<<i<<"] set to "<<R[i]<<endl;
	      
		  //arglist[0] = 30.; //1 is for simple chi^2. For multiparameter errors this needs to be increased. 30 adds ~ 1 min.
		  //gMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //Set the ERRordef or UP.
	      
		  //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //gMinuit->mnparm(j, Form("Qich%d",j+1), Qich[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are Amroun's Qi.
		      gMinuit->mnparm(j, Form("Qich%d",j+1), Qich_best[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are best fit Qi values so far.
		    }
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //gMinuit->mnparm(ngaus+j, Form("Qim%d",j+1), Qim[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are Amroun's Qi.
		      gMinuit->mnparm(ngaus+j, Form("Qim%d",j+1), Qim_best[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are best fit Qi values so far.
		    }
	      
		  arglist[0] = 10000.;//Max calls. 50000.
		  arglist[1] = 0.1;//Tolerance for convergance. 1 seems to give the same results as 0.1.
		  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
		  if(improve == 1)
		    {
		      gMinuit->mnimpr();
		    }
		  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		  cout<<"- Updated R["<<i<<"] = "<<R[i]<<"   New Chi^2 = "<<amin<<endl;
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //Set Qi values equal to the fitted parameters.
		      gMinuit->GetParameter(j,Qich[j],Qicherr[j]);
		      gMinuit->GetParameter(ngaus+j,Qim[j],Qimerr[j]);
		      cout<<"R["<<j<<"] = "<<R[j]<<"   Qich["<<j<<"] = "<<Qich[j]<<"   Qim["<<j<<"] = "<<Qim[j]<<endl;
		    }
		  if(amin<R_best_chi2)
		    {
		      for(Int_t j=0;j<ngaus;j++)
			{
			  Qich_best[j] = Qich[j];
			  Qim_best[j] = Qim[j];
			}
		      R_best[i] = R[i]; //Store best Ri value thus far.
		      R_best_chi2 = amin; //Store the best updated chi2 so far.
		      chi2_better = 1;
		    }
		  else
		    {
		      R[i] = R_best[i]; //No improvement so reset Ri to the best previous value.
		      cout<<"No improvement found. Setting R["<<i<<"] to best value = "<<R_best[i]<<".   With best Chi^2 = "<<R_best_chi2<<"."<<endl;
		      chi2_better = 0;
		    }
		} 
	      //Check for better Ri values above the initial Ri value.
	      chi2_better = 1;
	      R[i] = R_init[i];  //Start from initial Ri and move up now.
	      while(chi2_better == 1)
		{
		  R[i] = R[i] + 0.1;
		  cout<<"R["<<i<<"] set to "<<R[i]<<endl;

		  //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //gMinuit->mnparm(j, Form("Qich%d",j+1), Qich[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are Amroun's Qi.
		      gMinuit->mnparm(j, Form("Qich%d",j+1), Qich_best[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are best fit Qi values so far.
		    }
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //gMinuit->mnparm(ngaus+j, Form("Qim%d",j+1), Qim[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are Amroun's Qi.
		      gMinuit->mnparm(ngaus+j, Form("Qim%d",j+1), Qim_best[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are best fit Qi values so far.
		    }

		  arglist[0] = 10000.;//Max calls. 50000.
		  arglist[1] = 0.1;//Tolerance for convergance. 1 seems to give the same results as 0.1.
		  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
		  if(improve == 1)
		    {
		      gMinuit->mnimpr();
		    }
		  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		  cout<<"+ Updated R["<<i<<"] = "<<R[i]<<"   New Chi^2 = "<<amin<<endl;
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //Set Qi values equal to the fitted parameters.
		      gMinuit->GetParameter(j,Qich[j],Qicherr[j]);
		      gMinuit->GetParameter(ngaus+j,Qim[j],Qimerr[j]);
		      cout<<"R["<<j<<"] = "<<R[j]<<"   Qich["<<j<<"] = "<<Qich[j]<<"   Qim["<<j<<"] = "<<Qim[j]<<endl;
		    }
		  if(amin<R_best_chi2)
		    {
		      for(Int_t j=0;j<ngaus;j++)
			{
			  Qich_best[j] = Qich[j];
			  Qim_best[j] = Qim[j];
			}
		      R_best[i] = R[i]; //Store best Ri value thus far.
		      R_best_chi2 = amin; //Store the best updated chi2 so far.
		      chi2_better = 1;
		    }
		  else
		    {
		      R[i] = R_best[i]; //No improvement so reset Ri to the best previous value.
		      cout<<"No improvement found. Setting R["<<i<<"] to final best value = "<<R_best[i]<<".   With best Chi^2 = "<<R_best_chi2<<"."<<endl;
		      chi2_better = 0;
		    }
		} 
	    }
	  //Final fit using all of the Ri_best values.
	  cout<<"Preforming final fit with optimized R[i] values."<<endl;
	  for(Int_t j=0;j<ngaus;j++)
	    {
	      cout<<"R["<<j<<"] = "<<R[j]<<endl;
	    }
	  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
	  if(improve == 1)
	    {
	      gMinuit->mnimpr();
	    }
	  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	}//End optimize Ri.

      //Print Chi^2 and residual for each data point.
      if(showplots == 1)
	{ 
	  for(Int_t i=0;i<(nlines-skip);i++)
	    {
	      cout<<"Chi2["<<i<<"] = "<<Chi2[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   XSfit(E0,theta,par)["<<i<<"] = "<<xsfit[i]<<"   XSexp/XSfit = "<<sigexp[i]/xsfit[i]<<"   residual["<<i<<"] = "<<residual[i]<<endl;
	    }
	}

      Double_t maxchi2 = 0.;
      Double_t total_chi2 = 0.;
      //Double_t maxQ2 = 0.; //Moved to global.
      if(showplots == 1)
	{
	  //Create a plot of Chi^2 vs. theta.
	  TCanvas* c1=new TCanvas("c1");
	  c1->SetGrid();

	  for(Int_t i=0;i<(nlines-skip);i++)
	    {
	      //Store max chi2.
	      if(Chi2[i]>maxchi2)
		{
		  maxchi2 = Chi2[i];
		  cout<<"maxchi2 = "<<maxchi2<<endl;
		}
	      //Sum all the individual points' chi2 (should match fcn output value).
	      total_chi2 = total_chi2 + Chi2[i];
	    }
	  cout<<"Total Chi^2 = "<<total_chi2<<endl;

	  //Find max Q^2 for things like FB fitting.
	  for(Int_t i=0;i<(datapts);i++)
	    {
	      //Store max chi2.
	      if(Q2[i]>maxQ2)
		{
		  maxQ2 = Q2[i];
		  cout<<"maxQ2 = "<<maxQ2<<endl;
		}
	    }
      
	  TH2D *hchi = new TH2D("hchi","\Chi^{2} vs. Scattering Angle" , 161, 0., 160., maxchi2+21,0., maxchi2+20);
	  for(Int_t i=0;i<(nlines-skip);i++)
	    {
	      //hchi->Fill(theta[i],Chi2[i]);
	    }
	  hchi->SetMarkerStyle(20);
	  hchi->SetMarkerSize(1);
	  gStyle->SetOptStat(0);
	  hchi->Draw();

	  //Plot 59 Amroun data points. Removed two points without energy listed.
	  for (Int_t i=0;i<Amroun_pts;i++) 
	    {
	      m1 = new TMarker(theta[i], Chi2[i], 20);
	      m1->SetMarkerColor(2);
	      m1->SetMarkerSize(1);
	      m1->Draw();
	    }

	  //Plot 118 Collard 1965 (Amroun ref 5) data points.
	  for (Int_t i=Amroun_pts;i<(Amroun_pts+Collard_pts);i++) 
	    {
	      m2 = new TMarker(theta[i], Chi2[i], 20);
	      m2->SetMarkerColor(4);
	      m2->SetMarkerSize(1);
	      m2->Draw();
	    }

	  //Plot 22 Szlata 1977 (Amroun ref 8) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts);i<(Amroun_pts+Collard_pts+Szlata_pts);i++) 
	    {
	      m3 = new TMarker(theta[i], Chi2[i], 20);
	      m3->SetMarkerColor(3);
	      m3->SetMarkerSize(1);
	      m3->Draw();
	    }

	  //Plot 27 Dunn 1983 (Amroun ref 10) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i++) 
	    {
	      m4 = new TMarker(theta[i], Chi2[i], 20);
	      m4->SetMarkerColor(6);
	      m4->SetMarkerSize(1);
	      m4->Draw();
	    }

	  //Plot 16 (skipping 2) JLab (Camsonne) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i++) 
	    {
	      m5 = new TMarker(theta[i], Chi2[i], 20);
	      m5->SetMarkerColor(1);
	      m5->SetMarkerSize(1);
	      m5->Draw();
	    }

	  //Plot 5 Nakagawa 2001 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i++) 
	    {
	      m6 = new TMarker(theta[i], Chi2[i], 20);
	      m6->SetMarkerColor(7);
	      m6->SetMarkerSize(1);
	      m6->Draw();
	    }

	  //Plot my data point.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i++) 
	    {
	      m7 = new TMarker(theta[i], Chi2[i], 20);
	      m7->SetMarkerColor(kOrange+7);
	      m7->SetMarkerSize(1);
	      m7->Draw();
	    }

	  //Plot 11 Arnold 1978 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts);i++) 
	    {
	      m8 = new TMarker(theta[i], Chi2[i], 20);
	      m8->SetMarkerColor(kGreen+2);
	      m8->SetMarkerSize(1);
	      m8->Draw();
	    }

	  //auto legend1 = new TLegend(0.1,0.7,0.48,0.9); //Places legend in upper left corner of histogram.
	  TLegend *legend1;
	  legend1 = new TLegend(0.62,0.7,0.9,0.9); //Places legend in upper right corner of histogram.
	  legend1->AddEntry(m2,"Collard 1965","p");
	  legend1->AddEntry(m3,"Szlata 1977","p");
	  legend1->AddEntry(m8,"Arnold 1978","p");
	  legend1->AddEntry(m4,"Dunn 1983","p");
	  legend1->AddEntry(m1,"Amroun 1994","p");
	  legend1->AddEntry(m6,"Nakagawa 2001","p");
	  legend1->AddEntry(m5,"Camsonne 2016","p");
	  legend1->AddEntry(m7,"Current Analysis, Ye 2018","p");
	  legend1->Draw();

	  //Save the canvas as a .png file and a .C file.
	  c1->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Chi2_vs_Angle.png");
	  c1->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Chi2_vs_Angle.C");

	  //Create a plot of Chi^2 vs. Q^2.
	  TCanvas* cQ2=new TCanvas("cQ2");
	  cQ2->SetGrid();

	  TH2D *hQ2 = new TH2D("hQ2","Representative Fit #chi^{2} vs. Q^{2}" , datapts+1, 0., maxQ2+2, 500, 0., maxchi2+10);
	  for(Int_t i=0;i<datapts;i++)
	    {
	      //hQ2->Fill(Q2[i],Chi2[i]);
	    }
	  hQ2->SetMarkerStyle(20);
	  hQ2->SetMarkerSize(1);
	  gStyle->SetOptStat(0);
	  hQ2->Draw("p");

	  //Plot 59 Amroun data points. Removed two points without energy listed.
	  for (Int_t i=0;i<Amroun_pts;i++) 
	    {
	      TMarker *m1 = new TMarker(Q2[i], Chi2[i], 20);
	      m1->SetMarkerColor(2);
	      m1->SetMarkerSize(1);
	      m1->Draw();
	    }

	  //Plot 118 Collard 1965 (Amroun ref 5) data points.
	  for (Int_t i=Amroun_pts;i<(Amroun_pts+Collard_pts);i++) 
	    {
	      TMarker *m2 = new TMarker(Q2[i], Chi2[i], 20);
	      m2->SetMarkerColor(4);
	      m2->SetMarkerSize(1);
	      m2->Draw();
	    }

	  //Plot 22 Szlata 1977 (Amroun ref 8) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts);i<(Amroun_pts+Collard_pts+Szlata_pts);i++) 
	    {
	      TMarker *m3 = new TMarker(Q2[i], Chi2[i], 20);
	      m3->SetMarkerColor(3);
	      m3->SetMarkerSize(1);
	      m3->Draw();
	    }

	  //Plot 27 Dunn 1983 (Amroun ref 10) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i++) 
	    {
	      TMarker *m4 = new TMarker(Q2[i], Chi2[i], 20);
	      m4->SetMarkerColor(6);
	      m4->SetMarkerSize(1);
	      m4->Draw();
	    }

	  //Plot 16 (skipping 2) JLab data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i++) 
	    {
	      TMarker *m5 = new TMarker(Q2[i], Chi2[i], 20);
	      m5->SetMarkerColor(1);
	      m5->SetMarkerSize(1);
	      m5->Draw();
	    }

	  //Plot 5 Nakagawa 2001 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i++) 
	    {
	      TMarker *m6 = new TMarker(Q2[i], Chi2[i], 20);
	      m6->SetMarkerColor(7);
	      m6->SetMarkerSize(1);
	      m6->Draw();
	    }

	  //Plot my data point.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i++) 
	    {
	      TMarker *m7 = new TMarker(Q2[i], Chi2[i], 20);
	      m7->SetMarkerColor(kOrange+7);
	      m7->SetMarkerSize(1);
	      m7->Draw();
	    }

	  //Plot 11 Arnold 1978 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts);i++) 
	    {
	      TMarker *m8 = new TMarker(Q2[i], Chi2[i], 20);
	      m8->SetMarkerColor(kGreen+2);
	      m8->SetMarkerSize(1);
	      m8->Draw();
	    }

	  TLegend *legend2;
	  legend2 = new TLegend(0.62,0.7,0.9,0.9); //Places legend in upper right corner of histogram.
	  legend2->AddEntry(m2,"Collard 1965","p");
	  legend2->AddEntry(m3,"Szlata 1977","p");
	  legend2->AddEntry(m8,"Arnold 1978","p");
	  legend2->AddEntry(m4,"Dunn 1983","p");
	  legend2->AddEntry(m1,"Amroun 1994","p");
	  legend2->AddEntry(m6,"Nakagawa 2001","p");
	  legend2->AddEntry(m5,"Camsonne 2016","p");
	  legend2->AddEntry(m7,"Current Analysis, Ye 2018","p");
	  legend2->Draw();

	  //Save the canvas as a .png file and a .C file.
	  cQ2->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Chi2_vs_Q2.png");
	  cQ2->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Chi2_vs_Q2.C");
      
	  //Create a plot of XSexp/XSfit.
	  TCanvas* cxsfit_theta=new TCanvas("cxsfit_theta");
	  cxsfit_theta->SetGrid();
      
	  Double_t maxratio = 0.;
	  for(Int_t i=0;i<(nlines-skip);i++)
	    {
	      if(fabs(sigexp[i]/xsfit[i])>maxratio)
		{
		  maxratio = fabs(sigexp[i]/xsfit[i]);
		  cout<<"Max sigexp/sigfit = "<<maxratio<<endl;
		}
	    }
      
	  TH2D *hxsfit = new TH2D("hxsfit","Ratio of Experimental XS to XS from Fit vs. Scattering Angle" , 100, 0., 180., 100, 0., fabs(maxratio)+0.5);
	  for(Int_t i=0;i<(nlines-skip);i++)
	    {
	      //hxsfit->Fill(theta[i],sigexp[i]/xsfit[i]);
	    }
	  hxsfit->SetMarkerStyle(20);
	  hxsfit->SetMarkerSize(1);
	  gStyle->SetOptStat(0);
	  hxsfit->Draw();

	  //Plot 59 Amroun data points. Removed two points without energy listed.
	  for (Int_t i=0;i<Amroun_pts;i++) 
	    {
	      TMarker *m1 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m1->SetMarkerColor(2);
	      m1->SetMarkerSize(1);
	      m1->Draw();
	    }

	  //Plot 118 Collard 1965 (Amroun ref 5) data points.
	  for (Int_t i=Amroun_pts;i<(Amroun_pts+Collard_pts);i++) 
	    {
	      TMarker *m2 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m2->SetMarkerColor(4);
	      m2->SetMarkerSize(1);
	      m2->Draw();
	    }

	  //Plot 22 Szlata 1977 (Amroun ref 8) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts);i<(Amroun_pts+Collard_pts+Szlata_pts);i++) 
	    {
	      TMarker *m3 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m3->SetMarkerColor(3);
	      m3->SetMarkerSize(1);
	      m3->Draw();
	    }

	  //Plot 27 Dunn 1983 (Amroun ref 10) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i++) 
	    {
	      TMarker *m4 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m4->SetMarkerColor(6);
	      m4->SetMarkerSize(1);
	      m4->Draw();
	    }

	  //Plot 16 (skipping 2) JLab data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i++) 
	    {
	      TMarker *m5 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m5->SetMarkerColor(1);
	      m5->SetMarkerSize(1);
	      m5->Draw();
	    }

	  //Plot 5 Nakagawa 2001 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i++) 
	    {
	      TMarker *m6 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m6->SetMarkerColor(7);
	      m6->SetMarkerSize(1);
	      m6->Draw();
	    }

	  //Plot my data point.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i++) 
	    {
	      TMarker *m7 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m7->SetMarkerColor(kOrange+7);
	      m7->SetMarkerSize(1);
	      m7->Draw();
	    }

	  //Plot 11 Arnold 1978 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts);i++) 
	    {
	      TMarker *m8 = new TMarker(theta[i], sigexp[i]/xsfit[i], 20);
	      m8->SetMarkerColor(kGreen+2);
	      m8->SetMarkerSize(1);
	      m8->Draw();
	    }

	  TLegend *legend3;
	  legend3 = new TLegend(0.62,0.7,0.9,0.9); //Places legend in upper right corner of histogram.
	  legend3->AddEntry(m2,"Collard 1965","p");
	  legend3->AddEntry(m3,"Szlata 1977","p");
	  legend3->AddEntry(m8,"Arnold 1978","p");
	  legend3->AddEntry(m4,"Dunn 1983","p");
	  legend3->AddEntry(m1,"Amroun 1994","p");
	  legend3->AddEntry(m6,"Nakagawa 2001","p");
	  legend3->AddEntry(m5,"Camsonne 2016","p");
	  legend3->AddEntry(m7,"Current Analysis, Ye 2018","p");
	  legend3->Draw();

	  //Save the canvas as a .png file and a .C file.
	  cxsfit_theta->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/XSExp_Over_XSFit_vs_Theta.png");
	  cxsfit_theta->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/XSExp_Over_XSFit_vs_Theta.C");

	  //Create a plot of XSexp/XSfit.
	  TCanvas* cxsfit_q=new TCanvas("cxsfit_q");
	  cxsfit_q->SetGrid();

	  //Plot ratio of fit to experiment vs Q^2 instead of theta.
	  TH2D *hxsfitQ2 = new TH2D("hxsfitQ2","Ratio of Experimental XS to XS from Fit vs. q" , 1000, 0., 10., 100, 0., fabs(maxratio)+0.5);
	  for(Int_t i=0;i<(datapts);i++)
	    {
	      //hxsfit->Fill(theta[i],sigexp[i]/xsfit[i]);
	    }
	  hxsfitQ2->SetMarkerStyle(20);
	  hxsfitQ2->SetMarkerSize(1);
	  gStyle->SetOptStat(0);
	  hxsfitQ2->Draw();

	  //Plot 59 Amroun data points.
	  for (Int_t i=0;i<Amroun_pts;i++) 
	    {
	      //TMarker *m1 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m1 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m1->SetMarkerColor(2);
	      m1->SetMarkerSize(1);
	      m1->Draw();
	    }

	  //Plot 118 Collard 1965 (Amroun ref 5) data points.
	  for (Int_t i=Amroun_pts;i<(Amroun_pts+Collard_pts);i++) 
	    {
	      //TMarker *m2 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m2 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m2->SetMarkerColor(4);
	      m2->SetMarkerSize(1);
	      m2->Draw();
	    }

	  //Plot 22 Szlata 1977 (Amroun ref 8) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts);i<(Amroun_pts+Collard_pts+Szlata_pts);i++) 
	    {
	      //TMarker *m3 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m3 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m3->SetMarkerColor(3);
	      m3->SetMarkerSize(1);
	      m3->Draw();
	    }

	  //Plot 27 Dunn 1983 (Amroun ref 10) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i++) 
	    {
	      //TMarker *m4 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m4 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m4->SetMarkerColor(6);
	      m4->SetMarkerSize(1);
	      m4->Draw();
	    }

	  //Plot 16 (skipping 2) JLab data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i++) 
	    {
	      //TMarker *m5 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m5 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m5->SetMarkerColor(1);
	      m5->SetMarkerSize(1);
	      m5->Draw();
	    }

	  //Plot 5 Nakagawa 2001 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i++) 
	    {
	      //TMarker *m6 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m6 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m6->SetMarkerColor(7);
	      m6->SetMarkerSize(1);
	      m6->Draw();
	    }

	  //Plot my data point.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i++) 
	    {
	      //TMarker *m7 = new TMarker(Q2[i], sigexp[i]/xsfit[i], 20);
	      TMarker *m7 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m7->SetMarkerColor(kOrange+7);
	      m7->SetMarkerSize(1);
	      m7->Draw();
	    }

	  //Plot 11 Arnold 1978 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts);i++) 
	    {
	      TMarker *m8 = new TMarker(pow(Q2[i],0.5), sigexp[i]/xsfit[i], 20);
	      m8->SetMarkerColor(kGreen+2);
	      m8->SetMarkerSize(1);
	      m8->Draw();
	    }

	  TLegend *legend4;
	  legend4 = new TLegend(0.62,0.7,0.9,0.9); //Places legend in upper right corner of histogram.
	  legend4->AddEntry(m2,"Collard 1965","p");
	  legend4->AddEntry(m3,"Szlata 1977","p");
	  legend4->AddEntry(m8,"Arnold 1978","p");
	  legend4->AddEntry(m4,"Dunn 1983","p");
	  legend4->AddEntry(m1,"Amroun 1994","p");
	  legend4->AddEntry(m6,"Nakagawa 2001","p");
	  legend4->AddEntry(m5,"Camsonne 2016","p");
	  legend4->AddEntry(m7,"Current Analysis, Ye 2018","p");
	  legend4->Draw();

	  //Save the canvas as a .png file and a .C file.
	  cxsfit_q->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/XSExp_Over_XSFit_vs_Q.png");
	  cxsfit_q->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/XSExp_Over_XSFit_vs_Q.C");

	  //Plot residual of fit to experiment vs Q.
	  TCanvas* cresidual=new TCanvas("cresidual");
	  cresidual->SetGrid();

	  //TH2D *hxsresidualQ2 = new TH2D("hxsresidualQ2","Residual of Experimental XS to XS from SOG Fit vs. q" , 1000, 0., pow(maxQ2,0.5)+0.5, 100, -3., 3.);
	  TH2D *hxsresidualQ2 = new TH2D("hxsresidualQ2","Residual of Representative Fit vs. Q^{2}" , 1000, 0., maxQ2+2., 100, -10., 10.);
	  for(Int_t i=0;i<(datapts);i++)
	    {
	      //hxsresidualQ2->Fill(theta[i],sigexp[i]/xsfit[i]);
	    }
	  hxsresidualQ2->SetMarkerStyle(20);
	  hxsresidualQ2->SetMarkerSize(1);
	  gStyle->SetOptStat(0);
	  hxsresidualQ2->Draw();

	  //Plot 59 Amroun data points.
	  for (Int_t i=0;i<Amroun_pts;i++) 
	    {
	      TMarker *m1 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m1 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m1->SetMarkerColor(2);
	      m1->SetMarkerSize(1);
	      m1->Draw();
	    }

	  //Plot 118 Collard 1965 (Amroun ref 5) data points.
	  for (Int_t i=Amroun_pts;i<(Amroun_pts+Collard_pts);i++) 
	    {
	      TMarker *m2 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m2 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m2->SetMarkerColor(4);
	      m2->SetMarkerSize(1);
	      m2->Draw();
	    }

	  //Plot 22 Szlata 1977 (Amroun ref 8) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts);i<(Amroun_pts+Collard_pts+Szlata_pts);i++) 
	    {
	      TMarker *m3 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m3 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m3->SetMarkerColor(3);
	      m3->SetMarkerSize(1);
	      m3->Draw();
	    }

	  //Plot 27 Dunn 1983 (Amroun ref 10) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i++) 
	    {
	      TMarker *m4 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m4 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m4->SetMarkerColor(6);
	      m4->SetMarkerSize(1);
	      m4->Draw();
	    }

	  //Plot 16 (skipping 2) JLab data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i++) 
	    {
	      TMarker *m5 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m5 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m5->SetMarkerColor(1);
	      m5->SetMarkerSize(1);
	      m5->Draw();
	    }

	  //Plot 5 Nakagawa 2001 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i++) 
	    {
	      TMarker *m6 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m6 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m6->SetMarkerColor(7);
	      m6->SetMarkerSize(1);
	      m6->Draw();
	    }

	  //Plot my data point.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i++) 
	    {
	      TMarker *m7 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m7 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m7->SetMarkerColor(kOrange+7);
	      m7->SetMarkerSize(1);
	      m7->Draw();
	    }

	  //Plot 11 Arnold 1978 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts);i++) 
	    {
	      TMarker *m8 = new TMarker(Q2[i], residual[i], 20);
	      //TMarker *m8 = new TMarker(pow(Q2[i],0.5), residual[i], 20);
	      m8->SetMarkerColor(kGreen+2);
	      m8->SetMarkerSize(1);
	      m8->Draw();
	    }

	  TLegend *legend5;
	  legend5 = new TLegend(0.1,0.7,0.48,0.9); //Places legend in upper left corner of histogram.
	  //legend5 = new TLegend(0.62,0.7,0.9,0.9); //Places legend in upper right corner of histogram.
	  legend5->AddEntry(m2,"Collard 1965","p");
	  legend5->AddEntry(m3,"Szlata 1977","p");
	  legend5->AddEntry(m8,"Arnold 1978","p");
	  legend5->AddEntry(m4,"Dunn 1983","p");
	  legend5->AddEntry(m1,"Amroun 1994","p");
	  legend5->AddEntry(m6,"Nakagawa 2001","p");
	  legend5->AddEntry(m5,"Camsonne 2016","p");
	  legend5->AddEntry(m7,"Current Analysis, Ye 2018","p");
	  legend5->Draw();

	  //Save the canvas as a .png file and a .C file.
	  cresidual->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Residual_vs_Q.png");
	  cresidual->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Residual_vs_Q.C");

	}//End showplots.  
 
      if(fitvars == 0)
	{
	  Qichtot = 0; //Reset to zero in case looping over several fits.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      //Qich[i] = fxs0->GetParameter(i);
	      //Qich[i] = gMinuit->GetParameter(i,1.,1.);
	      gMinuit->GetParameter(i,Qich[i],Qicherr[i]);
	      Qichtot = Qichtot + Qich[i];
	      if(showplots == 1)
		{
		  cout<<"Qich["<<i<<"] = "<<Qich[i]<<"   Qicherr["<<i<<"] = "<<Qicherr[i]<<endl;
		}
	    }
	  if(showplots == 1)
	    {
	      cout<<"Qichtot = "<<Qichtot<<endl;
	    }
	  /*
	    for(Int_t i=0;i<ngaus;i++)
	    {
	    Qich[i] = Qich[i]/Qichtot;
	    cout<<"Qich["<<i<<"] = "<<Qich[i]<<endl;
	    }*/

	  Qimtot = 0; //Reset to zero in case looping over several fits.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      //Qim[i] = fxs0->GetParameter(ngaus+i);
	      gMinuit->GetParameter(ngaus+i,Qim[i],Qimerr[i]);
	      Qimtot = Qimtot + Qim[i];
	      if(showplots == 1)
		{
		  cout<<"Qim["<<i<<"] = "<<Qim[i]<<"   Qimerr["<<i<<"] = "<<Qimerr[i]<<endl;
		}
	    }
	  if(showplots == 1)
	    {
	      cout<<"Qimtot = "<<Qimtot<<endl;
	    }     

	  if(showplots == 1)
	    {
	      for(Int_t i=0;i<ngaus;i++)
		{
		  cout<<"R["<<i<<"] = "<<R[i]<<endl;
		}
	 
	      cout<<"Gamma = "<<Gamma<<endl;
	    }
	}

      //Overrides fitted Qi. 
      if(Amroun_Qi == 1)
	{
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      Qich[i] = Qich_Amroun[i];
	      Qim[i] = Qim_Amroun[i];
	    }
	}


      //Fill text file with fcn (chi2), Qichtot, Qimtot, Ri, Qich, Qim.
      output<<amin<<" "<<amin/(datapts-2*ngaus-1)<<" "<<datapts*TMath::Log(amin/datapts)+TMath::Log(datapts)*ngaus<<" "<<datapts*TMath::Log(amin/datapts)+2*ngaus<<" "<<Qichtot<<" "<<Qimtot<<" "<<R[0]<<" "<<R[1]<<" "<<R[2]<<" "<<R[3]<<" "<<R[4]<<" "<<R[5]<<" "<<R[6]<<" "<<R[7]<<" "<<R[8]<<" "<<R[9]<<" "<<R[10]<<" "<<R[11]<<" "<<Qich[0]<<" "<<Qich[1]<<" "<<Qich[2]<<" "<<Qich[3]<<" "<<Qich[4]<<" "<<Qich[5]<<" "<<Qich[6]<<" "<<Qich[7]<<" "<<Qich[8]<<" "<<Qich[9]<<" "<<Qich[10]<<" "<<Qich[11]<<" "<<Qim[0]<<" "<<Qim[1]<<" "<<Qim[2]<<" "<<Qim[3]<<" "<<Qim[4]<<" "<<Qim[5]<<" "<<Qim[6]<<" "<<Qim[7]<<" "<<Qim[8]<<" "<<Qim[9]<<" "<<Qim[10]<<" "<<Qim[11]<<endl;


    }//End loop over minimization for different Ri values or bootstrapping.
  //Close output file.
  output.close();

  //print_fit();
  /*
  if(showgaus == 1)
    {
      //Plot individual Gaussians with their fit parameters. 
      TF1 *g0 = new TF1("g0", fitg, yminFF, ymaxFF+20,2.);
      g0->SetParameters(Qich[0],R[0]);
      g0->SetLineColor(1);
      g0->SetNpx(npdraw);
      g0->Draw("csame");
      TF1 *g1 = new TF1("g1", fitg, yminFF, ymaxFF+20,2.);
      g1->SetParameters(Qich[1],R[1]);
      g1->SetLineColor(2);
      g1->SetNpx(npdraw);
      g1->Draw("cSame");
      TF1 *g2 = new TF1("g2", fitg, yminFF, ymaxFF+20,2.);
      g2->SetParameters(Qich[2],R[2]);
      g2->SetLineColor(3);
      g2->SetNpx(npdraw);
      g2->Draw("cSame");
      TF1 *g3 = new TF1("g3", fitg, yminFF, ymaxFF+20,2.);
      g3->SetParameters(Qich[3],R[3]);
      g3->SetLineColor(4);
      g3->SetNpx(npdraw);
      g3->Draw("cSame");
      TF1 *g4 = new TF1("g4", fitg, yminFF, ymaxFF+20,2.);
      g4->SetParameters(Qich[4],R[4]);
      g4->SetLineColor(5);
      g4->SetNpx(npdraw);
      g4->Draw("cSame");
      TF1 *g5 = new TF1("g5", fitg, yminFF, ymaxFF+20,2.);
      g5->SetParameters(Qich[5],R[5]);
      g5->SetLineColor(6);
      g5->SetNpx(npdraw);
      g5->Draw("cSame");
      TF1 *g6 = new TF1("g6", fitg, yminFF, ymaxFF+20,2.);
      g6->SetParameters(Qich[6],R[6]);
      g6->SetLineColor(7);
      g6->SetNpx(npdraw);
      g6->Draw("cSame");
      TF1 *g7 = new TF1("g7", fitg, yminFF, ymaxFF+20,2.);
      g7->SetParameters(Qich[7],R[7]);
      g7->SetLineColor(8);
      g7->SetNpx(npdraw);
      g7->Draw("cSame");
      TF1 *g8 = new TF1("g8", fitg, yminFF, ymaxFF+20,2.);
      g8->SetParameters(Qich[8],R[8]);
      g8->SetLineColor(9);
      g8->SetNpx(npdraw);
      g8->Draw("cSame");
      TF1 *g9 = new TF1("g9", fitg, yminFF, ymaxFF+20,2.);
      g9->SetParameters(Qich[9],R[9]);
      g9->SetLineColor(46);
      g9->SetNpx(npdraw);
      g9->Draw("cSame");
      TF1 *g10 = new TF1("g10", fitg, yminFF, ymaxFF+20,2.);
      g10->SetParameters(Qich[10],R[10]);
      g10->SetLineColor(11);
      g10->SetNpx(npdraw);
      g10->Draw("cSame");
      TF1 *g11 = new TF1("g11", fitg, yminFF, ymaxFF+20,2.);
      g11->SetParameters(Qich[11],R[11]);
      g11->SetLineColor(12);
      g11->SetNpx(npdraw);
      g11->Draw("cSame");

      cout<<"Integral g0 = "<<g0->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g1 = "<<g1->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g2 = "<<g2->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g3 = "<<g3->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g4 = "<<g4->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g5 = "<<g5->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g6 = "<<g6->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g7 = "<<g7->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g8 = "<<g8->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g9 = "<<g9->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g10 = "<<g10->Integral(yminFF,ymaxFF)<<endl;
      cout<<"Integral g11 = "<<g11->Integral(yminFF,ymaxFF)<<endl;
    }
  */
  if(showplots == 1)
    {
      TCanvas* c2=new TCanvas("c2");
      c2->SetGrid();
      c2->SetLogy();
      c2->SetTitle("Charge Form Factor");

      TF1 *fChFF = new TF1("fChFF",ChFF_Q2,yminFF,ymaxFF+54,1);
      //TF1 *fChFF = new TF1("fChFF",ChFF,yminFF,ymaxFF,1);
      //cout<<fChFF->Eval(0.000001)<<"!!!!!!!!"<<endl;
      cout<<"Qichtot = "<<Qichtot<<endl;
      fChFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
      fChFF->Draw("L");
      //c2->SetTitle("Charge Form Factor");
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

      //Save the canvas as a .png file and a .C file.
      c2->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/ChFF.png");
      c2->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/ChFF.C");

      //Plot the magnetic FF using the parameters determined by the SOG fit above. 
      //Make a new canvas to plot data.
      TCanvas* c4=new TCanvas("c4");
      c4->SetGrid();
      c4->SetLogy();
     
      TF1 *fMFF = new TF1("fMFF",MFF_Q2,yminFF,ymaxFF+54,1);
      //cout<<fMFF->Eval(0.000001)<<"!!!!!!!!"<<endl;
      cout<<"Qimtot = "<<Qimtot<<endl;
      fMFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
      fMFF->Draw("L");
      c4->SetTitle("He3 Magnetic Form Factor");
      //fChFF->SetTitle("C12 Charge Form Factor","#Q^2 (#fm^-2)","#F_{Ch}(q)");
      //fMFF->GetHistogram()->SetTitle("^{3}He Magnetic Form Factor");
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

      //Being used to test plotting FF world data. Moved to Plot_FFs.C since it's faster.
      /*
      Float_t xtest1[1],ytest1[1],xtest2[1],ytest2[1];
      xtest1[0] = 30.0;
      ytest1[0] = 0.0001;
      xtest2[0] = 40.0;
      ytest2[0] = 0.0001;
      
      TGraph *gr1 = new TGraph (1, xtest1, ytest1); 
      gr1->SetMarkerColorAlpha(kBlue, 0.35);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerSize(5);
      //gr1->SetLineColorAlpha(kBlue, 0.35);
      //gr1->SetLineWidth(2);
      gr1->Draw("same p");

      TGraph *gr2 = new TGraph (1, xtest2, ytest2); 
      gr2->SetMarkerColorAlpha(kRed, 0.35);
      gr2->SetMarkerStyle(20);
      gr2->SetMarkerSize(5);
      gr2->Draw("same p");

      TLegend *test_leg;
      test_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
      test_leg->AddEntry("fMFF","New ^{3}He |F_{m}(Q^{2})| Fit","l");
      test_leg->AddEntry(gr1,"Test 1","p");
      test_leg->AddEntry(gr2,"Test 2","p");
      test_leg->Draw();
      */

      /*
      TGraph *gr1 = new TGraph (2);
      for(Int_t i=0;i<1;i++)
	{
	  gr1->SetPoint(i,xtest[0],ytest[0]);
	}
      for(Int_t i=1;i<2;i++)
	{
	  gr1->SetPoint(i,xtest[1],ytest[1]);
	}
      */

      //Save the canvas as a .png file and a .C file.
      c4->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/MFF.png");
      c4->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/MFF.C");
     
      cout<<"Chi^2 (amin) = "<<amin<<"   Reduced Chi^2 = "<<amin<<"/("<<datapts<<" - "<<2*ngaus<<" - 1) = "<<amin/(datapts-2*ngaus-1)<<endl;
      cout<<"BIC = "<<datapts<<"*ln("<<amin<<"/"<<datapts<<")+"<<ngaus<<"*ln("<<datapts<<") = "<<datapts*TMath::Log(amin/datapts)+TMath::Log(datapts)*ngaus<<"   AIC = "<<datapts<<"*ln("<<amin<<"/"<<datapts<<")+"<<2*ngaus<<" = "<<datapts*TMath::Log(amin/datapts)+2*ngaus<<endl;

      //Now draw both FFs on the same plot for the GRC poster. 
      TCanvas* c5=new TCanvas("c5");
      c5->SetGrid();
      /*
      c5->Divide(1,3);
      c5->cd(1)->SetLogy();
      c5->cd(1);
      fChFF->Draw("L");
      */
      //cout<<fChFF->Eval(35)<<endl;
      //Now add Amroun's fit to the plot for comparison.
      /*for(Int_t i=0;i<ngaus;i++)
	{
	Qich[i] = Qich_Amroun[i];
	}*/
      //Qich[0] = 2;
      TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,yminFF,ymaxFF+54,1);
      //cout<<fChFF_Amroun->Eval(35)<<endl;
      fChFF_Amroun->SetNpx(npdraw);
      fChFF_Amroun->SetLineColor(4);
      //fChFF_Amroun->Draw("L same");
      TLegend *ChFF_leg;
      ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
      ChFF_leg->AddEntry("fChFF","New ^{3}He |F_{ch}(Q^{2})| Fit","l");
      ChFF_leg->AddEntry("fChFF_Amroun","Fit from Amroun et al.","l");
      //ChFF_leg->Draw();
     
      /*
      c5->cd(2);
      c5->cd(2)->SetLogy();
      //cout<<fMFF->Eval(35)<<endl;
      fMFF->Draw("L");
      */     

      //Now add Amroun's fit to the plot for comparison.
      /*for(Int_t i=0;i<ngaus;i++)
	{
	Qim[i] = Qim_Amroun[i];
	}*/
      TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,ymaxFF+54,1);
      //cout<<fMFF_Amroun->Eval(30)<<endl;
      fMFF_Amroun->SetNpx(npdraw);
      fMFF_Amroun->SetLineColor(4);
      //fMFF_Amroun->Draw("L same");
      TLegend *MFF_leg;
      MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
      MFF_leg->AddEntry("fMFF","New ^{3}He |F_{m}(Q^{2})| Fit","l");
      MFF_leg->AddEntry("fMFF_Amroun","Fit from Amroun et al.","l");
      //MFF_leg->Draw();

      //c5->cd(3);
      TH1F *h1 = new TH1F("h1", "^{3}He Cross Section World Data Distribution", 70, yminFF, ymaxFF+64);
      if(bootstrap == 0)
	{
	  for(Int_t i=0;i<datapts;i++)
	    {
	      h1->Fill(Q2[i]);
	    }
	}
      if(bootstrap == 1)
	{
	  for(Int_t i=0;i<datapts;i++)
	    {
	      h1->Fill(Q2_bs[i]);
	    }
	}
      h1->GetXaxis()->SetTitle("Q^{2} (fm^{-2})");
      h1->GetXaxis()->CenterTitle(true);
      h1->GetXaxis()->SetLabelSize(0.05);
      h1->GetXaxis()->SetTitleSize(0.06);
      h1->GetXaxis()->SetTitleOffset(0.75);
      h1->GetYaxis()->SetTitle("Number of Measurements");
      h1->GetYaxis()->CenterTitle(true);
      h1->GetYaxis()->SetLabelSize(0.05);
      h1->GetYaxis()->SetTitleSize(0.06);
      h1->GetYaxis()->SetTitleOffset(0.75);
      h1->SetFillColor(4);
      h1->Draw();

      //Save the canvas as a .png file and a .C file.
      c5->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FF_Comparison.png");
      c5->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FF_Comparison.C");
    }//End showplots.

  //Plot the charge density from I. Sick. 
  TCanvas* crho=new TCanvas("crho");
  crho->SetGrid();
 
  TF1 *frho_ch = new TF1("frho_ch",rho_ch,0.,4.,1);
  frho_ch->SetNpx(npdraw);
  frho_ch->SetTitle("^{3}He Charge Charge Density");
  frho_ch->GetHistogram()->GetXaxis()->SetTitle("r (fm)");
  frho_ch->GetHistogram()->GetXaxis()->CenterTitle(true);
  frho_ch->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  frho_ch->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  frho_ch->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);
  frho_ch->GetHistogram()->GetYaxis()->SetTitle("e/fm^{3}");
  frho_ch->GetHistogram()->GetYaxis()->CenterTitle(true);
  frho_ch->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  frho_ch->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  frho_ch->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  frho_ch->SetLineColor(2);
  frho_ch->Draw();

  if(showgaus == 1)
    {
      //Plot individual Gaussians with their fit parameters. 
      TF1 *g0 = new TF1("g0", fitg_rho, yminFF, ymaxFF+20,2.);
      g0->SetParameters(Qich[0],R[0]);
      g0->SetLineColor(1);
      g0->SetNpx(npdraw);
      g0->Draw("csame");
      TF1 *g1 = new TF1("g1", fitg_rho, yminFF, ymaxFF+20,2.);
      g1->SetParameters(Qich[1],R[1]);
      g1->SetLineColor(30);
      g1->SetNpx(npdraw);
      g1->Draw("cSame");
      TF1 *g2 = new TF1("g2", fitg_rho, yminFF, ymaxFF+20,2.);
      g2->SetParameters(Qich[2],R[2]);
      g2->SetLineColor(3);
      g2->SetNpx(npdraw);
      g2->Draw("cSame");
      TF1 *g3 = new TF1("g3", fitg_rho, yminFF, ymaxFF+20,2.);
      g3->SetParameters(Qich[3],R[3]);
      g3->SetLineColor(4);
      g3->SetNpx(npdraw);
      g3->Draw("cSame");
      TF1 *g4 = new TF1("g4", fitg_rho, yminFF, ymaxFF+20,2.);
      g4->SetParameters(Qich[4],R[4]);
      g4->SetLineColor(5);
      g4->SetNpx(npdraw);
      g4->Draw("cSame");
      TF1 *g5 = new TF1("g5", fitg_rho, yminFF, ymaxFF+20,2.);
      g5->SetParameters(Qich[5],R[5]);
      g5->SetLineColor(6);
      g5->SetNpx(npdraw);
      g5->Draw("cSame");
      TF1 *g6 = new TF1("g6", fitg_rho, yminFF, ymaxFF+20,2.);
      g6->SetParameters(Qich[6],R[6]);
      g6->SetLineColor(7);
      g6->SetNpx(npdraw);
      g6->Draw("cSame");
      TF1 *g7 = new TF1("g7", fitg_rho, yminFF, ymaxFF+20,2.);
      g7->SetParameters(Qich[7],R[7]);
      g7->SetLineColor(8);
      g7->SetNpx(npdraw);
      g7->Draw("cSame");
      TF1 *g8 = new TF1("g8", fitg_rho, yminFF, ymaxFF+20,2.);
      g8->SetParameters(Qich[8],R[8]);
      g8->SetLineColor(9);
      g8->SetNpx(npdraw);
      g8->Draw("cSame");
      TF1 *g9 = new TF1("g9", fitg_rho, yminFF, ymaxFF+20,2.);
      g9->SetParameters(Qich[9],R[9]);
      g9->SetLineColor(46);
      g9->SetNpx(npdraw);
      g9->Draw("cSame");
      TF1 *g10 = new TF1("g10", fitg_rho, yminFF, ymaxFF+20,2.);
      g10->SetParameters(Qich[10],R[10]);
      g10->SetLineColor(11);
      g10->SetNpx(npdraw);
      g10->Draw("cSame");
      TF1 *g11 = new TF1("g11", fitg_rho, yminFF, ymaxFF+20,2.);
      g11->SetParameters(Qich[11],R[11]);
      g11->SetLineColor(12);
      g11->SetNpx(npdraw);
      g11->Draw("cSame");

      Double_t int_rho = frho_ch->Integral(yminFF,ymaxFF); //Integral of rho (no normalization).
      cout<<"Integral Rho = "<<int_rho<<endl;
      cout<<"Integral g0 = "<<g0->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g0->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g1 = "<<g1->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g1->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g2 = "<<g2->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g2->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g3 = "<<g3->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g3->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g4 = "<<g4->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g4->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g5 = "<<g5->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g5->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g6 = "<<g6->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g6->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g7 = "<<g7->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g7->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g8 = "<<g8->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g8->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g9 = "<<g9->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g9->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g10 = "<<g10->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g10->Integral(yminFF,ymaxFF)/int_rho<<endl;
      cout<<"Integral g11 = "<<g11->Integral(yminFF,ymaxFF)<<"   Fraction of charge density held by this Gaussian = "<<g11->Integral(yminFF,ymaxFF)/int_rho<<endl;

    }

  TLegend *rho_leg;
  rho_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  rho_leg->AddEntry("frho_ch","New ^{3}He Charge Density","l");
  rho_leg->AddEntry("g0","Gaussian 1","l");
  rho_leg->AddEntry("g1","Gaussian 2","l");
  rho_leg->AddEntry("g2","Gaussian 3","l");
  rho_leg->AddEntry("g3","Gaussian 4","l");
  rho_leg->AddEntry("g4","Gaussian 5","l");
  rho_leg->AddEntry("g5","Gaussian 6","l");
  rho_leg->AddEntry("g6","Gaussian 7","l");
  rho_leg->AddEntry("g7","Gaussian 8","l");
  rho_leg->AddEntry("g8","Gaussian 9","l");
  rho_leg->AddEntry("g9","Gaussian 10","l");
  rho_leg->AddEntry("g10","Gaussian 11","l");
  rho_leg->AddEntry("g11","Gaussian 12","l");
  rho_leg->Draw();

  TF1 *frho_ch_int = new TF1("frho_ch_int",rho_ch_int,0.,4.,1);
  frho_ch_int->SetNpx(npdraw);
  frho_ch_int->SetLineColor(4);
  //frho_ch_int->Draw("same");

  //Save the canvas as a .png file and a .C file.
  crho->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Charge_Density.png");
  crho->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/Charge_Density.C");

  TF1 *frho_rms = new TF1("frho_rms",rho_rms,0.,4.,1);

  cout<<"Integral of frho_ch_int = "<<frho_ch_int->Integral(0.0,10.)<<"e."<<endl;
  //Calculate rms radius of nucleus. rms radius = int(rho(r) * r^2 dr).
  //cout<<"rms radius = "<<frho_ch_int->Integral(0.0,10.)*e/(4*pi)<<endl;
  cout<<"Integral of frho_ch = "<<frho_ch->Integral(0.0,10.)<<"e."<<endl;
  //cout<<"rms radius = "<<pow(frho_rms->Integral(0.0,10.)*e,0.5)<<endl;
  //cout<<"rms radius no division = "<<pow( frho_rms->Integral(0.0,10.) ,0.5)<<endl; //I. Sick.
  cout<<"rms radius = "<<pow( frho_rms->Integral(0.0,10.)/frho_ch_int->Integral(0.0,10.) ,0.5)<<endl; //I think this is the correct rms formulation.
  //cout<<""<<frho_ch->Integral(0.,10.)<<endl;
   
  //Define ChFF_Deriv to be the charge FF and then find its derivative at 0 to get charge radius,
  cout<<"ChFF_Deriv(~0) = "<<ChFF_Deriv(0.0000001)<<endl;
  double x0 = 0.002;
  ROOT::Math::Functor1D f1D(&ChFF_Deriv);
 
  ROOT::Math::RichardsonDerivator rd;
  rd.SetFunction(f1D);
  cout<<"First Derivative:   "<<rd.Derivative1(x0)<<"   -6*dFC(0)/dQ^2 = "<<-6*rd.Derivative1(x0)<<"   rms radius = "<<pow(-6*rd.Derivative1(x0),0.5)<<endl;
  //std::cout << "Second Derivative:  " << rd.Derivative2(x0) << std::endl;
  //std::cout << "Third Derivative:   " << rd.Derivative3(x0) << std::endl;
 

  if(useFB == 1)
    {
      //Initiate Minuit for minimization of Fourrier-Bessel fit.
      gSystem->Load("libMathMore");            //Needed to use cyl_bessel_j() function.
      TMinuit *gMinuit_FB = new TMinuit(2*nFB);  //initialize TMinuit with a maximum number of parameters.
      gMinuit_FB->SetFCN(fcn_FB);
     
      Double_t arglist_FB[10];
      Int_t ierflg_FB = 0;
     
      arglist_FB[0] = 1.;
      gMinuit_FB->mnexcm("SET ERR", arglist_FB ,1,ierflg_FB);
     
      //Set step sizes.
      static Double_t stepsize_FB[4] = {0.1 , 0.1 , 0.01 , 0.001};
     
      //Set starting guesses for parameters. (Use Retzlaff's FB parameters.)
      //FC
      for(Int_t i=0;i<nFB;i++)
	{
	  gMinuit_FB->mnparm(i, Form("av%d",i+1), av[i], stepsize_FB[3], 0.,1.,ierflg_FB);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.001,Qich[i]+0.001,ierflg);
	}
      if(useFB_FM == 1)
	{
	  //FM (Guess that Retzlaff's FC FB parameters are close enough.) Doesn't seem to really work.
	  for(Int_t i=nFB;i<2*nFB;i++)
	    {
	      gMinuit_FB->mnparm(i, Form("av%d",i+1), av[i-nFB], stepsize_FB[3], 0.,1.,ierflg_FB);
	      //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.001,Qich[i]+0.001,ierflg);
	    }
	}
     
      // Now ready for minimization step
      arglist_FB[0] = 500.;//50000.
      arglist_FB[1] = 1.;
      //cout<<"Sup1"<<endl;
      gMinuit_FB->mnexcm("MIGRAD", arglist_FB , 2, ierflg_FB);
      //cout<<"Sup2"<<endl;
      // Print results
      Double_t amin_FB,edm_FB,errdef_FB;
      Int_t nvpar_FB,nparx_FB,icstat_FB;
      gMinuit_FB->mnstat(amin_FB, edm_FB, errdef_FB, nvpar_FB, nparx_FB, icstat_FB);
      //gMinuit_FB->mnprin(3,amin_FB);
     
      //Set av and averr to the fitted parameters.
      if(useFB_FM == 0)
	{
	  for(Int_t i=0;i<nFB;i++)
	    {
	      gMinuit_FB->GetParameter(i,av[i],averr[i]);
	    }
	}
      if(useFB_FM == 1)
	{
	  for(Int_t i=0;i<2*nFB;i++)
	    {
	      gMinuit_FB->GetParameter(i,av[i],averr[i]);
	    }
	}

      if(showplots == 1)
	{ 
	  for(Int_t i=0;i<(nlines-skip);i++)
	    {
	      cout<<"Chi2_FB["<<i<<"] = "<<Chi2_FB[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   FB(E0,theta,par)["<<i<<"] = "<<FBfit[i]<<"   XSexp/XSfit = "<<sigexp[i]/FBfit[i]<<"   residual_FB["<<i<<"] = "<<residual_FB[i]<<endl;
	    }

	  //Plot residual of FB fit to experiment vs Q.
	  TCanvas* cresidual_FB=new TCanvas("cresidual_FB");
	  cresidual_FB->SetGrid();
	 
	  TH2D *hxsresidualQ2_FB = new TH2D("hxsresidualQ2_FB","Residual of Experimental XS to XS from FB Fit vs. q" , 1000, 0., pow(maxQ2,0.5)+0.5, 100, 0., 12.);
	  for(Int_t i=0;i<(datapts);i++)
	    {
	      //hxsresidualQ2->Fill(theta[i],sigexp[i]/xsfit[i]);
	    }
	  hxsresidualQ2_FB->SetMarkerStyle(20);
	  hxsresidualQ2_FB->SetMarkerSize(1);
	  gStyle->SetOptStat(0);
	  hxsresidualQ2_FB->Draw();

	  //Plot 59 Amroun data points.
	  for (Int_t i=0;i<Amroun_pts;i++) 
	    {
	      //TMarker *m1 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m1 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m1->SetMarkerColor(2);
	      m1->SetMarkerSize(1);
	      m1->Draw();
	    }
	 
	  //Plot 118 Collard 1965 (Amroun ref 5) data points.
	  for (Int_t i=Amroun_pts;i<(Amroun_pts+Collard_pts);i++) 
	    {
	      //TMarker *m2 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m2 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m2->SetMarkerColor(4);
	      m2->SetMarkerSize(1);
	      m2->Draw();
	    }
	 
	  //Plot 22 Szlata 1977 (Amroun ref 8) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts);i<(Amroun_pts+Collard_pts+Szlata_pts);i++) 
	    {
	      //TMarker *m3 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m3 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m3->SetMarkerColor(3);
	      m3->SetMarkerSize(1);
	      m3->Draw();
	    }
	 
	  //Plot 27 Dunn 1983 (Amroun ref 10) data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i++) 
	    {
	      //TMarker *m4 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m4 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m4->SetMarkerColor(6);
	      m4->SetMarkerSize(1);
	      m4->Draw();
	    }
	 
	  //Plot 16 (skipping 2) JLab data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i++) 
	    {
	      //TMarker *m5 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m5 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m5->SetMarkerColor(1);
	      m5->SetMarkerSize(1);
	      m5->Draw();
	    }
	 
	  //Plot 5 Nakagawa 2001 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i++) 
	    {
	      //TMarker *m6 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m6 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m6->SetMarkerColor(7);
	      m6->SetMarkerSize(1);
	      m6->Draw();
	    }
	 
	  //Plot my data point.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i++) 
	    {
	      //TMarker *m7 = new TMarker(Q2[i], residual[i], 20);
	      TMarker *m7 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m7->SetMarkerColor(kOrange+7);
	      m7->SetMarkerSize(1);
	      m7->Draw();
	    }

	  //Plot 11 Arnold 1978 data points.
	  for (Int_t i=(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts);i<(Amroun_pts+Collard_pts+Szlata_pts+Dunn_pts+Camsonne_pts+Nakagawa_pts+my_pts+Arnold_pts);i++) 
	    {
	      TMarker *m8 = new TMarker(pow(Q2[i],0.5), residual_FB[i], 20);
	      m8->SetMarkerColor(kGreen+2);
	      m8->SetMarkerSize(1);
	      m8->Draw();
	    }
	 
	  TLegend *legend6;
	  legend6 = new TLegend(0.62,0.7,0.9,0.9); //Places legend in upper right corner of histogram.
	  legend6->AddEntry(m2,"Collard 1965","p");
	  legend6->AddEntry(m3,"Szlata 1977","p");
	  legend6->AddEntry(m8,"Arnold 1978","p");
	  legend6->AddEntry(m4,"Dunn 1983","p");
	  legend6->AddEntry(m1,"Amroun 1994","p");
	  legend6->AddEntry(m6,"Nakagawa 2001","p");
	  legend6->AddEntry(m5,"Camsonne 2016","p");
	  legend6->AddEntry(m7,"Current Analysis, Ye 2018","p");
	  legend6->Draw();

	  //Save the canvas as a .png file and a .C file.
	  cresidual_FB->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_Residual_vs_Q.png");
	  cresidual_FB->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_Residual_vs_Q.C");

	  //Plot the FB FC. 
	  TCanvas* cFB_FC=new TCanvas("cFB_FC");
	  cFB_FC->SetGrid();
	  cFB_FC->SetLogy();

	  TF1 *fFB_FC = new TF1("fFB_FC",FB_FC,yminFF,ymaxFF+54,1);
	  fFB_FC->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
	  fFB_FC->Draw("L");

	  //Save the canvas as a .png file and a .C file.
	  cFB_FC->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_FC.png");
	  cFB_FC->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_FC.C");

	  if(useFB_FM == 1)
	    {
	      //Plot the FB FM. 
	      TCanvas* cFB_FM=new TCanvas("cFB_FM");
	      cFB_FM->SetGrid();
	      cFB_FM->SetLogy();
	     
	      TF1 *fFB_FM = new TF1("fFB_FM",FB_FM,yminFF,ymaxFF+54,1);
	      fFB_FM->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
	      fFB_FM->Draw("L");
	      //Save the canvas as a .png file and a .C file.
	      cFB_FM->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_FM.png");
	      cFB_FM->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_FM.C");
	    }

	  //Plot the FB charge density. 
	  TCanvas* cFB_rho=new TCanvas("cFB_rho");
	  cFB_rho->SetGrid();
	  
	  TF1 *fFB_rho = new TF1("fFB_rho",rho_ch_FB,yminFF,6.,1);
	  fFB_rho->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
	  fFB_rho->Draw("L");
	  //Save the canvas as a .png file and a .C file.
	  cFB_rho->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_rho.png");
	  cFB_rho->SaveAs("/home/skbarcus/Tritium/Analysis/SOG/Output/FB_rho.C");

	  cout<<"Integral of rho = "<<fFB_rho->Integral(0.0,10.)<<" e."<<endl;

	}//End showplots.

      //Calculate FB rms charge radius.
      TF1 *fFB_rho_ch_int = new TF1("fFB_rho_ch_int",FB_rho_ch_int,0.,10.,1);
      TF1 *fFB_rho_rms = new TF1("fFB_rho_rms",FB_rho_rms,0.,10.,1);
      //cout<<"FB rms charge radius = "<<pow( fFB_rho_rms->Integral(0.0,10.) ,0.5)<<endl;
      cout<<"FB rms charge radius = "<<pow( fFB_rho_rms->Integral(0.0,10.)/fFB_rho_ch_int->Integral(0.0,10.) ,0.5)<<endl; //I think this is the correct rms formulation.
    }//End use FB.
 
  st->Stop();
  cout<<"*********************************************"<<endl;
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}

int main()
{
  //Double_t xyz = 5.5;
  Global_Fit_3He_SOG();
}
