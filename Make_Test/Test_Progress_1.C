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
Double_t e = 1.60217662E-19;             //Electron charge C.
Double_t alpha = 0.0072973525664;//1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.

Int_t loops = 1;
const Int_t datapts = 246+11;//248
Int_t userand = 10;                       //0 = use predetermined Ri from Amroun. 1 = use random Ri in generated in a range around Amroun's. 2 = use random Ri, ngaus=12, generated in increments of 0.1 with larger possible spacing at greater radii. 3 = use predetermined Ri for the purposes of trying to tune the fit by hand. 4 = ngaus=8. 5 = ngaus=9. 6 = ngaus=10. 7 = ngaus=11.
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 0;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t fft = 0;                           //0 = don't use FFT to try to get a charge radii. 1 = do use FFT to extract a charge radii.
Int_t Amroun_Qi = 0;                     //1 = Override fitted Qi and use Amroun's values.
Int_t showplots = 1;
Int_t useFB = 0;                         //Turn on Fourier Bessel fit.
Int_t useFB_GM = 1;                      //0 = Turn on Fourier Bessel fit just for GE. 1 = Turn on Fourier Bessel fit attempting GE and GM.
Int_t improve = 0;                       //1 = use mnimpr() to check for other minima around the one MIGRAD finds.
Int_t MINOS = 0;                         //1 = use MINOS to calculate parameter errors. With ERRordef=30, npar=24, 10000 calls took about 1.5 hours and gave results only slightly different from intial parameter errors given. Several pars were hitting limits. 
Int_t optimize_Ri = 0;                   //1 = Have code loop over each Ri value shifting it 0.1 higher and 0.1 lower until chi2 stops improving.  
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
Int_t ngaus_Amroun = 12;
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
Double_t Gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
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
Float_t Q2[datapts];

Int_t Amroun_pts = 57;                 //Dropped two points with no energy value give.
Int_t Collard_pts = 118;
Int_t Szlata_pts = 22;
Int_t Dunn_pts = 27;
Int_t Camsonne_pts = 16;               //Dropped two points with crazy Chi^2 values. Should reevaluate eventually.
Int_t Nakagawa_pts = 5;
Int_t my_pts = 1;
Int_t Arnold_pts = 11;                 //These XSs had to be calculated from A^1/2 function.

Double_t m = 2.;
//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
//Double_t R[12] = {0.2,0.7,1.3,1.5,2.1,2.8,3.6,4.2,5.2,0.,0.,0.};
//Double_t R[12] = {0.1,0.6,1.,1.5,2.1,2.4,3.,3.7,4.4,4.7,0.,0.};//9/12/18 pretty good 2.
//Double_t R[12] = {0.1,0.6,1.,1.5,2.0,2.4,3.1,3.9,4.5,4.8,0.,0.}; //9/12/18 pretty good 1. Also these Ri work well with Qi 2. Very smooth Fm.
Double_t R[12] = {0.1,0.5,0.9,1.3,1.5,1.9,2.3,2.8,3.1,3.8,4.3,5.};//Best n=12 from Amroun starting.
//Double_t R[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_init[12] = {};
Double_t R_best[12] = {};
Double_t R_best_chi2 = 0;
Double_t Qich[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};//Amroun
Double_t Qim[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
//Double_t Qich[12] = {0.0289116,0.176012,0.227652,0.18408,0.186425,0.093576,0.0329847,0.052855,0.016552,0.00285043,0.00615456,1.58614E-11};
//Double_t Qim[12] = {0.0585454,0.160715,0.222426,0.156211,0.191486,0.125172,0.00162158,0.0602476,0.0228372,6.94389E-13,0.0192538,1.03119E-11};

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
//Double_t Qich[12] = {0.0411639,0.236983,0.204183,0.276745,0.0977119,0.0724304,0.0541865,0.0207177,1.26197e-07,0.00411213,0.,0.};//9/12/18 pretty good 2.
//Double_t Qim[12] = {0.07112,0.213554,0.205701,0.238441,0.162279,0.017472,0.0869983,0.0239365,0.0295892,5.58285e-10,0.,0.};
//Double_t Qich[12] = {0.0411535,0.237047,0.203676,0.277975,0.092553,0.0821181,0.0562533,0.0157997,2.61438e-10,0.00176867,0.,0.};//9/12/18 pretty good 1.
//Double_t Qim[12] = {0.0692335,0.221858,0.189445,0.254013,0.127296,0.0508983,0.0726007,0.0160226,0.0224739,4.38656e-10,0.,0.}

Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
Double_t Qich_best[12] = {};
Double_t Qim_best[12] = {};
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

  val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + (pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + pow(tan(theta*deg2rad/2),2))*pow(fitm,2.) ); //magnetic moment for C12 is 0 -> no mag part of XS.
  //cout<<"XS = "<<val<<endl;
  return val;
}

Double_t FB(float E0, float theta, Double_t *par)
{
  Double_t val = 0.;

  Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
  Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
  Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
  Double_t FB_GE_sum = 0.;
  Double_t FB_GE_temp = 0.;
  Double_t FB_GM_sum = 0.;
  Double_t FB_GM_temp = 0.;
  Double_t R_FB = 5.;  //fm
  Double_t mottxs = 0.;
  Double_t tau = 0;
  Double_t rho = E0/MtHe3;

  //Calculate Mott XS.
  mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/GeV2fm;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  //cout<<"MottXS = "<<mottxs<<endl;
  //Calculate tau.
  tau = Q2eff/(4*pow(MtHe3,2.)*GeV2fm);

  //Calculate GE.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_GE_temp = ( -4 * par[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_GE_sum = FB_GE_sum + FB_GE_temp;
    }

  //Calculate GM.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_GM_temp = ( -4 * par[i-1+nFB] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_GM_sum = FB_GM_sum + FB_GM_temp;
    }

  //val = FB_sum;
  if(useFB_GM == 0)
    {
      //val = pow(Z,2.) * mottxs * pow(FB_GE_sum,2.)/(1+tau); //Just GE (GM=0). This was for the recoil proton being measure.
      val = mottxs * pow(FB_GE_sum,2.)/(1+tau);
    }
  if(useFB_GM == 1)
    {
      //val = pow(Z,2.) * mottxs * (   (  pow(FB_GE_sum,2.)+tau*pow(FB_GM_sum,2.)  )/(1+tau) + 2*tau*pow(FB_GM_sum,2.)*pow(1/(1+rho),2.)*pow(1/tan(theta*deg2rad),2.)  );  //Try to fit GE and GM. This was for the recoil proton being measure.
      val = mottxs * (   (  pow(FB_GE_sum,2.)+tau*pow(FB_GM_sum,2.)  )/(1+tau) + 2*tau*pow(FB_GM_sum,2.)*pow(tan(theta*deg2rad/2.),2.)  );
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
  for(Int_t i=0;i<datapts;i++) 
    //for(Int_t i=0;i<226;i++) 
    {
      delta  = (sigexp[i]-XS(E0[i],theta[i],par))/uncertainty[i];
      chisq += delta*delta;
      Chi2[i] = delta*delta;
      residual[i] = (sigexp[i] - XS(E0[i],theta[i],par))/sigexp[i]; 
      xsfit[i] = XS(E0[i],theta[i],par);
      //cout<<"xsfit["<<i<<"] = "<<xsfit[i]<<endl;
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
      residual_FB[i] = (sigexp[i] - FB(E0[i],theta[i],par))/sigexp[i]; 
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
  fChFF1->GetHistogram()->GetYaxis()->SetTitle("|F_{ch}(q^{2})|");
  fChFF1->GetHistogram()->GetYaxis()->CenterTitle(true);
  fChFF1->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  fChFF1->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  fChFF1->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
  fChFF1->GetHistogram()->GetXaxis()->SetTitle("q^{2} (fm^{-2})");
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

  c_FF->SaveAs("Test_Print.png");
  //c1->SaveAs(Form("./%i/%s/pictures_png/%i_%s_wire%i_wire%i.png",run,plane,run,plane,wire1,wire2));
*/
}

Double_t test(Double_t X)
{
  cout<<pi<<endl;
  Double_t val = X;
  return val;
}

void Test()
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
 output<<"FCN    Qichtot    Qimtot  R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch    Q0m    Q1m    Q2m    Q3m    Q4m    Q5m    Q6m    Q7m    Q8m    Q9m    Q10m    Q11m"<<endl;



  /*
  Double_t test(float E0, float theta, Double_t *par)
  {
    Double_t val = 0;
    return val;
  }
  */
  cout<<test(muHe3)<<endl;
  cout<<"Hello there!"<<endl;
}

int main()
{
  //Double_t xyz = 5.5;
  Test();
}
