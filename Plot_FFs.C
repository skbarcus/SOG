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

Int_t Rep_Fit = 122;    //Fit# from Multifit requires -1. From Plot_FFs or Chi2_Sorted Fit# is just the number.
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
Int_t ngaus_Amroun = 12;                        //Number of Gaussians used to fit data from Amroun.
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
Double_t gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
Double_t theta = 0.;//21.04;
Double_t theta_cor = 0.;                //Theta that corrects for the Q^2eff adjustment. Basically when we plot the XS and FFs the Q2[0] is really Q^2eff if we don't use this theta_cor. This variable is for the slightly smaller theta representing the real scattering angle.
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
Double_t R[12] = {0.3,0.7,0.9,1.4,1.7,2.2,2.9,3.6,4.3,4.9};//Min Chi^2 of 357 (#556 combined). 
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit

//Double_t Qich[12] = {0.0440183,0.116665,0.202577,0.26934,0.0690628,0.179559,0.0854789,0.0318623,0.00963141,9.4369e-16,0.,0.};//9/15/18 51.
//Double_t Qim[12] = {0.0725889,0.0926902,0.209459,0.231037,0.079899,0.188591,0.0836548,0.0440054,0.0235449,2.42131e-10,0.,0.};
//Double_t Qich[12] = {0.0366999,0.207852,0.0914645,0.143132,0.189135,0.193206,0.0976287,0.0374533,0.0115873,7.80417E-10,0.,0.};//48
//Double_t Qim[12] = {0.0589855,0.265678,2.06176E-10,0.248498,0.068945,0.209543,0.0894644,0.0417275,0.0185915,4.6805E-10,0.,0.};
//Double_t Qich[12] = {0.0877489,0.157569,0.251965,0.153987,0.18569,0.114409,0.042258,0.0141493,5.25852e-12,0.,0.,0.};//9/14/18 42
//Double_t Qim[12] = {0.151938,0.0687805,0.318107,0.0627664,0.234607,0.100287,0.0562909,0.00777628,0.018992,0.,0.,0.};
//Double_t Qich[12] = {0.087498,0.157869,0.245378,0.101457,0.223892,0.123962,0.0483213,0.014929,0.00494579,0.,0.,0.};//9/14/18 40
//Double_t Qim[12] = {0.148784,0.0789473,0.304254,0.026821,0.248647,0.120223,0.05251,0.0205607,0.0136972,0.,0.,0.};
//Double_t Qich[12] = {0.0973238,0.238908,6.91669E-14,0.287191,0.219123,0.107724,0.0415777,0.0141707,0.0010956,0.,0.,0.};//38
//Double_t Qim[12] = {0.156104,0.129574,0.314365,1.9168E-13,0.262833,0.127868,0.098301,0.0779688,0.0382781,0.,0.,0.};
//Double_t Qich[12] = {0.0346018,0.186348,0.296016,0.208675,0.141902,0.0958882,0.0407531,0.000755174,0.00273895,0.,0.,0.};//34
//Double_t Qim[12] = {0.0897322,0.129928,0.299251,0.198647,0.131083,0.1176,0.0376243,0.0275592,0.00978282,0.,0.,0.};
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
//Double_t Qich[12] = {0.109056,0.322187,2.53868E-10,0.271415,0.140634,0.110921,0.0422077,0.00991945,9.37812E-10,0.,0.,0.};//23
//Double_t Qim[12] = {0.168311,0.172808,0.166683,0.189821,0.155563,0.169979,0.103117,0.111622,0.0574435,0.,0.,0.};
//Double_t Qich[12] = {0.110029,0.313268,0.188735,0.195501,0.116327,0.0664477,0.0154455,0.000856117,0.,0.,0.,0.};//18
//Double_t Qim[12] = {0.15363,0.241968,0.227636,0.141837,0.150786,0.0857793,0.0873129,0.0398393,0.,0.,0.,0.};
//Double_t Qich[12] = {0.0649212,0.308528,0.299338,0.187702,0.102192,0.0379261,0.00712276,2.74347e-12,0.,0.,0.,0.};//11
//Double_t Qim[12] = {0.0997782,0.274768,0.285966,0.191182,0.109338,0.0645759,0.0435615,0.0126547,0.,0.,0.,0.};
//Double_t Qich[12] = {0.0784469,0.247165,0.406019,0.120177,0.137968,4.57535E-11,0.0200847,2.63439E-9,0.,0.,0.,0.};//7
//Double_t Qim[12] = {0.0770148,0.298502,0.282963,0.175066,0.0769078,0.0381075,0.0899692,0.0675,0.,0.,0.,0.};
//Double_t Qich[12] = {0.0896211,0.234994,0.367953,0.193344,0.0460444,0.0763543,1.05355e-12,1.20741E-10,0.,0.,0.,0.};//5
//Double_t Qim[12] = {0.0738571,0.287688,0.219872,0.125082,0.000139601,7.65055e-13,0.0676226,0.0808108,0.,0.,0.,0.};
//Double_t Qich[12] = {0.0433495,0.189163,0.153897,0.236151,0.140139,0.155516,0.0690597,0.0159103,0.00495952,2.61974E-12,0.,0.};//4.
//Double_t Qim[12] = {0.0711952,0.177756,0.146683,0.219108,0.142134,0.180037,0.0992368,0.0517247,0.0363374,0.0161585,0.,0.};
//Double_t Qich[12] = {0.0419664,0.220685,0.160807,0.293756,0.0525773,0.147678,0.0700733,0.0169877,0.00368665,6.48132E-12,0.,0.};//9/12/18 3.
//Double_t Qim[12] = {0.0701842,0.207855,0.149425,0.271277,0.064722,0.153486,0.083882,0.0352331,0.0250242,1.03941E-8,0.,0.};
//Double_t Qich[12] = {0.0411639,0.236983,0.204183,0.276745,0.0977119,0.0724304,0.0541865,0.0207177,1.26197e-07,0.00411213,0.,0.};//2
//Double_t Qim[12] = {0.07112,0.213554,0.205701,0.238441,0.162279,0.017472,0.0869983,0.0239365,0.0295892,5.58285e-10,0.,0.};
//Double_t Qich[12] = {0.0411535,0.237047,0.203676,0.277975,0.092553,0.0821181,0.0562533,0.0157997,2.61438e-10,0.00176867,0.,0.};//9/12/18 pretty good 1.
//Double_t Qim[12] = {0.0692335,0.221858,0.189445,0.254013,0.127296,0.0508983,0.0726007,0.0160226,0.0224739,4.38656E-10,0.,0.};
//Double_t Qim[12] = {0.059577,0.0966408,0.24088,0.227914,0.123278,0.160815,0.0806045,0.00262875,0.0494118,0.0140003,0.,0.};
//Double_t Qich[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
//Double_t Qim[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};

//Double_t Qich[12] = {0.0511376,0.212372,0.255953,0.235852,0.104617,0.11099,0.0458289,0.0164548,0.00625,0.000121054};   //Avg. Ri 357 fits x^2<505.
//Double_t Qim[12] = {0.106411,0.21999,0.219267,0.210835,0.118368,0.120021,0.0755027,0.0153253,0.0224693,0.0127534};     //Avg. Ri 357 fits x^2<505.
Double_t Qich[12] = {0.101632,0.170317,0.165379,0.250085,0.0841009,0.140841,0.0673786,0.0214391,0.00704491,1.00142e-13};   //Min Chi^2 of 357 (#556 combined).
Double_t Qim[12] = {0.156834,0.0554965,0.250939,0.18088,0.120042,0.136881,0.0686315,0.0252427,0.0147734,2.08722e-13};     //Min Chi^2 of 357 (#556 combined).
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

const Int_t size = 1000;
Float_t Fit[size], Chi2[size],rChi2[size],BIC[size],AIC[size],Qichtot[size],Qimtot[size],R0[size],R1[size],R2[size],R3[size],R4[size],R5[size],R6[size],R7[size],R8[size],R9[size],R10[size],R11[size],Q0ch[size],Q1ch[size],Q2ch[size],Q3ch[size],Q4ch[size],Q5ch[size],Q6ch[size],Q7ch[size],Q8ch[size],Q9ch[size],Q10ch[size],Q11ch[size],Q0m[size],Q1m[size],Q2m[size],Q3m[size],Q4m[size],Q5m[size],Q6m[size],Q7m[size],Q8m[size],Q9m[size],Q10m[size],Q11m[size];
Float_t Fittemp,Chi2temp,rChi2temp,BICtemp,AICtemp,Qichtottemp,Qimtottemp,R0temp,R1temp,R2temp,R3temp,R4temp,R5temp,R6temp,R7temp,R8temp,R9temp,R10temp,R11temp,Q0chtemp,Q1chtemp,Q2chtemp,Q3chtemp,Q4chtemp,Q5chtemp,Q6chtemp,Q7chtemp,Q8chtemp,Q9chtemp,Q10chtemp,Q11chtemp,Q0mtemp,Q1mtemp,Q2mtemp,Q3mtemp,Q4mtemp,Q5mtemp,Q6mtemp,Q7mtemp,Q8mtemp,Q9mtemp,Q10mtemp,Q11mtemp;
Float_t Rmulti[size][12];
Float_t Qichmulti[size][12];
Float_t Qimmulti[size][12];

void Plot_FFs() 
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

  //Set SOG parameters to the representative fit's parameters.
  for(Int_t i=0;i<ngaus;i++)
    {
      R[i] = Rmulti[Rep_Fit][i];
      Qich[i] = Qichmulti[Rep_Fit][i];
      Qim[i] = Qimmulti[Rep_Fit][i];
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

  TF1 *fChFF_Amroun = new TF1("fChFF_Amroun",ChFF_Q2_Amroun,yminFF,ymaxFF+54,1);
  //cout<<fChFF_Amroun->Eval(35)<<endl;
  fChFF_Amroun->SetNpx(npdraw);
  fChFF_Amroun->SetLineColor(4);
  fChFF_Amroun->Draw("L same");
  auto ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
  ChFF_leg->AddEntry("fChFF","New ^{3}He |F_{ch}(q^{2})| Fit","l");
  ChFF_leg->AddEntry("fChFF_Amroun","^{3}He |F_{ch}(q^{2})| Fit from Amroun et al. [4]","l");
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

  TF1 *fMFF = new TF1("fMFF",MFF_Q2,0.,70.,1);
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

  TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,ymaxFF+54,1);
  //cout<<fMFF_Amroun->Eval(30)<<endl;
  fMFF_Amroun->SetNpx(npdraw);
  fMFF_Amroun->SetLineColor(4);
  fMFF_Amroun->Draw("L same");
  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  MFF_leg->AddEntry("fMFF","New ^{3}He |F_{m}(q^{2})| Fit","l");
  MFF_leg->AddEntry("fMFF_Amroun","^{3}He |F_{m}(q^{2})| Fit from Amroun et al. [4]","l");
  MFF_leg->Draw();

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
    mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta/2.0),4.0)))*pow(cos(theta/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

    //Calculate XS from FFs.
    val = mottxs * (1./eta) * ( (Q2[0]/q2_3)*pow(ChFF_Q2(Q2,par),2.) + (pow(muHe3,2.0)*Q2[0]/(2*pow(MtHe3,2)*GeV2fm))*(0.5*Q2[0]/q2_3 + pow(tan(theta/2),2))*pow(MFF_Q2(Q2,par),2.) ); 
    
    return val;
    
    //return MFF_Q2_Amroun(Q2,par) + ChFF_Q2_Amroun(Q2,par);
  }

  TF1 *fxs = new TF1("fxs",xs,yminFF,ymaxFF+54,1);
  //fMFF_Amroun->Draw("L");
  fxs->SetNpx(npdraw);
  fxs->Draw("L");
  fxs->SetTitle("^{3}He Cross Section at 3.356 Gev and 21.04#circ");
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
    mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta/2.0),4.0)))*pow(cos(theta/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

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

  auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
  MFF_leg->AddEntry("fxs","New ^{3}He Cross Section","l");
  MFF_leg->AddEntry("fxs_Amroun","^{3}He Cross Section from Amroun et al.","l");
  MFF_leg->Draw();

  //Plot my data point.
  m1 = new TMarker(34.0981, 1.33459E-10, 20);
  m1->SetMarkerColor(kOrange+7);
  m1->SetMarkerSize(1);
  m1->Draw();

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
  cout<<"Amroun Mott = "<<(  (pow(Z,2.)*((E0/(1.0+2.0*E0*pow(sin((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),2.0)/MtHe3))/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),4.0)))*pow(cos((2*TMath::ASin(  pow( (1/((4*pow(E0,2.)*GeV2fm/test_point)-(2.*E0/MtHe3))) , 0.5 )  ))/2.0),2.0)  ) * 1.0/25.7<<endl;
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
