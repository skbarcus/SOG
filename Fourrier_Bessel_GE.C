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

Double_t pi = TMath::Pi();
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t e = 1.60217662E-19;             //Electron charge C.
Double_t alpha = 0.0072973525664;//1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.

Int_t loops = 1;
const Int_t datapts = 6;//246
Int_t userand = 3;                       //0 = use predetermined Ri from Amroun. 1 = use random Ri in generated in a range around Amroun's. 2 = use random Ri generated in increments of 0.1 with larger possible spacing at greater radii. 3 = use predetermined Ri for the purposes of trying to tune the fit by hand.
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 0;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t fft = 0;                           //0 = don't use FFT to try to get a charge radii. 1 = do use FFT to extract a charge radii.
Int_t Amroun_Qi = 0;                     //1 = Override fitted Qi and use Amroun's values.
Int_t showplots = 1;
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
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
Float_t Ef_arr[datapts];
Float_t Q2_arr[datapts];
Float_t Q_arr[datapts];
Float_t Q2eff_arr[datapts];
Float_t Qeff_arr[datapts];
Float_t mottxs_arr[datapts];
Float_t tau_arr[datapts];
Float_t sig_red[datapts];                //Sig_exp/sig_mott * (1+tau)/Z^2
Float_t ratio[datapts];
Float_t ratio_retzlaff[datapts];

Double_t m = 2.;
//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
Double_t R[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t Qich[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
Double_t av[12] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t av_retzlaff[12] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t averr[12] ={};
Double_t Qicherr[12]={}; 
Double_t Qimerr[12]={};
Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};
Double_t Chi2_FB[datapts]={};
Double_t residual_FB[datapts]={};
Double_t FBfit[datapts]={};
Double_t FBfit_retzlaff[datapts]={};

Double_t FB(float E0, float theta, Double_t *par)
{
  Double_t val = 0.;

  Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
  Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
  Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
  Double_t FB_sum = 0.;
  Double_t FB_temp = 0.;
  Double_t R_FB = 5.;  //fm
  Double_t mottxs = 0.;
  Double_t tau = 0;

  //Calculate Mott XS.
  mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/GeV2fm;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  //cout<<"MottXS = "<<mottxs<<endl;
  //Calculate tau.
  tau = Q2eff/(4*pow(MtHe3,2.)*GeV2fm);

  //Calculate Ge.
  for(Int_t i=1; i<(nFB+1); i++)
    {
      FB_temp = ( -4 * par[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      //FB_temp = ( -4 * par[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( i*pi/R_FB * ROOT::Math::cyl_bessel_j(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
      FB_sum = FB_sum + FB_temp;
    }

  //val = FB_sum;
  //val = pow(Z,2.) * mottxs * pow(FB_sum,2.)/(1+tau); //This was for measuring proton recoil.
  val = mottxs * pow(FB_sum,2.)/(1+tau);
  return val;
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
      residual_FB[i] = sigexp[i] - FB(E0[i],theta[i],par); 
      FBfit[i] = FB(E0[i],theta[i],par);
      //cout<<"FBfit["<<i<<"] = "<<FBfit[i]<<endl;
    }
  f = chisq;
}

void Fourrier_Bessel_GE()
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Make a new canvas to plot data.
  //TCanvas* c1=new TCanvas("c1");
  //c1->SetGrid();

  //Read in data from text file.
  //Open file.

  if(usedifmin == 0)
    {
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_640.txt","r");
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data.txt","r");
    }

  if(usedifmin == 1)
    {
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_640.txt","r");
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Dunn_1983.txt","r");
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Dunn_1983_Low_Q2.txt","r");
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

	Q2_arr[nlines-skip] = 4 * E0[nlines-skip] * (E0[nlines-skip]/(1.0+2.0*E0[nlines-skip]*pow(sin(theta[nlines-skip]*deg2rad/2.0),2.0)/MtHe3)) * pow(sin(theta[nlines-skip]*deg2rad/2.0),2.) * GeV2fm;
	Q_arr[nlines-skip] = pow(Q2_arr[nlines-skip],0.5);
	Ef_arr[nlines-skip] = E0[nlines-skip]/(1.0+2.0*E0[nlines-skip]*pow(sin(theta[nlines-skip]*deg2rad/2.0),2.0)/MtHe3);
	Q2eff_arr[nlines-skip] = pow( pow(Q2_arr[nlines-skip],0.5) * (1.0+(1.5*Z*alpha)/(E0[nlines-skip]*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);
	Qeff_arr[nlines-skip] = pow(Q2eff_arr[nlines-skip],2.);
	mottxs_arr[nlines-skip] = (  (pow(Z,2.)*(Ef_arr[nlines-skip]/E0[nlines-skip])) * (pow(alpha,2.0)/(4.0*pow(E0[nlines-skip],2.0)*pow(sin(theta[nlines-skip]*deg2rad/2.0),4.0)))*pow(cos(theta[nlines-skip]*deg2rad/2.0),2.0)  ) * 1.0/GeV2fm;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7. 
	tau_arr[nlines-skip] = Q2_arr[nlines-skip]/(4*pow(MtHe3,2.)*GeV2fm);
	sig_red[nlines-skip] = sigexp[nlines-skip]/mottxs_arr[nlines-skip] * (1+tau_arr[nlines-skip]);

	//cout<<"MottXS = "<<mottxs_arr[nlines-skip]<<endl;
	cout<<"Q2_arr = "<<Q2_arr[nlines-skip]<<"   Q_arr = "<<Q_arr[nlines-skip]<<endl;
	nlines++;
      }
  }

  //Print the data read from the file.
  if(showplots == 1)
    { 
      for(int i=0; i<(nlines-skip); i++)
	{
	  cout<<"E0["<<i<<"] = "<<E0[i]<<"   theta["<<i<<"] = "<<theta[i]<<"   Q^2_arr["<<i<<"] = "<<Q2_arr[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   uncertainty["<<i<<"] = "<<uncertainty[i]<<endl;
	}
      
      cout<<"Number of lines = "<<nlines<<endl;
    }
  fclose(fp);

 //Initiate Minuit for minimization of Fourrier-Bessel fit.
 TMinuit *gMinuit_FB = new TMinuit(nFB);  //initialize TMinuit with a maximum of 24 params
 gMinuit_FB->SetFCN(fcn_FB);
 Double_t arglist_FB[10];
 Int_t ierflg_FB = 0;
 
 arglist_FB[0] = 1.;
 gMinuit_FB->mnexcm("SET ERR", arglist_FB ,1,ierflg_FB);

 //Set step sizes.
 static Double_t stepsize_FB[4] = {0.001 , 0.1 , 0.01 , 0.001};
 gSystem->Load("libMathMore"); //Needed to use cyl_bessel_j() function. I have no earthly idea why this needs to be located here but it does.
 //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
 for(Int_t i=0;i<nFB;i++)
   {
     gMinuit_FB->mnparm(i, Form("av%d",i+1), av[i], stepsize_FB[0], 0.,1.,ierflg_FB);
     //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.001,Qich[i]+0.001,ierflg);
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
 
 
 if(showplots == 1)
   { 
     for(Int_t i=0;i<(nlines-skip);i++)
       {
	 cout<<"Chi2_FB["<<i<<"] = "<<Chi2_FB[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   FB_XS(E0,theta,par)["<<i<<"] = "<<FBfit[i]<<"   XSexp/FB_XS_fit = "<<sigexp[i]/FBfit[i]<<endl;//"   residual["<<i<<"] = "<<residual[i]<<endl;
       }
   }
 
 //Fill av[i] with the fitted coefficients. 
 for(Int_t i=0;i<nFB;i++)
   {
     gMinuit->GetParameter(i,av[i],averr[i]);
     cout<<"av["<<i<<"] = "<<av[i]<<"   averr["<<i<<"] = "<<averr[i]<<endl;
   }
 
 //Plot FB fit of GE from Retzlaff 1984.
 Double_t FB_Q2(Double_t *Q2, Double_t *par)
 {
   Double_t val = 0.;
   Double_t FB_sum = 0.;
   Double_t FB_temp = 0.;
   Double_t R_FB = 5.;  //fm
   Double_t mottxs = 0.;
   Double_t tau = 0;
   /*
   //Calculate Mott XS.
   mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
   
   //Calculate tau.
   tau = Q2[0]/(4*pow(MtHe3,2.));
   */

   //Calculate Ge.
   for(Int_t i=1; i<(nFB+1); i++)
     {
       FB_temp = ( -4 * av[i-1] * sin( pow(Q2[0],0.5) * R_FB ) ) / ( pow(Q2[0],0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2[0] - pow(i*pi/R_FB,2.))  );
       FB_sum = FB_sum + FB_temp;
     }
   
   //val = FB_sum;                 //If want just GE.
   val = FB_sum * FB_sum;        //If want GE^2.
   return val;
 }

 Double_t FB_Q_retzlaff(Double_t *Q, Double_t *par)
 {
   Double_t val = 0.;
   Double_t FB_sum = 0.;
   Double_t FB_temp = 0.;
   Double_t R_FB = 5.;  //fm

   //Calculate Ge.
   for(Int_t i=1; i<(nFB+1); i++)
     {
       FB_temp = ( -4 * av_retzlaff[i-1] * sin( Q[0] * R_FB ) ) / ( Q[0] * i * ROOT::Math::sph_bessel(1,i*pi) * (pow(Q[0],2.) - pow(i*pi/R_FB,2.))  );
       FB_sum = FB_sum + FB_temp;
       //FB_temp = 0;
     }
   
   val = FB_sum;                 //If want just GE.
   //val = FB_sum * FB_sum;        //If want GE^2.
   return val;
 }

 Double_t FB_Q2_retzlaff(Double_t *Q2, Double_t *par)
 {
   Double_t val = 0.;
   Double_t FB_sum = 0.;
   Double_t FB_temp = 0.;
   Double_t R_FB = 5.;  //fm

   //Calculate Ge.
   for(Int_t i=1; i<(nFB+1); i++)
     {
       FB_temp = ( -4 * av_retzlaff[i-1] * sin( pow(Q2[0],0.5) * R_FB ) ) / ( pow(Q2[0],0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2[0] - pow(i*pi/R_FB,2.))  );
       FB_sum = FB_sum + FB_temp;
       //FB_temp = 0;
     }
   
   val = FB_sum;                 //If want just GE.
   //val = FB_sum * FB_sum;        //If want GE^2.
   return val;
 }

 TCanvas* c1=new TCanvas("c1");
 c1->SetGrid();
 //c1->SetLogy();

 graph = new TGraph(nlines-skip,Q2_arr,sig_red);
 //Draw the new TGraph called graph on the canvas. 
 graph->GetXaxis()->SetLimits(0.,8.);
 graph->SetMinimum(0);
 graph->SetMaximum(1.);
 //graph->SetLineWidth(3);
 //graph->SetLineColor(4);
 graph->SetFillColor(0);
 graph->SetMarkerStyle(20);
 graph->SetTitle("Main Title; X Axis Title; Y Axis Title");
 c1->Update();
 graph->Draw("ap");
 c1->Update();
 //graph_expected.SetFillColor(kYellow);
 //graph_expected.DrawClone("E3AL"); // E3 draws the band
 //leg.AddEntry(graph,"Experimental Cross Section Data");
  

 TF1 *FB_func = new TF1("FB_func",FB_Q2,yminFF,ymaxFF+54,1);
 FB_func->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
 FB_func->Draw("same L");

 //Plot FB fit of rho from Retzlaff 1984.
 Double_t rho_r(Double_t *r, Double_t *par)
 {
   Double_t val = 0.;
   Double_t rho_sum = 0.;
   Double_t rho_temp = 0.;
   Double_t R_FB = 5.;  //fm
   for(Int_t i=1; i<(nFB+1); i++)
     {
       
       rho_temp = av[i-1] * ROOT::Math::sph_bessel(0,(i*pi/R_FB)*r[0]);
       rho_sum = rho_sum + rho_temp;
       //cout<<i<<"   rho_temp = "<<rho_temp<<"   rho_sum = "<<rho_sum<<endl;
     }
   val = rho_sum;
   return val;
 }

 TCanvas* c2=new TCanvas("c2");
 c2->SetGrid();

 TF1 *rho_func = new TF1("rho_func",rho_r,0.,5.,1);
 rho_func->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
 rho_func->Draw("L");
 cout<<rho_func->Eval(5.)<<endl;

 for(Int_t i=0; i<datapts; i++)
   {
     //cout<<"MottXS = "<<mottxs_arr[i]<<endl;
     cout<<"sig_red["<<i<<"] = "<<sig_red[i]<<"   GE^2["<<i<<"] = "<<FB_func->Eval(Q2_arr[i])<<endl;
   }
 //Check that Root plots the Bessel functions correctly.
 /*
 TCanvas* c_bessel=new TCanvas("c3");
 c_bessel->SetGrid();
 TF1 *j0 = new TF1("j0","ROOT::Math::cyl_bessel_j(0,x)",0.,5.);
 //TF1 *j0 = new TF1("j0","sin(x)",0.,5.);
 TF1 *j1 = new TF1("j1","ROOT::Math::cyl_bessel_j(1,x)",0.,5.);
 j0->Draw();
 j1->Draw("same");
 */

 //Calculate Retzlaff's XSs from t=his fits.
 Double_t FB_retzlaff(float E0, float theta)
 {
   Double_t val = 0.;
   
   Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
   Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
   Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
   Double_t FB_sum = 0.;
   Double_t FB_temp = 0.;
   Double_t R_FB = 5.;  //fm
   Double_t mottxs = 0.;
   Double_t tau = 0;
   
   
   //Calculate Mott XS.
   mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/GeV2fm;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
   //cout<<"MottXS = "<<mottxs<<endl;
   //Calculate tau.
   tau = Q2eff/(4*pow(MtHe3,2.)*GeV2fm);

   //Calculate Ge.
   for(Int_t i=1; i<(nFB+1); i++)
     {
       FB_temp = ( -4 * av_retzlaff[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( pow(Q2eff,0.5) * i * ROOT::Math::sph_bessel(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
       //FB_temp = ( -4 * par[i-1] * sin( pow(Q2eff,0.5) * R_FB ) ) / ( i*pi/R_FB * ROOT::Math::cyl_bessel_j(1,i*pi) * (Q2eff - pow(i*pi/R_FB,2.))  );
       FB_sum = FB_sum + FB_temp;
     }
   
   //val = FB_sum;
   //val = pow(Z,2.) * mottxs * pow(FB_sum,2.)/(1+tau); //This was for measuring proton recoil.
   val = mottxs * pow(FB_sum,2.)/(1+tau);
   return val;
 }

 for(Int_t i=0;i<datapts;i++)
   {
     FBfit_retzlaff[i] = FB_retzlaff(E0[i],theta[i]);
     ratio_retzlaff[i] = sigexp[i]/FBfit_retzlaff[i];
     cout<<"FBfit_retzlaff["<<i<<"] = "<<FBfit_retzlaff[i]<<"   ratio_retzlaff["<<i<<"] = "<<ratio_retzlaff[i]<<endl;
   }
 /*
 for(Int_t i=0;i<datapts;i++)
   {
     cout<<"ratio_retzlaff["<<i<<"] = "<<ratio_retzlaff[i]<<endl;
   }
 */
 TCanvas* c3=new TCanvas("c3");
 c3->SetGrid();

 graph1 = new TGraph(nlines-skip,Q_arr,ratio_retzlaff);
 //Draw the new TGraph called graph on the canvas. 
 graph1->GetXaxis()->SetLimits(0.,2.);
 graph1->SetMinimum(0.6);
 graph1->SetMaximum(1.4);
 //grap1->SetLineWidth(3);
 //graph1->SetLineColor(4);
 graph1->SetFillColor(0);
 graph1->SetMarkerStyle(20);
 graph1->SetTitle("Ratio of experimental XS to FB XS Fit from Retzlaff; q; Ratio of X to FB XS Fit");
 graph1->Draw("ap");

 TCanvas* c4=new TCanvas("c4");
 c4->SetGrid();
 /*
 TH1 *hratio = new TH1D("hratio", "hratio", 1000, 0., 5.);
 for(Int_t i=0;i<datapts;i++)
   {
     hratio->Fill(sigexp[i]/FBfit[i]);
   }
 hratio->Draw();
 */
 for(Int_t i=0;i<datapts;i++)
   {
     ratio[i] = sigexp[i]/FBfit[i];
   }

 graph2 = new TGraph(nlines-skip,Q_arr,ratio);
 //Draw the new TGraph called graph on the canvas. 
 graph2->GetXaxis()->SetLimits(0.,2.);
 graph2->SetMinimum(0.6);
 graph2->SetMaximum(1.4);
 //graph2->SetLineWidth(3);
 //graph2->SetLineColor(4);
 graph2->SetFillColor(0);
 graph2->SetMarkerStyle(20);
 graph2->SetTitle("Ratio of experimental XS to FB XS Fit; q; Ratio of X to FB XS Fit");
 graph2->Draw("ap");

 TCanvas* c4=new TCanvas("c5");
 c5->SetGrid();
 c5->SetLogy();
 TF1 *FB_func_retzlaff = new TF1("FB_func_retzlaff",FB_Q2_retzlaff,0.,60.,1);
 FB_func_retzlaff->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
 FB_func_retzlaff->Draw("L");

 st->Stop();
 cout<<"*********************************************"<<endl;
 cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
