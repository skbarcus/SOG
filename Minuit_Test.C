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

Double_t pi = 3.141592654;
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t alpha = 1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.

Int_t userand = 1;
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 0;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
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
Float_t theta[100];                     //Angle in degrees.
Float_t qeff[100];                      //q effective in fm^-1.
Float_t sigexp[100];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[100];
Float_t E0[100];

Double_t m = 2.;
//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
Double_t R[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t Qich[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};




  //Make a function for the XS using the SOG parameterization that can be minimized to fit the measured cross sections.
  //This XS function fits only the Qi values.
  Double_t XS(float E0, float theta, Double_t *par)
  {//cout<<"!!!!!!!!!!!!Yo!!!!!!!!!!!!"<<endl;
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
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=6 A=12
                
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
    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 
	//Fit just the Qi values using predetermined R[i] values.
	sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
	
       fitch =  fitch + sumchtemp;
      }
    //}
    //fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
    fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
   
    //if(par[ngaus+0]+par[ngaus+1]+par[ngaus+2]+par[ngaus+3]+par[ngaus+4]+par[ngaus+5]+par[ngaus+6]+par[ngaus+7]+par[ngaus+8]+par[ngaus+9]+par[ngaus+10]+par[ngaus+11] == 1.)
    //{
    //Define SOG for magnetic FF.
    for(Int_t i=0; i<ngaus; i++)
      {
	//Fit just the Qi values using predetermined R[i] values.
	summtemp = (par[ngaus+i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );	

	fitm = fitm + summtemp;
      }
    //}
    fitm = fitm * exp(-0.25*Q2eff*pow(gamma,2.0));   //For some reason had fabs(fitm).

    val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2))*(0.5*Q2eff/q2_3 + pow(tan(theta*deg2rad/2),2))*pow(fitm,2.) ); //magnetic moment for C12 is 0 -> no mag part of XS.
    return val;
  }




 Double_t Test(float E0, float theta, Double_t *par)
  {
    //cout<<"!!!!!!!!!!!!Yo!!!!!!!!!!!!"<<endl;
    Double_t val = 0.;
    return val;
  }

  //Create a Chi^2 function to minimize. 
  void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
  {
    const Int_t nbins = 59;
    //Int_t i;
    //cout<<"Yo1111111111"<<endl;
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    //cout<<"Yo2222222222"<<endl;
    for(Int_t i=0;i<nbins; i++) 
      {
	//cout<<"Yo333333333"<<endl;
	delta  = (sigexp[i]-XS(E0[i],theta[i],par))/uncertainty[i];
	chisq += delta*delta;
	//cout<<"Yo444444444"<<endl;
      }
    //cout<<"Yo555555555"<<endl;
    f = chisq;
    //cout<<"Yo666666666"<<endl;
  }

void Minuit_Test() 
{
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
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data.txt","r");
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
	nlines++;
      }
  }

//Print the data read from the file. 
  for(int i=0; i<59; i++)
    {
      cout<<"E0["<<i<<"] = "<<E0[i]<<"   theta["<<i<<"] = "<<theta[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   uncertainty["<<i<<"] = "<<uncertainty[i]<<endl;
    }

  cout<<"Number of lines = "<<nlines<<endl;
  fclose(fp);

  //Initiate minimization procedure. 
  TMinuit *gMinuit = new TMinuit(24);  //initialize TMinuit with a maximum of 24 params
  gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1.;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

   //Set step sizes.
   static Double_t step[4] = {0.001 , 0.1 , 0.01 , 0.001};

   //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
   for(Int_t i=0;i<ngaus;i++)
     {
       gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], step[0], 0,0,ierflg);
     }
   for(Int_t i=0;i<ngaus;i++)
     {
       gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], step[0], 0,0,ierflg);
     }
   
   // Now ready for minimization step
   arglist[0] = 50000.;
   arglist[1] = 1.;
   cout<<"Sup1"<<endl;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   cout<<"Sup2"<<endl;
   // Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);


}
