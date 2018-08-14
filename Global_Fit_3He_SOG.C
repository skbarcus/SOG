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

Int_t loops = 1;
const Int_t datapts = 248;
Int_t userand = 2;
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 0;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t fft = 0;                           //0 = don't use FFT to try to get a charge radii. 1 = do use FFT to extract a charge radii.
Int_t showplots = 1;
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
Float_t theta[1000];                     //Angle in degrees.
Float_t qeff[1000];                      //q effective in fm^-1.
Float_t sigexp[1000];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[1000];
Float_t E0[1000];

Double_t m = 2.;
//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
Double_t R[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t Qich[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};
Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};
Double_t Qicherr[12]={}; 
Double_t Qimerr[12]={};
Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};

  //Make a function for the XS using the SOG parameterization that can be minimized to fit the measured cross sections.
  //This XS function fits only the Qi values.
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
	sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
	
       fitch =  fitch + sumchtemp;
       //cout<<"fitch["<<i<<"] = "<<fitch<<endl;
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
	//cout<<"fitm["<<i<<"] = "<<fitm<<endl;
      }
    //}
    fitm = fitm * exp(-0.25*Q2eff*pow(gamma,2.0));   //For some reason had fabs(fitm).
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

  //Create a Chi^2 function to minimize. 
  void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
  {
    const Int_t nbins = datapts;//177
    //Int_t i;
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    Double_t res;
    for(Int_t i=0;i<nbins; i++) 
      {
	delta  = (sigexp[i]-XS(E0[i],theta[i],par))/uncertainty[i];
	chisq += delta*delta;
	Chi2[i] = delta*delta;
	residual[i] = sigexp[i] - XS(E0[i],theta[i],par); 
	xsfit[i] = XS(E0[i],theta[i],par);
	//cout<<"xsfit["<<i<<"] = "<<xsfit[i]<<endl;
      }
    f = chisq;
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

  if(usedifmin == 0)
    {
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_640.txt","r");
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data.txt","r");
    }

  if(usedifmin == 1)
    {
      //FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_640.txt","r");
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Amroun_3He_Data.txt","r");
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
	nlines++;
      }
  }

  //Print the data read from the file.
  if(showplots == 1)
    { 
      for(int i=0; i<(nlines-skip); i++)
	{
	  cout<<"E0["<<i<<"] = "<<E0[i]<<"   theta["<<i<<"] = "<<theta[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   uncertainty["<<i<<"] = "<<uncertainty[i]<<endl;
	}
      
      cout<<"Number of lines = "<<nlines<<endl;
    }
  fclose(fp);

  //Create an output file to store fit results.
 std::ofstream output ("Ri_Chi2.txt", std::ofstream::out);
 output<<"FCN    Qichtot    Qimtot  R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch    Q0m    Q1m    Q2m    Q3m    Q4m    Q5m    Q6m    Q7m    Q8m    Q9m    Q10m    Q11m"<<endl;

 //Begin loop over fit with different Ri values each time.
 for(Int_t q=0;q<loops;q++)
   {
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

  arglist[0] = 1.;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set step sizes.
  static Double_t stepsize[4] = {0.001 , 0.1 , 0.01 , 0.001};
  
  //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
  for(Int_t i=0;i<ngaus;i++)
    {
      gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], 0.,1.,ierflg);
      //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.001,Qich[i]+0.001,ierflg);
      //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.000000000001,Qich[i]+0.000000000001,ierflg);
      //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.05,Qich[i]+0.05,ierflg);
    }
  for(Int_t i=0;i<ngaus;i++)
    {
      gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], 0.,1.,ierflg);
      //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], Qim[i]-0.001,Qim[i]+0.001,ierflg);
      //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], Qim[i]-0.000000000001,Qim[i]+0.000000000001,ierflg);
      //gMinuit->mnparm(ngaus+i, Form("Qim%d",i+1), Qim[i], stepsize[0], Qim[i]-0.05,Qim[i]+0.05,ierflg);
    }
  
  // Now ready for minimization step
  arglist[0] = 500.;//50000.
  arglist[1] = 1.;
  //cout<<"Sup1"<<endl;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  //cout<<"Sup2"<<endl;
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  
  //Print Chi^2 and residual for each data point.
  if(showplots == 1)
    { 
      for(Int_t i=0;i<(nlines-skip);i++)
	{
	  cout<<"Chi2["<<i<<"] = "<<Chi2[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   XSfit(E0,theta,par)["<<i<<"] = "<<xsfit[i]<<"   XSexp/XSfit = "<<sigexp[i]/xsfit[i]<<endl;//"   residual["<<i<<"] = "<<residual[i]<<endl;
	}
    }

  Double_t maxchi2 = 0.;
  Double_t total_chi2 = 0.;
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
      
      TH2D *hchi = new TH2D("hchi","\Chi^{2} vs. Scattering Angle" , 161, 0., 160., maxchi2+21,0., maxchi2+20);
      for(Int_t i=0;i<(nlines-skip);i++)
	{
	  //hchi->SetBinContent(E0[i],Chi2[i],1.);
      hchi->Fill(theta[i],Chi2[i]);
	}
      hchi->SetMarkerStyle(20);
      hchi->SetMarkerSize(1);
      hchi->Draw();
      
      //Create a plot of XSexp/XSfit.
      TCanvas* cxsfit=new TCanvas("cxsfit");
      cxsfit->SetGrid();
      
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
      hxsfit->Fill(theta[i],sigexp[i]/xsfit[i]);
	}
      hxsfit->SetMarkerStyle(20);
      hxsfit->SetMarkerSize(1);
      hxsfit->Draw();
    }//End showplots.  
  
  Double_t xs2(Double_t *angle, Double_t *par)
  {
    Double_t val = 0.;
    Double_t mottxs = 0.;
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;
    Double_t fitm = 0.;
    Double_t summtemp = 0.;
    
    Ef = E0/(1.0+2.0*E0*pow(sin(angle[0]*deg2rad/2.0),2.0)/MtHe3);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle[0]*deg2rad/2.0),2.0) * GeV2fm;
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=6 A=12
    
    Double_t W = E0 - Ef;
    //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
    Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    Double_t eta = 1.0 + Q2eff/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.
    
    //Calculate Mott XS.
    mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle[0]*deg2rad/2.0),4.0)))*pow(cos(angle[0]*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    
    //Double_t sumchtemp1 = 0.;
    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 
	//cout<<"R["<<i<<"] = "<<R[i]<<endl;
	/*
	  if(fitvars == 0)
	  {
	    //I can't figure out if shit in here influences the fit even if the if statement is never entered. It seems to influence the fit at least in that it can break it.
	    cout<<"Sup?"<<endl;
	//Fit just the Qi values using predetermined R[i] values.
	    //sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
	    //sumchtemp = par[i]*R[i]*gamma*Q2eff+pow(par[i],3.)+10000000.+exp(par[i])/(1+par[i]);
	    //sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) );
	    //sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0)))*1000000.+Q2eff+cos(pow(Q2eff,0.5)*R[i]);
	    //sumchtemp = cos(pow(Q2eff,0.5)*R[i]);
	  }*/
	//cout<<"*********"<<endl;
	//cout<<"sumchtemp = "<<sumchtemp<<endl;
	//Fit R[i] values and Qi.
	//sumchtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );

	//Fit R[i] values, Qi coefficients, and gamma.
	sumchtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(par[3*ngaus],2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(par[3*ngaus],2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );
	//cout<<"sumchtemp = "<<sumchtemp<<endl;
	//cout<<"*********"<<endl;
	//Use 'known' C12 coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
	
       fitch =  fitch + sumchtemp;
      }
    
    //fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
    fitch =  fitch * exp(-0.25*Q2eff*pow(par[3*ngaus],2.0));
   
    
    //Define SOG for magnetic FF.
    for(Int_t i=0; i<ngaus; i++)
      {
	//Fit just the Qi values using predetermined R[i] values.
	//summtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
       
	//Fit R[i] values and Qi.
	//summtemp = (par[2*ngaus+i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );

	//Fit R[i] values, Qi coefficients, and gamma.
	summtemp = (par[2*ngaus+i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(par[3*ngaus],2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(par[3*ngaus],2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );	
	
	//Use 'known' He3 coefficients and R[i] values. 
	//summtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );

	fitm = fitm + summtemp;
      }
    
    fitm = fitm * exp(-0.25*Q2eff*pow(par[3*ngaus],2.0));   //For some reason had fabs(fitm).


    val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2))*(0.5*Q2eff/q2_3 + pow(tan(angle[0]*deg2rad/2),2))*pow(fitm,2.) ); //magnetic moment for C12 is 0 -> no mag part of XS.

    return val;
  }
 




 //Calculate XS fitting only the Qi.
 Double_t xs0(Double_t *angle, Double_t *par)
  {
    Double_t val = 0.;
    Double_t mottxs = 0.;
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;
    Double_t fitm = 0.;
    Double_t summtemp = 0.;

    Ef = E0/(1.0+2.0*E0*pow(sin(angle[0]*deg2rad/2.0),2.0)/MtHe3);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle[0]*deg2rad/2.0),2.0) * GeV2fm;
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
    mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle[0]*deg2rad/2.0),4.0)))*pow(cos(angle[0]*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    
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

    val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2))*(0.5*Q2eff/q2_3 + pow(tan(angle[0]*deg2rad/2),2))*pow(fitm,2.) ); //magnetic moment for C12 is 0 -> no mag part of XS.

    return val;
  }

//Calculate XS fitting the Qi and Ri, but not gamma.
 Double_t xs1(Double_t *angle, Double_t *par)
  {
    Double_t val = 0.;
    Double_t mottxs = 0.;
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;
    Double_t fitm = 0.;
    Double_t summtemp = 0.;

    Ef = E0/(1.0+2.0*E0*pow(sin(angle[0]*deg2rad/2.0),2.0)/MtHe3);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle[0]*deg2rad/2.0),2.0) * GeV2fm;
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=6 A=12
                
    Double_t W = E0 - Ef;
    //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
    Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    Double_t eta = 1.0 + Q2eff/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.

    //Calculate Mott XS.
    mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle[0]*deg2rad/2.0),4.0)))*pow(cos(angle[0]*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    
    //Double_t sumchtemp1 = 0.;
    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 
	//Fit R[i] values and Qi.
	sumchtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );
	
       fitch =  fitch + sumchtemp;
      }
    
    fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
    
    //Define SOG for magnetic FF.
    for(Int_t i=0; i<ngaus; i++)
      {
	//Fit R[i] values and Qi.
	summtemp = (par[2*ngaus+i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );

	fitm = fitm + summtemp;
      }
    
    fitm = fitm * exp(-0.25*Q2eff*pow(gamma,2.0));   //For some reason had fabs(fitm).

    val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2))*(0.5*Q2eff/q2_3 + pow(tan(angle[0]*deg2rad/2),2))*pow(fitm,2.) ); //magnetic moment for C12 is 0 -> no mag part of XS.

    return val;
  }

 TF1 *fxs0 = new TF1("fxs0", xs0, ymin, ymax, ngaus*2);
 TF1 *fxs1 = new TF1("fxs1", xs1, ymin, ymax, ngaus*3);
 TF1 *fxs2 = new TF1("fxs2", xs2, ymin, ymax, ngaus*3+1);
 /*
 if(userand == 1)
   {
     //Generate random R[i] values. 
     Double_t d = 0.49;
     Double_t step = 0.5;
     gRandom->SetSeed(0);                    //Sets new random seed.
     TF1 *rand = new TF1("rand","x",0.,.01);
     R[0] = rand->GetRandom();
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
 */
 /*
 for(Int_t i=0;i<ngaus;i++)
   {
     cout<<"R["<<i<<"] = "<<R[i]<<endl;
   }
 */

 //Plot the charge FF using the parameters determined by the Minuit SOG fit above. 
 //Make a new canvas to plot data.
 if(showplots == 1)
   {
     TCanvas* c2=new TCanvas("c2");
     c2->SetGrid();
     c2->SetLogy();
   }
 //Set Qi and R[i] = to the corresponding fit parameters.

 Double_t  Qichtot = 0.;
 Double_t  Qimtot = 0.;
 
 if(fitvars == 0)
   {
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
	 
	 cout<<"Gamma = "<<gamma<<endl;
       }
   }

 if(fitvars == 1)
   {
     for(Int_t i=0;i<ngaus;i++)
       {
	 Qich[i] = fxs1->GetParameter(i);
	 cout<<"Qich["<<i<<"] = "<<Qich[i]<<endl;
       }
     
     for(Int_t i=0;i<ngaus;i++)
       {
	 Qim[i] = fxs1->GetParameter(2*ngaus+i);
	 cout<<"Qim["<<i<<"] = "<<Qim[i]<<endl;
       }
     
     for(Int_t i=0;i<ngaus;i++)
       {
	 R[i] = fxs1->GetParameter(ngaus+i);
	 cout<<"R["<<i<<"] = "<<R[i]<<endl;
       }
     
     cout<<"Gamma = "<<gamma<<endl;
   }

 if(fitvars == 2)
   {
     for(Int_t i=0;i<ngaus;i++)
       {
	 Qich[i] = fxs2->GetParameter(i);
	 cout<<"Qich["<<i<<"] = "<<Qich[i]<<endl;
       }
     
     for(Int_t i=0;i<ngaus;i++)
       {
	 Qim[i] = fxs2->GetParameter(2*ngaus+i);
	 cout<<"Qim["<<i<<"] = "<<Qim[i]<<endl;
       }
     
     for(Int_t i=0;i<ngaus;i++)
       {
	 R[i] = fxs2->GetParameter(ngaus+i);
	 cout<<"R["<<i<<"] = "<<R[i]<<endl;
       }
     
     cout<<"Gamma = "<<fxs2->GetParameter(3*ngaus)<<endl;
     //If fit gamma too need to redefine gamma to the fit parameter.
     gamma = fxs2->GetParameter(3*ngaus);
   }

 //Fill text file with fcn (chi2), Qichtot, Qimtot, Ri, Qich, Qim.
 output<<amin<<" "<<Qichtot<<" "<<Qimtot<<" "<<R[0]<<" "<<R[1]<<" "<<R[2]<<" "<<R[3]<<" "<<R[4]<<" "<<R[5]<<" "<<R[6]<<" "<<R[7]<<" "<<R[8]<<" "<<R[9]<<" "<<R[10]<<" "<<R[11]<<" "<<Qich[0]<<" "<<Qich[1]<<" "<<Qich[2]<<" "<<Qich[3]<<" "<<Qich[4]<<" "<<Qich[5]<<" "<<Qich[6]<<" "<<Qich[7]<<" "<<Qich[8]<<" "<<Qich[9]<<" "<<Qich[10]<<" "<<Qich[11]<<" "<<Qim[0]<<" "<<Qim[1]<<" "<<Qim[2]<<" "<<Qim[3]<<" "<<Qim[4]<<" "<<Qim[5]<<" "<<Qim[6]<<" "<<Qim[7]<<" "<<Qim[8]<<" "<<Qim[9]<<" "<<Qim[10]<<" "<<Qim[11]<<endl;


}//End loop over minimization for different Ri values.
 //Close output file.
 output.close();

 if(showplots == 1)
   {
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
	sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
	
	fitch = fitch + sumchtemp;
      }
    //Convert to fm. Not sure I need to do this.
    //fitch = fitch * exp(-0.25*pow(Q[0]*pow(GeV2fm,0.5),2.)*pow(gamma,2.0));
    fitch = fitch * exp(-0.25*pow(Q[0],2.)*pow(gamma,2.0));
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
	sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
	fitch = fitch + sumchtemp;
      }
    //Convert to fm. Not sure I need to do this.
    //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
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
    for(Int_t i=0; i<ngaus; i++)
      { 	
	//Use SOG fit for C12 Qi coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Convert to fm. Not sure I need to do this.
	//sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0]*GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0]*GeV2fm,0.5)*R[i])/(pow(Q2[0]*GeV2fm,0.5)*R[i])) );
	sumchtemp = (Qich_Amroun[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
	fitch = fitch + sumchtemp;
      }
    //Convert to fm. Not sure I need to do this.
    //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
    fitch = fitch * exp(-0.25*Q2[0]*pow(gamma,2.0));
    fitch = fabs(fitch);
    return fitch;
 }

 TF1 *fChFF = new TF1("fChFF",ChFF_Q2,yminFF,ymaxFF+54,1);
 //TF1 *fChFF = new TF1("fChFF",ChFF,yminFF,ymaxFF,1);
 cout<<fChFF->Eval(0.000001)<<"!!!!!!!!"<<endl;
 fChFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
 fChFF->Draw("L");
 c2->SetTitle("Charge Form Factor");
 //fChFF->SetTitle("C12 Charge Form Factor","#Q^2 (#fm^-2)","#F_{Ch}(q)");
 fChFF->SetTitle("^{3}He Charge Form Factor");
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
}//End showplots.

 Double_t fitg(Double_t *Q, Double_t *par)
 {
   Double_t val = 0.;

   val = (par[0]/(1.0+2.0*pow(par[1],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*par[1]) + (2.0*pow(par[1],2.0)/pow(gamma,2.0)) * (sin(Q[0]*par[1])/(Q[0]*par[1])) );

   val = val * exp(-0.25*pow(Q[0],2.)*pow(gamma,2.0));

    return val;
 }

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
 

 /*
 //Use ingot's R[i] and Qi values to check my FF plot.
 Double_t sick(Double_t *Q, Double_t *par)
 {
 Double_t fitch = 0.;
	Double_t sumchtemp = 0.;
	gamma = .8;
	Double_t Rsick[12] = {0.00000001,0.4,1.,1.3,1.7,2.3,2.7,3.5,4.3,5.4,6.7};//,7.9,9.1};
	Double_t Qisick[12] = {0.01669, 0.050325, 0.128621, 0.180515, 0.219097, 0.0278416, 0.058779, 0.057817, 0.007739, 0.002001, 0.000007};

    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 
	sumchtemp = (Qisick[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
	
	fitch = fitch + sumchtemp;
      }
    
    fitch = fitch * exp(-0.25*pow(Q[0],2.)*pow(gamma,2.0));
    fitch = fabs(fitch);
    return fitch;
 }*/

 /*
 Double_t test(Double_t *z, Double_t *par)
 {
   Double_t val = 0.;
   val = z[0]*Qi[1]*R[1];
   return val;
 }

 TF1 *ftest = new TF1("ftest",test, 0., 30.,1);
 cout<<ftest->Eval(10.)<<"!!!!!!!!"<<endl;
 ftest->Draw("same");
 */




 if(fft == 1)
   {

  //Fill a new histogram with data using the fit function. This will then be inverse Fourier transformed to obtain the charge distribution.
  TH1 *hChFF = new TH1D("hChFF", "hChFF", n+1, yminFF, ymaxFF);
  //Fill the histogram with function values
  for (Int_t i=0; i<=n; i++)
    {
      Double_t x = 0.;
      x = yminFF + (Double_t(i)/(n))*(range+0.);//range;  //Only fill bins up to max of fitted range.
      //cout<<x<<endl;
      if(x==0.)
	{
	  x = 0.0000001; //Inf/NaN at zero.
	}
      
      /*
      if(x<=ymaxFF && x<truncate)
	{
	  hChFF->SetBinContent(i+1, fChFF->Eval(x));
	}
      
      if(x>ymaxFF)
	{
	  //hChFF->SetBinContent(i+1, (0.8*cos(12.*x-3.1)+0.5)*exp(-x*.1));
	  //hChFF->SetBinContent(i+1, 0.);
	}
      */

      hChFF->SetBinContent(i+1, fChFF->Eval(x));
      
      hChFF->GetEntries();
    }
  //hfit->SetFillColor(17);
  //hChFF->Draw("same");

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry("fChFF","C12 Charge Form Factor","l");
  legend->AddEntry(hChFF,"C12 Charge Form Factor Histogram","f");
  legend->AddEntry("g0","Gaussian 0","l");
  legend->AddEntry("g1","Gaussian 1","l");
  legend->AddEntry("g2","Gaussian 2","l");
  legend->AddEntry("g3","Gaussian 3","l");
  legend->AddEntry("g4","Gaussian 4","l");
  legend->AddEntry("g5","Gaussian 5","l");
  legend->AddEntry("g6","Gaussian 6","l");
  legend->AddEntry("g7","Gaussian 7","l");
  legend->AddEntry("g8","Gaussian 8","l");
  legend->AddEntry("g9","Gaussian 9","l");
  legend->AddEntry("g10","Gaussian 10","l");
  legend->AddEntry("g11","Gaussian 11","l");
  //legend->Draw();
   
  //Inverse Fourier transform the C12 charge FF to get the charge distribution of C12.
  //Create arrays for real and complex inputs for inverse FFT. 
  Double_t *re_full = new Double_t[n+1];
  Double_t *im_full = new Double_t[n+1];
  
  //Fill the real and complex arrays. The complex array is all zeros since we have real data. The real data is from the histogram binned to the chare form factor dervived from the SOG fit. 
  for(Int_t i=0;i<n+1;i++)
    {
      re_full[i] = hChFF->GetBinContent(i+1);
      im_full[i] = 0;
      //cout<<"re_full["<<i<<"] = "<<re_full[i]<<endl;
    }

  //Make a new canvas to plot data.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  //c3->SetLogy();

  TVirtualFFT *iFFT = TVirtualFFT::FFT(1, &ndim, "C2R M K");   //Re and Mag look the same for C2R which makes sense. C2CBACKWARD mag definitely different from C2R and looks wrong. Same for C2CBackward re although first min x position looks ok. Stick to C2R mag it appears correct. 
  iFFT->SetPointsComplex(re_full,im_full);
  iFFT->Transform();
  TH1 *hcharge = 0;
  //TH1 *hcharge = new TH1D("hcharge", "hcharge", nback*10.+1, ymin, ymax);
  //Let's look at the output
  hcharge = TH1::TransformHisto(iFFT,hcharge,"mag");     //Not totally sure if we want re or mag.
  hcharge->SetTitle("C12 Charge Distribution");
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
  TH1 *Ch_dist = new TH1D("Ch_dist", "Ch_dist", n+1, -(n+1)/(ymaxFF*2.), (n+1)/(ymaxFF*2.));   

  Int_t inflection = (n)/2.;
  //cout<<"inflection = "<<inflectionback<<endl;
  
  //Negative Fourier frequencies.
  for(Int_t i=0;i<n+2;i++)
    {
      if(i>inflection+1)
	{
	  Ch_dist->SetBinContent(i-1-inflection,hcharge->GetBinContent(i)/(1./(range/(n+1.))));//range/(n+1.))));
	  //cout<<"i - inflectionback - 1 = "<<i<<" - "<<inflectionback<<"-1 = "<<i-inflectionback-1<<endl;
	}
    }
  
  //Positive Fourier frequencies.
  for(Int_t i=0;i<n;i++)
    {
      if(i<=inflection)
	{
	  Ch_dist->SetBinContent(i+inflection+1,hcharge->GetBinContent(i+1)/(1./(range/(n+1.))));
	  //cout<<i+inflection+1<<endl;
	}
    }

  Ch_dist->SetTitle("He3 Charge Distribution");
  Ch_dist->GetXaxis()->SetTitle("r");
  Ch_dist->GetYaxis()->SetTitle("#rho(r)");
  Ch_dist->Draw("");

  /*
  //Now (forward) Fourier transform the H3 charge distribution to recover the H3 charge FF. 

  //Make a new canvas to plot data.
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();
  
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
  c4->SetLogy();
  hH3chFFback->SetTitle("H3 Charge Form Factor Transformed back (FF->iFFT(FF)->FFT(iFFT(FF))=FF)");
  hH3chFFback->GetXaxis()->SetTitle("Q^{2} [fm^{-2}]");
  hH3chFFback->GetYaxis()->SetTitle("Fch(Q)");
  hH3chFFback->Draw("");
  */
   }//End fft.

 if(showplots == 1)
   {
     //Plot the magnetic FF using the parameters determined by the SOG fit above. 
     //Make a new canvas to plot data.
     TCanvas* c4=new TCanvas("c4");
     c4->SetGrid();
     c4->SetLogy();
     
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
	   
	   summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
	   
	   fitm = fitm + summtemp;
	 }
       
       fitm = fitm * exp(-0.25*pow(Q[0],2.)*pow(gamma,2.0));
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
	   
	   summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
	   fitm = fitm + summtemp;
	 }
       
       fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
       fitm = fabs(fitm);
       return fitm;
     }
     
     //Reset Ri to equal Amroun's values.
     for(Int_t i=0;i<12;i++)
       {
	 R[i] = R_Amroun[i];
       }
     
     //Define Amroun's FF. Not sure why I can't just change Qi in MFF_Q2.
     Double_t MFF_Q2_Amroun(Double_t *Q2, Double_t *par)
     {
       Double_t fitm = 0.;
       Double_t summtemp = 0.;
       
       //Define SOG for magnetic FF.
       for(Int_t i=0; i<ngaus; i++)
	 { 	
	   //Use SOG fit for C12 Qi coefficients and R[i] values. 
	   //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
	   summtemp = (Qim_Amroun[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	   
	   fitm = fitm + summtemp;
	 }
       
       fitm = fitm * exp(-0.25*Q2[0]*pow(gamma,2.0));
       fitm = fabs(fitm);
       return fitm;
     }
     
     TF1 *fMFF = new TF1("fMFF",MFF_Q2,yminFF,ymaxFF+54,1);
     cout<<fMFF->Eval(0.000001)<<"!!!!!!!!"<<endl;
     fMFF->SetNpx(npdraw);   //Sets number of points to use when drawing the function. 
     fMFF->Draw("L");
     c4->SetTitle("He3 Magnetic Form Factor");
     //fChFF->SetTitle("C12 Charge Form Factor","#Q^2 (#fm^-2)","#F_{Ch}(q)");
     //fMFF->GetHistogram()->SetTitle("^{3}He Magnetic Form Factor");
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
     
     //Now draw both FFs on the same plot for the GRC poster. 
     TCanvas* c5=new TCanvas("c5");
     c5->SetGrid();
     c5->Divide(1,2);
     c5->cd(1)->SetLogy();
     c5->cd(1);
     fChFF->Draw("L");
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
     fChFF_Amroun->Draw("L same");
     auto ChFF_leg = new TLegend(0.49,0.64,0.9,0.9); //(0.1,0.7,0.48,0.9)
     ChFF_leg->AddEntry("fChFF","New ^{3}He |F_{ch}(q^{2})| Fit","l");
     ChFF_leg->AddEntry("fChFF_Amroun","^{3}He |F_{ch}(q^{2})| Fit from Amroun et al. [4]","l");
     ChFF_leg->Draw();
     
     c5->cd(2);
     c5->cd(2)->SetLogy();
     //cout<<fMFF->Eval(35)<<endl;
     fMFF->Draw("L");
     
     //Now add Amroun's fit to the plot for comparison.
     /*for(Int_t i=0;i<ngaus;i++)
       {
       Qim[i] = Qim_Amroun[i];
       }*/
     TF1 *fMFF_Amroun = new TF1("fMFF_Amroun",MFF_Q2_Amroun,yminFF,ymaxFF+54,1);
     //cout<<fMFF_Amroun->Eval(30)<<endl;
     fMFF_Amroun->SetNpx(npdraw);
     fMFF_Amroun->SetLineColor(4);
     fMFF_Amroun->Draw("L same");
     auto MFF_leg = new TLegend(0.49,0.65,0.9,0.9); //(0.1,0.7,0.48,0.9)
     MFF_leg->AddEntry("fMFF","New ^{3}He |F_{m}(q^{2})| Fit","l");
     MFF_leg->AddEntry("fMFF_Amroun","^{3}He |F_{m}(q^{2})| Fit from Amroun et al. [4]","l");
     MFF_leg->Draw();
   }//End showplots.
 
  st->Stop();
  cout<<"*********************************************"<<endl;
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
/*
int main()
{
  Global_Fit_3He_SOG();
}
*/
