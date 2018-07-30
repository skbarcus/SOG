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

Int_t data_set = 0;        //Select which set of data points to use. Make code more clever later.

//Sets of data from Alex's 3He high Q^2 paper.
Double_t XS_Q2_55_1[2] = {2.77E-13,3.27E-15};        //fm^2/sr
Double_t uncertainty_Q2_55_1[2] = {0.39E-13,0.13E-15};                //fm^2/sr
Double_t E0_Q2_55_1[2] = {3.304, 0.9893};            //GeV
Double_t theta_Q2_55_1[2] = {27.24, 140.31};         //degrees
Double_t XS_Q2_60_8[2] = {2.14E-14,1.13E-15};        //fm^2/sr
Double_t uncertainty_Q2_60_8[2] = {0.72E-14,0.80E-15};                //fm^2/sr
Double_t E0_Q2_60_8[2] = {3.304, 1.052};             //GeV
Double_t theta_Q2_60_8[2] = {28.86, 140.51};         //degrees

Double_t XS_Q2_24_7[2] = {2.29E-9,2.80E-11};        //fm^2/sr
Double_t uncertainty_Q2_24_7[2] = {0.12E-9,0.20E-11};                //fm^2/sr
Double_t E0_Q2_24_7[2] = {3.304, 0.7391};             //GeV
Double_t theta_Q2_24_7[2] = {17.52, 97.78};         //degrees

Double_t XS_Q2_30_2[3] = {5.16E-10,3.95E-12,1.51E-12};        //fm^2/sr
Double_t uncertainty_Q2_30_2[3] = {0.29E-10,0.38E-12,0.19E-12};                //fm^2/sr
Double_t E0_Q2_30_2[3] = {3.304, 0.7391, 0.6878};             //GeV
Double_t theta_Q2_30_2[3] = {19.5, 118.99, 139.99};         //degrees

Double_t pi = 3.141592654;
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1.0/0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
//Double_t angle = 75.31*deg2rad;                   //Scattering angle [rad].
//Double_t theta = 75.31;

Double_t alpha = 1.0/137.0;              //Fine structure constant.
//Double_t E0 = 0.0;               //Initial electron energy [GeV].
Double_t Ef = 1.0;                       //Final energy of the electron after scattering.
Double_t MtH3 = 3.0160492*0.9315;       //Mass of trinucleon (H3 or He3) [GeV].
Double_t MtHe3 = 3.0160293*0.9315;
Double_t Q2 = 0.;                       //fm^-2

Double_t XSr[3] = {};
Double_t epsilon[3] = {};
Double_t uncertainty[3] = {};
Double_t tau = 0.;

void Rosenbluth_Separation() 
{

  if(data_set == 0)
    {
      Q2 = 55.1;             
    }
  if(data_set == 1)
    {
      Q2 = 60.8;
    }
  if(data_set == 2)
    {
      Q2 = 24.7;
    }
  if(data_set == 3)
    {
      Q2 = 30.2;
    }
  
  //Calculate tau.
  tau = Q2/(4.*pow(MtHe3,2.)*GeV2fm);     //convert to fm^-2
  cout<<"tau = "<<tau<<endl;

  //Calculate Mott cross section and remember to multiply by Z^2 and the recoil factor Ef/Ei. 
  Double_t MottXS(Double_t E0, Double_t theta)
  {
    Double_t mottxs = 0.; 
    Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
    
    mottxs = (  (pow(2,2)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

    cout<<"E0 = "<<E0<<"   theta = "<<theta<<endl;
    cout<<"Ef = "<<Ef<<"   mottxs = "<<mottxs<<endl;

    return mottxs;
  }
  
  Double_t Calc_Epsilon(Double_t theta)
  {
    Double_t val = 0.; 
    val = 1/(1+2*(1+tau)*pow(tan(theta*deg2rad/2.0),2.));
    return val;
  }
  
  Double_t XS_reduced(Double_t XSexp, Double_t E0, Double_t theta)
  {
    Double_t val = 0.; 
    val = ( XSexp / MottXS(E0,theta) ) * (1+tau);
    return val;
  }
  
  for(Int_t i=0.;i<3;i++)
    {
      if(data_set == 0)
	{
	  XSr[i] = XS_reduced(XS_Q2_55_1[i],E0_Q2_55_1[i],theta_Q2_55_1[i]);
	  epsilon[i] = Calc_Epsilon(theta_Q2_55_1[i]); 
	  uncertainty[i] = uncertainty_Q2_55_1[i]/MottXS(E0_Q2_55_1[i],theta_Q2_55_1[i]) * (1+tau);
	}

      if(data_set == 1)
	{
	  XSr[i] = XS_reduced(XS_Q2_60_8[i],E0_Q2_60_8[i],theta_Q2_60_8[i]);
	  epsilon[i] = Calc_Epsilon(theta_Q2_60_8[i]); 
	  uncertainty[i] = uncertainty_Q2_60_8[i]/MottXS(E0_Q2_60_8[i],theta_Q2_60_8[i]) * (1+tau);
	}

      if(data_set == 2)
	{
	  XSr[i] = XS_reduced(XS_Q2_24_7[i],E0_Q2_24_7[i],theta_Q2_24_7[i]);
	  epsilon[i] = Calc_Epsilon(theta_Q2_24_7[i]); 
	  uncertainty[i] = uncertainty_Q2_24_7[i]/MottXS(E0_Q2_24_7[i],theta_Q2_24_7[i]) * (1+tau);
	}

      if(data_set == 3)
	{
	  XSr[i] = XS_reduced(XS_Q2_30_2[i],E0_Q2_30_2[i],theta_Q2_30_2[i]);
	  epsilon[i] = Calc_Epsilon(theta_Q2_30_2[i]); 
	  uncertainty[i] = uncertainty_Q2_30_2[i]/MottXS(E0_Q2_30_2[i],theta_Q2_30_2[i]) * (1+tau);
	}
      cout<<"XSr["<<i<<"] = "<<XSr[i]<<"   epsilon["<<i<<"] = "<<epsilon[i]<<endl;
    }

  //Make Rosenbluth separation plot.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  TGraphErrors *graph = new TGraphErrors(2,epsilon,XSr,0,uncertainty);
  //Draw the new TGraph called graph on the canvas. 
  graph->Draw("");
  graph->SetLineWidth(1);
  graph->SetLineColor(4);
  graph->SetFillColor(0);
  graph->SetMarkerColor(1);
  graph->SetMarkerSize(0.4);
  graph->SetMarkerStyle(20);
  graph->SetTitle("Rosenbluth Separation Plot; #\epsilon; #\sigma_{r}");
  
}
