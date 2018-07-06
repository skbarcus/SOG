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

void Sieve_Aperture_Mott() 
{

  /*
  //Make a new canvas to plot data.
  TCanvas* c=new TCanvas("c");
  c->SetGrid();
  gROOT->Reset();
  */

  Int_t i = 0;
  Int_t j = 0;
  Double_t step = 0.0;
  Double_t pi = 3.141592654;
  Double_t deg2rad = pi/180.0;
  Double_t GeV2fm = 1.0/0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
  Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
  Double_t C = 299792458.0;                //Speed of light [m/s]. 
  Double_t fudge = 1.0;                 //Fudge factor for matching experimental data if I'm off by 4pi somehow from steradians.     

  Double_t lhrshorz = 0.03;                //Horizontal acceptance of LHRS [rad].
  Double_t lhrsangle = 12.5;               //LHRS angle [degrees].
  Double_t angle = 0.0;                    //Angle from target to point on sieve plate aperture [rad].
  Double_t dtarg2sieve = 1.109;            //Distance from target to sieve [m].
  Double_t dlhrshorz = 0.0;                //Width of LHRS acceptance on the sieve face [m].
  Double_t wquad = 0.043;                  //Width of the quadralateral sieve aperture [m].
  Double_t dquad = -0.0215;                //Distance from center of the quadralateral sieve aperture [m].
  Double_t hquad = 0.0;                    //Height of quadralateral aperture [m]. 
  Double_t quadangle = 0.0;                //Angular acceptance over the sieve aperture [rad]. 

  Double_t alpha = 1.0/137.0;              //Fine structure constant.
  Double_t E0 = 1.1;                       //Initial electron energy [GeV].
  Double_t Ef = 0.0;                       //Final energy of the electron after scattering.
  Double_t EfH3 = 0.0;
  Double_t EfHe3 = 0.0;
  Double_t Efproton = 0.0;
  Double_t mottxsection = 0.0;             //Mott cross section [1/GeV^2].

  Double_t Q2 = 0.0;                       //Q^2 in GeV.
  Double_t Q2H3 = 0.0;
  Double_t Q2He3 = 0.0;
  Double_t Q2effH3 = 0.0;                   //Q^2 with Coulomb corrections.
  Double_t Q2effHe3 = 0.0;
  Double_t Q2proton = 0.0;
  Double_t Q2effproton = 0.0;
  Double_t protondipole = 0.0;             //Reduced cross section for proton dipole. 

  Int_t npoints = 1000;
  Double_t x[10000];
  Double_t y[10000];

  Int_t nmottpoints = 1000;
  Double_t xmott[10000];
  Double_t ymott[10000];

  //Int_t ndipolepoints = 43;
  Double_t xdipole[10000];
  Double_t ydipole[10000];

  //Int_t nmottdipolepoints = 43;
  Double_t xmottdipole[10000];
  Double_t ymottdipole[10000];

  Int_t nmottwidepoints = 10000;
  Double_t xmottwide[10000];
  Double_t ymottwide[10000];

  Double_t xmottdipolewide[10000];
  Double_t ymottdipolewide[10000];

  Double_t Fch = 0.0;                    //Charge form factor.
  Double_t Fm = 0.0;                     //Magnetic form factor. 
  Double_t gamma = 0.8*pow(2.0/3.0,0.5);                  //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
  Double_t R[12] = {0.1, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9, 3.4, 4.0, 4.6, 5.2};  //Radii [fm].
  Double_t QH3ch[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121, 0.0, 0.0};
  Double_t QH3m[12] = {0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077, 0.0, 0.0};
  Double_t QHe3ch[12] = {0.027614, 0.170847, 0.219805, 0.170486, 0.134453, 0.100953, 0.074310, 0.053970, 0.023689, 0.017502, 0.002034, 0.004338};
  Double_t QHe3m[12] = {0.059785, 0.138368, 0.281326, 0.000037, 0.289808, 0.019056, 0.114825, 0.042296, 0.028345, 0.018312, 0.007843, 0.0};
  Double_t sumH3ch = 0.0;
  Double_t sumH3m = 0.0;
  Double_t sumHe3ch = 0.0;
  Double_t sumHe3m = 0.0;
  Double_t sumH3chtemp = 0.0;
  Double_t sumH3mtemp = 0.0;
  Double_t sumHe3chtemp = 0.0;
  Double_t sumHe3mtemp = 0.0;

  Double_t MtH3 = 3.0160492*0.9315;       //Mass of trinucleon (H3 or He3) [GeV].
  Double_t MtHe3 = 3.0160293*0.9315;
  Double_t Mtproton = 0.938272;           //Proton mass (for proton dipole) [GeV].      
  Double_t etaH3 = 0.0;                   //eta = 1+Q^2/(4*MT^2).
  Double_t etaHe3 = 0.0;
  Double_t muH3 = 2.9788*(3.0/1.0); //2.793-2*1.913 is too naive.    //Magnetic moment of trinucleon (H3 or He3). NIST: http://physics.nist.gov/cgi-bin/cuu/Results?search_for=magnet+moment   //MCEEP Code for H3 and He3 eleastic FFs has magnetic moments multiplied by 3.0/Z. I don't know why but it works. Maybe it's a factor of A/Z?
  Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.
  Double_t q2_3H3 = 0.0;                  //Momentum transfer squared three vector.
  Double_t q2_3He3 = 0.0;
  Double_t wH3 = 0.0;                     //Omega = E0 - E'.
  Double_t wHe3 = 0.0; 
  Double_t H3min = 1.0;                   //First minimum of cross section in fm^2.
  Double_t He3min = 1.0;
  Double_t H3minQ2 = 0.0;                 //Q^2 of first minimum of cross section in fm^-2.
  Double_t He3minQ2 = 0.0;
  Double_t rH3 = 0.0;                //Radius of H3 or He3 assmuing sphere of homogenous charge. 
  Double_t rHe3 = 0.0;

  Int_t nFFpoints = 1000;
  Double_t xFH3ch[10000];
  Double_t yFH3ch[10000];
  Double_t xFH3m[10000];
  Double_t yFH3m[10000];
  Double_t xFHe3ch[10000];
  Double_t yFHe3ch[10000];
  Double_t xFHe3m[10000];
  Double_t yFHe3m[10000];
  Double_t xcrosssectionH3[10000];
  Double_t ycrosssectionH3[10000];
  Double_t xcrosssectionHe3[10000];
  Double_t ycrosssectionHe3[10000];
  Double_t mottcrosssection[10000];
  Double_t mottcrosssectionH3[10000];
  Double_t mottcrosssectionHe3[10000];

  Double_t xH3crosssection[10000];
  Double_t yH3crosssection[10000];
  Double_t xHe3crosssection[10000];
  Double_t yHe3crosssection[10000];
  
  //Some experimental cross section data.
  //Tritium Experimental Data

  Double_t E189H3Q2[11] = {35.0, 40.0, 45.0, 50.0, 55.0, 65.0, 80.0, 85.0, 90.0, 95.0, 155.0};   //Angles to calculate Q^2 for this data.
  Double_t E189H3xsection[11] = {1.194*pow(10.0,-2.0), 6.030*pow(10.0,-3.0), 3.262*pow(10.0,-3.0), 1.958*pow(10.0,-3.0), 1.169*pow(10.0,-3.0), 4.900*pow(10.0,-4.0), 1.585*pow(10.0,-4.0), 1.013*pow(10.0,-4.0), 7.870*pow(10.0,-5.0), 5.325*pow(10.0,-5.0), 4.462*pow(10.0,-6.0)};   //xsections in mb/sr.
  Double_t E189H3Q2error[11] = {};
  Double_t E189H3xsectionerror[11] = {3.5, 3.4, 3.4, 3.5, 3.5, 3.4, 3.5, 3.4, 3.5, 3.4, 3.5};   //Percent error on cross section measurement. 

  Double_t E508H3Q2[18] = {30.0, 40.0, 45.0, 50.0, 54.0, 55.0, 60.0, 155.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 82.0, 90.0, 100.0, 155.0};   //Angles to calculate Q^2 for this data.
  Double_t E508H3xsection[18] = {8.670*pow(10.0,-4.0), 1.030*pow(10.0,-4.0), 3.752*pow(10.0,-5.0), 1.440*pow(10.0,-5.0), 6.789*pow(10.0,-6.0), 5.630*pow(10.0,-6.0), 2.328*pow(10.0,-6.0), 6.751*pow(10.0,-11.0), 1.031*pow(10.0,-4.0), 3.820*pow(10.0,-5.0), 1.420*pow(10.0,-5.0), 5.536*pow(10.0,-6.0), 2.202*pow(10.0,-6.0), 9.513*pow(10.0,-7.0), 6.949*pow(10.0,-8.0), 2.114*pow(10.0,-8.0), 6.417*pow(10.0,-9.0), 7.831*pow(10.0,-11.0)};   //xsections in mb/sr.
  Double_t E508H3Q2error[18] = {};
  Double_t E508H3xsectionerror[18] = {1.5, 1.5, 1.6, 1.5, 1.6, 1.6, 1.6, 23.2, 1.6, 1.6, 1.6, 1.8, 2.4, 2.5, 3.5, 4.1, 5.7, 21.4};   //Percent error on cross section measurement. 

  Double_t E683H3Q2[9] = {35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 66.0, 155.0, 104.0};   //Angles to calculate Q^2 for this data.
  Double_t E683H3xsection[9] = {4.244*pow(10.0,-5.0), 1.054*pow(10.0,-5.0), 2.814*pow(10.0,-6.0), 7.933*pow(10.0,-7.0), 2.275*pow(10.0,-7.0), 6.653*pow(10.0,-8.0), 1.703*pow(10.0,-8.0), 1.185*pow(10.0,-11.0), 5.600*pow(10.0,-11.0)};   //xsections in mb/sr.
  Double_t E683H3Q2error[9] = {};
  Double_t E683H3xsectionerror[9] = {1.8, 1.6, 1.5, 1.8, 2.5, 2.1, 3.2, 51.8, 50.0};   //Percent error on cross section measurement. 

  //Helium Experimental Data
  
  Double_t E314He3Q2[8] = {30.0, 35.0 ,40.0 ,45.0 ,50.0 ,55.0 ,60.0 ,155.0};   //Angles to calculate Q^2 for this data.
  Double_t E314He3xsection[8] = {1.923*pow(10.0,-2.0), 7.973*pow(10.0,-3.0), 3.486*pow(10.0,-3.0), 1.581*pow(10.0,-3.0), 7.511*pow(10.0,-4.0), 3.622*pow(10.0,-4.0), 1.835*pow(10.0,-4.0), 3.401*pow(10.0,-8.0)};   //xsections in mb/sr.
  Double_t E314He3Q2error[8] = {};
  Double_t E314He3xsectionerror[8] = {3.3, 3.2, 3.3, 3.3, 3.3, 3.3, 3.2, 4.1};     //Percent error on cross section measurement. 

  Double_t E413He3Q2[6] = {30.0, 35.0 ,40.0 ,45.0 ,50.0 ,155.0};   //Angles to calculate Q^2 for this data.
  Double_t E413He3xsection[6] = {6.547*pow(10.0,-3.0), 2.314*pow(10.0,-3.0), 8.687*pow(10.0,-4.0), 3.336*pow(10.0,-4.0), 1.320*pow(10.0,-4.0), 6.422*pow(10.0,-10.0)};   //xsections in mb/sr.
  Double_t E413He3Q2error[6] = {};
  Double_t E413He3xsectionerror[6] = {3.2, 3.3, 3.2, 3.5, 3.2, 7.8};     //Percent error on cross section measurement.
  
  Double_t E640He3Q2[16] = {45.0, 50.0, 54.0, 58.0, 62.0, 64.0, 66.0, 70.0, 74.0, 78.0, 82.0, 85.0, 94.0, 102.0, 112.0, 123.0};   //Angles to calculate Q^2 for this data.
  Double_t E640He3xsection[16] = {6.769*pow(10.0,-6.0), 1.510*pow(10.0,-6.0), 4.244*pow(10.0,-7.0), 1.116*pow(10.0,-7.0), 2.792*pow(10.0,-8.0), 1.476*pow(10.0,-8.0), 9.059*pow(10.0,-9.0), 4.876*pow(10.0,-9.0), 4.403*pow(10.0,-9.0), 3.735*pow(10.0,-9.0), 3.091*pow(10.0,-9.0), 2.623*pow(10.0,-9.0), 1.290*pow(10.0,-9.0), 7.627*pow(10.0,10.0), 2.780*pow(10.0,-10.0), 8.129*pow(10.0,-11.0),};   //xsections in mb/sr.
  Double_t E640He3Q2error[16] = {};
  Double_t E640He3xsectionerror[16] = {3.8, 3.9, 6.8, 7.2, 4.8, 6.4, 6.4, 6.5, 9.1, 9.5, 8.1, 7.3, 8.4, 13.9, 15.1, 29.7};     //Percent error on cross section measurement.

  //Tritium Experimental Data

  //Convert H3 E189 units to my units.
  if(E0 == 0.1892)
    {
      for(i=0; i<11; i++)
	{
	  E189H3Q2[i] = 4 * 0.1892 * (0.1892/(1.0+2.0*0.1892*pow(sin(E189H3Q2[i]*deg2rad/2.0),2.0)/MtH3)) * pow(sin(E189H3Q2[i]*deg2rad/2.0),2.0) * GeV2fm;             //Calculate Q^2 in fm^-2
	  E189H3xsection[i] = (E189H3xsection[i]/10.0)*fudge;         //Convert mb/sr to fm^2.
	  
	  E189H3xsectionerror[i] = (E189H3xsectionerror[i]/100.0)*E189H3xsection[i];  //Calculate size of error bars.
	  
	  cout<<"E189H3Q2 = "<<E189H3Q2[i]<<"   E189H3xsection = "<<E189H3xsection[i]<<endl;
	}
    }

  //Convert H3 E508 units to my units.
  if(E0 == 0.5084)
    {
      for(i=0; i<18; i++)
	{
	  E508H3Q2[i] = 4 * 0.5084 * (0.5084/(1.0+2.0*0.5084*pow(sin(E508H3Q2[i]*deg2rad/2.0),2.0)/MtH3)) * pow(sin(E508H3Q2[i]*deg2rad/2.0),2.0) * GeV2fm;             //Calculate Q^2 in fm^-2
	  E508H3xsection[i] = (E508H3xsection[i]/10.0)*fudge;         //Convert mb/sr to fm^2.
	  
	  E508H3xsectionerror[i] = (E508H3xsectionerror[i]/100.0)*E508H3xsection[i];  //Calculate size of error bars.
	  
	  cout<<"E508H3Q2 = "<<E508H3Q2[i]<<"   E508H3xsection = "<<E508H3xsection[i]<<endl;
	}
    }

  //Convert H3 E683 units to my units.
  if(E0 == 0.6837)
    {
      for(i=0; i<9; i++)
	{
	  E683H3Q2[i] = 4 * 0.6837 * (0.6837/(1.0+2.0*0.6837*pow(sin(E683H3Q2[i]*deg2rad/2.0),2.0)/MtH3)) * pow(sin(E683H3Q2[i]*deg2rad/2.0),2.0) * GeV2fm;             //Calculate Q^2 in fm^-2
	  E683H3xsection[i] = (E683H3xsection[i]/10.0)*fudge;         //Convert mb/sr to fm^2.
	  
	  E683H3xsectionerror[i] = (E683H3xsectionerror[i]/100.0)*E683H3xsection[i];  //Calculate size of error bars.
	  
	  cout<<"E683H3Q2 = "<<E683H3Q2[i]<<"   E683H3xsection = "<<E683H3xsection[i]<<endl;
	}
    }



  //Helium Experimental Data
  
  //Convert He3 E314 units to my units. 
  if(E0 == 0.3144)
    {
      for(i=0; i<8; i++)
	{
	  //EfHe3 = 0.3144/(1.0+2.0*0.3144*pow(sin(E314Q2[i]*deg2rad/2.0),2.0)/MtHe3);         //Calculate Ef for angle.
	  E314He3Q2[i] = 4 * 0.3144 * (0.3144/(1.0+2.0*0.3144*pow(sin(E314He3Q2[i]*deg2rad/2.0),2.0)/MtHe3)) * pow(sin(E314He3Q2[i]*deg2rad/2.0),2.0) * GeV2fm;             //Calculate Q^2 in fm^-2
	  E314He3xsection[i] = (E314He3xsection[i]/10.0)*fudge;         //Convert mb/sr to fm^2.
	  
	  E314He3xsectionerror[i] = (E314He3xsectionerror[i]/100.0)*E314He3xsection[i];  //Calculate size of error bars.
	  
	  cout<<"E314He3Q2 = "<<E314He3Q2[i]<<"   E314He3xsection = "<<E314He3xsection[i]<<endl;
	}
    }
  
  //Convert He3 E640 units to my units.
  if(E0 == 0.64)
    {
      for(i=0; i<8; i++)
	{
	  E640He3Q2[i] = 4 * 0.640 * (0.640/(1.0+2.0*0.640*pow(sin(E640He3Q2[i]*deg2rad/2.0),2.0)/MtHe3)) * pow(sin(E640He3Q2[i]*deg2rad/2.0),2.0) * GeV2fm;             //Calculate Q^2 in fm^-2
	  E640He3xsection[i] = (E640He3xsection[i]/10.0)*fudge;         //Convert mb/sr to fm^2.
	  
	  E640He3xsectionerror[i] = (E640He3xsectionerror[i]/100.0)*E640He3xsection[i];  //Calculate size of error bars.
	  
	  cout<<"E640He3Q2 = "<<E640He3Q2[i]<<"   E640He3xsection = "<<E640He3xsection[i]<<endl;
	}
    }

  //Convert He3 E413 units to my units.
  if(E0 == 0.4133)
    {
      for(i=0; i<6; i++)
	{
	  E413He3Q2[i] = 4 * 0.4133 * (0.4133/(1.0+2.0*0.4133*pow(sin(E413He3Q2[i]*deg2rad/2.0),2.0)/MtHe3)) * pow(sin(E413He3Q2[i]*deg2rad/2.0),2.0) * GeV2fm;             //Calculate Q^2 in fm^-2
	  E413He3xsection[i] = (E413He3xsection[i]/10.0)*fudge;         //Convert mb/sr to fm^2.
	  
	  E413He3xsectionerror[i] = (E413He3xsectionerror[i]/100.0)*E413He3xsection[i];  //Calculate size of error bars.
	  
	  cout<<"E413He3Q2 = "<<E413He3Q2[i]<<"   E413He3xsection = "<<E413He3xsection[i]<<endl;
	}
    }
  

  //Calculate width of LHRS acceptance on the sieve face [m].
  dlhrshorz = 2.0*dtarg2sieve*sin(0.03);
  //cout<<"dlhrshorz = "<<dlhrshorz<<endl;

  //Calculate angular acceptance over the sieve aperture [rad].
  quadangle = (wquad/dlhrshorz)*0.03;
  //cout<<"quadangle = "<<quadangle<<endl;

  //Calculate Mott cross section for LHRS angle (test).
  mottxsection = (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((lhrsangle*deg2rad)/2.0),4.0)))*pow(cos((lhrsangle*deg2rad)/2.0),2.0);
  //cout<<"Mott cross section = "<<mottxsection<<endl;

  //Make a new canvas to plot data.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  gROOT->Reset();

  for(i=0; i<nmottwidepoints; i++)
    {
      //Starting on the left side of the sieve aperture track across to the right.
      step = i;    //Convert i to double. 
      angle = 0 + (step/nmottwidepoints)*180.0*deg2rad;
      mottxsection = (  (pow(alpha,2)/(4.0*pow(E0,2)*pow(sin((angle)/2),4)))*pow(cos((angle)/2),2)  ) * 1.0/25.7;   //Now in fm^2.
      Q2 = (  4.0*pow(E0,2.0)*pow(sin(angle/2.0),2.0)  ) * GeV2fm;    //fm^-2
      protondipole = pow((1+(Q2*1.0/GeV2fm)/0.71),-2.0);              //Unitless.

      //cout<<"angle = "<<angle<<"   dquad = "<<dquad<<"   hquad = "<<hquad<<endl;
      //cout<<"Mott Cross Section = "<<mottxsection<<endl;

      xmottwide[i] = Q2;
      ymottwide[i] = mottxsection;

      xmottdipolewide[i] = Q2;
      ymottdipolewide[i] = mottxsection*protondipole;

      //cout<<"x["<<i<<"] = "<<x[i]<<"   y["<<i<<"] = "<<y[i]<<endl;
      //dquad = dquad + 0.001;
    }

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph4 = new TGraph(nmottwidepoints,xmottwide,ymottwide);
  //Draw the new TGraph called graph on the canvas. 
  graph4->Draw();
  
  //Set X axis
  //c1->SetLogy();
  //graph4->GetXaxis()->SetLimits(0.0,1.5);
  //graph4->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph4->SetMinimum(pow(10.0,-4.0));
  graph4->SetMaximum(1.0);
  graph4->SetLineWidth(3);
  graph4->SetLineColor(2);
  graph4->SetFillColor(0);
  //gPad->SetLogx(1);
  //gPad->SetLogy(1);
  c1->SetLogy();
  graph4->GetXaxis()->SetRangeUser(0.0,pi); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph4->SetTitle("Mott Cross Section with and without Proton Dipole FF; Q^2 [fm^-2]; Cross Section [fm^2]");
  //graph_expected.SetFillColor(kYellow);
  //graph_expected.DrawClone("E3AL"); // E3 draws the band

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph5 = new TGraph(nmottwidepoints,xmottdipolewide,ymottdipolewide);
  //Draw the new TGraph called graph on the canvas. 
  graph5->Draw("same");
  graph5->SetLineWidth(3);
  graph5->SetLineColor(4);
  graph5->SetFillColor(0);

  // Draw the Legend
  TLegend leg1(0.9,.7,.56,.9,"");
  leg1.SetFillColor(0);
  leg1.AddEntry(graph4,"Mott Cross Section");
  leg1.AddEntry(graph5,"Mott Cross Section * Proton Dipole FF");
  leg1.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  //Calculate charge FFs.
  for(i=0; i<nFFpoints; i++)
    {
      step = i;    //Convert i to double. 
      angle = 0.00000001 + (step/nFFpoints)*180.0*deg2rad;    //Note: 0.00000001 offset is just to avoid Q^2=0 at angle=0.
      EfH3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtH3);
      EfHe3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtHe3);
      Q2H3 = 4.0*E0*EfH3*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2He3 = 4.0*E0*EfHe3*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2effH3 = pow( pow(Q2H3,0.5) * (1.0+(1.5*1*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);   //Z=1
      Q2effHe3 = pow( pow(Q2He3,0.5) * (1.0+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);  //Z=2

      //Calculate the sum part of the SOG paramaterization for H3 charge FF.
      for(j=0; j<12; j++)
	{
	  sumH3chtemp = (QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effH3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effH3,0.5)*R[j])/(pow(Q2effH3,0.5)*R[j])) );
	  //sumH3chtemp = QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0));   //Fine.
	  //sumH3chtemp = cos(pow(Q2,0.5)*R[j]);   //Fine.
	  //sumH3chtemp = 2.0*pow(R[j],2.0)/pow(gamma,2.0);   //Fine.
	  //sumH3chtemp = (sin(pow(Q2,0.5)*R[j])/(pow(Q2,0.5)*R[j]));   //NAN. Issue is at angle = 0 -> Q = 0.
	  sumH3ch = sumH3ch + sumH3chtemp;
	  //cout<<"sumH3ch = "<<sumH3ch<<endl;
	}

      //Calculate the sum part of the SOG paramaterization for He3 charge FF.
      for(j=0; j<12; j++)
	{
	  sumHe3chtemp = (QHe3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effHe3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effHe3,0.5)*R[j])/(pow(Q2effHe3,0.5)*R[j])) );
	  sumHe3ch = sumHe3ch + sumHe3chtemp;
	  //cout<<"sumH3ch = "<<sumH3ch<<endl;
	}
      
      xFH3ch[i] = Q2H3;
      yFH3ch[i] = exp(-(1.0/4.0)*Q2effH3*pow(gamma,2.0))*fabs(sumH3ch);

      xFHe3ch[i] = Q2He3;
      yFHe3ch[i] = exp(-(1.0/4.0)*Q2effHe3*pow(gamma,2.0))*fabs(sumHe3ch);
      
      sumH3ch = 0.0;      //Reset sumH3ch.
      sumHe3ch = 0.0;     //Reset sumHe3ch.

      //cout<<"Q^2 = "<<xFH3ch[i]<<"   FH3ch = "<<yFH3ch[i]<<endl;

      /*
      if(Q2H3<25.0)
	{
	  if(yFH3ch[i]<H3min && Q2H3<15.0)
	    {
	      H3min = yFH3ch[i];
	      H3minQ2 = Q2H3;
	    }

	  if(yFHe3ch[i]<He3min && Q2H3<13.0)
	    {
	      He3min = yFHe3ch[i];
	      He3minQ2 = Q2He3;
	    }
	}
      */

    }

  /*
  rH3 = (4.5*hbar)/pow(((H3minQ2*0.0389*pow(10.0,9.0))/pow(C,1.0)),0.5);
  rHe3 = (4.5*hbar)/pow(((He3minQ2*0.0389*pow(10.0,9.0))/pow(C,1.0)),0.5);
  */

  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  gROOT->Reset();
  c2->Divide(2,1);

  c2->cd(1);
  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph6 = new TGraph(nFFpoints,xFH3ch,yFH3ch);
  //Draw the new TGraph called graph on the canvas. 
  graph6->Draw();
  
  //Set X axis
  graph6->GetXaxis()->SetLimits(0.0,25.0);
  //graph6->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph6->SetMinimum(pow(10.0,-4.0));
  //graph6->SetMinimum(-1.0);
  graph6->SetMaximum(1.0);
  //graph6->GetYaxis()->SetUserRange(-1.0,1.0);
  graph6->SetLineWidth(3);
  graph6->SetLineColor(2);
  graph6->SetFillColor(0);
  c2->cd(1)->SetLogy();
  //graph6->GetXaxis()->SetRangeUser(0.0,25.0); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph6->SetTitle("H3 Charge FF; Q^2 [fm^-2]; H3 Fch(q)");

  // Draw the Legend
  TLegend leg2(0.9,.7,.56,.9,"");
  leg2.SetFillColor(0);
  leg2.AddEntry(graph6,"H3 Charge FF");
  leg2.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  c2->cd(2);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph7 = new TGraph(nFFpoints,xFHe3ch,yFHe3ch);
  //Draw the new TGraph called graph on the canvas. 
  graph7->Draw();
  
  //Set X axis
  graph7->GetXaxis()->SetLimits(0.0,25.0);
  //graph6->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph7->SetMinimum(pow(10.0,-4.0));
  //graph6->SetMinimum(-1.0);
  graph7->SetMaximum(1.0);
  //graph6->GetYaxis()->SetUserRange(-1.0,1.0);
  graph7->SetLineWidth(3);
  graph7->SetLineColor(2);
  graph7->SetFillColor(0);
  c2->cd(2)->SetLogy();
  //graph6->GetXaxis()->SetRangeUser(0.0,25.0); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph7->SetTitle("He3 Charge FF; Q^2 [fm^-2]; He3 Fch(q)");

  // Draw the Legend
  TLegend leg2(0.9,.7,.56,.9,"");
  leg2.SetFillColor(0);
  leg2.AddEntry(graph7,"He3 Charge FF");
  leg2.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  //Calculate magnetic FFs.
  for(i=0; i<nFFpoints; i++)
    {
      step = i;    //Convert i to double. 
      angle = 0.00000001 + (step/nFFpoints)*180.0*deg2rad;    //Note: 0.00000001 offset is just to avoid Q^2=0 at angle=0.
      //cout<<angle*(1.0/deg2rad)<<endl;
      EfH3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtH3);
      EfHe3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtHe3);
      Q2H3 = 4.0*E0*EfH3*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2He3 = 4.0*E0*EfHe3*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2effH3 = pow( pow(Q2H3,0.5) * (1.0+(1.5*1*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);   //Z=1
      Q2effHe3 = pow( pow(Q2He3,0.5) * (1.0+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);  //Z=2

      //Calculate the sum part of the SOG paramaterization for H3 magnetic FF.
      for(j=0; j<12; j++)
	{
	  sumH3mtemp = (QH3m[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effH3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effH3,0.5)*R[j])/(pow(Q2effH3,0.5)*R[j])) );
	  //sumH3chtemp = QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0));   //Fine.
	  //sumH3chtemp = cos(pow(Q2,0.5)*R[j]);   //Fine.
	  //sumH3chtemp = 2.0*pow(R[j],2.0)/pow(gamma,2.0);   //Fine.
	  //sumH3chtemp = (sin(pow(Q2,0.5)*R[j])/(pow(Q2,0.5)*R[j]));   //NAN. Issue is at angle = 0 -> Q = 0.
	  sumH3m = sumH3m + sumH3mtemp;
	  //cout<<"sumH3m = "<<sumH3m<<endl;
	}

      //Calculate the sum part of the SOG paramaterization for He3 magnetic FF.
      for(j=0; j<12; j++)
	{
	  sumHe3mtemp = (QHe3m[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effHe3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effHe3,0.5)*R[j])/(pow(Q2effHe3,0.5)*R[j])) );
	  sumHe3m = sumHe3m + sumHe3mtemp;
	  //cout<<"sumH3m = "<<sumH3m<<endl;
	}
      
      //Calculate FFs.
      xFH3m[i] = Q2H3;
      yFH3m[i] = exp(-(1.0/4.0)*Q2effH3*pow(gamma,2.0))*fabs(sumH3m);

      xFHe3m[i] = Q2He3;
      yFHe3m[i] = exp(-(1.0/4.0)*Q2effHe3*pow(gamma,2.0))*fabs(sumHe3m);
      
      sumH3m = 0.0;      //Reset sumH3ch.
      sumHe3m = 0.0;     //Reset sumHe3ch.

      //cout<<"Q^2 = "<<xFH3m[i]<<"   FH3m = "<<yFH3m[i]<<endl;

      //Calculation of H3 and He3 cross sections.
      wH3 = E0 - EfH3;                  
      wHe3 = E0 - EfHe3;
      //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
      q2_3H3 = fabs(  pow(wH3,2.0)*GeV2fm - Q2effH3  );                  //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
      q2_3He3 = fabs(  pow(wHe3,2.0)*GeV2fm - Q2effHe3  );
      etaH3 = 1.0 + Q2effH3/(4.0*pow(MtH3,2.0)*GeV2fm);        //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.
      etaHe3 = 1.0 + Q2effHe3/(4.0*pow(MtHe3,2.0)*GeV2fm); 

      //mottcrosssection[i] = (  (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle/2.0),4.0)))*pow(cos(angle/2.0),2.0)  ) * 1.0/25.7;   //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
      //Calculate Mott cross section and remember to multiply by Z^2 and the recoil factor Ef/Ei. 
      mottcrosssectionH3[i] = (  (pow(1,2)*(EfH3/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle/2.0),4.0)))*pow(cos(angle/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
      mottcrosssectionHe3[i] = (  (pow(2,2)*(EfHe3/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle/2.0),4.0)))*pow(cos(angle/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

      xcrosssectionH3[i] = Q2H3;
      ycrosssectionH3[i] = fabs( mottcrosssectionH3[i]*(1.0/etaH3) * (  (Q2effH3/q2_3H3)*pow(yFH3ch[i],2.0) + (pow(muH3,2.0)*Q2effH3/(2.0*pow(MtH3,2.0)*GeV2fm)) * (Q2effH3/(2.0*q2_3H3) + pow(tan(angle/2.0),2.0)) * pow(yFH3m[i],2.0)  ) );
      xcrosssectionHe3[i] = Q2He3;
      ycrosssectionHe3[i] = fabs( mottcrosssectionHe3[i]*(1.0/etaHe3) * (  (Q2effHe3/q2_3He3)*pow(yFHe3ch[i],2.0) + (pow(muHe3,2.0)*Q2effHe3/(2.0*pow(MtHe3,2.0)*GeV2fm)) * (Q2effHe3/(2.0*q2_3He3) + pow(tan(angle/2.0),2.0)) * pow(yFHe3m[i],2.0)  ) );

      /* 
      if(Q2H3<40.0)
	{
	  if(ycrosssectionH3[i]<H3min && Q2H3<17.0)
	    {
	      H3min = ycrosssectionH3[i];
	      H3minQ2 = Q2H3;
	    }

	  if(ycrosssectionHe3[i]<He3min && Q2H3<13.0)
	    {
	      He3min = ycrosssectionHe3[i];
	      He3minQ2 = Q2He3;
	    }
	}
      */

      //cout<<"q2_3H3 = "<<q2_3H3<<"   q2_3He3 = "<<q2_3He3<<endl;
      //cout<<"Q2H3 = "<<Q2H3<<"   Q2effH3 = "<<Q2effH3<<"   Q2He3 = "<<Q2He3<<"   Q2effHe3 = "<<Q2effHe3<<endl;
      //cout<<"H3min = "<<H3min<<"   H3minQ2 = "<<H3minQ2<<"   He3min = "<<He3min<<"   He3minQ2 = "<<He3minQ2<<endl;
      //cout<<"FH3ch^2 = "<<pow(yFH3ch[i],2.0)<<"   FHe3ch^2 = "<<pow(yFHe3ch[i],2.0)<<"   FH3m^2 = "<<pow(yFH3m[i],2.0)<<"   FHe3m^2 = "<<pow(yFHe3m[i],2.0)<<endl;
      //cout<<"wH3 = "<<wH3<<"   wHe3 = "<<wHe3<<"   q2_3H3 = "<<q2_3H3<<"   q2_3He3 = "<<q2_3He3<<"   etaH3 = "<<etaH3<<"   etaHe3 = "<<etaHe3<<endl;
      //cout<<"Q^2 = "<<Q2H3<<"   Mott Cross Section = "<<mottcrosssection[i]<<"   H3 Cross Section = "<<ycrosssectionH3[i]<<"   He3 Cross Section = "<<ycrosssectionHe3[i]<<endl;
    }

  //Calculate radius of H3 and He3 using the first minimum of their cross sections and assuming that they are homogenous spheres. Povh, Rith pg. 64 Examples of Form Factors. 
  
  /*
  rH3 = (4.5*hbar)/pow(((H3minQ2*0.0389*pow(10.0,9.0))/pow(C,1.0)),0.5);
  rHe3 = (4.5*hbar)/pow(((He3minQ2*0.0389*pow(10.0,9.0))/pow(C,1.0)),0.5);

  cout<<"H3min = "<<H3min<<" fm^2   H3minQ2 = "<<H3minQ2<<" fm^-2   He3min = "<<He3min<<" fm^2   He3minQ2 = "<<He3minQ2<<" fm^-2"<<endl;
  cout<<"rH3 = "<<rH3<<" m   rHe3 =  "<<rHe3<<" m"<<endl;
  */

 //Make a new canvas to plot data.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  gROOT->Reset();
  c3->Divide(2,1);

  c3->cd(1);
  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph8 = new TGraph(nFFpoints,xFH3m,yFH3m);
  //Draw the new TGraph called graph on the canvas. 
  graph8->Draw();
  
  //Set X axis
  graph8->GetXaxis()->SetLimits(0.0,40.0);
  //graph6->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph8->SetMinimum(pow(10.0,-4.0));
  //graph6->SetMinimum(-1.0);
  graph8->SetMaximum(1.0);
  //graph6->GetYaxis()->SetUserRange(-1.0,1.0);
  graph8->SetLineWidth(3);
  graph8->SetLineColor(2);
  graph8->SetFillColor(0);
  c3->cd(1)->SetLogy();
  //graph6->GetXaxis()->SetRangeUser(0.0,25.0); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph8->SetTitle("H3 Magnetic FF; Q^2 [fm^-2]; H3 Fm(q)");

  // Draw the Legend
  TLegend leg3(0.9,.7,.56,.9,"");
  leg3.SetFillColor(0);
  leg3.AddEntry(graph8,"H3 Magnetic FF");
  leg3.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  c3->cd(2);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph9 = new TGraph(nFFpoints,xFHe3m,yFHe3m);
  //Draw the new TGraph called graph on the canvas. 
  graph9->Draw();
  
  //Set X axis
  graph9->GetXaxis()->SetLimits(0.0,40.0);
  //graph6->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph9->SetMinimum(pow(10.0,-4.0));
  //graph6->SetMinimum(-1.0);
  graph9->SetMaximum(1.0);
  //graph6->GetYaxis()->SetUserRange(-1.0,1.0);
  graph9->SetLineWidth(3);
  graph9->SetLineColor(2);
  graph9->SetFillColor(0);
  c3->cd(2)->SetLogy();
  //graph6->GetXaxis()->SetRangeUser(0.0,25.0); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph9->SetTitle("He3 Magnetic FF; Q^2 [fm^-2]; H3 Fm(q)");

  // Draw the Legend
  TLegend leg4(0.9,.7,.56,.9,"");
  leg4.SetFillColor(0);
  leg4.AddEntry(graph7,"He3 Magnetic FF");
  leg4.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.


//Make a new canvas to plot data.
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();
  gROOT->Reset();
  c4->Divide(2,1);

  c4->cd(1);
  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph10 = new TGraph(nFFpoints,xcrosssectionH3,ycrosssectionH3);
  //Draw the new TGraph called graph on the canvas. 
  graph10->Draw();
  
  //Set X axis
  graph10->GetXaxis()->SetLimits(0.0,40.0);
  //graph6->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph10->SetMinimum(pow(10.0,-13.0));
  //graph6->SetMinimum(-1.0);
  graph10->SetMaximum(pow(10.0,-1.0));
  //graph6->GetYaxis()->SetUserRange(-1.0,1.0);
  graph10->SetLineWidth(3);
  graph10->SetLineColor(2);
  graph10->SetFillColor(0);
  c4->cd(1)->SetLogy();
  //graph6->GetXaxis()->SetRangeUser(0.0,25.0); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph10->SetTitle("H3 Cross Section; Q^2 [fm^-2]; dsig/domega [fm^2/sr]");

  //Plot tritium experimental data.

  if(E0 == 0.1892)
    {
      //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
      graph14 = new TGraphErrors(11,E189H3Q2,E189H3xsection,E189H3Q2error,E189H3xsectionerror);
      graph14->Draw("psame");
      graph14->SetMarkerStyle(8);
      graph14->SetMarkerSize(0.75);
    }  

  if(E0 == 0.5084)
    {
      //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
      graph14 = new TGraphErrors(18,E508H3Q2,E508H3xsection,E508H3Q2error,E508H3xsectionerror);
      graph14->Draw("psame");
      graph14->SetMarkerStyle(8);
      graph14->SetMarkerSize(0.75);
    }  

  if(E0 == 0.6837)
    {
      //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
      graph14 = new TGraphErrors(9,E683H3Q2,E683H3xsection,E683H3Q2error,E683H3xsectionerror);
      graph14->Draw("psame");
      graph14->SetMarkerStyle(8);
      graph14->SetMarkerSize(0.75);
    }  


  // Draw the Legend
  TLegend leg5(0.9,.7,.56,.9,"");
  leg5.SetFillColor(0);
  leg5.AddEntry(graph10,"H3 Cross Section");
  leg5.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  c4->cd(2);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph11 = new TGraph(nFFpoints,xcrosssectionHe3,ycrosssectionHe3);
  //Draw the new TGraph called graph on the canvas. 
  graph11->Draw();
  
  //Set X axis
  graph11->GetXaxis()->SetLimits(0.0,40.0);
  //graph6->GetXaxis()->SetRange(0.0,1.5);
  //Set Y axis Min and Max (not sure why different from X).
  graph11->SetMinimum(pow(10.0,-13.0));
  //graph6->SetMinimum(-1.0);
  graph11->SetMaximum(pow(10.0,-1.0));
  //graph6->GetYaxis()->SetUserRange(-1.0,1.0);
  graph11->SetLineWidth(3);
  graph11->SetLineColor(2);
  graph11->SetFillColor(0);
  c4->cd(2)->SetLogy();
  //graph6->GetXaxis()->SetRangeUser(0.0,25.0); //Need SetRangeUser to set x axis if y axis is log for some reason.
  graph11->SetTitle("He3 Cross Section; Q^2 [fm^-2]; dsig/domega [fm^2/sr]");

  /*
  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph14 = new TGraph(8,E314He3Q2,E314He3xsection);
  graph14->Draw("psame");
  graph14->SetMarkerStyle(8);
  graph14->SetMarkerSize(0.75);
*/

  //Plot helium experimental data.

  if(E0 == 0.3144)
    {
      //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
      graph15 = new TGraphErrors(8,E314He3Q2,E314He3xsection,E314He3Q2error,E314He3xsectionerror);
      graph15->Draw("psame");
      graph15->SetMarkerStyle(8);
      graph15->SetMarkerSize(0.75);
    }

  if(E0 == 0.640)
    {
      //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
      graph15 = new TGraphErrors(16,E640He3Q2,E640He3xsection,E640He3Q2error,E640He3xsectionerror);
      graph15->Draw("psame");
      graph15->SetMarkerStyle(8);
      graph15->SetMarkerSize(0.75);
    }  

  if(E0 == 0.4133)
    {
      //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
      graph15 = new TGraphErrors(6,E413He3Q2,E413He3xsection,E413He3Q2error,E413He3xsectionerror);
      graph15->Draw("psame");
      graph15->SetMarkerStyle(8);
      graph15->SetMarkerSize(0.75);
    }  
  
  // Draw the Legend
  TLegend leg6(0.9,.7,.56,.9,"");
  leg6.SetFillColor(0);
  leg6.AddEntry(graph11,"He3 Cross Section");
  //leg6.AddEntry(graph15,"Experimental Data E=0.3144 GeV");
  leg6.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.


  //Make a new canvas to plot data.
  TCanvas* c=new TCanvas("c");
  c->SetGrid();
  gROOT->Reset();

  //Calculate the Mott cross section periodically across the sieve aperture and weight the result by the aperture height to see if rate is constant over aperture. 

  for(i=0; i<npoints; i++)
    {
      //Starting on the left side of the sieve aperture track across to the right.
      step = i;    //Convert i to double. 
      angle = lhrsangle*deg2rad-0.5*quadangle+(step/npoints)*quadangle;
      dquad = (angle-lhrsangle*deg2rad)*(0.06653/0.03)*-1.0; //Multiply by -1 to make higher rates on right side of sieve plate from beam's prespective. 
      hquad = (-0.008/0.043)*(dquad)+0.01;
      mottxsection = (  (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((angle)/2.0),4.0)))*pow(cos((angle)/2.0),2.0)  ) * 1.0/25.7;   //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.


      EfH3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtH3);
      EfHe3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtHe3);
      Q2H3 = 4.0*E0*EfH3*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2He3 = 4.0*E0*EfHe3*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2effH3 = pow( pow(Q2H3,0.5) * (1.0+(1.5*1*(1/137.0))/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);   //Z=1
      Q2effHe3 = pow( pow(Q2He3,0.5) * (1.0+(1.5*2*(1/137.0))/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);  //Z=2


      //Q2 = 4.0*E0*EfH3*pow(sin(angle/2.0),2.0) * GeV2fm;              //fm^-2

      Efproton = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/Mtproton);
      Q2proton = 4.0*E0*Efproton*pow(sin(angle/2.0),2.0) * GeV2fm;
      Q2effproton = pow( pow(Q2proton,0.5) * (1.0+(1.5*1*(1/137.0))/(E0*pow(GeV2fm,0.5)*1.12*pow(1.0,1.0/3.0))) ,2.0);
      protondipole = pow((1+(Q2proton*1.0/GeV2fm)/0.71),-2.0);                  //Unitless.




      //Calculate the sum part of the SOG paramaterization for H3 charge FF.
      for(j=0; j<12; j++)
	{
	  sumH3chtemp = (QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effH3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effH3,0.5)*R[j])/(pow(Q2effH3,0.5)*R[j])) );
	  //sumH3chtemp = QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0));   //Fine.
	  //sumH3chtemp = cos(pow(Q2,0.5)*R[j]);   //Fine.
	  //sumH3chtemp = 2.0*pow(R[j],2.0)/pow(gamma,2.0);   //Fine.
	  //sumH3chtemp = (sin(pow(Q2,0.5)*R[j])/(pow(Q2,0.5)*R[j]));   //NAN. Issue is at angle = 0 -> Q = 0.
	  sumH3ch = sumH3ch + sumH3chtemp;
	  //cout<<"sumH3ch = "<<sumH3ch<<endl;
	}

      //Calculate the sum part of the SOG paramaterization for He3 charge FF.
      for(j=0; j<12; j++)
	{
	  sumHe3chtemp = (QHe3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effHe3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effHe3,0.5)*R[j])/(pow(Q2effHe3,0.5)*R[j])) );
	  sumHe3ch = sumHe3ch + sumHe3chtemp;
	  //cout<<"sumH3ch = "<<sumH3ch<<endl;
	}
      
      xFH3ch[i] = Q2H3;
      yFH3ch[i] = exp(-(1.0/4.0)*Q2effH3*pow(gamma,2.0))*fabs(sumH3ch);

      xFHe3ch[i] = Q2He3;
      yFHe3ch[i] = exp(-(1.0/4.0)*Q2effHe3*pow(gamma,2.0))*fabs(sumHe3ch);
      
      sumH3ch = 0.0;      //Reset sumH3ch.
      sumHe3ch = 0.0;     //Reset sumHe3ch.

      //Calculate the sum part of the SOG paramaterization for H3 magnetic FF.
      for(j=0; j<12; j++)
	{
	  sumH3mtemp = (QH3m[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effH3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effH3,0.5)*R[j])/(pow(Q2effH3,0.5)*R[j])) );
	  //sumH3chtemp = QH3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0));   //Fine.
	  //sumH3chtemp = cos(pow(Q2,0.5)*R[j]);   //Fine.
	  //sumH3chtemp = 2.0*pow(R[j],2.0)/pow(gamma,2.0);   //Fine.
	  //sumH3chtemp = (sin(pow(Q2,0.5)*R[j])/(pow(Q2,0.5)*R[j]));   //NAN. Issue is at angle = 0 -> Q = 0.
	  sumH3m = sumH3m + sumH3mtemp;
	  //cout<<"sumH3m = "<<sumH3m<<endl;
	}

      //Calculate the sum part of the SOG paramaterization for He3 magnetic FF.
      for(j=0; j<12; j++)
	{
	  sumHe3mtemp = (QHe3m[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effHe3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effHe3,0.5)*R[j])/(pow(Q2effHe3,0.5)*R[j])) );
	  sumHe3m = sumHe3m + sumHe3mtemp;
	  //cout<<"sumH3m = "<<sumH3m<<endl;
	}
      
      //Calculate FFs.
      xFH3m[i] = Q2H3;
      yFH3m[i] = exp(-(1.0/4.0)*Q2effH3*pow(gamma,2.0))*fabs(sumH3m);

      xFHe3m[i] = Q2He3;
      yFHe3m[i] = exp(-(1.0/4.0)*Q2effHe3*pow(gamma,2.0))*fabs(sumHe3m);
      
      sumH3m = 0.0;      //Reset sumH3ch.
      sumHe3m = 0.0;     //Reset sumHe3ch.

      //cout<<"Q^2 = "<<xFH3m[i]<<"   FH3m = "<<yFH3m[i]<<endl;

      //Calculation of H3 and He3 cross sections.
      //wH3 = (Q2H3*1.0/GeV2fm)/(2.0*MtH3);          //Put Q^2 back to GeV^2 to match mass then switch back to fm^-2 in q2_3 calculations.
      //wHe3 = (Q2H3*1.0/GeV2fm)/(2.0*MtHe3);
      wH3 = E0 - EfH3;                  
      wHe3 = E0 - EfHe3;
      q2_3H3 = fabs(  pow(wH3,2.0)*GeV2fm - Q2effH3  );         //Everything back in fm^-2.
      q2_3He3 = fabs(  pow(wHe3,2.0)*GeV2fm - Q2effHe3 );
      etaH3 = 1.0 + Q2effH3/(4.0*pow(MtH3,2.0)*GeV2fm);    //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.
      etaHe3 = 1.0 + Q2effHe3/(4.0*pow(MtHe3,2.0)*GeV2fm); 

      mottcrosssectionH3[i] = ( (pow(1,2)*(EfH3/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((angle)/2.0),4.0)))*pow(cos((angle)/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
      mottcrosssectionHe3[i] = ( (pow(2,2)*(EfHe3/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((angle)/2.0),4.0)))*pow(cos((angle)/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

      xcrosssectionH3[i] = Q2H3;
      ycrosssectionH3[i] = fabs( mottcrosssectionH3[i]*(1.0/etaH3) * (  (Q2effH3/q2_3H3)*pow(yFH3ch[i],2.0) + (pow(muH3,2.0)*Q2effH3/(2.0*pow(MtH3,2.0))) * (Q2effH3/(2.0*q2_3H3) + pow(tan(angle/2.0),2.0)) * pow(yFH3m[i],2.0)  ) );
      xcrosssectionHe3[i] = Q2He3;
      ycrosssectionHe3[i] = fabs( mottcrosssectionHe3[i]*(1.0/etaHe3) * (  (Q2effHe3/q2_3He3)*pow(yFHe3ch[i],2.0) + (pow(muHe3,2.0)*Q2effHe3/(2.0*pow(MtHe3,2.0))) * (Q2effHe3/(2.0*q2_3He3) + pow(tan(angle/2.0),2.0)) * pow(yFHe3m[i],2.0)  ) );


      //cout<<"ycrosssectionH3 = "<<ycrosssectionH3[i]<<"   ycrosssectionHe3 = "<<ycrosssectionHe3[i]<<endl;
      //cout<<"angle = "<<angle<<"   dquad = "<<dquad<<"   hquad = "<<hquad<<endl;
      //cout<<"Mott Cross Section = "<<mottxsection<<endl;
      //cout<<"Q^2 = "<<Q2<<"   protondipole = "<<protondipole<<endl;

      x[i] = dquad;
      y[i] = mottxsection*hquad*0.008;

      xdipole[i] = dquad;
      ydipole[i] = mottxsection*hquad*protondipole*0.008;

      xH3crosssection[i] = dquad;
      yH3crosssection[i] = mottcrosssectionH3[i]*hquad*ycrosssectionH3[i];

      xHe3crosssection[i] = dquad;
      yHe3crosssection[i] = mottcrosssectionHe3[i]*hquad*ycrosssectionHe3[i];

      //cout<<"x["<<i<<"] = "<<x[i]<<"   y["<<i<<"] = "<<y[i]<<endl;
      //dquad = dquad + 0.001;
      //cout<<"Q^2 = "<<Q2<<"   dquad = "<<dquad<<"   yH3crosssection = "<<yH3crosssection[i]<<"   yHe3crosssection = "<<yHe3crosssection[i]<<endl;
    }

  for(i=0; i<nmottpoints; i++)
    {
      //Starting on the left side of the sieve aperture track across to the right.
      step = i;    //Convert i to double. 
      angle = lhrsangle*deg2rad-0.5*quadangle+(step/nmottpoints)*quadangle;
      dquad = (angle-lhrsangle*deg2rad)*(0.06653/0.03)*-1.0; //Multiply by -1 to make higher rates on right side of sieve plate from beam's prespective. 
      hquad = (-0.008/0.043)*(dquad)+0.01;
      mottxsection = (  (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin((angle)/2.0),4.0)))*pow(cos((angle)/2.0),2.0)  ) * 1.0/25.7;   //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
      //Q2 = 4.0*pow(E0,2.0)*pow(sin(angle/2.0),2.0) * GeV2fm;             //fm^-2
      //protondipole = pow((1+(Q2*1.0/GeV2fm)/0.71),-2.0);                 //Unitless.


      Efproton = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/Mtproton);
      Q2proton = 4.0*E0*Efproton*pow(sin(angle/2.0),2.0) * GeV2fm;         //[fm^-2]
      Q2effproton = pow( pow(Q2proton,0.5) * (1.0+(1.5*1*(1/137.0))/(E0*pow(GeV2fm,0.5)*1.12*pow(1.0,1.0/3.0))) ,2.0);
      protondipole = pow((1+(Q2proton*1.0/GeV2fm)/0.71),-2.0);                  //Unitless.



      //cout<<"angle = "<<angle<<"   dquad = "<<dquad<<"   hquad = "<<hquad<<endl;
      //cout<<"Mott Cross Section = "<<mottxsection<<endl;

      xmott[i] = dquad;
      ymott[i] = mottxsection*0.00008;

      xmottdipole[i] = dquad;
      ymottdipole[i] = mottxsection*protondipole*0.00008;

      //cout<<"ymott = "<<ymott[i]<<"   ymottdipole = "<<ymottdipole[i]<<endl;
      //cout<<"x["<<i<<"] = "<<x[i]<<"   y["<<i<<"] = "<<y[i]<<endl;
      //dquad = dquad + 0.001;
    }

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph = new TGraph(npoints,x,y);
  //Draw the new TGraph called graph on the canvas. 
  graph->Draw();
  
  //Set X axis
  graph->GetXaxis()->SetLimits(-0.022,0.022);
  //Set Y axis Min and Max (not sure why different from X).
  graph->SetMinimum(0);
  graph->SetMaximum(0.60*pow(10.0,-6.0));
  graph->SetLineWidth(3);
  graph->SetLineColor(4);
  graph->SetFillColor(0);
  //gPad->SetLogx(1);
  //gPad->SetLogy(1);
  graph->SetTitle("Sieve Aperture Rates; Distance from Horizontal Center of Sieve Aperture [m]; Cross Section [fm^2]");
  //graph_expected.SetFillColor(kYellow);
  //graph_expected.DrawClone("E3AL"); // E3 draws the band

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph1 = new TGraph(nmottpoints,xmott,ymott);
  //Draw the new TGraph called graph on the canvas. 
  graph1->Draw("same");
  graph1->SetLineWidth(3);
  graph1->SetLineColor(2);
  graph1->SetFillColor(0);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph2 = new TGraph(npoints,xdipole,ydipole);
  //Draw the new TGraph called graph on the canvas. 
  graph2->Draw("same");
  graph2->SetLineWidth(3);
  graph2->SetLineColor(1);
  graph2->SetFillColor(0);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph3 = new TGraph(nmottpoints,xmottdipole,ymottdipole);
  //Draw the new TGraph called graph on the canvas. 
  graph3->Draw("same");
  graph3->SetLineWidth(3);
  graph3->SetLineColor(3);
  graph3->SetFillColor(0);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph12 = new TGraph(npoints,xH3crosssection,yH3crosssection);
  //Draw the new TGraph called graph on the canvas. 
  graph12->Draw("same");
  graph12->SetLineWidth(3);
  graph12->SetLineColor(7);
  graph12->SetFillColor(0);

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph13 = new TGraph(npoints,xHe3crosssection,yHe3crosssection);
  //Draw the new TGraph called graph on the canvas. 
  graph13->Draw("same");
  graph13->SetLineWidth(3);
  graph13->SetLineColor(6);
  graph13->SetFillColor(0);

  // Draw the Legend
  TLegend leg(0.9,.7,.56,.9,"");
  leg.SetFillColor(0);
  leg.AddEntry(graph,"Mott Cross Section * Height Sieve Aperture * 8*10^-3");
  leg.AddEntry(graph1,"Mott Cross Section * 8*10^-5");
  leg.AddEntry(graph2,"Mott Cross Section * Proton Dipole * Height Sieve Aperture * 8*10^-3");
  leg.AddEntry(graph3,"Mott Cross Section * Proton Dipole * 8*10^-5");
  leg.AddEntry(graph12,"Mott Cross Section * H3 Cross Section * Height Sieve Aperture");
  leg.AddEntry(graph13,"Mott Cross Section * He3 Cross Section * Height Sieve Aperture");
  leg.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

}

