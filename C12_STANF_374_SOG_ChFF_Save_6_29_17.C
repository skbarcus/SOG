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

Double_t pi = 3.141592654;
Double_t deg2rad = pi/180.0;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t alpha = 1.0/137.0;              //Fine structure constant.

Int_t userand = 0;
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 0;
//Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = No Fit use preset R[i] and Qi.
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 11;                        //Number of Gaussians used to fit data.
Double_t Z = 6.;                         //Atomic number C12.
Double_t A = 12.;                        //Mass number C12.
Double_t MtC = 12.0*0.9315;              //Mass of C12 in GeV.
Double_t gamma = 0.935959;//0.8/pow(3./2.,0.5);//1.2*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
Double_t E0 = 0.3745;                    //Initial e- energy GeV.
Double_t Ef = 0.;                        //Final e- energy GeV.
Double_t ymin = 30.;//30
Double_t ymax = 100.;//100
Double_t yminFF = 0.;//30
Double_t ymaxFF = 6.;
Double_t range = fabs(ymax - ymin);
Int_t n = 10000;
Int_t ndim = n+1;
Double_t truncate = 25.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 2.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char* str[1000];                          //Variable to read lines of the data file.
Float_t thetatemp,qefftemp,sigexptemp,uncertaintytemp;
Float_t theta[100];                     //Angle in degrees.
Float_t qeff[100];                      //q effective in fm^-1.
Float_t sigexp[100];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[100];

Double_t m = 2.;
//Double_t R[12] = {0.1*m, 0.5*m, 0.9*m, 1.3*m, 1.6*m, 2.0*m, 2.4*m, 2.9*m, 3.4*m, 4.0*m, 4.6*m, 5.2*m};  //Radii [fm].
Double_t R[11] = {0.00000001,0.4,1.,1.3,1.7,2.3,2.7,3.5,4.3,5.4,6.7};//,7.9,9.1};
Double_t Qi[11] = {0.01669, 0.050325, 0.128621, 0.180515, 0.219097, 0.0278416, 0.058779, 0.057817, 0.007739, 0.002001, 0.000007};

void C12_STANF_374_SOG() 
{
  //Make a new canvas to plot data.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  //Read in data from text file.
  //Open file.

  if(usedifmin == 0)
    {
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/STANF_374_No_Min.txt","r");
    }

  if(usedifmin == 1)
    {
      FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/STANF_374.txt","r");
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
	ncols = fscanf(fp,"%f %f %f %f",&thetatemp, &qefftemp, &sigexptemp, &uncertaintytemp);
	if (ncols < 0) break;    
	//cout<<thetatemp<<"   "<<qefftemp<<"   "<<sigexptemp<<"   "<<uncertaintytemp<<endl;
	theta[nlines-skip] = thetatemp;
	qeff[nlines-skip] = qefftemp;
	sigexp[nlines-skip] = sigexptemp;
        uncertainty[nlines-skip] = uncertaintytemp;
	//Fill histograms with x and y data.
	//h1->Fill(x);
	//h2->Fill(x,y);
	//h2->SetMarkerSize(5);
	//Fill ntuple with x and y data.
	//ntuple->Fill(x,y);
	//Count the number of lines in the file. 
	nlines++;
      }
  }

//Print the data read from the file. 
  for(int i=0; i<37; i++)
    {
      cout<<"theta["<<i<<"] = "<<theta[i]<<"   qeff["<<i<<"] = "<<qeff[i]<<"   sigexp["<<i<<"] = "<<sigexp[i]<<"   uncertainty["<<i<<"] = "<<uncertainty[i]<<endl;
    }

  cout<<"Number of lines = "<<nlines<<endl;
  fclose(fp);

 //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  TGraphErrors *graph = new TGraphErrors(nlines-skip,theta,sigexp,0,uncertainty);
  //Draw the new TGraph called graph on the canvas. 
  graph->Draw("");
  c1->SetLogy();
  //Set X axis
  //graph->GetXaxis()->SetLimits(-12,12);
  //Set Y axis Min and Max (not sure why different from X).
  //graph->SetMinimum(0);
  //graph->SetMaximum(120);
  graph->SetLineWidth(1);
  graph->SetLineColor(4);
  graph->SetFillColor(0);
  graph->SetMarkerColor(1);
  graph->SetMarkerSize(0.4);
  graph->SetMarkerStyle(20);
  graph->SetTitle("C12 Cross Section; Angle (degrees); #\sigma_{exp}");
  //graph_expected.SetFillColor(kYellow);
  //graph_expected.DrawClone("E3AL"); // E3 draws the band

  // Draw the Legend
  TLegend leg(0.9,.7,.56,.9,"Legend Title");
  leg.SetFillColor(0);
  leg.AddEntry(graph,"Curve Name");
  leg.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.



  /*
  //Define fit function for H3 charge FF as a function of Q^2.
  Double_t fitch(Double_t *angle, Double_t *par)
  {
    Double_t fitval = 0.;
    Double_t sumH3chtemp = 0.;
    Ef = E0/(1.0+2.0*E0*pow(sin(angle[0]*deg2rad/2.0),2.0)/MtC);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle[0]*deg2rad/2.0),2.0) * GeV2fm;
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*6.*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(12.0,1.0/3.0))) ,2.0);   //A=6 Z=12
    
    for(Int_t i=0; i<ngaus; i++)
      {
	//sumH3chtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Fit R[i] values and Qi.
	sumH3chtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );

	fitval = fitval + sumH3chtemp;
      }
    
    fitval = fabs( fitval ) * exp(-0.25*Q2eff*pow(gamma,2.0));
    return fitval;
  }
  */
  
  //TF1 *H3chFF = new TF1("H3chFF",fitch, ymin, ymax,npar);
  //H3chFF->SetParLimits(0,0.00001,100.);
  //c1->SetLogy();
  
  /*
  H3chFF->SetParameter(0,1.);
  H3chFF->SetParameter(1,1.);
  H3chFF->SetParameter(2,1.);
  H3chFF->SetParameter(3,1.);
  H3chFF->SetParameter(4,1.);
  H3chFF->SetParameter(5,1.);
  H3chFF->SetParameter(6,1.);
  H3chFF->SetParameter(7,1.);
  H3chFF->SetParameter(8,1.);
  H3chFF->SetParameter(9,1.);
  H3chFF->SetParameter(10,1.);
  H3chFF->SetParameter(11,1.);
*/
  /*
  H3chFF->SetParameter(12,R[0]);
  H3chFF->SetParameter(13,R[1]);
  H3chFF->SetParameter(14,R[2]);
  H3chFF->SetParameter(15,R[3]);
  H3chFF->SetParameter(16,R[4]);
  H3chFF->SetParameter(17,R[5]);
  H3chFF->SetParameter(18,R[6]);
  H3chFF->SetParameter(19,R[7]);
  H3chFF->SetParameter(20,R[8]);
  H3chFF->SetParameter(21,R[9]);
  H3chFF->SetParameter(22,R[10]);
  H3chFF->SetParameter(23,R[11]);
  */
  /*
  H3chFF->SetParLimits(0,0.000001,10000);
  H3chFF->SetParLimits(1,0.000001,10000);
  H3chFF->SetParLimits(2,0.000001,10000);
  H3chFF->SetParLimits(3,0.000001,10000);
  H3chFF->SetParLimits(4,0.000001,10000);
  H3chFF->SetParLimits(5,0.000001,10000);
  H3chFF->SetParLimits(6,0.000001,10000);
  H3chFF->SetParLimits(7,0.000001,10000);
  H3chFF->SetParLimits(8,0.000001,10000);
  H3chFF->SetParLimits(9,0.000001,10000);
  H3chFF->SetParLimits(10,0.000001,10000);
  H3chFF->SetParLimits(11,0.000001,10000);

  H3chFF->SetParLimits(12,0.000001,10000);
  H3chFF->SetParLimits(13,0.000001,10000);
  H3chFF->SetParLimits(14,0.000001,10000);
  H3chFF->SetParLimits(15,0.000001,10000);
  H3chFF->SetParLimits(16,0.000001,10000);
  H3chFF->SetParLimits(17,0.000001,10000);
  H3chFF->SetParLimits(18,0.000001,10000);
  H3chFF->SetParLimits(19,0.000001,10000);
  H3chFF->SetParLimits(20,0.000001,10000);
  H3chFF->SetParLimits(21,0.000001,10000);
  H3chFF->SetParLimits(22,0.000001,10000);
  H3chFF->SetParLimits(23,0.000001,10000);
*/
  /*
  graph->Fit(H3chFF,"0");
  graph->Draw("same");
  H3chFF->SetLineColor(3);
  H3chFF->Draw("same");
  cout<<H3chFF->Eval(42.)<<endl;

  
 //Define fit function for H3 magnetic FF as a function of Q^2.
  Double_t fitm(Double_t *Q2, Double_t *par)
  {
    Double_t fitval = 0.;
    Double_t sumH3mtemp = 0.;
    
    for(Int_t i=0; i<ngaus; i++)
      {
	//sumH3chtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	//Fit R[i] values and Qi.
	sumH3mtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*par[ngaus+i])/(pow(Q2[0],0.5)*par[ngaus+i])) );

	fitval = fitval + sumH3mtemp;
      }
    
    fitval = fabs( fitval ) * exp(-0.25*Q2[0]*pow(gamma,2.0));
    return fitval;
  }
 */ 
  
  //TF1 *H3mFF = new TF1("H3mFF",fitm, ymin, ymax,npar);
  //H3chFF->SetParLimits(0,0.00001,100.);
  //c1->SetLogy();
  
  /*
  H3mFF->SetParameter(0,1.);
  H3mFF->SetParameter(1,1.);
  H3mFF->SetParameter(2,1.);
  H3mFF->SetParameter(3,1.);
  H3mFF->SetParameter(4,1.);
  H3mFF->SetParameter(5,1.);
  H3mFF->SetParameter(6,1.);
  H3mFF->SetParameter(7,1.);
  H3mFF->SetParameter(8,1.);
  H3mFF->SetParameter(9,1.);
  H3mFF->SetParameter(10,1.);
  H3mFF->SetParameter(11,1.);
*/
  /*
  H3mFF->SetParameter(12,R[0]);
  H3mFF->SetParameter(13,R[1]);
  H3mFF->SetParameter(14,R[2]);
  H3mFF->SetParameter(15,R[3]);
  H3mFF->SetParameter(16,R[4]);
  H3mFF->SetParameter(17,R[5]);
  H3mFF->SetParameter(18,R[6]);
  H3mFF->SetParameter(19,R[7]);
  H3mFF->SetParameter(20,R[8]);
  H3mFF->SetParameter(21,R[9]);
  H3mFF->SetParameter(22,R[10]);
  H3mFF->SetParameter(23,R[11]);
  */
  /*
  H3mFF->SetParLimits(0,0.000001,10000);
  H3mFF->SetParLimits(1,0.000001,10000);
  H3mFF->SetParLimits(2,0.000001,10000);
  H3mFF->SetParLimits(3,0.000001,10000);
  H3mFF->SetParLimits(4,0.000001,10000);
  H3mFF->SetParLimits(5,0.000001,10000);
  H3mFF->SetParLimits(6,0.000001,10000);
  H3mFF->SetParLimits(7,0.000001,10000);
  H3mFF->SetParLimits(8,0.000001,10000);
  H3mFF->SetParLimits(9,0.000001,10000);
  H3mFF->SetParLimits(10,0.000001,10000);
  H3mFF->SetParLimits(11,0.000001,10000);

  H3mFF->SetParLimits(12,0.000001,10000);
  H3mFF->SetParLimits(13,0.000001,10000);
  H3mFF->SetParLimits(14,0.000001,10000);
  H3mFF->SetParLimits(15,0.000001,10000);
  H3mFF->SetParLimits(16,0.000001,10000);
  H3mFF->SetParLimits(17,0.000001,10000);
  H3mFF->SetParLimits(18,0.000001,10000);
  H3mFF->SetParLimits(19,0.000001,10000);
  H3mFF->SetParLimits(20,0.000001,10000);
  H3mFF->SetParLimits(21,0.000001,10000);
  H3mFF->SetParLimits(22,0.000001,10000);
  H3mFF->SetParLimits(23,0.000001,10000);
*/

  /*
  graph->Fit(H3mFF,"0");
  //graph->Draw("same");
  H3mFF->SetLineColor(4);
  H3mFF->Draw("same");
  cout<<H3mFF->Eval(42.)<<endl;
*/



  /*
  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  
Double_t mottxs(Double_t *angle2, Double_t *par)
  {
    Double_t val = 0.; 
    Ef = E0/(1.0+2.0*E0*pow(sin(angle2[0]*deg2rad/2.0),2.0)/MtC);

    val = (  (pow(6.,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle2[0]*deg2rad/2.0),4.0)))*pow(cos(angle2[0]*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.

    return val;
  }

  TF1 *fmottxs = new TF1("fmottxs",mottxs, ymin, ymax,1);
  cout<<"!!! = "<<fmottxs->Eval(40.)<<endl;
  //fmottxs->Draw();
*/

 Double_t xs(Double_t *angle, Double_t *par)
  {
    Double_t val = 0.;
    Double_t mottxs = 0.;
    Double_t fitch = 0.;
    Double_t sumchtemp = 0.;
    Double_t fitm = 0.;
    Double_t summtemp = 0.;

    Ef = E0/(1.0+2.0*E0*pow(sin(angle[0]*deg2rad/2.0),2.0)/MtC);
    Double_t Q2 = 4.0*E0*Ef*pow(sin(angle[0]*deg2rad/2.0),2.0) * GeV2fm;
    Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=6 A=12
                
    Double_t W = E0 - Ef;
    //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtC);
    Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    Double_t eta = 1.0 + Q2eff/(4.0*pow(MtC,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.

    //Calculate Mott XS.
    mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle[0]*deg2rad/2.0),4.0)))*pow(cos(angle[0]*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    

    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 
	//Fit just the Qi values using predetermined R[i] values.
	//sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
       
	//Fit R[i] values and Qi.
	//sumchtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );

	//Fit R[i] values, Qi coefficients, and gamma.
	sumchtemp = (par[i]/(1.0+2.0*pow(par[ngaus+i],2.0)/pow(par[2*ngaus],2.0))) * ( cos(pow(Q2eff,0.5)*par[ngaus+i]) + (2.0*pow(par[ngaus+i],2.0)/pow(par[2*ngaus],2.0)) * (sin(pow(Q2eff,0.5)*par[ngaus+i])/(pow(Q2eff,0.5)*par[ngaus+i])) );
	
	
	//Use 'known' C12 coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
	
       fitch =  fitch + sumchtemp;
      }
    
    //fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
    fitch =  fitch * exp(-0.25*Q2eff*pow(par[2*ngaus],2.0));
   
    
    //Define SOG for magnetic FF.
    for(Int_t i=0; i<ngaus; i++)
      {
	//summtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
	//Fit R[i] values and Qi.
	summtemp = (par[2*ngaus +i]/(1.0+2.0*pow(par[3*ngaus+i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2eff,0.5)*par[3*ngaus+i]) + (2.0*pow(par[3*ngaus+i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2eff,0.5)*par[3*ngaus+i])/(pow(Q2eff,0.5)*par[3*ngaus+i])) );

	fitm = fitm + summtemp;
      }
    
    fitm = fabs( fitm ) * exp(-0.25*Q2eff*pow(gamma,2.0));

    val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) ); //magnetic moment for C12 is 0 -> no mag part of XS.

    return val;
  }
 
 //c2->SetLogy(); 
 graph->Draw();
 TF1 *fxs = new TF1("fxs",xs, ymin, ymax,ngaus*2+1);//ngaus*2);

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

 for(Int_t i=0;i<ngaus;i++)
   {
     cout<<"R["<<i<<"] = "<<R[i]<<endl;
   }
 
 /*
 fxs->SetParameter(0,1.);
 fxs->SetParameter(1,1.);
 fxs->SetParameter(2,1.);
 fxs->SetParameter(3,1.);
 fxs->SetParameter(4,1.);
 fxs->SetParameter(5,1.);
 fxs->SetParameter(6,1.);
 fxs->SetParameter(7,1.);
 fxs->SetParameter(8,1.);
 fxs->SetParameter(9,1.);
 fxs->SetParameter(10,1.);
 fxs->SetParameter(11,1.);
*/
 /*
 fxs->SetParameter(0,0.001);
 fxs->SetParameter(1,0.001);
 fxs->SetParameter(2,0.001);
 fxs->SetParameter(3,0.001);
 fxs->SetParameter(4,0.001);
 fxs->SetParameter(5,0.001);
 fxs->SetParameter(6,0.001);
 fxs->SetParameter(7,0.001);
 fxs->SetParameter(8,0.001);
 fxs->SetParameter(9,0.001);
 fxs->SetParameter(10,0.001);
 */
 //fxs->SetParameter(11,0.001);
 //fxs->SetParameter(12,0.001);

 fxs->SetParameter(0,0.01669);
 fxs->SetParameter(1,0.050325);
 fxs->SetParameter(2,0.128621);
 fxs->SetParameter(3,0.180515);
 fxs->SetParameter(4,0.219097);
 fxs->SetParameter(5,0.0278416);

 fxs->SetParameter(6,0.00000001);
 fxs->SetParameter(7,0.00000001);
 fxs->SetParameter(8,0.00000001);
 fxs->SetParameter(9,0.00000001);
 fxs->SetParameter(10,0.00000001);
 /*
 fxs->SetParameter(6,0.058779);
 fxs->SetParameter(7,0.057817);
 fxs->SetParameter(8,0.007739);
 fxs->SetParameter(9,0.002001);
 fxs->SetParameter(10,0.000007);
 */
 /*
 fxs->SetParameter(12,R[0]);
 fxs->SetParameter(13,R[1]);
 fxs->SetParameter(14,R[2]);
 fxs->SetParameter(15,R[3]);
 fxs->SetParameter(16,R[4]);
 fxs->SetParameter(17,R[5]);
 fxs->SetParameter(18,R[6]);
 fxs->SetParameter(19,R[7]);
 fxs->SetParameter(20,R[8]);
 fxs->SetParameter(21,R[9]);
 fxs->SetParameter(22,R[10]);
 fxs->SetParameter(23,R[11]);
*/

 fxs->SetParameter(ngaus,0.00000001);
 fxs->SetParameter(ngaus+1,R[1]);
 fxs->SetParameter(ngaus+2,R[2]);
 fxs->SetParameter(ngaus+3,R[3]);
 fxs->SetParameter(ngaus+4,R[4]);
 fxs->SetParameter(ngaus+5,R[5]);
 fxs->SetParameter(ngaus+6,R[6]);
 fxs->SetParameter(ngaus+7,R[7]);
 fxs->SetParameter(ngaus+8,R[8]);
 fxs->SetParameter(ngaus+9,R[9]);
 fxs->SetParameter(ngaus+10,R[10]);

 fxs->SetParameter(2*ngaus,0.8);
 //fxs->SetParameter(ngaus+11,R[11]);
 //fxs->SetParameter(ngaus+12,R[12]);

 if(userand == 0)
   {
     Double_t d = 0.1001;//0.5001//0.1001
     
     fxs->SetParLimits(11,0.00000001,0.1001);//0.5
     fxs->SetParLimits(12,R[1]-d,R[1]+d);
     fxs->SetParLimits(13,R[2]-d,R[2]+d);
     fxs->SetParLimits(14,R[3]-d,R[3]+d);
     fxs->SetParLimits(15,R[4]-d,R[4]+d);
     fxs->SetParLimits(16,R[5]-d,R[5]+d);
     fxs->SetParLimits(17,R[6]-d,R[6]+d);
     fxs->SetParLimits(18,R[7]-d,R[7]+d);
     fxs->SetParLimits(19,R[8]-d,R[8]+d);
     fxs->SetParLimits(20,R[9]-d,R[9]+d);
     fxs->SetParLimits(21,R[10]-d,R[10]+d);
   } 

 if(userand == 1)
   {
     Double_t dx = 0.5001;//0.5001
     
     fxs->SetParLimits(11,0.00000001,0.001);//0.5
     fxs->SetParLimits(12,R[1]-dx,R[1]+dx);
     fxs->SetParLimits(13,R[2]-dx,R[2]+dx);
     fxs->SetParLimits(14,R[3]-dx,R[3]+dx);
     fxs->SetParLimits(15,R[4]-dx,R[4]+dx);
     fxs->SetParLimits(16,R[5]-dx,R[5]+dx);
     fxs->SetParLimits(17,R[6]-dx,R[6]+dx);
     fxs->SetParLimits(18,R[7]-dx,R[7]+dx);
     fxs->SetParLimits(19,R[8]-dx,R[8]+dx);
     fxs->SetParLimits(20,R[9]-dx,R[9]+dx);
     fxs->SetParLimits(21,R[10]-dx,R[10]+dx);
   } 

 // Set stat options
 //gStyle->SetStatY(0.9);                
 // Set y-position (fraction of pad size)
 //gStyle->SetStatX(0.9);                
 // Set x-position (fraction of pad size)
 gStyle->SetStatW(0.15);                
 // Set width of stat-box (fraction of pad size)
 gStyle->SetStatH(0.025);

 graph->Fit(fxs,"");
 gStyle->SetOptFit(100);
 fxs->Draw("same");
 
 Double_t chi2 = fxs->GetChisquare();
 cout<<"Chi^2 = "<<chi2<<endl;

 
 //Plot the charge FF using the parameters determined by the SOG fit above. 
 //Make a new canvas to plot data.
 TCanvas* c2=new TCanvas("c2");
 c2->SetGrid();
 c2->SetLogy();

 //Set Qi and R[i] = to the corresponding ffit parameters.
 for(Int_t i=0;i<ngaus;i++)
   {
     Qi[i] = fxs->GetParameter(i);
     cout<<"Qi["<<i<<"] = "<<Qi[i]<<endl;
   }

 for(Int_t i=0;i<ngaus;i++)
   {
     R[i] = fxs->GetParameter(i+ngaus);
     cout<<"R["<<i<<"] = "<<R[i]<<endl;
   }

 cout<<"Gamma = "<<fxs->GetParameter(2*ngaus)<<endl;

 Double_t ChFF(Double_t *Q, Double_t *par)
 {
   Double_t fitch = 0.;
   Double_t sumchtemp = 0.;

    //Define SOG for charge FF.
    for(Int_t i=0; i<ngaus; i++)
      { 	
	//Use SOG fit for C12 Qi coefficients and R[i] values. 
	//sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

	sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
	
	fitch = fitch + sumchtemp;
      }
    
    fitch = fitch * exp(-0.25*pow(Q[0],2.)*pow(gamma,2.0));
    fitch = fabs(fitch);
    return fitch;
 }

 TF1 *fChFF = new TF1("fChFF",ChFF,yminFF,ymaxFF,1);
 cout<<fChFF->Eval(0.000001)<<"!!!!!!!!"<<endl;
 fChFF->Draw("C");
 c2->SetTitle("C12 Charge Form Factor");
 //fChFF->SetTitle("C12 Charge Form Factor","#Q^2 (#fm^-2)","#F_{Ch}(q)");
 fChFF->GetHistogram()->GetYaxis()->SetTitle("|F_{Ch}(q)|");
 fChFF->GetHistogram()->GetXaxis()->SetTitle("q (fm^{-1})");



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
     TF1 *g0 = new TF1("g0", fitg, yminFF, ymaxFF,2.);
     g0->SetParameters(Qi[0],R[0]);
     g0->SetLineColor(1);
     g0->Draw("csame");
     TF1 *g1 = new TF1("g1", fitg, yminFF, ymaxFF,2.);
     g1->SetParameters(Qi[1],R[1]);
     g1->SetLineColor(2);
     g1->Draw("cSame");
     TF1 *g2 = new TF1("g2", fitg, yminFF, ymaxFF,2.);
     g2->SetParameters(Qi[2],R[2]);
     g2->SetLineColor(3);
     g2->Draw("cSame");
     TF1 *g3 = new TF1("g3", fitg, yminFF, ymaxFF,2.);
     g3->SetParameters(Qi[3],R[3]);
     g3->SetLineColor(4);
     g3->Draw("cSame");
     TF1 *g4 = new TF1("g4", fitg, yminFF, ymaxFF,2.);
     g4->SetParameters(Qi[4],R[4]);
     g4->SetLineColor(5);
     g4->Draw("cSame");
     TF1 *g5 = new TF1("g5", fitg, yminFF, ymaxFF,2.);
     g5->SetParameters(Qi[5],R[5]);
     g5->SetLineColor(6);
     g5->Draw("cSame");
     TF1 *g6 = new TF1("g6", fitg, yminFF, ymaxFF,2.);
     g6->SetParameters(Qi[6],R[6]);
     g6->SetLineColor(7);
     g6->Draw("cSame");
     TF1 *g7 = new TF1("g7", fitg, yminFF, ymaxFF,2.);
     g7->SetParameters(Qi[7],R[7]);
     g7->SetLineColor(8);
     g7->Draw("cSame");
     TF1 *g8 = new TF1("g8", fitg, yminFF, ymaxFF,2.);
     g8->SetParameters(Qi[8],R[8]);
     g8->SetLineColor(9);
     g8->Draw("cSame");
     TF1 *g9 = new TF1("g9", fitg, yminFF, ymaxFF,2.);
     g9->SetParameters(Qi[9],R[9]);
     g9->SetLineColor(10);
     g9->Draw("cSame");
     TF1 *g10 = new TF1("g10", fitg, yminFF, ymaxFF,2.);
     g10->SetParameters(Qi[10],R[10]);
     g10->SetLineColor(11);
     g10->Draw("cSame");
     TF1 *g11 = new TF1("g11", fitg, yminFF, ymaxFF,2.);
     g11->SetParameters(Qi[11],R[11]);
     g11->SetLineColor(12);
     g11->Draw("cSame");
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

  /*
  //Fill a new histogram with data using the fit function. This will then be inverse Fourier transformed to view the charge distribution.
  TH1 *hfit = new TH1D("hH3chFF", "hH3chFF", n+1, ymin, ymax);
  //Fill the histogram with function values
  for (Int_t i=0; i<=n; i++)
    {
      Double_t x = 0.;
      x = ymin + (Double_t(i)/(n))*(range+0.);//range;  //Only fill bins up to max of fitted range.
      //cout<<x<<endl;
      if(x==0.)
	{
	  x = 0.0000001; //Inf/NaN at zero.
	}
      
      if(x<=ymax && x<truncate)
	{
	  hH3chFF->SetBinContent(i+1, H3chFF->Eval(x));
	}
      
      if(x>ymax)
	{
	  //hfit->SetBinContent(i+1, (0.8*cos(12.*x-3.1)+0.5)*exp(-x*.1));
	  //hH3chFF->SetBinContent(i+1, 0.);
	}
      
      hH3chFF->GetEntries();
    }
  //hfit->SetFillColor(17);
  hH3chFF->Draw("same");
  
  //Inverse Fourier transform the H3 charge FF to get the charge distribution of H3.
  //Create arrays for real and complex imnputs for inverse FFT. 
  Double_t *re_full = new Double_t[n+1];
  Double_t *im_full = new Double_t[n+1];
  
  //Fill the real and complex arrays. The complex array is all zeros since we have real data. The real data is from the histo fit. 
  for(Int_t i=0;i<n+1;i++)
    {
      re_full[i] = hH3chFF->GetBinContent(i+1);
      im_full[i] = 0;
      //cout<<"re_full["<<i<<"] = "<<re_full[i]<<endl;
    }

  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();

  TVirtualFFT *iFFT = TVirtualFFT::FFT(1, &ndim, "C2R M K");   //Re and Mag look the same for C2R which makes sense. C2CBACKWARD mag definitely different from C2R and looks wrong. Same for C2CBackward re although first min x position looks ok. Stick to C2R mag it appears correct. 
  iFFT->SetPointsComplex(re_full,im_full);
  iFFT->Transform();
  TH1 *hcharge = 0;
  //TH1 *hcharge = new TH1D("hcharge", "hcharge", nback*10.+1, ymin, ymax);
  //Let's look at the output
  hcharge = TH1::TransformHisto(iFFT,hcharge,"mag");     //Not totally sure if we want re or mag.
  hcharge->SetTitle("H3 Charge Distribution");
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
  TH1 *H3ch_dist = new TH1D("H3ch_dist", "H3ch_dist", n+1, -(n+1)/(ymax*2.), (n+1)/(ymax*2.));   

  Int_t inflection = (n)/2.;
  //cout<<"inflection = "<<inflectionback<<endl;
  
  //Negative Fourier frequencies.
  for(Int_t i=0;i<n+2;i++)
    {
      if(i>inflection+1)
	{
	  H3ch_dist->SetBinContent(i-1-inflection,hcharge->GetBinContent(i)/(1./(range/(n+1.))));//range/(n+1.))));
	  //cout<<"i - inflectionback - 1 = "<<i<<" - "<<inflectionback<<"-1 = "<<i-inflectionback-1<<endl;
	}
    }
  
  //Positive Fourier frequencies.
  for(Int_t i=0;i<n;i++)
    {
      if(i<=inflection)
	{
	  H3ch_dist->SetBinContent(i+inflection+1,hcharge->GetBinContent(i+1)/(1./(range/(n+1.))));
	  //cout<<i+inflection+1<<endl;
	}
    }
  c2->SetLogy();
  H3ch_dist->SetTitle("H3 Charge Distribution");
  H3ch_dist->GetXaxis()->SetTitle("r");
  H3ch_dist->GetYaxis()->SetTitle("#rho(r)");
  H3ch_dist->Draw("");


  //Now (forward) Fourier transform the H3 charge distribution to recover the H3 charge FF. 

  //Make a new canvas to plot data.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  
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
  c3->SetLogy();
  hH3chFFback->SetTitle("H3 Charge Form Factor Transformed back (FF->iFFT(FF)->FFT(iFFT(FF))=FF)");
  hH3chFFback->GetXaxis()->SetTitle("Q^{2} [fm^{-2}]");
  hH3chFFback->GetYaxis()->SetTitle("Fch(Q)");
  hH3chFFback->Draw("");
  */
}
