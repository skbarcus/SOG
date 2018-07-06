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

void Calculate_XS() 
{
  Double_t pi = 3.141592654;
  Double_t deg2rad = pi/180.0;
  Double_t GeV2fm = 1.0/0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
  Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
  Double_t C = 299792458.0;                //Speed of light [m/s]. 
  Double_t angle = 21.0*deg2rad;                   //Scattering angle [rad].

  Double_t alpha = 1.0/137.0;              //Fine structure constant.
  Double_t E0 = 3.356;                     //Initial electron energy [GeV].
  Double_t Ef = 0.0;                       //Final energy of the electron after scattering.
  Double_t EfH3 = 0.0;
  Double_t EfHe3 = 0.0;
  Double_t Efproton = 0.0;
  Double_t mottxsection = 0.0;             //Mott cross section [1/GeV^2].

  Double_t Q2 = 0.0;                       //Q^2 in GeV.
  Double_t Q2H3 = 0.0;
  Double_t Q2He3 = 0.0;
  Double_t Q2effH3 = 0.0;                  //Q^2 with Coulomb corrections.
  Double_t Q2effHe3 = 0.0;

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
  Double_t Fch3H, Fch3He, Fm3H, Fm3He;    //Form factors.
  Double_t mottxs3H, mottxs3He;           //Mott XS.
  Double_t xs3H, xs3He;                   //Absolute XS.
  Double_t xsch3He, xsm3He;

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

  //Calculate charge FFs.
  EfH3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtH3);
  EfHe3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtHe3);
  Q2H3 = 4.0*E0*EfH3*pow(sin(angle/2.0),2.0) * GeV2fm;
  Q2He3 = 4.0*E0*EfHe3*pow(sin(angle/2.0),2.0) * GeV2fm;
  Q2effH3 = pow( pow(Q2H3,0.5) * (1.0+(1.5*1*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);   //Z=1
  Q2effHe3 = pow( pow(Q2He3,0.5) * (1.0+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);  //Z=2
  
  cout<<"3He: Ef = "<<EfHe3<<" GeV,   Q^2 3He = "<<Q2He3<<" fm^-2,   Qeff^2 = "<<Q2effHe3<<" fm^-2"<<endl;

  //Calculate the sum part of the SOG paramaterization for H3 charge FF.
  for(Int_t j=0; j<12; j++)
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
  for(Int_t j=0; j<12; j++)
    {
      sumHe3chtemp = (QHe3ch[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effHe3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effHe3,0.5)*R[j])/(pow(Q2effHe3,0.5)*R[j])) );
      sumHe3ch = sumHe3ch + sumHe3chtemp;
      //cout<<"sumH3ch = "<<sumH3ch<<endl;
    }
  
  //xFH3ch[i] = Q2H3;
  Fch3H = exp(-(1.0/4.0)*Q2effH3*pow(gamma,2.0))*fabs(sumH3ch);
  
  //xFHe3ch[i] = Q2He3;
  Fch3He = exp(-(1.0/4.0)*Q2effHe3*pow(gamma,2.0))*fabs(sumHe3ch);
  
  cout<<"Fch3He = "<<Fch3He<<endl;

  sumH3ch = 0.0;      //Reset sumH3ch.
  sumHe3ch = 0.0;     //Reset sumHe3ch.
  
  
  //Calculate magnetic FFs.
  EfH3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtH3);
  EfHe3 = E0/(1.0+2.0*E0*pow(sin(angle/2.0),2.0)/MtHe3);
  Q2H3 = 4.0*E0*EfH3*pow(sin(angle/2.0),2.0) * GeV2fm;
  Q2He3 = 4.0*E0*EfHe3*pow(sin(angle/2.0),2.0) * GeV2fm;
  Q2effH3 = pow( pow(Q2H3,0.5) * (1.0+(1.5*1*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);   //Z=1
  Q2effHe3 = pow( pow(Q2He3,0.5) * (1.0+(1.5*2*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(3.0,1.0/3.0))) ,2.0);  //Z=2
  
  //Calculate the sum part of the SOG paramaterization for H3 magnetic FF.
  for(Int_t j=0; j<12; j++)
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
  for(Int_t j=0; j<12; j++)
    {
      sumHe3mtemp = (QHe3m[j]/(1.0+2.0*pow(R[j],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2effHe3,0.5)*R[j]) + (2.0*pow(R[j],2.0)/pow(gamma,2.0)) * (sin(pow(Q2effHe3,0.5)*R[j])/(pow(Q2effHe3,0.5)*R[j])) );
      sumHe3m = sumHe3m + sumHe3mtemp;
      //cout<<"sumH3m = "<<sumH3m<<endl;
    }
  
  //Calculate FFs.
  //xFH3m[i] = Q2H3;
  Fm3H = exp(-(1.0/4.0)*Q2effH3*pow(gamma,2.0))*fabs(sumH3m);
  
  //xFHe3m[i] = Q2He3;
  Fm3He = exp(-(1.0/4.0)*Q2effHe3*pow(gamma,2.0))*fabs(sumHe3m);

  cout<<"Fm3He = "<<Fm3He<<endl;
  
  sumH3m = 0.0;      //Reset sumH3ch.
  sumHe3m = 0.0;     //Reset sumHe3ch.
  
  //Calculation of H3 and He3 cross sections.
  wH3 = E0 - EfH3;                  
  wHe3 = E0 - EfHe3;
  //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
  q2_3H3 = fabs(  pow(wH3,2.0)*GeV2fm - Q2effH3  );                  //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
  q2_3He3 = fabs(  pow(wHe3,2.0)*GeV2fm - Q2effHe3  );
  etaH3 = 1.0 + Q2effH3/(4.0*pow(MtH3,2.0)*GeV2fm);        //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.
  etaHe3 = 1.0 + Q2effHe3/(4.0*pow(MtHe3,2.0)*GeV2fm); 
  
  //Calculate Mott cross section and remember to multiply by Z^2 and the recoil factor Ef/Ei. 
  mottxs3H = (  (pow(1,2)*(EfH3/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle/2.0),4.0)))*pow(cos(angle/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  mottxs3He = (  (pow(2,2)*(EfHe3/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(angle/2.0),4.0)))*pow(cos(angle/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  cout<<"3He: Mott XS = "<<mottxs3He<<" fm^2"<<endl;
  
  //xcrosssectionH3[i] = Q2H3;
  xs3H = fabs( mottxs3H*(1.0/etaH3) * (  (Q2effH3/q2_3H3)*pow(Fch3H,2.0) + (pow(muH3,2.0)*Q2effH3/(2.0*pow(MtH3,2.0)*GeV2fm)) * (Q2effH3/(2.0*q2_3H3) + pow(tan(angle/2.0),2.0)) * pow(Fm3H,2.0)  ) );
  //xcrosssectionHe3[i] = Q2He3;
  xs3He = fabs( mottxs3He*(1.0/etaHe3) * (  (Q2effHe3/q2_3He3)*pow(Fch3He,2.0) + (pow(muHe3,2.0)*Q2effHe3/(2.0*pow(MtHe3,2.0)*GeV2fm)) * (Q2effHe3/(2.0*q2_3He3) + pow(tan(angle/2.0),2.0)) * pow(Fm3He,2.0)  ) );

  //Calculate the charge and magnetic contributions to the XS.
  xsch3He = fabs( mottxs3He*(1.0/etaHe3) * (  (Q2effHe3/q2_3He3)*pow(Fch3He,2.0)  ) );
  xsm3He = fabs( mottxs3He*(1.0/etaHe3) * (  (pow(muHe3,2.0)*Q2effHe3/(2.0*pow(MtHe3,2.0)*GeV2fm)) * (Q2effHe3/(2.0*q2_3He3) + pow(tan(angle/2.0),2.0)) * pow(Fm3He,2.0)  ) );

		 cout<<"3He: Charge contribution to XS = "<<xsch3He<<" fm^2/sr.   Magnetic contribution to XS = "<<xsm3He<<" fm^2/sr. Percent magnetic contribution = "<<100.*xsm3He/(xsch3He+xsm3He)<<"%."<<endl;
		 cout<<"3He: XS = "<<xs3He<<" fm^2/sr = "<<xs3He*pow(10,4)<<" ub/sr"<<endl;
}
