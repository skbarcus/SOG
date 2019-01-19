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

Int_t target = 1;        //0 = 3He, 1 = 3H
Int_t ngaus = 8;
Double_t maxchi2 = 603;
const Int_t size = 3000;
Int_t skip = 1.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char* str[1000];                          //Variable to read lines of the data file.
Float_t Chi2[size],rChi2[size],BIC[size],AIC[size],Qichtot[size],Qimtot[size],R0[size],R1[size],R2[size],R3[size],R4[size],R5[size],R6[size],R7[size],R8[size],R9[size],R10[size],R11[size],R12[size],Q0ch[size],Q1ch[size],Q2ch[size],Q3ch[size],Q4ch[size],Q5ch[size],Q6ch[size],Q7ch[size],Q8ch[size],Q9ch[size],Q10ch[size],Q11ch[size],Q12ch[size],Q0m[size],Q1m[size],Q2m[size],Q3m[size],Q4m[size],Q5m[size],Q6m[size],Q7m[size],Q8m[size],Q9m[size],Q10m[size],Q11m[size],Q12m[size];
Float_t Chi2temp,rChi2temp,BICtemp,AICtemp,Qichtottemp,Qimtottemp,R0temp,R1temp,R2temp,R3temp,R4temp,R5temp,R6temp,R7temp,R8temp,R9temp,R10temp,R11temp,R12temp,Q0chtemp,Q1chtemp,Q2chtemp,Q3chtemp,Q4chtemp,Q5chtemp,Q6chtemp,Q7chtemp,Q8chtemp,Q9chtemp,Q10chtemp,Q11chtemp,Q12chtemp,Q0mtemp,Q1mtemp,Q2mtemp,Q3mtemp,Q4mtemp,Q5mtemp,Q6mtemp,Q7mtemp,Q8mtemp,Q9mtemp,Q10mtemp,Q11mtemp,Q12mtemp;
Int_t repeats_temp = 0;
Int_t Repeats[size];
Double_t Top_Ten[10][3];  //Stores the number of times a fit repeated followed by the Chi^2 of that fit followed by the first fit number with that Chi^2.
Int_t max_repeat = 0;  //Stores the largest number of repeat fits (same Ri and Chi^2).
Int_t replace = 0;

void Repeat_Fit_Analysis()
{
  FILE *fp;
  if(target == 0)
    {
      //3He
      fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=12_1352_12_22_2018.txt","r");//Final values.
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_Final_n=12_Short_12_22_2018.txt","r");
    }
  
  if(target == 1)
    {
      //3H
      fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=8_2600_12_22_2018.txt","r");//Final values.
      //fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Ri_Fits_3H_Final_n=8_Short_12_22_2018.txt","r");
    }

  //Read in data.
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
	  
	  //Reading in the normal datafile with max ngaus=12.
	  ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp, &R3temp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &Q0chtemp, &Q1chtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp, &Q11chtemp);
	  ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp, &Q8mtemp, &Q9mtemp);
	  ncols = ncols + fscanf(fp,"%f %f", &Q10mtemp, &Q11mtemp);
	  
	  //If using n=13 need to turn on reading the n=13 parameter columns.
	  /*
	    ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Chi2temp, &rChi2temp, &BICtemp, &AICtemp, &Qichtottemp, &Qimtottemp, &R0temp, &R1temp, &R2temp, &R3temp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &R4temp, &R5temp, &R6temp, &R7temp, &R8temp, &R9temp, &R10temp, &R11temp, &R12temp, &Q0chtemp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q1chtemp, &Q2chtemp, &Q3chtemp, &Q4chtemp, &Q5chtemp, &Q6chtemp, &Q7chtemp, &Q8chtemp, &Q9chtemp, &Q10chtemp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &Q11chtemp, &Q12chtemp, &Q0mtemp, &Q1mtemp, &Q2mtemp, &Q3mtemp, &Q4mtemp, &Q5mtemp, &Q6mtemp, &Q7mtemp);
	    ncols = ncols + fscanf(fp,"%f %f %f %f %f", &Q8mtemp, &Q9mtemp, &Q10mtemp, &Q11mtemp, &Q12mtemp);
	  */
	  
	  //cout<<"ncols = "<<ncols<<endl;
	  if (ncols < 0) break;    
	  
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
	  //R12[nlines-skip] = R12temp;
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
	  //Q12ch[nlines-skip] = Q12chtemp;
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
	  //Q12m[nlines-skip] = Q12mtemp;
	  /*
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
	    //Rmulti[nlines-skip][12] = R12temp;
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
	    //Qichmulti[nlines-skip][12] = Q12chtemp;
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
	    //Qimmulti[nlines-skip][12] = Q12mtemp;
	    */
	  //cout<<"!!! Chi2["<<nlines-skip<<"]"<<Chi2[nlines-skip]<<endl;

	  nlines++;
	}
    }

  //Print the data read from the file.
  for(int i=0; i<(nlines-skip); i++)
    {
      cout<<"Fit# = "<<i<<"  Chi2 = "<<Chi2[i]<<"  rChi2 = "<<rChi2[i]<<"  BIC = "<<BIC[i]<<"  AIC = "<<AIC[i]<<"  Qichtot = "<<Qichtot[i]<<"  Qimtot = "<<Qimtot[i]<<endl;
      cout<<endl;
      cout<<"R[0] = "<<R0[i]<<"  R[1] = "<<R1[i]<<"  R[2] = "<<R2[i]<<"  R[3] = "<<R3[i]<<"  R[4] = "<<R4[i]<<"  R[5] = "<<R5[i]<<"  R[6] = "<<R6[i]<<"  R[7] = "<<R7[i]<<"  R[8] = "<<R8[i]<<"  R[9] = "<<R9[i]<<"  R[10] = "<<R10[i]<<"  R[11] = "<<R11[i]<<endl;
      cout<<endl;
      cout<<"Qch[0] = "<<Q0ch[i]<<"  Qch[1] = "<<Q1ch[i]<<"  Qch[2] = "<<Q2ch[i]<<"  Qch[3] = "<<Q3ch[i]<<"  Qch[4] = "<<Q4ch[i]<<"  Qch[5] = "<<Q5ch[i]<<"  Qch[6] = "<<Q6ch[i]<<"  Qch[7] = "<<Q7ch[i]<<"  Qch[8] = "<<Q8ch[i]<<"  Qch[9] = "<<Q9ch[i]<<"  Qch[10] = "<<Q10ch[i]<<"  Qch[11] = "<<Q11ch[i]<<endl;
      cout<<endl;
      cout<<"Qm[0] = "<<Q0m[i]<<"  Qm[1] = "<<Q1m[i]<<"  Qm[2] = "<<Q2m[i]<<"  Qm[3] = "<<Q3m[i]<<"  Qm[4] = "<<Q4m[i]<<"  Qm[5] = "<<Q5m[i]<<"  Qm[6] = "<<Q6m[i]<<"  Qm[7] = "<<Q7m[i]<<"  Qm[8] = "<<Q8m[i]<<"  Qm[9] = "<<Q9m[i]<<"  Qm[10] = "<<Q10m[i]<<"  Qm[11] = "<<Q11m[i]<<endl;
      cout<<"-----------------------------------------------------"<<endl;
    }  
  cout<<"Number of lines = "<<nlines<<endl;

  fclose(fp);

  //Loop over the chi^2 values and check how many are equal to on another. 

  for(Int_t i=0;i<(nlines-1);i++)
    {
      repeats_temp = 0;//Reset repeat counter.
      for(Int_t j=0;j<(nlines-1);j++)
	{
	  if(Chi2[i]==Chi2[j] && j!=i)
	    {
	      repeats_temp++;
	    }
	}
      Repeats[i] = repeats_temp;
    }

  for(Int_t i=0;i<(nlines-1);i++)
    {
      if(Chi2[i]<maxchi2)
	{
	  //cout<<"Repeats["<<i<<"] = "<<Repeats[i]<<endl;

	  //Find largest value.
	  if(Repeats[i]>max_repeat)
	    {
	      for(Int_t j=9;j>0;j--)
		{
		  Top_Ten[j][0] = Top_Ten[j-1][0];
		  Top_Ten[j][1] = Top_Ten[j-1][1];
		  Top_Ten[j][2] = Top_Ten[j-1][2];
		}

	      Top_Ten[0][0] = Repeats[i];
	      Top_Ten[0][1] = Chi2[i];
	      Top_Ten[0][2] = i;
	      max_repeat = Repeats[i];
	    }
      
	  //Check if value is larger than any of the top 10. If so add it in descending order and move all other indices below it down by 1.
	  for(Int_t j=1;j<9;j++)
	    {
	      replace = 0;
	      if(Repeats[i]<max_repeat && Repeats[i]>Top_Ten[j][0] && Repeats[i]!=Top_Ten[1][0] && Repeats[i]!=Top_Ten[2][0] && Repeats[i]!=Top_Ten[3][0] && Repeats[i]!=Top_Ten[4][0] && Repeats[i]!=Top_Ten[5][0] && Repeats[i]!=Top_Ten[6][0] && Repeats[i]!=Top_Ten[7][0] && Repeats[i]!=Top_Ten[8][0] && Repeats[i]!=Top_Ten[9][0])//Messy but it works.
		{
		  for(Int_t k=9;k>j;k--)
		    {
		      Top_Ten[k][0] = Top_Ten[k-1][0];
		      Top_Ten[k][1] = Top_Ten[k-1][1];
		      Top_Ten[k][2] = Top_Ten[k-1][2];
		    }
		  replace = 1;
		  Top_Ten[j][0] = Repeats[i];
		  Top_Ten[j][1] = Chi2[i];
		  Top_Ten[j][2] = i;

		  //cout<<i<<endl;
		  for(Int_t z=0;z<10;z++)
		    {
		      //cout<<"Top_Ten["<<z<<"][0] = "<<Top_Ten[z][0]<<"   Top_Ten["<<z<<"][1] = "<<Top_Ten[z][1]<<endl;
		    }
		}
	      //Top_Ten[j][0] = Repeats[i];
	      //Top_Ten[j][1] = Chi2[i];
	      if(replace==1)
		{
		  break;
		}
	    }  
	} 
    }

  for(Int_t i=0;i<10;i++)
    {
      cout<<"Top_Ten["<<i<<"][0] = "<<Top_Ten[i][0]<<"   Top_Ten["<<i<<"][1] = "<<Top_Ten[i][1]<<"   Top_Ten["<<i<<"][2] = "<<Top_Ten[i][2]<<endl;
    }
}

