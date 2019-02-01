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

//Adding theory curves from Marcucci 2016. 4 theory predictions per target (2) per form factor (2) = 16 lines.
const Int_t size1 = 200;
Int_t skip1; 
Int_t nlines1;
char* str1[1000];
Float_t x[size1],y[size1];
Float_t x_temp,y_temp;

void Q_to_Q2() 
{
  //Create an output file to store a single fit result in.
  std::ofstream output ("3He_Fch_Conventional_Q2.txt", std::ofstream::out);

  //Fch conventional approach Marcucci 2016.
  FILE *fp1;

  fp1 = fopen("/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_Conventional.txt","r");
    
  //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines1 < skip1)
	{
	  fgets(str1,1000,fp1);
	  nlines1++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols1 = fscanf(fp1,"%f %f", &x_temp, &y_temp);
	  
	  if (ncols1 < 0) break;    
	  
	  x[nlines1-skip1] = x_temp * x_temp;
	  y[nlines1-skip1] = y_temp;
	  cout<<"x["<<nlines1-skip1<<"] = "<<x[nlines1-skip1]<<"   y["<<nlines1-skip1<<"] = "<<y[nlines1-skip1]<<endl;

	  //Print the squared Q value and magnitude to the new file.
	  output<<x[nlines1-skip1]<<"    "<<y[nlines1-skip1]<<endl;

	  nlines1++;
	}
    }
  fclose(fp1);
  cout<<"nlines1 = "<<nlines1<<endl;

}
