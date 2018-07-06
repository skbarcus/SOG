#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>

void Read_Example_4Col_12C_Names() {

  //Make a new canvas to plot data.
  TCanvas* c=new TCanvas("c");
  c->SetGrid();
  gROOT->Reset();
  
  //Open first file (blank disc).
  FILE *fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/STANF_374.txt","r");
  //Open file to write output to.
  std::ofstream output ("junk.txt", std::ofstream::out);
  
  Float_t thetatemp,qefftemp,sigexptemp,uncertaintytemp;          //Temporary variables to hold data from data file.
  Int_t ncols;          //Set how many columns of data we have in the data file.
  Int_t nlines = 0;     //Counts number of lines in the data file. 
  Int_t skip = 2;       //Gives number of lines to skip at top of data file. 
  char* str[1000];       //Variable to read lines of the data file.
  Float_t theta[39];    //Arrays to hold X and Y data.
  Float_t qeff[39];
  Float_t sigexp[39];
  Float_t uncertainty[39];
  
  //Create a new root file.
  TFile *f = new TFile("basic.root","RECREATE");
  //Create histograms.
  TH1F *h1 = new TH1F("h1","x distribution",23,-11,11);
  TH2F *h2 = new TH2F("h2","x vs. y distribution",2300,-11,11,111,0,110);
  //Create ntuple. 
  TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y");
  
  //Read in blank disc data.
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
	cout<<thetatemp<<"   "<<qefftemp<<"   "<<sigexptemp<<"   "<<uncertaintytemp<<endl;
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
  //Print number of lines with data.
  printf(" found %d points\n",nlines - skip);
  
  //Close data file. 
  fclose(fp);
  f->Write();
  
}

