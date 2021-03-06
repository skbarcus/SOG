#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>

void Read_Example() {

  //Make a new canvas to plot data.
  TCanvas* c=new TCanvas("c");
  c->SetGrid();
  gROOT->Reset();
  
  //Open first file (blank disc).
  FILE *fp = fopen("/home/skbarcus/Grinch/WLS_Graphs/Absorption_Tests/spectra_HR4C0611_1-140/Read_Example.txt","r");
  //Open file to write output to.
  std::ofstream output ("output.txt", std::ofstream::out);
  
  Float_t x,y;          //Temporary variables to hold data from data file.
  Int_t ncols;          //Set how many columns of data we have in the data file.
  Int_t nlines = 0;     //Counts number of lines in the data file. 
  Int_t skip = 5;       //Gives number of lines to skip at top of data file. 
  char* str[100];       //Variable to read lines of the data file.
  double xval[21];    //Arrays to hold X and Y data.
  double yval[21];
  
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
	fgets(str,100,fp);
	nlines++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(fp,"%f %f",&x, &y);
	if (ncols < 0) break;    
	xval[nlines-skip] = x;
	yval[nlines-skip] = y;
	//Fill histograms with x and y data.
	h1->Fill(x);
	h2->Fill(x,y);
	//h2->SetMarkerSize(5);
	//Fill ntuple with x and y data.
	ntuple->Fill(x,y);
	//Count the number of lines in the file. 
	nlines++;
      }
  }
  
  //Print the data read from the file. 
  for(int i=0; i<21; i++)
    {
      cout<<"xval["<<i<<"]= "<<xval[i]<<" yval["<<i<<"]= "<<yval[i]<<endl;
    }
  //Print number of lines with data.
  printf(" found %d points\n",nlines - skip);
  
  //Close data file. 
  fclose(fp);
  f->Write();
  
  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph = new TGraph(nlines-skip,xval,yval);
  //Draw the new TGraph called graph on the canvas. 
  graph->Draw();
  
  //Set X axis
  graph->GetXaxis()->SetLimits(-12,12);
  //Set Y axis Min and Max (not sure why different from X).
  graph->SetMinimum(0);
  graph->SetMaximum(120);
  graph->SetLineWidth(3);
  graph->SetLineColor(4);
  graph->SetFillColor(0);
  graph->SetTitle("Main Title; X Axis Title; Y Axis Title");
  //graph_expected.SetFillColor(kYellow);
  //graph_expected.DrawClone("E3AL"); // E3 draws the band

  // Draw the Legend
  TLegend leg(0.9,.7,.56,.9,"Legend Title");
  leg.SetFillColor(0);
  leg.AddEntry(graph,"Curve Name");
  leg.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  //Write absorption data to text file.
  //Writes two columns of data. Takes input data and multiplies it by 2.
  for(i=0; i<(nlines - skip); i++)
    {
	  output<<xval[i]*2<<"   "<<yval[i]*2<<endl;  
    }
  output.close();
  
}

