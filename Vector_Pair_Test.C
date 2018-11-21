
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

//root -l Vector_Pair_Test.C+ works for executing this code.
Int_t skip = 0.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file.
Int_t ncols;                             //Set how many columns of data we have in the data file.
char str[1000];                          //Variable to read lines of the data file.//No pointer now if I want to compile.
const Int_t size = 1000;
Float_t Fit[size], Chi2[size];
Float_t Fittemp,Chi2temp;

void Vector_Pair_Test() 
{
  using namespace std;

  FILE *fp;
  fp = fopen("/home/skbarcus/Tritium/Analysis/SOG/Chi2_Sorted.txt","r");

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
	ncols = fscanf(fp,"%f %f", &Fittemp, &Chi2temp);
	  //cout<<"ncols = "<<ncols<<endl;
	  if (ncols < 0) break;    
	  
	  Fit[nlines-skip] = Fittemp;               //Only if using representative fit text file. 
	  Chi2[nlines-skip] = Chi2temp;
	  
	  nlines++;
	}
    }

  vector< pair <Double_t,Int_t> > Test;
  pair<Double_t,Int_t> pair;

  /*
  pair.first = 2;
  pair.second = 1;
  Test.push_back(pair);
  pair.first = 5;
  pair.second = 2;
  Test.push_back(pair);
  */

  for(Int_t i = 0;i<nlines-skip;i++)
    {
      pair.first = Chi2[i];
      pair.second = Fit[i];
      Test.push_back(pair);
    }

  sort(Test.begin(),Test.end());

  for(Int_t i = 0;i<nlines-skip;i++)
    {
      cout<<Test[i].first<<"   "<<Test[i].second<<endl;
    }
}

/*
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>

typedef std::pair<int, int> pair;

int Vector_Pair_Test()
{
	std::vector<pair> v = { { 1, 2 }, { 6, 4 }, { 3, 4 }, { 6, 1 } };

	// sorts pairs in increasing order of their first value
	std::sort(v.begin(), v.end());

	// sorts pairs in decreasing order of their first value
	// std::sort(v.begin(), v.end(), std::greater<>());

	//for (const pair &p: v) {
	//	std::cout << '{' << p.first << ',' << p.second << '}' << '\n';
	//}

	return 0;
}
*/
