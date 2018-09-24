//Double_t x[10] = {0,1,2,3,4,5,6,7,8,9};

Double_t Test(Double_t *x, Double_t *par)
{
  Double_t val = 0.;

  val = x[0]+par[0];

  return val;
}

void Test_Multi_Plot() 
{
  const Int_t nfunc = 10;
  TF1 **functions = new TF1*[nfunc];
  for (Int_t i=0;i<nfunc;i++) 
    {
      char fname[20];
      sprintf(fname,"f%d",i);
      //functions[i] = new TF1(fname,"gaus",-10,10);
      //functions[i]->SetParameters(1,0,0.5+i/2.);
      functions[i] = new TF1(fname,Test,-10,10,1);
      functions[i]->SetParameter(0,i);
      //functions[i]->SetParameters(1.);
      if (i == 0) functions[i]->Draw();
      else        functions[i]->Draw("lsame");
    }  
}
