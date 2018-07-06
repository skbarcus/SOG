
// Simple macro showing capabilities of triple slider
//Authors: Bertrand Bellenot, Ilka Antcheva
//Modified by Scott Barcus Feb 2017 to examine how cuts on the Cherenkov spectrum influence the shower and preshower.

#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGTripleSlider.h"

enum ETestCommandIdentifiers {
   HId1,
   HId2,
   HId3,
   HCId1,
   HCId2,
   HSId1
};

class TTripleSliderDemo : public TGMainFrame {

private:
   TRootEmbeddedCanvas *fCanvas;
   TGLayoutHints       *fLcan;
   TF1                 *fFitFcn;
   TH1                 *histo1;


   TGHorizontalFrame   *fHframe0, *fHframe1, *fHframe2;
   TGLayoutHints       *fBly, *fBfly1, *fBfly2, *fBfly3;
   TGTripleHSlider     *fHslider1;
   TGTextEntry         *fTeh1, *fTeh2, *fTeh3;
   TGTextBuffer        *fTbh1, *fTbh2, *fTbh3;
   TGCheckButton       *fCheck1, *fCheck2;

public:
   TTripleSliderDemo();
   virtual ~TTripleSliderDemo();

   void CloseWindow();
   void DoText(const char *text);
   void DoSlider();
   void HandleButtons();

   ClassDef(TTripleSliderDemo, 0)
};

//______________________________________________________________________________
TTripleSliderDemo::TTripleSliderDemo() : TGMainFrame(gClient->GetRoot(), 100, 100)
{

   char buf[32];
   SetCleanup(kDeepCleanup);
   // Create an embedded canvas and add to the main frame, centered in x and y
   // and with 30 pixel margins all around
   fCanvas = new TRootEmbeddedCanvas("Canvas", this, 1200, 600);
   
   TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.0,0.96);
   pad1->Draw();
   pad1->Divide(2,1);
   pad1->cd(1);
   pad1->cd(1)->SetLogy();

   fLcan = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
   AddFrame(fCanvas, fLcan);
   fCanvas->GetCanvas()->SetFillColor(33);
   fCanvas->GetCanvas()->SetFrameFillColor(41);
   fCanvas->GetCanvas()->SetBorderMode(0);
   fCanvas->GetCanvas()->SetGrid();
   fCanvas->GetCanvas()->SetLogy();

   fHframe0 = new TGHorizontalFrame(this, 0, 0, 0);

   fCheck1 = new TGCheckButton(fHframe0, "&Constrained", HCId1);
   fCheck2 = new TGCheckButton(fHframe0, "&Relative", HCId2);
   fCheck1->SetState(kButtonUp);
   fCheck2->SetState(kButtonUp);
   fCheck1->SetToolTipText("Pointer position constrained to slider sides");
   fCheck2->SetToolTipText("Pointer position relative to slider position");

   fHframe0->Resize(200, 50);

   fHframe1 = new TGHorizontalFrame(this, 0, 0, 0);

   fHslider1 = new TGTripleHSlider(fHframe1, 190, kDoubleScaleBoth, HSId1,
                                   kHorizontalFrame,
                                   GetDefaultFrameBackground(),
                                   kFALSE, kFALSE, kFALSE, kFALSE);
   fHslider1->Connect("PointerPositionChanged()", "TTripleSliderDemo",
                      this, "DoSlider()");
   fHslider1->Connect("PositionChanged()", "TTripleSliderDemo",
                      this, "DoSlider()");
   //Set slider range.
   fHslider1->SetRange(0.0,3000.0);

   fHframe1->Resize(200, 25);

   fHframe2 = new TGHorizontalFrame(this, 0, 0, 0);

   fTeh1 = new TGTextEntry(fHframe2, fTbh1 = new TGTextBuffer(5), HId1);
   fTeh2 = new TGTextEntry(fHframe2, fTbh2 = new TGTextBuffer(5), HId2);
   fTeh3 = new TGTextEntry(fHframe2, fTbh3 = new TGTextBuffer(5), HId3);

   fTeh1->SetToolTipText("Minimum (left) Value of Slider");
   fTeh2->SetToolTipText("Pointer Position Value");
   fTeh3->SetToolTipText("Maximum (right) Value of Slider");

   fTbh1->AddText(0, "0.0");
   fTbh2->AddText(0, "0.0");
   fTbh3->AddText(0, "0.0");

   fTeh1->Connect("TextChanged(char*)", "TTripleSliderDemo", this,
                  "DoText(char*)");
   fTeh2->Connect("TextChanged(char*)", "TTripleSliderDemo", this,
                  "DoText(char*)");
   fTeh3->Connect("TextChanged(char*)", "TTripleSliderDemo", this,
                  "DoText(char*)");

   fCheck1->Connect("Clicked()", "TTripleSliderDemo", this,
                    "HandleButtons()");
   fCheck2->Connect("Clicked()", "TTripleSliderDemo", this,
                    "HandleButtons()");

   fHframe2->Resize(100, 25);

   //--- layout for buttons: top align, equally expand horizontally
   fBly = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5);

   //--- layout for the frame: place at bottom, right aligned
   fBfly1 = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 5, 5);
   fBfly2 = new TGLayoutHints(kLHintsTop | kLHintsLeft,    5, 5, 5, 5);
   fBfly3 = new TGLayoutHints(kLHintsTop | kLHintsRight,   5, 5, 5, 5);

   fHframe0->AddFrame(fCheck1, fBfly2);
   fHframe0->AddFrame(fCheck2, fBfly2);
   fHframe1->AddFrame(fHslider1, fBly);
   fHframe2->AddFrame(fTeh1, fBfly2);
   fHframe2->AddFrame(fTeh2, fBfly1);
   fHframe2->AddFrame(fTeh3, fBfly3);

   AddFrame(fHframe0, fBly);
   AddFrame(fHframe1, fBly);
   AddFrame(fHframe2, fBly);

   // Set main frame name, map sub windows (buttons), initialize layout
   // algorithm via Resize() and map main frame
   SetWindowName("Cherenkov Cuts on Shower and Preshower");
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();

   //Set slider positions.
   fHslider1->SetPosition(0.0,1500);
   fHslider1->SetPointerPosition(0.0);

   //Print slider positions to boxes.
   sprintf(buf, "%.3f", fHslider1->GetMinPosition());
   fTbh1->Clear();
   fTbh1->AddText(0, buf);
   sprintf(buf, "%.3f", fHslider1->GetPointerPosition());
   fTbh2->Clear();
   fTbh2->AddText(0, buf);
   sprintf(buf, "%.3f", fHslider1->GetMaxPosition());
   fTbh3->Clear();
   fTbh3->AddText(0, buf);

   //Draw Cherenkov Spectrum with a line representing the cut.
   Double_t ymax = 10000;
   T->Draw("L.cer.asum_c>>h1(500,0,3000)");
   h1->GetYaxis()->SetTitle("Counts");
   h1->GetXaxis()->SetTitle("L.cer.asum_c");
   h1->SetMaximum(ymax);
   h1->GetYaxis()->SetTitleOffset(1.3);

   Int_t cermin = 0;
   TLine *line1 = new TLine(cermin, 0, cermin, ymax);
   line1->SetLineColor(kRed);
   line1->Draw("SAME");

   pad1->cd(2);
   
   //Draw shower and preshower with no cut yet applied.
   T->Draw("L.prl1.e:L.prl2.e>>h2(500,0,1500,500,0,1500)","","colz");
   h2->GetYaxis()->SetTitle("L.prl1.e");
   h2->GetXaxis()->SetTitle("L.prl2.e");
   h2->GetYaxis()->SetTitleOffset(1.55);

}

//______________________________________________________________________________
TTripleSliderDemo::~TTripleSliderDemo()
{
   // Clean up

   Cleanup();
}

//______________________________________________________________________________
void TTripleSliderDemo::CloseWindow()
{
   // Called when window is closed via the window manager.

   delete this;
}

//______________________________________________________________________________
void TTripleSliderDemo::DoText(const char * /*text*/)
{
   // Handle text entry widgets.

   TGTextEntry *te = (TGTextEntry *) gTQSender;
   Int_t id = te->WidgetId();

   switch (id) {
      case HId1:
         fHslider1->SetPosition(atof(fTbh1->GetString()),
                                fHslider1->GetMaxPosition());
         break;
      case HId2:
         fHslider1->SetPointerPosition(atof(fTbh2->GetString()));
         break;
      case HId3:
         fHslider1->SetPosition(fHslider1->GetMinPosition(),
                                atof(fTbh1->GetString()));
         break;
      default:
         break;
   }
   fCanvas->GetCanvas()->Modified();
   fCanvas->GetCanvas()->Update();
}

//______________________________________________________________________________
void TTripleSliderDemo::DoSlider()
{
   // Handle slider widgets.
   char buf[32];

   //Update text boxes.
   sprintf(buf, "%.3f", fHslider1->GetMinPosition());
   fTbh1->Clear();
   fTbh1->AddText(0, buf);
   fTeh1->SetCursorPosition(fTeh1->GetCursorPosition());
   fTeh1->Deselect();
   gClient->NeedRedraw(fTeh1);

   sprintf(buf, "%.3f", fHslider1->GetPointerPosition());
   fTbh2->Clear();
   fTbh2->AddText(0, buf);
   fTeh2->SetCursorPosition(fTeh2->GetCursorPosition());
   fTeh2->Deselect();
   gClient->NeedRedraw(fTeh2);

   sprintf(buf, "%.3f", fHslider1->GetMaxPosition());
   fTbh3->Clear();
   fTbh3->AddText(0, buf);
   fTeh3->SetCursorPosition(fTeh3->GetCursorPosition());
   fTeh3->Deselect();
   gClient->NeedRedraw(fTeh3);


   Double_t ymax = 10000;

   pad1->cd(1);
   //Draw the Cherenkov spectrum.
   T->Draw("L.cer.asum_c>>h1(500,0,3000)");
   h1->GetYaxis()->SetTitle("Counts");
   h1->GetXaxis()->SetTitle("L.cer.asum_c");
   h1->SetMaximum(ymax);
   h1->GetYaxis()->SetTitleOffset(1.3);

   Int_t cermin = 300;
   Double_t slidervalue = 0;
   slidervalue = fHslider1->GetPointerPosition();

   //Draw a line representing the Cherenkov cut and tie it to the slider.
   TLine *line1 = new TLine(cermin, 0, cermin, ymax);
   line1->SetLineColor(kRed);
   line1->DrawLine(slidervalue,0,slidervalue,ymax);
   
   pad1->cd(2);

   //Create a cut on the Cherenkov spectrum and tie it to the current value of the slider.
   TCut ct1 = Form("L.cer.asum_c>%f",slidervalue);  

   //Draw the shower and preshower with the Cherenkov cut applied.
   T->Draw("L.prl1.e:L.prl2.e>>h2(500,0,1500,500,0,1500)",ct1,"colz");
   h2->GetYaxis()->SetTitle("L.prl1.e");
   h2->GetXaxis()->SetTitle("L.prl2.e");
   h2->GetYaxis()->SetTitleOffset(1.55);

   //Set slider min and max to adjust the X axis range of second histogram.
   h2->GetXaxis()->SetRangeUser(fHslider1->GetMinPosition(), fHslider1->GetMaxPosition());

   fCanvas->GetCanvas()->Modified();
   fCanvas->GetCanvas()->Update();
}

//______________________________________________________________________________
void TTripleSliderDemo::HandleButtons()
{
   // Handle different buttons.

   TGButton *btn = (TGButton *) gTQSender;
   Int_t id = btn->WidgetId();

   switch (id) {
      case HCId1:
         fHslider1->SetConstrained(fCheck1->GetState());
         break;
      case HCId2:
         fHslider1->SetRelative(fCheck2->GetState());
         break;
      default:
         break;
   }
}


void Cherenkov_Cut_11386()
{
  //Attach ROOT file with the data to be analyzed.
  TChain *T = new TChain("T");
  T->Add("./rootfiles/left_gmp_11386.root");
  T->Add("./rootfiles/left_gmp_11386_1.root");

  new TTripleSliderDemo(); 
}
