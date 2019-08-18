#include <iostream>

#include <TROOT.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>

const Double_t g_MinFP = 0.002;

void draw2oo3(){

  TCanvas* canvas = new TCanvas();
  TH1F* hFrame = new TH1F("hFrame","improvement with 2-out-of-3, 3-out-of-4 and 2-out-of-4 logic",
			  100,g_MinFP,1); 
  hFrame->SetStats(false);

  TF1* f_2oo3 = new TF1( "f_2oo3", "[0]*x*x+[1]*x", g_MinFP, 1 );
  f_2oo3->SetParameters(-2,3);
  f_2oo3->SetLineColor(1);
  TF1* f_3oo4 = new TF1( "f_3oo4", "[0]*x*x*x+[1]*x*x", g_MinFP, 1 );
  f_3oo4->SetParameters(-3,4);
  f_3oo4->SetLineColor(3);
  TF1* f_2oo4 = new TF1( "f_2oo4", "[0]*x*x*x+[1]*x*x+[2]*x", g_MinFP, 1 );
  f_2oo4->SetParameters(3,-6,8);
  f_2oo4->SetLineColor(4);

  hFrame->Draw();
  f_2oo3->Draw("same");
  f_3oo4->Draw("same");
  f_2oo4->Draw("same");
  hFrame->GetXaxis()->SetTitle("#it{f}");
  hFrame->GetYaxis()->SetTitle("#it{F/f}");
  hFrame->GetYaxis()->SetRangeUser(2e-4,8);
  canvas->SetLogx();
  canvas->SetLogy();

  TLegend* legend = new TLegend(0.15,0.9,0.4,0.6); 
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(f_2oo3,"2-out-of-3 logic","L");
  legend->AddEntry(f_3oo4,"3-out-of-4 logic","L");
  legend->AddEntry(f_2oo4,"2-out-of-4 logic","L");
  legend->Draw();

  return;
}
