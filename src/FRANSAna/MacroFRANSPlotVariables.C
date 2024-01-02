#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

vector<TPad*> buildpadcanvas(int nx, int ny){
    vector<TPad*> Tp;
    double x=0, y=1, dx=1./nx, dy=1./ny;
    TPad *pad = new TPad("","", 0, 0, 1, 1, -1, -1, -1);
    Tp.push_back(pad);
    for(int i=1; i<=nx; i++){
      y=1;
      for(int j=1; j<=ny; j++){
        //TPad *pad = new TPad("a", "a",x, y, x+dx,y+dy);
        //TPad pad("a", "a",x, y, x+dx,y+dy,1, 1, 2);
        //cout<<x<<" "<<y<<endl;
        TPad *pad = new TPad("","", x, y-dy, x+dx, y, -1, -1, -1);
        Tp.push_back(pad);
        y-=dy;
        // Bottom margin
        pad->SetBottomMargin(0.15);
        // Left margin
        pad->SetLeftMargin(0.15);
      }
      x+=dx;
      
    }
    for(int i=0; i<=nx*ny; i++){
      Tp.at(i)->Draw();
    }
    return Tp;
}

int MacroFRANSPlotVariables(std::string fInputFileName="", int view = 2, std::string fTreeDirName = "FRANSCheatedVx/", std::string fTreeName = "FRANSTree")
{

  //--------- Configuration Parameters
  bool fUseMultiplicityVars = false; //fUseMultiplicityVars = true;
  std::string fBGLabel = "";  //fBGLabel = "QE"; //fBGLabel = "RES";
  std::string fVw = std::to_string(view);

  //--------- Input file
  TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
  fFile->ls();
  TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

  // Axis offsets
  gStyle->SetTitleOffset(1.25, "Y");
  // Set divisions
  gStyle->SetNdivisions(505, "X");
  // Line widths
  gStyle->SetLineWidth(2);
  // Graphs line widths
  gStyle->SetHistLineWidth(3);
  // Draw variables for signal and background
  TCanvas* c = new TCanvas("c", "Variables", 1000, 700);
  std::vector<TPad*> Tp = buildpadcanvas(3,2);
  Tp.at(1)->cd();

  TH1F* hAlphaFrame = new TH1F("hAlphaFrame", ";#alpha;AU", 50, 0, 10000);
  hAlphaFrame->SetStats(0);
  TH1F* hAlphaSig = new TH1F("hAlphaSig", "Signal;#alpha", hAlphaFrame->GetNbinsX(), hAlphaFrame->GetXaxis()->GetXmin(), hAlphaFrame->GetXaxis()->GetXmax());
  TH1F* hAlphaBG = new TH1F("hAlphaBG", "Background;#alpha", hAlphaFrame->GetNbinsX(), hAlphaFrame->GetXaxis()->GetXmin(), hAlphaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fAlpha"+">>hAlphaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fAlpha"+">>hAlphaBG").c_str(), "IsSignal!=1");
  // Normalize
  hAlphaSig->Scale(1./hAlphaSig->Integral());
  hAlphaBG->Scale(1./hAlphaBG->Integral());
  // Set maximum
  hAlphaFrame->SetMaximum(std::max(hAlphaSig->GetMaximum(), hAlphaBG->GetMaximum())*1.2);
  // Set colors
  hAlphaSig->SetLineColor(kAzure+1);
  hAlphaBG->SetLineColor(kOrange+1);
  hAlphaFrame->Draw();
  hAlphaSig->Draw("hist same");
  hAlphaBG->Draw("hist same");
  // Legend
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hAlphaSig, "Signal", "l");
  legend->AddEntry(hAlphaBG, "Background", "l");
  legend->Draw("same");

  Tp.at(2)->cd();
  TH1F* hEtaFrame = new TH1F("hEtaFrame", ";#eta;AU", 50, 0, 20);
  hEtaFrame->SetStats(0);
  TH1F* hEtaSig = new TH1F("hEtaSig", "Signal;#eta", hEtaFrame->GetNbinsX(), hEtaFrame->GetXaxis()->GetXmin(), hEtaFrame->GetXaxis()->GetXmax());
  TH1F* hEtaBG = new TH1F("hEtaBG", "Background;#eta", hEtaFrame->GetNbinsX(), hEtaFrame->GetXaxis()->GetXmin(), hEtaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fEta"+">>hEtaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fEta"+">>hEtaBG").c_str(), "IsSignal!=1");
  // Normalize
  hEtaSig->Scale(1./hEtaSig->Integral());
  hEtaBG->Scale(1./hEtaBG->Integral());
  // Set maximum
  hEtaFrame->SetMaximum(std::max(hEtaSig->GetMaximum(), hEtaBG->GetMaximum())*1.2);
  // Set colors
  hEtaSig->SetLineColor(kAzure+1);
  hEtaBG->SetLineColor(kOrange+1);
  hEtaFrame->Draw();
  hEtaSig->Draw("hist same");
  hEtaBG->Draw("hist same");
  legend->Draw("same");


  Tp.at(3)->cd();
  TH1F* hDeltaFrame = new TH1F("hDeltaFrame", ";#Delta;AU", 50, 0, 1);
  hDeltaFrame->SetStats(0);
  TH1F* hDeltaSig = new TH1F("hDeltaSig", "Signal;#Delta", hDeltaFrame->GetNbinsX(), hDeltaFrame->GetXaxis()->GetXmin(), hDeltaFrame->GetXaxis()->GetXmax());
  TH1F* hDeltaBG = new TH1F("hDeltaBG", "Background;#Delta", hDeltaFrame->GetNbinsX(), hDeltaFrame->GetXaxis()->GetXmin(), hDeltaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fDelta"+">>hDeltaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fDelta"+">>hDeltaBG").c_str(), "IsSignal!=1");
  // Normalize
  hDeltaSig->Scale(1./hDeltaSig->Integral());
  hDeltaBG->Scale(1./hDeltaBG->Integral());
  // Set maximum
  hDeltaFrame->SetMaximum(std::max(hDeltaSig->GetMaximum(), hDeltaBG->GetMaximum())*1.2);
  // Set colors
  hDeltaSig->SetLineColor(kAzure+1);
  hDeltaBG->SetLineColor(kOrange+1);
  hDeltaFrame->Draw();
  hDeltaSig->Draw("hist same");
  hDeltaBG->Draw("hist same");
  legend->Draw("same");


  Tp.at(4)->cd();
  TH1F* hFitScoreFrame = new TH1F("hFitScoreFrame", ";Fit score;AU", 50, 0, 1);
  hFitScoreFrame->SetStats(0);
  TH1F* hFitScoreSig = new TH1F("hFitScoreSig", "Signal;Fit score", hFitScoreFrame->GetNbinsX(), hFitScoreFrame->GetXaxis()->GetXmin(), hFitScoreFrame->GetXaxis()->GetXmax());
  TH1F* hFitScoreBG = new TH1F("hFitScoreBG", "Background;Fit score", hFitScoreFrame->GetNbinsX(), hFitScoreFrame->GetXaxis()->GetXmin(), hFitScoreFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fFitScore"+">>hFitScoreSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fFitScore"+">>hFitScoreBG").c_str(), "IsSignal!=1");
  // Normalize
  hFitScoreSig->Scale(1./hFitScoreSig->Integral());
  hFitScoreBG->Scale(1./hFitScoreBG->Integral());
  // Set maximum
  hFitScoreFrame->SetMaximum(std::max(hFitScoreSig->GetMaximum(), hFitScoreBG->GetMaximum())*1.2);
  // Set colors
  hFitScoreSig->SetLineColor(kAzure+1);
  hFitScoreBG->SetLineColor(kOrange+1);
  hFitScoreFrame->Draw();
  hFitScoreSig->Draw("hist same");
  hFitScoreBG->Draw("hist same");
  legend->Draw("same");


  Tp.at(5)->cd();
  TH1F* hIotaFrame = new TH1F("hIotaFrame", ";#iota;AU", 20, 0, 10);
  hIotaFrame->SetStats(0);
  TH1F* hIotaSig = new TH1F("hIotaSig", "Signal;#iota", hIotaFrame->GetNbinsX(), hIotaFrame->GetXaxis()->GetXmin(), hIotaFrame->GetXaxis()->GetXmax());
  TH1F* hIotaBG = new TH1F("hIotaBG", "Background;#iota", hIotaFrame->GetNbinsX(), hIotaFrame->GetXaxis()->GetXmin(), hIotaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fIota"+">>hIotaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fIota"+">>hIotaBG").c_str(), "IsSignal!=1");
  // Normalize
  hIotaSig->Scale(1./hIotaSig->Integral());
  hIotaBG->Scale(1./hIotaBG->Integral());
  // Set maximum
  hIotaFrame->SetMaximum(std::max(hIotaSig->GetMaximum(), hIotaBG->GetMaximum())*1.2);
  // Set colors
  hIotaSig->SetLineColor(kAzure+1);
  hIotaBG->SetLineColor(kOrange+1);
  hIotaFrame->Draw();
  hIotaSig->Draw("hist same");
  hIotaBG->Draw("hist same");
  legend->Draw("same");


  c->Update();
  c->WaitPrimitive();
  // Save as pdf
  c->SaveAs(("FRANSPlotVariables_"+fVw+".pdf").c_str());

  return 0;
}
