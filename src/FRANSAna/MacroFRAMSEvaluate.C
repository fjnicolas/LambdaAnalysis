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

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int MacroFRAMSEvaluate()
{
  double fScoreCut = 0.15;
  //--------- Configuration Parameters
  // Use reco samples with reco vertex
  bool fUseReco=false; //fUseReco=true;

  bool fEvaluateSignal = true; fEvaluateSignal = false;

  //Background label name
  std::string fBGLabel = "";
  fBGLabel = "QE";
  //fBGLabel="Res";
  //fBGLabel="Co";

  // Input weights
  std::string fWeightsPath = "Weights/FRAMSSelectionTMVA_BDT.weights.xml";

  //--------- Input tree
  std::string tree_dirname = "framsTrue/";
  if(fUseReco) tree_dirname = "framsReco/";
  std::string tree_name = "FRAMSTree";
  if(fEvaluateSignal) tree_name+="S";
  //else tree_name+="BG";
  else tree_name+="S";

  std::string BGfilename = "InputFiles_v54/FRAMSAnaOutput_BG"+fBGLabel+".root";
  std::string Sfilename = "InputFiles_v54/FRAMSAnaOutput_S.root";

  // new files
  BGfilename = "InputFiles_v69/QE/analyzeItOutput_R1-1_SR1-1000.root";
  Sfilename = "InputFiles_v69/V0/analyzeItOutput_R1-1_SR1-10.root";

  std::string filename=Sfilename;
  if(!fEvaluateSignal) filename = BGfilename;


  TFile *fInputFile = new TFile(filename.c_str(),"READ");
  TTree *fTree = (TTree *)fInputFile->Get((tree_dirname+tree_name).c_str());

  // Create TMVA::Reader object
  TMVA::Reader *fTMVAReader = new TMVA::Reader();

  // Add variables
  Float_t fAlpha_C, fEta_C, fDelta_C, fFitScore_C;
  Float_t fGap, fProtonKE, fPionKE;

  fTMVAReader->AddVariable( "Alpha_C", &fAlpha_C );
  fTMVAReader->AddVariable( "Eta_C", &fEta_C );
  fTMVAReader->AddVariable( "Delta_C", &fDelta_C);
  fTMVAReader->AddVariable( "FitScore_C", &fFitScore_C );

  fTMVAReader->AddSpectator( "Gap", &fGap );
  fTMVAReader->AddSpectator( "ProtonKE", &fProtonKE );
  fTMVAReader->AddSpectator( "PionKE", &fPionKE );

  fTMVAReader->BookMVA( "FRAMS BDT",  fWeightsPath.c_str()  );


  Double_t Alpha_C, Eta_C, Delta_C, FitScore_C;
  Double_t Gap, ProtonKE, PionKE;
  fTree->SetBranchAddress("Alpha_C", &Alpha_C);
  fTree->SetBranchAddress( "Eta_C", &Eta_C );
  fTree->SetBranchAddress( "Delta_C", &Delta_C);
  fTree->SetBranchAddress( "FitScore_C", &FitScore_C );

  fTree->SetBranchAddress( "Gap", &Gap );
  fTree->SetBranchAddress( "ProtonKE", &ProtonKE );
  fTree->SetBranchAddress( "PionKE", &PionKE );

  //create one-dimensional TEfficiency object with fixed bin size
  TEfficiency* pEff = new TEfficiency("effS","Selection Efficiency;Gap [cm];",21,-1,20);
  bool passCut;

  TH2F hC("hC", ";Gap [cm]; Score", 21, -1, 20, 40, -1., 1.);

  for (Long64_t ievt=0; ievt<fTree->GetEntries();ievt++) {

    fTree->GetEntry(ievt);

    fAlpha_C=Alpha_C; fEta_C=Eta_C; fDelta_C=Delta_C; fFitScore_C=FitScore_C;
    fGap=Gap; fProtonKE=ProtonKE; fPionKE=PionKE;

    // retrieve the corresponding MVA output
    double score = fTMVAReader->EvaluateMVA( "FRAMS BDT" );

    std::cout<<ievt<<" "<<score<<std::endl;

    passCut = score > fScoreCut;
    pEff->Fill(passCut, Gap);

    hC.Fill(Gap, score);


  } // end of event loop

  delete fTMVAReader;

  TCanvas* c1 = new TCanvas("example","",600,400);
  c1->Divide(2,1);
  c1->SetFillStyle(1001);
  c1->SetFillColor(kWhite);
  c1->cd(1);
  pEff->Draw("AP");

  c1->cd(2);
  hC.Draw();

  c1->cd();

  c1->Update(); c1->WaitPrimitive();
  return 0;
}
