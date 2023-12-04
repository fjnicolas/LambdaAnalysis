#include "FRANSAnaCommandLineParser.h"

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

int MacroFRAMSTraining(std::string fInputFileName="", std::string fTreeDirName = "", std::string fTreeName = "FRAMSTree")
{

  //--------- Configuration Parameters
  bool fUseReco=false; fUseReco=true;
  bool fUseBestView=false; //fUseBestView=true;
  bool fUseMultiplicityVars = false; //fUseMultiplicityVars = true;
  
  // Apply quality cuts
  bool fUseQualityCut=false; //fUseQualityCut=true;
  TCut fCutS = ""; //"Gap>1.5 && ProtonKE>0.03 && PionKE>0.02";
  TCut fCutBG = "";

  fCutS = "NOrigins_C<3";
  fCutBG = "NOrigins_C<3";

  // Number of events for training
  int nTrain = 400;

  //Background label name
  std::string fBGLabel = "Inclusive";
  //fBGLabel = "QE";
  //BGLabel="RES";

  // Input file
  TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
  fFile->ls();
  TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
  TTree *fTreeS = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
  TTree *fTreeBG = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());


  //--------- Output file
  gSystem->Exec( "rm -rf TMVAResults" );
  gSystem->Exec( "mkdir TMVAResults" );
  std::string fOutputTMVAROOtFileName = "TMVAResults/TMTMVAResults.root";
  TFile* outputFile = TFile::Open( fOutputTMVAROOtFileName.c_str(), "CREATE" );
  
  
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained
  std::map<std::string,int> Use;
  // Rectangular cut optimisation
  Use["Cuts"]            = 0;
  // 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  // Linear Discriminant Analysis
  Use["Fisher"]          = 0;
  // Boosted Decision Trees
  Use["BDT"]             = 1;
  // Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0;


  // Define factory and data loader
  TMVA::Factory *factory = new TMVA::Factory( "FRAMSSelectionTMVA", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");


  //Add variables for MVA

  if(fUseMultiplicityVars){
    dataloader->AddVariable( "FRANSObj2.fNOrigins_C", "N_{C}", "", 'I' );
    dataloader->AddVariable( "FRANSObj2.fNOriginsM1_C", "N^{1}_{C}", "", 'I' );
    dataloader->AddVariable( "FRANSObj2.fNOriginsM2_C", "N^{2}_{C}", "", 'I' );
    //dataloader->AddVariable( "FRANSObj2.NOriginsM3_C", "N^{>3}_C", "", 'I' );
    dataloader->AddVariable( "FRANSObj2.fHitDensity_C", "d", "", 'D' );
  }


  dataloader->AddVariable( "FRANSObj2.fAlpha", "#alpha_{C}", "", 'D' );
  dataloader->AddVariable( "FRANSObj2.fEta", "#eta_{C}", "", 'D' );
  dataloader->AddVariable( "FRANSObj2.fDelta", "#Delta_{C}", "", 'D' );
  dataloader->AddVariable( "FRANSObj2.fFitScore", "r_{C}", "", 'D' );
  
  

  //Add spectator variables
  dataloader->AddSpectator( "FRANSObj2.fGap", "Gap", "", 'D' );
  dataloader->AddSpectator( "FRANSObj2.fProtonKE", "ProtonKE", "", 'D' );
  dataloader->AddSpectator( "FRANSObj2.fPionKE", "PionKE", "", 'D' );

  // You can add an arbitrary number of signal or background trees
  Double_t signalWeight     = 0.1;
  Double_t backgroundWeight = 1.0;
  dataloader->AddSignalTree    ( fTreeS,     signalWeight );
  dataloader->AddBackgroundTree( fTreeBG, backgroundWeight );

  // define signal and background cuts
  TCut signalCut      = "IsSignal==1";
  TCut bgCutInclusive = "IsSignal!=1";
  TCut bgCutQE = "IsSignal!=1 && IntMode==0";
  TCut bgCutRES = "IsSignal!=1  && IntMode==1";
  TCut bgCut = bgCutInclusive;
  if(fBGLabel=="QE") bgCut = bgCutQE;
  else if(fBGLabel=="RES") bgCut = bgCutRES;

  TCut mycuts = "";
  if(fUseQualityCut){
    signalCut = signalCut + fCutS ;
    bgCut = bgCut + fCutBG ;
  }

  std::string tmva_options = "nTrain_Signal="+std::to_string(nTrain)+":nTrain_Background="+to_string(nTrain);
  tmva_options+=":SplitMode=Random:NormMode=NumEvents:!V";
  dataloader->PrepareTrainingAndTestTree( signalCut, bgCut, tmva_options );


  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
                         "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
  // Boosted Decision Tree
  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                         "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();


  // Save the output
  outputFile->Close();

  delete factory;
  delete dataloader;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( fOutputTMVAROOtFileName.c_str() );

  return 0;
}
