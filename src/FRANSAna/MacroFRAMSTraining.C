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

int MacroFRAMSTraining(std::string fInputFileName="")
{

  //--------- Configuration Parameters
  // Use reco samples with reco vertex
  bool fUseReco=false; 
  fUseReco=true;

  // Apply truth gap cut
  bool fUseCut=false; //fUseCut=true;
  TCut fCutS = "Gap>1."; //"Gap>1.5 && ProtonKE>0.03 && PionKE>0.02";
  TCut fCutBG = "";

  bool fUseBestView=false; //fUseBestView=true;

  // Number of events for training
  int nTrain = 350;

  //Background label name
  std::string fBGLabel = "Inclusive";
  //fBGLabel = "QE";
  //BGLabel="RES";

  //--------- Input trees
  std::string tree_dirname = "";
  tree_dirname = "framsTrue/";
  if(fUseReco) tree_dirname = "framsReco/";
  tree_dirname = "";
  std::string tree_name = "FRAMSTree";

  /*std::string BGfilename = "InputFiles_v54/FRAMSAnaOutput_BG"+fBGLabel+".root";
  std::string Sfilename = "InputFiles_v54/FRAMSAnaOutput_S.root";
  // new files
  BGfilename = "InputFiles_v69/QE/analyzeItOutput_R1-1_SR1-1000.root";
  Sfilename = "InputFiles_v69/V0/analyzeItOutput_R1-1_SR1-10.root";
  TFile *fBG = new TFile(BGfilename.c_str(),"READ");
  //TTree *background = (TTree *)fBG->Get((tree_dirname+tree_name+"BG").c_str());
  TTree *background = (TTree *)fBG->Get((tree_dirname+tree_name+"S").c_str());
  TFile *fS = new TFile(Sfilename.c_str(),"READ");
  TTree *signalTree = (TTree *)fS->Get((tree_dirname+tree_name+"S").c_str());*/

  // Input file
  TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
  fFile->ls();
  TTree *fTree = (TTree *)fFile->Get((tree_dirname+tree_name).c_str());


  //--------- Output file
  bool fUseCosmics=(fBGLabel=="Co");
  std::string outputTMVAROOtFile = "TMVAResults/TMVAResults";
  if(fUseReco) outputTMVAROOtFile = outputTMVAROOtFile + "_reco";
  else outputTMVAROOtFile = outputTMVAROOtFile + "_true";
  if(fUseCut) outputTMVAROOtFile = outputTMVAROOtFile + "_cut";
  if(fUseCosmics)  outputTMVAROOtFile = outputTMVAROOtFile + "_cosmics";
  outputTMVAROOtFile = outputTMVAROOtFile + ".root";
  // ROOT file with TMVA output
  TString outfileName( outputTMVAROOtFile.c_str() );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained
  std::map<std::string,int> Use;
  // Rectangular cut optimisation
  Use["Cuts"]            = 0;
  // 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  // Linear Discriminant Analysis
  Use["Fisher"]          = 1;
  // Boosted Decision Trees
  Use["BDT"]             = 1;
  // Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0;


  // Define factory and data loader
  TMVA::Factory *factory = new TMVA::Factory( "FRAMSSelectionTMVA", outputFile,
                                             "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");


  //Add variables for MVA
  if(fUseBestView){
    dataloader->AddVariable( "Alpha", "#alpha", "", 'D' );
    dataloader->AddVariable( "Eta", "#eta", "", 'D' );
    dataloader->AddVariable( "Delta", "#Delta", "", 'D' );
    dataloader->AddVariable( "FitScore", "r", "", 'D' );
  }
  else{
    dataloader->AddVariable( "Alpha_C", "#alpha_{C}", "", 'D' );
    dataloader->AddVariable( "Eta_C", "#eta_{C}", "", 'D' );
    dataloader->AddVariable( "Delta_C", "#Delta_{C}", "", 'D' );
    dataloader->AddVariable( "FitScore_C", "r_{C}", "", 'D' );
  }
  

  //Add spectator variables
  dataloader->AddSpectator( "Gap", "Gap", "", 'D' );
  dataloader->AddSpectator( "ProtonKE", "ProtonKE", "", 'D' );
  dataloader->AddSpectator( "PionKE", "PionKE", "", 'D' );

  // You can add an arbitrary number of signal or background trees
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  /*dataloader->AddSignalTree    ( signalTree,     signalWeight );
  dataloader->AddBackgroundTree( background, backgroundWeight );*/
  dataloader->AddSignalTree    ( fTree,     signalWeight );
  dataloader->AddBackgroundTree( fTree, backgroundWeight );

  // define signal and background
  TCut signalCut      = "IsSignal==1";
  TCut bgCutInclusive = "IsSignal!=1";
  TCut bgCutQE = "IsSignal!=1 && IntMode==0";
  TCut bgCutRES = "IsSignal!=1  && IntMode==1";
  TCut bgCut = bgCutInclusive;
  if(fBGLabel=="QE") bgCut = bgCutQE;
  else if(fBGLabel=="RES") bgCut = bgCutRES;
  dataloader->PrepareTrainingAndTestTree( signalCut, bgCut, "SplitMode=random:!V" );

  // Prepare for training
  TCut mycuts = "";
  if(fUseCut) mycuts = fCutS;
  TCut mycutb = fCutBG;
  std::string tmva_options;
  tmva_options = "nTrain_Signal="+std::to_string(nTrain)+":nTrain_Background="+to_string(nTrain)+":SplitMode=Random:NormMode=NumEvents:!V";
  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, tmva_options );


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
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

  return 0;
}
