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

#include "CutEfficienciesDefinitions.C"

using namespace TMVA::Experimental;

//---------  Main function
void LambdaBDTAnalysis(std::string fInputFileName="", bool useBatchMode=false, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    // Number of events for training
    int nTrain = 4000;

    TCut fTruthInFV("TruthIsFiducial==1");
    TCut fTruthInAV("abs(NuvX)<200 && abs(NuvY)<200 && NuvZ>0 && NuvZ<500");

    // Batch mode
    useBatchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Output file
    gSystem->Exec( "rm -rf TMVAResults" );
    gSystem->Exec( "mkdir TMVAResults" );
    std::string fOutputTMVAROOtFileName = "TMVAResults/TMTMVAResults.root";
    TFile* outputFile = TFile::Open( fOutputTMVAROOtFileName.c_str(), "CREATE" );

    // Input file
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    // This loads the library
    TMVA::Tools::Instance();

    // Default MVA methods to be trained
    std::map<std::string,int> Use;
    // Rectangular cut optimisation
    Use["Cuts"]            = 1;
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


    dataloader->AddVariable( "AngleFRANSScore", "AngleFRANSScore", "", 'D' );
    //dataloader->AddVariable( "NUnOriginsMultGT3", "N^{1}_{C}", "", 'I' );
    dataloader->AddVariable( "NUnOrigins", "N^{2}_{C}", "", 'I' );
    //dataloader->AddVariable( "CRUMBSScore", "CRUMBS", "", 'D' );
    dataloader->AddVariable( "AngleDecayContainedDiff", "#alpha", "", 'D' );
    dataloader->AddVariable( "AngleLengthMainTrack", "Main track Length [cm]", "", 'I' );
    //dataloader->AddVariable( "AngleLengthTrack1", "Track 1 Length [cm]", "", 'I' );
    //dataloader->AddVariable( "AngleLengthTrack2", "Track 2 Length [cm]", "", 'I' );
    //dataloader->AddVariable( "AngleNHitsMainTrack", "Main track # hits", "", 'I' );
    //dataloader->AddVariable( "AngleNHitsTrack1", "Track 1 # hits", "", 'I' );
    //dataloader->AddVariable( "AngleNHitsTrack2", "Track 2 # hits", "", 'I' );
    dataloader->AddVariable( "AngleMinNHits", "# hits min", "", 'I' );
    dataloader->AddVariable( "NUnassociatedHits", "# unassociated hits", "", 'I' );
    dataloader->AddVariable( "FRANSScorePANDORA", "FRANS PANDORA", "", 'D' );
    //dataloader->AddVariable( "AngleDirtHits", "Dirt Hits", "", 'I' );
    //dataloader->AddVariable( "NShowers", "# showers", "", 'I' );
    dataloader->AddVariable( "NShowerHits", "# shower hits", "", 'I' );
    //dataloader->AddVariable( "AngleLongestIsMain", "LongestIsMain", "", 'I' );
    //dataloader->AddVariable( "ShowerEnergy", "Shower Energy", "", 'D' );
    //dataloader->AddVariable( "AngleOpeningAngle", "Opening Angle [º]", "", 'I' );
    //dataloader->AddVariable( "AngleCoveredArea", "Covered Area", "", 'D' );
    


    //dataloader->AddVariable( "AngleGap", "Gap", "", 'D' );
    //dataloader->AddVariable( "ShowerEnergy", "Shower Energy", "", 'D' );
    //dataloader->AddVariable( "FRANSScorePANDORA", "FRANS PANDORA", "", 'D' );
    

    // You can add an arbitrary number of signal or background trees
    std::string fMinimalCut = "RecoIsFiducial && NOriginsPairOneTwo>0 && NAngles>=1";

    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    dataloader->AddSignalTree    ( fTree,     signalWeight );
    dataloader->AddBackgroundTree( fTree, backgroundWeight );

    dataloader->AddCut(TCut(fMinimalCut.c_str()));
    // define signal and background cuts
    TCut signalCut = fTruthInAV + TCut("IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)");//+TCut(fMinimalCut.c_str());
    TCut bgCut = fTruthInAV && TCut("IntOrigin==1 && IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)");//+TCut(fMinimalCut.c_str());

    std::string tmva_options = "";//"nTrain_Signal="+std::to_string(nTrain)+":nTrain_Background="+to_string(nTrain);
    tmva_options+=":SplitMode=Random:NormMode=NumEvents:!V";
    dataloader->PrepareTrainingAndTestTree( signalCut, bgCut, tmva_options );

    std::cout<<"Training and test trees prepared\n";


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


}