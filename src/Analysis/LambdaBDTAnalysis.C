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

std::map<std::string, bool> extractLabels(const std::string& xmlDataFile) {
    std::map<std::string, bool> labelVector;

    std::ifstream file(xmlDataFile);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return labelVector;
    }

    std::string xmlData((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    size_t pos = 0;

    // Search for the beginning of the <Variable> tags
    while ((pos = xmlData.find("<Variable", pos)) != std::string::npos) {
        // Find the Label attribute within the current <Variable> tag
        size_t labelPos = xmlData.find("Label=\"", pos);
        if (labelPos != std::string::npos) {
            labelPos += 7; // Move past the "Label=\"" part
            // Find the closing quote of the Label attribute
            size_t endQuotePos = xmlData.find("\"", labelPos);
            if (endQuotePos != std::string::npos) {
                // Extract the Label value and add it to the vector
                std::string label = xmlData.substr(labelPos, endQuotePos - labelPos);
                labelVector[label] = true;
            }
        }

        // Move to the next position after the current <Variable> tag
        pos += 1;
    }

    return labelVector;
}

//---------  Main function
void LambdaBDTAnalysis(std::string fInputFileName="", double nTrainFrac = -1, bool useBatchMode=false, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
   
    TCut fTruthInFV("TruthIsFiducial==1");
    TCut fTruthInAV("TruthIsAV==1");

    // ---- Minimal cut
    TCut fMinimalCut("RecoIsFiducial && NAngles>=1 && AnglePassFit && AnglePassChargeFit");
    //TCut fSelCuts("FRANSScorePANDORA>=0.2 && AngleDecayContainedDiff<1 && NUnOrigins<1");

    // Batch mode
    useBatchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Output file
    gSystem->Exec( "rm -rf TMVAResultsTest" );
    gSystem->Exec( "mkdir TMVAResultsTest" );
    std::string fOutputTMVAROOtFileName = "TMVAResultsTest/TMVAResults.root";
    TFile* outputFile = TFile::Open( fOutputTMVAROOtFileName.c_str(), "CREATE" );

    // Input file
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

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
    // Random Forest
    Use["RandomForest"]    = 0;

    
    // Define factory and data loader
    TMVA::Factory *factory = new TMVA::Factory( "FRAMSSelectionTMVA", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");


    dataloader->AddVariable( "AngleFRANSScore", "AngleFRANSScore", "", 'D' );
    //dataloader->AddVariable( "NUnOriginsMultGT3", "N^{2}", "", 'I' );
    //dataloader->AddVariable( "NUnOrigins", "N", "", 'I' );
    //dataloader->AddVariable( "CRUMBSScore", "CRUMBS", "", 'D' );
    dataloader->AddVariable( "AngleDecayContainedDiff", "#alpha", "", 'D' );
    dataloader->AddVariable( "AngleLengthMainTrack", "Main track Length [cm]", "", 'I' );
    //dataloader->AddVariable( "AngleLengthTrack1", "Track 1 Length [cm]", "", 'I' );
    //dataloader->AddVariable( "AngleLengthTrack2", "Track 2 Length [cm]", "", 'I' );
    //dataloader->AddVariable( "AngleNHitsMainTrack", "Main track # hits", "", 'I' );
    //dataloader->AddVariable( "AngleNHitsTrack1", "Track 1 # hits", "", 'I' );
    //dataloader->AddVariable( "AngleNHitsTrack2", "Track 2 # hits", "", 'I' );
    //dataloader->AddVariable( "AngleMinNHits", "# hits min", "", 'I' );
    //dataloader->AddVariable( "NUnassociatedHits", "# unassociated hits", "", 'I' );
    //dataloader->AddVariable( "FRANSScorePANDORA", "FRANS PANDORA", "", 'D' );
    //dataloader->AddVariable( "AngleDirtHits", "Dirt Hits", "", 'I' );
    dataloader->AddVariable( "NShowers", "# showers", "", 'I' );
    //dataloader->AddVariable( "NShowerHits", "# shower hits", "", 'I' );
    //dataloader->AddVariable( "AngleLongestIsMain", "LongestIsMain", "", 'I' );
    //dataloader->AddVariable( "ShowerEnergy", "Shower Energy", "", 'D' );
    dataloader->AddVariable( "AngleOpeningAngle", "Opening Angle [º]", "", 'I' );
    //dataloader->AddVariable( "AngleCoveredArea", "Covered Area", "", 'D' );
    //dataloader->AddVariable( "AngleGap", "Gap", "", 'D' );
    //dataloader->AddVariable( "ShowerEnergy", "Shower Energy", "", 'D' );
    // Calorimetry
    //dataloader->AddVariable( "AngleTwoLinesChi2", "Two Lines Chi2", "", 'D' );
    //dataloader->AddVariable( "AnglePassFit", "PassFit", "", 'I' );
    //dataloader->AddVariable( "AnglePassChargeFit", "Pass Charge Fit", "", 'I' );
    //dataloader->AddVariable( "AngleBandOverlap", "B and Overlap", "", 'D' );
    
    dataloader->AddVariable( "AngleBandCrossHits", "Band Overlap", "", 'D' );
    
    //dataloader->AddVariable( "AngleChargeRatioFit", "Charge Ratio Fit", "", 'D' );
    //dataloader->AddVariable( "AngleChargeDifferenceFit", "Charge Difference Fit", "", 'D' );
    
    //dataloader->AddVariable( "AngleChargeRatioIntegral", "Charge Ratio Integral", "", 'D' );
    
    //dataloader->AddVariable( "AngleChargeDifferenceIntegral", "Charge Difference", "", 'D' );
    //dataloader->AddVariable( "AngleChargeRatioAverage", "Charge Ratio Average", "", 'D' );
    
    dataloader->AddVariable( "AngleVertexHitIntegralRatio", "Vertex Hit Integral Ratio", "", 'D' );
    
    //dataloader->AddVariable( "AngleTrackLengthRatio", "Track Length Ratio", "", 'D' );
    //dataloader->AddVariable( "AngleResidualRange1RMS", "Residuals Range 1 RMS", "", 'D' );
    //dataloader->AddVariable( "AngleResidualRange2RMS", "Residuals Range 2 RMS", "", 'D' );
    //dataloader->AddVariable( "AngleResidualRangeMinRMS", "Residuals Range Min RMS", "", 'D' );
    
    dataloader->AddVariable( "AngleNVertexHits", "N Vertex Hits", "", 'I' );
    dataloader->AddVariable( "AngleNBulkHits", "N Bulk Hits", "", 'I' );


    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    dataloader->AddSignalTree    ( fTree,     signalWeight );
    dataloader->AddBackgroundTree( fTree, backgroundWeight );
    dataloader->AddCut(fMinimalCut);

    // define signal and background cuts
    TCut signalCut = fTruthInFV + TCut("IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)");
    TCut bgCut = fTruthInFV && TCut("IntOrigin==1 && IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)");
    int nMaxSignal = fTree->Draw(">>selectedEntries", signalCut + fMinimalCut, "entrylist")-1;
    int nMaxBg = fTree->Draw(">>selectedEntries", bgCut + fMinimalCut, "entrylist")-1;
    std::cout<<"nMaxSignal: "<<nMaxSignal<<std::endl;
    std::cout<<"nMaxBg: "<<nMaxBg<<std::endl;

    std::string tmva_options = "";
    if(nTrainFrac<0) tmva_options = "nTrain_Signal="+std::to_string(nMaxSignal)+":nTrain_Background="+to_string(nMaxBg);
    else tmva_options = "nTrain_Signal="+std::to_string((int)nTrainFrac*nMaxSignal)+":nTrain_Background="+to_string((int)nTrainFrac*nMaxBg);
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
    if (Use["RandomForest"])  // Adaptive Boost
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "RandomForest",
                             "!H:!V:NTrees=2500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

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