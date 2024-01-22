#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "LambdaBDTAnalysis.C"
#include "LambdaAnalysis.C"


void LoadAnalysisMacros(){
    std::cout<<" Loaded ana macros\n";
    return;
}

void MacroEvaluateBDTAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTrees
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, 3.3e20, potScaling, potScalingSignal);
    std::cout<<"POT scaling: "<<potScaling<<" POT scaling signal: "<<potScalingSignal<<std::endl;
    
    RunEvaluateBDTAnalysis(fTree, fTreeHeader, "dataset/weights/FRANSSelectionTMVA_BDT.weights.xml", potScaling, potScalingSignal, 3.3e20);
}


void MacroBDTAnalysis(std::string fInputFileName="", double nTrainFrac = -1, bool useBatchMode=false, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree"){
    RunLambdaBDTAnalysis(fInputFileName, nTrainFrac, useBatchMode, fTreeDirName, fTreeName);
}


void MacroRunAndEvaluateBDTAnalysis(std::string fInputFileName="", std::string fInputFileNameTest="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree"){

    // First make the training
    RunLambdaBDTAnalysis(fInputFileName, -1, batchMode, fTreeDirName, fTreeName);

    // Then evaluate the BDT
    //--------- Input TTrees
    TFile *fFile = new TFile(fInputFileNameTest.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, 3.3e20, potScaling, potScalingSignal);
    std::cout<<"POT scaling: "<<potScaling<<" POT scaling signal: "<<potScalingSignal<<std::endl;
    
    RunEvaluateBDTAnalysis(fTree, fTreeHeader, "dataset/weights/FRANSSelectionTMVA_BDT.weights.xml", potScaling, potScalingSignal, 3.3e20);


}