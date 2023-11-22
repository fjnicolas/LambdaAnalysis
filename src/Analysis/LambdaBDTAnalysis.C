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

#include "../EventHandle/LambdaAnaTree.cpp"

#include "CutEfficienciesDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

using namespace TMVA::Experimental;

//---------  Main function
void LambdaBDTAnalysis(std::string fInputFileName="", bool useBatchMode=false, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    // Number of events for training
    int nTrain = 1000;

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


    //Add variables for MVA
    //dataloader->AddVariable( "RecoIsFiducial", "RecoIsFiducial", "", 'I' );
    //dataloader->AddVariable( "NOriginsPairOneTwo", "NOriginsPairOneTwo", "", 'I' );
    //dataloader->AddVariable( "NAngles", "NAngles", "", 'I' );
    

    dataloader->AddVariable( "AngleFRANSScore", "AngleFRANSScore", "", 'D' );
    dataloader->AddVariable( "NOriginsMultGT3", "N^{1}_{C}", "", 'I' );
    dataloader->AddVariable( "NOrigins", "N^{2}_{C}", "", 'I' );

    dataloader->AddVariable( "FRANSScorePANDORA", "FRANS PANDORA", "", 'D' );

    dataloader->AddVariable( "CRUMBSScore", "CRUMBS", "", 'D' );
    
    dataloader->AddVariable( "AngleGap", "Gap", "", 'D' );
    dataloader->AddVariable( "AngleDecayContainedDiff", "#alpha", "", 'D' );
    dataloader->AddVariable( "AngleNHits", "V # hits", "", 'I' );

    dataloader->AddVariable( "ShowerEnergy", "Shower Energy", "", 'D' );
    dataloader->AddVariable( "NShwTh100", "# showers > 100 MeV", "", 'D' );    
    

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

    std::string tmva_options = "";// "nTrain_Signal="+std::to_string(nTrain)+":nTrain_Background="+to_string(nTrain);
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




//---------  Main function
void LambdaEvaluateBDTAnalysis(std::string fInputFileName="", bool useBatchMode=false, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree"){

    TCut fTruthInFV("TruthIsFiducial==1");
    TCut fTruthInAV("abs(NuvX)<200 && abs(NuvY)<200 && NuvZ>0 && NuvZ<500");

    // Batch mode
    useBatchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    // Input file
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    // Input ana tree
    LambdaAnaTree fAnaTreeHandle(fTree, true);

    fAnaTreeHandle.PrintTree();
    // Create TMVA::Reader object
    TMVA::Reader *fTMVAReader = new TMVA::Reader();

    float AngleFRANSScore;
    fTMVAReader->AddVariable( "AngleFRANSScore", &AngleFRANSScore );
    float NOriginsMultGT3;
    fTMVAReader->AddVariable( "NOriginsMultGT3", &NOriginsMultGT3 );
    float NOrigins;
    fTMVAReader->AddVariable( "NOrigins", &NOrigins );
    float FRANSScorePANDORA;
    fTMVAReader->AddVariable( "FRANSScorePANDORA", &FRANSScorePANDORA );
    float AngleGap;
    float CRUMBSScore;
    fTMVAReader->AddVariable( "CRUMBSScore", &CRUMBSScore );
    fTMVAReader->AddVariable( "AngleGap", &AngleGap );
    float AngleDecayContainedDiff;
    fTMVAReader->AddVariable( "AngleDecayContainedDiff", &AngleDecayContainedDiff );
    float AngleNHits;
    fTMVAReader->AddVariable( "AngleNHits", &AngleNHits );
    float ShowerEnergy;
    fTMVAReader->AddVariable( "ShowerEnergy", &ShowerEnergy );
    float NShwTh100;
    fTMVAReader->AddVariable( "NShwTh100", &NShwTh100 );
    

    // Load the BDT
    fTMVAReader->BookMVA( "FRAMS BDT", "dataset/weights/FRAMSSelectionTMVA_BDT.weights.xml" );

    // Score histogram
    TH1F *hScoreS = new TH1F("hScoreS", "hScore", 100, -1, 1);
    TH1F *hScoreBG = new TH1F("hScoreBG", "hScore", 100, -1, 1);
    TH1F *hScoreCosmic = new TH1F("hScoreCosmic", "hScoreCosmic", 100, -1, 1);



    // loop over the TTree
    for (Long64_t ievt=0; ievt<fAnaTreeHandle.GetEntries();ievt++) {
        fAnaTreeHandle.GetEntry(ievt);

        
        // assign variables
        AngleFRANSScore = fAnaTreeHandle.fAngleFRANSScore;
        NOriginsMultGT3 = fAnaTreeHandle.fNOriginsMultGT3;
        NOrigins = fAnaTreeHandle.fNOrigins;
        FRANSScorePANDORA = fAnaTreeHandle.fFRANSScorePANDORA;
        AngleGap = fAnaTreeHandle.fAngleGap;
        AngleDecayContainedDiff = fAnaTreeHandle.fAngleDecayContainedDiff;
        AngleNHits = fAnaTreeHandle.fAngleNHits;
        ShowerEnergy = fAnaTreeHandle.fShowerEnergy;
        NShwTh100 = fAnaTreeHandle.fNShwTh100;
        CRUMBSScore = fAnaTreeHandle.fCRUMBSScore;


        // check active volume
        bool fTruthIsActive = std::abs(fAnaTreeHandle.fNuvX)<200 && std::abs(fAnaTreeHandle.fNuvY)<200 && fAnaTreeHandle.fNuvZ>0 && fAnaTreeHandle.fNuvZ<500;
    

        // retrieve the corresponding MVA output
        double score = fTMVAReader->EvaluateMVA( "FRAMS BDT" );

        // check minimal cut
        bool passCut = fAnaTreeHandle.fRecoIsFiducial && fAnaTreeHandle.fNOriginsPairOneTwo>0 && fAnaTreeHandle.fNAngles>=1;
        if(!passCut) score = -0.95;

        //std::cout<<ievt<<" "<<score<<std::endl;

        bool isSignal = fAnaTreeHandle.fIntOrigin==1 && fAnaTreeHandle.fIntDirt==0 && (fAnaTreeHandle.fIntNLambda>0 && fAnaTreeHandle.fIntMode==0 && std::abs(fAnaTreeHandle.fIntNuPDG)!=12);
        bool isNuBG = fAnaTreeHandle.fIntOrigin==1 && fAnaTreeHandle.fIntDirt==0 && !(fAnaTreeHandle.fIntNLambda>0 && fAnaTreeHandle.fIntMode==0 && std::abs(fAnaTreeHandle.fIntNuPDG)!=12);
        // Fill histogram for signals
        if(fTruthIsActive){
            if(isSignal)
                hScoreS->Fill(score);
            else if(isNuBG)
                hScoreBG->Fill(score);
        }
        else if(!fTruthIsActive && (fAnaTreeHandle.fIntOrigin==1 || fAnaTreeHandle.fIntDirt==1) ){
            hScoreCosmic->Fill(score);
        }


    } // end of event loop


    // Draw histogram
    TCanvas *cScore = new TCanvas("cScore", "cScore", 800, 600);
    hScoreBG->Draw();
    hScoreS->Draw("same");
    hScoreCosmic->Draw("same");
    // set colors
    hScoreBG->SetLineColor(kRed);
    hScoreS->SetLineColor(kBlue);
    hScoreCosmic->SetLineColor(kGreen);
    // set legend
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(hScoreBG,"Background","l");
    leg->AddEntry(hScoreS,"Signal","l");
    leg->AddEntry(hScoreCosmic,"Cosmic","l");
    leg->Draw();

    // scan score acis, get NSignals and NBg
    int nBins = hScoreS->GetNbinsX();
    double ScoreValues[nBins];
    double NSignals[nBins];
    double NBg[nBins];
    double NCo[nBins];
    double EffS[nBins];
    double EffBg[nBins];
    double EffCo[nBins];

    for(int i=0; i<nBins; i++){
        ScoreValues[i] = hScoreS->GetBinCenter(i+1);
        NSignals[i] = hScoreS->Integral(i+1, nBins);
        NBg[i] = hScoreBG->Integral(i+1, nBins);
        NCo[i] = hScoreCosmic->Integral(i+1, nBins);
        EffS[i] = NSignals[i]/hScoreS->Integral();
        EffBg[i] = NBg[i]/hScoreBG->Integral();
        EffCo[i] = NCo[i]/hScoreCosmic->Integral();
        std::cout<<i<<"-Score="<<ScoreValues[i]<<"Signal/BG: "<<NSignals[i]<<" "<<NBg[i]<<" "<<NCo[i] << " Eff:"<<EffS[i]<<" "<<EffBg[i]<<" "<<EffCo[i]<<std::endl;
    }

    // log scale
    cScore->SetLogy();

    cScore->Update();
    cScore->WaitPrimitive();

    // cout total entries in histograms
    std::cout<<"Total entries in signal histogram: "<<hScoreS->GetEntries()<<std::endl;
    std::cout<<"Total entries in background histogram: "<<hScoreBG->GetEntries()<<std::endl;

    cScore->SaveAs("OutputPlots/Score.pdf");
}