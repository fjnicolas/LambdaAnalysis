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
#include "CutEfficienciesLATeXInterface.C"

#include "../EventHandle/LambdaAnaTree.cpp"

using namespace TMVA::Experimental;


//---------  Main function
void EvaluateBDTAnalysis(TTree *fTree, TTree *fTreeHeader, std::string fWeightFilePath = "dataset/weights/FRAMSSelectionTMVA_BDT.weights.xml",
    double potNorm=1,
    double potNormSignal=1,
    double totalPOTNorm=1){

    TCut fTruthInFV("TruthIsFiducial==1");
    TCut fTruthInAV("TruthIsAV==1");


    // Input ana tree
    LambdaAnaTree fAnaTreeHandle(fTree, true);

    // Create TMVA::Reader object
    TMVA::Reader *fTMVAReader = new TMVA::Reader();

    // Declare variables
    float AngleFRANSScore;
    float NUnOriginsMultGT3;
    float NUnOrigins;
    float CRUMBSScore;
    float AngleDecayContainedDiff;
    float AngleLengthMainTrack;
    float AngleLengthTrack1;
    float AngleLengthTrack2;
    float AngleNHitsMainTrack;
    float AngleNHitsTrack1;
    float AngleNHitsTrack2;
    float AngleMinNHits;
    float NUnassociatedHits;
    float ShowerEnergy;
    float FRANSScorePANDORA;
    float AngleCoveredArea;
    float AngleDirtHits;
    float AngleOpeningAngle;
    float NShowers;
    float NShowerHits;
    float AngleLongestIsMain;
    // Calorimetry variables
    float AngleTwoLinesChi2;
    float AnglePassFit;
    float AnglePassChargeFit;
    float AngleBandOverlap;
    float AngleBandCrossHits;
    float AngleChargeRatioFit;
    float AngleChargeDifferenceFit;
    float AngleChargeRatioIntegral;
    float AngleChargeDifferenceIntegral;
    float AngleChargeRatioAverage;
    float AngleVertexHitIntegralRatio;
    float AngleTrackLengthRatio;
    float AngleResidualRange1RMS;
    float AngleResidualRange2RMS;
    float AngleResidualRangeMinRMS;
    float AngleNVertexHits;
    float AngleNBulkHits;


    


    // Add variables to the reader
    std::map<std::string, bool> fVarLabels = extractLabels(fWeightFilePath);
    std::cout << "Variables to use:" << std::endl;
    for (auto const& x : fVarLabels) {
        std::cout << x.first << std::endl;
    }
    if(fVarLabels["AngleFRANSScore"]==true) fTMVAReader->AddVariable( "AngleFRANSScore", &AngleFRANSScore );
    if(fVarLabels["NUnOriginsMultGT3"]==true) fTMVAReader->AddVariable( "NUnOriginsMultGT3", &NUnOriginsMultGT3 );
    if(fVarLabels["NUnOrigins"]==true) fTMVAReader->AddVariable( "NUnOrigins", &NUnOrigins );
    if(fVarLabels["CRUMBSScore"]==true) fTMVAReader->AddVariable( "CRUMBSScore", &CRUMBSScore );
    if(fVarLabels["AngleDecayContainedDiff"]==true) fTMVAReader->AddVariable( "AngleDecayContainedDiff", &AngleDecayContainedDiff );
    if(fVarLabels["AngleLengthMainTrack"]==true) fTMVAReader->AddVariable( "AngleLengthMainTrack", &AngleLengthMainTrack );
    if(fVarLabels["AngleLengthTrack1"]==true) fTMVAReader->AddVariable( "AngleLengthTrack1", &AngleLengthTrack1 );
    if(fVarLabels["AngleLengthTrack2"]==true) fTMVAReader->AddVariable( "AngleLengthTrack2", &AngleLengthTrack2 );
    if(fVarLabels["AngleNHitsMainTrack"]==true) fTMVAReader->AddVariable( "AngleNHitsMainTrack", &AngleNHitsMainTrack );
    if(fVarLabels["AngleNHitsTrack1"]==true) fTMVAReader->AddVariable( "AngleNHitsTrack1", &AngleNHitsTrack1 );
    if(fVarLabels["AngleNHitsTrack2"]==true) fTMVAReader->AddVariable( "AngleNHitsTrack2", &AngleNHitsTrack2 );
    if(fVarLabels["AngleMinNHits"]==true) fTMVAReader->AddVariable( "AngleMinNHits", &AngleMinNHits );
    if(fVarLabels["NUnassociatedHits"]==true) fTMVAReader->AddVariable( "NUnassociatedHits", &NUnassociatedHits );
    if(fVarLabels["ShowerEnergy"]==true) fTMVAReader->AddVariable( "ShowerEnergy", &ShowerEnergy );
    if(fVarLabels["FRANSScorePANDORA"]==true) fTMVAReader->AddVariable( "FRANSScorePANDORA", &FRANSScorePANDORA );
    if(fVarLabels["AngleCoveredArea"]==true) fTMVAReader->AddVariable( "AngleCoveredArea", &AngleCoveredArea );
    if(fVarLabels["AngleDirtHits"]==true) fTMVAReader->AddVariable( "AngleDirtHits", &AngleDirtHits );
    if(fVarLabels["NShowers"]==true) fTMVAReader->AddVariable( "NShowers", &NShowers );
    if(fVarLabels["NShowerHits"]==true) fTMVAReader->AddVariable( "NShowerHits", &NShowerHits );
    if(fVarLabels["AngleOpeningAngle"]==true) fTMVAReader->AddVariable( "AngleOpeningAngle", &AngleOpeningAngle );
    if(fVarLabels["AngleLongestIsMain"]==true) fTMVAReader->AddVariable( "AngleLongestIsMain", &AngleLongestIsMain );
    // Calo
    if(fVarLabels["AngleTwoLinesChi2"]==true) fTMVAReader->AddVariable( "AngleTwoLinesChi2", &AngleTwoLinesChi2 );
    if(fVarLabels["AnglePassFit"]==true) fTMVAReader->AddVariable( "AnglePassFit", &AnglePassFit );
    if(fVarLabels["AnglePassChargeFit"]==true) fTMVAReader->AddVariable( "AnglePassChargeFit", &AnglePassChargeFit );
    if(fVarLabels["AngleBandOverlap"]==true) fTMVAReader->AddVariable( "AngleBandOverlap", &AngleBandOverlap );
    if(fVarLabels["AngleBandCrossHits"]==true) fTMVAReader->AddVariable( "AngleBandCrossHits", &AngleBandCrossHits );
    if(fVarLabels["AngleChargeRatioFit"]==true) fTMVAReader->AddVariable( "AngleChargeRatioFit", &AngleChargeRatioFit );
    if(fVarLabels["AngleChargeDifferenceFit"]==true) fTMVAReader->AddVariable( "AngleChargeDifferenceFit", &AngleChargeDifferenceFit );
    if(fVarLabels["AngleChargeRatioIntegral"]==true) fTMVAReader->AddVariable( "AngleChargeRatioIntegral", &AngleChargeRatioIntegral );
    if(fVarLabels["AngleChargeDifferenceIntegral"]==true) fTMVAReader->AddVariable( "AngleChargeDifferenceIntegral", &AngleChargeDifferenceIntegral );
    if(fVarLabels["AngleChargeRatioAverage"]==true) fTMVAReader->AddVariable( "AngleChargeRatioAverage", &AngleChargeRatioAverage );
    if(fVarLabels["AngleVertexHitIntegralRatio"]==true) fTMVAReader->AddVariable( "AngleVertexHitIntegralRatio", &AngleVertexHitIntegralRatio );
    if(fVarLabels["AngleTrackLengthRatio"]==true) fTMVAReader->AddVariable( "AngleTrackLengthRatio", &AngleTrackLengthRatio );
    if(fVarLabels["AngleResidualRange1RMS"]==true) fTMVAReader->AddVariable( "AngleResidualRange1RMS", &AngleResidualRange1RMS );
    if(fVarLabels["AngleResidualRange2RMS"]==true) fTMVAReader->AddVariable( "AngleResidualRange2RMS", &AngleResidualRange2RMS );
    if(fVarLabels["AngleResidualRangeMinRMS"]==true) fTMVAReader->AddVariable( "AngleResidualRangeMinRMS", &AngleResidualRangeMinRMS );
    if(fVarLabels["AngleNVertexHits"]==true) fTMVAReader->AddVariable( "AngleNVertexHits", &AngleNVertexHits );
    if(fVarLabels["AngleNBulkHits"]==true) fTMVAReader->AddVariable( "AngleNBulkHits", &AngleNBulkHits );
    
    // Load the BDT
    fTMVAReader->BookMVA( "FRAMS BDT", fWeightFilePath.c_str() );

    // --- Set style
    // bottom margin
    gStyle->SetPadBottomMargin(0.15);
    // stats 0
    gStyle->SetOptStat(0);
    // plot line widths
    gStyle->SetLineWidth(2);
    // graph line width
    gStyle->SetHistLineWidth(2);
    

    // Score histogram
    TH1F *hScoreS = new TH1F("hScoreS", "; Classifier score; # entries", 100, -1, 1);
    TH1F *hScoreBG = new TH1F("hScoreBG", "; Classifier score; # entries", 100, -1, 1);
    TH1F *hScoreCosmic = new TH1F("hScoreCosmic", "; Classifier score; # entries", 100, -1, 1);

    // Save entries that pass cuts
    std::vector<int> passCutEventID;
    std::vector<std::string> passCutEventLabels;

    // loop over the TTree
    std::cout<<"Looping over the TTree"<<std::endl;
    std::cout<<"Number of entries: "<<fAnaTreeHandle.GetEntries()<<std::endl;
    for (int ievt=0; ievt<fAnaTreeHandle.GetEntries();ievt++) {
        
        fAnaTreeHandle.GetEntry(ievt);

        // assign variables
        AngleFRANSScore = fAnaTreeHandle.fAngleFRANSScore;
        FRANSScorePANDORA = fAnaTreeHandle.fFRANSScorePANDORA;
        NUnOriginsMultGT3 = fAnaTreeHandle.fNUnOriginsMultGT3;
        NUnOrigins = fAnaTreeHandle.fNUnOrigins;
        CRUMBSScore = fAnaTreeHandle.fCRUMBSScore;
        AngleDecayContainedDiff = fAnaTreeHandle.fAngleDecayContainedDiff;
        AngleLengthMainTrack = fAnaTreeHandle.fAngleLengthMainTrack;
        AngleLengthTrack1 = fAnaTreeHandle.fAngleLengthTrack1;
        AngleLengthTrack2 = fAnaTreeHandle.fAngleLengthTrack2;
        AngleNHitsMainTrack = fAnaTreeHandle.fAngleNHitsMainTrack;
        AngleNHitsTrack1 = fAnaTreeHandle.fAngleNHitsTrack1;
        AngleNHitsTrack2 = fAnaTreeHandle.fAngleNHitsTrack2;
        NUnassociatedHits = fAnaTreeHandle.fNUnassociatedHits;
        ShowerEnergy = fAnaTreeHandle.fShowerEnergy;
        AngleCoveredArea = fAnaTreeHandle.fAngleCoveredArea;
        AngleDirtHits = fAnaTreeHandle.fAngleDirtHits;
        NShowers = fAnaTreeHandle.fNShowers;
        NShowerHits = fAnaTreeHandle.fNShowerHits;
        AngleLongestIsMain = fAnaTreeHandle.fAngleLongestIsMain;
        // Calo vars
        AngleTwoLinesChi2 = fAnaTreeHandle.fAngleTwoLinesChi2;
        AnglePassFit = fAnaTreeHandle.fAnglePassFit;
        AnglePassChargeFit = fAnaTreeHandle.fAnglePassChargeFit;
        AngleBandOverlap = fAnaTreeHandle.fAngleBandOverlap;
        AngleBandCrossHits = fAnaTreeHandle.fAngleBandCrossHits;
        AngleChargeRatioFit = fAnaTreeHandle.fAngleChargeRatioFit;
        AngleChargeDifferenceFit = fAnaTreeHandle.fAngleChargeDifferenceFit;
        AngleChargeRatioIntegral = fAnaTreeHandle.fAngleChargeRatioIntegral;
        AngleChargeDifferenceIntegral = fAnaTreeHandle.fAngleChargeDifferenceIntegral;
        AngleChargeRatioAverage = fAnaTreeHandle.fAngleChargeRatioAverage;
        AngleVertexHitIntegralRatio = fAnaTreeHandle.fAngleVertexHitIntegralRatio;
        AngleTrackLengthRatio = fAnaTreeHandle.fAngleTrackLengthRatio;
        AngleResidualRange1RMS = fAnaTreeHandle.fAngleResidualRange1RMS;
        AngleResidualRange2RMS = fAnaTreeHandle.fAngleResidualRange2RMS;
        AngleResidualRangeMinRMS = fAnaTreeHandle.fAngleResidualRangeMinRMS;
        AngleNVertexHits = fAnaTreeHandle.fAngleNVertexHits;
        AngleNBulkHits = fAnaTreeHandle.fAngleNBulkHits;

        
        // check active volume
        bool fTruthIsActive = 1;// fAnaTreeHandle.fTruthIsAV;
        
        // retrieve the corresponding MVA output
        double score = fTMVAReader->EvaluateMVA( "FRAMS BDT" );
       

        // check minimal cut
        bool passCut = fAnaTreeHandle.fRecoIsFiducial && fAnaTreeHandle.fNOriginsPairOneTwo>0 && fAnaTreeHandle.fNAngles>=1 && fAnaTreeHandle.fAnglePassChargeFit==1;
        if(!passCut) score = -0.95;

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


        //sif(passCut && isSignal && NShowers>=2) {
        if(passCut && !isSignal && score>0.2){
            std::string label = std::to_string(fAnaTreeHandle.fRunID)+":"+std::to_string(fAnaTreeHandle.fSubrunID);
            std::cout<<"Pass cuts BG Event "<<label<<" score: "<<score<<" FRANSPANDORA:"<<FRANSScorePANDORA<<std::endl;
            passCutEventID.push_back(fAnaTreeHandle.fEventID);
            passCutEventLabels.push_back(label);
        }

      
    } // end of event loop

    std::cout<<" Plotting histograms"<<std::endl;
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
    TLegend *leg = new TLegend(0.7,0.7,0.95,0.9);
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


    // output file for efficiencies
    std::ofstream fEffFile;
    fEffFile.open("ScoreEfficiencies.txt");
    fEffFile<<"Score Signal Background Cosmic"<<std::endl;
    for(int i=0; i<nBins; i++){
        ScoreValues[i] = hScoreS->GetBinCenter(i+1);
        NSignals[i] = hScoreS->Integral(i+1, nBins);
        NBg[i] = hScoreBG->Integral(i+1, nBins);
        NCo[i] = hScoreCosmic->Integral(i+1, nBins);
        EffS[i] = NSignals[i]/hScoreS->Integral();
        EffBg[i] = NBg[i]/hScoreBG->Integral();
        EffCo[i] = NCo[i]/hScoreCosmic->Integral();
        fEffFile<<i<<"-Score="<<ScoreValues[i]<<"  Signal/BG/Cosmic: "<<NSignals[i]<<" "<<NBg[i]<<" "<<NCo[i] << " Eff:"<<EffS[i]<<" "<<EffBg[i]<<" "<<EffCo[i]<<"    "<<potNormSignal*NSignals[i]<<" "<<potNorm*NBg[i]<<" "<<potNorm*NCo[i]<<std::endl;
    }
    fEffFile.close();

    // log scale
    cScore->SetLogy();
    cScore->Update();

    // Notmalize NSignals to POT
    for(int i=0; i<nBins; i++){
        NSignals[i] = NSignals[i]*potNormSignal;
        NBg[i] = NBg[i]*potNorm;
        NCo[i] = NCo[i]*potNorm;
    }

    // Plot NSignals and NBg normalized to POT as a function of the score with TGraph
    TCanvas *cScoreNorm = new TCanvas("cScoreNorm", "cScoreNorm", 800, 600);
    TGraph *gNSignals = new TGraph(nBins, ScoreValues, NSignals);
    TGraph *gNBg = new TGraph(nBins, ScoreValues, NBg);
    TGraph *gNCo = new TGraph(nBins, ScoreValues, NCo);
    TGraph *gEffS = new TGraph(nBins, ScoreValues, EffS);
    TGraph *gEffBg = new TGraph(nBins, ScoreValues, EffBg);
    TGraph *gEffCo = new TGraph(nBins, ScoreValues, EffCo);

    gNSignals->SetLineColor(kBlue);
    gNBg->SetLineColor(kRed);
    gNCo->SetLineColor(kGreen);
    gEffS->SetLineColor(kBlue);
    gEffBg->SetLineColor(kRed);
    gEffCo->SetLineColor(kGreen);

    gNSignals->SetLineWidth(2);
    gNBg->SetLineWidth(2);
    gNCo->SetLineWidth(2);

    // Draw
    gNBg->Draw("AL");
    gNSignals->Draw("L");
    gNCo->Draw("L");


    // set legend
    TLegend *legNorm = new TLegend(0.7,0.7,0.95,0.9);
    legNorm->AddEntry(gNBg,"Background","l");
    legNorm->AddEntry(gNSignals,"Signal","l");
    legNorm->AddEntry(gNCo,"Cosmic+Dirt","l");
    // title total pot, scientific notation
    std::string title = "Total POT: ";
    title += std::to_string(totalPOTNorm/1e20);
    title += "e20";
    legNorm->SetHeader(title.c_str());
    legNorm->Draw();

    // axis labels
    gNBg->GetXaxis()->SetTitle("Classifier score");
    gNBg->GetYaxis()->SetTitle("# events");

    // range user
    gNBg->GetYaxis()->SetRangeUser(0.1, 70);
    gNBg->GetXaxis()->SetRangeUser(0.1, 0.35);

    // Wait
    cScoreNorm->Update();

    // Save canvas to pdf
    cScore->SaveAs("Score.pdf");
    cScoreNorm->SaveAs("ScoreNorm.pdf");
    // set log scale
    cScoreNorm->SetLogy();
    gNBg->GetYaxis()->SetRangeUser(0.1, 1000);
    cScoreNorm->Update();
    cScoreNorm->SaveAs("ScoreNormLog.pdf");

    // cout total entries in histograms
    std::cout<<"Total entries in signal histogram: "<<hScoreS->GetEntries()<<std::endl;
    std::cout<<"Total entries in background histogram: "<<hScoreBG->GetEntries()<<std::endl;


    // Read TreeHeader 
    unsigned int fSubRunId;
    fTreeHeader->SetBranchAddress("SubRunID", &fSubRunId);
    unsigned int fRunId;
    fTreeHeader->SetBranchAddress("RunID", &fRunId);
    std::string *fLArInputFileName = new std::string;
    fTreeHeader->SetBranchAddress("InputFileName", &fLArInputFileName);
    // Create map of run-subrun filename
    std::map<std::string, std::vector<std::string>> subrunFilenameMap;

    std::map<std::string, int> fileCounterMap;
    for(size_t i=0; i<fTreeHeader->GetEntries(); ++i){
        fTreeHeader->GetEntry(i);
            
        //if(fLArInputFileName->find("V0Lambda") != std::string::npos || fLArInputFileName->find("V0Overlay") != std::string::npos) continue;
        
        // check the string includes the substring "Inclusive"
        if(fLArInputFileName->find("Inclusive") == std::string::npos) continue;
    
        std::string runSubrunLabel = std::to_string(fRunId) + ":" + std::to_string(fSubRunId);

        // if the subrun is not in the map, add it
        if(subrunFilenameMap.find(runSubrunLabel) == subrunFilenameMap.end()){
            subrunFilenameMap[runSubrunLabel] = std::vector<std::string>();
        }

        // add the filename to the vector
        subrunFilenameMap[runSubrunLabel].push_back(*fLArInputFileName);

        // if the file is not in the map, add it
        if(fileCounterMap.find(*fLArInputFileName) == fileCounterMap.end()){
            fileCounterMap[*fLArInputFileName] = 1;
        }
        else{
            fileCounterMap[*fLArInputFileName]++;
        } 
    }
    
    std::cout<<"Number of files: "<<fileCounterMap.size()<<std::endl;
    // loop over pass cuts and print run, subrun and event in file
    std::ofstream fPassCutFile;
    fPassCutFile.open("handScanBDT.txt");
    for(int i=0; i<passCutEventLabels.size(); i++){
        std::string eventLabel = passCutEventLabels[i]; 

        std::cout<<"Event label: "<<eventLabel<<std::endl;       
        
        std::vector<std::string> filenames = subrunFilenameMap[eventLabel];
        std::cout<<"Number of files: "<<filenames.size()<<std::endl;
        for(auto& filename : filenames){
            std::cout<<"Filename: "<<filename<<std::endl;
            fPassCutFile<<eventLabel<<":"<<passCutEventID[i]<<" "<<filename<<std::endl;
        }
        
    }


}

void MacroEvaluateBDT(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, 3.3e20, potScaling, potScalingSignal);
    std::cout<<"POT scaling: "<<potScaling<<" POT scaling signal: "<<potScalingSignal<<std::endl;
    
    EvaluateBDTAnalysis(fTree, fTreeHeader, "dataset/weights/FRAMSSelectionTMVA_BDT.weights.xml", potScaling, potScalingSignal, 3.3e20);
}