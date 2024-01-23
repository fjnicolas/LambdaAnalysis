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

#include "../SObjects/LambdaAnaTree.cpp"

using namespace TMVA::Experimental;


void PrintHandScanEvents(TTree *treeHeader, std::string outputFilename, std::vector<int> highestScoreEventID, std::vector<std::string> highestScoreEventLabels, std::vector<double> highestScoreEventScores, std::string fileKeyLabel=""){
    
    // Read TreeHeader 
    unsigned int fSubRunId;
    treeHeader->SetBranchAddress("SubRunID", &fSubRunId);
    unsigned int fRunId;
    treeHeader->SetBranchAddress("RunID", &fRunId);
    std::string *fLArInputFileName = new std::string;
    treeHeader->SetBranchAddress("InputFileName", &fLArInputFileName);
    
    // Create map of run-subrun filename
    std::map<std::string, std::vector<std::string>> subrunFilenameMap;
    std::map<std::string, int> fileCounterMap;
    for(size_t i=0; i<treeHeader->GetEntries(); ++i){
        
        treeHeader->GetEntry(i);
            
        // check the string includes the substring "Inclusive"
        if(fileKeyLabel!=""){
            if(fLArInputFileName->find(fileKeyLabel) == std::string::npos) continue;
        }
    
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

    // output file    
    std::ofstream fPassCutFile;
    fPassCutFile.open(outputFilename);
    
    // order highestScore vectors by score value
    std::vector<int> highestScoreEventIDOrdered;
    std::vector<std::string> highestScoreEventLabelsOrdered;
    std::vector<double> highestScoreEventScoresOrdered;
    while(highestScoreEventID.size()>0){
        double highestScore = highestScoreEventScores[0];
        int highestScoreIndex = 0;
        for(int i=1; i<highestScoreEventScores.size(); i++){
            if(highestScoreEventScores[i]>highestScore){
                highestScore = highestScoreEventScores[i];
                highestScoreIndex = i;
            }
        }
        highestScoreEventIDOrdered.push_back(highestScoreEventID[highestScoreIndex]);
        highestScoreEventLabelsOrdered.push_back(highestScoreEventLabels[highestScoreIndex]);
        highestScoreEventScoresOrdered.push_back(highestScoreEventScores[highestScoreIndex]);
        highestScoreEventID.erase(highestScoreEventID.begin()+highestScoreIndex);
        highestScoreEventLabels.erase(highestScoreEventLabels.begin()+highestScoreIndex);
        highestScoreEventScores.erase(highestScoreEventScores.begin()+highestScoreIndex);
    }

    // print highest score events
    for(int i=0; i<highestScoreEventIDOrdered.size(); i++){
        std::cout<<"\nEvent "<<highestScoreEventIDOrdered[i]<<" "<<highestScoreEventLabelsOrdered[i]<<" score: "<<highestScoreEventScoresOrdered[i]<<std::endl;
        std::string eventLabel = highestScoreEventLabelsOrdered[i];
        std::vector<std::string> filenames = subrunFilenameMap[eventLabel];
        std::cout<<"Number of files: "<<filenames.size()<<std::endl;
        for(auto& filename : filenames){
            std::cout<<"Filename: "<<filename<<std::endl;
            fPassCutFile<<eventLabel<<":"<<highestScoreEventIDOrdered[i]<<" "<<filename<<std::endl;
        }
    }

    return;
}

std::map<std::string, int> fColorMap = {
    {"Signal", kBlue},
    {"Background", kRed},
    {"Cosmic", kGreen}
};



void PlotHistograms(std::map<std::string, TH1F*>& histogramMap, const std::string& outputDir, double potNormSignal, double potNorm, double totalPOTNorm) {
    std::cout << " Plotting histograms" << std::endl;

    // Parameters
    int fTargetNEvents = 20;

    // Create TCanvas
    TCanvas *cScore = new TCanvas("cScore", "cScore", 800, 600);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 0.5, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0.5, 1, 1);
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0., 0.5, 0.5);
    TPad *pad4 = new TPad("pad4", "pad4", 0.5, 0., 1, 0.5);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();

    // General gStyle
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);

    // Output file for efficiencies
    std::ofstream fEffFile;
    fEffFile.open(outputDir+"Efficiencies.txt");
    fEffFile<<"Score   Signal   Background   Cosmic"<<std::endl;

    // --- Score distribution ---
    pad1->cd();
    pad1->SetLogy();
    // Iterate through the histogram map
    for (auto const& hist : histogramMap) {
        // if first histogram, draw it
        if (hist == *histogramMap.begin()) {
            hist.second->Draw();
        }
        else {
            hist.second->Draw("same");
        }
        
        hist.second->SetLineColor(fColorMap[hist.first]);
    }

    // Set legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (auto const& hist : histogramMap)
        leg->AddEntry(hist.second, hist.first.c_str(), "l");
    leg->Draw("same");

    // --- Efficiencies ---
    int nBins = histogramMap["Signal"]->GetNbinsX();
    double ScoreValues[nBins];
    double NSignals[nBins];
    double NBg[nBins];
    double NCo[nBins];
    double EffS[nBins];
    double EffBg[nBins];
    double EffCo[nBins];

    // Fill
    for(int i=0; i<nBins; i++){
        ScoreValues[i] = histogramMap["Signal"]->GetBinCenter(i+1);
        
        double n;

        // signal
        n = histogramMap["Signal"]->Integral(i+1, nBins);
        NSignals[i] = n*potNormSignal;
        EffS[i] = n/histogramMap["Signal"]->Integral();

        // background
        n = histogramMap["Background"]->Integral(i+1, nBins);
        NBg[i] = n*potNorm;
        EffBg[i] = n/histogramMap["Background"]->Integral();

        // cosmic
        n = histogramMap["Cosmic"]->Integral(i+1, nBins);
        NCo[i] = n*potNorm;
        EffCo[i] = n/histogramMap["Cosmic"]->Integral();

        fEffFile<<i<<"-Score="<<ScoreValues[i]<<"  Signal/BG/Cosmic: "<<NSignals[i]<<" "<<NBg[i]<<" "<<NCo[i] << " Eff:"<<EffS[i]<<" "<<EffBg[i]<<" "<<EffCo[i]<<"    "<<potNormSignal*NSignals[i]<<" "<<potNorm*NBg[i]<<" "<<potNorm*NCo[i]<<std::endl;
    }

    pad2->cd();
    TGraph *gEffS = new TGraph(nBins, ScoreValues, EffS);
    TGraph *gEffBg = new TGraph(nBins, ScoreValues, EffBg);
    TGraph *gEffCo = new TGraph(nBins, ScoreValues, EffCo);
    gEffS->SetLineColor(fColorMap["Signal"]);
    gEffBg->SetLineColor(fColorMap["Background"]);
    gEffCo->SetLineColor(fColorMap["Cosmic"]);
    gEffS->Draw("AL");
    gEffBg->Draw("L");
    gEffCo->Draw("L");
    // axis labels
    gEffS->GetXaxis()->SetTitle("Classifier score");
    gEffS->GetYaxis()->SetTitle("efficiency");
    

    // title total pot, scientific notation
    std::string title = "POT: "+std::to_string(totalPOTNorm/1e20)+"e20";
    leg->SetHeader(title.c_str());
    leg->Draw("same");


    // --- # events POT normalized ---
    pad3->cd();
    TGraph *gNSignals = new TGraph(nBins, ScoreValues, NSignals);
    TGraph *gNBg = new TGraph(nBins, ScoreValues, NBg);
    TGraph *gNCo = new TGraph(nBins, ScoreValues, NCo);
    gNSignals->SetLineColor(fColorMap["Signal"]);
    gNBg->SetLineColor(fColorMap["Background"]);
    gNCo->SetLineColor(fColorMap["Cosmic"]);
    gNSignals->Draw("AL");
    gNBg->Draw("L");
    gNCo->Draw("L");
    leg->Draw("same");

    // axis labels
    gNSignals->GetXaxis()->SetTitle("Classifier score");
    gNSignals->GetYaxis()->SetTitle("# events");

    // Draw horizontal line at y=20
    TLine *l20 = new TLine(ScoreValues[0], fTargetNEvents, ScoreValues[nBins-1], fTargetNEvents);
    l20->SetLineColor(kBlack);
    l20->SetLineStyle(2);
    l20->Draw("same");

    // draw copy of pad3 in pad4, with log scale

    pad4->cd();
    pad4->SetLogy();
    gNSignals->Draw("AL");
    gNBg->Draw("L");
    gNCo->Draw("L");
    leg->Draw("same");
    l20->Draw("same");

    // Save canvas to pdf
    cScore->Update();
    cScore->SaveAs((outputDir + "BDTResponse.pdf").c_str());
    //pad1->SaveAs((outputDir + "ScorePad1.pdf").c_str());
    
    // Same plots, with zoom
    // range users
    gNSignals->GetXaxis()->SetRangeUser(0.1, 0.6);
    gNSignals->GetYaxis()->SetRangeUser(0.1, 500);
    gEffS->GetXaxis()->SetRangeUser(0.1, 0.6);
    gEffS->GetYaxis()->SetRangeUser(0., 0.4);
    // change l20 range
    l20->SetX1(0.1);
    l20->SetX2(0.6);
    cScore->Update();
    cScore->SaveAs((outputDir + "BDTResponseZoom.pdf").c_str());

    

    fEffFile.close();
}


//---------  Main function
void RunEvaluateBDTAnalysis(TTree *fTree, TTree *fTreeHeader, std::string fWeightFilePath = "dataset/weights/FRANSSelectionTMVA_BDT.weights.xml",
    double potNorm=1, double potNormSignal=1, double totalPOTNorm=1)
{

    TCut fTruthInFV("TruthIsFiducial==1");
    TCut fTruthInAV("TruthIsAV==1");

    int fNHighestScoreEvents = 5;

    // Output directory
    std::string fOutputDir = "OutputBDTAna/";
    gSystem->mkdir(fOutputDir.c_str(), kTRUE);

    // Input ana tree
    LambdaAnaTree fAnaTreeHandle(fTree, true);

    // Create TMVA::Reader object
    TMVA::Reader *fTMVAReader = new TMVA::Reader();

    // Declare variables
    float AngleFRANSScore;
    float NUnOriginsMultGT3;
    float NUnOrigins;
    float NOrigins;
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
    if(fVarLabels["NOrigins"]==true) fTMVAReader->AddVariable( "NOrigins", &NOrigins );
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
    if(fVarLabels["FRANSScorePANDORA"]==true) fTMVAReader->AddVariable( "FRANSScorePANDORA", &FRANSScorePANDORA );
    if(fVarLabels["AngleCoveredArea"]==true) fTMVAReader->AddVariable( "AngleCoveredArea", &AngleCoveredArea );
    if(fVarLabels["AngleDirtHits"]==true) fTMVAReader->AddVariable( "AngleDirtHits", &AngleDirtHits );
    if(fVarLabels["NShowers"]==true) fTMVAReader->AddVariable( "NShowers", &NShowers );
    if(fVarLabels["NShowerHits"]==true) fTMVAReader->AddVariable( "NShowerHits", &NShowerHits );
    if(fVarLabels["ShowerEnergy"]==true) fTMVAReader->AddVariable( "ShowerEnergy", &ShowerEnergy );
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
    fTMVAReader->BookMVA( "FRANS BDT", fWeightFilePath.c_str() );

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
    // Make previous a map
    std::map<std::string, TH1F*> hScoreMap{
        {"Signal", new TH1F("hScoreS", "; Classifier score; # entries", 100, -1, 1)},
        {"Background", new TH1F("hScoreBG", "; Classifier score; # entries", 100, -1, 1)},
        {"Cosmic", new TH1F("hScoreCosmic", "; Classifier score; # entries", 100, -1, 1)}
    };


    // Store highest BG score events
    std::vector<int> highestScoreEventID;
    std::vector<std::string> highestScoreEventLabels;
    std::vector<double> highestScoreEventScores;
    
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
        NOrigins = fAnaTreeHandle.fNOrigins;
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
        double score = fTMVAReader->EvaluateMVA( "FRANS BDT" );
       
        // check minimal cut
        bool minimalCut = fAnaTreeHandle.fRecoIsFiducial && fAnaTreeHandle.fNAngles>=1;
        bool passFitCuts = fAnaTreeHandle.fAnglePassChargeFit==1;//&& fAnaTreeHandle.fAnglePassFit==1;
        bool selCuts = fAnaTreeHandle.fFRANSScorePANDORA>0.2 && fAnaTreeHandle.fAngleDecayContainedDiff<1 && fAnaTreeHandle.fNUnOrigins<=0;
        bool passCut;

        passCut = minimalCut && passFitCuts;
        if(!passCut) score = -0.95;

        bool isSignal = fAnaTreeHandle.fIntOrigin==1 && fAnaTreeHandle.fIntDirt==0 && (fAnaTreeHandle.fIntNLambda>0 && fAnaTreeHandle.fIntMode==0 && std::abs(fAnaTreeHandle.fIntNuPDG)!=12);
        bool isNuBG = fAnaTreeHandle.fIntOrigin==1 && fAnaTreeHandle.fIntDirt==0 && !(fAnaTreeHandle.fIntNLambda>0 && fAnaTreeHandle.fIntMode==0 && std::abs(fAnaTreeHandle.fIntNuPDG)!=12);
        
        // Fill histogram for signals
        if(fTruthIsActive){
            if(isSignal)
                hScoreMap["Signal"]->Fill(score);
            else if(isNuBG)
                hScoreMap["Background"]->Fill(score);
        }
        else if(!fTruthIsActive && (fAnaTreeHandle.fIntOrigin==1 || fAnaTreeHandle.fIntDirt==1) ){
            hScoreMap["Cosmic"]->Fill(score);
        }

        // Fill highest score events
        if(passCut && !isSignal){
            //fAnaTreeHandle.PrintEventInfo();
            if(highestScoreEventID.size()<fNHighestScoreEvents){
                highestScoreEventID.push_back(fAnaTreeHandle.fEventID);
                std::string label = std::to_string(fAnaTreeHandle.fRunID)+":"+std::to_string(fAnaTreeHandle.fSubrunID);
                highestScoreEventLabels.push_back(label);
                highestScoreEventScores.push_back(score);
            }
            else{
                // check if score is higher than the lowest score in the vector
                double lowestScore = highestScoreEventScores[0];
                int lowestScoreIndex = 0;
                for(int i=1; i<fNHighestScoreEvents; i++){
                    if(highestScoreEventScores[i]<lowestScore){
                        lowestScore = highestScoreEventScores[i];
                        lowestScoreIndex = i;
                    }
                }
                // if score is higher than the lowest score, replace it
                if(score>lowestScore){
                    highestScoreEventID[lowestScoreIndex] = fAnaTreeHandle.fEventID;
                    std::string label = std::to_string(fAnaTreeHandle.fRunID)+":"+std::to_string(fAnaTreeHandle.fSubrunID);
                    highestScoreEventLabels[lowestScoreIndex] = label;
                    highestScoreEventScores[lowestScoreIndex] = score;
                }
            }
        }

      
    } // end of event loop

    PlotHistograms(hScoreMap, fOutputDir, potNormSignal, potNorm, totalPOTNorm);



    // --- Save BG events with highest score
    PrintHandScanEvents(fTreeHeader, fOutputDir+"HandScanEvents.txt", highestScoreEventID, highestScoreEventLabels, highestScoreEventScores, "Inclusive");
    
}