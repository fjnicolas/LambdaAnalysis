#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

//--------- Settings ---------
//--------- Scale POT
bool fScaleHistogramsToPOT = 0;
//--------- Plot all the variables
bool fPlotAllVars = 0;
//-------- POT normalization
double fPOTTotalNorm = 3.3e20;
//---------  LATeX output file
std::string fOutputFileName = "CutEfficiencies";
std::string fOutputFileNameNormalized = "CutEfficienciesNormalized";

//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = {
    {fTruthInAV+"IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Signal", true, "Signal"}
    ,{fTruthInAV+"IntOrigin==1 &&  IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "BG  #nu", false, "BG BNB"}
};

//---------  Phase space cuts
std::vector<PlotDef> fPhaseSpaceDefs = {};//psLambdaKinematics;

//---------  Cuts
std::vector<PlotDef> fCutDefs = cutDefsPID;


//---------  Load function
void LambdaAnalysis(){
    std::cout<<"LambdaAnalysis function loaded"<<std::endl;
    return;
}

//---------  Main function
void RunLambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTreePost")
{

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName+"/PhaseSpace").c_str());
    gSystem->Exec(("mkdir "+fOutputDirName+"/OtherDistributions").c_str());

    //--------- Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);


    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //----------------- POT normalization
    double potScalingBg = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, fPOTTotalNorm, potScalingBg, potScalingSignal, (fTreeDirName+"pottree").c_str());
    if(!fScaleHistogramsToPOT){
        potScalingBg = 1;
        potScalingSignal = 1;
    }
    std::cout<<"POT scaling: "<<potScalingBg<<" POT scaling signal: "<<potScalingSignal<<std::endl;

    //--------- Matrix to store the number of events
    std::vector< std::vector<int> > nEventsMatrix;
    std::vector<PlotDef> cutDefsForTable;

    //--------- Get vector with the cuts to accumulate
    std::vector<PlotDef> cutsToAccumulate;
    if(fPlotAllVars){
        for (size_t i = 0; i < fCutDefs.size(); ++i) {
            if(fCutDefs[i].GetAccumulateCut()){
                cutsToAccumulate.push_back(fCutDefs[i]);
            }
        }
    }

    //--------- Loop over the cuts
    std::vector<AnaPlot> anaPlots;
    TCut previousCut("");
    for (size_t i = 0; i < fCutDefs.size(); ++i) {
        
        // Get the current cut
        TCut currentCut = TCut(fCutDefs[i].GetCut());

        // Create the plot handle
        AnaPlot anaPlot(i, fCutDefs[i], sampleDefs, fPhaseSpaceDefs, cutsToAccumulate);

        // Draw the histograms
        anaPlot.DrawHistograms(fTree, previousCut, 0, potScalingSignal, potScalingBg);
        anaPlots.push_back(anaPlot);

        // Id the cut is to be accumulated
        if(fCutDefs[i].GetAccumulateCut()){
            // Store the previous cut
            previousCut = previousCut && currentCut;

            // Draw the histograms again
            anaPlot.DrawHistograms(fTree, previousCut, 1, potScalingSignal, potScalingBg);

            // Store the numbers for final efficiency table
            cutDefsForTable.push_back(fCutDefs[i]);
            nEventsMatrix.push_back({});
            size_t cutIndex = nEventsMatrix.size()-1;
            std::map<std::string, int> nEventsMap = anaPlot.GetCountsV();
            for(const auto& sample : sampleDefs){
                nEventsMatrix[cutIndex].push_back(nEventsMap[sample.GetLabelS()]);
            }

        }
    
    }

    // --- Loop over the samples and set NEvents
    for(auto& sample : sampleDefs){
        int nEvents = fTree->Draw( "", TCut(sample.GetVar())+fCounterCut, "goff");
        sample.SetNEvents(nEvents);
        std::cout<<"Sample: "<<sample.GetLabelS()<<" NEvents: "<<nEvents<<std::endl;
    }
    
    //--------- Create the LaTeX table
    GenerateAndCompileTeXTable(cutDefsForTable, sampleDefs, nEventsMatrix, fOutputFileName, "Cut efficiencies", fOutputDirName);

    //--------- Create the LaTeX table (POT normalized)
    GenerateAndCompileTeXTable(cutDefsForTable, sampleDefs, nEventsMatrix, fOutputFileNameNormalized, "Cut efficiencies", fOutputDirName, potScalingBg, potScalingSignal, fPOTTotalNorm);
   
    //--------- Output hand scans
    CreateHandScanList(fTree, fTreeHeader, previousCut, sampleDefs);

    return;
}



void RunCutLoopAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputCutsLoop";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    
    //---------- Set of cuts to loop over
    std::vector<PlotDef> loopCuts = cutDefsLoop;   
    // Minimal cut
    PlotDef minimalCut = minimalCutForLoop;
    
    // Loop over all possible orders using std::next_permutation
    int permutationId = 0;
    do {
        
        // Start Cut
        TCut currentCut("");

        // Vector of cuts including the minimal cut
        std::vector<PlotDef> loopCutsWithMinimal = loopCuts;
        loopCutsWithMinimal.insert(loopCutsWithMinimal.begin(), minimalCut);
        loopCutsWithMinimal.insert(loopCutsWithMinimal.begin(), startCutForLoop);

        //--------- Matrix to store the number of events
        std::vector< std::vector<int> > nEventsMatrix;
        nEventsMatrix.resize(loopCutsWithMinimal.size());
        for (size_t i = 0; i < loopCutsWithMinimal.size(); ++i) {
            nEventsMatrix[i].resize(sampleDefs.size());
        }

        for (const auto& cut : loopCutsWithMinimal) {
            std::cout << "Adding cut: " << cut.GetCut() << "\n";
            currentCut+=TCut(cut.GetCut());

            for(const auto& sample : sampleDefs){
                TCut currentSampleCut = TCut(sample.GetVar());

                int n = fTree->Draw( "", currentCut+currentSampleCut, "goff");
                nEventsMatrix[&cut - &loopCutsWithMinimal[0]][&sample - &sampleDefs[0]] = n;
            }
            
        }


        //--------- Create the LaTeX table
        GenerateAndCompileTeXTable(loopCutsWithMinimal, sampleDefs, nEventsMatrix, fOutputFileName+std::to_string(permutationId), "Cut efficiencies "+std::to_string(permutationId), fOutputDirName);

        permutationId++;


    } while (std::next_permutation(loopCuts.begin(), loopCuts.end()));

    return 0;
}