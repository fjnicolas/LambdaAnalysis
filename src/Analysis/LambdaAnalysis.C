#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

std::string fTruthInFV = "TruthIsFiducial==1 &&";
std::string fTruthInAV = "abs(NuvX)<200 && abs(NuvY)<200 && NuvZ>0 && NuvZ<500 && ";

//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = {
    {fTruthInAV+"IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Signal", true}
    ,{fTruthInAV+"IntOrigin==1 &&  IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Background  Nu", false}
    //,{"IntOrigin==2 || IntDirt==1", "Dirt+Cosmic", false}
    //,{"IntOrigin==1 && IntDirt==1", "Dirt", false}
    //,{"IntOrigin==2", "Cosmic", false}
    /*,{fTruthInFV+"IntNLambda==0 && IntMode==0 && abs(IntNuPDG)!=12", "QE", false}
    ,{fTruthInFV+"IntMode==1 && abs(IntNuPDG)!=12", "RES", false}
    ,{fTruthInFV+"IntMode==2 && abs(IntNuPDG)!=12", "DIS", false}
    ,{fTruthInFV+"(IntMode==3 || IntMode==10) && abs(IntNuPDG)!=12", "COH and MEC", false}
    ,{fTruthInFV+"abs(IntNuPDG)==12", "NuE", false} */
    //,{fTruthInFV+"IntNLambda==0 && IntMode==0 && abs(IntNuPDG)!=12", "QE", false}
    //,{fTruthInFV+"IntMode==1 && abs(IntNuPDG)!=12", "RES", false}
};

//-------- POT normalization
double fPOTTotalNorm = 3.3e20;

//---------  Phase space cuts
std::vector<PlotDef> phaseSpaceDefs = {};// phaseSpaceVars;

//---------  LATeX output file
std::string fOutputFileName = "CutEfficiencies";
std::string fOutputFileNameNormalized = "CutEfficienciesNormalized";

TCut fCounterCut = "SliceID==0";

std::vector<PlotDef> cutDefs = cutDefsTalk2;


//---------  Load function
void LambdaAnalysis(){
    std::cout<<"LambdaAnalysis function loaded"<<std::endl;
    return;
}

//---------  Main function
void RunLambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName+"/PhaseSpace").c_str());

    //--------- Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Scale POT
    bool fScaleHistogramsToPOT = false; //fScaleHistogramsToPOT = true;

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //----------------- POT normalization
    double potScalingBg = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, fPOTTotalNorm, potScalingBg, potScalingSignal, (fTreeDirName+"pottree").c_str());
    double potScaleBg = fScaleHistogramsToPOT? 1.*potScalingBg/fPOTTotalNorm: 1;
    double potScaleSignal = fScaleHistogramsToPOT? 1.*potScalingSignal/fPOTTotalNorm: 1;

    //--------- Matrix to store the number of events
    std::vector< std::vector<int> > nEventsMatrix;
    std::vector<PlotDef> cutDefsForTable;

    //--------- Loop over the cuts
    std::vector<AnaPlot> anaPlots;
    TCut previousCut("");
    for (size_t i = 0; i < cutDefs.size(); ++i) {

        TCut currentCut = TCut(cutDefs[i].GetCut());
        AnaPlot anaPlot(i, cutDefs[i], sampleDefs, phaseSpaceDefs);

        anaPlot.DrawHistograms(fTree, previousCut, 0, potScaleSignal, potScaleBg);
        anaPlots.push_back(anaPlot);

        if(cutDefs[i].GetAccumulateCut()){
            previousCut = previousCut && currentCut;
            anaPlot.DrawHistograms(fTree, previousCut, 1, potScaleSignal, potScaleBg);

            cutDefsForTable.push_back(cutDefs[i]);
            nEventsMatrix.push_back({});
            size_t cutIndex = nEventsMatrix.size()-1;
            std::map<std::string, int> nEventsMap = anaPlot.GetCountsV();
            for(const auto& sample : sampleDefs){
                nEventsMatrix[cutIndex].push_back(nEventsMap[sample.GetLabelS()]);
            }
        }
    
    }

    //loop over the samples and set NEvents
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