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


//---------  Main function
void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    //---------  Remove all *.pdf with gSystem
    gSystem->Exec("rm -rf OutputPlots");
    gSystem->Exec("mkdir OutputPlots");
    gSystem->Exec("mkdir OutputPlots/PhaseSpace");

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
    ReadPOT(fFile, fPOTTotalNorm, potScaling, potScalingSignal);


    //--------- Loop over the cuts
    std::vector<AnaPlot> anaPlots;
    TCut previousCut("");
    for (size_t i = 0; i < cutDefs.size(); ++i) {

        TCut currentCut = TCut(cutDefs[i].GetCut());
        AnaPlot anaPlot(i, cutDefs[i], sampleDefs, phaseSpaceDefs);

        anaPlot.DrawHistograms(fTree, previousCut);
        anaPlots.push_back(anaPlot);

        if(cutDefs[i].GetAccumulateCut()){
            previousCut = previousCut && currentCut;
            anaPlot.DrawHistograms(fTree, previousCut, 1);
        }
    }


    //loop over the samples and set NEvents
    for(auto& sample : sampleDefs){
        int nEvents = fTree->Draw( "", TCut(sample.GetVar())+fCounterCut, "goff");
        sample.SetNEvents(nEvents);
        std::cout<<"Sample: "<<sample.GetLabelS()<<" NEvents: "<<nEvents<<std::endl;
    }
    
    //--------- Create the LaTeX table
    GenerateAndCompileTeXTable(sampleDefs, anaPlots, 0, fOutputFileName, "Cut efficiencies");

    //--------- Create the LaTeX table (POT normalized)
    GenerateAndCompileTeXTable(sampleDefs, anaPlots, 0, fOutputFileNameNormalized, "Cut efficiencies", potScaling, potScalingSignal, fPOTTotalNorm);
   
    //--------- Output hand scans
    CreateHandScanList(fTree, fTreeHeader, previousCut, sampleDefs);
    

    return;
}


