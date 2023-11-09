#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"


//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = {
     {"IntNLambda>0 && IntMode==0 && abs(IntNuPDG!=12)", "Signal", true}
    ,{"IntNLambda==0 && IntMode==0 && abs(IntNuPDG!=12)", "QE", false}
    ,{"IntNLambda==0 && IntMode==1 && abs(IntNuPDG!=12)", "RES", false}
    ,{"IntNLambda==0 && IntMode==2 && abs(IntNuPDG!=12)", "DIS", false}
    ,{"IntNLambda==0 && (IntMode==3 || IntMode==10) && abs(IntNuPDG!=12)", "COH and MEC", false}
    ,{"abs(IntNuPDG==12)", "NuE", false}
};


//--------- Cut definitions
double fCutMinNAngles = 1;
double fCutFRANS = 0.15;
double fCutNOrigins = 4;
double fCutNOriginsM3 = 0;
double fCutFRANSPANDORA = 0.2;
double fCutNShw = 1;
double fCutShwEnergy = 135;

std::vector<PlotDef> cutDefs = {
     {"TruthIsFiducial",  "0==0",            CutType::kNone,   0, {0,2,2}, true, "Truth in FV",  "No \\ cut"}
    ,{"TruthIsFiducial",  "TruthIsFiducial", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    ,{"NShwTh100", "NShwTh100",       CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>100MeV)", "NShower(>100MeV)"}
    ,{"NShwTh75", "NShwTh75",         CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>75MeV)", "NShower(>75MeV)"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"MainShowerEnergy", "0==0",     CutType::kLeft, 0, {0, 500, 100},  false, "MainShowerEnergy [MeV]", "MainShowerEnergy"}
    ,{"MainShowerScore", "0==0",      CutType::kLeft, 0, {0, 0.55, 20 }, false, "MainShowerScore", "MainShowerScore"}

    ,{"NOriginsMult1", "0==0", CutType::kNone, 0, {0, 6, 6}, false, "# origins mult 1",  "\\# \\ origins \\ mult \\ 1"}
    ,{"NOriginsMult2", "0==0", CutType::kNone, 0, {0, 6, 6}, false, "# origins mult 2",  "\\# \\ origins \\ mult \\ 2"}
    
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}

    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS,      {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}
    ,{"NOriginsMultGT3", "NOriginsMultGT3", CutType::kLeft,  fCutNOriginsM3, {0, 5, 5}  , true, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NOrigins",        "NOrigins",        CutType::kLeft,  fCutNOrigins,   {0, 15, 15}, true, "# origins",        "\\#  \\ origins"}
    
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}

    ,{"AngleGap",                "0==0", CutType::kLeft,  20, {0, 30, 40}, false, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "0==0", CutType::kLeft,  10, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHits",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}

    ,{"NShwTh100", "NShwTh100",       CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>100MeV)", "NShower(>100MeV)"}
    ,{"NShwTh75", "NShwTh75",         CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>75MeV)", "NShower(>75MeV)"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"MainShowerEnergy", "0==0",     CutType::kLeft, 0, {0, 500, 100},  false, "MainShowerEnergy [MeV]", "MainShowerEnergy"}
    ,{"MainShowerScore", "0==0",      CutType::kLeft, 0, {0, 0.55, 20 }, false, "MainShowerScore", "MainShowerScore"}
                                                    
}; 


//---------  Main function
void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    //---------  LATeX output file
    std::string fOutputFileName = "CutEfficiencies";

    //---------  Remove all *.pdf with gSystem
    gSystem->Exec("rm -rf OutputPlots");
    gSystem->Exec("mkdir OutputPlots");

    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );


    //--------- Loop over the cuts
    std::vector<AnaPlot> anaPlots;
    TCut previousCut("");
    for (size_t i = 0; i < cutDefs.size(); ++i) {

        TCut currentCut = TCut(cutDefs[i].GetCut());
        AnaPlot anaPlot(i, cutDefs[i], sampleDefs);

        anaPlot.DrawHistograms(fTree, previousCut);
        anaPlots.push_back(anaPlot);

        if(cutDefs[i].GetAccumulateCut()){
            previousCut = previousCut && currentCut;
            anaPlot.DrawHistograms(fTree, previousCut, 1);
        }
    }


    //--------- Create the LaTeX table
    GenerateAndCompileTeXTable(sampleDefs, anaPlots, 1, fOutputFileName, "Cut efficiencies");
   
    //--------- Output hand scans
    CreateHandScanList(fTree, fTreeHeader, previousCut, sampleDefs);

    return;
}


