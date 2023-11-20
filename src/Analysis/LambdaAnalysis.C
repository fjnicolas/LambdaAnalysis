#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

std::string fTruthInFV = "TruthIsFiducial==1 &&";
std::string fTruthInAV = "abs(NuvX)<200 && abs(NuvY)<200 && NuvZ>0 && NuvZ<500 && ";

//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = {
    {fTruthInAV+"IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Signal", true}
    ,{fTruthInAV+"IntOrigin==1 && IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Background Nu", false}
    ,{"IntOrigin==1 && IntDirt==1", "Dirt", false}
    ,{"IntOrigin==2", "Cosmic", false}
    /*,{fTruthInFV+"IntNLambda==0 && IntMode==0 && abs(IntNuPDG)!=12", "QE", false}
    ,{fTruthInFV+"IntMode==1 && abs(IntNuPDG)!=12", "RES", false}
    ,{fTruthInFV+"IntMode==2 && abs(IntNuPDG)!=12", "DIS", false}
    ,{fTruthInFV+"(IntMode==3 || IntMode==10) && abs(IntNuPDG)!=12", "COH and MEC", false}
    ,{fTruthInFV+"abs(IntNuPDG)==12", "NuE", false} */
};

//-------- POT normalization
double fPOTTotalNorm = 3.3e20;

//---------  LATeX output file
std::string fOutputFileName = "CutEfficiencies";
std::string fOutputFileNameNormalized = "CutEfficienciesNormalized";

//--------- Cut definitions
double fCutCRUMBS = -0.2;
double fCutMinNAngles = 1;
double fCutFRANS = 0.15;
double fCutNOrigins = 4;
double fCutNOriginsM3 = 0;
double fCutFRANSPANDORA = 0.2;
double fCutNShw = 1;
double fCutShwEnergy = 135;

TCut fCounterCut = "SliceID==0";

std::vector<PlotDef> cutDefs = {
    {"TruthIsFiducial",  "0==0",            CutType::kNone,   0, {0,2,2}, true, "Truth in FV",  "No \\ cut"}
    ,{"CRUMBSScore",  "0==0",                CutType::kNone,   0, {-1,1,20}, false, "CRUMBSScore",  "CRUMBSScore"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    ,{"NOriginsMult1", "0==0", CutType::kNone, 0, {0, 6, 6}, false, "# origins mult 1",  "\\# \\ origins \\ mult \\ 1"}
    ,{"NOriginsMult2", "0==0", CutType::kNone, 0, {0, 6, 6}, false, "# origins mult 2",  "\\# \\ origins \\ mult \\ 2"}
    
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS,      {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}

    /*,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"AngleGap",                "AngleGap", CutType::kLeft, 20, {0, 30, 40}, true, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  6, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHitsTrack1>5 && AngleNHitsTrack2>5", "AngleNHitsTrack1>5 && AngleNHitsTrack2>5", CutType::kRight, 1, {0, 2, 2}, false, "min(track 1 # hits, track 2 # hits) ", "Min Track1/2 hits"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 20, {0, 200, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, true, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}*/

    ,{"NOriginsMultGT3", "NOriginsMultGT3", CutType::kLeft,  fCutNOriginsM3, {0, 5, 5}  , true, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NOrigins",        "NOrigins",        CutType::kLeft,  fCutNOrigins,   {0, 15, 15}, true, "# origins",        "\\#  \\ origins"}
    
    /*,{"AngleGap",                "AngleGap", CutType::kLeft,  20, {0, 30, 40}, true, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  7.5, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHits",              "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, true, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"AngleNHitsTrack1",  "AngleNHitsTrack1",  CutType::kRight, 5, {0, 200, 40}, true, "Angle track 1 # hits", "Angle \\ track \\ 1 \\# \\ hits"}
    ,{"AngleNHitsTrack2",  "AngleNHitsTrack2",  CutType::kRight, 5, {0, 200, 40}, true, "Angle track 2 # hits", "Angle \\ track \\ 2 \\# \\ hits"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 30, {0, 200, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}*/

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    
    ,{"AngleNHitsTrack1:AngleNHitsTrack2", "0==0", CutType::k2D, 10, {0, 50, 50}, false, "Test2D", "Test2D"}

    ,{"AngleGap",                "0==0", CutType::kLeft,  20, {0, 30, 40}, false, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "0==0", CutType::kLeft,  10, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHitsTrack1>5 && AngleNHitsTrack2>5", "AngleNHitsTrack1>5 && AngleNHitsTrack2>5", CutType::kRight, 1, {0, 2, 2}, false, "min(track 1 # hits, track 2 # hits) ", "Min Track1/2 hits"}
    ,{"NUnassociatedHits",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}


    ,{"NShwTh100", "NShwTh100",       CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>100MeV)", "NShower(>100MeV)"}
    ,{"NShwTh75", "NShwTh75",         CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>75MeV)", "NShower(>75MeV)"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"MainShowerEnergy", "0==0",     CutType::kLeft, 0, {0, 500, 100},  false, "MainShowerEnergy [MeV]", "MainShowerEnergy"}
    ,{"MainShowerScore", "0==0",      CutType::kLeft, 0, {0, 0.55, 20 }, false, "MainShowerScore", "MainShowerScore"}

    //,{"AngleNHits",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    //,{"AngleNHitsTrack1",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle track 1 # hits", "Angle \\ track \\ 1 \\# \\ hits"}
    //,{"AngleNHitsTrack2",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle track 2 # hits", "Angle \\ track \\ 2 \\# \\ hits"}
                                                    
}; 




//---------  Main function
void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

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
    

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, fPOTTotalNorm, potScaling, potScalingSignal);


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


