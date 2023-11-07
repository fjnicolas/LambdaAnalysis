#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

std::vector<SampleDef> sampleDefs = {
    {"IntNLambda>0 && IntMode==0", "Signal", true},
    {"IntNLambda==0 && IntMode==0", "QE", false},
    {"IntNLambda==0 && IntMode==1", "RES", false},
    {"IntNLambda==0 && IntMode==2", "DIS", false},
    {"IntNLambda==0 && (IntMode==3 || IntMode==10)", "COH and MEC", false}
};

std::string fCutFRANSPANDORA = "0.15";
std::string fCutFRANS = "0.15 ";
std::string fCutNOrigins = "4";
std::string fCutNOriginsM3 = "0";
std::string fCutNShw = "1";

std::vector<PlotDef> cutDefs = {
    {"", "TruthIsFiducial", "TruthIsFiducial>=0", 1, "Truth in FV",  "No \\ cut", CutType::kCenter, {0,2,2},  true},
    {"", "TruthIsFiducial", "TruthIsFiducial==1", 1, "Truth in FV", "Truth \\ in \\ FV", CutType::kCenter, {0,2,2}, true},
    {"", "RecoIsFiducial", "RecoIsFiducial==1", 1, "Reco in FV", "Reco \\ in \\ FV", CutType::kCenter, {0,2,2}, true},

    //{"", "FRANSScorePANDORA", "FRANSScorePANDORA>"+fCutFRANSPANDORA, 0.15, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA>"+fCutFRANSPANDORA, CutType::kRight,{-.5,.5,40}, true},

    {"", "NShwTh100", "NShwTh100<="+fCutNShw, 1, "# showers (>100MeV)", "NShower(>100MeV)<="+fCutNShw, CutType::kLeft, {0, 6, 6}, false, true},
    {"", "ShowerEnergy", "ShowerEnergy<=135", 135, "Shower Energy [MeV]", "ShowerEnergy<=135 \\ [MeV]", CutType::kLeft, {0, 300, 15}, false, true},
    {"", "MainShowerEnergy", "0==0", 0, "MainShowerEnergy [MeV]", "", CutType::kLeft, {0, 500, 100}, false, true},
    {"", "MainShowerScore", "0==0", 0, "MainShowerScore", "", CutType::kLeft, {0, 0.55, 20 }, false, true},


    {"", "NOriginsMult1", "0==0",  1, "# origins mult 1",  "\\# \\ origins \\ mult \\ 1", CutType::kRight, {0, 6, 6}, false},
    {"", "NOriginsMult2", "0==0", 1, "# origins mult 2",  "\\# \\ origins \\ mult \\ 2", CutType::kRight, {0, 6, 6}, false},
    
    {"", "NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", 1, "# origins mult 1>0, # origins mult2>0", "\\# \\ origins \\ mult \\ 1 >0, \\# \\ origins \\ mult \\ 2 > 0", CutType::kCenter, {0,2,2}, true},
    {"", "NAngles", "NAngles>0", 1, "# V", "\\# \\ V>0", CutType::kRight,{0,5,5}, true},

    {"", "AngleGap", "0==0", 0, "Gap [cm]", "Gap \\ [cm]", CutType::kLeft,{0, 50, 100}, false, true},
    {"", "AngleDecayContainedDiff", "0==0", 0, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle", CutType::kLeft, {0, 20, 40}, false, true},
    {"", "AngleNHits", "0==0", 0, "Angle # hits", "Angle \\ \\# \\ hits", CutType::kRight,{0, 250, 25}, false, true},

    
    {"", "AngleFRANSScore", "AngleFRANSScore>"+fCutFRANS, 0.15, "V FRANS score", "V \\ FRANS \\ score>"+fCutFRANS, CutType::kRight, {-.5,.5,40}, true, true},
    {"", "NOriginsMultGT3", "NOriginsMultGT3<="+fCutNOriginsM3, 1, "# origins mult 3",  "\\# \\ origins \\ mult \\ 3<="+fCutNOriginsM3, CutType::kLeft, {0, 5, 5}, true},
    {"", "NOrigins", "NOrigins<="+fCutNOrigins, 4, "# origins", "\\#  \\ origins  <="+fCutNOrigins, CutType::kLeft, {0,15,15}, true},
    
    {"", "FRANSScorePANDORA", "FRANSScorePANDORA>"+fCutFRANSPANDORA, 0.15, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA>"+fCutFRANSPANDORA, CutType::kRight,{-.5,.5,40}, true, true},

    {"", "AngleGap", "0==0", 0, "Gap [cm]", "Gap \\ [cm]", CutType::kLeft,{0, 50, 100}, false, true},
    {"", "AngleDecayContainedDiff", "0==0", 0, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle", CutType::kLeft, {0, 20, 40}, false, true},
    {"", "AngleNHits", "0==0", 0, "Angle # hits", "Angle \\ \\# \\ hits", CutType::kRight,{0, 250, 50}, false, true},

    //{"", "NShwTh75", "NShwTh75<="+fCutNShw, "# showers (>75MeV)", "NShower(>75MeV)<="+fCutNShw, CutType::kRight,{0,6, 6}, true},
    //{"", "ShowerEnergy", "ShowerEnergy<=135", "Shower energy [MeV]", "Shower \\ energy < 135 \\ [MeV]", CutType::kRight,{0,500,25}, true},
    
    //{"", "NShwTh75", "NShwTh75<="+fCutNShw, "# showers (>75MeV)", "NShower(>75MeV)<="+fCutNShw, CutType::kRight,{0,6, 6}, true},
                                                    
}; 





void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    //---------  Remove all *.pdf with gSystem
    gSystem->Exec("rm *.pdf");
    gSystem->Exec("rm OutputPlots.root");

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    fFile->ls();
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    std::vector<AnaPlot> anaPlots;

    // Cut defitions
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



    generateAndCompileTeXTable(sampleDefs, anaPlots, 1, "output.tex", "Cut efficiencies");



    // Create a 2D histogram using TTree::Draw
    std::cout<<" BAD BACKGROUNDS QE:\n";
    std::cout<<previousCut<<std::endl;
    TH2D* h2_all = new TH2D("h2_all", "2D Histogram", 1000, 1, 1001, 1000, 1, 1001);
    TCut debugCut = previousCut && ( TCut(sampleDefs[1].GetVar())  );
    fTree->Draw("EventID:SubrunID>>h2_all", debugCut);

    // Get the number of bins in the histogram
    int numBinsX = h2_all->GetNbinsX();
    int numBinsY = h2_all->GetNbinsY();
    for (int i = 1; i <= numBinsX; i++) {
        for (int j = 1; j <= numBinsY; j++) {
            double binValue = h2_all->GetBinContent(i, j);
            if(binValue>0){
                std::cout<<"Bin: "<<i<<", "<<j<<", value: "<<binValue<<std::endl;
            }
        }
    }

    debugCut = previousCut && ( TCut(sampleDefs[2].GetVar()) );
    fTree->Draw("EventID:SubrunID>>h2_all", debugCut);
    std::cout<<" BAD BACKGROUNDS RES:\n";
    // Get the number of bins in the histogram
    for (int i = 1; i <= numBinsX; i++) {
        for (int j = 1; j <= numBinsY; j++) {
            double binValue = h2_all->GetBinContent(i, j);
            if(binValue>0){
                std::cout<<"Bin: "<<i<<", "<<j<<", value: "<<binValue<<std::endl;
            }
        }
    }

    return;
}

