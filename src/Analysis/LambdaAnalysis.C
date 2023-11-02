#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesLATeXInterface.C"


std::vector<PlotDef> cutDefs1 = {
                                    {"", "TruthIsFiducial", "TruthIsFiducial>=0", "No \\ cut", "", "", true},
                                    {"", "TruthIsFiducial", "TruthIsFiducial==1", "Truth \\ in \\ FV", "", "", true},
                                    {"", "RecoIsFiducial", "RecoIsFiducial==1", "Reco \\ in \\ FV", "", "", true},
                                    {"", "NOriginsPairOneTwo", "NOriginsPairOneTwo>0", "\\# \\ 2-1 \\ origin \\ pairs > 0", "", "", true},
                                    {"", "NAngles", "NAngles>0", "\\# \\ V>0", "", "", true},
                                    {"", "AngleFRANSScore", "AngleFRANSScore>0.1", "V \\ FRANS \\ score>0.1", "", "", true},
                                    {"", "NOriginsMultGT3", "NOriginsMultGT3==0", "\\# \\ origins \\ mult \\ 3==0", "", "", true},
                                    {"", "NOrigins", "NOrigins<=5", "\\# \\ origins <=5", "", "", true},
                                    {"", "FRANSScorePANDORA", "FRANSScorePANDORA>0.15", "FRANS \\ score \\ PANDORA>0.15", "", "", true}
                                };       

std::vector<PlotDef> cutDefs2 = {
                                    {"", "TruthIsFiducial", "TruthIsFiducial>=0", "No \\ cut", "", "", true},
                                    {"", "TruthIsFiducial", "TruthIsFiducial==1", "Truth \\ in \\ FV", "", "", true},
                                    {"", "RecoIsFiducial", "RecoIsFiducial==1", "Reco \\ in \\ FV", "", "", true},
                                     {"", "FRANSScorePANDORA", "FRANSScorePANDORA>0.15", "FRANS \\ score \\ PANDORA>0.15", "", "", true},
                                    {"", "NOriginsPairOneTwo", "NOriginsPairOneTwo>0", "\\# \\ 2-1 \\ origin \\ pairs > 0", "", "", true},
                                    {"", "NAngles", "NAngles>0", "\\# \\ V>0", "", "", true},
                                    {"", "AngleFRANSScore", "AngleFRANSScore>0.1", "V \\ FRANS \\ score>0.1", "", "", true},
                                    {"", "NOriginsMultGT3", "NOriginsMultGT3==0", "\\# \\ origins \\ mult \\ 3==0", "", "", true},
                                    {"", "NOrigins", "NOrigins<=5", "\\# \\ origins <=5", "", "", true}
                                };                                     
                                   


std::vector<std::string> sampleDefs1 = { "IntNLambda>0 && IntMode==0", "IntNLambda==0 && IntMode==0", "IntNLambda==0 && IntMode==1" };
std::vector<std::string> sampleDefNames1 = { "Signal", "QE", "RES" };


std::vector<PlotDef> cutDefs = cutDefs2; 

std::vector<std::string> sampleDefs = sampleDefs1;
std::vector<std::string> sampleDefNames = sampleDefNames1;


void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    // Set batch mode
    gROOT->SetBatch(batchMode);

    //--------- Configuration Parameters
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    fFile->ls();
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    std::vector<std::vector<int>> histogramCounts(cutDefs.size(), std::vector<int>(sampleDefs.size(), 0));

    // Cut defitions
    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    TCut currentCut("");

    for (size_t i = 0; i < cutDefs.size(); ++i) {
        currentCut = currentCut && TCut(cutDefs[i].cut.c_str()); 
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            TCut sampelCut(sampleDefs[j].c_str());
            sampelCut = sampelCut && currentCut;
            TH1F *h1 = new TH1F("h1",cutDefs[i].cut.c_str(),2,0,2);
            fTree->Draw("RunID>>h1",sampelCut);
            h1->Draw();

            histogramCounts[i][j] = h1->GetEntries();

            c1->Update();
            c1->WaitPrimitive();

        }
        
        
    }

        for (size_t i = 0; i < cutDefs.size(); ++i) {
        std::cout << cutDefs[i].varLabel << "\t";
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            std::cout << histogramCounts[i][j] << "\t";
        }
        std::cout << std::endl;
    }


    generateAndCompileTeXTable(cutDefs, sampleDefNames, histogramCounts, "output.tex", "Cut efficiencies");



    return;
}

