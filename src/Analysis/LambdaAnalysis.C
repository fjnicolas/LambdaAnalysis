#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesLATeXInterface.C"


std::string fCutFRANSPANDORA = "0.15";
std::string fCutFRANS = "0.15 ";
std::string fCutNOrigins = "4";
std::string fCutNOriginsM3 = "0";
std::string fCutNShw = "1";

std::vector<PlotDef> cutDefs1 = {
                                    {"", "TruthIsFiducial", "TruthIsFiducial>=0", "Truth in FV", "No \\ cut", {0,2,2}, false, false},
                                    {"", "TruthIsFiducial", "TruthIsFiducial==1", "Truth in FV", "Truth \\ in \\ FV", {0,2,2}, true, false},
                                    {"", "RecoIsFiducial", "RecoIsFiducial==1", "Reco in FV", "Reco \\ in \\ FV", {0,2,2}, true, false},

                                    //{"", "FRANSScorePANDORA", "FRANSScorePANDORA>"+fCutFRANSPANDORA, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA>"+fCutFRANSPANDORA, {-.5,.5,40}, true},

                                    //{"", "NShwTh75", "!(NOrigins>4 && NShwTh75<="+fCutNShw+")", "# showers (>100MeV)", "NShower(>100MeV)<="+fCutNShw, {0,6, 6}, true},

                                    {"", "NOriginsMult1", "0==0", "# origins mult 1",  "\\# \\ origins \\ mult \\ 1", {0, 6, 6}, false, false},
                                    {"", "NOriginsMult2", "0==0", "# origins mult 2",  "\\# \\ origins \\ mult \\ 2", {0, 6, 6}, false, false},
                                    {"", "NOriginsMultGT3", "0==0", "# origins mult >=3",  "\\# \\ origins \\ mult \\ >= 3", {0, 6, 6}, false, false},

                                    {"", "NShwTh100", "NShwTh100<="+fCutNShw, "# showers (>100MeV)", "NShower(>100MeV)<="+fCutNShw, {0, 6, 6}, false},
                                    {"", "ShowerEnergy", "ShowerEnergy<=135", "Shower Energy [MeV]", "ShowerEnergy<=135 \\ [MeV]", {0, 300, 15}, false},

                                    {"", "NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", "# origins mult 1>0, # origins mult2>0", "\\# \\ origins \\ mult \\ 1 >0, \\# \\ origins \\ mult \\ 2 > 0", {0,2,2}, true},
                                    {"", "NAngles", "NAngles>0", "# V", "\\# \\ V>0", {0,5,5}, true},
                                    
                                    {"", "AngleFRANSScore", "AngleFRANSScore>"+fCutFRANS, "V FRANS score", "V \\ FRANS \\ score>"+fCutFRANS, {-.5,.5,40}, true},
                                    {"", "NOriginsMultGT3", "NOriginsMultGT3<="+fCutNOriginsM3, "# origins mult 3",  "\\# \\ origins \\ mult \\ 3<="+fCutNOriginsM3, {0, 5, 5}, true},
                                    {"", "NOrigins", "NOrigins<="+fCutNOrigins, "# origins", "\\#  \\ origins  <="+fCutNOrigins, {0,15,15}, true},
                                    
                                    {"", "NShwTh100", "NShwTh100<="+fCutNShw, "# showers (>100MeV)", "NShower(>100MeV)<="+fCutNShw, {0, 6, 6}, true},
                                    {"", "ShowerEnergy", "ShowerEnergy<=135", "Shower Energy [MeV]", "ShowerEnergy<=135 \\ [MeV]", {0, 300, 15}, true},

                                    {"", "FRANSScorePANDORA", "FRANSScorePANDORA>"+fCutFRANSPANDORA, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA>"+fCutFRANSPANDORA, {-.5,.5,40}, true},
                                
                                    //{"", "NShwTh75", "NShwTh75<="+fCutNShw, "# showers (>75MeV)", "NShower(>75MeV)<="+fCutNShw, {0,6, 6}, true},
                                    //{"", "ShowerEnergy", "ShowerEnergy<=135", "Shower energy [MeV]", "Shower \\ energy < 135 \\ [MeV]", {0,500,25}, true},
                                   
                                    //{"", "NShwTh75", "NShwTh75<="+fCutNShw, "# showers (>75MeV)", "NShower(>75MeV)<="+fCutNShw, {0,6, 6}, true},
                                                                
                                    
                                    
                                };       


std::vector<std::string> sampleDefs1 = { "IntNLambda>0 && IntMode==0", "IntNLambda==0 && IntMode==0", "IntNLambda==0 && IntMode==1" };
std::vector<std::string> sampleDefNames1 = { "Signal", "QE", "RES" };

std::vector<PlotDef> cutDefs = cutDefs1; 

std::vector<std::string> sampleDefs = sampleDefs1;
std::vector<std::string> sampleDefNames = sampleDefNames1;

std::vector<int> fColors = {kRed+1, kBlue-4, kGreen+4, kOrange+5};
std::vector<int> fLineStyle = {1, 2, 9, 3, 6, 10};

void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    // Set batch mode
    gROOT->SetBatch(batchMode);
    SetGStyle();

    // shell command remove all *.pdf with gSystem
    gSystem->Exec("rm *.pdf");

    //--------- Configuration Parameters
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    fFile->ls();
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    std::vector<std::vector<int>> histogramCounts(cutDefs.size(), std::vector<int>(sampleDefs.size(), 0));
    
    std::vector<std::vector<TH1F*>> histograms;
    for (size_t i = 0; i < cutDefs.size(); ++i) {
        std::vector<TH1F*> hVec;
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            std::string plotName = std::to_string(i)+"_"+std::to_string(j);
            TH1F *hAux = new TH1F(plotName.c_str(), plotName.c_str(), cutDefs[i].bins.nBins, cutDefs[i].bins.x1, cutDefs[i].bins.x2);
            hVec.push_back(hAux);
        }
        histograms.push_back(hVec);
    }

    // Cut defitions
    TCut previousCut("");

    for (size_t i = 0; i < cutDefs.size(); ++i) {

        TCut currentCut = TCut(cutDefs[i].cut.c_str()); 
        
        TCanvas *c1 = new TCanvas("c1","c1",800,600);
        TCanvas *c2 = new TCanvas("c2"," c2",800,600);
        std::string plotAxisLabels=";"+cutDefs[i].varLabel+"; # events";
        TH1F *hAux = new TH1F("hh", plotAxisLabels.c_str(), cutDefs[i].bins.nBins, cutDefs[i].bins.x1, cutDefs[i].bins.x2);
        hAux->Draw();

        // Add legend
        TLegend *legend = new TLegend(0.75, 0.75, 0.9, 0.9);

        if(cutDefs[i].log==true){
            hAux->SetMinimum(0.5);
            gStyle->SetOptLogy(1);
        }
        else{
            hAux->SetMinimum(0);
            gStyle->SetOptLogy(0);
        }

        int maxVal=0;
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            
            TCut sampelCut(sampleDefs[j].c_str());
            sampelCut = sampelCut;
            
            std::string plotName = std::to_string(i)+"_"+std::to_string(j);
            
            c2->cd();
            fTree->Draw( (cutDefs[i].var+">>"+plotName).c_str(), sampelCut && previousCut );
            std::cout<<"Plotting: "<<cutDefs[i].var<<std::endl;

            TH1F *hCount = new TH1F("hCount", "hCount", 2, 0, 2);
            fTree->Draw( (cutDefs[i].var+">>hCount").c_str(), sampelCut && previousCut && currentCut );
            
            c1->cd();
            histograms[i][j]->SetLineColor(fColors[j]);
            histograms[i][j]->SetLineWidth(3);
            histograms[i][j]->SetLineStyle(fLineStyle[j]);
            histograms[i][j]->SetStats(0);
            legend->AddEntry(histograms[i][j], sampleDefNames[j].c_str(), "l");
            
            histogramCounts[i][j] = hCount->GetEntries();
            if(histograms[i][j]->GetMaximum()>maxVal) maxVal = histograms[i][j]->GetMaximum();
        }



        hAux->GetYaxis()->SetRangeUser(0.5, maxVal*1.1);
        hAux->SetStats(0);
        // hAux in log scale if the cut is in log

        
        hAux->Draw();
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            histograms[i][j]->Draw("hist same");
        }

        // add legent transparency
        legend->SetFillColorAlpha(0, 0);
        
        legend->Draw("same");
        std::cout<<"Plotting: "<<cutDefs[i].log<<std::endl;
        

        c1->Update();
        c1->SaveAs(("plot"+std::to_string(i)+cutDefs[i].var+".pdf").c_str());
        delete c1;
        
        if(cutDefs[i].accumulateCut)
            previousCut = previousCut && currentCut;
    }

    for (size_t i = 0; i < cutDefs.size(); ++i) {
        std::cout << cutDefs[i].varLabel << "\t";
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            std::cout << histogramCounts[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    

    generateAndCompileTeXTable(cutDefs, sampleDefNames, histogramCounts, "output.tex", "Cut efficiencies");



    // Create a 2D histogram using TTree::Draw
    std::cout<<" BAD BACKGROUNDS QE:\n";
    std::cout<<previousCut<<std::endl;
    TH2D* h2_all = new TH2D("h2_all", "2D Histogram", 1000, 1, 1001, 1000, 1, 1001);
    TCut debugCut = previousCut && ( TCut(sampleDefs[1].c_str())  );
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

    debugCut = previousCut && ( TCut(sampleDefs[2].c_str()) );
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

