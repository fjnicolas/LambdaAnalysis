#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

struct PlotDef {

    std::string suffix;
    std::string var;
    std::string cut;
    std::string varLabel;
    std::string cutLabel;
    std::string bins;
    bool accumulateCut;

    PlotDef(const std::string& _suffix = "",
            const std::string& _var = "",
            const std::string& _cut = "",
            const std::string& _varLabel = "",
            const std::string& _cutLabel = "",
            const std::string& _bins = "",
            bool _accumulateCut = true)
        : suffix(_suffix),
          var(_var),
          cut(_cut),
          varLabel(_varLabel),
          cutLabel(_cutLabel),
          bins(_bins),
          accumulateCut(_accumulateCut) {
    }
};


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


void generateAndCompileTeXTable(
    const std::vector<PlotDef>& cutDefs,
    const std::vector<std::string>& sampleDefNames,
    const std::vector<std::vector<int>>& histogramCounts,
    const std::string& fileName,
    const std::string& tableCaption
) {
    // Open a file stream for writing the .tex file
    std::ofstream texFile(fileName);

    // Check if the file is open
    if (!texFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return;
    }

    // Write the LaTeX document preamble
    
    texFile << "\\documentclass{article}" << std::endl;
    texFile << "\\usepackage{graphicx}" << std::endl;
    texFile << "\\usepackage{pdflscape}"<<std::endl;
    texFile << "\\begin{document}" << std::endl;
    texFile << "\\begin{landscape}" << std::endl;


    // Write the LaTeX code for the table
    texFile << "\\begin{table}[h]" << std::endl;
    texFile << "\\centering" << std::endl;
    texFile << "\\begin{tabular}{|c|";
    for(size_t k = 0; k<sampleDefNames.size(); k++) texFile<<"c|";
    texFile << "}" << std::endl;
    texFile << "\\hline" << std::endl;

    // Write the header row
    texFile << " Cut & ";
    for(size_t k = 0; k<sampleDefNames.size(); k++){
        if(k==sampleDefNames.size()-1) texFile << sampleDefNames[k];
        else texFile << sampleDefNames[k] << " & ";
    }
    texFile << "\\\\ \\hline" << std::endl;

    // Write the data rows
    for (size_t i = 0; i < histogramCounts.size(); ++i) {
        texFile << "$ {\\rm " << cutDefs[i].varLabel << "}$" << " & ";
        for (size_t j = 0; j < histogramCounts[i].size(); ++j) {
            std::ostringstream streamObjEff;
            streamObjEff <<  std::fixed << std::setprecision(2) << 100.*histogramCounts[i][j]/histogramCounts[1][j];
            std::cout<<"Efficiency: "<<streamObjEff.str()<<" %"<<std::endl;
            // Get string from out
            if(j==sampleDefNames.size()-1) texFile << histogramCounts[i][j] <<" ("<< streamObjEff.str() <<" \\%)" << std::endl;
            else texFile << histogramCounts[i][j] <<" ("<< streamObjEff.str() <<" \\%)"  << " & ";
        }
        texFile << "\\\\ \\hline" << std::endl;
    }
 
    texFile << "\\end{tabular}" << std::endl;
    texFile << "\\caption{" << tableCaption << "}" << std::endl;
    texFile << "\\end{table}" << std::endl;

    texFile << "\\end{landscape}" << std::endl;

    // Close the LaTeX document
    texFile << "\\end{document}" << std::endl;

    // Close the file
    texFile.close();

    // Compile the .tex file
    std::string compileCommand = "pdflatex " + fileName;
    int compileResult = std::system(compileCommand.c_str());

    if (compileResult != 0) {
        std::cerr << "Failed to compile the .tex file." << std::endl;
    }
}


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

