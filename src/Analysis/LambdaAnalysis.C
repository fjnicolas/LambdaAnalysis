#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

std::vector<std::string> sampleDefs1 = { "IntNLambda>0 && IntMode==0", "IntNLambda==0 && IntMode==0", "IntNLambda==0 && IntMode==1" };
std::vector<std::string> sampleDefNames1 = { "Signal", "QE", "RES" };

std::vector<std::string> cutDefs1 = { "RunID>0", "NOriginsPairOneTwo>0", "NAngles>0", "AngleFRANSScore>0.1", "NOriginsMultGT3==0", "NOrigins<=6" };

std::vector<std::string> cutDefs2 = { "RunID>0", "FRANSScorePANDORA>0.1", "NOriginsPairOneTwo>0", "NAngles>0", "AngleFRANSScore>0.1", "NOriginsMultGT3==0", "NOrigins<=6" };

std::vector<std::string> cutDefs = cutDefs2; 

std::vector<std::string> sampleDefs = sampleDefs1;
std::vector<std::string> sampleDefNames = sampleDefNames1;


void generateAndCompileTeXTable(
    const std::vector<std::string>& cutDefs,
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
    texFile << "\\begin{document}" << std::endl;

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
        texFile << cutDefs[i] << " & ";
        for (size_t j = 0; j < histogramCounts[i].size(); ++j) {
            std::ostringstream streamObjEff;
            streamObjEff <<  std::fixed << std::setprecision(2) << 100.*histogramCounts[i][j]/histogramCounts[0][j];
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


void LambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "LambdaAnaTree/", std::string fTreeName = "LambdaAnaTree")
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
        currentCut = currentCut && TCut(cutDefs[i].c_str()); 
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            TCut sampelCut(sampleDefs[j].c_str());
            sampelCut = sampelCut && currentCut;
            TH1F *h1 = new TH1F("h1",cutDefs[i].c_str(),2,0,2);
            fTree->Draw("RunID>>h1",sampelCut);
            h1->Draw();

            histogramCounts[i][j] = h1->GetEntries();

            c1->Update();
            c1->WaitPrimitive();

        }
        
        
    }

        for (size_t i = 0; i < cutDefs.size(); ++i) {
        std::cout << cutDefs[i] << "\t";
        for (size_t j = 0; j < sampleDefs.size(); ++j) {
            std::cout << histogramCounts[i][j] << "\t";
        }
        std::cout << std::endl;
    }


    generateAndCompileTeXTable(cutDefs, sampleDefNames, histogramCounts, "output.tex", "Cut efficiencies");



    return;
}

