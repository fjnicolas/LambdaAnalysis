#ifndef CUTEFFICIENCIES_LATEXINTERFACE_H
#define CUTEFFICIENCIES_LATEXINTERFACE_H

void GenerateAndCompileTeXTable(
    const std::vector<SampleDef>& sampleDefs,
    const std::vector<AnaPlot>& anaPlots,
    size_t denominatorIndex,
    const std::string& fileName,
    const std::string& tableCaption
) {
    // Open a file stream for writing the .tex file
    std::ofstream texFile(fileName+".tex");

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
    for(size_t k = 0; k<sampleDefs.size(); k++) texFile<<"c|";
    texFile << "}" << std::endl;
    texFile << "\\hline" << std::endl;


    std::map<std::string, int> histogramCounts0 = anaPlots[denominatorIndex].GetCountsV();

    // Write the header row
    texFile << " Cut & ";
    int cont=0;
    
    for(auto& sampleName : histogramCounts0){
        if( cont == sampleDefs.size()-1 )
            texFile << sampleName.first;
        else
            texFile << sampleName.first << " & ";
        
        cont++;
    }

    texFile << "\\\\ \\hline" << std::endl;

    // Write the data rows
    for (size_t i = 0; i < anaPlots.size(); ++i) {
        if(anaPlots[i].GetPlotDef().GetAccumulateCut() == false) continue;
        texFile << "$ {\\rm " << anaPlots[i].GetPlotDef().GetCutLabel() << "}$" << " & ";
        std::map<std::string, int> histogramCounts = anaPlots[i].GetCountsV();
        
        cont=0;
        for(auto& sampleName : histogramCounts){
            std::ostringstream streamObjEff;
            streamObjEff <<  std::fixed << std::setprecision(2) << 100.*sampleName.second/histogramCounts0[sampleName.first];
            // Get string from out
            if(  cont == sampleDefs.size()-1 )
                texFile << sampleName.second <<" ("<< streamObjEff.str() <<" \\%)";
            else
                texFile << sampleName.second <<" ("<< streamObjEff.str() <<" \\%)"  << " & ";
            cont++;
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

    gSystem->Exec( ("mv "+fileName+".pdf OutputPlots").c_str() );
    gSystem->Exec( ("rm "+fileName+".{aux,log,tex,pdf}").c_str() );

}



void CreateHandScanList(TTree *fTree, TTree *fTreeHeader, TCut cut, std::vector<SampleDef> sampleDefs){

    //--------- Loop over the TreeHeader
    unsigned int fSubRunId;
    fTreeHeader->SetBranchAddress("SubRunID", &fSubRunId);
    std::string *fLArInputFileName = new std::string;
    fTreeHeader->SetBranchAddress("InputFileName", &fLArInputFileName);
    // Create map of subrun filename
    std::map<int, std::string> subrunFilenameMap;
    for(size_t i=0; i<fTreeHeader->GetEntries(); ++i){
        fTreeHeader->GetEntry(i);
        
        // check the string includes the substring "Inclusive"
        if(fLArInputFileName->find("Inclusive") == std::string::npos){
            continue;
        }

        subrunFilenameMap[fSubRunId] = *fLArInputFileName;
    }

    // Output 
    std::ofstream handScanEvents;
    handScanEvents.open("OutputPlots/handScanEvents.txt");

    // loop over signal types
    for(size_t k=0; k<sampleDefs.size(); k++){
        // Skip if it is the signal
        if(sampleDefs[k].IsSignal()) continue;
        
        std::cout<<" Sample: "<<sampleDefs[k].GetLabel()<<std::endl;
        
        // Create the histogram
        TH2D* h2_all = new TH2D("h2_all", "2D Histogram", 1000, 1, 1001, 1000, 1, 1001);
        
        // Create the cut
        TCut debugCut = cut && ( TCut(sampleDefs[k].GetVar())  );
        
        // Draw the histogram
        // Y axis is event ID
        // X axis is subrun ID
        fTree->Draw("EventID:SubrunID>>h2_all", debugCut);
        // Get the number of bins in the histogram
        int numBinsX = h2_all->GetNbinsX();
        int numBinsY = h2_all->GetNbinsY();
        for (int sr = 1; sr <= numBinsX; sr++) {
            for (int e = 1; e <= numBinsY; e++) {
                double binValue = h2_all->GetBinContent(sr, e);
                if(binValue>0){
                    std::cout<<"SubRunID: "<<sr<<" EventID: "<<e<<" BinValue: "<<binValue<<std::endl;
                    std::cout<<"Filename: "<<subrunFilenameMap[sr]<<std::endl;
                    handScanEvents<<"1:"<<sr<<":"<<e<<" "<<subrunFilenameMap[sr]<<std::endl;
                }
            }
        }        
    }

    handScanEvents.close();

    return;

}


#endif