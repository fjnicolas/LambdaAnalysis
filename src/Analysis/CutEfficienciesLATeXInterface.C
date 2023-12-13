#ifndef CUTEFFICIENCIES_LATEXINTERFACE_H
#define CUTEFFICIENCIES_LATEXINTERFACE_H

#define WriteRelativeEfficiencies 0

//--------- Create the LaTeX table (POT normalized)
void GenerateAndCompileTeXTable(
    const std::vector<SampleDef>& sampleDefs,
    const std::vector<AnaPlot>& anaPlots,
    size_t denominatorIndex,
    const std::string& fileName,
    const std::string& tableCaption,
    double potNorm=1,
    double potNormSignal=1,
    double totalPOTNorm=1
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
    for(size_t k = 0; k<sampleDefs.size(); k++)
        texFile<<"c|";
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

    std::map<std::string, double> potScalingMap;
    for(auto& sampleName : histogramCounts0){
        if(sampleName.first.find("Signal") != std::string::npos)
            potScalingMap[sampleName.first] = potNormSignal;
        else
            potScalingMap[sampleName.first] = potNorm;
    }

    // Write the data rows
    // last stored index
    size_t lastStoredIndex = 0;
    for (size_t i = 0; i < anaPlots.size(); ++i) {
        if(anaPlots[i].GetPlotDef().GetAccumulateCut() == false) continue;
        texFile << "$ {\\rm " << anaPlots[i].GetPlotDef().GetCutLabel() << "}$" << " & ";
        std::map<std::string, int> histogramCounts = anaPlots[i].GetCountsV();
        
        std::map<std::string, int> histogramCountsPre = anaPlots[lastStoredIndex].GetCountsV();

        cont=0;
        for(auto& sampleName : histogramCounts){

            double potScaling = potScalingMap[sampleName.first];
            
            std::ostringstream streamObjEff;
            streamObjEff << " (" << std::fixed << std::setprecision(2) << 100.*sampleName.second/histogramCounts0[sampleName.first];
            if(WriteRelativeEfficiencies==1){
                streamObjEff << " \\%) \\ $\\epsilon_{r}$= " << std::setprecision(0) << 100.*sampleName.second/histogramCountsPre[sampleName.first] << "\\%";
            }
            else{
                streamObjEff << " \\%)";
            }
            // Get string from out
            if(  cont == sampleDefs.size()-1 )
                texFile << potScaling*sampleName.second << streamObjEff.str();
            else
                texFile << potScaling*sampleName.second << streamObjEff.str() << " & ";
            cont++;
        }

        texFile << "\\\\ \\hline" << std::endl;

        lastStoredIndex = i;
    }
    



    texFile << "$ {\\rm " << "Final eff" << "}$" << " & ";
    
    cont=0;
    std::map<std::string, int> histogramCounts = anaPlots[ lastStoredIndex].GetCountsV();
    std::map<std::string, int> initialCounts;
    for(auto& sample:sampleDefs){
        initialCounts[sample.GetLabelS()] = sample.GetNEvents();
    }

    for(auto& sampleName : histogramCounts){

        double potScaling = potScalingMap[sampleName.first];

        std::ostringstream streamObjEff;
        streamObjEff <<  std::fixed << std::setprecision(2) << 100.*sampleName.second/initialCounts[sampleName.first];
        // Get string from out
        if(  cont == sampleDefs.size()-1 )
            texFile << potScaling*sampleName.second <<" ("<< streamObjEff.str() <<" \\%)";
        else
            texFile << potScaling*sampleName.second <<" ("<< streamObjEff.str() <<" \\%)"  << " & ";
        cont++;
    }


    texFile << "\\\\ \\hline" << std::endl;

    
    std::string newTableCaption = tableCaption;
    if(totalPOTNorm!=1){
        // total pot in scientific notation
        std::ostringstream streamObj;
        streamObj <<  std::scientific << std::setprecision(2) << totalPOTNorm;
        // Get string from out
        std::string totalPOTNormStr = streamObj.str();
        newTableCaption += " (Normalized to "+totalPOTNormStr+" POT)";
    }
    texFile << "\\end{tabular}" << std::endl;
    texFile << "\\caption{" << newTableCaption << "}" << std::endl;
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


//--------- Read the POT normalization
void ReadPOT(TFile *fFile, double totalPOTNorm, double& potScaling, double& potScalingSignal, std::string fTreeName = "originsAna/pottree"){
    
    // Read TreePOT
    TTree *fTreePOT = (TTree *)fFile->Get( fTreeName.c_str() );


    // -------- Get the accumulated POT
    double pot = 0;
    double averageintmode;
    bool inclusive;
    fTreePOT->SetBranchAddress("pot",&pot);
    fTreePOT->SetBranchAddress("averageintmode",&averageintmode);
    fTreePOT->SetBranchAddress("inclusive",&inclusive);
    double totalPOT = 0;
    double totalPOTSignal = 0;
    for (int i = 0; i < fTreePOT->GetEntries(); ++i) {
        fTreePOT->GetEntry(i);
        //std::cout<<"averageintmode: "<<averageintmode<<" inclusive: "<<inclusive<<" POT: "<<pot<<std::endl;
        if(averageintmode==0 && inclusive==0) totalPOTSignal+=pot;
        else if(pot<1e17) totalPOT+=pot;
    }
    
    potScaling = totalPOTNorm/totalPOT;
    potScalingSignal = totalPOTNorm/totalPOTSignal;
    std::cout<<" Accumulated POT: "<<totalPOT<<std::endl;
    std::cout<<" Accumulated POT Signal: "<<totalPOTSignal<<std::endl;
    std::cout<<" POT scaling factors: "<<potScaling<<" "<<potScalingSignal<<std::endl;

    return;   
}


// --------- Create output list with events that pass the cuts
void CreateHandScanList(TTree *fTree, TTree *fTreeHeader, TCut cut, std::vector<SampleDef> sampleDefs){

    //--------- Loop over the TreeHeader
    unsigned int fSubRunId;
    fTreeHeader->SetBranchAddress("SubRunID", &fSubRunId);
    unsigned int fRunId;
    fTreeHeader->SetBranchAddress("RunID", &fRunId);
    std::string *fLArInputFileName = new std::string;
    fTreeHeader->SetBranchAddress("InputFileName", &fLArInputFileName);
    // Create map of run-subrun filename
    std::map<std::string, std::vector<std::string>> subrunFilenameMap;

    std::map<std::string, int> fileCounterMap;

    for(size_t i=0; i<fTreeHeader->GetEntries(); ++i){
        fTreeHeader->GetEntry(i);
        
        // check the string includes the substring "Inclusive"
        if(fLArInputFileName->find("V0Lambda") != std::string::npos || fLArInputFileName->find("V0Overlay") != std::string::npos){
            continue;
        }
        
        
        std::string runSubrunLabel = std::to_string(fRunId) + ":" + std::to_string(fSubRunId);

        // if the subrun is not in the map, add it
        if(subrunFilenameMap.find(runSubrunLabel) == subrunFilenameMap.end()){
            subrunFilenameMap[runSubrunLabel] = std::vector<std::string>();
        }

        // add the filename to the vector
        subrunFilenameMap[runSubrunLabel].push_back(*fLArInputFileName);


        // if the file is not in the map, add it
        if(fileCounterMap.find(*fLArInputFileName) == fileCounterMap.end()){
            fileCounterMap[*fLArInputFileName] = 1;
        }
        else{
            fileCounterMap[*fLArInputFileName]++;
        }
          
    }

    // print counter map
    for(auto& fileCounter : fileCounterMap){
        std::cout << "File: " << fileCounter.first << " Count: " << fileCounter.second << std::endl;
    }
    // cout map size
    std::cout << "Map size: " << fileCounterMap.size() << std::endl;

    // Output 
    std::ofstream handScanEvents;
    handScanEvents.open("OutputPlots/handScanEvents.txt");

    // loop over signal types
    for (size_t k = 0; k < sampleDefs.size(); k++) {
        // Skip if it is the signal
        if (sampleDefs[k].IsSignal()) continue;

        std::cout << " Sample: " << sampleDefs[k].GetLabel() << std::endl;

        // Create the histogram
        TH3D* h3_all = new TH3D("h3_all", "3D Histogram", 100, 1, 101, 1100, 1, 1101, 1100, 1, 1101);

        // Create the cut
        TCut debugCut = cut && (TCut(sampleDefs[k].GetVar()));

        // Draw the histogram
        // Z axis is the run number, Y axis is event ID, and X axis is subrun ID
        fTree->Draw("EventID:SubrunID:RunID>>h3_all", debugCut);

        // Get the number of bins in the histogram
        int numBinsX = h3_all->GetNbinsX();
        int numBinsY = h3_all->GetNbinsY();
        int numBinsZ = h3_all->GetNbinsZ();

        for (int run = 1; run <= numBinsX; run++) {
            for (int sr = 1; sr <= numBinsY; sr++) {
                for (int e = 1; e <= numBinsZ; e++) {
                    double binValue = h3_all->GetBinContent(run, sr, e);
                    if (binValue > 0) {
                        std::cout << "RunNumber: " << run << " SubRunID: " << sr << " EventID: " << e << " BinValue: " << binValue << std::endl;
                        std::string runSubrunLabel = std::to_string(run) + ":" + std::to_string(sr);
                        for (size_t i = 0; i < subrunFilenameMap[runSubrunLabel].size(); i++){
                            std::cout << "   Filename: " << subrunFilenameMap[runSubrunLabel][i] << std::endl;
                            handScanEvents << run << ":" << sr << ":" << e << " " << subrunFilenameMap[runSubrunLabel][i] << std::endl;
                        }
                    }
                }
            }
        }
    }


    handScanEvents.close();

    return;

}

#endif