#ifndef CUTEFFICIENCIES_LATEXINTERFACE_H
#define CUTEFFICIENCIES_LATEXINTERFACE_H

#define WriteRelativeEfficiencies 0


//--------- Read the POT normalization
void ReadPOT(TFile *fFile, double totalPOTNorm, double& potScaling, double& potScalingSignal, std::string fTreeName = "originsAna/pottree"){
    
    // Read TreePOT
    std::cout<<"Reading POT tree "<<fTreeName<<std::endl;
    TTree *fTreePOT = (TTree *)fFile->Get( fTreeName.c_str() );
    if (!fTreePOT) {
        std::cerr << "Failed to read the POT tree." << std::endl;
        return;
    }
    
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
        
        // check the string includes the substring "V0Lambda"
        //if(fLArInputFileName->find("V0Lambda") != std::string::npos || fLArInputFileName->find("V0Overlay") != std::string::npos) continue;
        

        if(fLArInputFileName->find("Inclusive") == std::string::npos) continue;
        
        
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
    handScanEvents.open("OutputPlots/HandScanEvents.txt");

    // loop over signal types
    for (size_t k = 0; k < sampleDefs.size(); k++) {
        // Skip if it is the signal
        if (sampleDefs[k].IsSignal()) continue;

        std::cout << " Sample: " << sampleDefs[k].GetLabel() << std::endl;

        // Create the histogram
        TH3D* h3_all = new TH3D("h3_all", "3D Histogram", 100, 1, 101, 2100, 1, 2101, 2100, 1, 2101);

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

//--------- Generate and compile the LaTeX table
void GenerateAndCompileTeXTable(
    const std::vector<PlotDef> plotDefs,
    const std::vector<SampleDef> sampleDefs,
    const std::vector<std::vector<int>> matrixCounter,
    const std::string fileName,
    const std::string tableCaption,
    const std::string dirLabel,
    double potNorm=1,
    double potNormSignal=1,
    double totalPOTNorm=1,
    int denominatorIndex=0
) {
    //--- Open a file stream for writing the .tex file
    std::ofstream texFile(fileName+".tex");

    //--- Check if the file is open
    if (!texFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return;
    }

    //--- Denominator weights
    std::vector<int> counts0 = matrixCounter[denominatorIndex];

    //--- POT scaling map weights
    std::vector<double> potScalingMap;
    for(auto& sampleName : sampleDefs){
        if(sampleName.GetLabelS()=="Signal")
            potScalingMap.push_back( potNormSignal );
        else
            potScalingMap.push_back( potNorm );
    }

    //--- Write the LaTeX document preamble
    texFile << "\\documentclass{article}" << std::endl;
    texFile << "\\usepackage{graphicx}" << std::endl;
    texFile << "\\usepackage{pdflscape}"<<std::endl;
    texFile << "\\begin{document}" << std::endl;
    texFile << "\\begin{landscape}" << std::endl;


    //--- Write the LaTeX code for the table
    texFile << "\\begin{table}[h]" << std::endl;
    texFile << "\\centering" << std::endl;
    texFile << "\\begin{tabular}{|c|";
    for(size_t k = 0; k<sampleDefs.size(); k++)
        texFile<<"c|";
    texFile << "}" << std::endl;
    texFile << "\\hline" << std::endl;

    

    //--- Write the header row
    texFile << " Cut & ";
    int cont=0;  
    for(auto& sampleName : sampleDefs){
        if( cont == sampleDefs.size()-1 )
            texFile << sampleName.GetLatexLabel();
        else
            texFile << sampleName.GetLatexLabel() << " & ";
        
        cont++;
    }
    texFile << "\\\\ \\hline" << std::endl;


    // Write the data rows
    // last stored index
    size_t lastStoredIndex = 0;
    for (size_t i = 0; i < matrixCounter.size(); ++i) {

        texFile << "$ {\\rm " << plotDefs[i].GetCutLabelS() << "}$" << " & ";

        std::vector<int> counts = matrixCounter[i];
        std::vector<int> countsPre = matrixCounter[lastStoredIndex];

        cont=0;
        for (size_t j = 0; j < matrixCounter[i].size(); ++j) {
        //for(auto& sampleName : histogramCounts){

            double potScaling = potScalingMap[j];
            std::string sampleName = sampleDefs[j].GetLabelS();
            int counts = matrixCounter[i][j];
            
            std::ostringstream streamObjEff;
            streamObjEff << " (" << std::fixed << std::setprecision(3) << 100.*counts/counts0[j];
            if(WriteRelativeEfficiencies==1){
                streamObjEff << " \\%) \\ $\\epsilon_{r}$= " << std::setprecision(0) << 100.*counts/countsPre[j] << "\\%";
            }
            else{
                streamObjEff << " \\%)";
            }

            // Get string from out
            if(  cont == sampleDefs.size()-1 )
                texFile << potScaling*counts << streamObjEff.str();
            else
                texFile << potScaling*counts << streamObjEff.str() << " & ";
            cont++;
        }

        texFile << "\\\\ \\hline" << std::endl;

        lastStoredIndex = i;
    }
    


    /*
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
    */

    // --- Write the caption and end the table    
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


    //--- Close the LaTeX document
    texFile << "\\end{document}" << std::endl;

    //--- Close the file
    texFile.close();

    //--- Compile the .tex file
    std::string compileCommand = "pdflatex " + fileName;
    int compileResult = std::system(compileCommand.c_str());
    if (compileResult != 0) {
        std::cerr << "Failed to compile the .tex file." << std::endl;
    }
    gSystem->Exec( ("mv "+fileName+".pdf "+dirLabel).c_str() );
    gSystem->Exec( ("rm "+fileName+".{aux,log,tex,pdf}").c_str() );

    return;

}

#endif