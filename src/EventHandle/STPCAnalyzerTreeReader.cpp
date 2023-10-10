////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "STPCAnalyzerTreeReader.h"

MyTPCTreeReader::~MyTPCTreeReader(){
    // Clean up when the class is destroyed
    delete file;
}

int MyTPCTreeReader::NEntries(){return fNEntries;}

MyTPCTreeReader::MyTPCTreeReader(TString fileName, std::string treeName) {
    // Open the ROOT file
    file = new TFile(fileName);

    // Get the TTree from the file (assuming the TTree is named "myTree")
    tree = (TTree *)file->Get(treeName.c_str());

    fNEntries = tree->GetEntries();

    // Set branch addresses for event information
    tree->SetBranchAddress("EventID", &eventID);
    tree->SetBranchAddress("SubRunID", &subrunID);
    tree->SetBranchAddress("RunID", &runID);

    // Set branch addresses for true vertex information
    tree->SetBranchAddress("TrueVEnergy", &nuvE);
    tree->SetBranchAddress("TrueVt", &nuvT);
    tree->SetBranchAddress("TrueVx", &nuvX);
    tree->SetBranchAddress("TrueVy", &nuvY);
    tree->SetBranchAddress("TrueVz", &nuvZ);
    tree->SetBranchAddress("TrueVU", &nuvU);
    tree->SetBranchAddress("TrueVV", &nuvV);
    tree->SetBranchAddress("TrueVC", &nuvC);
    tree->SetBranchAddress("TrueVTimeTick", &nuvTimeTick);

    // Set branch addresses for reco vertex information
    tree->SetBranchAddress("RecoVx", &recnuvX);
    tree->SetBranchAddress("RecoVy", &recnuvY);
    tree->SetBranchAddress("RecoVz", &recnuvZ);
    tree->SetBranchAddress("RecoVU", &recnuvU);
    tree->SetBranchAddress("RecoVV", &recnuvV);
    tree->SetBranchAddress("RecoVC", &recnuvC);
    tree->SetBranchAddress("RecoVTimeTick", &recnuvTimeTick);

    // True PDGs
    tree->SetBranchAddress("TruePrimariesPDG", &truePrimeriesPDG);
    
    // Set branch addresses for hits information
    tree->SetBranchAddress("HitsView", &hitsView);
    tree->SetBranchAddress("HitsIntegral", &hitsIntegral);
    tree->SetBranchAddress("HitsPeakTime", &hitsPeakTime);
    tree->SetBranchAddress("HitsChannel", &hitsChannel);
    tree->SetBranchAddress("HitsRMS", &hitsRMS);
    tree->SetBranchAddress("HitsStartT", &hitsStartT);
    tree->SetBranchAddress("HitsEndT", &hitsEndT);
    tree->SetBranchAddress("HitsChi2", &hitsChi2);
    tree->SetBranchAddress("HitsNDF", &hitsNDF);
    tree->SetBranchAddress("HitsClusterID", &hitsClusterID);
}


bool MyTPCTreeReader::GetEntry(int entry) {
    if (entry < 0 || entry >= tree->GetEntries())
        return false;
    
    tree->GetEntry(entry);
    return true;
}