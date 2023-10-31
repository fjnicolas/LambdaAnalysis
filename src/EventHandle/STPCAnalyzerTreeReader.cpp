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


// Class to get the hits in a specific view
std::vector<SHit> GetHitsInView(
    int view,
    std::vector<int> *_X,
    std::vector<double> *_Y,
    std::vector<double> *_Int,
    std::vector<double> *_Wi,
    std::vector<double> *_ST,
    std::vector<double> *_ET,
    std::vector<int> *_View,
    std::vector<double> *_Chi2)
{

    // set variables
    std::vector<SHit> hitList;
    int nTotalHits = _X->size();

    // loop over the hits
    for (int i = 0; i < nTotalHits; i++) {

        // filter channels for the view        
        if ( _View->at(i)==view ) {
            SHit hit(-1, _X->at(i), _Y->at(i), _Wi->at(i), _Int->at(i), _ST->at(i), _ET->at(i), _Chi2->at(i));
            hitList.push_back(hit);
        }
    }

    return hitList;

}



// Class to read the TPCAnalyzer TTree
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

    tree->SetBranchAddress("IntMode", &intMode);
    tree->SetBranchAddress("IntNLambda", &intNLambda);

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



