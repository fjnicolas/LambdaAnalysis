////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_STPCANALYZERTREEREADER_H
#define TPC_STPCANALYZERTREEREADER_H

#include <vector>

// ROOT includess
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


class MyTPCTreeReader {
private:
    TFile *file;    // Pointer to the ROOT file containing the TTree
    TTree *tree;    // Pointer to the TTree object

    int fNEntries;
    
public:
    // Event information
    int eventID;
    int subrunID;
    int runID;
    
    // True Vertex information
    double nuvE, nuvT, nuvX, nuvY, nuvZ;
    int nuvU, nuvV, nuvC, nuvTimeTick;
    
    // Reco Vertex information
    double recnuvX;
    double recnuvY;
    double recnuvZ;
    int recnuvU;
    int recnuvV;
    int recnuvC;
    int recnuvTimeTick;

    // MC PDGs
    std::vector<int> *truePrimeriesPDG = new std::vector<int>;

    // Hits information
    std::vector<int> *hitsView = new std::vector<int>;
    std::vector<double> *hitsIntegral = new std::vector<double>;
    std::vector<double> *hitsPeakTime = new std::vector<double>;
    std::vector<int> *hitsChannel = new std::vector<int>;
    std::vector<double> *hitsRMS = new std::vector<double>;
    std::vector<double> *hitsStartT = new std::vector<double>;
    std::vector<double> *hitsEndT = new std::vector<double>;
    std::vector<double> *hitsChi2 = new std::vector<double>;
    std::vector<double> *hitsNDF = new std::vector<double>;
    std::vector<int> *hitsClusterID = new std::vector<int>;

    MyTPCTreeReader(TString fileName, std::string treeName);

    ~MyTPCTreeReader();

    bool GetEntry(int entry);

    int NEntries();

};



#endif // TPC_LINES_PARAMETERS_H