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

#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleHits.h"
#else
#include "TPCSimpleHits.h"
#endif



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
                    std::vector<double> *_Chi2);



// Class to read the TPCAnalyzer TTree
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

    int intMode;
    int intNLambda;
    
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