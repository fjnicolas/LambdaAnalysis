////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of LambdaAnaTree
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LAMBDAANATREE_H
#define TPC_LAMBDAANATREE_H

#include <vector>

// ROOT includess
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


// Class to read the TPCAnalyzer TTree
class LambdaAnaTree {
private:
    TTree *fTree;    // Pointer to the TTree object
    
public:

    // Constructors
    LambdaAnaTree();
    LambdaAnaTree(TTree* tree);
    void SetTree(TTree* tree);

    void InitializeTree();


    // Event information
    int fEventID;
    int fSubrunID;
    int fRunID;
    int fSliceID;
    
    // True variables
    int fIntMode;
    int fIntCCNC;
    int fIntNProtons;
    int fIntNNeutrons;
    int fIntNPi0;
    int fIntNPip;
    int fIntNPim;
    int fIntNMuonP;
    int fIntNMuonM;
    int fIntNElectronP;
    int fIntNElectronM;
    int fIntNLambda;

    // True Vertex information
    double fNuvE, fNuvT, fNuvX, fNuvY, fNuvZ;
    bool fTruthIsFiducial;

    // Reco Vertex information
    double fRecnuvX;
    double fRecnuvY;
    double fRecnuvZ;
    bool fRecoIsFiducial;

    // FRANS PANDORA information
    double fFRANSScorePANDORA;

    // Reco origins information
    int fNOrigins;
    int fNOriginsMult1;
    int fNOriginsMult2;
    int fNOriginsMultGT3;
    int fNOriginsPairOneTwo;

    // Angle information
    int fNAngles;
    double fAngleFRANSScore;
    int fAngleNHits;
    int fAngleMainTrackNHits;
    bool fAngleLongestIsMain;


    void FillTree() { fTree->Fill(); }
    void WriteTree() { fTree->Write(); }

    void ResetVars();

};



#endif // TPC_LINES_PARAMETERS_H