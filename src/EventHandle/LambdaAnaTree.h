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
    std::string fInputFileName;
    int fSliceID;
    
    // True variables
    int fIntOrigin;
    bool fIntDirt;
    int fIntMode;
    int fIntType;
    int fIntCCNC;
    int fIntNuPDG;
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

    // Lambda decay true information
    double fGap;
    double fLambdaKE;
    double fProtonKE;
    double fPionKE;

    // Cosmic rejection
    double fCRUMBSScore;

    // Reco Vertex information
    double fRecnuvX;
    double fRecnuvY;
    double fRecnuvZ;
    bool fRecoIsFiducial;

    // PANDORA information
    double fFRANSScorePANDORA;
    double fShowerEnergy;
    int fNShowers;
    int fNShwTh75;
    int fNShwTh100;
    double fMainShowerEnergy;
    double fMainShowerScore;
    std::vector<double> fShowerEnergyVect;
    std::vector<double> fShowerScoreVect;


    // Reco origins information
    int fNOrigins;
    int fNOriginsMult1;
    int fNOriginsMult2;
    int fNOriginsMultGT3;
    int fNOriginsPairOneTwo;

    // Angle information
    int fNAngles;
    double fAngleFRANSScore;
    double fAngleGap;
    int fAngleNHits;
    int fAngleNHitsTrack1;
    int fAngleNHitsTrack2;
    int fAngleNHitsMainTrack;
    double fAngleLengthTrack1;
    double fAngleLengthTrack2;
    double fAngleLengthMainTrack;
    double fAngleDecayContainedDiff;
    bool fAngleLongestIsMain;
    int fNUnassociatedHits;


    void FillTree() { fTree->Fill(); }
    void WriteTree() { fTree->Write(); }

    void ResetVars();

};



#endif // TPC_LINES_PARAMETERS_H