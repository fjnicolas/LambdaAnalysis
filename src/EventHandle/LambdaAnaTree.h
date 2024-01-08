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

// C++ includes
#include <iostream>
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
    LambdaAnaTree(TTree* tree, bool readMode=0);
    void SetTree(TTree* tree, bool readMode=0);

    template <typename T> void SetBranch(const char* name, T* variable, bool setAddress);

    void InitializeTree(bool readMode);

    // Event information
    int fEventID;
    int fSubrunID;
    int fRunID;
    std::string fInputFileName;
    std::string *fInputFileNameRead;
    int fSliceID;
    
    // True variables
    int fIntOrigin;
    bool fIntDirt;
    int fIntMode;
    int fIntType;
    int fIntCCNC;
    int fIntNuPDG;
    double fIntNuEnergy;
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
    bool fTruthIsAV;

    // Slice truth quality
    double fSliceCompleteness;
    double fSlicePurity;

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
    int fNVertexTracks;
    double fShowerEnergy;
    int fNShowers;
    int fNShowerHits;
    int fNShwTh75;
    int fNShwTh100;
    int fNShowersOutROI;
    double fMainShowerEnergy;
    double fMainShowerScore;
    std::vector<double> fShowerEnergyVect;
    std::vector<double> fShowerScoreVect;
    std::vector<int> fShowerNHitsVect;
    std::vector<int> fShowerOutROINHitsVect;


    // Reco origins information
    int fNOrigins;
    int fNOriginsMult1;
    int fNOriginsMult2;
    int fNOriginsMultGT3;
    int fNOriginsPairOneTwo;

    // Reco origins information (after V selection)
    int fNUnOrigins;
    int fNUnOriginsMult1;
    int fNUnOriginsMult2;
    int fNUnOriginsMultGT3;

    // Angle information
    int fNAngles;
    double fAngleFRANSScore;
    double fAngleGap;
    int fAngleNHits;
    int fAngleNHitsTrack1;
    int fAngleNHitsTrack2;
    int fAngleMinNHits;
    int fAngleNHitsMainTrack;
    double fAngleLengthTrack1;
    double fAngleLengthTrack2;
    double fAngleLengthMainTrack;
    double fAngleDecayContainedDiff;
    bool fAngleLongestIsMain;
    double fAngleCoveredArea;
    int fAngleDirtHits;
    double fAngleDirtHitsRatio;
    int fAngleDirtHitsWires;
    double fAngleDirtHitsWiresRatio;
    double fAngleOpeningAngle;

    int fNFreeHits;
    int fNUnassociatedHits;

    // Calorimetry
    bool fAnglePassFit;
    double fAngleTwoLinesChi2;

    int fAngleNVertexHits;
    int fAngleNBulkHits;
    double fAngleVertexHitIntegralRatio;
    double fAngleVertexHitIntegralDifference;
    double fAngleVertexHitIntegralRelativeDifference;

    double fAngleTrackLength1;
    double fAngleTrackLength2;
    double fAngleTrackLengthRatio;
    double fAngleResidualRange1RMS;
    double fAngleResidualRange2RMS;

    double fAngleChargeRatioAverage;
    double fAngleChargeDifferenceAverage;
    double fAngleChargeRelativeDifferenceAverage;

    bool fAnglePassChargeFit;
    double fAngleBandOverlap;
    double fAngleChargeRatioFit;
    double fAngleChargeDifferenceFit;
    double fAngleChargeRelativeDifferenceFit;



    void FillTree() { fTree->Fill(); }
    void WriteTree() { fTree->Write(); }
    void PrintTree() { fTree->Print();}

    int GetEntries() { return fTree->GetEntries(); }
    void GetEntry(int i);

    void ResetVars();

};



#endif // TPC_LINES_PARAMETERS_H