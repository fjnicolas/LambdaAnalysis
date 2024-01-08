////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of LambdaAnaTree
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "LambdaAnaTree.h"

LambdaAnaTree::LambdaAnaTree()
: fTree(nullptr)
{}

LambdaAnaTree::LambdaAnaTree(TTree* tree, bool readMode)
: fTree(tree),
  fInputFileNameRead(new std::string)
{
    InitializeTree(readMode);
}

void LambdaAnaTree::SetTree(TTree* tree, bool readMode)
{
    fTree = tree;
    InitializeTree(readMode);
}

template <typename T>
void LambdaAnaTree::SetBranch(const char* name, T* variable, bool setAddress) {
    if (setAddress) {
        fTree->SetBranchAddress(name, variable);
        std::cout << "Setting branch address for " << name << std::endl;
    } else {
        fTree->Branch(name, variable);
    }
}

// instance types
template void LambdaAnaTree::SetBranch(const char* name, int* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, float* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, double* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, bool* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, std::string* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, std::vector<float>* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, std::vector<double>* variable, bool setAddress);


void LambdaAnaTree::InitializeTree(bool readMode){
    // Set branch addresses for event information
    SetBranch("EventID", &fEventID, readMode);
    SetBranch("SubrunID", &fSubrunID, readMode);
    SetBranch("RunID", &fRunID, readMode);
    SetBranch("SliceID", &fSliceID, readMode);

    // Set branch addresses for true variables
    SetBranch("IntOrigin", &fIntOrigin, readMode);
    SetBranch("IntMode", &fIntMode, readMode);
    SetBranch("IntDirt", &fIntDirt, readMode);
    SetBranch("IntType", &fIntType, readMode);
    SetBranch("IntCCNC", &fIntCCNC, readMode);
    SetBranch("IntNuPDG", &fIntNuPDG, readMode);
    SetBranch("IntNuEnergy", &fIntNuEnergy, readMode);
    SetBranch("IntNProtons", &fIntNProtons, readMode);
    SetBranch("IntNNeutrons", &fIntNNeutrons, readMode);
    SetBranch("IntNPi0", &fIntNPi0, readMode);
    SetBranch("IntNPip", &fIntNPip, readMode);
    SetBranch("IntNPim", &fIntNPim, readMode);
    SetBranch("IntNMuonP", &fIntNMuonP, readMode);
    SetBranch("IntNMuonM", &fIntNMuonM, readMode);
    SetBranch("IntNElectronP", &fIntNElectronP, readMode);
    SetBranch("IntNElectronM", &fIntNElectronM, readMode);
    SetBranch("IntNLambda", &fIntNLambda, readMode);

    // Set branch addresses for true vertex information
    SetBranch("NuvE", &fNuvE, readMode);
    SetBranch("NuvT", &fNuvT, readMode);
    SetBranch("NuvX", &fNuvX, readMode);
    SetBranch("NuvY", &fNuvY, readMode);
    SetBranch("NuvZ", &fNuvZ, readMode);
    SetBranch("TruthIsFiducial", &fTruthIsFiducial, readMode);
    SetBranch("TruthIsAV", &fTruthIsAV, readMode);

    // Set branch addresses for slice information
    SetBranch("SliceCompleteness", &fSliceCompleteness, readMode);
    SetBranch("SlicePurity", &fSlicePurity, readMode);

    // Set branch addresses for lambda true information
    SetBranch("Gap", &fGap, readMode);
    SetBranch("LambdaKE", &fLambdaKE, readMode);
    SetBranch("ProtonKE", &fProtonKE, readMode);
    SetBranch("PionKE", &fPionKE, readMode);

    // Set branch addresses for cosmic rejection
    SetBranch("CRUMBSScore", &fCRUMBSScore, readMode);

    // Set branch addresses for reco vertex information
    SetBranch("RecnuvX", &fRecnuvX, readMode);
    SetBranch("RecnuvY", &fRecnuvY, readMode);
    SetBranch("RecnuvZ", &fRecnuvZ, readMode);
    SetBranch("RecoIsFiducial", &fRecoIsFiducial, readMode);

    // Set branch addresses for FRANS PANDORA information
    SetBranch("FRANSScorePANDORA", &fFRANSScorePANDORA, readMode);
    SetBranch("NVertexTracks", &fNVertexTracks, readMode);
    SetBranch("ShowerEnergy", &fShowerEnergy, readMode);
    SetBranch("NShowers", &fNShowers, readMode);
    SetBranch("NShowerHits", &fNShowerHits, readMode);
    SetBranch("NShwTh75", &fNShwTh75, readMode);
    SetBranch("NShwTh100", &fNShwTh100, readMode);
    SetBranch("MainShowerEnergy", &fMainShowerEnergy, readMode);
    SetBranch("MainShowerScore", &fMainShowerScore, readMode);
    if(!readMode){
        SetBranch("ShowerEnergyVect", &fShowerEnergyVect, readMode);
        SetBranch("ShowerScoreVect", &fShowerScoreVect, readMode);
        SetBranch("ShowerNHitsVect", &fShowerNHitsVect, readMode);
        SetBranch("ShowerOutROINHitsVect", &fShowerOutROINHitsVect, readMode);
    }
    
    // Set branch addresses for # origins information
    SetBranch("NOrigins", &fNOrigins, readMode);
    SetBranch("NOriginsMult1", &fNOriginsMult1, readMode);
    SetBranch("NOriginsMult2", &fNOriginsMult2, readMode);
    SetBranch("NOriginsMultGT3", &fNOriginsMultGT3, readMode);
    SetBranch("NOriginsPairOneTwo", &fNOriginsPairOneTwo, readMode);

    // Set branch addresses for # unassociated origins information
    SetBranch("NUnOrigins", &fNUnOrigins, readMode);
    SetBranch("NUnOriginsMult1", &fNUnOriginsMult1, readMode);
    SetBranch("NUnOriginsMult2", &fNUnOriginsMult2, readMode);
    SetBranch("NUnOriginsMultGT3", &fNUnOriginsMultGT3, readMode);

    // Set branch addresses for angle information
    SetBranch("NAngles", &fNAngles, readMode);
    SetBranch("AngleFRANSScore", &fAngleFRANSScore, readMode);
    SetBranch("AngleGap", &fAngleGap, readMode);
    SetBranch("AngleNHits", &fAngleNHits, readMode);
    SetBranch("AngleNHitsTrack1", &fAngleNHitsTrack1, readMode);
    SetBranch("AngleNHitsTrack2", &fAngleNHitsTrack2, readMode);
    SetBranch("AngleMinNHits", &fAngleMinNHits, readMode);
    SetBranch("AngleNHitsMainTrack", &fAngleNHitsMainTrack, readMode);
    SetBranch("AngleLengthTrack1", &fAngleLengthTrack1, readMode);
    SetBranch("AngleLengthTrack2", &fAngleLengthTrack2, readMode);
    SetBranch("AngleLengthMainTrack", &fAngleLengthMainTrack, readMode);
    SetBranch("AngleLongestIsMain", &fAngleLongestIsMain, readMode);
    SetBranch("AngleDecayContainedDiff", &fAngleDecayContainedDiff, readMode);
    SetBranch("AngleCoveredArea", &fAngleCoveredArea, readMode);
    SetBranch("AngleDirtHits", &fAngleDirtHits, readMode);
    SetBranch("AngleDirtHitsRatio", &fAngleDirtHitsRatio, readMode);
    SetBranch("AngleDirtHitsWires", &fAngleDirtHitsWires, readMode);
    SetBranch("AngleDirtHitsWiresRatio", &fAngleDirtHitsWiresRatio, readMode);
    SetBranch("AngleOpeningAngle", &fAngleOpeningAngle, readMode);

    SetBranch("NFreeHits", &fNFreeHits, readMode);
    SetBranch("NUnassociatedHits", &fNUnassociatedHits, readMode);

    // Set branch addresses for charge ratio information
    SetBranch("AnglePassFit", &fAnglePassFit, readMode);
    SetBranch("AngleTwoLinesChi2", &fAngleTwoLinesChi2, readMode);
    SetBranch("AngleNVertexHits", &fAngleNVertexHits, readMode);
    SetBranch("AngleNBulkHits", &fAngleNBulkHits, readMode);
    SetBranch("AngleVertexHitIntegralRatio", &fAngleVertexHitIntegralRatio, readMode);
    SetBranch("AngleVertexHitIntegralDifference", &fAngleVertexHitIntegralDifference, readMode);
    SetBranch("AngleVertexHitIntegralRelativeDifference", &fAngleVertexHitIntegralRelativeDifference, readMode);
    SetBranch("AngleTrackLength1", &fAngleTrackLength1, readMode);
    SetBranch("AngleTrackLength2", &fAngleTrackLength2, readMode);
    SetBranch("AngleTrackLengthRatio", &fAngleTrackLengthRatio, readMode);
    SetBranch("AngleResidualRange1RMS", &fAngleResidualRange1RMS, readMode);
    SetBranch("AngleResidualRange2RMS", &fAngleResidualRange2RMS, readMode);
    SetBranch("AngleChargeRatioAverage", &fAngleChargeRatioAverage, readMode);
    SetBranch("AngleChargeDifferenceAverage", &fAngleChargeDifferenceAverage, readMode);
    SetBranch("AngleChargeRelativeDifferenceAverage", &fAngleChargeRelativeDifferenceAverage, readMode);
    SetBranch("AnglePassChargeFit", &fAnglePassChargeFit, readMode);
    SetBranch("AngleBandOverlap", &fAngleBandOverlap, readMode);
    SetBranch("AngleChargeRatioFit", &fAngleChargeRatioFit, readMode);
    SetBranch("AngleChargeDifferenceFit", &fAngleChargeDifferenceFit, readMode);
    SetBranch("AngleChargeRelativeDifferenceFit", &fAngleChargeRelativeDifferenceFit, readMode);


    if(readMode==false)
        SetBranch("InputFileName", &fInputFileName, readMode);
    else{
        SetBranch("InputFileName", &fInputFileNameRead, readMode);
    }

}

void LambdaAnaTree::ResetVars(){
    fEventID = -999;
    fSubrunID = -999;
    fRunID = -999;
    fInputFileName = "";
    fSliceID = -999;

    fIntOrigin = -999;
    fIntDirt = false;
    fIntMode = -999;
    fIntType = -999;
    fIntCCNC = -999;
    fIntNuPDG = -999;
    fIntNuEnergy = -999;
    fIntNProtons = -999;
    fIntNNeutrons = -999;
    fIntNPi0 = -999;
    fIntNPip = -999;
    fIntNPim = -999;
    fIntNMuonP = -999;
    fIntNMuonM = -999;
    fIntNElectronP = -999;
    fIntNElectronM = -999;
    fIntNLambda = -999;

    fNuvE = -999;
    fNuvT = -999;
    fNuvX = -999;
    fNuvY = -999;
    fNuvZ = -999;
    fTruthIsFiducial = false;
    fTruthIsAV = false;

    fSliceCompleteness = -999;
    fSlicePurity = -999;

    fGap = -999;
    fLambdaKE = -999;
    fProtonKE = -999;
    fPionKE = -999;

    fCRUMBSScore = -999;

    fRecnuvX = -999;
    fRecnuvY = -999;
    fRecnuvZ = -999;
    fRecoIsFiducial = false;

    fFRANSScorePANDORA = -999;
    fNVertexTracks = -999;
    fShowerEnergy = -999;
    fNShwTh75 = -999;
    fNShwTh100 = -999;
    fNShowers = -999;
    fNShowerHits = -999;
    fNShowersOutROI = -999;
    fMainShowerEnergy = -999;
    fMainShowerScore = -999;
    fShowerEnergyVect.clear();
    fShowerScoreVect.clear();
    fShowerNHitsVect.clear();
    fShowerOutROINHitsVect.clear();

    fNOrigins = -999;
    fNOriginsMult1 = -999;
    fNOriginsMult2 = -999;
    fNOriginsMultGT3 = -999;
    fNOriginsPairOneTwo = -999;

    fNUnOrigins = -999;
    fNUnOriginsMult1 = -999;
    fNUnOriginsMult2 = -999;
    fNUnOriginsMultGT3 = -999;

    fNAngles = -999;
    fAngleFRANSScore = -999;
    fAngleGap = -999;
    fAngleNHits = -999;
    fAngleNHitsTrack1 = -999;
    fAngleNHitsTrack2 = -999;
    fAngleMinNHits = -999;
    fAngleNHitsMainTrack = -999;
    fAngleLengthTrack1 = -999;
    fAngleLengthTrack2 = -999;
    fAngleLengthMainTrack = -999;
    fAngleLongestIsMain = false;
    fAngleDecayContainedDiff = -999;
    fAngleCoveredArea = -999;
    fAngleDirtHits = -999;
    fAngleDirtHitsRatio = -999;
    fAngleDirtHitsWires = -999;
    fAngleDirtHitsWiresRatio = -999;
    fAngleOpeningAngle = -999;

    fAnglePassFit = false;
    fAngleTwoLinesChi2 = -999;
    fAngleNVertexHits = -999;
    fAngleNBulkHits = -999;
    fAngleVertexHitIntegralRatio = -999;
    fAngleVertexHitIntegralDifference = -999;
    fAngleVertexHitIntegralRelativeDifference = -999;
    fAngleTrackLength1 = -999;
    fAngleTrackLength2 = -999;
    fAngleTrackLengthRatio = -999;
    fAngleResidualRange1RMS = -999;
    fAngleResidualRange2RMS = -999;
    fAngleChargeRatioAverage = -999;
    fAngleChargeDifferenceAverage = -999;
    fAngleChargeRelativeDifferenceAverage = -999;
    fAnglePassChargeFit = false;
    fAngleBandOverlap = -999;
    fAngleChargeRatioFit = -999;
    fAngleChargeDifferenceFit = -999;
    fAngleChargeRelativeDifferenceFit = -999;


    fNFreeHits = -999;
    fNUnassociatedHits = -999;

}


void LambdaAnaTree::GetEntry(int i){
    fTree->GetEntry(i);
}