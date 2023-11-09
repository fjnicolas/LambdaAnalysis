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

LambdaAnaTree::LambdaAnaTree(TTree* tree)
: fTree(tree)
{
    InitializeTree();
}

void LambdaAnaTree::SetTree(TTree* tree)
{
    fTree = tree;
    InitializeTree();
}


void LambdaAnaTree::InitializeTree(){

    // Set branch addresses for event information
    fTree->Branch("EventID", &fEventID);
    fTree->Branch("SubrunID", &fSubrunID);
    fTree->Branch("RunID", &fRunID);
    fTree->Branch("InputFileName", &fInputFileName);
    fTree->Branch("SliceID", &fSliceID);

    // Set branch addresses for true variables
    fTree->Branch("IntMode", &fIntMode);
    fTree->Branch("IntType", &fIntType);
    fTree->Branch("IntCCNC", &fIntCCNC);
    fTree->Branch("IntNuPDG", &fIntNuPDG);
    fTree->Branch("IntNProtons", &fIntNProtons);
    fTree->Branch("IntNNeutrons", &fIntNNeutrons);
    fTree->Branch("IntNPi0", &fIntNPi0);
    fTree->Branch("IntNPip", &fIntNPip);
    fTree->Branch("IntNPim", &fIntNPim);
    fTree->Branch("IntNMuonP", &fIntNMuonP);
    fTree->Branch("IntNMuonM", &fIntNMuonM);
    fTree->Branch("IntNElectronP", &fIntNElectronP);
    fTree->Branch("IntNElectronM", &fIntNElectronM);
    fTree->Branch("IntNLambda", &fIntNLambda);

    // Set branch addresses for true vertex information
    fTree->Branch("NuvE", &fNuvE);
    fTree->Branch("NuvT", &fNuvT);
    fTree->Branch("NuvX", &fNuvX);
    fTree->Branch("NuvY", &fNuvY);
    fTree->Branch("NuvZ", &fNuvZ);
    fTree->Branch("TruthIsFiducial", &fTruthIsFiducial);

    // Set branch addresses for reco vertex information
    fTree->Branch("RecnuvX", &fRecnuvX);
    fTree->Branch("RecnuvY", &fRecnuvY);
    fTree->Branch("RecnuvZ", &fRecnuvZ);
    fTree->Branch("RecoIsFiducial", &fRecoIsFiducial);

    // Set branch addresses for FRANS PANDORA information
    fTree->Branch("FRANSScorePANDORA", &fFRANSScorePANDORA);
    fTree->Branch("ShowerEnergy", &fShowerEnergy);
    fTree->Branch("NShowers", &fNShowers);
    fTree->Branch("NShwTh75", &fNShwTh75);
    fTree->Branch("NShwTh100", &fNShwTh100);
    fTree->Branch("MainShowerEnergy", &fMainShowerEnergy);
    fTree->Branch("MainShowerScore", &fMainShowerScore);
    fTree->Branch("ShowerEnergyVect", &fShowerEnergyVect);
    fTree->Branch("ShowerScoreVect", &fShowerScoreVect);
    
    // Set branch addresses for FRANS PANDORA information
    fTree->Branch("NOrigins", &fNOrigins);
    fTree->Branch("NOriginsMult1", &fNOriginsMult1);
    fTree->Branch("NOriginsMult2", &fNOriginsMult2);
    fTree->Branch("NOriginsMultGT3", &fNOriginsMultGT3);
    fTree->Branch("NOriginsPairOneTwo", &fNOriginsPairOneTwo);

    // Set branch addresses for angle information
    fTree->Branch("NAngles", &fNAngles);
    fTree->Branch("AngleFRANSScore", &fAngleFRANSScore);
    fTree->Branch("AngleGap", &fAngleGap);
    fTree->Branch("AngleNHits", &fAngleNHits);
    fTree->Branch("AngleMainTrackNHits", &fAngleMainTrackNHits);
    fTree->Branch("AngleLongestIsMain", &fAngleLongestIsMain);
    fTree->Branch("AngleDecayContainedDiff", &fAngleDecayContainedDiff);

}

void LambdaAnaTree::ResetVars(){
    fEventID = -999;
    fSubrunID = -999;
    fRunID = -999;
    fInputFileName = "";
    fSliceID = -999;

    fIntMode = -999;
    fIntType = -999;
    fIntCCNC = -999;
    fIntNuPDG = -999;
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

    fRecnuvX = -999;
    fRecnuvY = -999;
    fRecnuvZ = -999;
    fRecoIsFiducial = false;

    fFRANSScorePANDORA = -999;
    fShowerEnergy = -999;
    fNShwTh75 = -999;
    fNShwTh100 = -999;
    fMainShowerEnergy = -999;
    fMainShowerScore = -999;
    fShowerEnergyVect.clear();
    fShowerScoreVect.clear();

    fNOrigins = -999;
    fNOriginsMult1 = -999;
    fNOriginsMult2 = -999;
    fNOriginsMultGT3 = -999;
    fNOriginsPairOneTwo = -999;

    fNAngles = -999;
    fAngleFRANSScore = -999;
    fAngleGap = -999;
    fAngleNHits = -999;
    fAngleMainTrackNHits = -999;
    fAngleLongestIsMain = false;
    fAngleDecayContainedDiff = -999;

}
