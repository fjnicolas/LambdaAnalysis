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
    fTree->Branch("SliceID", &fSliceID);

    // Set branch addresses for true variables
    fTree->Branch("IntMode", &fIntMode);
    fTree->Branch("IntCCNC", &fIntCCNC);
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

    // Set branch addresses for reco vertex information
    fTree->Branch("RecnuvX", &fRecnuvX);
    fTree->Branch("RecnuvY", &fRecnuvY);
    fTree->Branch("RecnuvZ", &fRecnuvZ);

    // Set branch addresses for FRANS PANDORA information
    fTree->Branch("FRANSScorePANDORA", &fFRANSScorePANDORA);
    
    // Set branch addresses for FRANS PANDORA information
    fTree->Branch("NOrigins", &fNOrigins);
    fTree->Branch("NOriginsMult1", &fNOriginsMult1);
    fTree->Branch("NOriginsMult2", &fNOriginsMult2);
    fTree->Branch("NOriginsMultGT3", &fNOriginsMultGT3);
    fTree->Branch("NOriginsPairOneTwo", &fNOriginsPairOneTwo);

    // Set branch addresses for angle information
    fTree->Branch("NAngles", &fNAngles);
    fTree->Branch("AngleFRANSScore", &fAngleFRANSScore);
    fTree->Branch("AngleNHits", &fAngleNHits);
    fTree->Branch("AngleMainTrackNHits", &fAngleMainTrackNHits);
    fTree->Branch("AngleLongestIsMain", &fAngleLongestIsMain);

}
