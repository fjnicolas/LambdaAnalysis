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


// Class to read the TPCAnalyzer TTree
LambdaAnaTree::~LambdaAnaTree(){
}

LambdaAnaTree::LambdaAnaTree(std::string treeName):
    tree ( new TTree(treeName.c_str(), "Example ROOT Tree") ) { 

    // Set branch addresses for event information
    tree->Branch("EventID", &fEventID);
    tree->Branch("SubrunID", &fSubrunID);
    tree->Branch("RunID", &fRunID);
    tree->Branch("SliceID", &fSliceID);

    // Set branch addresses for true variables
    tree->Branch("IntMode", &fIntMode);
    tree->Branch("IntCCNC", &fIntCCNC);
    tree->Branch("IntNProtons", &fIntNProtons);
    tree->Branch("IntNNeutrons", &fIntNNeutrons);
    tree->Branch("IntNPi0", &fIntNPi0);
    tree->Branch("IntNPip", &fIntNPip);
    tree->Branch("IntNPim", &fIntNPim);
    tree->Branch("IntNMuonP", &fIntNMuonP);
    tree->Branch("IntNMuonM", &fIntNMuonM);
    tree->Branch("IntNElectronP", &fIntNElectronP);
    tree->Branch("IntNElectronM", &fIntNElectronM);
    tree->Branch("IntNLambda", &fIntNLambda);

    // Set branch addresses for true vertex information
    tree->Branch("NuvE", &fNuvE);
    tree->Branch("NuvT", &fNuvT);
    tree->Branch("NuvX", &fNuvX);
    tree->Branch("NuvY", &fNuvY);
    tree->Branch("NuvZ", &fNuvZ);

    // Set branch addresses for reco vertex information
    tree->Branch("RecnuvX", &fRecnuvX);
    tree->Branch("RecnuvY", &fRecnuvY);
    tree->Branch("RecnuvZ", &fRecnuvZ);

    // Set branch addresses for FRANS PANDORA information
    tree->Branch("FRANSScorePANDORA", &fFRANSScorePANDORA);
    
    // Set branch addresses for FRANS PANDORA information
    tree->Branch("NOrigins", &fNOrigins);
    tree->Branch("NOriginsMult1", &fNOriginsMult1);
    tree->Branch("NOriginsMult2", &fNOriginsMult2);
    tree->Branch("NOriginsMultGT3", &fNOriginsMultGT3);
    tree->Branch("NOriginsPairOneTwo", &fNOriginsPairOneTwo);

    // Set branch addresses for angle information
    tree->Branch("NAngles", &fNAngles);
    tree->Branch("AngleFRANSScore", &fAngleFRANSScore);
    tree->Branch("AngleNHits", &fAngleNHits);
    tree->Branch("AngleMainTrackNHits", &fAngleMainTrackNHits);
    tree->Branch("AngleLongestIsMain", &fAngleLongestIsMain);

}
