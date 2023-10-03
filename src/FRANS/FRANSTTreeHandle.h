///////////////////////////////////////////////////////////////////////
/// File: ChargeDensity.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////
#ifndef SBND_FRANSTTREEHANDLE_H
#define SBND_FRANSTTREEHANDLE_H

#include <TTree.h>

class FRANSTTree {
public:
    FRANSTTree(TTree* tree)
        : fTree(tree)
    {
        InitializeTree();
    }

    void FillData(int view, int runID, int subRunID, int eventID, int isSignal,
                  double delta, double eta, double fitScore,
                  double alpha, double omega, double tau, double iota)
    {
        fRunID = runID;
        fSubRunID = subRunID;
        fEventID = eventID;

        fNNuInts = 0;
        fIntType = -1;
        fIntMode = -1;
        fGap = -1.;
        fIsSignal = isSignal;
        fProtonKE=0;
        fPionKE=0;
        fLambdaKE=0;
        
        if(view==0){
            fDelta_U = delta;
            fEta_U = eta;
            fFitScore_U = fitScore;
            fAlpha_U = alpha;
            fOmega_U = omega;
            fTau_U = tau;
            fIota_U = iota;
        }
        else if(view==1){
            fDelta_V = delta;
            fEta_V = eta;
            fFitScore_V = fitScore;
            fAlpha_V = alpha;
            fOmega_V = omega;
            fTau_V = tau;
            fIota_V = iota;
        }
        else if(view==2){
            fDelta_C = delta;
            fEta_C = eta;
            fFitScore_C = fitScore;
            fAlpha_C = alpha;
            fOmega_C = omega;
            fTau_C = tau;
            fIota_C = iota;
        }
        
    }

    void FillTree(){ fTree->Fill(); }

private:
    void InitializeTree()
    {
        fTree->Branch("RunID", &fRunID, "RunID/I");
        fTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
        fTree->Branch("EventID", &fEventID, "EventID/I");

        fTree->Branch("NNuInts", &fNNuInts, "NNuInts/I");
        fTree->Branch("IntMode", &fIntMode, "IntMode/I");
        fTree->Branch("IntType", &fIntType, "IntType/I");
        fTree->Branch("IsSignal", &fIsSignal, "IsSignal/I");
        fTree->Branch("Gap", &fGap, "Gap/D");
        fTree->Branch("LambdaKE", &fLambdaKE, "LambdaKE/D");
        fTree->Branch("ProtonKE", &fProtonKE, "ProtonKE/D");
        fTree->Branch("PionKE", &fPionKE, "PionKE/D");

        // TTree for BDT training
        fTree->Branch("Delta_U", &fDelta_U, "Delta_U/D");
        fTree->Branch("Eta_U", &fEta_U, "Eta_U/D");
        fTree->Branch("FitScore_U", &fFitScore_U, "FitScore_U/D");
        fTree->Branch("Alpha_U", &fAlpha_U, "Alpha_U/D");
        fTree->Branch("Omega_U", &fOmega_U, "Omega_U/D");
        fTree->Branch("Tau_U", &fTau_U, "Tau_U/D");
        fTree->Branch("Iota_U", &fIota_U, "Iota_U/D");
        
        fTree->Branch("Delta_V", &fDelta_V, "Delta_V/D");
        fTree->Branch("Eta_V", &fEta_V, "Eta_V/D");
        fTree->Branch("FitScore_V", &fFitScore_V, "FitScore_V/D");
        fTree->Branch("Alpha_V", &fAlpha_V, "Alpha_V/D");
        fTree->Branch("Omega_V", &fOmega_V, "Omega_V/D");
        fTree->Branch("Tau_V", &fTau_V, "Tau_V/D");
        fTree->Branch("Iota_V", &fIota_V, "Iota_V/D");
       
        fTree->Branch("Delta_C", &fDelta_C, "Delta_C/D");
        fTree->Branch("Eta_C", &fEta_C, "Eta_C/D");
        fTree->Branch("FitScore_C", &fFitScore_C, "FitScore_C/D");
        fTree->Branch("Alpha_C", &fAlpha_C, "Alpha_C/D");
        fTree->Branch("Omega_C", &fOmega_C, "Omega_C/D");
        fTree->Branch("Tau_C", &fTau_C, "Tau_C/D");
        fTree->Branch("Iota_C", &fIota_C, "Iota_C/D");
    }

private:
    TTree* fTree;
    std::vector<int> fViews;

    int fRunID;
    int fSubRunID;
    int fEventID;


    int fNNuInts;
    int fIntMode;
    int fIntType;
    double fGap;
    double fLambdaKE;
    double fProtonKE;
    double fPionKE;
    int fIsSignal;

    double fDelta_U, fDelta_V, fDelta_C;
    double fEta_U, fEta_V, fEta_C;
    double fFitScore_U, fFitScore_V, fFitScore_C;
    double fAlpha_U, fAlpha_V, fAlpha_C;
    double fOmega_U, fOmega_V, fOmega_C;
    double fTau_U, fTau_V, fTau_C;
    double fIota_U, fIota_V, fIota_C;
    double fScore_U, fScore_V, fScore_C;
};



#endif
