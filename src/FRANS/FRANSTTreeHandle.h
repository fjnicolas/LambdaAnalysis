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
                  double alpha, double omega, double tau, double iota,
                  int nOrigins=0, int nOriginsM1=0, int nOriginsM2=0, int nOriginsM3=0, double hitDensity=0)
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
            fNOrigins_U = nOrigins;
            fNOriginsM1_U = nOriginsM1;
            fNOriginsM2_U = nOriginsM2;
            fNOriginsM3_U = nOriginsM3;
            fHitDensity_U = hitDensity;
        }
        else if(view==1){
            fDelta_V = delta;
            fEta_V = eta;
            fFitScore_V = fitScore;
            fAlpha_V = alpha;
            fOmega_V = omega;
            fTau_V = tau;
            fIota_V = iota;
            fNOrigins_V = nOrigins;
            fNOriginsM1_V = nOriginsM1;
            fNOriginsM2_V = nOriginsM2;
            fNOriginsM3_V = nOriginsM3;
            fHitDensity_V = hitDensity;
        }
        else if(view==2){
            fDelta_C = delta;
            fEta_C = eta;
            fFitScore_C = fitScore;
            fAlpha_C = alpha;
            fOmega_C = omega;
            fTau_C = tau;
            fIota_C = iota;
            fNOrigins_C = nOrigins;
            fNOriginsM1_C = nOriginsM1;
            fNOriginsM2_C = nOriginsM2;
            fNOriginsM3_C = nOriginsM3;
            fHitDensity_C = hitDensity;
        }
        else if(view==-1){
            fDelta = delta;
            fEta = eta;
            fFitScore = fitScore;
            fAlpha = alpha;
            fOmega = omega;
            fTau = tau;
            fIota = iota;
            fNOrigins = nOrigins;
            fNOriginsM1 = nOriginsM1;
            fNOriginsM2 = nOriginsM2;
            fNOriginsM3 = nOriginsM3;
            fHitDensity = hitDensity;
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
        fTree->Branch("NOrigins_U", &fNOrigins_U, "fNOrigins_U/I");
        fTree->Branch("NOriginsM1_U", &fNOriginsM1_U, "NOriginsM1_U/I");
        fTree->Branch("NOriginsM2_U", &fNOriginsM2_U, "NOriginsM2_U/I");
        fTree->Branch("NOriginsM3_U", &fNOriginsM3_U, "NOriginsM3_U/I");
        fTree->Branch("HitDensity_U", &fHitDensity_U, "HitDensity_U/D");

        fTree->Branch("Delta_V", &fDelta_V, "Delta_V/D");
        fTree->Branch("Eta_V", &fEta_V, "Eta_V/D");
        fTree->Branch("FitScore_V", &fFitScore_V, "FitScore_V/D");
        fTree->Branch("Alpha_V", &fAlpha_V, "Alpha_V/D");
        fTree->Branch("Omega_V", &fOmega_V, "Omega_V/D");
        fTree->Branch("Tau_V", &fTau_V, "Tau_V/D");
        fTree->Branch("Iota_V", &fIota_V, "Iota_V/D");
        fTree->Branch("NOrigins_V", &fNOrigins_V, "NOrigins_V/I");
        fTree->Branch("NOriginsM1_V", &fNOriginsM1_V, "NOriginsM1_V/I");
        fTree->Branch("NOriginsM2_V", &fNOriginsM2_V, "NOriginsM2_V/I");
        fTree->Branch("NOriginsM3_V", &fNOriginsM3_V, "NOriginsM3_V/I");
        fTree->Branch("HitDensity_V", &fHitDensity_V, "HitDensity_V/D");
       
        fTree->Branch("Delta_C", &fDelta_C, "Delta_C/D");
        fTree->Branch("Eta_C", &fEta_C, "Eta_C/D");
        fTree->Branch("FitScore_C", &fFitScore_C, "FitScore_C/D");
        fTree->Branch("Alpha_C", &fAlpha_C, "Alpha_C/D");
        fTree->Branch("Omega_C", &fOmega_C, "Omega_C/D");
        fTree->Branch("Tau_C", &fTau_C, "Tau_C/D");
        fTree->Branch("Iota_C", &fIota_C, "Iota_C/D");
        fTree->Branch("NOrigins_C", &fNOrigins_C, "NOrigins_C/I");
        fTree->Branch("NOriginsM1_C", &fNOriginsM1_C, "NOriginsM1_C/I");
        fTree->Branch("NOriginsM2_C", &fNOriginsM2_C, "NOriginsM2_C/I");
        fTree->Branch("NOriginsM3_C", &fNOriginsM3_C, "NOriginsM3_C/I");
        fTree->Branch("HitDensity_C", &fHitDensity_C, "HitDensity_C/D");

        fTree->Branch("Delta", &fDelta, "Delta/D");
        fTree->Branch("Eta", &fEta, "Eta/D");
        fTree->Branch("FitScore", &fFitScore, "FitScore/D");
        fTree->Branch("Alpha", &fAlpha, "Alpha/D");
        fTree->Branch("Omega", &fOmega, "Omega/D");
        fTree->Branch("Tau", &fTau, "Tau/D");
        fTree->Branch("Iota", &fIota, "Iota/D");
        fTree->Branch("NOrigins", &fNOrigins, "NOrigins/I");
        fTree->Branch("NOriginsM1", &fNOriginsM1, "NOriginsM1/I");
        fTree->Branch("NOriginsM2", &fNOriginsM2, "NOriginsM2/I");
        fTree->Branch("NOriginsM3", &fNOriginsM3, "NOriginsM3/I");
        fTree->Branch("HitDensity", &fHitDensity, "HitDensity/D");
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

    double fDelta_U, fDelta_V, fDelta_C, fDelta;  // Add fDelta
    double fEta_U, fEta_V, fEta_C, fEta;  // Add fEta
    double fFitScore_U, fFitScore_V, fFitScore_C, fFitScore;  // Add fFitScore
    double fAlpha_U, fAlpha_V, fAlpha_C, fAlpha;  // Add fAlpha
    double fOmega_U, fOmega_V, fOmega_C, fOmega;  // Add fOmega
    double fTau_U, fTau_V, fTau_C, fTau;  // Add fTau
    double fIota_U, fIota_V, fIota_C, fIota;  // Add fIota
    double fScore_U, fScore_V, fScore_C, fScore;  // Add fScore

    double fHitDensity_U, fHitDensity_V, fHitDensity_C, fHitDensity;
    int fNOrigins_U, fNOriginsM1_U, fNOriginsM2_U, fNOriginsM3_U;
    int fNOrigins_V, fNOriginsM1_V, fNOriginsM2_V, fNOriginsM3_V;
    int fNOrigins_C, fNOriginsM1_C, fNOriginsM2_C, fNOriginsM3_C;
    int fNOrigins, fNOriginsM1, fNOriginsM2, fNOriginsM3;


};



#endif
