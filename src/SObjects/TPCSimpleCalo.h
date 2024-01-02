////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////


#ifndef TPC_SIMPLE_CALO_H
#define TPC_SIMPLE_CALO_H


#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include <TVectorD.h>
#include <TMatrixD.h>


#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleHits.h"
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleTriangles.h"
#include "sbndcode/LambdaAnalysis/src/FRANS/ChargeDensity.h"
#include "sbndcode/LambdaAnalysis/src/FRANS/ChargeDensityPset.h"
#else
#include "TPCSimpleHits.h"
#include "TPCSimpleTriangles.h"
#include "ChargeDensity.h"
#include "ChargeDensityPset.h"
#endif

// --- Drawing options ---
const double fMaxMarkerSize = 5;
const double fAlpha = 0.6;
const double fColor1 = kAzure+7;
const double fColor2 = kOrange+8;

class SCalo {
public:
    // Constructor to initialize the collection
    SCalo(const std::vector<SHit>& points, double angle=0);

    // Method to display the calculated path lengths
    void Display();

    std::vector<double> GetResidualRange() const {return fResidualRange;};
    std::vector<double> GetDepositedEnergy() const {return fDepositedEnergy;};
    double GetTrackLength() const {return fTrackLength;};

private:

    // Conversion factors
    double fHitIntegralToEnergy;
    double fStepXToLength;
    double fStepYToLength;
    
    // Track legnth
    double fTrackLength;
    // TrackOrientationAngle
    double fCosTrackAngle;
    double fTrackAngle;
    // Define the hit points and path lengths as private members
    std::vector<SHit> fHitList;
    std::vector<double> fPathLengths;
    std::vector<double> fDepositedEnergy;
    std::vector<double> fResidualRange;

    // Function to calculate the distance between two hit points
    double CalculateDistance(const SHit& p1, const SHit& p2);

    // Method to calculate and store the path lengths
    void CalculatePathLengths();

    // Function to calculate residual range and deposited energy
    void CalculateResidualRange();
};



class STriangleCalo {
    
    private:
        STriangle fTriangle;
        SCalo fCalo1;
        SCalo fCalo2;

        double fChargeRatioAverage;
        double fChargeDifferenceAverage;
        double fChargeRelativeDifferenceAverage;

        double fChargeRatioFit;
        double fChargeDifferenceFit;
        double fChargeRelativeDifferenceFit;

        double fVertexHitIntegralRatio;
        double fVertexHitIntegralDifference;
        double fVertexHitIntegralRelativeDifference;

        bool fPassFit;
        double fTwoLinesChi2;

        void MakeEnergyLossVsResidualRangePlot(SCalo calo1, SCalo calo2, TPad *pad);

    public:
        STriangleCalo(STriangle triangle);

        // Function to display dEdx
        void CreateEnergyLossVsResidualRangePlot();
        // Triangle joint fit analysis
        void JointFitAnalysis(unsigned int maxHits, double widthTol, bool useHitError, double& fitSlope1, double &fitSlope2, ChargeDensity & chargeDensityAlgo);
        void JointFitAnalysisFisher();

        // Return functions
        double ChargeRatioAverage() {return fChargeRatioAverage;};
        double ChargeDifferenceAverage() {return fChargeDifferenceAverage;};
        double ChargeRelativeDifferenceAverage() {return fChargeRelativeDifferenceAverage;};
        double ChargeRatioFit() {return fChargeRatioFit;};
        double ChargeDifferenceFit() {return fChargeDifferenceFit;};
        double ChargeRelativeDifferenceFit() {return fChargeRelativeDifferenceFit;};
        double VertexHitIntegralRatio() {return fVertexHitIntegralRatio;};
        double VertexHitIntegralDifference() {return fVertexHitIntegralDifference;};
        double VertexHitIntegralRelativeDifference() {return fVertexHitIntegralRelativeDifference;};
        bool PassFit() {return fPassFit;};
        double TwoLinesChi2() {return fTwoLinesChi2;};

};



#endif // TPC_SIMPLE_CALO_H