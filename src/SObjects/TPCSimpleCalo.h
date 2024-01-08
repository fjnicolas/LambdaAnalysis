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
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include <TVectorD.h>
#include <TMatrixD.h>
#include "TVirtualFitter.h"


#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleHits.h"
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleTriangles.h"
#else
#include "TPCSimpleHits.h"
#include "TPCSimpleTriangles.h"
#endif

// --- Drawing options ---
const double fMaxMarkerSize = 5;
const double fAlpha = 0.6;
const double fColor1 = kAzure+7;
const double fColor2 = kOrange+8;

class SCalo {
public:
    // Constructor to initialize the collection
    SCalo(const std::vector<SHit>& points = {}, double angle=0);

    // Method to display the calculated path lengths
    void Display();

    std::vector<double> GetResidualRange() const {return fResidualRange;};
    std::vector<double> GetDepositedEnergy() const {return fDepositedEnergy;};
    double GetTrackLength() const {return fTrackLength;};
    double GetResidualRangeRMS() const {return fResidualRangeRMS;};

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
    // Residual range RMS
    double fResidualRangeRMS;

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

        std::vector<SHit> fAllHits, fHitsTrack1, fHitsTrack2, fVertexHits, fResidualHits;
        // Fit parameters
        bool fPassFit;
        double fM1;
        double fN1;
        double fM2;
        double fXB;
        double fN2;

        // Range
        double fXVertex;
        double fYVertex;
        double fXFitVtx;
        double fYFitVtx;
        double fMinX;
        double fMaxX;
        double fMinY;
        double fMaxY;
        double fMaxIntegral;

        double fChargeRatioAverage;
        double fChargeDifferenceAverage;
        double fChargeRelativeDifferenceAverage;

        bool fPassChargeFit;
        double fBandOverlap;
        double fChargeRatioFit;
        double fChargeDifferenceFit;
        double fChargeRelativeDifferenceFit;

        int fNVertexHits;
        int fNBulkHits;
        double fVertexHitIntegralRatio;
        double fVertexHitIntegralDifference;
        double fVertexHitIntegralRelativeDifference;

        double fTrackLength1;
        double fTrackLength2;
        double fTrackLengthRatio;
        double fResidualRange1RMS;
        double fResidualRange2RMS;

        double fTwoLinesChi2;

        void MakeEnergyLossVsResidualRangePlot(SCalo calo1, SCalo calo2, TPad *pad);

    public:
        STriangleCalo(STriangle triangle);

        // Function to display dEdx
        void CreateEnergyLossVsResidualRangePlot();
        // Triangle joint fit analysis
        void JointFitAnalysis(unsigned int maxHits, double widthTol, bool useHitError);
        void JointFitAnalysisFisher();
        void Display(TCanvas *c1);

        // Return functions
        double ChargeRatioAverage() {return fChargeRatioAverage;};
        double ChargeDifferenceAverage() {return fChargeDifferenceAverage;};
        double ChargeRelativeDifferenceAverage() {return fChargeRelativeDifferenceAverage;};
        bool PassFit() {return fPassFit;};
        double ChargeRatioFit() {return fChargeRatioFit;};
        double ChargeDifferenceFit() {return fChargeDifferenceFit;};
        double ChargeRelativeDifferenceFit() {return fChargeRelativeDifferenceFit;};
        int NVertexHits() {return fNVertexHits;};
        int NBulkHits() {return fNBulkHits;};
        double VertexHitIntegralRatio() {return fVertexHitIntegralRatio;};
        double VertexHitIntegralDifference() {return fVertexHitIntegralDifference;};
        double VertexHitIntegralRelativeDifference() {return fVertexHitIntegralRelativeDifference;};
        bool PassChargeFit() {return fPassChargeFit;};
        double BandOverlap() {return fBandOverlap;};
        double TwoLinesChi2() {return fTwoLinesChi2;};
        double TrackLengthRatio() {return fTrackLengthRatio;};
        double TrackLength1() {return fTrackLength1;};
        double TrackLength2() {return fTrackLength2;};
        double ResidualRange1RMS() {return fResidualRange1RMS;};
        double ResidualRange2RMS() {return fResidualRange2RMS;};

};



#endif // TPC_SIMPLE_CALO_H