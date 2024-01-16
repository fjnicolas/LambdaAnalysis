////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_ALGO_H
#define TPC_LINES_ALGO_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <unordered_set>
#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleHits.h"
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleClusters.h"
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleTriangles.h"
#include "sbndcode/LambdaAnalysis/src/SObjects/TPCSimpleEvents.h"

#include "sbndcode/LambdaAnalysis/src/TPCLinesAlgo/TPCLinesParameters.h"
#include "sbndcode/LambdaAnalysis/src/TPCLinesAlgo/TPCLinesHough.h"
#include "sbndcode/LambdaAnalysis/src/TPCLinesAlgo/TPCLinesTrackFinder.h"
#include "sbndcode/LambdaAnalysis/src/TPCLinesAlgo/TPCLinesVertexFinder.h"
#include "sbndcode/LambdaAnalysis/src/TPCLinesAlgo/TPCLinesDisplay.h"
#include "sbndcode/LambdaAnalysis/src/TPCLinesAlgo/TPCLinesDirectionRecoUtils.h"

#else
#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleTriangles.h"
#include "TPCSimpleEvents.h"

#include "TPCLinesParameters.h"
#include "TPCLinesHough.h"
#include "TPCLinesTrackFinder.h"
#include "TPCLinesVertexFinder.h"
#include "TPCLinesDisplay.h"
#include "TPCLinesDirectionRecoUtils.h"
#endif





class TPCLinesAlgo {
    private:
        // Algorithm parameters
        TPCLinesAlgoPsetType fTPCLinesPset;

        // Input variables
        size_t fNTotalHits;
        std::vector<SHit> fHitList;
        std::vector<SHit> fHitListOutROI;
        std::vector<SHit> fUnmatchedHits;
        std::map<int, int> fClusterIdCounter;
        SVertex fVertex;
        SVertex fVertexTrue;
        int fMinX;
        double fMinY;
        int fMaxX;

        std::vector<SLinearCluster> fFinalTrackList;
        SLinearCluster fMainDirection;
        std::vector<STriangle> fAngleList;
        std::vector<SOrigin> fOrigins;

        // Hough algorithm
        TPCLinesHough fHoughAlgo;

        // Track finder and clustering algorithm
        TPCLinesTrackFinder fTrackFinder;

        // Vertex finder algorithm
        TPCLinesVertexFinder fVertexFinder;

        // Display app directory
        TPCLinesDisplay fDisplay;
        std::string fDisplayAppPath;

        std::vector<SHit> RemoveIsolatedHits(std::vector<SHit> hitListForHough, std::vector<SHit>& discardedHits, double maxD, int minNeighbours);

        std::vector<SLinearCluster> MergeIsolatedHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, double dCleaning1D, std::vector<SHit> & discardedHits, double dTh = 3);

        std::vector<SLinearCluster> MergeOutOfROIHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList);

        std::vector<SLinearCluster> CreateOutOfROIClusters(std::vector<SLinearCluster> recoTrackList);

        std::vector<SHit> GetHitsInCluster(int clusterId);

        SEvent fRecoEvent;
        SPoint fMainVertex;

    public:
        // --- Constructor ---
        TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset);
        
        // --- Function to set the input variables ---
        bool SetHitList(int view,
                        std::vector<int>& vertex,
                        std::vector<int>& vertexTrue,
                        std::vector<SHit> hits);
        
        // --- Main analysis function ---
        void AnaView(std::string eventLabel);

        // --- Getters ---
        // Get average desnity of hits per wire
        double GetAverageHitDensity();
        
        // Get best view (best average chi2)
        int GetBestView(std::vector<int> *_Ch, std::vector<double> *_Chi2);
        
        // Get X-Y shifts
        int ShiftX(){ return fMinX; };
        double ShiftY(){ return fMinY; };
        
        // Get N input hits
        int GetNInputHits() {return fNTotalHits; };

        // Get reco event
        SEvent GetRecoEvent() { return fRecoEvent; };
        SPoint GetMainVertex() { return fMainVertex; };

        // --- Display function ---
        void Display(std::string name, TCanvas *canvas = nullptr);

        
};

#endif // TPC_SIMPLE_LINES_H