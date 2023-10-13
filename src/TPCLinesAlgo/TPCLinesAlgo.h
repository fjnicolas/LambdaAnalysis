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


/*#include "SObjects/TPCLinesParameters.cpp"
#include "SObjects/TPCSimpleHits.h"
#include "SObjects/TPCSimpleClusters.h"
#include "SObjects/TPCSimpleTriangles.h"
#include "SObjects/TPCSimpleEvents.h"*/
#include "TPCLinesParameters.h"
#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleTriangles.h"
#include "TPCSimpleEvents.h"

#include "TPCLinesHough.h"
#include "TPCLinesTrackFinder.h"
#include "TPCLinesVertexFinder.h"
#include "TPCLinesDisplay.h"
#include "TPCLinesDirectionRecoUtils.h"


class TPCLinesAlgo {
    private:
        // Algorithm parameters
        TPCLinesAlgoPsetType fTPCLinesPset;

        // Input variables
        size_t fNTotalHits;
        std::vector<SHit> fHitList;
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

        std::vector<SLinearCluster> MergeIsolatedHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, double dCleaning1D, double dTh = 3);

        SEvent fRecoEvent;
        SPoint fMainVertex;

    public:
        // constructor
        TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset);
        
        // Function to set the input variables
        bool SetHitList(int view,
                            std::vector<int>& vertex,
                            std::vector<int>& vertexTrue,
                            std::vector<SHit> hits);
        
        // Function to analyze the view
        void AnaView(std::string eventLabel);

        // Function to get the average desnity of hits per wire
        double GetAverageHitDensity();

        // Function to get the best view
        int GetBestView(std::vector<int> *_Ch, std::vector<double> *_Chi2);

        int ShiftX(){ return fMinX; };
        double ShiftY(){ return fMinY; };

        SEvent GetRecoEvent() { return fRecoEvent; };
        SPoint GetMainVertex() { return fMainVertex; };
        // Display
        void Display(std::string labelEv, std::string label);
};

#endif // TPC_SIMPLE_LINES_H