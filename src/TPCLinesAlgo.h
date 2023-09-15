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


#include "SObjects/TPCLinesParameters.cpp"
#include "SObjects/TPCSimpleHits.h"
#include "SObjects/TPCSimpleClusters.h"
#include "SObjects/TPCSimpleTriangles.h"
#include "SObjects/TPCSimpleEvents.h"
#include "TPCLinesHough.h"
#include "TPCLinesTrackFinder.h"
#include "TPCLinesVertexFinder.h"
#include "TPCLinesDisplay.cpp"


class TPCLinesAlgo {
    private:
        // Algorithm parameters
        TPCLinesAlgoPsetType fTPCLinesPset;

        // Channel boundaries
        std::map<std::string, std::vector<int>> fChB =  {
            {"U0", {0, 1983}},
            {"V0", {1984, 3967}},
            {"C0", {3968, 5631}},
            {"U1", {5632, 7631}},
            {"V1", {7632, 9599}},
            {"C1", {9600, 11263}}
        };

        // Input variables
        size_t fNTotalHits;
        std::vector<SHit> fHitList;
        SVertex fVertex;
        SVertex fVertexTrue;
        int fMaxX;

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

    public:
        // constructor
        TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, std::string displayPath="");
        
        // Function to set the input variables
        void SetHitList(std::string view,
                        std::vector<int>& vertex,
                        std::vector<int>& vertexTrue,
                        std::vector<int> *_X,
                        std::vector<double> *_Y,
                        std::vector<double> *_Int,
                        std::vector<double> *_Wi,
                        std::vector<double> *_ST,
                        std::vector<double> *_ET,
                        std::string eventLabel="");
        
        // Function to analyze the view
        SEvent AnaView(std::string eventLabel);

        // Function to get the average desnity of hits per wire
        double GetAverageHitDensity();

        // Display
        void Display();                            
};

#endif // TPC_SIMPLE_LINES_H