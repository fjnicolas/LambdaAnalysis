
////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesVertexFinder
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPCLINES_VERTEXFINDER_H
#define TPCLINES_VERTEXFINDER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleLines.h"
#include "TPCSimpleTriangles.h"
#include "TPCLinesParameters.h"
#include "TPCLinesDistanceUtils.h"
#include "TPCSimpleEvents.h"
#include "TPCLinesDirectionRecoUtils.h"


class TPCLinesVertexFinder {
    private:
        // configuration parameters
        VertexFinderAlgorithmPsetType fTPCLinesVertexFinderPset;

        double GetAngle360(double x, double y);
        
        bool TrackTriangleJunctionConatined(SLinearCluster track, STriangle tri);
        
        int GetNHitsBetweenJunction(SLinearCluster track, STriangle tri, std::vector<SLinearCluster> trackList, SPoint intP, std::vector<int> track1ListIx, 
        std::vector<int> track2ListIx, double tol);
        
        bool areVectorsEqual(const std::vector<int>& vec1, const std::vector<int>& vec2);
        
        SPoint GetTracksIntersection(SLinearCluster track1, SLinearCluster track2, double dMax, bool useEdgeSlopes = true, bool useFit = false);
        
        int GetHitsContainedInLineEquation(LineEquation trackEq, std::vector<SHit> hitList, float tol = 1.0);
        
        std::vector<SHit> GetMutualHitsContainedInHypo(SLinearCluster track1, SLinearCluster track2,  SPoint intP, int nHits, float tol = 1.0);
        
        SPoint check_arrow_line_intersection(float Ax, float Ay, float Dx, float Dy, float line_slope, float line_intercept);
        
        SPoint GetTrackssEuationOppositePoint(SLinearCluster track, std::vector<SLinearCluster> trackList, SPoint p);

        SLinearCluster GetMainDirection(std::vector<SLinearCluster> & mainDirTrackList, 
                                std::vector<SLinearCluster>& freeTrackList,
                                std::vector<std::vector<SLinearCluster>> parallelTracks,
                                bool decideMainTrack,
                                int verbose);

        std::vector<SLinearCluster> GetCollinearTracks(SLinearCluster mainTrack, std::vector<SLinearCluster> trackList);

        bool LambdaDecayKinematicCheck(STriangle Triangle, SLinearCluster MainDirection, SLinearCluster track1, SLinearCluster track2, std::vector<SLinearCluster> FreeTracksList);


    public:
        // constructor
        TPCLinesVertexFinder(VertexFinderAlgorithmPsetType tpcLinesVertexFinderPset);

        // main function
        void GetAngleVertices(std::vector<SLinearCluster> trackList, std::vector<STriangle>& vertexList, std::vector<SPoint> &originList, SLinearCluster &mainDirection);

        std::vector<SOrigin> GetOrigins(std::vector<SLinearCluster> tracksList, SPoint ballVertex);

        
};

#endif // TPC_SIMPLE_LINES_H
