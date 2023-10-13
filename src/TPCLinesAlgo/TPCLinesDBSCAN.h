////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesDBSCAN.cpp
//
// \brief Definition of TPCLinesDBSCAN
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_DBSCAN_H
#define TPC_LINES_DBSCAN_H

#include <iostream>
#include <vector>
#include <cmath>

#include "TPCSimpleHits.h"

class DBSCAN {
public:
    DBSCAN(double epsilon, int minPts);

    void setDistanceFunction(double (*distFunc)( SHit&,  SHit&));

    void fit( std::vector<SHit>& points);

    std::vector<int>& getClusterAssignment();

private:
    double epsilon;
    int minPts;
    double (*distanceFunction)( SHit&,  SHit&);
    std::vector<int> clusterAssignment;

    std::vector<int> regionQuery( std::vector<SHit>& points,  SHit& p);
    void expandCluster( std::vector<SHit>& points, int pointIdx, int clusterIdx);
};

// --------------- Distance functions for DBSCAN
double DBSCANHitEuclidianDistance( SHit& p1,  SHit& p2);
double DBSCANHitManhattanDistance( SHit& p1,  SHit& p2);
double DBSCANHitWidthDistance( SHit& p1, SHit& p2);
double DBSCANHitOverlapDistance( SHit& p1, SHit& p2);


#endif