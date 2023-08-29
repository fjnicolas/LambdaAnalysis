////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesDBSCAN.cpp
//
// \brief Definition of TPCLinesDBSCAN
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <cmath>

#include "TPCSimpleHits.h"

class DBSCAN {
public:
    DBSCAN(double epsilon, int minPts) : epsilon(epsilon), minPts(minPts) {}

    void setDistanceFunction(double (*distFunc)( SHit&,  SHit&)) {
        distanceFunction = distFunc;
    }

    void fit( std::vector<SHit>& points);

     std::vector<int>& getClusterAssignment()  {
        return clusterAssignment;
    }

private:
    double epsilon;
    int minPts;
    double (*distanceFunction)( SHit&,  SHit&);
    std::vector<int> clusterAssignment;

    std::vector<int> regionQuery( std::vector<SHit>& points,  SHit& p);
    void expandCluster( std::vector<SHit>& points, int pointIdx, int clusterIdx);
};

double euclideanDistance( SHit& p1,  SHit& p2) {
    return std::sqrt((p1.X() - p2.X()) * (p1.X() - p2.X()) + (p1.Y() - p2.Y()) * (p1.Y() - p2.Y()));
}

double manhattanDistance( SHit& p1,  SHit& p2) {
    return std::abs(p1.X() - p2.X()) + std::abs(p1.Y() - p2.Y());
}

double GetpClusterDistanceW( SHit& p1, SHit& p2) {
    double d0 = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
    double dp = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y() + p2.Width(), 2));
    double dm = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y() - p2.Width(), 2));
    return std::min(d0, std::min(dp, dm));
}

double GetpClusterDistanceOverlap( SHit& p1, SHit& p2) {
    double y_range1_min = p1.Y() - p1.Width();
    double y_range1_max = p1.Y() + p1.Width();
    double y_range2_min = p2.Y() - p2.Width();
    double y_range2_max = p2.Y() + p2.Width();
    bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
    double dX = p1.X() - p2.X();
    double dY = 1;
    if (!overlap) {
        dY = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
    }
    return std::sqrt(std::pow(dX, 2) + std::pow(dY, 2));
}

void DBSCAN::fit( std::vector<SHit>& points) {
    clusterAssignment.assign(points.size(), 0);  // 0 indicates unassigned
    
    int clusterIdx = 0;
    for (int i = 0; i < points.size(); ++i) {
        if (clusterAssignment[i] == 0) {  // Not assigned to any cluster
            std::vector<int> neighbors = regionQuery(points, points[i]);
            if (neighbors.size() < minPts) {
                clusterAssignment[i] = -1;  // Mark as noise
            } else {
                ++clusterIdx;
                expandCluster(points, i, clusterIdx);
            }
        }
    }
}

std::vector<int> DBSCAN::regionQuery( std::vector<SHit>& points,  SHit& p) {
    std::vector<int> neighbors;
    for (int i = 0; i < points.size(); ++i) {
        if (distanceFunction(points[i], p) <= epsilon) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void DBSCAN::expandCluster( std::vector<SHit>& points, int pointIdx, int clusterIdx) {
    std::vector<int> seeds = regionQuery(points, points[pointIdx]);
    if (seeds.size() < minPts) {
        clusterAssignment[pointIdx] = -1;  // Mark as noise
        return;
    }
    
    clusterAssignment[pointIdx] = clusterIdx;
    
    for (int i = 0; i < seeds.size(); ++i) {
        int currentIdx = seeds[i];
        if (clusterAssignment[currentIdx] == 0) {  // Not assigned to any cluster
            clusterAssignment[currentIdx] = clusterIdx;
            std::vector<int> currentSHitNeighbors = regionQuery(points, points[currentIdx]);
            if (currentSHitNeighbors.size() >= minPts) {
                seeds.insert(seeds.end(), currentSHitNeighbors.begin(), currentSHitNeighbors.end());
            }
        }
        if (clusterAssignment[currentIdx] == -1) {  // Noise point
            clusterAssignment[currentIdx] = clusterIdx;
        }
    }
}
