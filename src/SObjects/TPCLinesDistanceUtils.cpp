////////////////////////////////////////////////////////////////////////////
//
// \file DistanceUtils.h
//
// \brief Definition of fistsance functions
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef DISTANCE_UTILS_H
#define DISTANCE_UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"

namespace TPCLinesDistanceUtils{
    // Mean and Std Dev function for std::vector
    double CalculateMean( std::vector<double>& values) {

        double sum = 0;
        for ( double& value : values) {
            sum += value;
        }
        return sum / values.size();
    }

    double CalculateStdDev( std::vector<double>& values) {

        double mean = CalculateMean(values);
        double sumSqDiff = 0;
        for ( double& value : values) {
            double diff = value - mean;
            sumSqDiff += diff * diff;
        }
        return std::sqrt(sumSqDiff / values.size());
    }


    template <typename T1, typename T2>
    double GetHitDistance(const T1& hit1, const T2& hit2) {
        double xdif = hit1.X() - hit2.X();
        double ydif = hit1.Y() - hit2.Y();

        return std::sqrt(xdif * xdif + ydif * ydif);
    }


    // Hit distance funcitons
    double GetPointDistance( SPoint p1,  SPoint p2) {

        double xdif = p1.X() - p2.X();
        double ydif = p1.Y() - p2.Y();

        return std::sqrt( xdif*xdif + ydif*ydif);
    }

    template <typename T1, typename T2>
    int GetHitDistance1D(const T1& hit1, const T2& hit2) {
        return std::abs(hit1.X() - hit2.X());
    }

    template <typename T1, typename T2>
    double GetHitDistanceW(const T1& hit1, const T2& hit2) {
        
        // distance in X
        double dX = std::pow(hit1.X() - hit2.X(), 2);
        
        // distance between centers
        double d0 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y(), 2));
        
        // dist of hit 1 to width of hit 2
        double dp1 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() + hit2.Width(), 2));
        double dm1 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() - hit2.Width(), 2));
        double d1 = std::min(dp1, dm1);
        
        // dist of hit 2  to width of hit 1
        double dp2 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() + hit1.Width(), 2));
        double dm2 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() - hit1.Width(), 2));
        double d2 = std::min(dp2, dm2);
        
        // return min
        return std::min(d0, (d1 + d2) / 2.);
    }

    template <typename T1, typename T2>
    double GetHitDistanceOverlap(const T1& hit1, const T2& hit2) {
        double y_range1_min = hit1.Y() - hit1.Width();
        double y_range1_max = hit1.Y() + hit1.Width();
        double y_range2_min = hit2.Y() - hit2.Width();
        double y_range2_max = hit2.Y() + hit2.Width();
        bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
        double dX = std::pow(hit1.X() - hit2.X(), 2);
        double dY = 1;
        if (!overlap) {
            dY = std::pow(hit1.Y() - hit2.Y(), 2);
        }
        return std::sqrt(dX + dY);
    }


    bool HitWidthOverlap(SHit& hit1, SHit& hit2) {
        std::pair<double, double> y_range1 = {hit1.Y() - hit1.Width(), hit1.Y() + hit1.Width()};
        std::pair<double, double> y_range2 = {hit2.Y() - hit2.Width(), hit2.Y() + hit2.Width()};
        
        return y_range1.first <= y_range2.second && y_range1.second >= y_range2.first;
    }



    double GetClusterMinDistance(SCluster cluster1, SCluster cluster2) {
        double minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = std::sqrt(std::pow(hit1.X() - hit2.X(), 2) + std::pow(hit1.Y() - hit2.Y(), 2));
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }

    double GetClusterConnectedness(SCluster cluster1, SCluster cluster2) {
        double minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = GetHitDistanceW(hit1, hit2);
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }

    double GetClusterConnectednessOverlap(SCluster cluster1, SCluster cluster2) {
        double minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = GetHitDistanceOverlap(hit1, hit2);
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }


    double GetpClusterDistanceW(SHit p1, SHit p2) {
        double d0 = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
        double dp = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y() + p2.Width(), 2));
        double dm = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y() - p2.Width(), 2));
        return std::min(d0, std::min(dp, dm));
    }


    double GetpClusterDistanceOverlap(SHit p1, SHit p2) {
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

}

#endif // TPC_SIMPLE_CLUSTERS_H