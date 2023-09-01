////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleClusters.h
//
// \brief Definition of SimpleTPCClusters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_CLUSTERS_H
#define TPC_SIMPLE_CLUSTERS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"

class SCluster {
    private:
        int fNHits = 0;

        std::vector<SHit> fHitList;
        std::vector<double> fConnectednessV;
        std::vector<double> fCompactnessV;
        std::vector<double> fConnectedness1DV;
        std::vector<double> fWidths;
        
        double fCompactness;
        double fCompactnessRMS;
        double fConnectedness;
        double fConnectednessRMS;
        double fConnectedness1D;
        double fConnectedness1DRMS;
        double fAverageWidth;

        static constexpr double DefaultMaxDistance = std::numeric_limits<double>::max();
    
    public:
        SCluster(std::vector<SHit> hitList={});

        int NHits(){return fNHits;}
        void AddHit(SHit& hit);
        
        // return methods
        std::vector<SHit> GetHits(){return fHitList;};
        std::vector<double> GetConnectednessV(){return fConnectednessV;};
        double GetConnectedness(){return fConnectedness;};
        double GetConnectednessRMS(){return fConnectednessRMS;};
        std::vector<double> GetCompactnessV(){return fCompactnessV;};
        double GetCompactness(){return fCompactness;};
        double GetCompactnessRMS(){return fCompactnessRMS;};
        double GetAverageWidth(){return fAverageWidth;};

        // get hit distance to clusters method
        template <typename T> double GetMinDistanceToCluster(const T& h);
        template <typename T> double GetMinDistanceToClusterW(const T& h);
        template <typename T> double GetMinDistanceToClusterOverlap(const T& h);
        template <typename T> double GetMinDistanceToCluster1D(const T& h);
        template <typename T> std::pair<SHit, double>  GetClosestHitToPoint(const T& h);
        

        friend std::ostream& operator<<(std::ostream& out, SCluster & c);

};



// Function to detect outliers using IQR method
std::vector<int> detect_outliers_iqr(std::vector<double> data, double threshold = 1.5) {
    // Calculate the first and third quartiles (25th and 75th percentiles)
    std::sort(data.begin(), data.end());
    double q1 = data[data.size() / 4];
    double q3 = data[data.size() * 3 / 4];

    // Calculate the IQR
    double iqr = q3 - q1;

    // Calculate the lower and upper bounds for outliers
    double lower_bound = q1 - threshold * iqr;
    double upper_bound = q3 + threshold * iqr;

    // Find the outliers in the data
    std::vector<int> outlier_indices;
    for (int i = 0; i < data.size(); i++) {
        if (data[i] < lower_bound || data[i] > upper_bound) {
            outlier_indices.push_back(i);
        }
    }

    return outlier_indices;
}



// Class to represent a linear cluster
class SLinearCluster {
    private:
        int fId;
        SCluster fHitCluster;

        static constexpr double DefaultMax = std::numeric_limits<double>::max();
   
        float fMinX = DefaultMax;
        float fMinY = DefaultMax;
        float fMaxX = 0;
        float fMaxY = 0;
        float fYAtMaxX = 0;
        float fYAtMinX = 0;
        float fMeanX = 0;
        float fMeanY = 0;

        LineEquation fTrackEquation;
        LineEquation fTrackEquationStart;
        LineEquation fTrackEquationEnd;
        
        SCluster fResidualHitCluster;
        SCluster fMainHitCluster;

        SPoint fStartPoint;
        SPoint fEndPoint;
        SPoint fCoMPoint;
        bool fHasStartEndPoints;
        bool fHasResidualHits;


    public:
        
        // constructor
        SLinearCluster(std::vector<SHit> hitList={});

        // return methods
        int GetId() {return fId;}
        std::vector<SHit> GetHits() {return fHitCluster.GetHits();}
        int NHits() {return fHitCluster.NHits();}
        double GetIntegral();

        // Getter functions
        float GetMinX()  { return fMinX; }
        float GetMinY()  { return fMinY; }
        float GetMaxX()  { return fMaxX; }
        float GetMaxY()  { return fMaxY; }
        float GetMeanX()  { return fMeanX; }
        float GetMeanY()  { return fMeanY; }
        float GetYatMinX()  { return fYAtMinX; }
        float GetYatMaxX()  { return fYAtMaxX; }
        SPoint GetStartPoint()  { return fStartPoint; }
        SPoint GetEndPoint()  { return fEndPoint; }
        SPoint GetCoMPoint()  { return fCoMPoint; }

        LineEquation GetTrackEquation() { return fTrackEquation; }
        LineEquation GetTrackEquationStart() { return fTrackEquationStart; }
        LineEquation GetTrackEquationEnd() { return fTrackEquationEnd; }

        SCluster GetHitCluster(){return fHitCluster;};
    
        std::vector<double> GetConnectednessV(){return fHitCluster.GetConnectednessV();};
        double GetConnectedness(){return fHitCluster.GetConnectedness();};
        std::vector<double> GetCompactnessV(){return fHitCluster.GetCompactnessV();};
        double GetCompactness(){return fHitCluster.GetCompactness();};
        double GetAverageWidth(){return fHitCluster.GetAverageWidth();};
    
        void AssignId(int id) {fId = id;}
        void AddHit(SHit hit) {fHitCluster.AddHit(hit);}

        void FillResidualHits();

};



#endif // TPC_SIMPLE_CLUSTERS_H