////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleClusters.cpp
//
// \brief Definition of SimpleTPCClusters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleClusters.h"

SCluster::SCluster(std::vector<SHit> hitList){

    fNHits = hitList.size();

    // sort hits by X coordinate
    std::sort(hitList.begin(), hitList.end(), [](SHit& a, SHit& b) {
        return a.X() < b.X();
    });
    if (fNHits > 0) {

        for(auto &hit: hitList){
            double dminComp = DefaultMaxDistance;
            double dminConn = DefaultMaxDistance;
            double dminX = DefaultMaxDistance;
            
            for(auto &hit2: hitList){
                if (&hit == &hit2) continue;

                // distance between centers               
                double d = TPCLinesDistanceUtils::GetHitDistance(hit, hit2);
                if(d<dminComp) dminComp=d;

                // calculate connectedeness
                d = TPCLinesDistanceUtils::GetHitDistanceW(hit, hit2);
                //d = GetHitDistanceOverlap(hit, hit2)
                if(d<dminConn) dminConn=d;

                // calculate connectedness 1D
                d = TPCLinesDistanceUtils::GetHitDistance1D(hit, hit2);
                if(d<dminX) dminX=d;
            }

            SHit newHit = hit;
            newHit.SetHitConnectivity(dminComp, dminConn, dminX);

            fCompactnessV.push_back(dminComp);
            fConnectednessV.push_back(dminConn);
            fConnectedness1DV.push_back(dminX);
            fWidths.push_back(newHit.Width());
            fHitList.push_back(newHit);
        }
            
        fCompactness = TPCLinesDistanceUtils::CalculateMean(fCompactnessV);
        fCompactnessRMS = TPCLinesDistanceUtils::CalculateStdDev(fCompactnessV);
        fConnectedness = TPCLinesDistanceUtils::CalculateMean(fConnectednessV);
        fConnectednessRMS = TPCLinesDistanceUtils::CalculateStdDev(fConnectednessV);
        fConnectedness1D = TPCLinesDistanceUtils::CalculateMean(fConnectedness1DV);
        fConnectedness1DRMS = TPCLinesDistanceUtils::CalculateStdDev(fConnectedness1DV);
        fAverageWidth = TPCLinesDistanceUtils::CalculateMean(fWidths);
    }
    else {
        fCompactness = 0;
        fCompactnessRMS = 0;
        fConnectedness = 0;
        fConnectednessRMS = 0;
        fConnectedness1D = 0;
        fConnectedness1DRMS = 0;
        fAverageWidth = 0;
    }
}

void SCluster::AddHit(SHit& hit) {
    fHitList.push_back(hit);
    fNHits++;
}

std::ostream& operator<<(std::ostream& out, SCluster & c)
{
    out <<  " --- Created Simple Cluster:\n"
    << "  Comp=" << c.fCompactness << " pm " << c.fCompactnessRMS
    << "  Conn=" << c.fConnectedness << " pm " << c.fConnectednessRMS
    << "  Conn1D=" <<c.fConnectedness1D << " pm " << c.fConnectedness1DRMS;

    return out;
}

template <typename T>
double SCluster::GetMinDistanceToCluster(const T& h) {
    double minD = 1e6;
    for (auto& hit : fHitList) {
        double d =TPCLinesDistanceUtils::GetHitDistance(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
double SCluster::GetMinDistanceToClusterW(const T& h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = TPCLinesDistanceUtils::GetHitDistanceW(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
double SCluster::GetMinDistanceToClusterOverlap(const T& h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = TPCLinesDistanceUtils::GetHitDistanceOverlap(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
double SCluster::GetMinDistanceToCluster1D(const T& h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = TPCLinesDistanceUtils::GetHitDistance1D(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
std::pair<SHit, double> SCluster::GetClosestHitToPoint(const T& h) {
    double minD = 1e4;
    size_t closestHitIx = 0;

    for (size_t ix =0; ix<fHitList.size(); ix++){
        SHit hit = fHitList[ix];
        double d = TPCLinesDistanceUtils::GetHitDistance(hit, h);
        if (d < minD){
            minD = d;
            closestHitIx = ix;
        }
    }

    std::pair<SHit, double> result(fHitList[closestHitIx], minD);
    return result;
}



//------------------------------  SLinearCluster

SLinearCluster::SLinearCluster(std::vector<SHit> hitList){

    fId = -1;
    fHitCluster = SCluster(hitList);

    fMinX = DefaultMax;
    fMinY = DefaultMax;
    fMaxX = 0;
    fMaxY = 0;
    fYAtMaxX = 0;
    fYAtMinX = 0;
    fMeanX = 0;
    fMeanY = 0;

    for (
        auto &hit : hitList) {
        if (hit.X() < fMinX) {
            fMinX = hit.X();
            fYAtMinX = hit.Y();
        }
        if (hit.X() > fMaxX) {
            fMaxX = hit.X();
            fYAtMaxX = hit.Y();
        }
        if (hit.Y() < fMinY) fMinY = hit.Y();
        if (hit.Y() > fMaxY) fMaxY = hit.Y();
        fMeanX += hit.X();
        fMeanY += hit.Y();
    }

    if (!hitList.empty()) {
        fMeanX /= hitList.size();
        fMeanY /= hitList.size();
        fCoMPoint = SPoint(fMeanX, fMeanY);
    }

    fStartPoint = SPoint(fMinX, fYAtMinX);
    fEndPoint = SPoint(fMaxX, fYAtMaxX);

    // fill track equations
    if(hitList.size()>2){
        TPCLinesPCA pcaAlgo;
        fTrackEquation = pcaAlgo.PerformPCA2D(hitList);
        if(hitList.size()>6){
            fTrackEquationStart = pcaAlgo.PerformPCA2DThreshold(hitList, 0.5 );
            fTrackEquationEnd = pcaAlgo.PerformPCA2DThreshold(hitList, 0.5, true);
        }
        else{
            fTrackEquationStart = fTrackEquation;
            fTrackEquationEnd = fTrackEquation;
        }
    }
    
    // members ro fill later
    fHasResidualHits = false;
    fHasStartEndPoints = false;
    
}


double SLinearCluster::GetIntegral(){
    double w = 0;
    for (SHit hit : GetHits()) {
        w += hit.Integral();
    }
    return w;
}


std::vector<int> SLinearCluster::detect_outliers_iqr(std::vector<float> data, float threshold) {
    std::vector<int> outlier_indices;

    // Calculate the first and third quartiles (25th and 75th percentiles)
    size_t n = data.size();
    std::vector<float> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());

    size_t q1_index = static_cast<size_t>(n * 0.25);
    size_t q3_index = static_cast<size_t>(n * 0.75);
    float q1 = sorted_data[q1_index];
    float q3 = sorted_data[q3_index];

    // Calculate the IQR
    float iqr = q3 - q1;

    // Calculate the lower and upper bounds for outliers
    float lower_bound = q1 - threshold * iqr;
    float upper_bound = q3 + threshold * iqr;

    // Find the outliers in the data
    for (size_t i = 0; i < n; ++i) {
        if (data[i] < lower_bound || data[i] > upper_bound) {
            outlier_indices.push_back(i);
        }
    }

    return outlier_indices;
}


// Function to detect outliers using IQR method
std::vector<int> SLinearCluster::detect_outliers_iqr2(std::vector<double> data, double threshold) {
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
    for (size_t i = 0; i < data.size(); i++) {
        if (data[i] < lower_bound || data[i] > upper_bound) {
            outlier_indices.push_back(i);
        }
    }

    return outlier_indices;
}

void SLinearCluster::FillResidualHits() {
    fHasResidualHits = true;
    //std::cout << "\n\n +-+-+-+-+-+-+- Filling residual hits +-+-+-+-+-+-+-" << std::endl;

    if (NHits() > 0) {
        // Get average residual
        std::vector<float> dV;
        for (SHit hit : GetHits()) {
            float d = fTrackEquation.GetDistance(SPoint(hit.X(), hit.Y()));
            dV.push_back(d);
        }

        std::vector<int> outliersIx = detect_outliers_iqr(dV, 1.5);
        // std::vector<int> outliersIx = detect_outliers_zscore(dV, 3);

        std::vector<SHit> resHitList;
        std::vector<SHit> mainHitList;
        for (int ix = 0; ix < NHits(); ++ix) {
            if (std::find(outliersIx.begin(), outliersIx.end(), ix) != outliersIx.end()) {
                resHitList.push_back(fHitCluster.GetHits()[ix]);
            } else {
                mainHitList.push_back(fHitCluster.GetHits()[ix]);
            }
        }

        fResidualHitCluster = SCluster(resHitList);
        fMainHitCluster = SCluster(mainHitList);

        TPCLinesPCA pcaAlgo;

        fTrackEquation= pcaAlgo.PerformPCA2D(mainHitList);
        if (mainHitList.size() > 6) {
            fTrackEquationStart = pcaAlgo.PerformPCA2DThreshold(mainHitList, 0.5 );
            fTrackEquationEnd = pcaAlgo.PerformPCA2DThreshold(mainHitList, 0.5, true);
        } else {
            fTrackEquationStart = fTrackEquation;
            fTrackEquationEnd = fTrackEquation;
        }
    }
}

// instances
template double SCluster::GetMinDistanceToCluster(const SHit& h);
template double SCluster::GetMinDistanceToCluster(const SPoint& h);
template double SCluster::GetMinDistanceToClusterW(const SHit& h);
template double SCluster::GetMinDistanceToClusterOverlap(const SHit& h);
template double SCluster::GetMinDistanceToCluster1D(const SHit& h);
template double SCluster::GetMinDistanceToCluster1D(const SPoint& h);
template std::pair<SHit, double>  SCluster::GetClosestHitToPoint(const SHit& h);
template std::pair<SHit, double>  SCluster::GetClosestHitToPoint(const SPoint& h);