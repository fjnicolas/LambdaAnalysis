////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesPCA.cpp
//
// \brief Definition of TPCLinesPCA
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesPCA.h"

double mean(const std::vector<SHit>& data, double (SHit::*attribute)() const) {
    double sum = 0.0;
    for (const SHit& s : data) {
        sum += (s.*attribute)();
    }
    return sum / data.size();
}

double pearsonCorrelation(const std::vector<SHit>& data) {
    double meanX = 0.0;
    double meanY = 0.0;
    double sumXY = 0.0;
    double sumX2 = 0.0;
    double sumY2 = 0.0;

    for (const SHit& s : data) {
        meanX+=s.X();
        meanY+=s.Y();
    }
    meanX/=data.size();
    meanY/=data.size();

    for (const SHit& s : data) {
        double x = s.X();
        double y = s.Y();
        sumXY += (x - meanX) * (y - meanY);
        sumX2 += std::pow(x - meanX, 2);
        sumY2 += std::pow(y - meanY, 2);
    }

    return sumXY / (std::sqrt(sumX2) * std::sqrt(sumY2));
}


LineEquation TPCLinesPCA::PerformPCA2D(std::vector<SHit>& data) {

    // first center the data
    double meanX = 0.0, meanY = 0.0;
    for (SHit& h : data) {
        meanX += h.X();
        meanY += h.Y();
    }
    meanX /= data.size();
    meanY /= data.size();

    // get the covariance and fill centered vectors
    double covXX = 0.0, covYY = 0.0, covXY = 0.0;
    for (SHit& h : data) {
        double x = h.X()-meanX;
        double y = h.Y()-meanY;
        covXX += x * x;
        covYY += y * y;
        covXY += x * y;
    }
    covXX /= data.size();
    covYY /= data.size();
    covXY /= data.size();

    // Total variance (sum of covariances)
    double total_variance = covXX + covYY;  

    double eigenvalue1 = (covXX + covYY + std::sqrt((covXX - covYY) * (covXX - covYY) + 4 * covXY * covXY)) / 2.0;
    double eigenvectorX1 = eigenvalue1 - covYY;
    double eigenvectorY1 = covXY;

    double explained_variance = eigenvalue1 / total_variance;

    double correlation = pearsonCorrelation(data);

    /* double eigenvalue2 = (covXX + covYY - std::sqrt((covXX - covYY) * (covXX - covYY) + 4 * covXY * covXY)) / 2.0;
    double eigenvectorX2 = eigenvalue2 - covYY;
    double eigenvectorY2 = covXY;*/


    // Calculate slopes (m) and intercepts (b) for lines
    double slope1 = eigenvectorY1 / eigenvectorX1;
    double yIntercept1 = meanY - slope1 * meanX; // Since the line passes through the mean
    /*std::cout << "Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << std::endl;
    std::cout << "Eigenvectors: (" << eigenvectorX1 << ", " << eigenvectorY1 << "), (" << eigenvectorX2 << ", " << eigenvectorY2 << ")" << std::endl;
    double slope2 = eigenvectorY2 / eigenvectorX2;
    double yIntercept2 = 0.0; // Since the line passes through the origin (0, 0)*/
    /*std::cout << "Slope-Intercept Form of Lines:" << std::endl;
    std::cout << "Line 1: y = " << slope1 << "x + " << yIntercept1 << std::endl;*/

    return LineEquation(slope1, yIntercept1, correlation);
}


LineEquation TPCLinesPCA::PerformPCA2DThreshold(std::vector<SHit> &hitList, float th, bool reverse) {
    if(th > 1) th = 1;

    std::vector<SHit> points = hitList;

    // Order points by x
    if(reverse==true)
        std::sort(points.begin(), points.end(), [](SHit& a, SHit& b) {return a.X() > b.X();});
    else
        std::sort(points.begin(), points.end(), [](SHit& a, SHit& b) {return a.X() < b.X();});
    
    // Get first selected point
    points.resize(static_cast<size_t>(th * points.size()));

    return PerformPCA2D(points);
}
