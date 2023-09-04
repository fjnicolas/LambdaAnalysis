////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesHough
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include <numeric>

#include "TPCLinesHough.h"

// Constructor
TPCLinesHough::TPCLinesHough(HoughAlgorithmPsetType tpcLinesHoughPset):
    fTPCLinesHoughPset(tpcLinesHoughPset),
    fDisplay(TPCLinesDisplay(tpcLinesHoughPset.Verbose>0))
{}


std::vector<SHit> TPCLinesHough::GetHitsInBall(std::vector<SHit> hitList, SVertex vertex, double d) {
    std::vector<SHit> hitsInBall;
    for (int ix = 0; ix < hitList.size(); ix++) {
        // distance: accounts different units for X and Y directions 
        if (std::sqrt((0.3 * 0.3) * std::pow(hitList[ix].X() - vertex.X(), 2) + (0.08 * 0.08) * std::pow(hitList[ix].Y() - vertex.Y(), 2)) < d) {
            hitsInBall.push_back(hitList[ix]);
        }
    }
    return hitsInBall;
}

LineEquation TPCLinesHough::ComputeLineEquationFromHough(double th, SPoint p) {
    double m = std::tan(th);
    double n = p.Y() - p.X() * m;
    return LineEquation(m, n);
}

double TPCLinesHough::HoughWeightDistance(LineEquation line, std::vector<SHit> hitList) {
    double weight = 0;
    for (auto hit : hitList) {
        double d = line.GetDistance(SPoint{hit.X(), hit.Y()});
        if (d < 1) {
            weight += 1;
        } else {
            weight += 1. / d;
        }
    }
    return weight;
}

double TPCLinesHough::HoughWeightAverageDistance(LineEquation line, std::vector<SHit> hitList) {
    std::vector<double> D;
    for (auto hit : hitList) {
        SPoint p{hit.X(), hit.Y()};
        SPoint pProj = line.GetLineClosestPoint(p);
        double d = std::sqrt(std::pow(p.X() - pProj.X(), 2) + std::pow(p.Y() - pProj.Y(), 2));
        D.push_back(d);
    }
    float mean = std::accumulate(D.begin(), D.end(), 0)/D.size();
    return 1. / mean;
}


double TPCLinesHough::GetLinearR2(std::vector<double> hypoV, std::vector<double> dataV) {
    if (dataV.size() > 1) {
        double SSres = 0.0;
        double SStot = 0.0;
        double meanData = std::accumulate(dataV.begin(), dataV.end(), 0)/dataV.size();
        for (int k = 0; k < dataV.size(); k++) {
            SSres += std::pow(hypoV[k] - dataV[k], 2);
            SStot += std::pow(dataV[k] - meanData, 2);
        }
        return 1. - SSres / SStot;
    } else {
        return 0;
    }
}


double TPCLinesHough::HoughWeightLinearR2(LineEquation line, std::vector<SHit> hitList) {
    std::vector<double> HypoYV;
    std::vector<double> DataYV;

    for (auto hit : hitList) {
        HypoYV.push_back(line.Slope() * hit.X() + line.Intercept());
        DataYV.push_back(hit.Y());
    }

    return GetLinearR2(HypoYV, DataYV);
}

HoughLine TPCLinesHough::GetBestHoughLine(std::vector<SHit> hitList, SVertex vertex)//, double thetaRes, double maxRadius, double maxDistance, bool displayDirectionHypo = false) {
{

    // hits near the vertex, used to make the Hough line hyposthesis
    std::vector<SHit> hitPivotList = GetHitsInBall(hitList, vertex, fTPCLinesHoughPset.MaxRadiusLineHypothesis);
    HoughLine bestHoughLine;

    std::cout<<"  Pivot hits: "<<hitPivotList.size()<<std::endl;

    // theta step (in radians)
    double thetaStep = M_PI*fTPCLinesHoughPset.ThetaRes/180;

    // loop over the pivot hits
    for (auto pivotHit : hitPivotList) {
       
        SPoint P0( pivotHit.X(), pivotHit.Y() );

        // make different hough hypiosthis for each angle
        for (double theta = 0; theta < M_PI; theta +=thetaStep) {

            LineEquation hypoLine = ComputeLineEquationFromHough(theta, P0);
            
            // get hits closer to the Hough hypothesis
            std::vector<SHit> hitHoughList;
            for (auto hit : hitList) {
                if (hit.Id() == pivotHit.Id()) continue;
    
                double d_clo = hypoLine.GetDistance(SPoint{hit.X(), hit.Y()});

                if (d_clo < fTPCLinesHoughPset.MaxDistanceTube) {
                    hitHoughList.push_back(hit);
                }
            }

            // get a score to the hough line
            double hypoLineWeight = HoughWeightDistance(hypoLine, hitHoughList);

            if (hypoLineWeight > bestHoughLine.Score()) {
                bestHoughLine.SetLineEquation ( hypoLine );
                bestHoughLine.SetScore( hypoLineWeight ); 
                bestHoughLine.SetNHits( hitHoughList.size() );

                if(fTPCLinesHoughPset.Verbose>=3) fDisplay.Show("New Hough direction", hitList, hypoLine, hitHoughList);
            }
        }
    }

    return bestHoughLine;
}