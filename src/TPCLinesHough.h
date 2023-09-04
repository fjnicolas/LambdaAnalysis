////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesHoughs
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_HOUGH_H
#define TPC_LINES_HOUGH_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>


#include "SObjects/TPCLinesParameters.cpp"
#include "SObjects/TPCSimpleHits.h"
#include "SObjects/TPCSimpleLines.h"
#include "TPCLinesDisplay.cpp"


class TPCLinesHough {
    private:
        HoughAlgorithmPsetType fTPCLinesHoughPset;

        std::vector<SHit> GetHitsInBall(std::vector<SHit> hitList, SVertex vertex, double d);
        LineEquation ComputeLineEquationFromHough(double th, SPoint p);

        // weight function
        double HoughWeightDistance(LineEquation line, std::vector<SHit> hitList);

        // alternative weight functions
        double HoughWeightAverageDistance(LineEquation line, std::vector<SHit> hitList);
        double GetLinearR2(std::vector<double> hypoV, std::vector<double> dataV);
        double HoughWeightLinearR2(LineEquation line, std::vector<SHit> hitList);

        // display app directory
        TPCLinesDisplay fDisplay;
        std::string fDisplayAppPath;


    public:
        // constructor
        TPCLinesHough(HoughAlgorithmPsetType tpcLinesHoughPset);

        // main function
        HoughLine GetBestHoughLine(std::vector<SHit> hitList, SVertex vertex);
        
};

#endif // TPC_SIMPLE_LINES_H