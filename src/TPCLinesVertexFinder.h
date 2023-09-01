
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

#include "SObjects/TPCSimpleHits.h"
#include "SObjects/TPCSimpleLines.h"
#include "SObjects/TPCSimpleTriangles.h"
#include "SObjects/TPCLinesDistanceUtils.cpp"
#include "TPCLinesDirectionRecoUtils.cpp"



class TPCLinesVertexFinder {
    private:



    public:
        // constructor
        TPCLinesVertexFinder(){};

        // main function
        void GetOrigins(std::vector<SLinearCluster> trackList, std::vector<STriangle>& vertexList, std::vector<SPoint> &originList, SLinearCluster &mainDirection);

        
};

#endif // TPC_SIMPLE_LINES_H
