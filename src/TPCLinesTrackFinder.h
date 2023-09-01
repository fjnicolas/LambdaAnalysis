////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesTrackFinder.h
//
// \brief Definition of TPCLinesTrackFinder
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_TRACKFINDER_H
#define TPC_LINES_TRACKFINDER_H


#include "SObjects/TPCLinesParameters.cpp"
#include "SObjects/TPCSimpleHits.h"
#include "SObjects/TPCSimpleLines.h"
#include "SObjects/TPCSimpleClusters.h"
#include "TPCLinesDisplay.cpp"


class TPCLinesTrackFinder {
    private:
        TrackFinderAlgorithmPsetType fTPCLinesTrackFinderPset;
        
        double GetHitLinearDensity(std::vector<SHit> data);
        
        LineEquation PerformPCA2D(std::vector<SHit>& data);
        
        void GetHitsInTube(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD);
        
        void GetHitsInTubeSingleWire(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD);
        
        std::vector<std::pair<int, int>> GetListFrequency(std::vector<int> vList);
        
        std::pair<double, double> ComputeConnectivityMode(std::vector<double> v, int minClusterHits, double weightWidth);

        std::vector<SLinearCluster> Get2DClusters(std::vector<SHit> & hits, double epsilon, int minPts, std::string option="");

        std::vector<SLinearCluster> CaptureMissingHits(std::vector<SLinearCluster> trackV, std::vector<SHit> hitList);

        std::vector<SLinearCluster> MakeTrack(std::vector<SLinearCluster> linearClusterList);

        // Display app directory
        TPCLinesDisplay fDisplay;


    public:
        // constructor
        TPCLinesTrackFinder(TrackFinderAlgorithmPsetType tpcLinesTrackFinderPset);

        // main function
        std::vector<SLinearCluster> ReconstructTracksFromHoughDirection(std::vector<SHit> &hitList, LineEquation houghLine, int trkIter);
        
};

#endif // TPC_SIMPLE_LINES_H