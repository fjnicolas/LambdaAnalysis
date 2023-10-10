////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_PARAMETERS_H
#define TPC_LINES_PARAMETERS_H

#include <vector>

// ROOT includes
#include "TH2F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"

#include "TPCSimpleEvents.h"

std::vector<TPad*> buildpadcanvas(int nx, int ny);

struct TrackFinderAlgorithmPsetType {
    int MaxDTube;
    double MaxDCluster;
    bool SingleWireMode;
    int MinClusterHits;
    double DCleaning;
    double ClusterCompletenessCut;
    double ClusterAngleCut;
    bool CaptureMissingHits;
    int MinTrackHits;
    float HitDensityThreshold;
    int Verbose;

    // constructor
    TrackFinderAlgorithmPsetType(){};

    TrackFinderAlgorithmPsetType(
        int _maxDTube,
        double _maxDCluster,
        bool _singleWireMode,
        int _minClusterHits,
        double _dCleaning,
        double _clusterCompletenessCut,
        double _clusterAngleCut,
        bool _captureMissingHits,
        int _minTrackHits,
        float _hitDensityThreshold,
        int _verbose) : 
        MaxDTube(_maxDTube),
        MaxDCluster(_maxDCluster),
        SingleWireMode(_singleWireMode),
        MinClusterHits(_minClusterHits), 
        DCleaning(_dCleaning),
        ClusterCompletenessCut(_clusterCompletenessCut),
        ClusterAngleCut(_clusterAngleCut), 
        CaptureMissingHits(_captureMissingHits),
        MinTrackHits(_minTrackHits),
        HitDensityThreshold(_hitDensityThreshold),
        Verbose(_verbose)
    {}

    void Print() const {
        std::cout << "MaxDTube: " << MaxDTube << std::endl;
        std::cout << "MaxDCluster: " << MaxDCluster << std::endl;
        std::cout << "SingleWireMode: " << (SingleWireMode ? "true" : "false") << std::endl;
        std::cout << "MinClusterHits: " << MinClusterHits << std::endl;
        std::cout << "DCleaning: " << DCleaning << std::endl;
        std::cout << "ClusterCompletenessCut: " << ClusterCompletenessCut << std::endl;
        std::cout << "ClusterAngleCut: " << ClusterAngleCut << std::endl;
        std::cout << "CaptureMissingHits: " << (CaptureMissingHits ? "true" : "false") << std::endl;
        std::cout << "MinTrackHits: " << MinTrackHits << std::endl;
        std::cout << "HitDensityThreshold: " << HitDensityThreshold << std::endl;
        std::cout << "Verbose: " << Verbose << std::endl;
    }
};


struct VertexFinderAlgorithmPsetType {
 
    double MaxDistToEdge = 3;
    bool RefineVertexIntersection = true;
    bool UseEdgesDiscard = true;
    float MaxTrackFractionInMain = 0.75;
    bool DecideMainTrack = false;
    bool AddCollinearLines = false;
    int Verbose;

    // constructor
    VertexFinderAlgorithmPsetType(){};
    
    VertexFinderAlgorithmPsetType(
        double _maxDistToEdge,
        bool _refineVertexIntersection,
        bool _useEdgesDiscard,
        float _maxTrackFractionInMain,
        bool _decideMainTrack,
        bool _addCollinearLines,
        int _verbose): 
        MaxDistToEdge(_maxDistToEdge),
        RefineVertexIntersection(_refineVertexIntersection),
        UseEdgesDiscard(_useEdgesDiscard),
        MaxTrackFractionInMain(_maxTrackFractionInMain),
        DecideMainTrack(_decideMainTrack),
        AddCollinearLines(_addCollinearLines),
        Verbose(_verbose)
    {}

    void Print() const {
        std::cout << "MaxDistToEdge: " << MaxDistToEdge << std::endl;
        std::cout << "RefineVertexIntersection: " << (RefineVertexIntersection ? "true" : "false") << std::endl;
        std::cout << "UseEdgesDiscard: " << (UseEdgesDiscard ? "true" : "false") << std::endl;
        std::cout << "MaxTrackFractionInMain: " << MaxTrackFractionInMain << std::endl;
        std::cout << "DecideMainTrack: " << (DecideMainTrack ? "true" : "false") << std::endl;
        std::cout << "AddCollinearLines: " << (AddCollinearLines ? "true" : "false") << std::endl;
        std::cout << "Verbose: " << Verbose << std::endl;
    }
};


struct HoughAlgorithmPsetType {

    double MaxRadiusLineHypothesis;
    double ThetaRes;
    double MaxDistanceTube;
    int MinHoughHits;
    int Verbose;

    // constructor
    HoughAlgorithmPsetType(){};

    HoughAlgorithmPsetType(
        double _maxRadiusLineHypothesis,
        double _thetaRes,
        double _maxDistanceTube,
        int _minHoughHits,
        int _verbose) : 
        MaxRadiusLineHypothesis(_maxRadiusLineHypothesis),
        ThetaRes(_thetaRes),
        MaxDistanceTube(_maxDistanceTube), 
        MinHoughHits(_minHoughHits),
        Verbose(_verbose)
    {}

    void Print() const {
        std::cout << "MaxRadiusLineHypothesis: " << MaxRadiusLineHypothesis << std::endl;
        std::cout << "ThetaRes: " << ThetaRes << std::endl;
        std::cout << "MaxDistanceTube: " << MaxDistanceTube << std::endl;
        std::cout << "MinHoughHits: " << MinHoughHits << std::endl;
        std::cout << "Verbose: " << Verbose << std::endl;
    }
};


struct TPCLinesAlgoPsetType{

    double MaxRadius;
    double DriftConversion;
    int MaxHoughTracks;
    int MinTrackHits;
    bool RemoveIsolatedHits;
    double MaxNeighbourDistance;
    int MinNeighboursHits;
    int VertexAlgorithm;
    std::string View;
    std::string OutputPath;
    int Verbose;
    int DebugMode;
    HoughAlgorithmPsetType HoughAlgorithmPset;
    TrackFinderAlgorithmPsetType TrackFinderAlgorithmPset;
    VertexFinderAlgorithmPsetType VertexFinderAlgorithmPset;
    
    // constructor
    TPCLinesAlgoPsetType(){};

    TPCLinesAlgoPsetType(
        double _maxRadius,
        double _driftConversion,
        int _maxHoughTracks,
        int _minTrackHits,
        bool _removeIsolatedHits,
        double _maxNeighbourDistance,
        int _minNeighboursHits,
        int _vertexAlgorithm,
        std::string _view,
        std::string _outputPath,
        int _verbose,
        int _debugMode,
        HoughAlgorithmPsetType _houghAlgorithmPset,
        TrackFinderAlgorithmPsetType _trackFinderPset,
        VertexFinderAlgorithmPsetType _vertexFinderPset) :
        MaxRadius(_maxRadius),
        DriftConversion(_driftConversion),
        MaxHoughTracks(_maxHoughTracks),
        MinTrackHits(_minTrackHits),
        RemoveIsolatedHits(_removeIsolatedHits),
        MaxNeighbourDistance(_maxNeighbourDistance), 
        MinNeighboursHits(_minNeighboursHits),
        VertexAlgorithm(_vertexAlgorithm),
        View(_view),
        OutputPath(_outputPath),
        Verbose(_verbose),
        DebugMode(_debugMode),
        HoughAlgorithmPset(_houghAlgorithmPset),
        TrackFinderAlgorithmPset(_trackFinderPset),
        VertexFinderAlgorithmPset(_vertexFinderPset)
    {}


    void Print() {
        std::cout << "MaxRadius: " <<  MaxRadius << std::endl;
        std::cout << "DriftConversion: " <<  DriftConversion << std::endl;
        std::cout << "MaxHoughTracks: " <<  MaxHoughTracks << std::endl;
        std::cout << "MinTrackHits: " <<  MinTrackHits << std::endl;
        std::cout << "RemoveIsolatedHits: " <<  RemoveIsolatedHits << std::endl;
        std::cout << "MaxNeighbourDistance: " <<  MaxNeighbourDistance << std::endl;
        std::cout << "MinNeighboursHits: " <<  MinNeighboursHits << std::endl;
        std::cout << "Verbose: " <<  Verbose << std::endl;
        std::cout << "DebugMode: " <<  DebugMode << std::endl;

        HoughAlgorithmPset.Print();
        TrackFinderAlgorithmPset.Print();
        VertexFinderAlgorithmPset.Print();
    }
};



// Function to print the struct's values
void printStruct(const TPCLinesAlgoPsetType&  );

#endif // TPC_LINES_PARAMETERS_H