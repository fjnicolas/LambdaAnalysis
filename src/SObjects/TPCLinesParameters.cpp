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
    int Verbose;

    // constructor
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
        Verbose(_verbose)
    {}
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
};


struct HoughAlgorithmPsetType {

    double MaxRadiusLineHypothesis;
    double ThetaRes;
    double MaxDistanceTube;
    int MinHoughHits;
    int Verbose;

    // constructor
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
};


struct TPCLinesAlgoPsetType{

    double MaxRadius;
    double DriftConversion;
    int MaxHoughTracks;
    int MinTrackHits;
    bool RemoveIsolatedHits;
    double MaxNeighbourDistance;
    int MinNeighboursHits;
    int Verbose;
    int DebugMode;
    HoughAlgorithmPsetType HoughAlgorithmPset;
    TrackFinderAlgorithmPsetType TrackFinderAlgorithmPset;
    VertexFinderAlgorithmPsetType VertexFinderAlgorithmPset;
    
    // constructor
    TPCLinesAlgoPsetType(
        double _maxRadius,
        double _driftConversion,
        int _maxHoughTracks,
        int _minTrackHits,
        bool _removeIsolatedHits,
        double _maxNeighbourDistance,
        int _minNeighboursHits,
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
        Verbose(_verbose),
        DebugMode(_debugMode),
        HoughAlgorithmPset(_houghAlgorithmPset),
        TrackFinderAlgorithmPset(_trackFinderPset),
        VertexFinderAlgorithmPset(_vertexFinderPset)
    {}
};


#endif // TPC_LINES_PARAMETERS_H