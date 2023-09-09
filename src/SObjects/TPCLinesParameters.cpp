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
    float HitDensityThreshold;
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


class SEventId {
    public:
        SEventId(int r=0, int sr=0, int e=0):
            fRun(r),
            fSubRun(sr),
            fEvent(e),
            fEventLabel(std::to_string(r)+"_"+std::to_string(sr)+"_"+std::to_string(e))
            {};


        SEventId(const SEventId& other) : SEventId(other.fRun, other.fSubRun, other.fEvent) {}

        std::string Label(){return fEventLabel;};
        int Run(){return fRun;};
        int SubRun(){return fSubRun;};
        int Event(){return fEvent;};

        friend std::ostream& operator<<(std::ostream& os, const SEventId& statusPrinter) {
        os << " ***** R=" << statusPrinter.fRun << " SR=" << statusPrinter.fSubRun << " E=" << statusPrinter.fEvent << " *****" << std::endl;
        return os;
    }

    private:
        int fRun, fSubRun, fEvent;
        std::string fEventLabel;
};

class SEventSelection {
    public:
        SEventSelection(SEventId ev):
            fEvent(ev),
            fNSlices(0),
            fNSelected(0),
            fNNotSelected(0),
            fNSkipped(0)
            {};

        SEventSelection(){
            fEvent = SEventId();
            fNSlices=0;
            fNSelected=0;
            fNNotSelected=0;
            fNSkipped=0;
        }
        
        SEventSelection(const SEventSelection& other) : SEventSelection(other.fEvent) {}

        SEventId EventId(){return fEvent;};

        void AddSelected(){
            fNSlices++;
            fNSelected++;
        };

        void AddSkipped(){
            fNSlices++;
            fNSkipped++;
        }

        void AddNotSelected(){
            fNSlices++;
            fNNotSelected++;
        }

        int NSlices(){return fNSlices;};
        int NSelected(){return fNSelected;};
        int NNotSelected(){return fNNotSelected;};
        int NSkipped(){return fNSkipped;};

    private:
        SEventId fEvent;
        int fNSlices;
        int fNSelected;
        int fNNotSelected;
        int fNSkipped;
};


class EfficiencyCalculator {
public:
    EfficiencyCalculator():
        fEventList({}),
        nEvents(0),
        nEventsSkipped(0),
        nProcessedEvents(0),
        nEventsSelected(0)
        {};

    void AddEvent(SEventId ev){
        if(fEventList.find(ev.Label()) == fEventList.end()){
            SEventSelection sEv(ev);
            fEventList[ev.Label()] = SEventSelection(sEv);
        }
    }

    void UpdateSkipped(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddSkipped();
    }

    void UpdateSelected(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddSelected();
    }

    void UpdateNotSelected(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddNotSelected();

    }

    void UpdateValues(){
        nEvents = fEventList.size();
        nEventsSkipped=0;
        nProcessedEvents=0;
        nEventsSelected=0;
        for(auto &pair:fEventList){
            if(pair.second.NSkipped()==pair.second.NSlices()){
                nEventsSkipped++;
            }
            else{
                nProcessedEvents++;
                if(pair.second.NSelected()>0){
                    nEventsSelected++;
                }
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os, EfficiencyCalculator& statusPrinter) {
        statusPrinter.UpdateValues();
        os << " ********** Final status-" << std::endl;
        os << " ... NTotalEvents=" << statusPrinter.nEvents << " NSkipped=" << statusPrinter.nEventsSkipped << std::endl;
        os << " ... NProcessed=" << statusPrinter.nProcessedEvents << " NSelected=" << statusPrinter.nEventsSelected;
        os << " efficiency=" << static_cast<double>(statusPrinter.nEventsSelected) / statusPrinter.nProcessedEvents;
        os << " efficiency all=" << static_cast<double>(statusPrinter.nEventsSelected) / statusPrinter.nEvents << std::endl;
        return os;
    }

private:
    std::map<std::string, SEventSelection> fEventList;
    int nEvents;
    int nEventsSkipped;
    int nProcessedEvents;
    int nEventsSelected;
};


#endif // TPC_LINES_PARAMETERS_H