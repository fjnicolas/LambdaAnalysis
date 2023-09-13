////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleEvents.h
//
// \brief Definition of SimpleEvent
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_EVENT_H
#define TPC_SIMPLE_EVENT_H

#include <iostream>
#include <string>
#include <iosfwd>


#include "TPCSimpleLines.h"
#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleTriangles.h"



class SOrigin {
    private:
        SPoint fVertex;
        std::vector<SLinearCluster> fTrackList;
        int fMultiplicity;
        bool fEdgeOrigin;
        
        
    public:
        SOrigin(SPoint p, std::vector<SLinearCluster> tracks, bool isEdge);
        
        SPoint GetPoint(){ return fVertex;};
        std::vector<SLinearCluster> GetTracks(){ return fTrackList; };
        int Multiplicity(){return fMultiplicity;};
        void AddTrack(SLinearCluster track, SPoint p); 
        bool HasTrackIndex(int ix);

        friend std::ostream& operator<<(std::ostream& out, SOrigin const& ori);

};



class SEvent {
    private:
        std::vector<STriangle> fAngleList;
        std::vector<SLinearCluster> fTrackList;
        std::vector<SOrigin> fOriginList;
    
        
    public:
        SEvent(std::vector<SOrigin> origins);

        std::vector<SOrigin> GetOrigins(){return fOriginList;};
        int GetNOrigins(){ return fOriginList.size();};

        int GetNOriginsMult(int mult);
        int GetNOriginsMultGt(int mult);
        
};

#endif // TPC_SIMPLE_HITS_H