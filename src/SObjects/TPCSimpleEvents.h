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
        int fNHits;
        bool fEdgeOrigin;
        
        
        
    public:
        SOrigin(SPoint p, std::vector<SLinearCluster> tracks, bool isEdge);
        
        SPoint GetPoint(){ return fVertex;};
        std::vector<SLinearCluster> GetTracks(){ return fTrackList; };
        int Multiplicity(){return fMultiplicity;};
        void AddTrack(SLinearCluster track, SPoint p); 
        bool HasTrackIndex(int ix);
        int NHits(){return fNHits;};
        SLinearCluster GetTrackEntry(size_t ix){ return fTrackList[ix];};

        //friend std::ostream& operator<<(std::ostream& out, SOrigin const& ori);

        friend std::ostream& operator<<(std::ostream& out, SOrigin & ori)
        {
            out << " SOrigin -- Vertex=(" << ori.GetPoint().X() << ", " <<ori.GetPoint().Y() << ")"<< " Mult="<<ori.Multiplicity()<< " Tracks= ";
            for(auto & trk: ori.GetTracks()){
                out<<trk.GetId()<<" ";
            }
            out<<"  NHits"<<ori.NHits()<<std::endl;
            return out;
        }
};



class SEvent {
    private:
        std::vector<SLinearCluster> fTrackList;
        std::vector<SOrigin> fOriginList;
        std::vector<STriangle> fAngleList;
        std::vector<SOrigin> fAssociatedOrigins;
        double fHitDensity;
    
        
    public:
        SEvent(std::vector<SLinearCluster> tracks={}, std::vector<SOrigin> origins={}, std::vector<STriangle> angles = {}, std::vector<SOrigin> associatedOrigins = {}, double hitDensity=0);

        std::vector<SOrigin> GetOrigins(){return fOriginList;};
        int GetNOrigins(){ return fOriginList.size();};
        int GetNOriginsMult(int mult);
        int GetNOriginsMultGt(int mult);

        int GetNAngles(){ return fAngleList.size();};

        std::vector<STriangle> GetAngleList(){return fAngleList;};
        std::vector<SOrigin> GetAssociatedOrigins(){return fAssociatedOrigins;};

        double HitDensity(){return fHitDensity;};
        
};

#endif // TPC_SIMPLE_HITS_H