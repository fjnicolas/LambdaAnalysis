////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleEvents.cpp
//
// \brief Definition of SimpleEvents
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleEvents.h"

SOrigin::SOrigin(SPoint p, std::vector<SLinearCluster> tracks, bool isEdge) :
    fVertex(p),
    fTrackList(tracks),
    fMultiplicity(tracks.size()),
    fNHits(0),
    fEdgeOrigin(isEdge)
{
    for(SLinearCluster & trk:tracks){
        fNHits+=trk.NHits();
    }
}

void SOrigin::AddTrack(SLinearCluster track, SPoint p){
    fTrackList.push_back(track);
    fMultiplicity++;
    SPoint newPoint = SPoint( (fVertex.X()+p.X())/2., (fVertex.Y()+p.Y())/2.  );
    fVertex = newPoint;
    fNHits+=track.NHits();
}

bool SOrigin::HasTrackIndex(int ix){
    bool hasIx=false;
    for(SLinearCluster &trk: fTrackList){
        if(ix==trk.GetId()) hasIx=true;
    }
    return hasIx;
}


SEvent::SEvent(std::vector<SLinearCluster> tracks, std::vector<SOrigin> origins, std::vector<STriangle> angles,  std::vector<SOrigin> associatedOrigins, double hitDensity):
    fTrackList(tracks),
    fOriginList(origins),
    fAngleList(angles),
    fAssociatedOrigins(associatedOrigins),
    fHitDensity(hitDensity)
{}


int SEvent::GetNOriginsMult(int mult){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(ori.Multiplicity()==mult){
            n++;
        }
    }
    return n;
}

int SEvent::GetNOriginsMultGt(int mult){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(ori.Multiplicity()>=mult){
            n++;
        }
    }
    return n;
}