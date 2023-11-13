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

SOrigin::SOrigin(SPoint p, std::vector<SLinearCluster> tracks, bool isEdge, double yError) :
    fVertex(p),
    fYError(yError),
    fTrackList(tracks),
    fMultiplicity(tracks.size()),
    fNHits(0),
    fEdgeOrigin(isEdge)
{
    for(SLinearCluster & trk:tracks){
        fNHits+=trk.NHits();
    }
}

void SOrigin::AddTrack(SLinearCluster track, SPoint p, double yError){
    fTrackList.push_back(track);
    fMultiplicity++;
    SPoint newPoint = SPoint( (fVertex.X()+p.X())/2., (fVertex.Y()+p.Y())/2.  );
    fVertex = newPoint;
    fYError = std::sqrt( std::pow(fYError,2) + std::pow(yError,2) );
    fNHits+=track.NHits();
}

bool SOrigin::HasTrackIndex(int ix){
    bool hasIx=false;
    for(SLinearCluster &trk: fTrackList){
        if(ix==trk.GetId()) hasIx=true;
    }
    return hasIx;
}

double SOrigin::TotalCharge(){
    double totalCharge = 0;
    for(SLinearCluster &trk: fTrackList){
        for(SHit &hit:trk.GetHits()){
            totalCharge+=hit.Integral();
        }
    }
    return totalCharge;
}


SEvent::SEvent(std::vector<SLinearCluster> tracks, std::vector<SOrigin> origins, std::vector<STriangle> angles,  std::vector<SOrigin> associatedOrigins, double hitDensity):
    fTrackList(tracks),
    fOriginList(origins),
    fAngleList(angles),
    fAssociatedOrigins(associatedOrigins),
    fHitDensity(hitDensity)
{}

int SEvent::GetNOrigins(){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(!ori.IsEdgeOrigin()) continue;
        n++;
    }
    return n;
}

int SEvent::GetNOriginsMult(int mult){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(!ori.IsEdgeOrigin()) continue;
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

int SEvent::GetNOriginsMultGt(int mult, int id1, int id2){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(ori.Multiplicity()>=mult){
            if( !(ori.HasTrackIndex(id1) || ori.HasTrackIndex(id2))){
                n++;
            }
        }
    }
    return n;
}

int SEvent::NHits(){
    int n=0;
    for(SLinearCluster & trk: fTrackList){
        n+=trk.NHits();
    }
    return n;
}