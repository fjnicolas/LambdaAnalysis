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
    fEdgeOrigin(isEdge)
{};

void SOrigin::AddTrack(SLinearCluster track, SPoint p){
    fTrackList.push_back(track);
    fMultiplicity++;
    SPoint newPoint = SPoint( (fVertex.X()+p.X())/2., (fVertex.Y()+p.Y())/2.  );
    fVertex = newPoint;
}

bool SOrigin::HasTrackIndex(int ix){
    bool hasIx=false;
    for(SLinearCluster &trk: fTrackList){
        if(ix==trk.GetId()) hasIx=true;
    }
    return hasIx;
}

std::ostream& operator<<(std::ostream& out, SOrigin & ori)
{
    out << " SOrigin -- Vertex=(" << ori.GetPoint().X() << ", " <<ori.GetPoint().Y() << ")"<< " Mult="<<ori.Multiplicity()<< " Tracks= ";
    for(auto & trk: ori.GetTracks()){
        out<<trk.GetId()<<" ";
    }
    out<<std::endl;
    return out;
}



SEvent::SEvent(std::vector<SOrigin> origins):
    fOriginList(origins)
{};


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