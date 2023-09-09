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
#include "TPCSimpleTriangles.h"



class SOrigin {
    private:
        SPoint fVertex;
        bool fEdgeOrigin;
        int fMultiplicity;
    
        
    public:
        SOrigin(SPoint p, bool isEdge, int mult);
    
        
};



class SEvent {
    private:
        std::vector<STriangle> fAngleList;
        std::vector<SLinearCluster> fTrackList;
        std::vector<SOrigin> fOriginList;
    
        
    public:
        SEvent(std::vector<SOrigin> origins);

        int GetNOrigins(){ return fOriginList.size();};
    
        
};

#endif // TPC_SIMPLE_HITS_H