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

SOrigin::SOrigin(SPoint p, bool isEdge, int mult) :
    fVertex(p),
    fEdgeOrigin(isEdge),
    fMultiplicity(mult)
{};


SEvent::SEvent(std::vector<SOrigin> origins):
    fOriginList(origins)
{};
