////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleTriangle.h
//
// \brief Definition of SimpleTriangle
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_TRIANGLE_H
#define TPC_SIMPLE_TRIANGLE_H

#include <iostream>
#include <string>
#include <iosfwd>


#include "TPCSimpleLines.h"
#include "TPCSimpleHits.h"



class STriangle {
    private:
        SPoint fMainVertex;
        SHit fMainVertexHit;
        SPoint fVertexB;
        SPoint fVertexC;
        SPoint fMidPoint;
        SPoint fDirectorVector;
        LineEquation fDirection;
        LineEquation fMomentumHypo1;
        LineEquation fMomentumHypo2;
        
    public:
        STriangle(SPoint main_vertex, SPoint vertex_b, SPoint vertex_c, SHit mainhit, double weight_b=1, double weight_c=1);
        
        SPoint GetMainVertex() const {
            return fMainVertex;
        }

        SHit GetMainVertexHit() const {
            return fMainVertexHit;
        }

        SPoint GetVertexB() const {
            return fVertexB;
        }

        SPoint GetVertexC() const {
            return fVertexC;
        }

        SPoint GetMidPoint() const {
            return fMidPoint;
        }

        SPoint GetDirectorVector() const {
            return fDirectorVector;
        }

        LineEquation GetDirection() const {
            return fDirection;
        }

        LineEquation GetMomentumHypo1() const {
            return fMomentumHypo1;
        }

        LineEquation GetMomentumHypo2() const {
            return fMomentumHypo2;
        }

        
};

#endif // TPC_SIMPLE_HITS_H