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
#include <map>

#include "TPCSimpleLines.h"
#include "TPCSimpleClusters.h"
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

        SLinearCluster fTrack1;
        SLinearCluster fTrack2;
        SLinearCluster fMainTrack;

        double fOpeningAngle;
        
    public:
        STriangle(SPoint main_vertex, SPoint vertex_b, SPoint vertex_c, SHit mainhit, SLinearCluster track2, SLinearCluster track1, SLinearCluster mainTrack, double weight_b=1, double weight_c=1);

        SLinearCluster GetTrack1() const {
            return fTrack1;
        }

        SLinearCluster GetTrack2() const {
            return fTrack2;
        }

        SLinearCluster GetMainTrack() const {
            return fMainTrack;
        }
        
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

        double GetOpeningAngle() const {
            return fOpeningAngle;
        }

        double GetSideLenghtB() const {
            return std::hypot( fMainVertex.X() - fVertexB.X(), fMainVertex.Y() - fVertexB.Y() );
        }

        double GetSideLenghtC() const {
            return std::hypot( fMainVertex.X() - fVertexC.X(), fMainVertex.Y() - fVertexC.Y() );
        }

        double GetOppositeSideLenght() const {
            return std::hypot( fVertexB.X() - fVertexC.X(), fVertexB.Y() - fVertexC.Y() );
        }

        

        double ComputeCoveredArea(std::vector<SHit> triangleHits, double widthTol);

        
};

#endif // TPC_SIMPLE_HITS_H