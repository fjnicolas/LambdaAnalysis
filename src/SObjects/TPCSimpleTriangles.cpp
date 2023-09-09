////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleTriangle.cpp
//
// \brief Definition of SimpleTriangle
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleTriangles.h"

double get_momentum(double ke, double m) {
    return ke * std::sqrt(1 + 2. * m / ke);
}

STriangle::STriangle(SPoint main_vertex, SPoint vertex_b, SPoint vertex_c, SHit mainhit, double weight_b, double weight_c)
{
    fMainVertex = main_vertex;
    fMainVertexHit = mainhit;
    
    if (vertex_b.X() < vertex_c.X()) {
        fVertexB = vertex_b;
        fVertexC = vertex_c;
    } else {
        fVertexB = vertex_c;
        fVertexC = vertex_b;
    }
    
    double xM = (fVertexB.X() + fVertexC.X()) / 2;
    double yM = (fVertexB.Y() + fVertexC.Y()) / 2;
    fMidPoint = SPoint(xM, yM);
    //std::cout << "MidPoint " << xM << " " << yM << std::endl;
    fDirectorVector = SPoint( fMainVertex.X() - xM, fMainVertex.Y() - yM );
    double slope = (yM - fMainVertex.Y()) / (xM - fMainVertex.X());
    double intercept = fMainVertex.Y() - slope * fMainVertex.X();
    fDirection = LineEquation(slope, intercept);

    double slopeB = (fVertexB.Y() - fMainVertex.Y()) / (fVertexB.X() - fMainVertex.X());
    double slopeC = (fVertexC.Y() - fMainVertex.Y()) / (fVertexC.X() - fMainVertex.X());

    double MProton = 938;
    double MPion = 130;

    double fConFactor = (1 / 0.0201293) * 23.6e-6;

    //std::cout << "Weights before conversion: " << weight_b << ", " << weight_c << std::endl;
    //std::cout << "Weights after conversion: " << weight_b * fConFactor << ", " << weight_c * fConFactor << std::endl;

    // Hypothesis 1
    double slope1 = (MProton * slopeB + MPion * slopeC) / (MProton + MPion);
    double intercept1 = fMainVertex.Y() - slope1 * fMainVertex.X();
    fMomentumHypo1 = LineEquation(slope1, intercept1);

    // Hypothesis 2
    double slope2 = (MPion * slopeB + MProton * slopeC) / (MProton + MPion);
    double intercept2 = fMainVertex.Y() - slope2 * fMainVertex.X();
    fMomentumHypo2 = LineEquation(slope2, intercept2);

    // Calculate momenta
    weight_b *= fConFactor;
    weight_c *= fConFactor;

    double Pb = get_momentum(weight_b, MProton);
    double Pc = get_momentum(weight_c, MPion);
    //std::cout << "Momentum weights hypo 1: " << Pb << ", " << Pc << std::endl;

    slope1 = (Pb * slopeB + Pc * slopeC) / (Pb + Pc);
    intercept1 = fMainVertex.Y() - slope1 * fMainVertex.X();
    fMomentumHypo1 = LineEquation(slope1, intercept1);


    Pb = get_momentum(weight_b, MPion);
    Pc = get_momentum(weight_c, MProton);
    //std::cout << "Momentum weights hypo 2: " << Pb << ", " << Pc << std::endl;

    slope2 = (Pb * slopeB + Pc * slopeC) / (Pb + Pc);
    intercept2 = fMainVertex.Y() - slope2 * fMainVertex.X(

    );
    fMomentumHypo2 = LineEquation(slope2, intercept2);

}

