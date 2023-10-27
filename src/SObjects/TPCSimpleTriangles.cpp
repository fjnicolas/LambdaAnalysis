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

STriangle::STriangle(SPoint main_vertex, SPoint vertex_b, SPoint vertex_c, SHit mainhit, SLinearCluster track2, SLinearCluster track1, SLinearCluster mainTrack, double weight_b, double weight_c)
{
    fTrack1 = track1;
    fTrack2 = track2;
    fMainTrack = mainTrack;
    
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

    // Calculate the opening angle in radians
    fOpeningAngle = std::abs(180 * ( std::atan((slope2 - slope1) / (1 + slope1 * slope2)) ) / M_PI);

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

/*double ComputeCoveredArea2(std::vector<SHit> triangleHits, double widthTol) {
    // Determine the line equation of the opposite side
    SPoint P1 = (fVertexB.X() < fVertexC.X()) ? fVertexB : fVertexC;
    SPoint P2 = (fVertexB.X() > fVertexC.X()) ? fVertexB : fVertexC;

    double m = (P2.Y() - P1.Y()) / (P2.X() - P1.X());
    double n = P1.Y() - m * P1.X();

    int N = 0;
    int NCovered = 0;

    for (int x = static_cast<int>(P1.X()); x < static_cast<int>(P2.X()); x++) {
        int y = static_cast<int>(m * x + n);
        N++;

        for (const SHit &hit : triangleHits) {
            if (x == static_cast<int>(hit.X()) &&
                y >= hit.Y() - widthTol*hit.Width() && y <= hit.Y() + widthTol*hit.Width()) {
                NCovered++;
                break; // No need to continue checking this x coordinate.
            }
        }
    }

    double coverageFraction = static_cast<double>(NCovered) / N;

    return coverageFraction;
}*/

struct Point {
    double x, y;
};

struct Line {
    double x, y, yWidth;
};

double interpolateX(const Point& p1, const Point& p2, double y) {
    double x1 = p1.x, y1 = p1.y;
    double x2 = p2.x, y2 = p2.y;

    if (y1 == y2) {
        return std::min(x1, x2);
    }

    return x1 + (y - y1) * (x2 - x1) / (y2 - y1);
}

double fractionOfLineInTriangle(const std::vector<Point>& triangle, const Line& line) {
    double topY = line.y - line.yWidth / 2;
    double bottomY = line.y + line.yWidth / 2;
    double totalLength = 0;

    for (double y = topY; y <= bottomY; y += 1.0) {
        double x1 = interpolateX(triangle[0], triangle[1], y);
        double x2 = interpolateX(triangle[1], triangle[2], y);
        double x3 = interpolateX(triangle[2], triangle[0], y);

        double min_x = std::min(std::min(x1, x2), x3);
        double max_x = std::max(std::max(x1, x2), x3);

        // Check if the line segment overlaps with the triangle for the given y-coordinate.
        if (min_x <= line.x && max_x >= line.x) {
            // Calculate the length of the line segment within the triangle.
            double segmentLength = std::min(max_x, line.x + line.yWidth / 2) - std::max(min_x, line.x - line.yWidth / 2);
            totalLength += segmentLength;
        }
    }

    return totalLength / line.yWidth;
}


double STriangle::ComputeCoveredArea(std::vector<SHit> triangleHits, double widthTol) {
   
    std::vector<Point> triangle = { {fMainVertex.X(), fMainVertex.Y()}, {fVertexB.X(), fVertexB.Y()}, {fVertexC.X(), fVertexC.Y()} };
    std::vector<Line> verticalLines = {{1.0, 3.0}, {2.0, 2.0}, {3.0, 2.0}};

    for (SHit &hit : triangleHits){
        Line line = {hit.X(), hit.Y(), hit.Width()};
        double fraction = fractionOfLineInTriangle(triangle, line);
    }

    //sdouble coverageFraction = static_cast<double>(NCovered) / N;

    return 1;
}
