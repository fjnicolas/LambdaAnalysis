////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleLines.cpp
//
// \brief Definition of TPCSimpleLines
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleLines.h"


// Class functions for LineEquation class 
LineEquation::LineEquation(float slope, float intercept):
    fM(slope),
    fN(intercept)
{}

SPoint LineEquation::GetLineClosestPoint(double a, double b, double c, SPoint P) {
    
    // line in form ax+by+c=0
    float x = (b * (b * P.X() - a * P.Y()) - a * c) / (a * a + b * b);
    float y = (a * (-b * P.X() + a * P.Y()) - b * c) / (a * a + b * b);
    
    return SPoint(x, y);
}


SPoint LineEquation::GetLineClosestPoint(SPoint P) {
    
    return GetLineClosestPoint(-fM, 1, -fN, P);
}


float LineEquation::GetDistance(SPoint p) {
    SPoint pProj = GetLineClosestPoint(-fM, 1, -fN, p);
    return std::sqrt((p.X() - pProj.X()) * (p.X() - pProj.X()) + (p.Y() - pProj.Y()) * (p.Y() - pProj.Y()));
}

float LineEquation::Evaluate(SPoint p) {
    return fM * p.X() + fN;
}


float LineEquation::EvaluateX(double x) {
    return fM * x + fN;
}


// Class functions for HoughEquation class 
HoughLine::HoughLine(LineEquation line, float score, int nHits):
    fEquation(line),
    fScore(score),
    fNHits(nHits)
{}