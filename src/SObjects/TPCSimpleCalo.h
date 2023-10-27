////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////


#ifndef TPC_SIMPLE_CALO_H
#define TPC_SIMPLE_CALO_H


#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "TPCSimpleHits.h"


class SCalo {
public:
    // Constructor to initialize the collection
    SCalo(const std::vector<SHit>& points);

    // Method to display the calculated path lengths
    void Display();

    std::vector<double> GetResidualRange() const {return fResidualRange;};
    std::vector<double> GetDepositedEnergy() const {return fDepositedEnergy;};
    double GetTrackLength() const {return fTrackLength;};

private:

    // Conversion factors
    double fHitIntegralToEnergy;
    double fStepXToLength;
    double fStepYToLength;

    // Track legnth
    double fTrackLength;
    // Define the hit points and path lengths as private members
    std::vector<SHit> fHitList;
    std::vector<double> fPathLengths;
    std::vector<double> fDepositedEnergy;
    std::vector<double> fResidualRange;

    // Function to calculate the distance between two hit points
    double CalculateDistance(const SHit& p1, const SHit& p2);

    // Method to calculate and store the path lengths
    void CalculatePathLengths();

    // Function to calculate residual range and deposited energy
    void CalculateResidualRange();
};


void CreateEnergyLossVsResidualRangePlot(const std::vector<SCalo>& caloObjects);

#endif // TPC_SIMPLE_CALO_H