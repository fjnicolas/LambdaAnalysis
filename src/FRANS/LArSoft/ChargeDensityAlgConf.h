///////////////////////////////////////////////////////////////////////
/// File: ChargeDensityAlgConf.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////
#include <iostream>


#include "fhiclcpp/types/Atom.h"


#ifndef SBND_CHARGEDENSITYALGCONF_H
#define SBND_CHARGEDENSITYALGCONF_H

namespace ChargeDensityConf{


  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<bool> applyRawSmoothing {
      Name("ApplyRawSmoothing"),
      Comment("Apply exponential average smoothing to Zeta(rho)")
    };

    fhicl::Atom<bool> applySmoothing {
      Name("ApplySmoothing"),
      Comment("Apply exponential average smoothing to Zeta(rho) derivative")
    };

    fhicl::Atom<bool> applyCumulativeSmoothing {
      Name("ApplyCumulativeSmoothing"),
      Comment("Apply exponential average smoothing to Zeta(rho) cumulative")
    };

    fhicl::Atom<unsigned int> nDriftPack {
      Name("NDriftPack"),
      Comment("Number of packed time ticks")
    };

    fhicl::Atom<unsigned int> nWirePack {
      Name("NWirePack"),
      Comment("Number of packed wires")
    };

    fhicl::Atom<float> expoAvSmoothPar {
      Name("ExpoAvSmoothPar"),
      Comment("Expo average smoothing paramter")
    };

    fhicl::Atom<int> unAvNeighbours{
      Name("UnAvNeighbours"),
      Comment("Neighbours used for unweigh smoothing")
    };

    fhicl::Atom<double> cumulativeCut {
      Name("CumulativeCut"),
      Comment("Cumulative vut")
    };

    fhicl::Atom<int> slidingWindowN {
      Name("SlidingWindowN"),
      Comment("Neighbours used in sliding window")
    };

    fhicl::Atom<int> maxRadius {
      Name("MaxRadius"),
      Comment("Maximum vall radius")
    };

    fhicl::Atom<int> nSamplesBeginSlope {
      Name("NSamplesBeginSlope"),
      Comment("Binsused to compute start slope")
    };

    fhicl::Atom<bool> debugMode {
      Name("DebugMode"),
      Comment("Debug mode")
    };

    fhicl::Atom<bool> calculateScore {
      Name("CalculateScore"),
      Comment("Calculate final BDT score")
    };

    fhicl::Atom<std::string> tMVAFilename {
      Name("TMVAFilename"),
      Comment("BDT input file")
    };

  };    //struct Config


  struct FRAMSPsetType {
    bool applyRawSmoothing;
    bool applySmoothing;
    bool applyCumulativeSmoothing;
    unsigned int nDriftPack;
    unsigned int nWirePack;
    float expoAvSmoothPar;
    int unAvNeighbours;
    double cumulativeCut;
    int slidingWindowN;
    int maxRadius;
    int nSamplesBeginSlope;
    bool debugMode;
    bool calculateScore;
    std::string tMVAFilename;

    // Constructor to initialize the members
    FRAMSPsetType(bool rawSmoothing, bool smoothing, bool cumulativeSmoothing,
           unsigned int driftPack, unsigned int wirePack, float smoothParam,
           int unweightedNeighbors, double cut, int windowN, int radius,
           int samplesSlope, bool debug, bool calculate, const std::string& filename)
        : applyRawSmoothing(rawSmoothing),
          applySmoothing(smoothing),
          applyCumulativeSmoothing(cumulativeSmoothing),
          nDriftPack(driftPack),
          nWirePack(wirePack),
          expoAvSmoothPar(smoothParam),
          unAvNeighbours(unweightedNeighbors),
          cumulativeCut(cut),
          slidingWindowN(windowN),
          maxRadius(radius),
          nSamplesBeginSlope(samplesSlope),
          debugMode(debug),
          calculateScore(calculate),
          tMVAFilename(filename) {
    }
};

}

#endif
