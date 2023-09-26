///////////////////////////////////////////////////////////////////////
/// File: ChargeDensityAlgConf.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////

#ifndef SBND_CHARGEDENSITYPSET_H
#define SBND_CHARGEDENSITYPSET_H

struct FRAMSPsetType {
  bool ApplyRawSmoothing;
  bool ApplySmoothing;
  bool ApplyCumulativeSmoothing;
  unsigned int NDriftPack;
  unsigned int NWirePack;
  float ExpoAvSmoothPar;
  int UnAvNeighbours;
  double CumulativeCut;
  int SlidingWindowN;
  int MaxRadius;
  int NSamplesBeginSlope;
  int Verbose;
  bool CalculateScore;
  std::string TMVAFilename;

  // Constructor to initialize the members
  FRAMSPsetType(bool rawSmoothing, bool smoothing, bool cumulativeSmoothing,
          unsigned int driftPack, unsigned int wirePack, float smoothParam,
          int unweightedNeighbors, double cut, int windowN, int radius,
          int samplesSlope, int verbose, bool calculate, const std::string& filename)
      : ApplyRawSmoothing(rawSmoothing),
        ApplySmoothing(smoothing),
        ApplyCumulativeSmoothing(cumulativeSmoothing),
        NDriftPack(driftPack),
        NWirePack(wirePack),
        ExpoAvSmoothPar(smoothParam),
        UnAvNeighbours(unweightedNeighbors),
        CumulativeCut(cut),
        SlidingWindowN(windowN),
        MaxRadius(radius),
        NSamplesBeginSlope(samplesSlope),
        Verbose(verbose),
        CalculateScore(calculate),
        TMVAFilename(filename) {
  }
};



#endif
