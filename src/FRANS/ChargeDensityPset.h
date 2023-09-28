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
  int NSamplesBeginSlope;
  int MaxRadius;
  int Verbose;
  bool CalculateScore;
  std::string TMVAFilename;
  std::string OutputPath;

  // Constructor to initialize the members
  FRAMSPsetType(bool rawSmoothing=false,
                bool smoothing=false,
                bool cumulativeSmoothing=false,
                unsigned int driftPack=4,
                unsigned int wirePack=1,
                float smoothParam=0.3,
                int unweightedNeighbors=1,
                double cut=0.8,
                int windowN=3,
                int samplesSlope=3,
                int radius=70,
                int verbose=0,
                bool calculate=false,
                const std::string& filename="",
                const std::string& outputpath="" )
                  : ApplyRawSmoothing(rawSmoothing),
                    ApplySmoothing(smoothing),
                    ApplyCumulativeSmoothing(cumulativeSmoothing),
                    NDriftPack(driftPack),
                    NWirePack(wirePack),
                    ExpoAvSmoothPar(smoothParam),
                    UnAvNeighbours(unweightedNeighbors),
                    CumulativeCut(cut),
                    SlidingWindowN(windowN),
                    NSamplesBeginSlope(samplesSlope),
                    MaxRadius(radius),
                    Verbose(verbose),
                    CalculateScore(calculate),
                    TMVAFilename(filename),
                    OutputPath(outputpath)
  {}
};



#endif
