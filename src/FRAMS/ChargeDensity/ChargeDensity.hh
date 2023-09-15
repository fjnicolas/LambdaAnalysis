///////////////////////////////////////////////////////////////////////
/// File: ChargeDensity.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////
#include <vector>
#include <limits>
#include <numeric>
#include <map>
#include <iostream>

#include "TGraphErrors.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art_root_io/TFileDirectory.h"

#include "lardataobj/RecoBase/Hit.h"

#include "sbndcode/FRAMS/ChargeDensity/VertexView.hh"
#include "sbndcode/FRAMS/ChargeDensity/ChargeDensityAlgConf.hh"
#include "TMVA/Reader.h"

#define DefaultMaxZSize 2000

#ifndef SBND_CHARGEDENSITY_H
#define SBND_CHARGEDENSITY_H

class ChargeDensity{

  public:
    ChargeDensity(ChargeDensity const&) = delete;
    ChargeDensity(ChargeDensity&&) = delete;
    ChargeDensity & operator=(ChargeDensity const&) = delete;
    ChargeDensity & operator=(ChargeDensity &&) = delete;

    ChargeDensity(ChargeDensityConf::Config const& config, unsigned int view, unsigned int tpc);

    void Fill(std::vector<art::Ptr<recob::Hit>> hitsVect, VertexView vertex);
    void HelloWorld();
    void Save2ROOT(art::TFileDirectory tfdir, std::string name);

    double Delta(){return fDelta;}
    double Eta(){return fEta;}
    double FitScore(){return fFitScore;}
    double Alpha(){return fAlpha;}
    double Omega(){return fOmega;}
    double Tau(){return fTau;}
    double Iota(){return fIota;}
    double Score(){return fScore;}

  private:

    // configuration parameters
    bool fApplyRawSmoothing;
    bool fApplySmoothing;
    bool fApplyCumulativeSmoothing;
    unsigned int fNDriftPack;
    unsigned int fNWirePack;
    float fExpoAvSmoothPar;
    int fUnAvNeighbours;
    double fCumulativeCut;
    int fSlidingWindowN;
    int fMaxRadius;
    int fNSamplesBeginSlope;
    bool fDebugMode;
    bool fCalculateScore;
    std::string fTMVAFilename;
    //ChargeDensityConf::ConfParameters_t fPset;

    // view information
    unsigned int fView;
    unsigned int fTPC;

    // input vertex
    VertexView fVertex;

    // vectors storing the charge profile
    std::vector<double> fZ;
    std::vector<double> fRho;
    std::vector<double> fZCum;
    std::vector<double> fZCumStart;
    std::vector<double> fZCumDer;
    std::vector<double> fZCumDerErr;

    double fNormUnAvSmooth;
    double fZetaNorm;
    double fMaxCumulative;

    // analysis parameters
    float fDelta;
    float fEta;
    float fFitScore;
    float fAlpha;
    float fOmega;
    float fTau;
    float fIota;
    float fScore;

    float fGap;
    float fProtonKE;
    float fPionKE;

    TMVA::Reader fTMVAReader;

    double GetDistance(int x, double y, int x0, double y0);
    void UpdateMetrics();
    void ApplyExpoAvSmoothing(std::vector<double>& wf);
    void ApplyUnAvSmoothing(std::vector<double>& wf);
    void SlidingWindow(std::vector<double>& wf);
    void FillCumulative();
    void Rescale(std::vector<double>& wf, double kappa);

};


#endif
