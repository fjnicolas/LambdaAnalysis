////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesDisplay
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_DISPLAY_H
#define TPC_LINES_DISPLAY_H

#include <vector>

// ROOT includes
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TLegend.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
//#include "THistPainter.h"
#include <TROOT.h>

#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleTriangles.h"
#include "TPCSimpleEvents.h"


class TPCLinesDisplay {
    private:
        void DrawHitScatter(std::vector<SHit> hitV, TLegend& leg, std::string label, int color, int style, double size, double errorAlpha);
        void DrawLinearCluster(SLinearCluster cluster, TLegend& leg, std::string label, int color, double size=1.1, int style=4);
        void DrawLine(LineEquation line, double xmin, double xmax, TLegend& leg, std::string label, int color, int style);
        TH2F GetFrame(std::vector<SHit> hitsV, std::string label);
        void DrawTriangle(STriangle tri, TLegend& leg, std::string label, int colorP, int color, double alpha);
        void DrawVertex(SVertex vertex, TLegend& leg, std::string label, int color, int marker, double alpha);
        void SetStyle();

        std::vector<int> fColors;
        std::vector<int> fColorsOrigins;
        double fLegendFontSize=0.12;
        std::string fOutputPath;
    public:
        TPCLinesDisplay(bool show, std::string outputPath="Plots/");

        void Show(
            std::string eventLabel,
            std::vector<SHit> allHitsV,
            LineEquation houghLine,
            std::vector<SHit> selectedHitsV={},
            std::vector<SLinearCluster> clustersV={},
            SLinearCluster mainDirection = SLinearCluster(std::vector<SHit> {}),
            std::vector<STriangle> originAngles = {}, 
            SVertex recoVertex = SVertex(),
            SVertex trueVertex = SVertex(),
            std::vector<SOrigin> origins = {} );
};

#endif // TPC_SIMPLE_LINES_H