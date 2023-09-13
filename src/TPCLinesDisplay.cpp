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
#include "THistPainter.h"
#include <TROOT.h>


#include "SObjects/TPCSimpleHits.h"
#include "SObjects/TPCSimpleClusters.h"
#include "SObjects/TPCSimpleTriangles.h"
#include "SObjects/TPCSimpleEvents.h"


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
        TPCLinesDisplay(bool show, std::string outputPath="Plots/"):
            fOutputPath(outputPath)
        {
           fColors = {kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5};

           fColorsOrigins = {40, 42, 46, 30, 35, 40, 42, 46, 30, 35};

           if(show==false){
                gROOT->SetBatch(true);
           }

           // clean the directory
            gSystem->Exec(("rm -rf "+fOutputPath).c_str());
           
        }

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

void TPCLinesDisplay::SetStyle(){
    //TITLES SIZE AND FONT
    gStyle->SetPalette(112,0);
    gStyle->SetTitleFont(132, "TXYZ");
    gStyle->SetTitleSize(0.05, "TXYZ");

    gStyle->SetTitleFont(132, "titleFont"); 
    
    // Off stats
    gStyle->SetOptStat(0);     // Turn off statistics box

    
    //LABELS SIZE AND FONT
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    gStyle->SetTitleYOffset (1.4);
    
    //AXIS OFFSETS AND SIZES
    /*gStyle->SetTitleXOffset (1.);
    gStyle->SetTitleXSize (0.05);
    gStyle->SetTitleYOffset (1.);
    gStyle->SetTitleYSize (0.05);*/


}

void TPCLinesDisplay::DrawVertex(SVertex vertex, TLegend& leg, std::string label, int color, int marker, double alpha){
    TGraph *pointGraph = new TGraph();
    pointGraph->SetPoint(0, vertex.Point().X(), vertex.Point().Y()); // Set the point's coordinates
    // Set marker style and size for the point
    pointGraph->SetMarkerStyle(marker);
    pointGraph->SetMarkerSize(2.);
    pointGraph->SetMarkerColorAlpha(color, alpha);

    if(label!="")
        leg.AddEntry(pointGraph, label.c_str(), "p");

    // Draw the point
    pointGraph->Draw("P");
}

TH2F TPCLinesDisplay::GetFrame(std::vector<SHit> hitsV, std::string label){

    std::vector<double> x, y, err;
    for(size_t ix=0; ix<hitsV.size(); ix++){
        x.push_back(hitsV[ix].X());
        y.push_back(hitsV[ix].Y());
        err.push_back(hitsV[ix].Width());
    }
    double minX = *std::min_element(x.begin(), x.end());
    double maxX = *std::max_element(x.begin(), x.end());
    double minY = *std::min_element(y.begin(), y.end());
    double maxY = *std::max_element(y.begin(), y.end());
    double overX = 0.1*(maxX-minX);
    double overY = 0.1*(maxY-minY);

    TH2F hFrame("frame", (label+";WireId;TimeTick [0.5 #mu s]").c_str(), 200, minX-overX, maxX+overX, 200, minY-overY, maxY+overY);
    hFrame.SetStats(0);
    return hFrame;
}

void TPCLinesDisplay::DrawLine(LineEquation line, double xmin, double xmax, TLegend& leg, std::string label, int color, int style){
    // Draw a horizontal line from x = 1 to x = 4 at y = 2
    double y1 = line.EvaluateX(xmin);
    double y2 = line.EvaluateX(xmax);
    TLine* horizontalLine = new TLine(xmin, y1, xmax, y2);
    horizontalLine->SetLineColor(color); // Set line color
    horizontalLine->SetLineStyle(style); // Set line color
    horizontalLine->SetLineWidth(2);     // Set line width
    horizontalLine->Draw();

    if(label!="")
        leg.AddEntry(horizontalLine, label.c_str(), "l");
}

void TPCLinesDisplay::DrawHitScatter(std::vector<SHit> hitsV, TLegend& leg, std::string label, int color, int style, double size, double errorAlpha){

    if(hitsV.size()>0){
        std::vector<double> x, y, err;
        for(size_t ix=0; ix<hitsV.size(); ix++){
            x.push_back(hitsV[ix].X());
            y.push_back(hitsV[ix].Y());
            err.push_back(hitsV[ix].Width());
        }

        TGraphErrors *g = new TGraphErrors(x.size(),&x[0],&y[0], nullptr, &err[0]); 
        g->SetMarkerColorAlpha(color, 0.5);
        g->SetMarkerStyle(style);
        g->SetMarkerSize(size);
        g->SetLineColorAlpha(kGray, errorAlpha);
        g->Draw("p");
        
        if(label!="")
            leg.AddEntry(g, label.c_str(), "p");
    }

    return;
}

void TPCLinesDisplay::DrawTriangle(STriangle tri, TLegend& leg, std::string label, int colorP, int color, double alpha){ 
    // Define the triangle's vertices
    Double_t x[3] = {tri.GetMainVertex().X(), tri.GetVertexB().X(), tri.GetVertexC().X()};
    Double_t y[3] = {tri.GetMainVertex().Y(), tri.GetVertexB().Y(), tri.GetVertexC().Y()};

    // Create a polyline representing the triangle
    TPolyLine *triangle = new TPolyLine(3, x, y, "F");
    triangle->SetFillColorAlpha(color, alpha);  // Set fill color

    TGraph *pointGraph = new TGraph();
    pointGraph->SetPoint(0, tri.GetMainVertex().X(), tri.GetMainVertex().Y()); // Set the point's coordinates
    // Set marker style and size for the point
    pointGraph->SetMarkerStyle(29);
    pointGraph->SetMarkerSize(2.2);
    pointGraph->SetMarkerColor(colorP);

    if(label!="")
        leg.AddEntry(pointGraph, label.c_str(), "p");


    // Draw the point
    pointGraph->Draw("P"); // A: Axis, P: Po

    // Draw the triangle
    triangle->Draw("F");
}

void TPCLinesDisplay::DrawLinearCluster(SLinearCluster cluster, TLegend& leg, std::string label, int color, double size, int style){
    if(cluster.NHits()>0){
        std::vector<SHit> hits = cluster.GetHits();
        std::vector<double> x, y;
        for(size_t ix=0; ix<hits.size(); ix++){
            x.push_back(hits[ix].X());
            y.push_back(hits[ix].Y());
        }

        TGraph *g = new TGraph(x.size(),&x[0],&y[0]); 
        g->SetMarkerColorAlpha(color, 0.5);
        g->SetMarkerStyle(style);;
        g->SetMarkerSize(size);
        g->SetLineWidth(20);
        //g->SetLineColorAlpha(kGray, errorAlpha);
        g->Draw("p");

        DrawLine(cluster.GetTrackEquation(), cluster.GetMinX(), cluster.GetMaxX(), leg, "", color, kSolid);
        DrawLine(cluster.GetTrackEquationStart(), cluster.GetMinX(), cluster.GetMaxX(), leg, "", color, kDashed);
        DrawLine(cluster.GetTrackEquationEnd(), cluster.GetMinX(), cluster.GetMaxX(), leg, "", color, kDashed);

        leg.AddEntry(g, ( "Cluster  "+std::to_string(cluster.GetId())).c_str(), "p");
    }

    
    return;
}



void TPCLinesDisplay::Show(
    std::string eventLabel,
    std::vector<SHit> allHitsV,
    LineEquation houghLine,
    std::vector<SHit> selectedHitsV,
    std::vector<SLinearCluster> clustersV,
    SLinearCluster mainDirection,
    std::vector<STriangle> originAngles,
    SVertex recoVertex,
    SVertex trueVertex,
    std::vector<SOrigin> origins)
{
    
    SetStyle();

    

    if(allHitsV.size()==0){return;}
    TCanvas c("c", eventLabel.c_str(), 0, 0, 1000, 800);

    TPad *pad1 = new TPad("pad1", "Graph Pad", 0.0, 0., .8, 1.0);
    TPad *pad2 = new TPad("pad2", "Legend Pad", 0.75, 0., 1.0, 1.);

    pad1->SetLeftMargin(0.15);
    pad2->SetLeftMargin(-0.2);
    pad1->Draw();
    pad2->Draw();
    
    TLegend legend(0., 0.1, 0.99, 0.9); // (x1, y1, x2, y2)

    pad1->cd();

    // general frame
    TH2F hFrame=GetFrame(allHitsV, eventLabel);
    hFrame.Draw();

    // Main direction
    if(mainDirection.NHits()>0){
        // selected hits scatter
        std::vector<SHit> mainDirHits = mainDirection.GetHits();
        DrawHitScatter(mainDirHits, legend, "MainLine", kBlack, 8, 1.8, 0);   
    }

    // all hits scatter
    DrawHitScatter(allHitsV, legend, "AllHits", 65, 8, 1, 0.6);

    if(selectedHitsV.size()!=0){
        // selected hits scatter
        DrawHitScatter(selectedHitsV, legend, "HoughHits", kRed, 24, 1.1, 0);

        // hough line
        DrawLine(houghLine, hFrame.GetXaxis()->GetXmin(), hFrame.GetXaxis()->GetXmax(), legend, "Hough line", kBlack, kDashed);
    }
   
    // linear clusters
    if(clustersV.size()>0){

        for(size_t cIx=0; cIx<clustersV.size(); cIx++){
            DrawLinearCluster(clustersV[cIx], legend, "Cluster", fColors[cIx], 1., 8);
        }
    }

    // triangles
    for(size_t oIx=0; oIx<originAngles.size(); oIx++){
        std::cout<<"Drawing origin angle "<<oIx<<std::endl;
        DrawTriangle(originAngles[oIx], legend, "Origin "+std::to_string(oIx), fColorsOrigins[oIx], 90, 0.5);
    }


    // triangles
    for(size_t oIx=0; oIx<origins.size(); oIx++){
        std::cout<<"Drawing origin "<<oIx<<std::endl;
        DrawVertex(origins[oIx].GetPoint(), legend, "Origin "+std::to_string(oIx), fColorsOrigins[oIx], 90, 0.5);
    }


    // draw PANDORA vertex
    if(recoVertex.Point().X()!=-1 && recoVertex.Point().Y()!=-1){
        DrawVertex(recoVertex, legend, "PANDORA vtx", kPink-9, 33, 1);
    }

    // draw PANDORA vertex
    if(trueVertex.Point().X()!=-1 && trueVertex.Point().Y()!=-1){
        DrawVertex(trueVertex, legend, "True vtx", kPink-12, 33, 1);
    }
    
    
    // Create a legend
    pad2->cd();
    legend.SetBorderSize(0);
    legend.SetTextSize(fLegendFontSize); 
    legend.Draw();


    // Check if the directory exists, create it if not
    gSystem->Exec(("pwd "+fOutputPath).c_str());
    if (!gSystem->OpenDirectory(fOutputPath.c_str())) {
        gSystem->Exec(("mkdir "+fOutputPath).c_str());
        gSystem->Exec(("mkdir "+fOutputPath+"/rootfiles").c_str());
    }  
    TFile* rootFile = new TFile((fOutputPath + "/rootfiles/"+eventLabel).c_str(), "RECREATE");
    c.Write();
    rootFile->Close();
    c.SaveAs((fOutputPath + "/" + eventLabel+".pdf").c_str());

    c.cd();

    
    c.Update();
    c.WaitPrimitive(); 





    return;

}

#endif // TPC_SIMPLE_LINES_H