////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleCalo.h"

#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"

// (1 / 0.0201293) e-/ADC*time_ticks x  23.6e-6 MeV/e-  x)

// Constructor to initialize the collection
SCalo::SCalo(const std::vector<SHit>& hits)
    : fHitIntegralToEnergy( (1 / 0.0201293) *  23.6e-6 ),
    fStepXToLength(0.3),
    fStepYToLength(0.075),
    fTrackLength(0),
    fHitList(hits)
{
    CalculateResidualRange();
}


// Function to calculate the distance between two hit points
double SCalo::CalculateDistance(const SHit& p1, const SHit& p2) {
    double dx = fStepXToLength * ( p2.X() - p1.X() );
    double dy = fStepYToLength * ( p2.Y() - p1.Y() );
    return std::sqrt(dx * dx + dy * dy);
}

// Method to calculate and store the path lengths
void SCalo::CalculatePathLengths() {
    fPathLengths.clear();

    // first hit
    fPathLengths.push_back( CalculateDistance(fHitList[0], fHitList[1]) );

    for (size_t i = 1; i < fHitList.size()-1; ++i) {
        double distance1 = CalculateDistance(fHitList[i - 1], fHitList[i]);
        double distance2 = CalculateDistance(fHitList[i], fHitList[i + 1]);
        double halfDistance = 0.5 * (distance1 + distance2);
        
        fPathLengths.push_back(halfDistance);
    }

    // last hit
    fPathLengths.push_back( CalculateDistance(fHitList[fHitList.size()-1], fHitList[fHitList.size()-2]) );

}

void SCalo::CalculateResidualRange() {
    
    CalculatePathLengths();
    
    fDepositedEnergy.clear();
    fResidualRange.clear();

    fTrackLength = 0;
    fTrackLength = std::accumulate(fPathLengths.begin(), fPathLengths.end(), fTrackLength);

    double pathLength = 0;
    for (size_t i = 0; i < fHitList.size(); ++i) {
        fResidualRange.push_back( fTrackLength - pathLength);
        fDepositedEnergy.push_back(fHitList[i].Integral() * fHitIntegralToEnergy / fPathLengths[i]);
        pathLength+=fPathLengths[i];
    }

}


void SCalo::Display(){
    std::cout<<"Displaying calorimetry...track length: "<<fTrackLength<<"\n";

    for (size_t i = 0; i < fPathLengths.size(); ++i) {
        std::cout << " Residual range / dEdx: " << fResidualRange[i] << " / " << fDepositedEnergy[i] << std::endl;
    }

}

void CreateEnergyLossVsResidualRangePlot(const std::vector<SCalo>& caloObjects) {
    TCanvas* canvas = new TCanvas("energy_loss_vs_residual_range", "Energy Loss vs. Residual Range", 800, 600);
    canvas->Divide(1,1);
    std::vector<int> fColorsOrigins = {40, 42, 46, 30, 35, 40, 42, 46, 30, 35};

    canvas->cd(1);

    const std::vector<double> residualRange1 = caloObjects[0].GetResidualRange();
    const std::vector<double> depositedEnergy1 = caloObjects[0].GetDepositedEnergy();
    TGraph* graph1 = new TGraph(residualRange1.size(), &residualRange1[0], &depositedEnergy1[0]);
    graph1->SetTitle("Track 1");
    graph1->SetLineColor(kRed+2);
    graph1->SetMarkerColor(kRed+2);
    graph1->SetMarkerStyle(20);

    const std::vector<double> residualRange2 = caloObjects[1].GetResidualRange();
    const std::vector<double> depositedEnergy2 = caloObjects[1].GetDepositedEnergy();
    TGraph* graph2 = new TGraph(residualRange2.size(), &residualRange2[0], &depositedEnergy2[0]);
    graph2->SetTitle("Track 2");
    graph2->SetLineColor(kBlue+2);
    graph2->SetMarkerColor(kBlue+2);
    graph2->SetMarkerStyle(20);

    double maxX = std::max(caloObjects[0].GetTrackLength(), caloObjects[1].GetTrackLength());
    TH2F hFrame("hFrame", ";Residual range [cm];dEdx [MeV/cm]", 200, 0, maxX, 100, 0, 20);
    hFrame.Draw();
    graph1->Draw("lp same");
    graph2->Draw("lp same");

   
    canvas->cd();

    canvas->Update();
    canvas->SaveAs("energy_loss_vs_residual_range.pdf");
}