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

#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "Fit/Fitter.h"
#include <Math/Functor.h>
#include "TLine.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TPaveText.h"

// Calorimetry conversion factor (ADC to charge)
// (1 / 0.0201293) e-/ADC*time_ticks x  23.6e-6 MeV/e-  x)



// Constructor to initialize the collection
SCalo::SCalo(const std::vector<SHit>& hits, double trackAngle)
    : fHitIntegralToEnergy( (1 / 0.0201293) *  23.6e-6 ),
    fStepXToLength(0.3),
    fStepYToLength(0.075),
    fTrackLength(0),
    fCosTrackAngle( std::abs(std::cos(trackAngle)) ),
    fTrackAngle(trackAngle),
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
        fDepositedEnergy.push_back( (fHitList[i].Integral() * fHitIntegralToEnergy) * fCosTrackAngle ); // / fPathLengths[i]);
        pathLength+=fPathLengths[i];
    }

}


void SCalo::Display(){
    std::cout<<"Displaying calorimetry...track length: "<<fTrackLength<<"\n";

    for (size_t i = 0; i < fPathLengths.size(); ++i) {
        std::cout << " Residual range / dEdx: " << fResidualRange[i] << " / " << fDepositedEnergy[i] << std::endl;
    }

}



// ------ STriangleCalo class ------

// --- Constructor ---
STriangleCalo::STriangleCalo(STriangle triangle):
    fTriangle(triangle),
    fCalo1(triangle.GetTrack1().GetHits()),
    fCalo2(triangle.GetTrack2().GetHits())
{
}


// --- Function to display dEdx ---
void STriangleCalo::MakeEnergyLossVsResidualRangePlot(SCalo fCalo1, SCalo fCalo2, TPad *pad) {
    pad->cd();

    // Get vectors
    const std::vector<double> residualRange1 = fCalo1.GetResidualRange();
    const std::vector<double> depositedEnergy1 = fCalo1.GetDepositedEnergy();
    const std::vector<double> residualRange2 = fCalo2.GetResidualRange();
    const std::vector<double> depositedEnergy2 = fCalo2.GetDepositedEnergy();

    // Define graph 1
    TGraph* graph1 = new TGraph(residualRange1.size(), &residualRange1[0], &depositedEnergy1[0]);
    graph1->SetTitle("Track 1");
    graph1->SetLineColor(fColor1);
    graph1->SetMarkerColor(fColor1);
    graph1->SetMarkerStyle(20);

    // Define graph 2
    TGraph* graph2 = new TGraph(residualRange2.size(), &residualRange2[0], &depositedEnergy2[0]);
    graph2->SetTitle("Track 2");
    graph2->SetLineColor(fColor2);
    graph2->SetMarkerColor(fColor2);
    graph2->SetMarkerStyle(20);

    // Frame
    double maxX = std::max(fCalo1.GetTrackLength(), fCalo2.GetTrackLength());
    double maxY = std::max(*std::max_element(depositedEnergy1.begin(), depositedEnergy1.end()), *std::max_element(depositedEnergy2.begin(), depositedEnergy2.end()));
    TH2F *hFrame = new TH2F("hFrame", ";Residual range [cm];Charge [AU]", 200, 0, maxX, 100, 0, maxY);
    hFrame->GetYaxis()->SetTitleOffset(.8);
    hFrame->GetXaxis()->SetTitleOffset(.8);
    hFrame->SetStats(0);
    hFrame->Draw();
    graph1->Draw("lp same");
    graph2->Draw("lp same");

    // TLegend
    TLegend *leg = new TLegend(0.7, 0.5, 0.85, 0.65);
    leg->AddEntry(graph1, "Track1", "lp");
    leg->AddEntry(graph2, "Track2", "lp");
    leg->Draw("same");

    // Fit the two graphs to two constants
    TF1 *fit1 = new TF1("fit1", "[0]", 0, maxX);
    fit1->SetLineColor(fColor1);
    fit1->SetLineWidth(2);
    graph1->Fit(fit1, "R");
    TF1 *fit2 = new TF1("fit2", "[0]", 0, maxX);
    fit2->SetLineColor(fColor2);
    fit2->SetLineWidth(2);
    graph2->Fit(fit2, "R");

    std::cout<<" dEdX fit done!\n";

    // Draw the fitted functions
    fit1->Draw("same");
    fit2->Draw("same");

    double charge1Fit = std::min(fit1->GetParameter(0), fit2->GetParameter(0));
    double charge2Fit = std::max(fit1->GetParameter(0), fit2->GetParameter(0));

    fChargeRatioFit = charge1Fit / charge2Fit;
    fChargeDifferenceFit = charge2Fit - charge1Fit;
    fChargeRelativeDifferenceFit = fChargeDifferenceFit / charge1Fit;

    double charge1 = std::accumulate(depositedEnergy1.begin(), depositedEnergy1.end(), 0.0) / depositedEnergy1.size();
    double charge2 = std::accumulate(depositedEnergy2.begin(), depositedEnergy2.end(), 0.0) / depositedEnergy2.size();

    double charge1Average = std::min(charge1, charge2);
    double charge2Average = std::max(charge1, charge2);

    fChargeRatioAverage = charge1Average / charge2Average;
    fChargeDifferenceAverage = charge2Average - charge1Average;
    fChargeRelativeDifferenceAverage = fChargeDifferenceAverage / charge1Average;

}


void STriangleCalo::CreateEnergyLossVsResidualRangePlot() {
    TCanvas *c1 = new TCanvas("canvasCalo","Calorimetry", 600,1200);
    c1->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.02,0.98,0.98);
    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, pad1);
    delete c1;
    delete pad1;
}

// --- Broken line fit ---
Double_t BrokenLine(Double_t *x,Double_t *par) {
  
  Double_t arg = 0.;
  if(x[0]<par[0])
    arg = par[1]+par[2]*(x[0]-par[0]);
  else
    arg = par[1]+par[3]*(x[0]-par[0]);
  
  return arg; 
}


// --- Fisher discriminant 2D ---
TVectorD FisherDiscriminant2D(const std::vector<SHit>& sample1, const std::vector<SHit>& sample2) {
    // Calculate means for each sample
    double mean1X = 0.0, mean1Y = 0.0;
    double mean2X = 0.0, mean2Y = 0.0;

    for (const auto& point : sample1) {
        mean1X += point.X();
        mean1Y += point.Y();
    }

    for (const auto& point : sample2) {
        mean2X += point.X();
        mean2Y += point.Y();
    }

    mean1X /= sample1.size();
    mean1Y /= sample1.size();

    mean2X /= sample2.size();
    mean2Y /= sample2.size();

    std::cout<<mean1X<<" "<<mean1Y<<" "<<mean2X<<" "<<mean2Y<<std::endl;

    // Calculate covariance matrices for each sample
    TMatrixD cov1(2, 2), cov2(2, 2);

    for (const auto& point : sample1) {
        cov1(0, 0) += (point.X() - mean1X) * (point.X() - mean1X);
        cov1(0, 1) += (point.X() - mean1X) * (point.Y() - mean1Y);
        cov1(1, 0) += (point.Y() - mean1Y) * (point.X() - mean1X);
        cov1(1, 1) += (point.Y() - mean1Y) * (point.Y() - mean1Y);
    }

    for (const auto& point : sample2) {
        cov2(0, 0) += (point.X() - mean2X) * (point.X() - mean2X);
        cov2(0, 1) += (point.X() - mean2X) * (point.Y() - mean2Y);
        cov2(1, 0) += (point.Y() - mean2Y) * (point.X() - mean2X);
        cov2(1, 1) += (point.Y() - mean2Y) * (point.Y() - mean2Y);
    }

    cov1 *= 1.0 / (sample1.size() - 1);
    cov2 *= 1.0 / (sample2.size() - 1);

    std::cout<<" Cov "<<cov1(0,0)<<" "<<cov1(0,1)<<" "<<cov1(1,0)<<" "<<cov1(1,1)<<std::endl;

    TMatrixD cov1Inv (cov1.Invert());

    std::cout<<" CovInv "<<cov1Inv(0,0)<<" "<<cov1Inv(0,1)<<" "<<cov1Inv(1,0)<<" "<<cov1Inv(1,1)<<std::endl;
    TVectorD menDiffVector(0, 1, mean1X - mean2X, mean1Y - mean2Y, "END");
    std::cout<<" menDiffVector "<<menDiffVector(0)<<" "<<menDiffVector(1)<<std::endl;
    
    // Calculate the Fisher's discriminant direction
    TVectorD w = cov1Inv * menDiffVector;

    // Print the discriminant direction
    std::cout << "Fisher's Discriminant Direction: (" << w[0] << ", " << w[1] << ")" << std::endl;

    return w;
}


// --- JointFit function ---
void STriangleCalo::JointFitAnalysisFisher(){

    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example", 800, 0, 800,1200);
    c1->SetGrid();
    // Create TPads
    c1->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.98,0.48);
    pad1->Draw();
    pad2->Draw();

    size_t nhits1 = fTriangle.GetNHitsTrack1();
    size_t nhits2 = fTriangle.GetNHitsTrack2();
    size_t maxNHits = std::min(nhits1, nhits2);
    maxNHits = std::min(nhits1, nhits2);
    //maxNHits = 20;

    // x-y values from hits
    size_t hitCont = 0;
    std::vector<double> xV1, yV1, xV2, yV2;
    for(SHit &h:fTriangle.GetTrack1().GetHits()){
        if(hitCont>maxNHits) continue;
        xV1.push_back(h.X());
        yV1.push_back(h.Y());
        hitCont++;
        
    }
    hitCont=0;
    for(SHit &h:fTriangle.GetTrack2().GetHits()){
        if(hitCont>maxNHits) continue;
        xV2.push_back(h.X());
        yV2.push_back(h.Y());
        hitCont++;
    }

    double xVertex = fTriangle.GetMainVertex().X();
    double yVertex = fTriangle.GetMainVertex().Y();

    // Fisher direcrtion
    TVectorD fWVector = FisherDiscriminant2D(fTriangle.GetTrack1().GetHits(), fTriangle.GetTrack2().GetHits());

    // opposite side
    double x1 = fTriangle.GetVertexB().X();
    double y1 = fTriangle.GetVertexB().Y();
    double x2 = fTriangle.GetVertexC().X();
    double y2 = fTriangle.GetVertexC().Y();

    double maxD = 0;
    for(size_t k=0; k<xV1.size(); k++){
        double d = std::hypot( xV1[k]-xVertex, yV1[k]-yVertex);
        if(d>maxD){
            x1 = xV1[k];
            y1 = yV1[k];
            maxD = d;
        }
    }

    maxD = 0;
    for(size_t k=0; k<xV2.size(); k++){
        double d = std::hypot( xV2[k]-xVertex, yV2[k]-yVertex);
        if(d>maxD){
            x2 = xV2[k];
            y2 = yV2[k];
            maxD = d;
        }
    }

    std::cout<<" CHEEEEECK \n";
    std::cout<<" Opposited coordinates: "<<x1<<":"<<y1<<"   "<<x2<<":"<<y2<<std::endl;

    double m = (y2-y1)/(x2-x1);

    double rawSlope1 = fTriangle.GetTrack1().GetTrackEquation().Slope();
    double rawSlope2 = fTriangle.GetTrack2().GetTrackEquation().Slope();

    std::vector<double> x, y;
    x.insert(x.end(), xV1.begin(), xV1.end());
    x.insert(x.end(), xV2.begin(), xV2.end());
    y.insert(y.end(), yV1.begin(), yV1.end());
    y.insert(y.end(), yV2.begin(), yV2.end());

    // project the points
    for(size_t i=0; i<x.size(); i++){
        x[i] = x[i]-xVertex;
        y[i] = y[i]-yVertex;
    }

    size_t n = x.size();


    // project into the line with points (17.3, 12.4) and (13.7, 6.65)
    // get the line equation
    // y = m*x + q
    pad1->cd();
    // create a graph
    TGraph *grRaw = new TGraph(n, &x[0], &y[0]);
    grRaw->SetMarkerColor(4);
    grRaw->SetMarkerStyle(21);
    grRaw->Draw("AP");
    // axis labels
    grRaw->GetXaxis()->SetTitle("Wire ID");
    grRaw->GetYaxis()->SetTitle("Time Tick [0.5#mu s]");
    // axis label offset
    grRaw->GetYaxis()->SetTitleOffset(.7);
    grRaw->GetXaxis()->SetTitleOffset(.83);
    grRaw->SetTitle("View (raw)");
    
    // angle
    double angle = atan(m);
    std::cout << "m: " << m << std::endl;
    std::cout << "angle: " << 180 * angle / TMath::Pi() << std::endl;
    std::cout << "opening angle: " << fTriangle.GetOpeningAngle() << std::endl;
    double isocAngle = ( 180 - fTriangle.GetOpeningAngle() ) / 2;
    std::cout << "isocAngle: " << isocAngle << std::endl;

    // convert the ray slopes
    rawSlope1 = sin(-angle)+rawSlope1*cos(-angle) / (cos(-angle)-rawSlope1*sin(-angle));
    rawSlope2 = sin(-angle)+rawSlope2*cos(-angle) / (cos(-angle)-rawSlope2*sin(-angle));
   
    // distances to the main vertex
    double d1 = std::hypot( x1, y1);
    double d2 = std::hypot( x2, y2);

    // min track angle
    double minTrackAngle1 = std::acos( x1/std::hypot(x1, y1) );
    double minTrackAngle2 = std::acos( x2/std::hypot(x2, y2) );
    std::cout<<" Side track angle: "<<180*minTrackAngle1/TMath::Pi()<<" "<<180*minTrackAngle2/TMath::Pi()<<std::endl;

    // Set angle to the isoceles triangle
    double minTrackAngle = (d1>d2)? minTrackAngle1:minTrackAngle2;
    double newSlope =  180 - 180*minTrackAngle/TMath::Pi() - isocAngle;
    newSlope = TMath::Pi() * newSlope / 180;
    angle = newSlope;
    std::cout<<" New slope: "<< newSlope << "angle "<<angle<<std::endl;

    // Set angle to the Fisher Direction
    angle  = fWVector(1)/fWVector(0);

    // project the points
    for(size_t i=0; i<x.size(); i++){
        double x0 = x[i];
        double y0 = y[i];
        double xproj = x0*cos(angle) + y0*sin(angle);
        double yproj = -x0*sin(angle) + y0*cos(angle);
        x[i] = xproj;
        y[i] = yproj;
    }
    // get min and max x values
    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    xmin -= 0.1 * (xmax-xmin);
    xmax += 0.1 * (xmax-xmin);
    std::cout << "xmin: " << xmin << std::endl;
    std::cout << "xmax: " << xmax << std::endl;


    // create a graph
    pad2->cd();
    
    TGraph *gr = new TGraph(n, &x[0], &y[0]);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetTitle("Wire ID");
    gr->GetYaxis()->SetTitle("Time Tick [0.5#mu s]");
    // axis label offset
    gr->GetYaxis()->SetTitleOffset(.7);
    gr->GetXaxis()->SetTitleOffset(.83);

    gr->SetTitle("Rotated view");
    gr->Draw("AP");
    c1->Update();


    // Fit to the brokenLine function
    TF1 *brokenLine = new TF1("brokenLine",BrokenLine, xmin, xmax, (int)(xmax-xmin));
    brokenLine->SetParameters(0, 0, rawSlope1, rawSlope2);
    brokenLine->SetParNames("breakPoint","y0","slope1","slope2");
    brokenLine->SetLineColor(kRed);
    brokenLine->SetLineWidth(2);
    gr->Fit(brokenLine,"R");

    // Draw the fitted function
    brokenLine->Draw("same");

    pad1->cd();
    // Undo the transformation
    double breakPoint = brokenLine->GetParameter(0);
    double y0 = brokenLine->GetParameter(1);
    double slope1 = brokenLine->GetParameter(2);
    double slope2 = brokenLine->GetParameter(3);
    
    // Convert the parameters
    breakPoint = breakPoint*cos(angle) - y0*sin(angle);
    y0 = breakPoint*sin(angle) + y0*cos(angle);
    
    slope1 = ( sin(angle)+slope1*cos(angle) ) / (cos(angle)-slope1*sin(angle));

    slope2 = ( sin(angle)+slope2*cos(angle) )/ (cos(angle)-slope2*sin(angle));

    std::cout<<" Break point: "<<breakPoint<<" y0: "<<y0<<" slope1: "<<slope1<<" slope2: "<<slope2<<std::endl;


    // Draw the fitted function
    xmin = *std::min_element(x.begin(), x.end());
    xmax = *std::max_element(x.begin(), x.end());
    TF1 *line1 = new TF1("line1","[0]+[1]*x", xmin, xmax);
    double q1 = y0 - slope1*breakPoint;
    line1->SetParameters(q1, slope1);
    line1->SetLineColor(kOrange+7);
    line1->SetLineWidth(2);
    line1->Draw("same");
    TF1 *line2 = new TF1("line2","[0]+[1]*x", xmin, xmax);
    double q2 = y0 - slope2*breakPoint;
    line2->SetParameters(q2, slope2);
    line2->SetLineColor(kGreen+2);
    line2->SetLineWidth(2);
    line2->Draw("same");

    // Draw m/q line
    TF1 *line3 = new TF1("line3","[0]+[1]*x", xmin, xmax);
    double q = y1-yVertex - m*(x1-xVertex);
    line3->SetParameters(q, m);
    line3->SetLineColor(kPink+2);
    line3->SetLineWidth(2);
    line3->Draw("same");

    
    c1->Update();

    // Print the results
    std::cout << "Fit results:" << std::endl;
    std::cout << "breakPoint: " << brokenLine->GetParameter(0) << std::endl;


}


// --- Sort hits by distance ---
void SortSHitVectorByDistance(std::vector<SHit>& hits, double externalX, double externalY) {
    // Lambda expression for comparing SHit objects based on distance
    auto distanceComparator = [externalX, externalY](const SHit& a, const SHit& b) {
        // Function to calculate the distance between two points
        auto calculateDistance = [](const SHit& point, double externalX, double externalY) {
            double deltaX = point.X() - externalX;
            double deltaY = point.Y() - externalY;
            return std::sqrt(deltaX * deltaX + deltaY * deltaY);
        };

        double distanceA = calculateDistance(a, externalX, externalY);
        double distanceB = calculateDistance(b, externalX, externalY);
        return distanceA < distanceB;
    };

    // Sort the vector using the distanceComparator lambda
    std::sort(hits.begin(), hits.end(), distanceComparator);
}


// --- JointFit analysis---
void STriangleCalo::JointFitAnalysis(unsigned int maxHits, double widthTol, bool useHitError, double& fitSlope1, double& fitSlope2, ChargeDensity & chargeDensityAlgo){

    // --- Create the canvas ---
    TCanvas *c1 = new TCanvas("canvasCalo","Calorimetry", 600,1200);
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.98,0.48);
    pad1->Draw();
    pad2->Draw();

    // --- Get the hits and the vertex ---
    // Vertex
    double xVertex = fTriangle.GetMainVertex().X();
    double yVertex = fTriangle.GetMainVertex().Y();

    // All hits
    std::vector<SHit> hitList = fTriangle.GetAllHits();

    // Hit coordinates fior the first maxHits hits
    std::vector<double> xV, yV, widthV;
    size_t n = 0;
    for(SHit &h:fTriangle.GetTrack1().GetHits()){
        if(n>=maxHits) break;
        xV.push_back(h.X()-xVertex);
        yV.push_back(h.Y()-yVertex);
        widthV.push_back(h.Width());
        n++;
    }
    n=0;
    for(SHit &h:fTriangle.GetTrack2().GetHits()){
        if(n>=maxHits) break;
        xV.push_back(h.X()-xVertex);
        yV.push_back(h.Y()-yVertex);
        widthV.push_back(h.Width());
        n++;
    }

    // --- Create the TGraph ---
    n = xV.size();
    TGraphErrors* gr = new TGraphErrors(n, &xV[0], &yV[0], 0, &widthV[0]);
    gr->SetMarkerStyle(3);

    // --- Chi2 minimizer function for two lines with common vertex ---
    auto Chi2TwoLines = [&](const double *par) {
        // minimisation function computing the sum of squares of residuals
        int np = gr->GetN();
        double f = 0;
        double *x = gr->GetX();
        double *y = gr->GetY();
        double *error = gr->GetEY();

        // loop over the data points
        for (int i=0;i<np;i++) {
            // get values
            double xH = x[i];
            double yH = y[i];
            double errorH =  useHitError? error[i]:1;

            // get residuals to the two lines
            double dr1 = ( yH - ( par[0]*xH + par[1] ) ) / errorH ;
            double dr2 = ( yH - ( par[2]*xH + (par[0]-par[2])*par[3] + par[1] ) ) / errorH;
            
            // add contribution to the chi2 (minimum residual to the two lines)
            double dr = std::min(dr1*dr1,dr2*dr2);
            f += dr*dr;
        }

        return f;
    };


    // --- Fit the graph to the function ---
    ROOT::Math::Functor fcn(Chi2TwoLines, 4);
    ROOT::Fit::Fitter  fitter;
    // Initializations
    double pStart[4] = {fTriangle.GetTrack1().GetTrackEquation().Slope(),
                        0,
                        fTriangle.GetTrack2().GetTrackEquation().Slope(),
                        0};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("m1");
    fitter.Config().ParSettings(1).SetName("n1");
    fitter.Config().ParSettings(2).SetName("m2");
    fitter.Config().ParSettings(3).SetName("xB");

    // Do the fit
    bool ok = fitter.FitFCN();
    if (!ok) {
        fPassFit = false;
        Error("line3Dfit","Line3D Fit failed");
    }
    else{
        fPassFit = true;
    }

    // Fit results
    const ROOT::Fit::FitResult & result = fitter.Result();

    // --- Draw the graph and the fit ---
    pad1->cd();
    // Get frame values
    double x1 = *std::min_element(xV.begin(), xV.end());
    double x2 = *std::max_element(xV.begin(), xV.end());
    double yMin = *std::min_element(yV.begin(), yV.end());
    double yMax = *std::max_element(yV.begin(), yV.end());
    x1-=0.1*(x2-x1);
    x2+=0.1*(x2-x1);
    yMax+=0.1*(yMax-yMin);
    yMin-=0.1*(yMax-yMin);
    TH2F * hFrame = new TH2F("hFrame", ";Wire ID;Time Tick [0.5 #mu s]", 200, x1, x2, 200, yMin, yMax);
    hFrame->GetYaxis()->SetTitleOffset(.8);
    hFrame->GetXaxis()->SetTitleOffset(.8);
    hFrame->SetStats(0);


    // Get parameters
    double m1 = result.Parameter(0);
    double n1 = result.Parameter(1);
    double m2 = result.Parameter(2);
    double xB = result.Parameter(3);
    double n2 = (m1-m2)*xB+n1;

    // Line1
    double y1 = m1*x1+n1;
    double y2 = m1*x2+n1;
    TLine *line1 = new TLine(x1,y1,x2,y2);
    line1->SetLineColor(fColor1);
    line1->SetLineWidth(4); 

    // Line2
    y1 = m2*x1+n2;
    y2 = m2*x2+n2;
    TLine *line2 = new TLine(x1,y1,x2,y2);
    line2->SetLineColor(fColor2);
    line2->SetLineWidth(4);

    // Fitted vertex
    double xFitVtx = result.Parameter(3);
    double yFitVtx = m1*xFitVtx+n1;
    TMarker *marker = new TMarker(xFitVtx, yFitVtx, 29);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(2);

    // TLegend
    TLegend *leg = new TLegend(0.7, 0.5, 0.85, 0.65);
    leg->AddEntry(gr,"Data","p");
    leg->AddEntry(line1, "Fitted line 1", "l");
    leg->AddEntry(line2, "Fitted line 2", "l");
    leg->AddEntry(marker, "Fitted vertex", "p");

    // Create text with the fit results
    TPaveText *pt = new TPaveText(0.7, 0.7, 0.85, 0.85, "NDC");
    pt->AddText(Form("Fit result #chi^{2} = %.2f\n",result.MinFcnValue()));
    pt->AddText(Form("m1 = %.2f #pm %.2f",result.Parameter(0),result.ParError(0)));
    pt->AddText(Form("n1 = %.2f #pm %.2f",result.Parameter(1),result.ParError(1)));
    pt->AddText(Form("m2 = %.2f #pm %.2f",result.Parameter(2),result.ParError(2)));
    pt->AddText(Form("xB = %.2f #pm %.2f",result.Parameter(3),result.ParError(3)));



    // Draw
    hFrame->Draw();
    gr->Draw("p same");
    line1->Draw("same");
    line2->Draw("same");
    marker->Draw("same");
    leg->Draw("same");
    pt->Draw("same");

    // --- Assing each hit to a track ---
    std::vector<SHit> hitsTrack1, hitsTrack2;
    double maxIntegral = 0;
    for(SHit &h:hitList){

        // shifted values 
        double x = h.X()-xVertex;
        double y = h.Y()-yVertex;

        // values from fitted slopes
        double yH1 = m1*x+n1;
        double yH2 = m2*x+n2;

        // closest tracks
        bool isTrack1 = std::abs( y-yH1 )<widthTol*h.Width();
        bool isTrack2 = std::abs( y-yH2 )<widthTol*h.Width();

        // if shared hit, split the integral in two hits
        if(isTrack1 && isTrack2){
            SHit sharedHit(-1, h.X(), h.Y(), h.Width(), 0.5*h.Integral());
            hitsTrack1.push_back(sharedHit);
            hitsTrack2.push_back(sharedHit);
        }
        else if(isTrack1){
            hitsTrack1.push_back(h);
        }
        else if(isTrack2){
            hitsTrack2.push_back(h);
        }
        else{
            // if not width compatible, add to the closest
            double d1 = std::abs(y-yH1);
            double d2 = std::abs(y-yH2);
            d1<d2? hitsTrack1.push_back(h):hitsTrack2.push_back(h);
        }

        // update maximum integral
        if(h.Integral()>maxIntegral) maxIntegral = h.Integral();
    }

    // --- Draw xV and yV markers after assignment ---
    for(SHit &h:hitsTrack1){
        double size = (h.Integral()/maxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-xVertex, h.Y()-yVertex, 21);
        marker->SetMarkerColorAlpha(fColor1, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
    }

    for(SHit &h:hitsTrack2){
        double size = (h.Integral()/maxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-xVertex, h.Y()-yVertex, 20);
        marker->SetMarkerColorAlpha(fColor2, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
    }

    // --- Sort by distance to fitted vertex ---
    SortSHitVectorByDistance(hitsTrack1, xVertex+xFitVtx, yVertex+yFitVtx);
    SortSHitVectorByDistance(hitsTrack2, xVertex+xFitVtx, yVertex+yFitVtx);

    // --- Calorimetry objects ---
    SCalo calo1(hitsTrack1, std::atan(m1));
    SCalo calo2(hitsTrack2, std::atan(m2));
    MakeEnergyLossVsResidualRangePlot(calo1, calo2, pad2);
    
    
    c1->Update();
    c1->cd();

    // Return values
    fitSlope1 = m1;
    fitSlope2 = m2;
    fTwoLinesChi2 = result.MinFcnValue();

    // vertex charge vs track charge
    double vertexHitsIntegral = 0;
    double bulkHitsIntegral = 0;
    int nVertexHits = 0;
    int nBulkHits = 0;
    // loop over the hits
    for(SHit &h:hitList){
        double d = std::hypot( h.X() - xVertex - xFitVtx, h.Y() - yVertex - yFitVtx );
        if(d<2.5){
            vertexHitsIntegral+=h.Integral();
            nVertexHits++;
        }
        else if(d>2.5 && d <10){
            bulkHitsIntegral+=h.Integral();
            nBulkHits++;
        }
    }
    std::cout<<"CHECK 1\n";

    double vertexCharge = vertexHitsIntegral/nVertexHits;
    double bulkCharge = bulkHitsIntegral/nBulkHits;
    
    fVertexHitIntegralRatio = vertexCharge/bulkCharge;
    fVertexHitIntegralDifference = vertexCharge-bulkCharge;
    fVertexHitIntegralRelativeDifference = fVertexHitIntegralDifference/bulkCharge;

    /*SVertex fVertexFit = SVertex( SPoint( xFitVtx+xVertex, yFitVtx+yVertex ), "");
    chargeDensityAlgo.Fill(hitList, fVertexFit);
    TCanvas *cDisplay = new TCanvas( "FinalRecoFRANS", "FinalRecoFRANS", 700, 0, 900, 1200);
    chargeDensityAlgo.Display(cDisplay);
    double eta = chargeDensityAlgo.Eta();
    std::cout<<"Eta: "<<eta<<std::endl;*/
    std::cout<<"CHECK 1\n";
    delete c1;
    //delete pad1;
    //delete pad2;
    std::cout<<"CHECK 1\n";
    return;

}