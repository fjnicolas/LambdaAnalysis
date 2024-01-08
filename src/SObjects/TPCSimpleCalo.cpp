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
#include "TStyle.h"

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

    if(fHitList.size()>1){
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
    else{
        fPathLengths.push_back(0);
    }

}

void SCalo::CalculateResidualRange() {
    
    CalculatePathLengths();
    
    fDepositedEnergy.clear();
    fResidualRange.clear();

    fTrackLength = 0;
    std::cout<<" fHitList.size() "<<fHitList.size()<<" PP: "<<fPathLengths.size()<<std::endl;
    if(fPathLengths.size()>1){
        fTrackLength = std::accumulate(fPathLengths.begin(), fPathLengths.end(), fTrackLength);

        double pathLength = 0;
        for (size_t i = 0; i < fHitList.size(); ++i) {
            fResidualRange.push_back( fTrackLength - pathLength);
            fDepositedEnergy.push_back( (fHitList[i].Integral() * fHitIntegralToEnergy) / fPathLengths[i]);
            //fDepositedEnergy.push_back( (fHitList[i].Integral() * fHitIntegralToEnergy) * fCosTrackAngle );
            pathLength+=fPathLengths[i];
        }

        // Check last points
        double minE = *std::min_element(fDepositedEnergy.begin(), fDepositedEnergy.end());
        std::cout << " minE " << minE << " " << fDepositedEnergy.front() << " " << fDepositedEnergy.back() << std::endl;

        if(fResidualRange.size()>3){
            while (!fDepositedEnergy.empty() && fDepositedEnergy.back() <= minE) {
                fDepositedEnergy.pop_back();
                fResidualRange.pop_back();
                fPathLengths.pop_back();
                fHitList.pop_back();
                minE = *std::min_element(fDepositedEnergy.begin(), fDepositedEnergy.end());
            }
        }

        double lastResidualRange = fResidualRange.back();
        std::cout << " lastResidualRange " << lastResidualRange << std::endl;
        // reshift 
        for(size_t i=0; i<fResidualRange.size(); i++){
            fResidualRange[i] = fResidualRange[i] - lastResidualRange;
        }
        fTrackLength = fTrackLength - lastResidualRange;
    }

    if(fResidualRange.size()>=2){
        
        std::vector<double> rangeDiff;
        for(size_t i=0; i<fResidualRange.size()-1; i++){
            rangeDiff.push_back( fResidualRange[i+1] - fResidualRange[i] );
        }
        double sq_sum = std::inner_product(rangeDiff.begin(), rangeDiff.end(), rangeDiff.begin(), 0.0);
        double mean = std::accumulate(rangeDiff.begin(), rangeDiff.end(), 0.0) / rangeDiff.size();
        fResidualRangeRMS = std::sqrt(sq_sum / rangeDiff.size() - mean * mean);

    }
    else{
        fResidualRangeRMS = -1;
    }



    

}


void SCalo::Display(){
    std::cout<<"Displaying calorimetry...track length: "<<fTrackLength<<"\n";

    for (size_t i = 0; i < fPathLengths.size(); ++i) {
        std::cout << " Residual range / dEdx: " << fResidualRange[i] << " / " << fDepositedEnergy[i] << std::endl;
    }

}


// ------ STriangleCalo class ------

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


// --- Constructor ---
STriangleCalo::STriangleCalo(STriangle triangle):
    fTriangle(triangle)
{
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

// --- Function to check overlap between two contours plots ---
double CheckOverlapRegion(TH1F *h1, TH1F *h2, double x1Max, double x2Max){
    double overlapArea = 0;

    double x1Min = h1->GetXaxis()->GetXmin();
    double x2Min = h2->GetXaxis()->GetXmin();
    
    double x1BinWidth = h1->GetXaxis()->GetBinWidth(1);
    double x2BinWidth = h2->GetXaxis()->GetBinWidth(1);

    double xMax = std::min(x1Max, x2Max);
    // First get total area of the two histograms in the range min-max
    // Now get the overlap
    double x1 = x1Min;
    double x2 = x2Min;    

    double area1 = 0;
    double area2 = 0;
    while(x1<xMax){
        area1+=2*h1->GetBinError(h1->FindBin(x1))*x1BinWidth;
        x1+=x1BinWidth;
    }
    while(x2<xMax){
        area2+=2*h2->GetBinError(h2->FindBin(x2))*x2BinWidth;
        x2+=x2BinWidth;
    }

    // Now get the overlap
    x1 = x1Min;
    x2 = x2Min;    
    while(x1<xMax && x2<xMax){
        // Get y values with errors
        double y1Error = h1->GetBinError(h1->FindBin(x1));
        double y1 = h1->GetBinContent(h1->FindBin(x1));
        double y2Error = h2->GetBinError(h2->FindBin(x2));
        double y2 = h2->GetBinContent(h2->FindBin(x2));

        double y1Low = y1-y1Error;
        double y1High = y1+y1Error;
        double y2Low = y2-y2Error;
        double y2High = y2+y2Error;

        // Calculate the overlap integral
        double overlap = ( std::min(y1High, y2High) - std::max(y1Low, y2Low) ) * x1BinWidth;
        if (overlap < 0.0) overlap = 0.0;
        //std::cout<<" x1: "<<x1<<" x2: "<<x2<<" Y1: "<<y1Low<<":"<<y1High<<" Y2: "<<y2Low<<":"<<y2High<<" overlap: "<<overlap<<std::endl;

        x1+=x1BinWidth;
        x2+=x2BinWidth;
        overlapArea+=2*overlap; // overlap counts to the two ranges
    }
    std::cout<<" overlapArea: "<<overlapArea<<" area1: "<<area1<<" area2: "<<area2<<std::endl;
    return overlapArea/(area1+area2);
}

// --- Function to display dEdx ---
void STriangleCalo::MakeEnergyLossVsResidualRangePlot(SCalo calo1, SCalo calo2, TPad *pad) {
    
    if(pad!=nullptr){
        pad->cd();
        gStyle->SetOptStat(0);
    }
    double fCondifenceLevel = 0.95;
    
    // Get vectors
    const std::vector<double> residualRange1 = calo1.GetResidualRange();
    const std::vector<double> depositedEnergy1 = calo1.GetDepositedEnergy();
    const std::vector<double> residualRange2 = calo2.GetResidualRange();
    const std::vector<double> depositedEnergy2 = calo2.GetDepositedEnergy();

    // Define graph 1
    TGraph* graph1 = new TGraph(residualRange1.size(), &residualRange1[0], &depositedEnergy1[0]);
    graph1->SetTitle("Track 1");
    graph1->SetLineColor(fColor1);
    graph1->SetLineStyle(kDashed);
    graph1->SetMarkerColor(fColor1);
    graph1->SetMarkerStyle(20);

    // Define graph 2
    TGraph* graph2 = new TGraph(residualRange2.size(), &residualRange2[0], &depositedEnergy2[0]);
    graph2->SetTitle("Track 2");
    graph2->SetLineColor(fColor2);
    graph2->SetMarkerColor(fColor2);
    graph2->SetLineStyle(kDashed);
    graph2->SetMarkerStyle(20);

    // Frame
    double maxX = std::max(calo1.GetTrackLength(), calo2.GetTrackLength());   
    double maxY1 = 0, maxY2 = 0;
    if(depositedEnergy1.size()>1)
        maxY1 = *std::max_element(depositedEnergy1.begin(), depositedEnergy1.end());
    if(depositedEnergy2.size()>1)
        maxY2 = *std::max_element(depositedEnergy2.begin(), depositedEnergy2.end()); 
    double maxY = std::max(maxY1, maxY2);
    
    TH2F *hFrame = new TH2F("hFrame", ";Residual range [cm];dQ/dx [AU]", 200, 0, 1.1*maxX, 100, 0, 1.1*maxY);
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

    TVirtualFitter* fitter = TVirtualFitter::GetFitter();
    
    // Fit the two graphs to two constants
    bool passFit1 = false;
    bool passFit1Exp = false;
    bool passFit2 = false;
    bool passFit2Exp = false;
    double maxLength = std::max(calo1.GetTrackLength(), calo2.GetTrackLength());
    TH1F *hint1 = new TH1F("hint1", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    TH1F *hint1Exp = new TH1F("hint1Exp", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    TH1F *hint2 = new TH1F("hint2", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    TH1F *hint2Exp = new TH1F("hint2Exp", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    hint1->SetFillColorAlpha(fColor1,0.5);
    hint1Exp->SetFillColorAlpha(fColor1,0.5);
    hint2->SetFillColorAlpha(fColor2,0.5);
    hint2Exp->SetFillColorAlpha(fColor2,0.5);

    TF1 *fit1 = new TF1("fit1", "[0]", 0, maxX);
    fit1->SetLineColor(fColor1);
    fit1->SetLineWidth(2);
    graph1->Fit(fit1, "R");
    fitter = TVirtualFitter::GetFitter();
    if(fitter){
        fitter->GetConfidenceIntervals(hint1, fCondifenceLevel);
        passFit1 = true;
        std::cout<<"Pass fit 1"<<std::endl;
    }
    
    TF1 *fit1Exp = new TF1("fit1Exp", "[0]*exp(x*[1])+[2]", 0, calo1.GetTrackLength());
    fit1Exp->SetLineColor(fColor1);
    fit1Exp->SetLineWidth(2);
    fit1Exp->SetParameters(maxY, -1, fit1->GetParameter(0));
    fit1Exp->SetParLimits(1, -1e6, 0);
    fit1Exp->SetParLimits(2, 0, 1e6);
    graph1->Fit(fit1Exp, "BR");    
    fitter = TVirtualFitter::GetFitter();
    if(fitter){
        fitter->GetConfidenceIntervals(hint1Exp, fCondifenceLevel);
        passFit1Exp = true;
        std::cout<<"Pass fit 1 exp"<<std::endl;
    }
    

    TF1 *fit2 = new TF1("fit2", "[0]", 0, calo2.GetTrackLength());
    fit2->SetLineColor(fColor2);
    fit2->SetLineWidth(2);
    graph2->Fit(fit2, "R");
    fitter = TVirtualFitter::GetFitter();
    if(fitter){
        fitter->GetConfidenceIntervals(hint2, fCondifenceLevel);
        passFit2 = true;
        std::cout<<"Pass fit 2"<<std::endl;
    }
    TF1 *fit2Exp = new TF1("fit2Exp", "[0]*exp(x*[1])+[2]", 0, calo2.GetTrackLength());
    fit2Exp->SetLineColor(fColor2);
    fit2Exp->SetLineWidth(2);
    fit2Exp->SetParameters(maxY, -1, fit2->GetParameter(0));
    fit2Exp->SetParLimits(1, -1e6, 0);
    fit2Exp->SetParLimits(2, 0, 1e6);
    graph2->Fit(fit2Exp, "BR");
    fitter = TVirtualFitter::GetFitter();
    if(fitter){
        fitter->GetConfidenceIntervals(hint2Exp, fCondifenceLevel);
        passFit2Exp = true;
        std::cout<<"Pass fit 2 exp"<<std::endl;
    }

    // Draw the fitted functions
    fit1->Draw("same");
    fit2->Draw("same");
    fit1Exp->Draw("same");
    fit2Exp->Draw("same");
    double p0Fit_1 = 1;
    double p0Fit_2 = 1;
    if(fit1Exp->GetParameter(0)>0 && fit1Exp->GetParameter(1)<0 && hint1Exp->Integral()>0){
        hint1Exp->Draw("e3 same");
        p0Fit_1 = fit1Exp->GetParameter(2);
    }
    else{
        hint1->Draw("e3 same");
        p0Fit_1 = fit1->GetParameter(0);
    }
    if(fit2Exp->GetParameter(0)>0 && fit2Exp->GetParameter(1)<0 && hint2Exp->Integral()>0){
        hint2Exp->Draw("e3 same");
        p0Fit_2 = fit2Exp->GetParameter(2);
    }
    else{
        hint2->Draw("e3 same");
        p0Fit_2 = fit2->GetParameter(0);
    }

    double charge1=0, charge2=0;
    if(depositedEnergy1.size()>1)
        charge1 = std::accumulate(depositedEnergy1.begin(), depositedEnergy1.end(), 0.0) / depositedEnergy1.size();
    if(depositedEnergy2.size()>1)
        charge2 = std::accumulate(depositedEnergy2.begin(), depositedEnergy2.end(), 0.0) / depositedEnergy2.size();
    double charge1Average = std::max(charge1, charge2);
    double charge2Average = std::min(charge1, charge2);
    fChargeRatioAverage = charge1Average / charge2Average;
    fChargeDifferenceAverage = charge2Average - charge1Average;
    fChargeRelativeDifferenceAverage = fChargeDifferenceAverage / charge1Average;


    bool pass1 = passFit1 || passFit1Exp;
    bool pass2 = passFit2 || passFit2Exp;

    fPassChargeFit = pass1 && pass2;
    if(fPassChargeFit){
        double charge1Fit = std::max(p0Fit_1, p0Fit_2);
        double charge2Fit = std::min(p0Fit_1, p0Fit_2);

        fChargeRatioFit = charge1Fit / charge2Fit;
        fChargeDifferenceFit = charge2Fit - charge1Fit;
        fChargeRelativeDifferenceFit = fChargeDifferenceFit / charge1Fit;

        // Check the overlap
        double overlap = CheckOverlapRegion(hint1Exp, hint2Exp, calo1.GetTrackLength(), calo2.GetTrackLength());
        std::cout<<"Overlap: "<<overlap<<std::endl;
        fBandOverlap = overlap;
    }
    else{
        fChargeRatioFit = 1;
        fChargeDifferenceFit = 0;
        fChargeRelativeDifferenceFit = 0;
        fBandOverlap = 1;
    }

}


void STriangleCalo::CreateEnergyLossVsResidualRangePlot() {
    TCanvas *c1 = new TCanvas("canvasCalo","Calorimetry", 600,1200);
    c1->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.02,0.98,0.98);
    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, pad1);
    delete c1;
    delete pad1;
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

// --- Get mean and std deviation of vector of hits ---
void GetStatisticsHitVector(std::vector<SHit> hits, double & mean, double & stddev){
    double sum = 0;
    double sum2 = 0;
    for(SHit & h:hits){
        sum+=h.Integral();
        sum2+=h.Integral() * h.Integral();
    }
    sum/=hits.size();
    sum2/=hits.size();
    mean = sum;
    stddev = std::sqrt(sum2-sum*sum);

    return;

}

// --- Right/left hits with respect to vertex ---
std::vector<SHit> GetSideHits(std::vector<SHit> hitV, double xVertex, std::vector<SHit> & residualHits){

    std::vector<SHit> vLeft, vRight;
    for(SHit & h:hitV){
        if(h.X()>=xVertex) vRight.push_back(h);
        else  vLeft.push_back(h);
    }

    if(vLeft.size()>vRight.size()) {
        residualHits.clear();
        residualHits = vRight;
        return vLeft;
    }
    else{
        residualHits.clear();
        residualHits = vLeft;
        return vRight;
    }
}



// --- JointFit analysis---
void STriangleCalo::JointFitAnalysis(unsigned int maxHits, double widthTol, bool useHitError){

    // --- Get the hits and the vertex ---
    // Vertex
    fXVertex = fTriangle.GetMainVertex().X();
    fYVertex = fTriangle.GetMainVertex().Y();

    // Hit coordinates for the first maxHits hits
    maxHits = std::min( maxHits, (unsigned int)std::min(fTriangle.GetNHitsTrack1(), fTriangle.GetNHitsTrack2()) );
    std::vector<double> xV, yV, widthV;
    size_t n = 0;
    for(SHit &h:fTriangle.GetTrack1().GetHits()){
        if(n>=maxHits) break;
        xV.push_back(h.X()-fXVertex);
        yV.push_back(h.Y()-fYVertex);
        widthV.push_back(h.Width());
        n++;
    }
    n=0;
    for(SHit &h:fTriangle.GetTrack2().GetHits()){
        if(n>=maxHits) break;
        xV.push_back(h.X()-fXVertex);
        yV.push_back(h.Y()-fYVertex);
        widthV.push_back(h.Width());
        n++;
    }

    // --- Get frame values ---
    fMinX = *std::min_element(xV.begin(), xV.end());
    fMaxX = *std::max_element(xV.begin(), xV.end());
    fMinY = *std::min_element(yV.begin(), yV.end());
    fMaxY = *std::max_element(yV.begin(), yV.end());

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

    // --- Fit results ---
    const ROOT::Fit::FitResult & result = fitter.Result();
    fM1 = result.Parameter(0);
    fN1 = result.Parameter(1);
    fM2 = result.Parameter(2);
    fXB = result.Parameter(3);
    fN2 = (fM1-fM2)*fXB+fN1;
    fXFitVtx = fXB;
    fYFitVtx = fM1*fXFitVtx+fN1;
    result.Print(std::cout);

    // --- Assing each hit to a track ---
    fMaxIntegral = 0;
    fHitsTrack1.clear();
    fHitsTrack2.clear();
    fVertexHits.clear();
    std::vector<SHit> hitList = fTriangle.GetAllHits();
    fAllHits = hitList;
    for(SHit &h:hitList){

        // shifted values 
        double x = h.X()-fXVertex;
        double y = h.Y()-fYVertex;

        if(x<fMinX) fMinX = x;
        if(y<fMinY) fMinY = y;
        if(x>fMaxX) fMaxX = x;
        if(y>fMaxY) fMaxY = y;

        // values from fitted slopes
        double yH1 = fM1*x+fN1;
        double yH2 = fM2*x+fN2;

        // closest tracks
        bool isTrack1 = std::abs( y-yH1 )<widthTol*h.Width();
        bool isTrack2 = std::abs( y-yH2 )<widthTol*h.Width();

        // if shared hit, split the integral in two hits
        if(isTrack1 && isTrack2){
            //SHit sharedHit(-1, h.X(), h.Y(), h.Width(), 0.5*h.Integral());
            //fHitsTrack1.push_back(sharedHit);
            //fHitsTrack2.push_back(sharedHit);
            fVertexHits.push_back(h);
        }
        else if(isTrack1){
            fHitsTrack1.push_back(h);
        }
        else if(isTrack2){
            fHitsTrack2.push_back(h);
        }
        else{
            // if not width compatible, add to the closest
            double d1 = std::abs(y-yH1);
            double d2 = std::abs(y-yH2);
            d1<d2? fHitsTrack1.push_back(h):fHitsTrack2.push_back(h);
        }

        // update maximum integral
        if(h.Integral()>fMaxIntegral) fMaxIntegral = h.Integral();
    }

    // --- Remove outliers hits ---
    std::vector<SHit> resHits;
    fHitsTrack1 = GetSideHits(fHitsTrack1, fXVertex+fXFitVtx, resHits);
    if(resHits.size()>0){
        fResidualHits.insert(fResidualHits.end(), resHits.begin(), resHits.end());
        resHits.clear();
    }
    fHitsTrack2 = GetSideHits(fHitsTrack2, fXVertex+fXFitVtx, resHits);
    if(resHits.size()>0){
        fResidualHits.insert(fResidualHits.end(), resHits.begin(), resHits.end());
        resHits.clear();
    }
    // Average track 1 hit intetgral
    double meanIntegral1 = 0, stdDev1 = 0;
    double meanIntegral2 = 0, stdDev2 = 0;
    GetStatisticsHitVector(fHitsTrack1, meanIntegral1, stdDev1);
    GetStatisticsHitVector(fHitsTrack2, meanIntegral2, stdDev2);

    double meanIntegralShared = 0, stdDevShared = 0;
    GetStatisticsHitVector(fVertexHits, meanIntegralShared, stdDevShared);
    std::vector<SHit> auxVertexHits = fVertexHits;
    fVertexHits.clear();
    for(SHit & h:auxVertexHits){
        fVertexHits.push_back(h);
    }

    // --- Calorimetry objects ---
    // Sort by distance to fitted vertex ---
    SortSHitVectorByDistance(fHitsTrack1, fXVertex+fXFitVtx, fYVertex+fYFitVtx);
    SortSHitVectorByDistance(fHitsTrack2, fXVertex+fXFitVtx, fYVertex+fYFitVtx);
    SortSHitVectorByDistance(fVertexHits, fXVertex+fXFitVtx, fYVertex+fYFitVtx);

    //fHitsTrack1.pop_back(); fHitsTrack2.pop_back();
    // Create the objects
    SCalo calo1(fHitsTrack1, std::atan(fM1));
    SCalo calo2(fHitsTrack2, std::atan(fM2));
    fCalo1 = calo1;
    fCalo2 = calo2;
    std::cout<<" Set calo objects\n";
    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, nullptr);
   
    // --- Return values ---
    fTwoLinesChi2 = result.MinFcnValue();

    fMinX-=0.1*(fMaxX-fMinX);
    fMaxX+=0.1*(fMaxX-fMinX);
    fMinY-=0.1*(fMaxY-fMinY);
    fMaxY+=0.1*(fMaxY-fMinY);

    // vertex charge vs track charge
    SortSHitVectorByDistance(hitList, fXVertex+fXFitVtx, fYVertex+fYFitVtx);
    double vertexHitsIntegral = 0;
    double bulkHitsIntegral = 0;
    int nVertexHits = 0;
    int nBulkHits = 0;
    // loop over the hits
    for(SHit &h:fVertexHits){
        vertexHitsIntegral+=h.Integral();
    }
    nVertexHits = fVertexHits.size();
    int nCounter=0;
    for(SHit &h:fHitsTrack1){
        if(nCounter>(int)fHitsTrack1.size()/2) continue;
        bulkHitsIntegral+=h.Integral();
        nBulkHits++;
        nCounter++;
    }
    nCounter=0;
    for(SHit &h:fHitsTrack2){
        if(nCounter>(int)fHitsTrack2.size()/2) continue;
        bulkHitsIntegral+=h.Integral();
        nBulkHits++;
        nCounter++;
    }

    double vertexCharge = vertexHitsIntegral/nVertexHits;
    double bulkCharge = bulkHitsIntegral/nBulkHits;
    
    fNVertexHits = nVertexHits;
    fNBulkHits = nBulkHits;
    if(nVertexHits>0 and nBulkHits>0){
        fVertexHitIntegralRatio = vertexCharge/bulkCharge;
        fVertexHitIntegralDifference = vertexCharge-bulkCharge;
        fVertexHitIntegralRelativeDifference = fVertexHitIntegralDifference/bulkCharge;
    }
    else{
        fVertexHitIntegralRatio = 0;
        fVertexHitIntegralDifference = 0;
        fVertexHitIntegralRelativeDifference = 0;
    }

    double trackLength1 = calo1.GetTrackLength();
    double trackLength2 = calo2.GetTrackLength();
    fTrackLength1 = std::max(trackLength1, trackLength2);
    fTrackLength2 = std::min(trackLength1, trackLength2);
    fTrackLengthRatio = fTrackLength1/fTrackLength2;
    fResidualRange1RMS = calo1.GetResidualRangeRMS();
    fResidualRange2RMS = calo2.GetResidualRangeRMS();

    std::cout << "NVertexHits: " << fNVertexHits <<std::endl;
    std::cout << "NBulkHits: " << fNBulkHits <<std::endl;
    std::cout << "VertexHitIntegralRatio: "<<fVertexHitIntegralRatio<<std::endl;
    std::cout << "VertexHitIntegralDifference: "<<fVertexHitIntegralDifference<<std::endl;
    std::cout << "VertexHitIntegralRelDifference: "<<fVertexHitIntegralRelativeDifference<<std::endl;
    std::cout << "--- TrackLengthVariables: "<<std::endl;
    std::cout << "TrackLengthRatio: "<<fTrackLengthRatio<<std::endl;
    std::cout << "Track 1 length: "<<fTrackLength1<<std::endl;
    std::cout << "Track 2 length: "<<fTrackLength2<<std::endl;
    std::cout << "Track 1 residual range RMS: "<<fResidualRange1RMS<<std::endl;
    std::cout << "Track 2 residual range RMS: "<<fResidualRange2RMS<<std::endl;
    std::cout << "--- Charge variables: "<<std::endl;
    std::cout << "Charge ratio (average): " << fChargeRatioAverage << std::endl;
    std::cout << "Charge difference (average): " << fChargeDifferenceAverage << std::endl;
    std::cout << "Charge relative difference (average): " << fChargeRelativeDifferenceAverage << std::endl;
    std::cout << "Charge ratio (fit) variables: "<<std::endl;
    std::cout << "Pass charge fit: " << fPassChargeFit << std::endl;
    std::cout << "Band overlap: " << fBandOverlap << std::endl;
    std::cout << "Charge ratio (fit): " << fChargeRatioFit << std::endl;
    std::cout << "Charge difference (fit): " << fChargeDifferenceFit << std::endl;
    std::cout << "Charge relative difference (fit): " << fChargeRelativeDifferenceFit << std::endl;
    


    return;

}


// --- Display function ---
void STriangleCalo::Display(TCanvas *c1){
    // --- Create the canvas ---
    if(c1==nullptr) {
        c1 = new TCanvas("canvasCalo","Calorimetry", 600,1200);
    }
    c1->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.98,0.48);
    pad1->Draw();
    pad2->Draw();

    // --- Draw the graph and the fit ---
    pad1->cd();
    TH2F * hFrame = new TH2F("hFrame", ";Wire ID;Time Tick [0.5 #mu s]", 200, fMinX, fMaxX, 200, fMinY, fMaxY);
    hFrame->GetYaxis()->SetTitleOffset(.8);
    hFrame->GetXaxis()->SetTitleOffset(.8);
    hFrame->SetStats(0);

    // Line1
    double y1 = fM1*fMinX+fN1;
    double y2 = fM1*fMaxX+fN1;
    TLine *line1 = new TLine(fMinX,y1,fMaxX,y2);
    line1->SetLineColor(fColor1);
    line1->SetLineWidth(4); 

    // Line2
    y1 = fM2*fMinX+fN2;
    y2 = fM2*fMaxX+fN2;
    TLine *line2 = new TLine(fMinX,y1,fMaxX,y2);
    line2->SetLineColor(fColor2);
    line2->SetLineWidth(4);

    // Fitted vertex
    TMarker *marker = new TMarker(fXB, fM1*fXB+fN1, 29);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(2);

    // TLegend
    TLegend *leg = new TLegend(0.7, 0.4, 0.85, 0.65);
    leg->AddEntry(line1, "Fitted line 1", "l");
    leg->AddEntry(line2, "Fitted line 2", "l");
    leg->AddEntry(marker, "Fitted vertex", "p");

    // Create text with the fit results
    TPaveText *pt = new TPaveText(0.7, 0.7, 0.85, 0.85, "NDC");
    pt->AddText(Form("Fit result #chi^{2} = %.2f\n",fTwoLinesChi2));
    pt->AddText(Form("m1 = %.2f #pm %.2f",fM1, 0.));
    pt->AddText(Form("n1 = %.2f #pm %.2f",fN1, 0.));
    pt->AddText(Form("m2 = %.2f #pm %.2f",fM2, 0.));
    pt->AddText(Form("xB = %.2f #pm %.2f",fXB, 0.));

    // Draw
    hFrame->Draw();
    line1->Draw("same");
    line2->Draw("same");
    marker->Draw("same");
    pt->Draw("same");
    
    std::cout<<" Drawing calorimetry...track length: "<<fHitsTrack1.size()<<"\n";
    // --- Draw xV and yV markers after assignment ---

    for(SHit &h:fAllHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kGray, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
    }

    for(SHit &h:fResidualHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kGreen, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
    }

    for(SHit &h:fHitsTrack1){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 22);
        marker->SetMarkerColorAlpha(fColor1, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        if (&h == &fHitsTrack1.front()) leg->AddEntry(marker, "Track 1 hits", "p");
    }

    for(SHit &h:fHitsTrack2){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 23);
        marker->SetMarkerColorAlpha(fColor2, 1);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        if (&h == &fHitsTrack2.front()) leg->AddEntry(marker, "Track 2 hits",  "p");
    }

    for(SHit &h:fVertexHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kViolet-6, 1);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        if (&h == &fVertexHits.front()) leg->AddEntry(marker, "Vertex hits",  "p");
    }

    leg->Draw("same");

    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, pad2); 



    c1->Update();
    c1->cd();

    return;
}