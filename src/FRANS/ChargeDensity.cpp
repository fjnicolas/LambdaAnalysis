#include "ChargeDensity.h"


ChargeDensity::ChargeDensity(FRANSPsetType const& config, int view)
  : fFRANSPset(config),
    fView(view)
{

  fNormUnAvSmooth=1./(2*fFRANSPset.UnAvNeighbours+1);

  // default values
  fDelta = 0.5;
  fEta = 1.;
  fFitScore = 1.;
  fAlpha = 10e3;
  fOmega = -1e4;
  fTau = 0.;
  fIota = 10e3;
  fScore = -1e4;

  fNHits = 0;
  fMeanChi2 = 0;

  // Right now we only use the collection view
  if(fFRANSPset.CalculateScore){
    if(fFRANSPset.UseAlpha)
      fTMVAReader.AddVariable( "FRANSObj"+std::to_string(view)+".fAlpha", &fAlpha );
    fTMVAReader.AddVariable( "FRANSObj"+std::to_string(view)+".fEta", &fEta );
    fTMVAReader.AddVariable( "FRANSObj"+std::to_string(view)+".fDelta", &fDelta);
    fTMVAReader.AddVariable( "FRANSObj"+std::to_string(view)+".fFitScore", &fFitScore );

    fTMVAReader.AddSpectator( "Gap", &fGap );
    fTMVAReader.AddSpectator( "ProtonKE", &fProtonKE );
    fTMVAReader.AddSpectator( "PionKE", &fPionKE );

    fTMVAReader.BookMVA( "FRANS BDT",  fFRANSPset.TMVAFilename.c_str()  );

  }
}


// ----------- Distance function -----------
double ChargeDensity::GetDistance(int x, double y, int x0, double y0) {
    const double deltaX = static_cast<double>(x - x0) / fFRANSPset.NWirePack;
    const double deltaY = (y - y0) / fFRANSPset.NDriftPack;
    return std::hypot(deltaX, deltaY);
}


// ----------- Distance function wires -----------
int ChargeDensity::GetDistanceWires(int x, int x0) {
    return static_cast<int> ( std::abs(x - x0) / fFRANSPset.NWirePack );
}


void ChargeDensity::Rescale(std::vector<double>& wf, double kappa){
  std::transform(wf.begin(), wf.end(), wf.begin(), [&kappa](double& z){return z*kappa;});
}


//---- Expo average smoothing
void ChargeDensity::ApplyExpoAvSmoothing(std::vector<double>& wf){
  std::transform (std::next(wf.begin(), 1), wf.end(), wf.begin(), std::next(wf.begin(), 1),
    [&](double _x, double _y) { return  fFRANSPset.ExpoAvSmoothPar*_x+ (1. - fFRANSPset.ExpoAvSmoothPar)*_y; }  );
}


//---- Unweighted average smoothing
void ChargeDensity::ApplyUnAvSmoothing(std::vector<double>& wf){
  std::vector<double> wf_aux(wf.begin(), wf.end());
  for(size_t bin=fFRANSPset.UnAvNeighbours; bin<wf.size()-fFRANSPset.UnAvNeighbours; bin++){
    double sum=0.;
    for(size_t nbin=bin-fFRANSPset.UnAvNeighbours; nbin<=bin+fFRANSPset.UnAvNeighbours; nbin++)
      sum+=wf_aux[nbin];
    wf[bin]=sum*fNormUnAvSmooth;
  }
}


//---- Sliding window
void ChargeDensity::SlidingWindow(std::vector<double>& wf){

  fZCumDer.clear(); fZCumDer.reserve(fZCum.size());
  fZCumDerErr.clear(); fZCumDerErr.reserve(fZCum.size());

  std::vector<float> x; x.reserve(fFRANSPset.SlidingWindowN);
  for(int k=0; k<fFRANSPset.SlidingWindowN; k++) x.push_back(k);

  std::vector<std::vector<float>> Z;
  Z.push_back(x); Z.push_back(x);

  Linear_Regression <float,float> Reg_Class;

  for(size_t k=0; k<wf.size()-fFRANSPset.SlidingWindowN; k++){

    std::vector<float> y;
    for(size_t j=k; j<k+fFRANSPset.SlidingWindowN; j++){
      y.push_back(wf[j]);
    }

    Z[1]=y;

    Reg_Class.fit(Z);

    //std::cout<<" Slope: "<<k<<" "<<Reg_Class.Slope()<<std::endl;

    fZCumDer.push_back(Reg_Class.Slope());
    fZCumDerErr.push_back(Reg_Class.SlopeErr());
  }

}


void ChargeDensity::FillCumulative(){
  fZCum.clear();
  fZCum.reserve(fZ.size());

  fMaxCumulative = std::accumulate(fZ.begin(), fZ.end(), 0.);

  fZCum.push_back(fZ[0]);
  for(size_t ix=1; ix<fZ.size(); ix++){
    fZCum.push_back( fZCum[ix-1] + fZ[ix] );
  }

  //for(size_t ix=0; ix<fZCum.size(); ix++){std::cout<<ix<<" "<<fZCum[ix]<<std::endl;}

  size_t start = std::upper_bound(fZCum.begin(), fZCum.end(), fFRANSPset.CumulativeCut*fMaxCumulative) - fZCum.begin();

  fZ.resize(start);
  fZCum.resize(start);
}


// ----------- Update Metrics -----------
void ChargeDensity::UpdateMetrics(){

  // --- absolute normalization
  fZetaNorm = *std::max_element(fZ.begin(), fZ.end());

  if(fFRANSPset.Verbose>0) std::cout<<"Updating metrics...size is: "<<fZ.size()<<std::endl;
  if(fFRANSPset.Verbose>0) { std::cout<<"Debug 0 (pre raw smoothing): "; for(size_t l=0; l<10; l++) {std::cout<<fZ[l]<<" ";} }  
  if(fFRANSPset.Verbose>0) std::cout<<"\nZetaNorm: "<<fZetaNorm<<std::endl;

  // --- apply raw smoothing
  if(fFRANSPset.ApplyRawSmoothing) {
    // keep first point without smoothing
    std::vector<double> Zsmoothed(std::next(fZ.cbegin(), 1), fZ.cend());
    ApplyExpoAvSmoothing(Zsmoothed); //!!!!! Check this
    Zsmoothed.insert(Zsmoothed.begin(), fZ[0]);
    fZ.clear();
    fZ.insert( fZ.begin(), Zsmoothed.begin(), Zsmoothed.end());
  }

  double ZetaNormI = 1./( *std::max_element(fZ.begin(), fZ.end()) );
  Rescale(fZ, ZetaNormI);
  
  if(fFRANSPset.Verbose>0){
    std::cout<<"\nDebug 1 (after raw smoothing): ";
    for(size_t l=0; l<10; l++) {std::cout<<fZ[l]<<" ";}
  }

  // --- fill cumulative and resize vector keeping only ROI
  FillCumulative();

  // --- apply the sliding window algorithm
  if((int)fZ.size()>fFRANSPset.SlidingWindowN){
    SlidingWindow(fZCum);

    if(fFRANSPset.ApplySmoothing){
      ApplyExpoAvSmoothing(fZCumDer);
      ApplyUnAvSmoothing(fZCumDer);
    }

    Rescale(fZCum, 1./fMaxCumulative);
    if(fFRANSPset.Verbose>0) { std::cout<<"\n\nDebug 2 (cumulative): "; for(size_t l=0; l<10; l++) {std::cout<<fZCum[l]<<" ";} }

    fRho.clear();
    fRho.reserve(fZ.size());
    for(size_t k=0; k<fZ.size(); k++){
      fRho.push_back(k);
    }

    //--- Now we start with the parameter estimation
    // first get ROIs
    int L = std::min(fFRANSPset.MaxRadius, (int)fZCumDer.size());
    std::vector<double> Slopes(fZCumDer.cbegin(), std::next(fZCumDer.cbegin(), L));
    std::vector<double> SlopesErr(fZCumDerErr.cbegin(), std::next(fZCumDerErr.cbegin(), L));
    std::cout<<"L "<<L<<" "<<fFRANSPset.MaxRadius<<std::endl;
    if(fFRANSPset.Verbose>0) { std::cout<<"\n\nDebug 3 (slopes): "; for(size_t l=0; l<10; l++) {std::cout<<Slopes[l]<<" ";} }

    // get first big jump in the slope and maximum slope indexes
    int ix_firstJump=-1;
    size_t ix_maxSlope=0;
    double maxSlope=0;
    for(size_t k=0; k<Slopes.size(); k++){
      if(Slopes[k]>maxSlope){
        maxSlope=Slopes[k];
        ix_maxSlope=k;
        //std::cout<<" New max slope at "<<ix_maxSlope<<" "<<maxSlope<<std::endl;
      }
      if(SlopesErr[k]>0.02 && ix_firstJump==-1){
        ix_firstJump=(int)k;
      }
    }

    if(fFRANSPset.Verbose>0) std::cout<<"\n\nFirtJump @ "<<ix_firstJump<<" MaxSlope @ "<<ix_maxSlope<<std::endl;
    double startSlope = 1.;
    double startNHits = 1e4;
    if(fFRANSPset.Verbose>0) std::cout<<"\nStart Slope: ";
    if(ix_firstJump>0){
      startSlope = std::accumulate(Slopes.cbegin(), std::next(Slopes.cbegin(), ix_firstJump), 0.)/(ix_firstJump);
      for(int k=0; k<ix_firstJump+1; k++) std::cout<<k<<":"<<Slopes[k]<<"  ";
      startNHits = std::accumulate(fZCounter.cbegin(), std::next(fZCounter.cbegin(), ix_firstJump), 0.)/(ix_firstJump);
    }
    else{
      startSlope = std::accumulate(Slopes.cbegin(), std::next(Slopes.cbegin(), fFRANSPset.NSamplesBeginSlope), 0.)/(fFRANSPset.NSamplesBeginSlope);
      for(int k=0; k<fFRANSPset.NSamplesBeginSlope; k++) std::cout<<k<<":"<<Slopes[k]<<"  ";
      startNHits = std::accumulate(fZCounter.cbegin(), std::next(fZCounter.cbegin(), fFRANSPset.NSamplesBeginSlope), 0.)/(fFRANSPset.NSamplesBeginSlope);
    }

    if(fFRANSPset.Verbose>0) std::cout<<"\nMax Slope Index "<<ix_maxSlope<<" "<<maxSlope<<" "<<ix_firstJump<<" "<<startSlope<<std::endl;

    if(startSlope!=0) fEta = maxSlope/startSlope;
    if(fEta>100) fEta=100;
    fIota = startNHits;
    fOmega = startSlope;
    fAlpha = fOmega * fZetaNorm;

    if(fFRANSPset.Verbose>0) std::cout<<"\nStartSlope: "<<startSlope<<"  MaxSlope"<<maxSlope<<"  ZetaNorm:"<<fZetaNorm<<std::endl;

    fTau = ix_maxSlope;
    if(fFRANSPset.Verbose>0) std::cout<<"\nIx max slope: "<<ix_maxSlope<<std::endl;
    if(ix_maxSlope>5){
      std::vector<double> CumStart(fZCum.cbegin(), std::next(fZCum.cbegin(), ix_maxSlope+1));
      if(fFRANSPset.Verbose>0) std::cout<< "   l = "<<CumStart.size()<<std::endl;
      fZCumStart = CumStart;
      double maxCumStart =  *std::max_element(fZCumStart.begin(), fZCumStart.end());
      //for(size_t j=0; j<fZCumStart.size(); j++) {std::cout<<j<<" "<<fZCumStart[j]<<std::endl;}
      Rescale(fZCumStart, 1./maxCumStart);

      if(fFRANSPset.Verbose>0) { std::cout<<"\n\nDebug Delta 4: "; for(size_t l=0; l<fZCumStart.size(); l++) {std::cout<<fZCumStart[l]<<" ";} }

      fDelta = 1.*(std::upper_bound(fZCumStart.begin(), fZCumStart.end(), 0.5) - fZCumStart.begin() ) / (ix_maxSlope+1);

      std::vector<double> fFitX;
      for(size_t j=0; j<fZCumStart.size(); j++) {fFitX.push_back(j);}

      std::vector<std::vector<double>> ScoreReg;
      ScoreReg.push_back(fFitX);
      ScoreReg.push_back(fZCumStart);

      Linear_Regression <double, double> Reg_Class;
      Reg_Class.fit(ScoreReg);
      if(isnan(Reg_Class.R2())==false)
        fFitScore = Reg_Class.R2();
      if(fFRANSPset.Verbose>0) std::cout<<" FIT: "<<Reg_Class.R2()<<" "<<Reg_Class.Slope()<<std::endl;



    }
    else{
      fZCumStart.push_back(fZCum[0]);
    }


  }

  if(fFRANSPset.CalculateScore)
    fScore = fTMVAReader.EvaluateMVA( "FRANS BDT" );

  std::cout<<"Alpha: "<<fAlpha<<" Omega: "<<fOmega<<" Iota: "<<fIota<<" Eta: "<<fEta<<" Delta: "<<fDelta<<" FitScore: "<<fFitScore<<" BDTScore:"<<fScore<<std::endl;
}


// ----------- Gauss function -----------
double gaussian(double x, double mu, double sig) {
    return exp(-0.5 * pow((x - mu) / sig, 2)) / (sig * sqrt(2 * M_PI));
}


// ----------- Reset function -----------
void ChargeDensity::Reset() {
    // Reset vectors
    fZ.clear();
    fZ.resize(DefaultMaxZSize, 0);
    fZCounter.clear();
    fZCounter.resize(DefaultMaxZSize, 0);
    fRho.clear();
    fZCum.clear();
    fZCumStart.clear();
    fZCumDer.clear();
    fZCumDerErr.clear();
    fNHits = 0;
    fHitDensity = 0;
    fMeanChi2 = 0;

    // default values
    fDelta = 0.5;
    fEta = 1.;
    fFitScore = 1.;
    fAlpha = 10e3;
    fOmega = 1.;
    fTau = 0.;
    fScore = -1e4;
}


// ----------- Fill function -----------
void ChargeDensity::Fill(std::vector<SHit> hitsVect, SVertex vertex) {

    // Reset all vectors and variables
    Reset();
  
    fVertex = vertex;
    double vTimeTick = vertex.Y();
    int vCh = vertex.X();

    std::vector<SHit> Hits;
    double d_vertexhit = 1e6;
    SHit VertexHit;

    // Find the closest hit to the vertex channel and TT
    for (auto &hit : hitsVect) {
        Hits.push_back(hit);
        double d = GetDistance(hit.X(), hit.Y(), vCh, vTimeTick);
        if (d < d_vertexhit) {
            d_vertexhit = d;
            VertexHit = hit;
        }
    }

    // Redefine start vertex
    vTimeTick = VertexHit.Y();
    vCh = VertexHit.X();

    std::cout << "Refactored vertex: " << vCh << " " << vTimeTick << std::endl;

    

    // Fill the vectors
    int dMax = 0;
    for (auto &hit : Hits) {
        double d = GetDistance(hit.X(), hit.Y(), vCh, vTimeTick);
        int dWires = GetDistanceWires(hit.X(), vCh);
        if (d < DefaultMaxZSize) {
            if(fFRANSPset.UseHitWidth){

              // lower and upper time ticks of the hit
              size_t lower_bin = std::floor(hit.Y()-hit.Width());
              size_t uppper_bin = std::ceil(hit.Y()+hit.Width());

              std::vector<double> hitWeights;
              std::vector<double> hitY;
              for(size_t k = lower_bin; k<=uppper_bin; k++){
                hitWeights.push_back(gaussian(k, hit.Y(), hit.Width()));
                hitY.push_back(k+0.5);
              }

              double sumWeights = 0;
              sumWeights = std::accumulate(hitWeights.begin(), hitWeights.end(), 0.);
              double normFactor = 1./sumWeights;
              
              for(size_t k=0; k<hitWeights.size(); k++){
                double dBin = GetDistance(hit.X(), hitY[k], vCh, vTimeTick);
                if(dBin<DefaultMaxZSize) fZ[static_cast<int>(dBin)] += normFactor * hitWeights[k] * hit.Integral();
              }

              fZCounter[dWires] += 1;
            
            }
            else{
              fZ[static_cast<int>(d)] += hit.Integral();
              fZCounter[dWires] += 1;
            }
            
            fNHits++;
            fHitDensity++;
            fMeanChi2 += hit.Chi2();
            if (d > dMax) dMax = static_cast<int>(d);
        }
    }

    fHitDensity /= dMax;
    fMeanChi2 /= fNHits;
    fZ.resize(dMax);


    if(fZ.size()>1)
      UpdateMetrics();

}


// ----------- Display function -----------
void ChargeDensity::Display(TCanvas *c){

  gStyle->SetPalette(112,0);
  gStyle->SetTitleFont(132, "TXYZ");
  gStyle->SetTitleSize(0.05, "TXYZ");

  gStyle->SetTitleFont(132, "titleFont"); 
  
  // Off stats
  gStyle->SetOptStat(0); 

  //LABELS SIZE AND FONT
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetTitleYOffset (1.4);
  
  c->cd();

  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.66, .98, 0.99);
  pad1->Draw();
  
  TPad *pad2A = new TPad("pad2A", "pad2A", 0., .33, .49, .66);
  pad2A->Draw();
  TPad *pad2B = new TPad("pad2B", "pad2B", 0.49, .33, .98, .66);
  pad2B->Draw();

  TPad *pad3A = new TPad("pad3A", "pad3A", 0., 0., .49, 0.33);
  pad3A->Draw();
  TPad *pad3B = new TPad("pad3B", "pad3B", 0.49, 0., .98, 0.33);
  pad3B->Draw();

 

  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.12);
  pad2A->SetBottomMargin(0.15);
  pad2A->SetLeftMargin(0.2);
  pad2A->SetRightMargin(0.);
  pad2B->SetBottomMargin(0.15);
  pad3A->SetBottomMargin(0.15);
  pad3A->SetLeftMargin(0.2);
  pad3B->SetBottomMargin(0.15);
  pad3B->SetRightMargin(0.);


  gStyle->SetTitleFont(52, "TXYZ");
  gStyle->SetTitleSize(0.05);
  //LABELS SIZE AND FONT
  gStyle->SetLabelFont(60, "XYZ");
  gStyle->SetLabelSize(0.08, "XYZ");
  //AXIS OFFSETS AND SIZES
  gStyle->SetTitleXOffset (0.85);
  gStyle->SetTitleXSize (0.08);
  gStyle->SetTitleYOffset (0.65);
  gStyle->SetTitleYSize (0.08);


  pad1->cd();
  TGraph *gr = new TGraph(fRho.size(), &fRho[0], &fZ[0]);
  gr->SetTitle("");
  gr->GetHistogram()->GetYaxis()->SetTitle("#zeta(#rho) [ADCxTT]");
  gr->GetHistogram()->GetXaxis()->SetTitle("#rho");
  gr->Draw("ALP");

  double score = fTMVAReader.EvaluateMVA( "FRANS BDT" );
  TLegend* leg1 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg1->SetBorderSize(1); leg1->SetTextFont(62); leg1->SetTextSize(0.1);
  std::ostringstream legLabel1; legLabel1 << std::setprecision(2);
  legLabel1 << "Score=" << score;
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.075);
  leg1->AddEntry(gr, legLabel1.str().c_str(), "");
  leg1->Draw("same");


  pad2A->cd();
  TGraph *grCum = new TGraph(fRho.size(), &fRho[0], &fZCum[0]);
  grCum->SetTitle("");
  grCum->GetHistogram()->GetYaxis()->SetTitle("Z(#rho) [ADCxTT]");
  grCum->GetHistogram()->GetXaxis()->SetTitle("#rho");
  grCum->GetHistogram()->GetYaxis()->SetRangeUser(0, 1);
  grCum->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  grCum->Draw("alp");

  
  pad2B->cd();
  TGraph *grCumStart = new TGraph(fZCumStart.size(), &fRho[0], &fZCumStart[0]);
  grCumStart->SetTitle("");
  //grCumStart->GetHistogram()->GetYaxis()->SetTitle("Z(#rho) [ADCxTT]");
  grCumStart->GetHistogram()->GetXaxis()->SetTitle("#rho");
  grCumStart->GetHistogram()->GetYaxis()->SetRangeUser(0, 1);
  grCumStart->SetLineColor(kRed+2);
  grCumStart->Draw("alp");

  TLegend* leg2 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.075);
  leg2->SetHeader("Start");
  std::ostringstream legLabel2;
  legLabel2 << std::setprecision(3);
  legLabel2 << "#Delta=" << fDelta;
  leg2->AddEntry(grCumStart, legLabel2.str().c_str(), "");
  legLabel2.str("");
  legLabel2 << "FitSc="<<fFitScore;
  leg2->AddEntry(grCumStart, legLabel2.str().c_str(), "");
  leg2->Draw("same");


  pad3A->cd();
  TGraph *grDer = new TGraphErrors(fZCumDer.size(), &fRho[0], &fZCumDer[0], 0, &fZCumDerErr[0]);
  grDer->SetTitle("");
  grDer->GetHistogram()->GetYaxis()->SetTitle("Z'(#rho) [ADCxTT]");
  grDer->GetHistogram()->GetXaxis()->SetTitle("#rho");
  grDer->Draw("alp");

  TLegend* leg3 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextFont(62);
  leg3->SetTextSize(0.075);
  std::ostringstream legLabel3;
  legLabel3 << std::setprecision(3);
  legLabel3 << "#eta=" << fEta;
  leg3->AddEntry(grDer, legLabel3.str().c_str(), "");
  legLabel3.str("");
  legLabel3 << " #alpha="<<fAlpha;
  leg3->AddEntry(grDer, legLabel3.str().c_str(), "");
  leg3->Draw("same");

  pad3B->cd();
  TGraph *grCounter = new TGraph(fRho.size(), &fRho[0], &fZCounter[0]);
  grCounter->SetTitle("");
  grCounter->GetHistogram()->GetYaxis()->SetTitle("# hits(#rho)");
  grCounter->GetHistogram()->GetXaxis()->SetTitle("#rho");
  grCounter->Draw("ALP");

  TLegend* leg4 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextFont(62);
  leg4->SetTextSize(0.075);
  std::ostringstream legLabel4;
  legLabel4 << std::setprecision(3);
  legLabel4 << "#iota=" << fIota;
  leg4->AddEntry(grDer, legLabel4.str().c_str(), "");
  leg4->Draw("same");


  c->cd();
  c->Update();

  return;
}

FRANSObj ChargeDensity::GetFRANSResult(){
  FRANSObj FRANSResult(fView, fScore, fHitDensity, fMeanChi2,
                        fDelta, fEta, fFitScore,
                        fAlpha, fOmega, fTau, fIota, 
                        0, 0, 0, 0);
  return FRANSResult;                       
}
