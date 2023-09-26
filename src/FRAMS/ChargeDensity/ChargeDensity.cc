#include "ChargeDensity.hh"
#include "LeastSquares.hh"
#include <TROOT.h>
#include <TStyle.h>
#include <math.h>
#include "TLegend.h"
#include <sstream>
#include <iomanip>

void ChargeDensity::HelloWorld(){
  std::cout<<"Helllo World!\n";
}

double ChargeDensity::GetDistance(int x, double y, int x0, double y0){
  //std::cout<<x-x0<<" "<<y-y0<<std::endl;
  return std::sqrt( std::pow( ((double)(x-x0))/fNWirePack, 2) + std::pow( (y-y0)/fNDriftPack, 2) );
}

ChargeDensity::ChargeDensity(ChargeDensityConf::FRAMSPsetType const& config, unsigned int view, unsigned int tpc)
  : //fPset(config),
  fApplyRawSmoothing( config.applyRawSmoothing() ),
  fApplySmoothing( config.applySmoothing() ),
  fApplyCumulativeSmoothing( config.applyCumulativeSmoothing() ),
  fNDriftPack( config.nDriftPack() ),
  fNWirePack( config.nWirePack() ),
  fExpoAvSmoothPar( config.expoAvSmoothPar() ),
  fUnAvNeighbours( config.unAvNeighbours() ),
  fCumulativeCut( config.cumulativeCut() ),
  fSlidingWindowN( config.slidingWindowN() ),
  fMaxRadius( config.maxRadius() ),
  fNSamplesBeginSlope( config.nSamplesBeginSlope() ),
  fDebugMode( config.debugMode() ),
  fCalculateScore( config.calculateScore() ),
  fTMVAFilename( config.tMVAFilename() ),
  fView(view),
  fTPC(tpc)
  {

  fNormUnAvSmooth=1./(2*fUnAvNeighbours+1);

  // default values
  fDelta = 0.5;
  fEta = 1.;
  fFitScore = 1.;
  fAlpha = 10e3;
  fOmega = 1.;
  fTau = 0.;
  fScore = -1e4;

  // Right now we only use the collection view
  if(fCalculateScore){
    fTMVAReader.AddVariable( "Alpha_C", &fAlpha );
    fTMVAReader.AddVariable( "Eta_C", &fEta );
    fTMVAReader.AddVariable( "Delta_C", &fDelta);
    fTMVAReader.AddVariable( "FitScore_C", &fFitScore );

    fTMVAReader.AddSpectator( "Gap", &fGap );
    fTMVAReader.AddSpectator( "ProtonKE", &fProtonKE );
    fTMVAReader.AddSpectator( "PionKE", &fPionKE );

    std::string file_name;
    cet::search_path sp("FW_SEARCH_PATH");
    if ( !sp.find_file(fTMVAFilename, file_name) )
     throw cet::exception("FRMASFilter module") << "BDT file " <<
         fTMVAFilename << " not found in FW_SEARCH_PATH\n";
    std::cout<<" BDT weights file name: "<<file_name<<std::endl;

    fTMVAReader.BookMVA( "FRAMS BDT",  file_name.c_str()  );
  }
}


void ChargeDensity::Rescale(std::vector<double>& wf, double kappa){
  std::transform(wf.begin(), wf.end(), wf.begin(), [&kappa](auto& z){return z*kappa;});
}


void ChargeDensity::ApplyExpoAvSmoothing(std::vector<double>& wf){
  std::transform (std::next(wf.begin(), 1), wf.end(), wf.begin(), std::next(wf.begin(), 1),
    [&](double _x, double _y) { return  fExpoAvSmoothPar*_x+ (1. - fExpoAvSmoothPar)*_y; }  );
}


void ChargeDensity::ApplyUnAvSmoothing(std::vector<double>& wf){
  std::vector<double> wf_aux(wf.begin(), wf.end());
  for(size_t bin=fUnAvNeighbours; bin<wf.size()-fUnAvNeighbours; bin++){
    double sum=0.;
    for(size_t nbin=bin-fUnAvNeighbours; nbin<=bin+fUnAvNeighbours; nbin++)
      sum+=wf_aux[nbin];
    wf[bin]=sum*fNormUnAvSmooth;
  }
}


void ChargeDensity::SlidingWindow(std::vector<double>& wf){

  fZCumDer.clear(); fZCumDer.reserve(fZCum.size());
  fZCumDerErr.clear(); fZCumDerErr.reserve(fZCum.size());

  std::vector<float> x; x.reserve(fSlidingWindowN);
  for(int k=0; k<fSlidingWindowN; k++) x.push_back(k);

  std::vector<std::vector<float>> Z;
  Z.push_back(x); Z.push_back(x);

  Linear_Regression <float,float> Reg_Class;

  for(size_t k=0; k<wf.size()-fSlidingWindowN; k++){

    std::vector<float> y;
    for(size_t j=k; j<k+fSlidingWindowN; j++){
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

  size_t start = std::upper_bound(fZCum.begin(), fZCum.end(), fCumulativeCut*fMaxCumulative) - fZCum.begin();

  fZ.resize(start);
  fZCum.resize(start);
}


void ChargeDensity::UpdateMetrics(){

  if(fDebugMode) std::cout<<"Updating metrics...size is: "<<fZ.size()<<std::endl;

  if(fDebugMode) { std::cout<<"Debug 0: "; for(size_t l=0; l<10; l++) {std::cout<<fZ[l]<<" ";} }

  //for(size_t k=0; k<fZ.size(); k++){std::cout<<k<<" "<<fZ[k]<<std::endl;}

  fZetaNorm = *std::max_element(fZ.begin(), fZ.end());

  if(fDebugMode) std::cout<<"ZetaNorm: "<<fZetaNorm<<std::endl;
  if(fApplyRawSmoothing) {
    // keep first point without smoothing
    std::vector<double> Zsmoothed(std::next(fZ.cbegin(), 1), fZ.cend());
    ApplyExpoAvSmoothing(Zsmoothed); //!!!!! Check this
    Zsmoothed.insert(Zsmoothed.begin(), fZ[0]);
    fZ.clear();
    fZ.insert( fZ.begin(), Zsmoothed.begin(), Zsmoothed.end());
  }

  double ZetaNormI = 1./( *std::max_element(fZ.begin(), fZ.end()) );
  if(fDebugMode) std::cout<<"Finished raw smooth \n";
  Rescale(fZ, ZetaNormI);

  if(fDebugMode) {std::cout<<"\n\nDebug 1: "; for(size_t l=0; l<10; l++) {std::cout<<fZ[l]<<" ";} }

  if(fDebugMode) std::cout<<"Finished norm \n";
  //fill cumulative and resize vector keeping only ROI
  FillCumulative();

  if(fDebugMode) { std::cout<<"\n\nDebug Cum 2: "; for(size_t l=0; l<10; l++) {std::cout<<fZCum[l]<<" ";} }

  if((int)fZ.size()>fSlidingWindowN){
    SlidingWindow(fZCum);

    if(fDebugMode) std::cout<<"Finished SW \n";

    if(fApplySmoothing){
      ApplyExpoAvSmoothing(fZCumDer);
      ApplyUnAvSmoothing(fZCumDer);
      if(fDebugMode) std::cout<<"Finished SW smooth \n";
    }

    Rescale(fZCum, 1./fMaxCumulative);

    fRho.clear(); fRho.reserve(fZ.size());
    for(size_t k=0; k<fZ.size(); k++){
      fRho.push_back(k);
    }



    //--- Now we start with the parameter estimation
    // first get ROIs
    int L = std::min(fMaxRadius, (int)fZCumDer.size());
    std::vector<double> Slopes(fZCumDer.cbegin(), std::next(fZCumDer.cbegin(), L));
    std::vector<double> SlopesErr(fZCumDerErr.cbegin(), std::next(fZCumDerErr.cbegin(), L));

    if(fDebugMode) { std::cout<<"\n\nDebug Slopes 3: "; for(size_t l=0; l<10; l++) {std::cout<<Slopes[l]<<" ";} }

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

    if(fDebugMode) std::cout<<"\nFirtJump @ "<<ix_firstJump<<" MaxSlope @ "<<ix_maxSlope<<std::endl;
    double startSlope = 1.;
    if(fDebugMode) std::cout<<" Start Slope: ";
    if(ix_firstJump>0){
      startSlope = std::accumulate(Slopes.cbegin(), std::next(Slopes.cbegin(), ix_firstJump), 0.)/(ix_firstJump);
      for(int k=0; k<ix_firstJump+1; k++) std::cout<<k<<":"<<Slopes[k]<<"  ";
    }
    else{
      startSlope = std::accumulate(Slopes.cbegin(), std::next(Slopes.cbegin(), fNSamplesBeginSlope), 0.)/(fNSamplesBeginSlope);
      for(int k=0; k<fNSamplesBeginSlope; k++) std::cout<<k<<":"<<Slopes[k]<<"  ";
    }

    if(fDebugMode) std::cout<<" Max Slope Index "<<ix_maxSlope<<" "<<maxSlope<<" "<<ix_firstJump<<" "<<startSlope<<std::endl;

    if(startSlope!=0) fEta = maxSlope/startSlope;
    if(fEta>100) fEta=100;
    fOmega = startSlope;
    fAlpha = fOmega * fZetaNorm;

    if(fDebugMode) std::cout<<" StartSlope: "<<startSlope<<"  MaxSlope"<<maxSlope<<"  ZetaNorm:"<<fZetaNorm<<std::endl;

    fTau = ix_maxSlope;
    if(fDebugMode) std::cout<<" Ix max slope: "<<ix_maxSlope<<std::endl;
    if(ix_maxSlope>5){
      std::vector<double> CumStart(fZCum.cbegin(), std::next(fZCum.cbegin(), ix_maxSlope+1));
      if(fDebugMode) std::cout<< "   l = "<<CumStart.size()<<std::endl;
      fZCumStart = CumStart;
      double maxCumStart =  *std::max_element(fZCumStart.begin(), fZCumStart.end());
      //for(size_t j=0; j<fZCumStart.size(); j++) {std::cout<<j<<" "<<fZCumStart[j]<<std::endl;}
      Rescale(fZCumStart, 1./maxCumStart);

      if(fDebugMode) { std::cout<<"\n\nDebug Delta 4: "; for(size_t l=0; l<fZCumStart.size(); l++) {std::cout<<fZCumStart[l]<<" ";} }

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
      if(fDebugMode) std::cout<<" FIT: "<<Reg_Class.R2()<<" "<<Reg_Class.Slope()<<std::endl;



    }
    else{
      fZCumStart.push_back(fZCum[0]);
    }


  }

  if(fCalculateScore)
    fScore = fTMVAReader.EvaluateMVA( "FRAMS BDT" );

  std::cout<<"Alpha: "<<fAlpha<<" Omega: "<<fOmega<<" Eta: "<<fEta<<" Delta: "<<fDelta<<" FitScore: "<<fFitScore<<" BDTScore:"<<fScore<<std::endl;
}


void ChargeDensity::Fill(std::vector<art::Ptr<recob::Hit>> hitsVect, VertexView vertex){

  fVertex = vertex;

  double vTimeTick = vertex.TT;
  int vCh;
  if(fView==0) vCh = vertex.U;
  else if(fView==1) vCh = vertex.V;
  else vCh = vertex.C;


  std::vector<art::Ptr<recob::Hit>> Hits;
  double d_vertexhit=1e6;
  recob::Hit VertexHit;



  // Get growing start point (closest hit to the vertex channel and TT)
  for(auto &hit: hitsVect){
    if(hit->WireID().TPC==fTPC && hit->WireID().Plane==fView){

      //std::cout<<hit->View()<<" "<<hit->WireID().TPC<<std::endl;
      Hits.push_back(hit);

      double d = GetDistance(hit->Channel(), hit->PeakTime(), vCh, vTimeTick);
      //std::cout<<"d="<<d<<" "<<hit->Channel()<<" "<<hit->Integral()<<std::endl;
      if(d<d_vertexhit){
        d_vertexhit = d;
        VertexHit = *hit;
      }

    }

  }

  //Redefine start vertex
  vTimeTick = VertexHit.PeakTime();
  vCh = VertexHit.Channel();

  std::cout<<" Refactored vertex: "<<vCh<<" "<<vTimeTick<<std::endl;

  fZ.clear();
  fZ.resize(DefaultMaxZSize, 0);
  int dMax=0;
  for(auto &hit: Hits){
    double d = GetDistance(hit->Channel(), hit->PeakTime(), vCh, vTimeTick);
    if(d<DefaultMaxZSize){
      fZ[(int)d]+=hit->Integral();
      if(d>dMax) dMax = (int) d;
    }
  }

  fZ.resize(dMax);

  if(fZ.size()>1)
    UpdateMetrics();

}


void ChargeDensity::Save2ROOT(art::TFileDirectory tfdir, std::string name){
  TCanvas * c = tfdir.makeAndRegister< TCanvas >(("c_"+name+"_vw"+std::to_string(fView)).c_str(),"Charge Density");
  c->Divide(1,3);

  gStyle->SetTitleFont(52, "TXYZ");
  gStyle->SetTitleSize(0.05);
  //LABELS SIZE AND FONT
  gStyle->SetLabelFont(60, "XYZ");
  gStyle->SetLabelSize(0.08, "XYZ");
  //AXIS OFFSETS AND SIZES
  gStyle->SetTitleXOffset (0.60);
  gStyle->SetTitleXSize (0.08);
  gStyle->SetTitleYOffset (0.65);
  gStyle->SetTitleYSize (0.08);

  c->cd(1);
  TGraph *gr = new TGraph(fRho.size(), &fRho[0], &fZ[0]);
  gr->SetTitle("");
  gr->GetHistogram()->GetYaxis()->SetTitle("#zeta(#rho) [ADCxTT]");
  gr->GetHistogram()->GetXaxis()->SetTitle("#rho");
  gr->Draw("ALP");

  double score = fTMVAReader.EvaluateMVA( "FRAMS BDT" );
  TLegend* leg1 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg1->SetBorderSize(1); leg1->SetTextFont(62); leg1->SetTextSize(0.1);
  std::ostringstream legLabel1; legLabel1 << std::setprecision(3);
  legLabel1 << "Score=" << score;
  leg1->AddEntry(gr, legLabel1.str().c_str(), "");
  leg1->Draw("same");


  c->cd(2);
  TGraph *grCum = new TGraph(fRho.size(), &fRho[0], &fZCum[0]);
  TGraph *grCumStart = new TGraph(fZCumStart.size(), &fRho[0], &fZCumStart[0]);
  grCum->SetTitle("");
  grCum->GetHistogram()->GetYaxis()->SetTitle("Z(#rho) [ADCxTT]");
  grCum->GetHistogram()->GetXaxis()->SetTitle("#rho");
  grCum->GetHistogram()->GetYaxis()->SetRangeUser(0, 1);
  grCumStart->SetLineColor(kRed);
  grCum->Draw("alp");
  grCumStart->Draw("lp same");

  TLegend* leg2 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg2->SetBorderSize(1); leg2->SetTextFont(62); leg2->SetTextSize(0.1);
  std::ostringstream legLabel2; legLabel2 << std::setprecision(3);
  legLabel2 << "#Delta=" << fDelta<< " FitScore="<<fFitScore;
  leg2->AddEntry(grCum, legLabel2.str().c_str(), "");
  leg2->Draw("same");


  c->cd(3);
  TGraph *grDer = new TGraphErrors(fZCumDer.size(), &fRho[0], &fZCumDer[0], 0, &fZCumDerErr[0]);
  grDer->SetTitle("");
  grDer->GetHistogram()->GetYaxis()->SetTitle("Z'(#rho) [ADCxTT]");
  grDer->GetHistogram()->GetXaxis()->SetTitle("#rho");
  grDer->Draw("alp");

  TLegend* leg3 = new TLegend(0.5, 0.60, 0.85, 0.85);
  leg3->SetBorderSize(1); leg3->SetTextFont(62); leg3->SetTextSize(0.1);
  std::ostringstream legLabel3; legLabel3 << std::setprecision(3);
  legLabel3 << "#eta=" << fEta<< " #alpha="<<fAlpha;
  leg3->AddEntry(grDer, legLabel3.str().c_str(), "");
  leg3->Draw("same");




  c->cd();
  c->Update();
}
