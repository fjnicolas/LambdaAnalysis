////////////////////////////////////////////////////////////////////////
// Class:       FRAMSAna
// Plugin Type: Producer (art v3_05_01)
// File:        FRAMSAna_module.cc
//
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"


#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"



#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//#include "sbndcode/FRAMS/FRAMSObj/FRAMSObj.h"
#include "sbnobj/SBND/FRAMSObj/FRAMSObj.h"

#include "sbndcode/FRAMS/ChargeDensity/ChargeDensity.hh"
#include "sbndcode/FRAMS/ChargeDensity/ChargeDensityAlgConf.hh"
#include "sbndcode/HyperonAnalyzer/LambdaTruthManager/LambdaTruthManager.hh"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"

#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>

#define fXFidCut 190
#define fYFidCut 190
#define fZFidCut1 10
#define fZFidCut2 490

#define fDefaulNeutrinoID 99999

namespace sbnd {
  class FRAMSAna
  ;
}


class sbnd::FRAMSAna : public art::EDProducer {
public:

  struct Config {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;

      fhicl::Atom<std::string> MCTruthLabel {Name("MCTruthLabel")};
      fhicl::Atom<std::string> HitLabel {Name("HitLabel")};
      fhicl::Atom<std::string> VertexLabel {Name("VertexLabel")};
      fhicl::Atom<std::string> RecoLabel {Name("RecoLabel")};
      fhicl::Atom<bool> ApplyFiducialCut {Name("ApplyFiducialCut")};
      fhicl::Sequence<int> Views {Name("Views")};
      fhicl::Atom<float> NuScoreCut {Name("NuScoreCut")};
      fhicl::Atom<float> TrueVertexBallRadius {Name("TrueVertexBallRadius")};
      fhicl::Atom<bool> UseTrueVertex {Name("UseTrueVertex")};
      fhicl::Atom<bool> FillLambdaTrue {Name("FillLambdaTrue")};
      fhicl::Atom<bool> SaveChargeProfiles {Name("SaveChargeProfiles")};
      fhicl::Atom<std::string> TreeName {Name("TreeName")};

      fhicl::Table<ChargeDensityConf::Config> ChargeDensityAlg{Name("ChargeDensityAlg")};
    }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;

  explicit FRAMSAna(Parameters const& config);

  // Plugins should not be copied or assigned.
  FRAMSAna(FRAMSAna const&) = delete;
  FRAMSAna(FRAMSAna&&) = delete;
  FRAMSAna & operator=(FRAMSAna const&) = delete;
  FRAMSAna & operator=(FRAMSAna &&) = delete;

  // Required functions.
  void produce(art::Event &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  void resetVars();
  int VertexToDriftTick(double vt, double vx);
  float GetSliceNuScore(art::Ptr<recob::PFParticle> &nupfp,
    std::vector<art::Ptr<recob::PFParticle>> pfpVect,
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns);
  art::Ptr<recob::PFParticle> GetNeutrinoPFP(std::vector<art::Ptr<recob::PFParticle>> pfpVect);
  double GetDistance(double p1[], double p2[]);

  // handles
  art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
  art::Handle< std::vector<simb::MCParticle> > mcparticleHandle;
  std::vector<art::Ptr<simb::MCTruth>> mctruthVect;
  std::vector<art::Ptr<simb::MCParticle>> mcpVect;

  //............................Read Recob Slice
  ::art::Handle<std::vector<recob::Slice>> sliceHandle;
  //............................Read PFPs
  ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVect;
  //............................Read vertex
  ::art::Handle<std::vector<recob::Vertex>> vertexHandle;
  std::vector<art::Ptr<recob::Vertex>> vertexVect;
  //............................Read hits
  ::art::Handle<std::vector<recob::Hit>> hitsHandle;
  std::vector<art::Ptr<recob::Hit>> hitsVect;

  // configuration parameters
  std::string fMCTruthLabel;
  std::string fHitLabel;
  std::string fVertexLabel;
  std::string fRecoLabel;
  bool fApplyFiducialCut;
  std::string fTreeName;
  std::vector<int> fViews;
  float fNuScoreCut;
  float fTrueVertexBallRadius;
  bool fUseTrueVertex;
  bool fFillLambdaTrue;
  bool fSaveChargeProfiles;
  ChargeDensityConf::Config fChargeDensityAlg;

  //True variables
  std::vector<int> fTruePrimariesPDG;
  std::vector<double> fTruePrimariesE;

  int fNAnalyzedEvents;

  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  unsigned int fNChannels;
  unsigned int fReadoutWindow;
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;

  VertexView fTrueVertex;

  //TTree variables
  TTree* fTree;
  int fEventID, fRunID, fSubRunID;

  int fNNuInts;
  int fIntMode;
  int fIntType;

  double fGap;
  double fLambdaKE;
  double fProtonKE;
  double fPionKE;

  double fDelta_U, fDelta_V, fDelta_C;
  double fEta_U, fEta_V, fEta_C;
  double fFitScore_U, fFitScore_V, fFitScore_C;
  double fAlpha_U, fAlpha_V, fAlpha_C;
  double fOmega_U, fOmega_V, fOmega_C;
  double fTau_U, fTau_V, fTau_C;
  double fIota_U, fIota_V, fIota_C;
  double fScore_U, fScore_V, fScore_C;


  //Load TFileService serrvice
  art::ServiceHandle<art::TFileService> tfs;
};


sbnd::FRAMSAna::FRAMSAna(Parameters const& config)
  : EDProducer{config},
  fMCTruthLabel( config().MCTruthLabel() ),
  fHitLabel( config().HitLabel() ),
  fVertexLabel( config().VertexLabel() ),
  fRecoLabel( config().RecoLabel() ),
  fApplyFiducialCut( config().ApplyFiducialCut() ),
  fTreeName( config().TreeName() ),
  fViews( config().Views() ),
  fNuScoreCut( config().NuScoreCut() ),
  fTrueVertexBallRadius( config().TrueVertexBallRadius() ),
  fUseTrueVertex( config().UseTrueVertex() ),
  fFillLambdaTrue( config().FillLambdaTrue() ),
  fSaveChargeProfiles( config().SaveChargeProfiles() ),
  fChargeDensityAlg( config().ChargeDensityAlg()),
  fNChannels(fGeom->Nchannels())
  // More initializers here.
{

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

  fTriggerOffsetTPC = clockData.TriggerOffsetTPC(); //in us
  fTickPeriodTPC = clockData.TPCClock().TickPeriod(); //in us
  fReadoutWindow = detProp.ReadOutWindowSize();
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  constexpr geo::TPCID tpcid{0, 0};
  fWirePlanePosition = std::abs( fGeom->Plane(geo::PlaneID{tpcid, 1}).GetCenter().X() );


  std::cout<<"  - Read TPC clocks...  ReadOutWindowSize: "<<fReadoutWindow<<"  TriggerOffsetTPC: "<<fTriggerOffsetTPC;
  std::cout<<"  TickPeriodTPC: "<<fTickPeriodTPC<<std::endl;
  std::cout<<"  - Drift Velocity: "<<fDriftVelocity<<"  WirePlanePosition: "<<fWirePlanePosition<<std::endl;

  produces< std::vector<sbnd::FRAMSObj> >();
  produces<art::Assns<recob::Slice, sbnd::FRAMSObj>>();
}


void sbnd::FRAMSAna::produce(art::Event &e)
{
  std::cout<<"Running FRAMSAna---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<"\n";

  //............................Event General Info
  fNAnalyzedEvents++;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fEventID = e.id().event();

  std::unique_ptr<  std::vector<sbnd::FRAMSObj> > framsObj_V( new  std::vector<sbnd::FRAMSObj> );
  std::unique_ptr<  art::Assns<recob::Slice, sbnd::FRAMSObj> > framsObjAssns_V( new art::Assns<recob::Slice, sbnd::FRAMSObj> );
  //std::unique_ptr<  sbnd::FRAMSObj > framsObj_A(new  sbnd::FRAMSObj );


  //Reset tree variables
  resetVars();

  //............................Read Truth Objects
  e.getByLabel(fMCTruthLabel, mctruthHandle);
  art::fill_ptr_vector(mctruthVect, mctruthHandle);
  e.getByLabel("largeant", mcparticleHandle);
  art::fill_ptr_vector(mcpVect, mcparticleHandle);
  std::cout<<" --- Saving MCTruth\n";
  std::cout<<"     --- Number of MCTruths: "<<mctruthVect.size()<<"\n";

  //Reconstructed fVertex
  double fTrueVx=-1e3, fTrueVy=-1e3, fTrueVz=-1e3, fTrueVt=-1e6;

  for (auto const& truth : mctruthVect) {

    for (int p = 0; p < truth->NParticles(); p++){
      simb::MCParticle const& TruePart = truth->GetParticle(p);

      if( TruePart.StatusCode()==1 ){
        fTruePrimariesPDG.push_back( TruePart.PdgCode() );
        fTruePrimariesE.push_back( TruePart.E() );
      }

      if( TruePart.Mother()==-1 && ( abs(TruePart.PdgCode())==12 || abs(TruePart.PdgCode())==14 ) ){

        if(fApplyFiducialCut && std::abs(TruePart.EndX())<fXFidCut && std::abs(TruePart.EndY())<fYFidCut && TruePart.EndZ()>fZFidCut1 && TruePart.EndZ()<fZFidCut2){
          fNNuInts++;

          fTrueVx=TruePart.EndX();
          fTrueVy=TruePart.EndY();
          fTrueVz=TruePart.EndZ();
          fTrueVt=TruePart.T();

          fIntMode = truth->GetNeutrino().Mode();
          fIntType = truth->GetNeutrino().InteractionType();
        }

        // printout
        std::cout << "     **** Mode=" << truth->GetNeutrino().Mode() <<"  IntType="<<truth->GetNeutrino().InteractionType()
        <<" CCNC=" << truth->GetNeutrino().CCNC()<<" Origin: " << truth->Origin() << " Target=" << truth->GetNeutrino().Target() << std::endl;
        std::cout<<"  -- Vertex: "<<fTrueVx<<" "<<fTrueVy<<" "<<fTrueVz<<std::endl;
        //std::cout<<"   - VertexWire: "<<fTrueVU<<" "<<fTrueVV<<" "<<fTrueVC<<" "<<fTrueVTimeTick<<std::endl;
      }
    }

  }

  std::cout<<" Number of nus in TPC fiducial volume: "<<fNNuInts<<std::endl;

  if(fFillLambdaTrue){
    LambdaTruthManager lambdaMgr(mctruthVect, mcpVect);
    fGap = lambdaMgr.Gap();
    fLambdaKE = lambdaMgr.LambdaKE();
    fProtonKE = lambdaMgr.ProtonKE();
    fPionKE = lambdaMgr.PionKE();

    std::cout<<"Gap is "<<fGap<<std::endl;
  }
  else{
    fGap = -1;
    fLambdaKE = -1;
    fProtonKE = -1;
    fPionKE = -1;
  }



  // Only store vertex and analyze event if there's a single nu interaciton in the TPC
  // This need to be updated/improved to include events with multiple neutrino events

  // Store the fVertex view
  if(fNNuInts==1){
    geo::Point_t p={fTrueVx, fTrueVy, fTrueVz};
    std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(p))<<" "<<fGeom->FindTPCAtPosition(p).TPC<<std::endl;
    int fTrueVU, fTrueVV, fTrueVC, fTrueVTimeTick;
    if( fGeom->HasTPC(fGeom->FindTPCAtPosition(p)) ){
      unsigned int tpcID=fGeom->FindTPCAtPosition(p).TPC;
      fTrueVU=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 0));
      fTrueVV=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 1));
      fTrueVC=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 2));
      fTrueVTimeTick=VertexToDriftTick(fTrueVt, fTrueVx);
    }
    else{
      fTrueVU=-1;
      fTrueVV=-1;
      fTrueVC=-1;
      fTrueVTimeTick=-1;
    }
    fTrueVertex.reco = false;
    fTrueVertex.x = fTrueVx; fTrueVertex.y = fTrueVy; fTrueVertex.z = fTrueVz; fTrueVertex.t = fTrueVt;
    fTrueVertex.U = fTrueVU; fTrueVertex.V = fTrueVV; fTrueVertex.C = fTrueVC; fTrueVertex.TT = fTrueVTimeTick;
    fTrueVertex.status = fTrueVU!=-1;
  }

  bool inFidVolume=true;
  if(fApplyFiducialCut && (std::abs(fTrueVertex.x)>fXFidCut ||
    std::abs(fTrueVertex.y)>fYFidCut || fTrueVertex.z<fZFidCut1 || fTrueVertex.x>fZFidCut2) )
  {
    inFidVolume=false;
  }

  if(inFidVolume==true){

    e.getByLabel(fRecoLabel, sliceHandle);
    e.getByLabel(fRecoLabel, pfpHandle);
    e.getByLabel(fRecoLabel, vertexHandle);
    e.getByLabel(fHitLabel, hitsHandle);

    //Slice association for Hits
    art::FindManyP<recob::Hit> slice_hit_assns (sliceHandle, e, fRecoLabel);
    //Slice association for PFParticles
    art::FindManyP<recob::PFParticle> slice_pfp_assns (sliceHandle, e, fRecoLabel);
    //PFParticle association for Vertex
    art::FindManyP<recob::Vertex> pfp_vertex_assns (pfpHandle, e, fRecoLabel);
    //PDParticle associated metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns (pfpHandle, e, fRecoLabel);

    std::vector< art::Ptr<recob::Slice> > sliceVect;
    art::fill_ptr_vector(sliceVect, sliceHandle);
    std::cout<<" Number of slices: "<<sliceVect.size()<<std::endl;
    double xyzVertex[3]={fTrueVertex.x, fTrueVertex.y, fTrueVertex.z};
    //Store candidate slices indexes
    std::vector<size_t> sliceCandidatesV;
    std::vector<VertexView> sliceCandidatesVertex;

    for(size_t ix=0; ix<sliceVect.size(); ix++){
      auto & slice = sliceVect[ix];
      pfpVect = slice_pfp_assns.at(slice.key());
      std::cout<<"   *** PFParticle size:"<<pfpVect.size()<<std::endl;

      art::Ptr<recob::PFParticle> nupfp;
      float nu_score=GetSliceNuScore(nupfp, pfpVect, pfp_pfpmeta_assns);

      std::cout<<"Nu Score: "<<nu_score<<std::endl;

      if(nu_score<fNuScoreCut) continue;

      vertexVect = pfp_vertex_assns.at(nupfp.key());

      double slc_vertex[3];
      for(const art::Ptr<recob::Vertex> &ver : vertexVect){
        ver->XYZ(slc_vertex);
        //fpfpVertexPosition.push_back( xyz_vec );
        double chi2=ver->chi2(), chi2ndof=ver->chi2PerNdof();
        std::cout<<"  --VERTEX  ID="<<ver->ID()<<"  x,y,z="<<slc_vertex[0]<<","<<slc_vertex[1]<<","<<slc_vertex[2];
        std::cout<<" Chi2="<<chi2<<" Chi2/DoF="<<chi2ndof<<" Status:"<<ver->status()<<",\n";
      }

      // if using true vertex, only use slice if it's in radius
      if(fUseTrueVertex){
        double d = GetDistance( xyzVertex, slc_vertex);
        if(d<fTrueVertexBallRadius){
          sliceCandidatesV.push_back(ix);
          sliceCandidatesVertex.push_back(fTrueVertex);
        }
      }
      else{

        geo::Point_t p={slc_vertex[0], slc_vertex[1], slc_vertex[2]};

        if(fApplyFiducialCut && std::abs(p.X())<fXFidCut && std::abs(p.Y())<fYFidCut && p.Z()>fZFidCut1 && p.Z()<fZFidCut2){
          VertexView recoVertex;
          std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(p))<<" "<<fGeom->FindTPCAtPosition(p).TPC<<std::endl;
          int fRecoVU, fRecoVV, fRecoVC, fRecoVTimeTick;
          if( fGeom->HasTPC(fGeom->FindTPCAtPosition(p)) ){
            unsigned int tpcID=fGeom->FindTPCAtPosition(p).TPC;
            fRecoVU=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 0));
            fRecoVV=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 1));
            fRecoVC=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 2));
            fRecoVTimeTick=VertexToDriftTick(0., p.X());
          }
          else{
            fRecoVU=-1;
            fRecoVV=-1;
            fRecoVC=-1;
            fRecoVTimeTick=-1;
          }
          recoVertex.reco = true;
          recoVertex.x = slc_vertex[0]; recoVertex.y = slc_vertex[1]; recoVertex.z = slc_vertex[2]; recoVertex.t = 0;
          recoVertex.U = fRecoVU; recoVertex.V = fRecoVV; recoVertex.C = fRecoVC; recoVertex.TT = fRecoVTimeTick;
          recoVertex.status = fRecoVU!=-1;

          sliceCandidatesVertex.push_back(recoVertex);
          sliceCandidatesV.push_back(ix);

        }

      }
    }

    std::cout<<" Number of candidate slices: "<<sliceCandidatesV.size()<<std::endl;
    for(size_t ix=0; ix<sliceCandidatesV.size(); ix++){
      std::cout<<"------ Candidate slice: "<<sliceCandidatesV[ix]<<std::endl;

      auto & slice = sliceVect[ sliceCandidatesV[ix] ];

      hitsVect = slice_hit_assns.at(slice.key());
      std::cout<<" --- Read recob::Hit... we have "<<hitsVect.size()<<" in total\n";

      for(size_t k=0; k<fViews.size(); k++){
        std::cout<<"Setting view "<<fViews[k]<<std::endl;
        int tpc=0;

        VertexView v;
        if(fUseTrueVertex){
          // If using true vertex, use the MCTruth vertex previously stored
          std::cout<<"Using true vertex\n";
          v=fTrueVertex;
          if(fTrueVertex.x>0) tpc=1;
        }
        else{
          // If using reco vertex, use the Slice vertex
          std::cout<<"Using reco vertex\n";
          v=sliceCandidatesVertex[ix];
          if(v.x>0) tpc=1;
        }

        ChargeDensity chargeDensity(fChargeDensityAlg, fViews[k], tpc);
        chargeDensity.Fill(hitsVect, v);

        // saving metrics to root file
        if(fSaveChargeProfiles){
          std::string dir_name = std::to_string(fRunID)+"_"+std::to_string(fSubRunID)+"_"+std::to_string(fEventID);
          art::TFileDirectory tfdir  = tfs->mkdir( dir_name );
          chargeDensity.Save2ROOT(tfdir, dir_name);
        }

        if(fViews[k]==0){
          fDelta_U = chargeDensity.Delta();
          fEta_U = chargeDensity.Eta();
          fFitScore_U = chargeDensity.FitScore();
          fAlpha_U = chargeDensity.Alpha();
          fOmega_U = chargeDensity.Omega();
          fTau_U = chargeDensity.Tau();
          fIota_U = chargeDensity.Iota();
          fScore_U = chargeDensity.Score();
        }
        else if(fViews[k]==1){
          fDelta_V = chargeDensity.Delta();
          fEta_V = chargeDensity.Eta();
          fFitScore_V = chargeDensity.FitScore();
          fAlpha_V = chargeDensity.Alpha();
          fOmega_V = chargeDensity.Omega();
          fTau_V = chargeDensity.Tau();
          fIota_V = chargeDensity.Iota();
          fScore_V = chargeDensity.Score();
        }
        else if(fViews[k]==2){
          fDelta_C = chargeDensity.Delta();
          fEta_C = chargeDensity.Eta();
          fFitScore_C = chargeDensity.FitScore();
          fAlpha_C = chargeDensity.Alpha();
          fOmega_C = chargeDensity.Omega();
          fTau_C = chargeDensity.Tau();
          fIota_C = chargeDensity.Iota();
          fScore_C = chargeDensity.Score();
        }

      }

      bool is_signal = fGap>0 ? 1:0;
      sbnd::FRAMSObj framsObj(
        fDelta_U, fEta_U, fFitScore_U, fAlpha_U, fOmega_U, fTau_U, fIota_U,
        fDelta_V, fEta_V, fFitScore_V, fAlpha_V, fOmega_V, fTau_V, fIota_V,
        fDelta_C, fEta_C, fFitScore_C, fAlpha_C, fOmega_C, fTau_C, fIota_C,
        fGap, fProtonKE, fPionKE, fScore_U, fScore_V, fScore_C, is_signal);

      // Add objects and its association to the art event
      framsObj_V->emplace_back(std::move(framsObj));
      util::CreateAssn(*this, e, *framsObj_V, slice, *framsObjAssns_V);

      fTree->Fill();
    }
  }

  e.put( std::move(framsObj_V) );
  e.put( std::move(framsObjAssns_V) );
  //e.put( std::move(framsObj_A) );


}


int sbnd::FRAMSAna::VertexToDriftTick(double vt, double vx){
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}


void sbnd::FRAMSAna::resetVars()
{
  mctruthVect.clear();
  mcpVect.clear();

  pfpVect.clear();
  vertexVect.clear();
  hitsVect.clear();

  fNNuInts = 0;
  fIntType = -1;
  fIntMode = -1;

  fGap = -1.;

  fTrueVertex.reco = false;
  fTrueVertex.x = -1e3; fTrueVertex.y = -1e3; fTrueVertex.z = -1e3; fTrueVertex.t = -1e3;
  fTrueVertex.U = -1e3; fTrueVertex.V = -1e3; fTrueVertex.C = -1e3; fTrueVertex.TT = -1e3;

  fDelta_U=defVal2;  fEta_U=defVal;  fFitScore_U=defVal;  fAlpha_U=defVal3;  fOmega_U=defVal;  fTau_U=defVal4;  fIota_U=defVal;
  fDelta_V=defVal2;  fEta_V=defVal;  fFitScore_V=defVal;  fAlpha_V=defVal3;  fOmega_V=defVal;  fTau_V=defVal4;  fIota_V=defVal;
  fDelta_C=defVal2;  fEta_C=defVal;  fFitScore_C=defVal;  fAlpha_C=defVal3;  fOmega_C=defVal;  fTau_C=defVal4;  fIota_C=defVal;
}


void sbnd::FRAMSAna::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // TTree for BDT training
  fTree=tfs->make<TTree>(fTreeName.c_str(), "FRAMS Output Tree");
  fTree->Branch("RunID", &fRunID, "RunID/I");
  fTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");

  fTree->Branch("NNuInts", &fNNuInts, "NNuInts/I");
  fTree->Branch("IntMode", &fIntMode, "IntMode/I");
  fTree->Branch("IntType", &fIntType, "IntType/I");
  fTree->Branch("Gap", &fGap, "Gap/D");
  fTree->Branch("LambdaKE", &fLambdaKE, "LambdaKE/D");
  fTree->Branch("ProtonKE", &fProtonKE, "ProtonKE/D");
  fTree->Branch("PionKE", &fPionKE, "PionKE/D");

  for(size_t k=0; k<fViews.size(); k++){
    if(fViews[k]==0){
      fTree->Branch("Delta_U", &fDelta_U, "Delta_U/D");
      fTree->Branch("Eta_U", &fEta_U, "Eta_U/D");
      fTree->Branch("FitScore_U", &fFitScore_U, "FitScore_U/D");
      fTree->Branch("Alpha_U", &fAlpha_U, "Alpha_U/D");
      fTree->Branch("Omega_U", &fOmega_U, "Omega_U/D");
      fTree->Branch("Tau_U", &fTau_U, "Tau_U/D");
      fTree->Branch("Iota_U", &fIota_U, "Iota_U/D");
    }
    if(fViews[k]==1){
      fTree->Branch("Delta_V", &fDelta_V, "Delta_V/D");
      fTree->Branch("Eta_V", &fEta_V, "Eta_V/D");
      fTree->Branch("FitScore_V", &fFitScore_V, "FitScore_V/D");
      fTree->Branch("Alpha_V", &fAlpha_V, "Alpha_V/D");
      fTree->Branch("Omega_V", &fOmega_V, "Omega_V/D");
      fTree->Branch("Tau_V", &fTau_V, "Tau_V/D");
      fTree->Branch("Iota_V", &fIota_V, "Iota_V/D");
    }
    if(fViews[k]==2){
      fTree->Branch("Delta_C", &fDelta_C, "Delta_C/D");
      fTree->Branch("Eta_C", &fEta_C, "Eta_C/D");
      fTree->Branch("FitScore_C", &fFitScore_C, "FitScore_C/D");
      fTree->Branch("Alpha_C", &fAlpha_C, "Alpha_C/D");
      fTree->Branch("Omega_C", &fOmega_C, "Omega_C/D");
      fTree->Branch("Tau_C", &fTau_C, "Tau_C/D");
      fTree->Branch("Iota_C", &fIota_C, "Iota_C/D");
    }
  }
  fNAnalyzedEvents=0;
}

void sbnd::FRAMSAna::endJob(){}


float sbnd::FRAMSAna::GetSliceNuScore(art::Ptr<recob::PFParticle> &nupfp,
  std::vector<art::Ptr<recob::PFParticle>> pfpVect, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns){

  float nu_score = -1.;

  for(size_t ix = 0; ix<pfpVect.size(); ix++ ){
    const art::Ptr<recob::PFParticle> &pfp = pfpVect[ix];
    //std::cout<<"PFParticlePDG:"<<pfp->PdgCode()<<" Primary="<<pfp->IsPrimary()<< " NDaughters="<<pfp->NumDaughters()<<" ID="<<pfp->Self()<<std::endl;
    if(pfp->IsPrimary() && ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 )  ){
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfp_pfpmeta_assns.at(pfp.key());
      art::Ptr<larpandoraobj::PFParticleMetadata> nupfpMeta = pfpMetaVec[0];
      std::map<std::string, float> nupfp_PropMap = nupfpMeta->GetPropertiesMap();

      if ( nupfp_PropMap.find("NuScore") != nupfp_PropMap.end() ){
        //std::cout<<"Score: "<<nupfp_PropMap["NuScore"]<<std::endl;
        nu_score = nupfp_PropMap["NuScore"];
      }
      //for(auto &prop: nupfp_PropMap){std::cout<<prop.first<<"  "<<prop.second<<std::endl;}
      nupfp = pfp;
    }
  }

  return nu_score;
}

double sbnd::FRAMSAna::GetDistance(double p1[], double p2[]){
  return std::sqrt( (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2])  );
}


DEFINE_ART_MODULE(sbnd::FRAMSAna)
