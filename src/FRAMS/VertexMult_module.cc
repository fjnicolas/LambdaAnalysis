////////////////////////////////////////////////////////////////////////
// Class:       VertexMult
// Plugin Type: analyzer (art v3_05_01)
// File:        VertexMult_module.cc
//
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
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

#include "sbndcode/FRAMS/ChargeDensity/ChargeDensity.hh"
#include "sbndcode/FRAMS/ChargeDensity/ChargeDensityAlgConf.hh"

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
#define fMaxRadius 17

namespace test {
  class VertexMult
;
}


class test::VertexMult : public art::EDAnalyzer {
public:

  struct Config {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;

      fhicl::Atom<std::string> MCTruthLabel {Name("MCTruthLabel")};
      fhicl::Atom<std::string> HitLabel {Name("HitLabel")};
      fhicl::Atom<std::string> RecoLabel {Name("RecoLabel")};
      fhicl::Atom<bool> ApplyFiducialCut {Name("ApplyFiducialCut")};
      fhicl::Atom<std::string> TreeName {Name("TreeName")};

    }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit VertexMult(Parameters const& config);

  // Plugins should not be copied or assigned.
  VertexMult(VertexMult const&) = delete;
  VertexMult(VertexMult&&) = delete;
  VertexMult & operator=(VertexMult const&) = delete;
  VertexMult & operator=(VertexMult &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  void resetVars();
  int VertexToDriftTick(double vt, double vx);

  // configuration parameters
  std::string fMCTruthLabel;
  std::string fHitLabel;
  std::string fRecoLabel;
  bool fApplyFiducialCut;
  std::string fTreeName;


  TTree* fTree;
  int fEventID, fRunID, fSubRunID;

  //True variables
  std::vector<int> fTruePrimariesPDG;
  std::vector<double> fTruePrimariesE;
  double fTrueVx;
  double fTrueVy;
  double fTrueVz;
  double fTrueVt;
  int fTrueVU;
  int fTrueVV;
  int fTrueVC;
  int fTrueVTimeTick;
  double fTrueVEnergy;


  //Reconstructed vertex
  double fRecoVx;
  double fRecoVy;
  double fRecoVz;
  int fRecoVU;
  int fRecoVV;
  int fRecoVC;
  int fRecoVTimeTick;

  int fNAnalyzedEvents;
  int fNVertexInRadius;

  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  unsigned int fNChannels;
  unsigned int fReadoutWindow;
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;



  //Load TFileService serrvice
  art::ServiceHandle<art::TFileService> tfs;
};


test::VertexMult::VertexMult(Parameters const& config)
  : EDAnalyzer{config},
  fMCTruthLabel( config().MCTruthLabel() ),
  fHitLabel( config().HitLabel() ),
  fRecoLabel( config().RecoLabel() ),
  fApplyFiducialCut( config().ApplyFiducialCut() ),
  fTreeName( config().TreeName() ),
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

}


void test::VertexMult::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  std::cout<<"Running VertexMult---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<"\n";
  //auto const fClockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  //auto const fDetProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, fClockData);

  //............................Event General Info
  fNAnalyzedEvents++;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fEventID = e.id().event();
  //Reset tree variables
  resetVars();

  //............................Read Truth Objects
  art::Handle<std::vector<simb::MCTruth>> mctruths;
  e.getByLabel(fMCTruthLabel, mctruths);
  std::cout<<" --- Saving MCTruth\n";

  for (auto const& truth : *mctruths) {

    for (int p = 0; p < truth.NParticles(); p++){
      simb::MCParticle const& TruePart = truth.GetParticle(p);

      if( TruePart.StatusCode()==1 ){
        fTruePrimariesPDG.push_back( TruePart.PdgCode() );
        fTruePrimariesE.push_back( TruePart.E() );
      }

      if( TruePart.Mother()==-1 && ( abs(TruePart.PdgCode())==12 || abs(TruePart.PdgCode())==14 ) ){
        fTrueVx=TruePart.EndX();
        fTrueVy=TruePart.EndY();
        fTrueVz=TruePart.EndZ();
        fTrueVt=TruePart.T();
        fTrueVEnergy=TruePart.E();

        if(fApplyFiducialCut && std::abs(fTrueVx)<fXFidCut && std::abs(fTrueVy)<fYFidCut && fTrueVz>fZFidCut1 && fTrueVz<fZFidCut2){
          geo::Point_t p={fTrueVx, fTrueVy, fTrueVz};
          std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(p))<<" "<<fGeom->FindTPCAtPosition(p).TPC<<std::endl;


          if( fGeom->HasTPC(fGeom->FindTPCAtPosition(p)) ){
            unsigned int tpcID=fGeom->FindTPCAtPosition(p).TPC;
            fTrueVU=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 0));
            fTrueVV=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 1));
            fTrueVC=fGeom->NearestChannel(p, geo::PlaneID(0, tpcID, 2));
            fTrueVTimeTick=VertexToDriftTick(fTrueVt, fTrueVx);
          }
        }
        else{
          fTrueVU=-1;
          fTrueVV=-1;
          fTrueVC=-1;
          fTrueVTimeTick=-1;
        }

        std::cout<<"  -- Vertex: "<<fTrueVx<<" "<<fTrueVy<<" "<<fTrueVz<<std::endl;
        std::cout<<"   - VertexWire: "<<fTrueVU<<" "<<fTrueVV<<" "<<fTrueVC<<" "<<fTrueVTimeTick<<std::endl;
      }
    }
  }


  // store the vertex view
  VertexView vertex;
  vertex.reco = false;
  vertex.x = fTrueVx; vertex.y = fTrueVy; vertex.z = fTrueVz; vertex.t = fTrueVt;
  vertex.U = fTrueVU; vertex.V = fTrueVV; vertex.C = fTrueVC; vertex.TT = fTrueVTimeTick;

  int tpc=0;
  if(vertex.x>0) tpc=1;
  std::cout<<tpc<<std::endl;

  bool inFidVolume=true;
  if(fApplyFiducialCut && (std::abs(vertex.x)>fXFidCut ||
    std::abs(vertex.y)>fYFidCut || vertex.z<fZFidCut1 || vertex.x>fZFidCut2) )
  {
    inFidVolume=false;
  }

  if(inFidVolume==true){

    //............................Read Recob Slice
    ::art::Handle<std::vector<recob::Slice>> sliceHandle;
    e.getByLabel(fRecoLabel, sliceHandle);
    //............................Read PFPs
    ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    e.getByLabel(fRecoLabel, pfpHandle);
    //............................Read vertex
    ::art::Handle<std::vector<recob::Vertex>> vertexHandle;
    e.getByLabel(fRecoLabel, vertexHandle);
    //Vector for recob PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfpVect;
    //Vector for recob Vertex
    std::vector<art::Ptr<recob::Vertex>> vertexVect;
    //Slice association for PFParticles
    art::FindManyP<recob::PFParticle> slice_pfp_assns (sliceHandle, e, fRecoLabel);
    //Get Vertex Association
    art::FindManyP<recob::Vertex> vertexAssoc (pfpHandle, e, fRecoLabel);
    //PDParticle associated metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns(pfpHandle, e, fRecoLabel);

    // Loop over slices
    fNVertexInRadius=0;

    std::vector< art::Ptr<recob::Slice> > sliceVect;
    art::fill_ptr_vector(sliceVect, sliceHandle);
    std::cout<<" Number of slices: "<<sliceVect.size()<<std::endl;
    for(auto & slice:sliceVect){
      pfpVect = slice_pfp_assns.at(slice.key());
      std::cout<<"   *** PFParticle size:"<<pfpVect.size()<<std::endl;

      size_t neutrinoID = fDefaulNeutrinoID;
      size_t neutrinoID_ix = fDefaulNeutrinoID;
      double xyz_nu_vertex[3];
      for(size_t ix = 0; ix<pfpVect.size(); ix++ ){
        const art::Ptr<recob::PFParticle> &pfp = pfpVect[ix];
        std::cout<<"PFParticlePDG:"<<pfp->PdgCode()<<" Primary="<<pfp->IsPrimary()<< " NDaughters="<<pfp->NumDaughters()
        <<" ID="<<pfp->Self()<<std::endl;
        if( !( pfp->IsPrimary() && ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ) ) continue;
        neutrinoID = pfp->Self();
        neutrinoID_ix = ix;
        vertexVect = vertexAssoc.at(pfp.key());
        for(const art::Ptr<recob::Vertex> &ver : vertexVect){
          ver->XYZ(xyz_nu_vertex);
        }
      }
      std::vector<double> nu_vec(std::begin(xyz_nu_vertex), std::end(xyz_nu_vertex));


      std::cout<<"    ** NeutrinoID:"<<neutrinoID<<" IX="<<neutrinoID_ix<<"\n\n";

      if(neutrinoID != fDefaulNeutrinoID){
        const art::Ptr<recob::PFParticle> &nupfp = pfpVect[neutrinoID_ix];

        const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfp_pfpmeta_assns.at(nupfp.key());
        art::Ptr<larpandoraobj::PFParticleMetadata> nupfpMeta = pfpMetaVec[0];
        std::map<std::string, float> nupfp_PropMap = nupfpMeta->GetPropertiesMap();

        for(auto &prop: nupfp_PropMap){
          std::cout<<prop.first<<"  "<<prop.second<<std::endl;
        }

        // lets first print the neutrino daughters
        std::cout<<"  Neutrino PFParticle   Daughters="<<nupfp->NumDaughters()<<std::endl;
        for(int d = 0; d<nupfp->NumDaughters(); d++ ){
          for(const art::Ptr<recob::PFParticle> &dpfp : pfpVect){
            if( nupfp->Daughter(d) != dpfp->Self() ) continue;
            std::cout<<"     Daughter "<<d<<"  PDG="<<dpfp->PdgCode()<<" ID="<<dpfp->Self()<<" Daughters="<<dpfp->NumDaughters()<<std::endl;
          }
        }

        std::cout<<"\n\n";
        //PFParticle's vertex loop
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
          if(pfp->Self()==neutrinoID) continue;
          std::cout<<"     PFParticle: "<<pfp->Self()<<" PDG: "<<pfp->PdgCode()<<std::endl;
          vertexVect = vertexAssoc.at(pfp.key());
          double xyz_vertex[3];
          for(const art::Ptr<recob::Vertex> &ver : vertexVect){
            ver->XYZ(xyz_vertex);
            //fpfpVertexPosition.push_back( xyz_vec );
            double chi2=ver->chi2(), chi2ndof=ver->chi2PerNdof();
            std::cout<<"  --VERTEX  ID="<<ver->ID()<<"  x,y,z="<<xyz_vertex[0]<<","<<xyz_vertex[1]<<","<<xyz_vertex[2];
            std::cout<<" Chi2="<<chi2<<" Chi2/DoF="<<chi2ndof<<" Status:"<<ver->status()<<",\n";
          }
          std::vector<double> xyz_vec(std::begin(xyz_vertex), std::end(xyz_vertex));
          double d = std::sqrt( (xyz_vec[0]-nu_vec[0])*(xyz_vec[0]-nu_vec[0])+(xyz_vec[1]-nu_vec[1])*(xyz_vec[1]-nu_vec[1])+(xyz_vec[2]-nu_vec[2])*(xyz_vec[2]-nu_vec[2])  );
          if(d<fMaxRadius) fNVertexInRadius++;
        }
      }
    }

    std::cout<<"   HHHHHHH fNVertexInRadius="<<fNVertexInRadius<<std::endl;
    fTree->Fill();
  }


}


int test::VertexMult::VertexToDriftTick(double vt, double vx){
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}


void test::VertexMult::resetVars()
{
}


void test::VertexMult::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fTree=tfs->make<TTree>(fTreeName.c_str(), "FRAMS Output Tree");
  fTree->Branch("RunID", &fRunID, "RunID/I");
  fTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");
  fTree->Branch("NVertexInRadius", &fNVertexInRadius, "NVertexInRadius/I");

  fNAnalyzedEvents=0;
}

void test::VertexMult::endJob(){}

DEFINE_ART_MODULE(test::VertexMult)
