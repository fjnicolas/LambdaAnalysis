////////////////////////////////////////////////////////////////////////
// Class:       FRAMSFilter
// Plugin Type: filter (Unknown Unknown)
// File:        FRAMSFilter_module.cc
//
// Generated at Thu Oct 21 03:39:48 2021 by Francisco Nicolas-Arnaldos using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Geometry
#include "larcore/Geometry/Geometry.h"

//#include "sbndcode/FRAMS/FRAMSObj/FRAMSObj.h"
#include "sbnobj/SBND/FRAMSObj/FRAMSObj.h"

namespace sbnd {
  class FRAMSFilter;
}


class sbnd::FRAMSFilter : public art::EDFilter {
public:
  explicit FRAMSFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FRAMSFilter(FRAMSFilter const&) = delete;
  FRAMSFilter(FRAMSFilter&&) = delete;
  FRAMSFilter& operator=(FRAMSFilter const&) = delete;
  FRAMSFilter& operator=(FRAMSFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  std::string fFRAMSLabel;
  double fScoreCut;
  bool fSelectionMode;
  bool fUseGapCut;
  double fGapCut;

  // Declare member data here.
  bool fKeepEvent;
};


sbnd::FRAMSFilter::FRAMSFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
  fFRAMSLabel ( p.get<std::string>("FRAMSLabel") ),
  fScoreCut ( p.get<double>("ScoreCut") ),
  fSelectionMode ( p.get<bool>("SelectionMode") ),
  fUseGapCut ( p.get<bool>("UseGapCut") ),
  fGapCut ( p.get<double>("GapCut") )
{
}

bool sbnd::FRAMSFilter::filter(art::Event& e)
{
  fKeepEvent=true;
  //Read MCTruth
  art::Handle< std::vector<sbnd::FRAMSObj> > FRAMSObjHandle;
  e.getByLabel(fFRAMSLabel, FRAMSObjHandle);

  if(FRAMSObjHandle->size()==0) fKeepEvent = false;
  else{
    for (size_t n = 0; n < FRAMSObjHandle->size(); n++) {
      sbnd::FRAMSObj const&  framsObj = (*FRAMSObjHandle)[n];

      std::cout<<n<<" DeltaC: "<<framsObj.Delta_C<<" Score: "<<framsObj.ScoreBDT_C<<std::endl;

      // if selection mode, keep everything above BDT cut
      if(fSelectionMode) fKeepEvent =  (framsObj.ScoreBDT_C>fScoreCut);
      // else, we keep the misselections, deopending on the event being signal or BG
      else{
        // if is BG, we keep the event if the score is greater than the cut, aka misselected as Lambda
        if(framsObj.IsSignal==0){
          fKeepEvent =  (framsObj.ScoreBDT_C>fScoreCut);
        }
        // if is BG, we keep the event if the score is smaller than the cut, aka not selected
        else{
          fKeepEvent =  (framsObj.ScoreBDT_C<fScoreCut);
          if(fUseGapCut) fKeepEvent = (fKeepEvent && framsObj.Gap>fGapCut);
        }
      }

    }
  }

  std::cout<<"KeepEvent "<<fKeepEvent<<std::endl;

  return fKeepEvent;
}

void sbnd::FRAMSFilter::beginJob()
{
}

void sbnd::FRAMSFilter::endJob()
{
}

DEFINE_ART_MODULE(sbnd::FRAMSFilter)
