art_make_library(
    LIBRARY_NAME   sbndcode_FRAMS_ChargeDensity
    SOURCE
	    ChargeDensity.cc
      ChargeDensityAlgConf.hh
      LeastSquares.hh
      VertexView.hh
    LIBRARIES
      larcore::Geometry_Geometry_service
      lardataobj::RecoBase
      larsim::Simulation
      lardataobj::Simulation
      nusimdata::SimulationBase
      canvas::canvas
      cetlib::cetlib
      cetlib_except::cetlib_except
      art::Framework_Core
      art::Framework_Principal
      art::Framework_Services_Registry
      art_root_io::tfile_support ROOT::Core
      art_root_io::TFileService_service
      art::Persistency_Common
      art::Persistency_Provenance
      art::Utilities
      messagefacility::MF_MessageLogger
      fhiclcpp::fhiclcpp
      ${ROOT_BASIC_LIB_LIST}
      ROOT::Geom
      ROOT::XMLIO
      ROOT::Gdml
      ROOT::Core
      nug4::ParticleNavigation
      ROOT::TMVA
)


install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)