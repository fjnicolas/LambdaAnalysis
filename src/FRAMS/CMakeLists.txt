set(
MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    nusimdata::SimulationBase
    larsim::MCCheater_ParticleInventoryService_service
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}

    art_root_io::tfile_support
    art_root_io::TFileService_service

    ROOT::Tree
    ROOT::Core

    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::MCBase
    lardataobj::Simulation
    nusimdata::SimulationBase
    nug4::ParticleNavigation

    larsim::Utils
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_PhotonBackTrackerService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service

    sbndcode::FRAMS_ChargeDensity
    sbnobj::SBND_FRAMSObj
    ROOT::TMVA
)

cet_build_plugin(FRAMSAna art::Module SOURCE FRAMSAna_module.cc LIBRARIES ${MODULE_LIBRARIES} )
cet_build_plugin(FRAMSFilter art::Module SOURCE FRAMSFilter_module.cc LIBRARIES ${MODULE_LIBRARIES} )
cet_build_plugin(VertexMult art::Module SOURCE VertexMult_module.cc LIBRARIES ${MODULE_LIBRARIES} )



install_source()
install_fhicl()

add_subdirectory(job)
add_subdirectory(config)
add_subdirectory(ChargeDensity)
