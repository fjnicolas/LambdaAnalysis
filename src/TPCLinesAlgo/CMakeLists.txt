message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED SOBJECTS" )

    art_make_library(
        LIBRARY_NAME   sbndcode_TPCLines
        
        SOURCE
            TPCLinesAlgo.cpp            
            TPCLinesDisplay.cpp            
            TPCLinesHough.cpp            
            TPCLinesTrackFinder.cpp            
            TPCLinesVertexFinder.cpp            
            TPCLinesDBSCAN.cpp            
            TPCLinesDirectionRecoUtils.cpp            
            TPCLinesParameters.cpp
        LIBRARIES
            ROOT::Tree
            ROOT::Core

            sbndcode_SObjects
            sbndcode_SEventHandle
      )

      install_headers()
      install_source()
endif()


add_library(tpclines_library
    TPCLinesAlgo.cpp
    TPCLinesAlgo.h
    TPCLinesDisplay.cpp
    TPCLinesDisplay.h
    TPCLinesHough.cpp
    TPCLinesHough.h
    TPCLinesTrackFinder.cpp
    TPCLinesTrackFinder.h
    TPCLinesVertexFinder.cpp
    TPCLinesVertexFinder.h
    TPCLinesDBSCAN.cpp
    TPCLinesDBSCAN.h
    TPCLinesDirectionRecoUtils.cpp
    TPCLinesDirectionRecoUtils.h
    TPCLinesParameters.cpp
    TPCLinesParameters.h
)

# Include directories for header files
#target_include_directories(sobjects_library PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link math_library to string_library
target_link_libraries(tpclines_library PRIVATE ROOT::Hist ROOT::Graf ROOT::Gpad)
target_link_libraries(tpclines_library PRIVATE sobjects_library)

