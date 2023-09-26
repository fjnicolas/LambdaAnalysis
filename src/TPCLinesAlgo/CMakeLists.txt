message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

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
)

# Include directories for header files
#target_include_directories(sobjects_library PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link math_library to string_library
target_link_libraries(tpclines_library PRIVATE ROOT::Hist ROOT::Graf ROOT::Gpad)
target_link_libraries(tpclines_library PRIVATE sobjects_library sevents_library)

