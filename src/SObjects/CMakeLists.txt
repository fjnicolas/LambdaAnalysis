message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(sobjects_library
    TPCSimpleHits.cpp
    TPCSimpleHits.h
    TPCLinesHitDistanceUtils.cpp
    TPCLinesHitDistanceUtils.h
    TPCSimpleClusters.cpp
    TPCSimpleClusters.h
    TPCSimpleLines.cpp
    TPCSimpleLines.h
    TPCSimpleTriangles.cpp
    TPCSimpleTriangles.h
    TPCLinesPCA.cpp
    TPCLinesPCA.h
    TPCLinesParameters.cpp
    TPCLinesDistanceUtils.cpp
    TPCLinesDistanceUtils.h
)

add_library(sevents_library
    TPCSimpleEvents.cpp
    TPCSimpleEvents.h
)



# Include directories for header files
target_include_directories(sobjects_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(sobjects_library PUBLIC ROOT::Hist ROOT::Gpad)


# Include directories for header files
target_include_directories(sevents_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(sevents_library PUBLIC ROOT::Hist ROOT::Gpad)
target_link_libraries(sevents_library PRIVATE sobjects_library)