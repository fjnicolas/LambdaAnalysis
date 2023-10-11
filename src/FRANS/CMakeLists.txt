message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(frans_library
  LeastSquares.h
  ChargeDensity.cpp
  ChargeDensity.h
  ChargeDensityPset.h
  FRANSTTreeHandle.h
  VertexView.h
)


# Include directories for header files
target_include_directories(frans_library PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(frans_library PRIVATE sobjects_library)
target_link_libraries(frans_library PRIVATE ROOT::Hist ROOT::Gpad)
target_link_libraries(frans_library PUBLIC ROOT::TMVA)
