message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(frans_library
  LeastSquares.h
  ChargeDensity.cpp
  ChargeDensity.h
  ChargeDensityPset.h
  VertexView.h
)


# Include directories for header files
target_include_directories(frans_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(frans_library PUBLIC sobjects_library ROOT::Hist ROOT::Gpad ROOT::TMVA)