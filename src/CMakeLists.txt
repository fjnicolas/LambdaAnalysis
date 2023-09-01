# Find the ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

set(Objects_DIR "${CMAKE_BINARY_DIR}/SObjects")
message(STATUS "CMake: SObjects ${Objects_DIR}")
find_package(TPCSimpleHits CONFIG )
find_package(TPCSimpleLines CONFIG )
find_package(TPCSimpleClusters CONFIG )
find_package(TPCSimpleTriangles CONFIG )
message(STATUS "---> ${TPCSimpleHits_FOUND}")


# Compile the C++ code into a shared library
add_library(TPCLinesDistanceUtils TPCLinesDistanceUtils.cpp)
target_link_libraries (TPCLinesDistanceUtils LINK_SHARED TPCSimpleHits)

