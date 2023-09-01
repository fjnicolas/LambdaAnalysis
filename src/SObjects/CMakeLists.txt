message(STATUS "CMake: SObjects")

# Find the ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Compile the C++ code into a shared library
add_library(TPCSimpleHits TPCSimpleHits.cpp TPCSimpleHits.h)

target_include_directories(TPCSimpleHits PUBLIC "${CURRENT_SOURCE_DIR}")
install(TARGETS TPCSimpleHits EXPORT TTPCSimpleHits LIBRARY DESTINATION lib/)
#export(TARGETS TPCSimpleHits NAMESPACE FILE TTPCSimpleHits ${CMAKE_BINARY_DIR}/SObjects/TPCSimpleHits/Config.cmake)

add_library(SObjects::TPCSimpleHits ALIAS TPCSimpleHits)
file(WRITE ${CMAKE_BINARY_DIR}/SObjects/TPCSimpleHits.cmake "")

add_library(TPCSimpleLines TPCSimpleLines.cpp TPCSimpleLines.h)
target_link_libraries (TPCSimpleLines LINK_SHARED TPCSimpleHits)


add_library(TPCSimpleClusters TPCSimpleClusters.cpp TPCSimpleClusters.h)
target_link_libraries (TPCSimpleClusters LINK_SHARED TPCSimpleHits TPCSimpleLines TPCLinesPCA TPCLinesDistanceUtils)


add_library(TPCLinesParameters TPCLinesParameters.cpp)

add_library(TPCSimpleTriangles TPCSimpleTriangles.cpp TPCSimpleTriangles.h)
target_link_libraries (TPCSimpleTriangles LINK_SHARED TPCSimpleHits TPCSimpleLines)



