# find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Net)
find_package(Boost REQUIRED COMPONENTS system filesystem)

# Add subdirectories for libraries
add_subdirectory(EventHandle)
add_subdirectory(SObjects)
add_subdirectory(TPCLinesAlgo)
add_subdirectory(FRANS)


# Add main executable
add_executable(RunAlgoTPCLines RunAlgoTPCLines.cpp)
# Link the libraries to the main executable
target_include_directories(RunAlgoTPCLines PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/TPCLinesAlgo ${CMAKE_CURRENT_SOURCE_DIR}/FRANS)
target_link_libraries(RunAlgoTPCLines PUBLIC  sobjects_library frans_library  seventhandle_library tpclines_library ROOT::Tree)

# Add main executable
add_executable(RunAlgoFRANS RunAlgoFRANS.cpp)
# Link the libraries to the main executable
target_include_directories(RunAlgoFRANS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FRANS ${CMAKE_CURRENT_SOURCE_DIR}/TPCLinesAlgo)
target_link_libraries(RunAlgoFRANS PUBLIC seventhandle_library sobjects_library tpclines_library frans_library)