# find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Net)

# Add subdirectories for libraries
add_subdirectory(Parser)
add_subdirectory(SObjects)
add_subdirectory(TPCLinesAlgo)


# Add the main executable
add_executable(RunAlgoTPCLines RunAlgoTPCLines.cpp)

# Link the libraries to the main executable
target_include_directories(RunAlgoTPCLines PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/TPCLinesAlgo)
target_link_libraries(RunAlgoTPCLines PUBLIC sobjects_library tpclines_library sevents_library commandlineparser_library ROOT::Tree)