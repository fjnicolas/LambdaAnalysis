message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(seventhandle_library
    CommandLineParser.cpp
    CommandLineParser.h
    SEventHandle.cpp
    SEventHandle.h
    STPCAnalyzerTreeReader.cpp
    STPCAnalyzerTreeReader.h
    SParserConfig.h
)

# Include directories for header files
target_include_directories(seventhandle_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(seventhandle_library PUBLIC frans_library tpclines_library sobjects_library ROOT::Tree  ROOT::Tree Boost::filesystem)