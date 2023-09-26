message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(seventhandle_library
    SeventHandle.cpp
    SEventHandle.h
    STPCAnalyzerTreeReader.cpp
    STPCAnalyzerTreeReader.h
)

# Include directories for header files
target_include_directories(seventhandle_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(seventhandle_library PUBLIC sobjects_library ROOT::Tree)
