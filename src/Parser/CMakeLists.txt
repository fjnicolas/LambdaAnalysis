message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(commandlineparser_library
    CommandLineParser.cpp
    CommandLineParser.h
)


# Include directories for header files
target_include_directories(commandlineparser_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})