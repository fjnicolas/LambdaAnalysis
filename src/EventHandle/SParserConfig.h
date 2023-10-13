////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_PARSERCONFIG_PARAMETERS_H
#define TPC_PARSERCONFIG_PARAMETERS_H

#include "ChargeDensityPset.h"
#include "TPCLinesParameters.h"

#include <boost/property_tree/info_parser.hpp>
#include <filesystem>
#include <fstream>
#include <map>

namespace fs = std::filesystem;

std::string FindFile(const std::string& fileName, const char* envVarName="LAMBDAANA_SEARCHPATH") {
    
    const char* envVarValue = std::getenv(envVarName);

    std::string fullPath="";

    if (envVarValue!=nullptr) {
        
        std::istringstream paths(envVarValue);
        std::string path;
        while (std::getline(paths, path, ':')) { 
            fs::path directoryPath(path);

            if (fs::is_directory(directoryPath)) {
                for (const auto& entry : fs::directory_iterator(directoryPath)) {
                    
                    if (fs::is_regular_file(entry) && entry.path().filename() == fileName) {
                        fullPath = entry.path().string(); // Return the full path of the file.
                    }
                }
            }
        }
    }
    else{
        std::cerr << "PATH environmental variable not found." << std::endl;
        return "";
    }

    return fullPath;
}


boost::property_tree::ptree GetPropertyTreeFromFileName(std::string filename, std::string blockName){
    
    // Define the section names
    const std::string begin_prolog = "BEGIN_PROLOG";
    const std::string end_prolog = "END_PROLOG";

    // Create a property tree
    boost::property_tree::ptree pt;

    // Open the configuration file
    std::ifstream config_file(filename);
    if (config_file.is_open()) {

        std::string line;
        bool inside_prolog = false;
        bool inside_block = false;

        while (std::getline(config_file, line)) {
            if (line == begin_prolog){
                inside_prolog = true;
            }
            else if (line == end_prolog) {
                inside_prolog = false;
            }
            else if (line ==  blockName) {
                inside_block = true;
            }
            else if (inside_prolog && inside_block) {
                
                if (!line.empty() && line[0] != '#') {  // Ignore empty lines and comments
                    size_t separator_pos = line.find(":");
                    if (separator_pos != std::string::npos) {
                        std::string key = line.substr(0, separator_pos);
                        std::string value = line.substr(separator_pos + 1);
                        
                        // Remove spaces
                        key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
                        value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());

                        // Add the key-value pair to the property tree
                        pt.put(key, value);
                    }
                }
                
            }
        }

        config_file.close();

    }



    return pt;
}


FRAMSPsetType ReadFRANSPset(std::string filename, std::string blockName){

    // Create the parameter set
    FRAMSPsetType fPset;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPset.ApplyRawSmoothing = pt.get<bool>("ApplyRawSmoothing");
    fPset.ApplySmoothing = pt.get<bool>("ApplySmoothing");
    fPset.ApplyCumulativeSmoothing = pt.get<bool>("ApplyCumulativeSmoothing");
    fPset.NDriftPack = pt.get<unsigned int>("NDriftPack");
    fPset.NWirePack = pt.get<unsigned int>("NWirePack");
    fPset.ExpoAvSmoothPar = pt.get<double>("ExpoAvSmoothPar");
    fPset.UnAvNeighbours = pt.get<int>("UnAvNeighbours");
    fPset.CumulativeCut = pt.get<double>("CumulativeCut");
    fPset.SlidingWindowN = pt.get<int>("SlidingWindowN");
    fPset.NSamplesBeginSlope = pt.get<int>("NSamplesBeginSlope");
    fPset.MaxRadius = pt.get<int>("MaxRadius");
    fPset.Verbose = pt.get<int>("Verbose");
    fPset.CalculateScore = pt.get<bool>("CalculateScore");
    fPset.TMVAFilename = pt.get<std::string>("TMVAFilename");
    fPset.TMVAFilename = fPset.TMVAFilename.substr(1, fPset.TMVAFilename.size()-2);
    fPset.OutputPath = pt.get<std::string>("OutputPath");
    
    return fPset;
}


TrackFinderAlgorithmPsetType ReadTrackFinderAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    TrackFinderAlgorithmPsetType fPset;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPset.MaxDTube = pt.get<int>("MaxDTube");
    fPset.MaxDCluster = pt.get<double>("MaxDCluster");
    fPset.SingleWireMode = pt.get<bool>("SingleWireMode");
    fPset.MinClusterHits = pt.get<int>("MinClusterHits");
    fPset.DCleaning = pt.get<double>("DCleaning");
    fPset.ClusterCompletenessCut = pt.get<double>("ClusterCompletenessCut");
    fPset.ClusterAngleCut = pt.get<double>("ClusterAngleCut");
    fPset.CaptureMissingHits = pt.get<bool>("CaptureMissingHits");
    fPset.MinTrackHits = pt.get<int>("MinTrackHits");
    fPset.HitDensityThreshold = pt.get<float>("HitDensityThreshold");
    fPset.Verbose = pt.get<int>("Verbose");
    
    
    return fPset;
}


VertexFinderAlgorithmPsetType ReadVertexFinderAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    VertexFinderAlgorithmPsetType fPset;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPset.MaxDistToEdge = pt.get<double>("MaxDistToEdge");
    fPset.RefineVertexIntersection = pt.get<bool>("RefineVertexIntersection");
    fPset.UseEdgesDiscard = pt.get<bool>("UseEdgesDiscard");
    fPset.MaxTrackFractionInMain = pt.get<float>("MaxTrackFractionInMain");
    fPset.DecideMainTrack = pt.get<bool>("DecideMainTrack");
    fPset.AddCollinearLines = pt.get<bool>("AddCollinearLines");
    fPset.Verbose = pt.get<int>("Verbose");
    
    
    return fPset;
}

HoughAlgorithmPsetType ReadHoughAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    HoughAlgorithmPsetType fPset;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPset.MaxRadiusLineHypothesis = pt.get<double>("MaxRadiusLineHypothesis");
    fPset.ThetaRes = pt.get<double>("ThetaRes");
    fPset.MaxDistanceTube = pt.get<double>("MaxDistanceTube");
    fPset.MinHoughHits = pt.get<int>("MinHoughHits");
    fPset.Verbose = pt.get<int>("Verbose");
    
    
    return fPset;
}

TPCLinesAlgoPsetType ReadTPCLinesAlgoPset(std::string filename, std::string blockName){

    // Create the parameter set
    TPCLinesAlgoPsetType fPset;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPset.MaxRadius = pt.get<double>("MaxRadius");
    fPset.DriftConversion = pt.get<double>("DriftConversion");
    fPset.MaxHoughTracks = pt.get<int>("MaxHoughTracks");
    fPset.MinTrackHits = pt.get<int>("MinTrackHits");
    fPset.RemoveIsolatedHits = pt.get<bool>("RemoveIsolatedHits");
    fPset.MaxNeighbourDistance = pt.get<double>("MaxNeighbourDistance");
    fPset.MinNeighboursHits = pt.get<int>("MinNeighboursHits");
    fPset.VertexAlgorithm = pt.get<int>("VertexAlgorithm");
    fPset.View = pt.get<std::string>("View");
    fPset.OutputPath = pt.get<std::string>("OutputPath");
    fPset.Verbose = pt.get<int>("Verbose");
    fPset.DebugMode = pt.get<int>("DebugMode");
    
    
    return fPset;
}

#endif // TPC_LINES_PARAMETERS_H