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
#include <fstream>
#include <map>


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
                std::cout<<"INSIDE PROLOG\n";
                
                if (!line.empty() && line[0] != '#') {  // Ignore empty lines and comments
                    size_t separator_pos = line.find(":");
                    if (separator_pos != std::string::npos) {
                        std::string key = line.substr(0, separator_pos);
                        std::string value = line.substr(separator_pos + 1);
                        
                        // Remove spaces
                        key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
                        value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());

                        std::cout<<key<<" "<<value<<std::endl;

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

    // Create and fill the parameter set
    FRAMSPsetType fPsetFRANS;
    
    
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);


    fPsetFRANS.ApplyRawSmoothing = pt.get<bool>("ApplyRawSmoothing");
    fPsetFRANS.ApplySmoothing = pt.get<bool>("ApplySmoothing");
    fPsetFRANS.ApplyCumulativeSmoothing = pt.get<bool>("ApplyCumulativeSmoothing");
    fPsetFRANS.NDriftPack = pt.get<unsigned int>("NDriftPack");
    fPsetFRANS.NWirePack = pt.get<unsigned int>("NWirePack");
    fPsetFRANS.ExpoAvSmoothPar = pt.get<double>("ExpoAvSmoothPar");
    fPsetFRANS.UnAvNeighbours = pt.get<int>("UnAvNeighbours");
    fPsetFRANS.CumulativeCut = pt.get<double>("CumulativeCut");
    fPsetFRANS.SlidingWindowN = pt.get<int>("SlidingWindowN");
    fPsetFRANS.NSamplesBeginSlope = pt.get<int>("NSamplesBeginSlope");
    fPsetFRANS.MaxRadius = pt.get<int>("MaxRadius");
    fPsetFRANS.Verbose = pt.get<int>("Verbose");
    fPsetFRANS.CalculateScore = pt.get<bool>("CalculateScore");
    fPsetFRANS.TMVAFilename = pt.get<std::string>("TMVAFilename");
    fPsetFRANS.OutputPath = pt.get<std::string>("OutputPath");

    std::cout << "Key1: " << pt.get<bool>("ApplyRawSmoothing") << std::endl;

    std::cout << "Key2: " << pt.get<std::string>("TMVAFilename") << std::endl;
    
    return fPsetFRANS;

}



FRAMSPsetType ReadFRANSPset2(std::string filename){
   
    // Create a property tree to hold the parsed data
    boost::property_tree::ptree pt;

    // Create an instance of the struct and populate its fields
    FRAMSPsetType fPsetFRANS;
    std::string fFRAMSPsetProlog = "FRANS_";

    // Parse the FHiCL file
    try {
        // Open and read the configuration file
        std::ifstream config_file(filename);
        const char custom_separator = ':';
        boost::property_tree::read_info(config_file, pt);

        //boost::property_tree::ptree frans_block = pt.get_child(fFRAMSPsetProlog);

        std::cout<<"IIIIIIIN\n";

        fPsetFRANS.ApplyRawSmoothing = pt.get<bool>(fFRAMSPsetProlog+"ApplyRawSmoothing");
        std::cout<<"Apply "<< fPsetFRANS.ApplyRawSmoothing<<std::endl;;
        fPsetFRANS.ApplySmoothing = pt.get<bool>(fFRAMSPsetProlog+"ApplySmoothing");
        fPsetFRANS.ApplyCumulativeSmoothing = pt.get<bool>(fFRAMSPsetProlog+"ApplyCumulativeSmoothing");
        fPsetFRANS.NDriftPack = pt.get<unsigned int>(fFRAMSPsetProlog+"NDriftPack");
        fPsetFRANS.NWirePack = pt.get<unsigned int>(fFRAMSPsetProlog+"NWirePack");
        fPsetFRANS.ExpoAvSmoothPar = pt.get<double>(fFRAMSPsetProlog+"ExpoAvSmoothPar");
        fPsetFRANS.UnAvNeighbours = pt.get<int>(fFRAMSPsetProlog+"UnAvNeighbours");
        fPsetFRANS.CumulativeCut = pt.get<double>(fFRAMSPsetProlog+"CumulativeCut");
        fPsetFRANS.SlidingWindowN = pt.get<int>(fFRAMSPsetProlog+"SlidingWindowN");
        fPsetFRANS.NSamplesBeginSlope = pt.get<int>(fFRAMSPsetProlog+"NSamplesBeginSlope");
        fPsetFRANS.MaxRadius = pt.get<int>(fFRAMSPsetProlog+"MaxRadius");
        fPsetFRANS.Verbose = pt.get<int>(fFRAMSPsetProlog+"Verbose");
        fPsetFRANS.CalculateScore = pt.get<bool>(fFRAMSPsetProlog+"CalculateScore");
        fPsetFRANS.TMVAFilename = pt.get<std::string>(fFRAMSPsetProlog+"TMVAFilename");
        fPsetFRANS.OutputPath = pt.get<std::string>(fFRAMSPsetProlog+"OutputPath");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return fPsetFRANS;
}


TrackFinderAlgorithmPsetType ReadTrackFinderAlgorithmPset(std::string filename) {
    // Create a property tree to hold the parsed data
    boost::property_tree::ptree pt;

    // Create an instance of the struct and populate its fields
    TrackFinderAlgorithmPsetType trackFinderPset;
    std::string prolog = "TrackFinder:"; // Adjust this as needed

    // Parse the configuration file
    try {
        // Open and read the configuration file
        std::ifstream config_file(filename);
        boost::property_tree::read_info(config_file, pt);

        trackFinderPset.MaxDTube = pt.get<int>(prolog + "MaxDTube");
        trackFinderPset.MaxDCluster = pt.get<double>(prolog + "MaxDCluster");
        trackFinderPset.SingleWireMode = pt.get<bool>(prolog + "SingleWireMode");
        trackFinderPset.MinClusterHits = pt.get<int>(prolog + "MinClusterHits");
        trackFinderPset.DCleaning = pt.get<double>(prolog + "DCleaning");
        trackFinderPset.ClusterCompletenessCut = pt.get<double>(prolog + "ClusterCompletenessCut");
        trackFinderPset.ClusterAngleCut = pt.get<double>(prolog + "ClusterAngleCut");
        trackFinderPset.CaptureMissingHits = pt.get<bool>(prolog + "CaptureMissingHits");
        trackFinderPset.MinTrackHits = pt.get<int>(prolog + "MinTrackHits");
        trackFinderPset.HitDensityThreshold = pt.get<float>(prolog + "HitDensityThreshold");
        trackFinderPset.Verbose = pt.get<int>(prolog + "Verbose");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return trackFinderPset;
}


VertexFinderAlgorithmPsetType ReadVertexFinderAlgorithmPset(std::string filename) {
    // Create a property tree to hold the parsed data
    boost::property_tree::ptree pt;

    // Create an instance of the struct and populate its fields
    VertexFinderAlgorithmPsetType vertexFinderPset;
    std::string prolog = "VertexFinder:"; // Adjust this as needed

    // Parse the configuration file
    try {
        // Open and read the configuration file
        std::ifstream config_file(filename);
        boost::property_tree::read_info(config_file, pt);

        vertexFinderPset.MaxDistToEdge = pt.get<double>(prolog + "MaxDistToEdge");
        vertexFinderPset.RefineVertexIntersection = pt.get<bool>(prolog + "RefineVertexIntersection");
        vertexFinderPset.UseEdgesDiscard = pt.get<bool>(prolog + "UseEdgesDiscard");
        vertexFinderPset.MaxTrackFractionInMain = pt.get<float>(prolog + "MaxTrackFractionInMain");
        vertexFinderPset.DecideMainTrack = pt.get<bool>(prolog + "DecideMainTrack");
        vertexFinderPset.AddCollinearLines = pt.get<bool>(prolog + "AddCollinearLines");
        vertexFinderPset.Verbose = pt.get<int>(prolog + "Verbose");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return vertexFinderPset;
}



HoughAlgorithmPsetType ReadHoughAlgorithmPset(std::string filename) {
    // Create a property tree to hold the parsed data
    boost::property_tree::ptree pt;

    // Create an instance of the struct and populate its fields
    HoughAlgorithmPsetType houghAlgorithmPset;
    std::string prolog = "Hough:"; // Adjust this as needed

    // Parse the configuration file
    try {
        // Open and read the configuration file
        std::ifstream config_file(filename);
        boost::property_tree::read_info(config_file, pt);

        houghAlgorithmPset.MaxRadiusLineHypothesis = pt.get<double>(prolog + "MaxRadiusLineHypothesis");
        houghAlgorithmPset.ThetaRes = pt.get<double>(prolog + "ThetaRes");
        houghAlgorithmPset.MaxDistanceTube = pt.get<double>(prolog + "MaxDistanceTube");
        houghAlgorithmPset.MinHoughHits = pt.get<int>(prolog + "MinHoughHits");
        houghAlgorithmPset.Verbose = pt.get<int>(prolog + "Verbose");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return houghAlgorithmPset;
}


TPCLinesAlgoPsetType ReadTPCLinesAlgoPset(std::string filename) {
    // Create a property tree to hold the parsed data
    boost::property_tree::ptree pt;

    // Create an instance of the struct and populate its fields
    TPCLinesAlgoPsetType tpcLinesAlgoPset;
    std::string prolog = "TPCLines:"; // Adjust this as needed

    // Parse the configuration file
    try {
        // Open and read the configuration file
        std::ifstream config_file(filename);
        boost::property_tree::read_info(config_file, pt);

        tpcLinesAlgoPset.MaxRadius = pt.get<double>(prolog + "MaxRadius");
        tpcLinesAlgoPset.DriftConversion = pt.get<double>(prolog + "DriftConversion");
        tpcLinesAlgoPset.MaxHoughTracks = pt.get<int>(prolog + "MaxHoughTracks");
        tpcLinesAlgoPset.MinTrackHits = pt.get<int>(prolog + "MinTrackHits");
        tpcLinesAlgoPset.RemoveIsolatedHits = pt.get<bool>(prolog + "RemoveIsolatedHits");
        tpcLinesAlgoPset.MaxNeighbourDistance = pt.get<double>(prolog + "MaxNeighbourDistance");
        tpcLinesAlgoPset.MinNeighboursHits = pt.get<int>(prolog + "MinNeighboursHits");
        tpcLinesAlgoPset.VertexAlgorithm = pt.get<int>(prolog + "VertexAlgorithm");
        tpcLinesAlgoPset.View = pt.get<std::string>(prolog + "View");
        tpcLinesAlgoPset.OutputPath = pt.get<std::string>(prolog + "OutputPath");
        tpcLinesAlgoPset.Verbose = pt.get<int>(prolog + "Verbose");
        tpcLinesAlgoPset.DebugMode = pt.get<int>(prolog + "DebugMode");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return tpcLinesAlgoPset;
}


#endif // TPC_LINES_PARAMETERS_H