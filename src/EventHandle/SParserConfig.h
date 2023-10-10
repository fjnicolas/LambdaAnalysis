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

FRAMSPsetType ReadFRANSPset(std::string filename){
   
    // Create a property tree to hold the parsed data
    boost::property_tree::ptree pt;

    // Create an instance of the struct and populate its fields
    FRAMSPsetType fPsetFRANS;
    std::string fFRAMSPsetProlog = "FRANS:";

    // Parse the FHiCL file
    try {
        // Open and read the configuration file
        std::ifstream config_file(filename);
        boost::property_tree::read_info(config_file, pt);

        fPsetFRANS.ApplyRawSmoothing = pt.get<bool>(fFRAMSPsetProlog+"ApplyRawSmoothing");
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