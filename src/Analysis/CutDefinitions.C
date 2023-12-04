#include "CutEfficienciesDefinitions.C"


//--------- Cut definitions
double fCutCRUMBS = -0.2;
double fCutMinNAngles = 1;
double fCutFRANS = 0.15;
double fCutNOrigins = 4;
double fCutNOriginsM3 = 0;
double fCutFRANSPANDORA = 0.2;
double fCutNShw = 1;
double fCutShwEnergy = 135;


std::vector<PlotDef> cutDefsBase = {
    {"TruthIsFiducial",  "0==0",            CutType::kNone,   0, {0,2,2}, true, "Truth in FV",  "No \\ cut"}
    ,{"NVertexTracks", "NVertexTracks", CutType::kCenter, 3,      {0, 8, 8},    false, "# tracks", "NTracks"}
    
    ,{"CRUMBSScore",  "0==0",                CutType::kNone,   0, {-1,1,20}, false, "CRUMBSScore",  "CRUMBSScore"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    ,{"NOriginsMult1", "0==0", CutType::kNone, 0, {0, 6, 6}, false, "# origins mult 1",  "\\# \\ origins \\ mult \\ 1"}
    ,{"NOriginsMult2", "0==0", CutType::kNone, 0, {0, 6, 6}, false, "# origins mult 2",  "\\# \\ origins \\ mult \\ 2"}
    
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS,      {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}

    /*,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"AngleGap",                "AngleGap", CutType::kLeft, 20, {0, 30, 40}, true, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  6, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHitsTrack1>5 && AngleNHitsTrack2>5", "AngleNHitsTrack1>5 && AngleNHitsTrack2>5", CutType::kRight, 1, {0, 2, 2}, false, "min(track 1 # hits, track 2 # hits) ", "Min Track1/2 hits"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 20, {0, 200, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, true, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}*/

    ,{"NOriginsMultGT3", "NOriginsMultGT3", CutType::kLeft,  fCutNOriginsM3, {0, 5, 5}  , true, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NOrigins",        "NOrigins",        CutType::kLeft,  fCutNOrigins,   {0, 15, 15}, true, "# origins",        "\\#  \\ origins"}
    
    /*,{"AngleGap",                "AngleGap", CutType::kLeft,  20, {0, 30, 40}, true, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  7.5, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHits",              "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, true, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"AngleNHitsTrack1",  "AngleNHitsTrack1",  CutType::kRight, 5, {0, 200, 40}, true, "Angle track 1 # hits", "Angle \\ track \\ 1 \\# \\ hits"}
    ,{"AngleNHitsTrack2",  "AngleNHitsTrack2",  CutType::kRight, 5, {0, 200, 40}, true, "Angle track 2 # hits", "Angle \\ track \\ 2 \\# \\ hits"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 30, {0, 200, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}*/

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    
    ,{"AngleNHitsTrack1:AngleNHitsTrack2", "0==0", CutType::k2D, 10, {0, 50, 50}, false, "Test2D", "Test2D"}

    ,{"AngleGap",                "0==0", CutType::kLeft,  20, {0, 30, 40}, false, "Gap [cm]", "Gap \\ [cm]"}
    ,{"AngleDecayContainedDiff", "0==0", CutType::kLeft,  10, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleNHitsTrack1>5 && AngleNHitsTrack2>5", "AngleNHitsTrack1>5 && AngleNHitsTrack2>5", CutType::kRight, 1, {0, 2, 2}, false, "min(track 1 # hits, track 2 # hits) ", "Min Track1/2 hits"}
    ,{"NUnassociatedHits",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    


    ,{"NVertexTracks", "NVertexTracks", CutType::kCenter, 3,      {0, 8, 8},    false, "# tracks", "NTracks"}
    ,{"NShwTh100", "NShwTh100",       CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>100MeV)", "NShower(>100MeV)"}
    ,{"NShwTh75", "NShwTh75",         CutType::kLeft, fCutNShw,      {0, 6, 6},    false, "# showers (>75MeV)", "NShower(>75MeV)"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"MainShowerEnergy", "0==0",     CutType::kLeft, 0, {0, 500, 100},  false, "MainShowerEnergy [MeV]", "MainShowerEnergy"}
    ,{"MainShowerScore", "0==0",      CutType::kLeft, 0, {0, 0.55, 20 }, false, "MainShowerScore", "MainShowerScore"}

    //,{"AngleNHits",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    //,{"AngleNHitsTrack1",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle track 1 # hits", "Angle \\ track \\ 1 \\# \\ hits"}
    //,{"AngleNHitsTrack2",              "0==0", CutType::kRight, 10, {0, 200, 40}, false, "Angle track 2 # hits", "Angle \\ track \\ 2 \\# \\ hits"}
                                                    
}; 



std::vector<PlotDef> cutDefsOriginsDistributions = {
    {"TruthIsFiducial",  "0==0",            CutType::kNone,   0, {0,2,2}, true, "Truth in FV",  "No \\ cut"}
    
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    ,{"NOrigins",        "NOrigins",        CutType::kNone,  fCutNOrigins,   {0, 5,  5}, false, "# vertices",        "\\#  \\ vertices", true}
    ,{"NOriginsMult1",    "NOriginsMult1",        CutType::kNone,  fCutNOrigins,   {0, 5, 5}, false, "# vertices multiplicity 1", "\\#  \\ vertices \\ mult \\ 1", true}
    ,{"NOriginsMult2",    "NOriginsMult2",        CutType::kNone,  fCutNOrigins,   {0, 5, 5}, false, "# vertices multiplicity 2", "\\#  \\ vertices \\ mult \\ 2", true}  
    ,{"NOriginsMultGT3", "NOriginsMultGT3", CutType::kNone,  fCutNOriginsM3, {0, 5, 5}  , true, "# vertices mult > 3", "\\# \\ vertices \\ mult \\ \\ > \\ 3", true}
};


std::vector<PlotDef> cutDefs2 = {
    {"TruthIsAV",  "TruthIsAV", CutType::kCenter, 1, {0,2,2}, true, "Truth in AV",  "Truth \\ in \\ AV"}
    ,{"NVertexTracks", "NVertexTracks", CutType::kCenter, 3,      {0, 8, 8},    false, "# tracks", "NTracks"}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    ,{"NShowerHits",               "NShowerHits",             CutType::kLeft, 20 , {0,200,50}, false, "# shower hits", "\\# \\ shower \\ hits" }
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    
    /*,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, false, "Main track # hits", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, false, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    ,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, false, "track 2 # hits ", "Track2 \\ \\# \\ hits"}*/


    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS,      {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}
    ,{"NUnOriginsMultGT3", "NUnOriginsMultGT3", CutType::kLeft,  0, {0, 5, 5}, false, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeft,  0,   {0, 15, 15}, true,  "# origins",  "\\#  \\ origins"}
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    //,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 15, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    
    
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, true, "# showers", "\\# \\ showers"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  1, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, true, "Longest is main", "Longest \\ is \\ main"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, true, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    ,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, true, "track 2 # hits ", "Track2 \\ \\# \\ hits"}
    


    ,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, false, "Longest is main", "Longest \\ is \\ main"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"AngleNHits", "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"NFreeHits",         "NFreeHits",         CutType::kLeft,  30, {0, 200, 40}, false,  "# free hits",  "\\#  \\ free \\ hits"}

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}

    
                                                    
}; 




std::vector<PlotDef> cutDefsTalk = {
    //{"TruthIsAV",  "TruthIsAV", CutType::kCenter, 1, {0,2,2}, true, "Truth in AV",  "Truth \\ in \\ AV"}
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    //,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,100}, true, "FRANS score", "FRANS \\ score"}

    
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, false, "# V", "\\# \\ V" }
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft, 1, {0, 15, 30}, false, "#Delta#alpha [#circ]", "\\Delta\\alpha\\ [\\circ]"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeftInt,  0,   {0, 15, 15}, false,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}
    
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NShowerHits", "NShowerHits", CutType::kLeft, 20, {0, 200, 20}, false, "# shower hits", "\\# \\ shower \\ hits"}

    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}

    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NShowerHits", "NShowerHits", CutType::kLeft, 20, {0, 200, 20}, false, "# shower hits", "\\# \\ shower \\ hits"}

    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}

    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft, 1, {0, 15, 30}, true, "#Delta#alpha [#circ]", "\\Delta\\alpha\\ [\\circ]"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeftInt,  0,   {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}

    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS,      {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}
    
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"NShowers", "NShowers", CutType::kLeftInt, 1, {0, 10, 10}, true, "# showers", "\\# \\ showers"}
    
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, true, "Longest is main", "Longest \\ is \\ main"}
    
   
   
   
    ,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, true, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    ,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, true, "track 2 # hits ", "Track2 \\ \\# \\ hits"}
    


    ,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeftInt,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeftInt,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, false, "Longest is main", "Longest \\ is \\ main"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"AngleNHits", "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"NFreeHits",         "NFreeHits",         CutType::kLeftInt,  30, {0, 200, 40}, false,  "# free hits",  "\\#  \\ free \\ hits"}

    
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    
}; 