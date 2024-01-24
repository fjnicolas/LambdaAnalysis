#include "CutEfficienciesDefinitions.C"


//--------- Cut definitions
double fCutCRUMBS = -0.2;
double fCutMinNAngles = 1;
double fCutFRANS = 0.2;
double fCutNOrigins = 4;
double fCutNOriginsM3 = 0;
double fCutFRANSPANDORA = 0.2;
double fCutNShw = 1;
double fCutShwEnergy = 135;


std::vector<PlotDef> phaseSpaceVars = {
    {"NuvZ",  "0==0", CutType::kNone,  0, {-10,510,50}, false, "Z_{#nu} [cm]"}
    ,{"LambdaKE", "0==0", CutType::kNone,  0, {0, 0.4, 30}, false, "KE_{#Lambda} [GeV]"}
    ,{"NuvE", "0==0", CutType::kNone,  0, {0,3, 30}, false, "E_{#nu} [GeV]"}
    ,{"Gap", "0==0", CutType::kNone, 0, {0, 20, 30}, false, "Gap [cm]"}
};


std::vector<PlotDef> cutDefsTalk = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    //,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}

    
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, false, "# V", "\\# \\ V" }
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft, 1, {0, 15, 30}, false, "#Delta#alpha [#circ]", "\\Delta\\alpha\\ [\\circ]"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeftInt,  0,   {0, 15, 15}, false,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NShowerHits", "NShowerHits", CutType::kLeft, 20, {0, 200, 20}, false, "# shower hits", "\\# \\ shower \\ hits"}

    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}

    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NShowerHits", "NShowerHits", CutType::kLeft, 20, {0, 200, 20}, false, "# shower hits", "\\# \\ shower \\ hits"}
}; 


std::vector<PlotDef> cutDefsTalk2 = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    //,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}

    

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, false, "# V", "\\# \\ V" }
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft, 1, {0, 15, 30}, false, "#Delta#alpha [#circ]", "\\Delta\\alpha\\ [\\circ]"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeftInt,  0,   {0, 15, 15}, false,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 5, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"AngleMainTrackOverlap", "AngleMainTrackOverlap", CutType::kRight,  2, {0, 50, 50}, false, "OverlapWithMainTrack", "OverlapWithMainTrack"}

    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}
    //,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}
    //,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}


    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 200, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}


    
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}

    ,{"AngleMainTrackOverlap", "AngleMainTrackOverlap", CutType::kRight,  2, {0, 50, 50}, true, "OverlapWithMainTrack", "OverlapWithMainTrack"}

    ,{"AnglePassFit", "AnglePassFit", CutType::kRight, 1, {0, 2, 2}, false, "AnglePassFit", "AnglePassFit"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, false, "AnglePassChargeFit", "AnglePassChargeFit"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 5, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, false, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kRight, 1, {0, 3, 20}, false, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kRight, 1, {0, 3, 20}, false, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kRight, 1, {0, 3, 20}, false, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleChargeRatioIntegral", "AngleChargeRatioIntegral", CutType::kNone, 1, {0, 3, 20}, false, "AngleChargeRatioIntegral", "AngleChargeRatioIntegral"}
    ,{"AngleBandCrossHits", "AngleBandCrossHits", CutType::kNone, 1, {0, 2, 50}, false, "AngleBandCrossHits", "AngleBandCrossHits"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 100, 25}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleNVertexHits", "AngleNVertexHits", CutType::kNone, 1, {0, 30, 30}, false, "AngleNVertexHits", "AngleNVertexHits"}
    ,{"AngleNBulkHits", "AngleNBulkHits", CutType::kNone, 1, {0, 100, 20}, false, "AngleNBulkHits", "AngleNBulkHits"}

    


    /*,{"AnglePassFit", "AnglePassFit", CutType::kRight, 1, {0, 2, 2}, true, "AnglePassFit", "AnglePassFit"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, true, "AnglePassChargeFit", "AnglePassChargeFit"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 5, {0, 200, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}

    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kRight, 1, {0, 3, 20}, false, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kRight, 1, {0, 3, 20}, true, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleChargeRatioIntegral", "AngleChargeRatioIntegral", CutType::kNone, 1, {0, 3, 20}, false, "AngleChargeRatioIntegral", "AngleChargeRatioIntegral"}
    ,{"AngleBandCrossHits", "AngleBandCrossHits", CutType::kNone, 1, {0, 2, 50}, false, "AngleBandCrossHits", "AngleBandCrossHits"}*/
    
    
    
    
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
    
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    
    /*,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, false, "Main track # hits", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, false, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    ,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, false, "track 2 # hits ", "Track2 \\ \\# \\ hits"}*/

    
    ,{"AnglePassFit", "AnglePassFit", CutType::kNone, 1, {0, 2, 2}, false, "AnglePassFit", "AnglePassFit"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kNone, 1, {0, 2, 2}, false, "AnglePassChargeFit", "AnglePassChargeFit"}
    ,{"AngleTwoLinesChi2", "AngleTwoLinesChi2", CutType::kNone, 1, {0, 1000, 100}, false, "Chi2", "Chi2"}
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kNone, 1, {0, 3, 20}, false, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleVertexHitIntegralDifference", "AngleVertexHitIntegralDifference", CutType::kNone, 1, {-10, 10, 100}, false, "AngleVertexHitIntegralDifference", "AngleVertexHitIntegralDifference"}
    ,{"AngleTrackLengthRatio", "AngleTrackLengthRatio", CutType::kNone, 1, {0, 3, 20}, false, "AngleTrackLengthRatio", "AngleTrackLengthRatio"}
    ,{"AngleResidualRangeMinRMS", "AngleResidualRangeMinRMS", CutType::kNone, 1, {0, 2, 40}, false, "AngleResidualRangeMinRMS", "AngleResidualRangeMinRMS"}
    ,{"AngleChargeRatioFit", "AngleChargeRatioFit", CutType::kNone, 1, {0, 3, 20}, false, "AngleChargeRatioFit", "AngleChargeRatioFit"}
    ,{"AngleChargeDifferenceFit", "AngleChargeDifferenceFit", CutType::kNone, 1, {-10, 10, 100}, false, "AngleChargeDifferenceFit", "AngleChargeDifferenceFit"}
    ,{"AngleChargeRatioIntegral", "AngleChargeRatioIntegral", CutType::kNone, 1, {0, 3, 20}, false, "AngleChargeRatioIntegral", "AngleChargeRatioIntegral"}
    ,{"AngleChargeDifferenceIntegral", "AngleChargeDifferenceIntegral", CutType::kNone, 1, {-10, 10, 100}, false, "AngleChargeDifferenceIntegral", "AngleChargeDifferenceIntegral"}
    ,{"AngleNVertexHits", "AngleNVertexHits", CutType::kNone, 1, {0, 100, 20}, false, "AngleNVertexHits", "AngleNVertexHits"}
    ,{"AngleNBulkHits", "AngleNBulkHits", CutType::kNone, 1, {0, 100, 20}, false, "AngleNBulkHits", "AngleNBulkHits"}
    ,{"AngleBandOverlap", "AngleBandOverlap", CutType::kNone, 1, {0, 2, 50}, false, "AngleBandOverlap", "AngleBandOverlap"}
    ,{"AngleBandCrossHits", "AngleBandCrossHits", CutType::kNone, 1, {0, 2, 50}, false, "AngleBandCrossHits", "AngleBandCrossHits"}


    
    //,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS,      {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}

    ,{"AnglePassFit", "AnglePassFit", CutType::kRight, 1, {0, 2, 2}, true, "AnglePassFit", "AnglePassFit"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, true, "AnglePassChargeFit", "AnglePassChargeFit"}

    ,{"NUnOriginsMultGT3", "NUnOriginsMultGT3", CutType::kLeft,  0, {0, 5, 5}, false, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeft,  0,   {0, 15, 15}, true,  "# origins",  "\\#  \\ origins"}
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    //,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 15, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}

    
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, true, "# showers", "\\# \\ showers"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  1, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    //,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, true, "Longest is main", "Longest \\ is \\ main"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    //,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, true, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    //,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, true, "track 2 # hits ", "Track2 \\ \\# \\ hits"}

    
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kRight, 1.2, {0, 3, 20}, true, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleResidualRangeMinRMS", "AngleResidualRangeMinRMS", CutType::kLeft, 0.1, {0, 2, 40}, true, "AngleResidualRangeMinRMS", "AngleResidualRangeMinRMS"}
    ,{"AngleChargeRatioFit", "AngleChargeRatioFit", CutType::kRight, 1.25, {0, 3, 20}, true, "AngleChargeRatioFit", "AngleChargeRatioFit"}
    ,{"AngleChargeRatioIntegral", "AngleChargeRatioIntegral", CutType::kRight, 1.25, {0, 3, 20}, true, "AngleChargeRatioIntegral", "AngleChargeRatioIntegral"}
    ,{"AngleNVertexHits", "AngleNVertexHits", CutType::kRight, 0, {0, 100, 20}, true, "AngleNVertexHits", "AngleNVertexHits"}   
    ,{"AngleBandCrossHits", "AngleBandCrossHits", CutType::kLeft, 0.5, {0, 2, 50}, true, "AngleBandCrossHits", "AngleBandCrossHits"}


    /*,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, false, "Longest is main", "Longest \\ is \\ main"}
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    ,{"AngleNHits", "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"NFreeHits",         "NFreeHits",         CutType::kLeft,  30, {0, 200, 40}, false,  "# free hits",  "\\#  \\ free \\ hits"}*/

    

    
                                                    
}; 




