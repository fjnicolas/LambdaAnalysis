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


PlotDef startCutForLoop = {"TruthIsFiducial", "TruthIsFiducial", CutType::kCenter, 1,   {0, 2, 2}, true,  "fiducial volume",  "TruthFiducialVolume"};
PlotDef minimalCutForLoop = {"RecoIsFiducial==1 && NAngles>=1 && AngleFRANSScore>0.2 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(RecoIsFiducial==1 && NAngles>=1 && AngleFRANSScore>0.2 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "minimal cut",  "MinimalCut"};


std::vector<PlotDef> cutDefsTalk2 = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    // FRANS (all events)
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    
    // Topological cuts
    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}

    // FRANS
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}

    // CLEANESS
    ,{"AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=50", "(AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=50)", CutType::kCenter, 1, {0, 2, 2}, true, "PreCalorimetryCut", "PreCalorimetryCut"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "AnglePzSign", "AnglePzSign"}  
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleDirtHits<=10 && NUnassociatedHits<=10", "(AngleDirtHits<=10 && NUnassociatedHits<=10)", CutType::kCenter, 1, {0, 2, 2}, true, "# unnassociated hits<10", "unnassociated hits<10"}
    //,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    // After Cuts Distributions
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, false, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  25, {0, 50, 25}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, false, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, false, "AnglePassChargeFit", "AnglePassChargeFit"}
    
}; 


std::vector<PlotDef> cutDefsPID = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}

    //,{"NVertexTracks", "NVertexTracks", CutType::kLeft, 3, {0, 10, 10}, true, "# vertex tracks", "NVertexTracks"}
    
    // Topological cuts
    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=10", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=10)", CutType::kCenter,  1, {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}

    //,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, 0.2, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}


    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, false, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, false, "AnglePzSign", "AnglePzSign"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  25, {0, 50, 25}, false, "#Dirt Hits", "Dirt \\ Hits"}

    //PID
    ,{"NTracksLI", "NTracksLI", CutType::kRight, 1, {0, 5, 5}, false, "# tracks LI", "\\# \\ tracks \\ LI"}
    ,{"NTracksHI", "NTracksHI", CutType::kRight, 1, {0, 5, 5}, false, "# tracks HI", "\\# \\ tracks \\ HI"}
    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "NTracksLI>0 && NTracksHI>0", "NTracksLI>0 \\ \\&\\& \\ NTracksHI>0"}


    //,{"AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=50", "(AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=50)", CutType::kCenter, 1, {0, 2, 2}, true, "PreCalorimetryCut", "PreCalorimetryCut"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "AnglePzSign", "AnglePzSign"}  
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleDirtHits<=10 && NUnassociatedHits<=10", "(AngleDirtHits<=10 && NUnassociatedHits<=10)", CutType::kCenter, 1, {0, 2, 2}, true, "# unnassociated hits<10", "unnassociated hits<10"}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, true, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    //,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    //,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    
    // After Cuts Distributions
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, false, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  25, {0, 50, 25}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, false, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, false, "AnglePassChargeFit", "AnglePassChargeFit"}


    ,{"InvariantMass", "InvariantMass", CutType::kCenter, 1115, {1000, 1500, 10}, false, "Invariant Mass [MeV]", "Invariant \\ Mass"}

    
}; 



//--------- CUT REPORITORY DEFINITIONS ----------------
std::vector<PlotDef> cutDefsRepository = {

    // Generic
    {"TruthIsAV",  "TruthIsAV", CutType::kCenter, 1, {0,2,2}, true, "Truth in AV",  "Truth \\ in \\ AV"}
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    
    // # origins
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    ,{"NUnOriginsMultGT3", "NUnOriginsMultGT3", CutType::kLeft,  0, {0, 5, 5}, false, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeft,  0,   {0, 15, 15}, true,  "# origins",  "\\#  \\ origins"}

    // FRANS
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}
    

    // CLEANESS
    ,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, false, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    ,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, false, "track 2 # hits ", "Track2 \\ \\# \\ hits"}
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, false, "Longest is main", "Longest \\ is \\ main"}
    ,{"AngleNHits", "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"NFreeHits",         "NFreeHits",         CutType::kLeft,  30, {0, 200, 40}, false,  "# free hits",  "\\#  \\ free \\ hits"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 15, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    
    // MAIN TRACK
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, false, "Main track # hits", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleMainTrackOverlap", "AngleMainTrackOverlap", CutType::kRight,  2, {0, 50, 50}, false, "OverlapWithMainTrack", "OverlapWithMainTrack"}
    ,{"AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion", CutType::kLeft,  0.1, {0, 1, 10}, true, "AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion"}

    // 2D CALORIMETRY VARIABLES
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


    // SHOWERS
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    
    // 3D CALORIMETRY VARIABLES

    // TOPOLIGICAL CUTS
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}
    ,{"NVertexTracks", "NVertexTracks", CutType::kLeft, 3, {0, 8, 8}, false, "# vertex tracks", "NVertexTracks"}
};


// --------- CUTS TO PLOT ORIGINS DISTRIBUTIONS ------------
std::vector<PlotDef> cutDefsOriginsDistributions = {
    {"TruthIsFiducial",  "0==0",            CutType::kNone,   0, {0,2,2}, true, "Truth in FV",  "No \\ cut"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    ,{"NOrigins",        "NOrigins",        CutType::kNone,  fCutNOrigins,   {0, 5,  5}, false, "# vertices",        "\\#  \\ vertices", true}
    ,{"NOriginsMult1",    "NOriginsMult1",        CutType::kNone,  fCutNOrigins,   {0, 5, 5}, false, "# vertices multiplicity 1", "\\#  \\ vertices \\ mult \\ 1", true}
    ,{"NOriginsMult2",    "NOriginsMult2",        CutType::kNone,  fCutNOrigins,   {0, 5, 5}, false, "# vertices multiplicity 2", "\\#  \\ vertices \\ mult \\ 2", true}  
    ,{"NOriginsMultGT3", "NOriginsMultGT3", CutType::kNone,  fCutNOriginsM3, {0, 5, 5}  , true, "# vertices mult > 3", "\\# \\ vertices \\ mult \\ \\ > \\ 3", true}
};


// --------- CUTS FOR BRAZIL CM ------------
std::vector<PlotDef> cutDefsTalk = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

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


// --------- CUTS FOR LOOP MACRO ------------
std::vector<PlotDef> cutDefsLoop{
    {"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleDirtHits<=10 && NUnassociatedHits<=10", "(AngleDirtHits<=10 && NUnassociatedHits<=10)", CutType::kCenter, 1, {0, 2, 2}, true, "# unnassociated hits<10", "unnassociated hits<10"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "AnglePzSign", "AnglePzSign"}  
    ,{"AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=45", "(AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=45)", CutType::kCenter, 1, {0, 2, 2}, true, "PreCalorimetryCut", "PreCalorimetryCut"}
    //,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 45, {0, 360, 36}, true, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    //,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, true, "AnglePassChargeFit", "AnglePassChargeFit"}
};