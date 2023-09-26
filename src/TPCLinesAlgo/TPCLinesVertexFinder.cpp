////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesHough
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesVertexFinder.h"


// ------- constructor
TPCLinesVertexFinder::TPCLinesVertexFinder(VertexFinderAlgorithmPsetType tpcLinesVertexFinderPset):
    fTPCLinesVertexFinderPset(tpcLinesVertexFinderPset)
{}

double TPCLinesVertexFinder::GetAngle360(double x, double y) {
    double a = std::atan(std::abs(y / x)) * 180.0 / M_PI;

    if (x > 0 && y < 0) { // quadrant 4
        a = 360 - a;
    } else if (x < 0 && y < 0) { // quadrant 3
        a = 180 + a;
    } else if (x < 0 && y > 0) { // quadrant 2
        a = 180 - a;
    }

    return a;
}


bool TPCLinesVertexFinder::TrackTriangleJunctionConatined(SLinearCluster track, STriangle tri){
    
    SPoint p1(track.GetMinX(), track.GetYatMinX());
    double d1 = TPCLinesDistanceUtils::GetHitDistance(p1, tri.GetMainVertex());
    SPoint p2(track.GetMaxX(), track.GetYatMaxX());
    double d2 = TPCLinesDistanceUtils::GetHitDistance(p2, tri.GetMainVertex());
    SPoint mainTrackVertex = (d1 < d2) ? p1 : p2;

    double juntionDirection[] = {
        tri.GetMainVertex().X() - mainTrackVertex.X(),
        tri.GetMainVertex().Y() - mainTrackVertex.Y()
    };

    double VDir1[] = {
        tri.GetVertexB().X() - tri.GetMainVertex().X(),
        tri.GetVertexB().Y() - tri.GetMainVertex().Y()
    };

    double VDir2[] = {
        tri.GetVertexC().X() - tri.GetMainVertex().X(),
        tri.GetVertexC().Y() - tri.GetMainVertex().Y()
    };

    double juntionDirectionAngle = GetAngle360(juntionDirection[0], juntionDirection[1]);
    double VDir1Angle = GetAngle360(VDir1[0], VDir1[1]);
    double VDir2Angle = GetAngle360(VDir2[0], VDir2[1]);

    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "V ANGLES: " << juntionDirectionAngle << ", " << VDir1Angle << ", " << VDir2Angle << std::endl;


    double minAngle = std::min(VDir1Angle, VDir2Angle);
    double maxAngle = std::max(VDir1Angle, VDir2Angle);

    double fExtraAngle = 5.0;
    bool junctionContained = false;
    
    if (maxAngle - minAngle < 180) {
        junctionContained = (minAngle-fExtraAngle<=juntionDirectionAngle) && (juntionDirectionAngle<=maxAngle+fExtraAngle);
    }
    else {
        junctionContained = (juntionDirectionAngle<=minAngle+fExtraAngle) || (juntionDirectionAngle>=maxAngle-fExtraAngle);
    }

    if (!junctionContained) {
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "JUNCTION NOT CONTAINED" << std::endl;
    }
    else{
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "JUNCTION CONTAINED" << std::endl;
    }

    return junctionContained;
}


int TPCLinesVertexFinder::GetNHitsBetweenJunction(SLinearCluster track, STriangle tri, std::vector<SLinearCluster> trackList, SPoint intP, 
std::vector<int> track1ListIx, std::vector<int> track2ListIx, double tol){// std::vector<SLinearCluster> track2List, std::vector<SLinearCluster> track1List, std::vector<SLinearCluster> track2List, ){

    SPoint p1(track.GetMinX(), track.GetYatMinX());
    double d1 = TPCLinesDistanceUtils::GetHitDistance(p1, tri.GetMainVertex());
    SPoint p2(track.GetMaxX(), track.GetYatMaxX());
    double d2 = TPCLinesDistanceUtils::GetHitDistance(p2, tri.GetMainVertex());
    SPoint mainTrackVertex = (d1 < d2) ? p1 : p2;

    double juntionDirection[] = {
        tri.GetMainVertex().X() - mainTrackVertex.X(),
        tri.GetMainVertex().Y() - mainTrackVertex.Y()
    };

    double juntionEdges[2] = { std::min(mainTrackVertex.X(), tri.GetMainVertex().X()), std::max(mainTrackVertex.X(), tri.GetMainVertex().X()) };
    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Edges: " << juntionEdges[0] << ", " << juntionEdges[1] << std::endl;

    double juncSlope = juntionDirection[1] / juntionDirection[0];
    double juncIntercept = intP.Y() - juncSlope * intP.X();

    int nHitsInMiddle = 0;
    for (SLinearCluster & eTrack : trackList) {
        if( std::find(track1ListIx.begin(), track1ListIx.end(), eTrack.GetId())!=track1ListIx.end()) continue;
        if( std::find(track2ListIx.begin(), track2ListIx.end(), eTrack.GetId())!=track2ListIx.end()) continue;
        
        //if (eTrack.GetId() == track1.GetId() || eTrack.GetId() == track2.GetId()) continue;
        
        if ((eTrack.GetMinX() > juntionEdges[1] && eTrack.GetMaxX() > juntionEdges[1]) ||
            (eTrack.GetMinX() < juntionEdges[0] && eTrack.GetMaxX() < juntionEdges[0])) {
            continue;
        }
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "In the middle track " << eTrack.GetId() << std::endl;
        
        for (SHit& hit : eTrack.GetHits()) {
            double yHypo = juncSlope * hit.X() + juncIntercept;
            if (std::abs(yHypo - hit.Y()) < tol*hit.Width()) {
                nHitsInMiddle++;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   In middle hit: " << hit.X() << ", " << hit.Y() << std::endl;
            }
        }
    }

    return nHitsInMiddle;

}


bool TPCLinesVertexFinder::areVectorsEqual(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        if (vec1[i] != vec2[i]) {
            return false;
        }
    }

    return true;
}


SPoint TPCLinesVertexFinder::GetTracksIntersection(SLinearCluster track1, SLinearCluster track2, double dMax, bool useEdgeSlopes, bool useFit){
    SPoint intP(-1, -1);

    if (useFit) {
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Using fit not implemented\n";
        // Implement the useFit logic here.
        // The code related to evaluating track fits and finding intersects goes here.
    }
    else {
        LineEquation lineEq1 = track1.GetTrackEquation();
        LineEquation lineEq2 = track2.GetTrackEquation();
        double xInt = (lineEq1.Intercept() - lineEq2.Intercept()) / (lineEq2.Slope() - lineEq1.Slope());
        double yInt = lineEq1.Slope() * xInt + lineEq1.Intercept();

        intP = SPoint(xInt, yInt);

        if (useEdgeSlopes == true) {
            if (std::abs(intP.X() - track1.GetMinX()) < std::abs(intP.X() - track1.GetMaxX())) {
                lineEq1 = track1.GetTrackEquationStart();
            } else {
                lineEq1 = track1.GetTrackEquationEnd();
            }

            if (std::abs(intP.X() - track2.GetMinX()) < std::abs(intP.X() - track2.GetMaxX())) {
                lineEq2 = track2.GetTrackEquationStart();
            } else {
                lineEq2 = track2.GetTrackEquationEnd();
            }

            xInt = (lineEq1.Intercept() - lineEq2.Intercept()) / (lineEq2.Slope() - lineEq1.Slope());
            yInt = lineEq1.Slope() * xInt + lineEq1.Intercept();
            intP = SPoint(xInt, yInt);
        }
    }

    SHit intHit = SHit(intP.X(), intP.Y());
    double d1 = track1.GetHitCluster().GetMinDistanceToCluster(intHit);
    double d2 = track2.GetHitCluster().GetMinDistanceToCluster(intHit);

    if (d1 < dMax || d2 < dMax) {
        return intP;
    } else {
        return SPoint(-1, -1); // Return an appropriate "no intersection" value.
    }
}


int TPCLinesVertexFinder::GetHitsContainedInLineEquation(LineEquation trackEq, std::vector<SHit> hitList, float tol) {
    
    int nhits = 0;

    for (SHit& hit : hitList){
        float yHypo = trackEq.Slope() * hit.X() + trackEq.Intercept();
        if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
            nhits++;
            //if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << hit.X() << " " << hit.Y() << " " << yHypo << std::endl;
        }
    }

    return nhits;
}


// GetHitsContainedInHypo
// Input two tracks, get the NHits closer to a vertex
// Check how many hits of each track are contained in the line equation of the other track
std::vector<SHit> TPCLinesVertexFinder::GetMutualHitsContainedInHypo(SLinearCluster track1, SLinearCluster track2,  SPoint intP, int nHits, float tol) {
    
    // Get hits closest to the itneraction point
    int maxHits = std::min(nHits, track1.NHits());
    std::vector<SHit> track1Hits_ = track1.GetHits();
    std::vector<SHit> track1Hits;
    if (std::abs(intP.X() - track1.GetMinX()) < std::abs(intP.X() - track1.GetMaxX())) {
        track1Hits = std::vector<SHit>(track1Hits_.begin(), std::next(track1Hits_.begin(), maxHits) ); 
    }
    else {
        track1Hits = std::vector<SHit>(std::prev(track1Hits_.end(), maxHits), track1Hits_.end());
    }
    SLinearCluster trk1(track1Hits); 

    maxHits = std::min(nHits, track2.NHits());
    std::vector<SHit> track2Hits;
    std::vector<SHit> track2Hits_ = track2.GetHits();
    if (std::abs(intP.X() - track2.GetMinX()) < std::abs(intP.X() - track2.GetMaxX())) {
        track2Hits = std::vector<SHit>(track2Hits_.begin(), std::next(track2Hits_.begin(), maxHits) ); 
    }
    else {
        track2Hits = std::vector<SHit>(std::prev(track2Hits_.end(), maxHits), track2Hits_.end());
    }
    SLinearCluster trk2(track2Hits); 

   
    LineEquation track1Eq = trk1.GetTrackEquation();
    
    LineEquation track2Eq = trk2.GetTrackEquation();

    std::vector<SHit> containedHits;

    for (SHit& hit : track2.GetHits()) {
        float yHypo = track1Eq.Slope() * hit.X() + track1Eq.Intercept();
        if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
            containedHits.push_back(hit);
        }
    }

    for (SHit& hit : track1.GetHits()) {
        float yHypo = track2Eq.Slope() * hit.X() + track2Eq.Intercept();
        if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
            containedHits.push_back(hit);
        }
    }

    return containedHits;
}


SPoint TPCLinesVertexFinder::check_arrow_line_intersection(float Ax, float Ay, float Dx, float Dy, float line_slope, float line_intercept) {
    float m = line_slope, b = line_intercept;

    // Calculate t from the intersection equation
    float t = (m * Ax - Ay + b) / (Dy - m * Dx);

    float x = Ax + t * Dx;
    float y = Ay + t * Dy;

    SPoint intersection_point(x, y);

    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "ArrowDir intP " << intersection_point.X() << " " << intersection_point.Y() << std::endl;

    if ((x - Ax >= 0) == (Dx >= 0)) {
        return intersection_point;
    } else {
        return SPoint(-1, -1);
    }
}


SPoint TPCLinesVertexFinder::GetTrackssEuationOppositePoint(SLinearCluster track, std::vector<SLinearCluster> trackList, SPoint p){
    float minX1 = 1e3;
    float maxX1 = 0;
    for (SLinearCluster trk : trackList) {
        minX1 = std::min(minX1, trk.GetMinX());
        maxX1 = std::max(maxX1, trk.GetMaxX());
    }

    std::vector<float> Xpoints;
    Xpoints.push_back( std::min(p.X(), minX1) );
    Xpoints.push_back( std::max(p.X(), maxX1) );

    double m1 = track.GetTrackEquation().Slope();
    double n1 = track.GetTrackEquation().Intercept();

    std::vector<float> Ypoints;
    for (auto x : Xpoints) {
        Ypoints.push_back(m1 * x + n1);
    }

    SPoint pA(Xpoints[0], Ypoints[0]);
    SPoint pB(Xpoints[1], Ypoints[1]);

    double dA = TPCLinesDistanceUtils::GetHitDistance(pA, p);
    double dB = TPCLinesDistanceUtils::GetHitDistance(pB, p);


    return  (dA > dB) ? pA : pB;
}


std::vector<SLinearCluster> TPCLinesVertexFinder::GetCollinearTracks(SLinearCluster mainTrack, std::vector<SLinearCluster> trackList){

    std::vector<SLinearCluster> pointingTracksV;


    // Get directions pointing to the main track vertex
    for(auto &trk:trackList){
        // intersection point with the main track
        SPoint intP = GetTracksIntersection(trk, mainTrack, 10000);

        // check the intersection point is close enough to the mian track vertex
        double d = TPCLinesDistanceUtils::GetHitDistance(mainTrack.GetStartPoint(), intP);
        if(d<10){
            pointingTracksV.push_back(trk);
        }
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   Track "<<trk.GetId()<<" D="<<d<<" IntP: "<<intP;
    }

    // sort by minX
    std::sort(pointingTracksV.begin(), pointingTracksV.end(), [&](SLinearCluster& l1, SLinearCluster& l2) {return l1.GetMinX() < l2.GetMinX();} );    



    std::vector<std::vector<SLinearCluster>> groupedTracks;
    if(pointingTracksV.size()>0){
        // get consecutive segments
        groupedTracks.push_back({pointingTracksV[0]});
        SLinearCluster lastCluster = pointingTracksV[0];
        if(pointingTracksV.size()>=2){
            for (size_t i = 1; i < pointingTracksV.size(); ++i) {
                SLinearCluster currentCluster = pointingTracksV[i];
                
                // Consider the track if thje start and end are within 1 wire
                bool collinear = false;
                if ( std::abs(currentCluster.GetMinX() - lastCluster.GetMaxX()) <= 2){
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   Contiguous "<<currentCluster.GetId()<<" to "<<lastCluster.GetId()<<std::endl;
                    // check is the edge hit for both tracks is compatible with the connectedness
                    SPoint endHit = lastCluster.GetEndPoint();
                    SPoint startHit = currentCluster.GetStartPoint();

                    double d12 = TPCLinesDistanceUtils::GetHitDistance(endHit, startHit);

                    double averageComp = (currentCluster.GetCompactness()+lastCluster.GetCompactness())/2.;
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    AvComp: "<<averageComp<<" D12: "<<d12<<std::endl;
                    if(d12<5.*averageComp){
                        collinear=true;
                    }
                }

                if(collinear){
                    groupedTracks.back().push_back(currentCluster);
                }
                else{
                    groupedTracks.push_back({currentCluster});
                }

                lastCluster = currentCluster;

            }
         }
    }

    std::vector<SLinearCluster> finalTrackList;
    for (size_t i = 0; i < groupedTracks.size(); ++i) {
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Grouped tracks "<<i<<std::endl;
        std::vector<SHit> hits;
        for (size_t j = 0; j < groupedTracks[i].size(); ++j) {
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Trk: "<<groupedTracks[i][j].GetId()<<std::endl;
            std::vector<SHit> trk_hits = groupedTracks[i][j].GetHits();
            hits.insert(hits.end(), trk_hits.begin(), trk_hits.end());
        }

        int nHits = std::min(hits.size(), (size_t)5);

        std::sort(hits.begin(), hits.end(), [](SHit& a, SHit& b) {return a.X() < b.X();});
        std::vector<SHit> hitsStart = hits;
        hitsStart.resize(nHits);

        std::sort(hits.begin(), hits.end(), [](SHit& a, SHit& b) {return a.X() > b.X();});
        std::vector<SHit> hitsEnd = hits;
        //hitsEnd.resize(nHits);

        double avIntStart=0;
        double avIntEnd=0;
        for(int k=0; k<nHits; k++){
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"S:"<<hitsStart[k].Integral()<<hitsStart[k]<<" E:"<<hitsEnd[k].Integral()<<hitsEnd[k]<<std::endl;
            avIntStart+=hitsStart[k].Integral();
        }
        avIntStart/=nHits;

        for(size_t k=0; k<hitsEnd.size(); k++){
            avIntEnd+=hitsEnd[k].Integral();
        }
        
        avIntEnd/=hitsEnd.size();

        SLinearCluster finalTrack(hits);
        

        // get average hit signal at the beginning

        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<< "  Final track "<<finalTrack.GetMinX()<<" "<<finalTrack.GetMaxX()<<" Integral start/end "<<avIntStart<<" "<<avIntEnd<<std::endl;

        // add ig the energy deposition is larger at the beginning and the main track is larger
        if(avIntStart>1.25*avIntEnd && hits.size()){//<mainTrack.NHits()){
            finalTrackList.push_back(finalTrack);
        }
    }
  
    
                    
    return finalTrackList;
} 


SLinearCluster TPCLinesVertexFinder::GetMainDirection(std::vector<SLinearCluster> & mainDirTrackList, 
                                std::vector<SLinearCluster>& freeTracksList,
                                std::vector<std::vector<SLinearCluster>> parallelTracks,
                                bool decideMainTrack,
                                int verbose)
{
    
    // longest
    std::vector<SLinearCluster> LongDirTrackList, LongFreeTracksList;
    SLinearCluster LongDirection = TPCLinesDirectionUtils::GetMainDirectionLongest(parallelTracks, LongDirTrackList, LongFreeTracksList);
    // downstream
    std::vector<SLinearCluster> DownDirTrackList, DownFreeTracksList;
    SLinearCluster DownDirection = TPCLinesDirectionUtils::GetMainDirectionDownstream(parallelTracks, DownDirTrackList, DownFreeTracksList);

    std::vector<int> LongDirectionIndexes;
    if(verbose>=1) std::cout << "Longest track: ";
    for (SLinearCluster &track : LongDirTrackList) {
        if(verbose>=1) std::cout << track.GetId() << " ";
        LongDirectionIndexes.push_back(track.GetId());
    }
    if(verbose>=1) std::cout << std::endl;

    if(verbose>=1) std::cout << "Downstream track: ";
    std::vector<int> DownDirectionIndexes;
    for (SLinearCluster &track : DownDirTrackList) {
        if(verbose>=1) std::cout << track.GetId() << " ";
        DownDirectionIndexes.push_back(track.GetId());
    }
    if(verbose>=1) std::cout << std::endl;

    bool downIsLong = areVectorsEqual(DownDirectionIndexes, LongDirectionIndexes);
    bool useLargest = true;

    if(decideMainTrack==true){
        if (!downIsLong) {
            double connTol = 3 * (LongDirection.GetConnectedness() + DownDirection.GetConnectedness()) / 2;
            double conn = TPCLinesDistanceUtils::GetClusterConnectedness(LongDirection.GetHitCluster(), DownDirection.GetHitCluster());
            SPoint intP = GetTracksIntersection(LongDirection, DownDirection, 20, true);

            if(verbose>=1) std::cout << "Long/main direction intersection: " << intP << " Connectedness: " << conn << " ConnTol: " << connTol << std::endl;

            if (conn > connTol) {
                useLargest = false;
                if(verbose>=1) std::cout << "  longest and most downstream are not connected, using the most downstream" << std::endl;
            }
            else {
                // if they are consecutive, still use the most downstream
                if(DownDirection.GetMaxX()<=LongDirection.GetMinX()){
                    useLargest = false;
                    if(verbose>=1) std::cout << "  longest is consecutive to the downstream, using the most downstream" << std::endl;
                }
                if(verbose>=1) std::cout << " connected, using the longest" << std::endl;
            }
        }
        else {
            if(verbose>=1) std::cout << "Longest/most downstream directions are the same" << std::endl;
        }
    }

    SLinearCluster MainDirection;
    mainDirTrackList.clear();
    freeTracksList.clear();
    if (!useLargest) {
        if(verbose>=1) std::cout << "Using the most downstream track" << std::endl;
        MainDirection = DownDirection;
        mainDirTrackList = DownDirTrackList;
        freeTracksList = DownFreeTracksList;
    }
    else{
        MainDirection = LongDirection;
        mainDirTrackList = LongDirTrackList;    
        freeTracksList = LongFreeTracksList;
    }

    return MainDirection;
}


void TPCLinesVertexFinder::GetOrigins(std::vector<SLinearCluster> trackList, std::vector<STriangle>& vertexList, std::vector<SPoint> &originList, SLinearCluster &mainDirection){

    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" In Origin finder\n";    

    // ------- reset variables
    vertexList.clear();
    originList.clear();


    // ------ Get the short/vertex tracks
    std::map<int, std::vector<int>> shortToLongTrackDict;
    std::map<int, int> shortConnectionTrackDict;
    std::vector<SLinearCluster> shortTrackList = TPCLinesDirectionUtils::GetVertexTracks(trackList, shortToLongTrackDict, shortConnectionTrackDict, 6, 3, 5, fTPCLinesVertexFinderPset.Verbose); 
    // remove the short tracks
    std::vector<SLinearCluster> newTrackList;
    for (SLinearCluster& track : trackList) {
        if (shortToLongTrackDict.find(track.GetId()) == shortToLongTrackDict.end()) {
            newTrackList.push_back(track);
        }
    }


    // ------ Get the parallel tracks
    std::vector<std::vector<SLinearCluster>> parallelTracks = TPCLinesDirectionUtils::GetParallelTracks(newTrackList, -2, 15, 30, fTPCLinesVertexFinderPset.Verbose);

    std::vector<std::vector<int>> parallelTrackClusterIndexList;
    for(size_t ix = 0; ix<parallelTracks.size(); ix++){
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Parallel track cluster "<<ix<<": ";
        std::vector<int> _indexes;
        for(size_t jx = 0; jx<parallelTracks[ix].size(); jx++){
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<parallelTracks[ix][jx].GetId()<<" ";
            _indexes.push_back(parallelTracks[ix][jx].GetId());
        }
        parallelTrackClusterIndexList.push_back(_indexes);
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<std::endl;
    }


    // ------ Get main tracks
    std::vector<SLinearCluster> MainDirTrackList;
    std::vector<SLinearCluster> FreeTracksList;
    SLinearCluster MainDirection = GetMainDirection(MainDirTrackList, FreeTracksList, parallelTracks, fTPCLinesVertexFinderPset.DecideMainTrack, fTPCLinesVertexFinderPset.Verbose);


    // ------- Look for possible intersections
    // ------ Loop 1
    if (FreeTracksList.size() > 0) {
        for (size_t ix = 0; ix < FreeTracksList.size(); ++ix) {
            SLinearCluster track1 = FreeTracksList[ix];


            // Loop through other tracks
            // ------ Loop 2
            for (size_t jx = ix + 1; jx < FreeTracksList.size(); ++jx) {
                SLinearCluster track2 = FreeTracksList[jx];

                // ------ get the associated parallel tracks
                std::vector<SLinearCluster> track1List;
                std::vector<SLinearCluster> track2List;
                std::vector<int> track1ListIds;
                std::vector<int> track2ListIds;
                for (size_t ix = 0; ix < parallelTrackClusterIndexList.size(); ++ix) {
                    std::vector<int> parTrackCluster = parallelTrackClusterIndexList[ix];
                    if (std::find(parTrackCluster.begin(), parTrackCluster.end(), track1.GetId()) != parTrackCluster.end()) {
                        track1List = parallelTracks[ix];
                        for(SLinearCluster&trk:track1List){
                            track1ListIds.push_back(trk.GetId());
                        }
                    }
                    else if (std::find(parTrackCluster.begin(), parTrackCluster.end(), track2.GetId()) != parTrackCluster.end())
                    {
                        track2List = parallelTracks[ix];
                        for(SLinearCluster&trk:track2List){
                            track2ListIds.push_back(trk.GetId());
                        }
                    }
                }


                // ------ Skip if its in the same parallel cluster as track 1
                if (std::find(track1ListIds.begin(), track1ListIds.end(),
                            track2.GetId()) != track1ListIds.end()) {
                    continue;
                }


                // ------ Check if the tracks are connected
                // calculate connection based on the tracks connectedes
                float connTol = 6 * (track1.GetConnectedness() + track2.GetConnectedness()) / 2;
                float minConn = 1e3;
                for (SLinearCluster trk1 : track1List) {
                    for (SLinearCluster trk2 : track2List) {
                        float conn = TPCLinesDistanceUtils::GetClusterConnectedness(trk1.GetHitCluster(), trk2.GetHitCluster());
                        if (conn < minConn) {
                            minConn = conn;
                        }
                    }
                }

                //  calculate connection based on short tracks
                bool connectedThroughShortTrack = false;
                if(shortConnectionTrackDict.find(track1.GetId()) != shortConnectionTrackDict.end()){
                    if(shortConnectionTrackDict[track1.GetId()] == track2.GetId()){
                        connectedThroughShortTrack=true;
                    }
                }

                bool connected = (minConn < connTol) || connectedThroughShortTrack;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "\n  -- Potential intersection " << track1.GetId() << " " << track2.GetId() << " MinConn:" << minConn << " Tol:s" << connTol << std::endl;
                if (!connected) {
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" ... not connnected\n";
                    continue;
                }


                // ------ Look for the intersection points
                SPoint intP = GetTracksIntersection(track1, track2, 50, fTPCLinesVertexFinderPset.RefineVertexIntersection);
                if(intP.X()==-1 and intP.Y()==-1) continue;

                // Get closest hits and vertex hits
                std::pair<SHit, double> cloHit1Pair = track1.GetHitCluster().GetClosestHitToPoint(intP);
                SHit cloHit1 = cloHit1Pair.first;
                double dHit1 = cloHit1Pair.second;
                std::pair<SHit, double> cloHit2Pair = track2.GetHitCluster().GetClosestHitToPoint(intP);
                SHit cloHit2 = cloHit2Pair.first;
                double dHit2 = cloHit2Pair.second;
                SHit cloHit = (dHit1 < dHit2) ? cloHit1 : cloHit2;
                //double dHit = (dHit1 < dHit2) ? dHit1 : dHit2;

                // Study the edges of the tracks
                std::vector<SHit> vertexHits = GetMutualHitsContainedInHypo(track1, track2, intP, 5);
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "    NVERTEX HITS " << vertexHits.size() << std::endl;

                float dEdge1 = std::min(std::abs(cloHit1.X() - track1.GetMinX()), std::abs(cloHit1.X() - track1.GetMaxX()));
                float dEdge2 = std::min(std::abs(cloHit2.X() - track2.GetMinX()), std::abs(cloHit2.X() - track2.GetMaxX()));

                for (SHit& hit : vertexHits) {
                    float min1 = std::min(std::abs(hit.X() - track1.GetMinX()), std::abs(hit.X() - track1.GetMaxX()));
                    float min2 = std::min(std::abs(hit.X() - track2.GetMinX()), std::abs(hit.X() - track2.GetMaxX()));
                    
                    if (min1 < dEdge1) {
                        dEdge1 = min1;
                    }
                    if (min2 < dEdge2) {
                        dEdge2 = min2;
                    }
                }

                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Clo Hit1/2 " << cloHit1 << " " << cloHit2 << " DEdges: " << dEdge1 << " " << dEdge2 << std::endl;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Clo Hit " << cloHit << std::endl;

                if (fTPCLinesVertexFinderPset.UseEdgesDiscard && (dEdge2 >= fTPCLinesVertexFinderPset.MaxDistToEdge || dEdge1 >= fTPCLinesVertexFinderPset.MaxDistToEdge)) {
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Skipping...intersection is not with edge hits" << std::endl;
                    continue;
                }

                float minD = 1e3;
                SHit vertexHit(intP.X(), intP.Y());
                for (SHit& hit : vertexHits) {
                    float d = std::hypot(hit.X() - intP.X(), hit.Y() - intP.Y());
                    if (d < minD) {
                        minD = d;
                        vertexHit = hit;
                    }
                }

                //  put intersection point in the closest hit
                intP = SPoint(vertexHit.X(), vertexHit.Y());
                bool thereIsIntersectionVertex = (intP.X()!=-1 && intP.Y()!=-1);

                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << "      Vertex set to: " << intP << std::endl;
                    std::cout << "V TRACKS" << std::endl;
                    for (SLinearCluster& t : track1List) {
                        std::cout << t.GetId() << " ";
                    }
                    std::cout << std::endl;
                    for (SLinearCluster& t : track2List) {
                        std::cout << t.GetId() << " ";
                    }
                    std::cout << std::endl;
                }
                

                if(thereIsIntersectionVertex){
                    
                    // ------ Create the triangle object
                    SPoint vertex1 = GetTrackssEuationOppositePoint( track1, track1List, intP );
                    SPoint vertex2 = GetTrackssEuationOppositePoint( track2, track2List, intP );
                    STriangle triangle = STriangle(intP, vertex1, vertex2, cloHit, track1.GetIntegral(), track2.GetIntegral());

                    // CHECK 1: check if the triangle direction intersects the main driection
                    // get closes segment in the MainDirection to the intersection point
                    double minD = 1e3;
                    SLinearCluster mainTrackDir = MainDirection;

                    for (SLinearCluster mTrack : MainDirTrackList) {
                        double d = mTrack.GetHitCluster().GetMinDistanceToCluster(intP);
                        if (d < minD) {
                            minD = d;
                            mainTrackDir = mTrack;
                        }
                    }

                    SPoint start_point(triangle.GetMainVertex().X(), triangle.GetMainVertex().Y());
                    SPoint direction_vector = triangle.GetDirectorVector();
                    double line_slope = mainTrackDir.GetTrackEquation().Slope();
                    double line_intercept = mainTrackDir.GetTrackEquation().Intercept();
                    
                    SPoint intersection_point = check_arrow_line_intersection(start_point.X(), start_point.Y(),direction_vector.X(), direction_vector.Y(), line_slope, line_intercept);
                    bool triangleIntersects = (intersection_point.X()!=-1 && intersection_point.Y()!=-1);
                    
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Arrow line intersects the line equation at: (" << intersection_point.X() << ", " << intersection_point.Y() << ")" << std::endl;
                
                    
                    // CHECK 2: junction between the main direction and the triangle vertex
                    // is contained within the triangle
                    bool junctionContained = TrackTriangleJunctionConatined(MainDirection, triangle);


                    // CHECK 3: check the juntion doesn't cross other tracks
                    //int nHitsInMiddle = GetNHitsBetweenJunction(MainDirection, triangle, FreeTracksList, intP, track1ListIds, track2ListIds, 1.5);
                    int nHitsInMiddle = GetNHitsBetweenJunction(MainDirection, triangle, FreeTracksList, intP, {track1.GetId()}, {track2.GetId()}, 1.5);
                    bool passJunctionIsFree = (nHitsInMiddle <= 1);
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "NHits in middle" << nHitsInMiddle << "Pass?" <<passJunctionIsFree  << std::endl;


                    // CHECK 4: how many hits of triangle tracks are contained in the main direciton hypothesis
                    bool passAngleTracksNotInMain = true;
                    int nhits_track1_inMainDir = GetHitsContainedInLineEquation(MainDirection.GetTrackEquation(), track1.GetHits());
                    int nhits_track2_inMainDir = GetHitsContainedInLineEquation(MainDirection.GetTrackEquation(), track2.GetHits());
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" NHits of track "<<track1.GetId()<<" in main direction: "<<nhits_track1_inMainDir<<std::endl;
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" NHits of track "<<track2.GetId()<<" in main direction: "<<nhits_track2_inMainDir<<std::endl;
                    passAngleTracksNotInMain = (1.*nhits_track1_inMainDir/track1.NHits()<fTPCLinesVertexFinderPset.MaxTrackFractionInMain && 
                    1.*nhits_track2_inMainDir/track2.NHits()<fTPCLinesVertexFinderPset.MaxTrackFractionInMain);
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Pass AngleTracksNotInMain: "<<passAngleTracksNotInMain<<std::endl;
                    

                    // CHECK 5: triangle tracks start/end are not next to the main track edge hit
                    bool passTriangleEdgesNotInMain = true; 
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Triangle vertex BC: "<<triangle.GetVertexB().X()<<" "<<triangle.GetVertexB().Y()<<
                    " "<<triangle.GetVertexC().X()<<" "<<triangle.GetVertexC().Y()<<std::endl;
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" MainDirectionEdgeHits: "<<MainDirection.GetStartPoint().X()<<" "<<MainDirection.GetStartPoint().Y()<<" "<<MainDirection.GetEndPoint().X()<<" "<<MainDirection.GetEndPoint().Y()<<std::endl;
                    
                    double distanceBMainStart = TPCLinesDistanceUtils::GetHitDistance(triangle.GetVertexB(), MainDirection.GetStartPoint());
                    double distanceBMainEnd = TPCLinesDistanceUtils::GetHitDistance(triangle.GetVertexB(), MainDirection.GetEndPoint());
                    double distanceCMainStart = TPCLinesDistanceUtils::GetHitDistance(triangle.GetVertexC(), MainDirection.GetStartPoint());
                    double distanceCMainEnd = TPCLinesDistanceUtils::GetHitDistance(triangle.GetVertexC(), MainDirection.GetEndPoint());
                    // get the average compactness of all the tracks in the main direction
                    double meanComp = 0;
                    for(SLinearCluster & trk:MainDirTrackList){
                        meanComp+=trk.GetCompactness();
                    }
                    meanComp/=MainDirTrackList.size();
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" MainDir compactness:"<<meanComp<<std::endl;
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << " Distances B: "<<distanceBMainStart<<" "<<distanceBMainEnd<<std::endl;
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << " Distances C: "<<distanceCMainStart<<" "<<distanceCMainEnd<<std::endl;
                    
                    double _tol = 1.;

                    passTriangleEdgesNotInMain = distanceBMainStart>_tol*meanComp &&
                                                    distanceBMainEnd>_tol*meanComp &&
                                                    distanceCMainStart>_tol*meanComp &&
                                                    distanceCMainEnd>_tol*meanComp;

                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Pass edge triangle hits not in main track edges: "<<passTriangleEdgesNotInMain<<std::endl;
                    
                    if(triangleIntersects && junctionContained && passJunctionIsFree && passAngleTracksNotInMain && passTriangleEdgesNotInMain){
                        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"FOUND ORIGIN!\n";
                        vertexList.push_back(triangle);
                        originList.push_back(intP);
                        mainDirection = MainDirection;
                    }

                }

            }
        }
    }


    // Look for collinear clusters
    if(fTPCLinesVertexFinderPset.AddCollinearLines){
        std::vector<SLinearCluster> collinearTriangles = GetCollinearTracks(MainDirection, FreeTracksList);
        if(collinearTriangles.size()>0){
            for(SLinearCluster& colTrack:collinearTriangles){
                SHit intHit(colTrack.GetStartPoint().X(), colTrack.GetStartPoint().Y());
                STriangle triangle = STriangle(colTrack.GetStartPoint(), colTrack.GetEndPoint(), colTrack.GetEndPoint(), intHit);
                vertexList.push_back(triangle);
                originList.push_back(triangle.GetMainVertex());
            }
        }
    }
    

    return;
}


std::vector<SOrigin> TPCLinesVertexFinder::GetInterectionsInBall(std::vector<SLinearCluster> tracksList, SPoint ballVertex){

    std::vector<STriangle> triangleList;
    std::vector<SOrigin> originList;
    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" In vertex finder\n";

    std::vector<bool> usedTrack(tracksList.size(), false);
    std::vector<std::pair<int, int>> intersectingTracks;

    // ------- Look for possible intersections
    // ------ Loop 1
    if (tracksList.size() > 0) {
        for (size_t ix = 0; ix < tracksList.size(); ++ix) {
            SLinearCluster track1 = tracksList[ix];
            std::cout<<ix<<std::endl;
            // Loop through other tracks
            // ------ Loop 2
            for (size_t jx = ix+1; jx < tracksList.size(); ++jx) {
                std::cout<<jx<<std::endl;
                SLinearCluster track2 = tracksList[jx];


                // ------ Check if the tracks are connected
                // calculate connection based on the tracks connectedes
                float connTol = 20 * std::max(track1.GetConnectedness(), track2.GetConnectedness());
                float conn = TPCLinesDistanceUtils::GetClusterConnectedness(track1.GetHitCluster(), track2.GetHitCluster());
                bool connected = (conn < connTol); 

                if (!connected) {
                    continue;
                }

                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "\n----- Potential intersection " << track1.GetId() << " " << track2.GetId() << " Conn:" << conn << " Tol:" << connTol << std::endl;


                // ------ Look for the intersection points
                SPoint intP = GetTracksIntersection(track1, track2, 10, fTPCLinesVertexFinderPset.RefineVertexIntersection);
                if(intP.X()==-1 and intP.Y()==-1) continue;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   ** Intersection in: " << intP ;
                

                // Get the closest hit of each track to the intersection
                std::pair<SHit, double> cloHit1Pair = track1.GetHitCluster().GetClosestHitToPoint(intP);
                SHit cloHit1 = cloHit1Pair.first;
                double dHit1 = cloHit1Pair.second;
                std::pair<SHit, double> cloHit2Pair = track2.GetHitCluster().GetClosestHitToPoint(intP);
                SHit cloHit2 = cloHit2Pair.first;
                double dHit2 = cloHit2Pair.second;
                // Choose the closest hit of the two candidates
                SHit cloHit = (dHit1 < dHit2) ? cloHit1 : cloHit2;
                //double dHit = (dHit1 < dHit2) ? dHit1 : dHit2;

                // Study the edges of the tracks
                std::vector<SHit> vertexHits = GetMutualHitsContainedInHypo(track1, track2, intP, 5);
                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << "    NVERTEX HITS " << vertexHits.size() << std::endl;
                    for(SHit & h:vertexHits) std::cout<<h;
                }

                if(vertexHits.size()==0){
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   No vertex hits...skipping\n";
                    continue;
                }

                float dEdge1 = std::min(std::abs(cloHit1.X() - track1.GetMinX()), std::abs(cloHit1.X() - track1.GetMaxX()));
                float dEdge2 = std::min(std::abs(cloHit2.X() - track2.GetMinX()), std::abs(cloHit2.X() - track2.GetMaxX()));

                for (SHit& hit : vertexHits) {
                    float min1 = std::min(std::abs(hit.X() - track1.GetMinX()), std::abs(hit.X() - track1.GetMaxX()));
                    float min2 = std::min(std::abs(hit.X() - track2.GetMinX()), std::abs(hit.X() - track2.GetMaxX()));
                    
                    if (min1 < dEdge1) {
                        dEdge1 = min1;
                    }
                    if (min2 < dEdge2) {
                        dEdge2 = min2;
                    }
                }

                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   Clo Hit1/2:\n  " << cloHit1 << "  " << cloHit2 << " DEdges: " << dEdge1 << " " << dEdge2 << std::endl;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Clo Hit: " << cloHit;

                
                if (fTPCLinesVertexFinderPset.UseEdgesDiscard && (dEdge2 >= fTPCLinesVertexFinderPset.MaxDistToEdge || dEdge1 >= fTPCLinesVertexFinderPset.MaxDistToEdge)) {
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "Skipping...intersection is not with edge hits" << std::endl;
                    continue;
                }
                else{
                    std::cout<<"   Pass edges test!\n";
                }

                // check connectedness overlap
                double dOverlap = TPCLinesDistanceUtils::GetClusterConnectednessOverlap(track1.GetHitCluster(), track2.GetHitCluster());
                std::cout<<"  Connectedness Overlap: "<<dOverlap<<std::endl;
                if(dOverlap>3){
                    std::cout<<"    Do not overlap\n";
                    continue;
                }

                float minD = 1e3;
                SHit vertexHit(intP.X(), intP.Y());
                for (SHit& hit : vertexHits) {
                    float d = std::hypot(hit.X() - intP.X(), hit.Y() - intP.Y());
                    if (d < minD) {
                        minD = d;
                        vertexHit = hit;
                    }
                }

                //  put intersection point in the closest hit
                intP = SPoint(vertexHit.X(), vertexHit.Y());

                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << "      Vertex set to: " << intP;
                }
                
                // check distance to the PANDORA vertex
                float dIntPBallVertex = std::hypot( 0.3*(intP.X() - ballVertex.X()), 0.075*(intP.Y() - ballVertex.Y()) );
                if( dIntPBallVertex < 30 ){
                    // make the triangle
                    SPoint vertex1 = GetTrackssEuationOppositePoint( track1, {track1}, intP );
                    SPoint vertex2 = GetTrackssEuationOppositePoint( track2, {track2}, intP );
                    STriangle triangle = STriangle(intP, vertex1, vertex2, cloHit, track1.GetIntegral(), track2.GetIntegral());
                    triangleList.push_back( triangle );
                    
                    if(usedTrack[ix]==false && usedTrack[jx]==false){
                        intersectingTracks.push_back(std::make_pair(ix, jx));
                        SOrigin newOr = SOrigin(intP, {track1, track2}, true);
                        originList.push_back( newOr );
                        // mark as used tracks
                        usedTrack[ix]=true;
                        usedTrack[jx]=true;
                        std::cout<<"   END = Adding completely new origin..."<<newOr;
                    }
                    else if(usedTrack[ix]==false || usedTrack[jx]==false){
                        
                        SLinearCluster newTrack = (usedTrack[ix]==false) ? track1 : track2;
                        SLinearCluster previousTrack = (usedTrack[ix]==true) ? track1 : track2;

                        bool merged=false;
                        for(SOrigin & ori:originList){
                            if(ori.HasTrackIndex(previousTrack.GetId())==true){
                                float d = std::hypot(ori.GetPoint().X() - intP.X(), ori.GetPoint().Y() - intP.Y());
                                if(d<10){
                                    ori.AddTrack(newTrack, intP);
                                    std::cout<<"   END = Adding new track to origin "<<ori.GetPoint().X()<<ix<<" "<<jx<<std::endl;
                                    usedTrack[ix]=true;
                                    usedTrack[jx]=true;
                                    merged=true;
                                }
                            }
                        }

                        if(merged==false){
                            SOrigin newOr(intP, {track1, track2}, true);
                            std::cout<<"   END = Adding new origin..."<<newOr;
                            originList.push_back( newOr );
                            usedTrack[ix]=true;
                            usedTrack[jx]=true;
                        }

                    }
                    
                }
                
            }
        }
    }

    // add origins of unmatched tracks with edge near PANDORA vertex
    for(size_t k=0; k<tracksList.size(); k++){
    
        if(usedTrack[k]==true) continue;
        SLinearCluster track = tracksList[k];
        float d1 = std::hypot( 0.3*(track.GetStartPoint().X() - ballVertex.X()), 0.075*(track.GetStartPoint().Y() - ballVertex.Y()) );
        float d2 = std::hypot( 0.3*(track.GetEndPoint().X() - ballVertex.X()), 0.075*(track.GetEndPoint().Y() - ballVertex.Y()) );
        float d = std::min(d1, d2);
        if( d < 20 && track.NHits()>=5 ){
            SPoint edgePoint = (d1<d2)? track.GetStartPoint() : track.GetEndPoint();
            originList.push_back( SOrigin(edgePoint, {track}, true) );  
        }

    }


    // sort by NHits
    std::sort(originList.begin(), originList.end(), [](SOrigin& obj1, SOrigin& obj2) {
        return obj1.NHits() > obj2.NHits();
    });


    return originList;
}