////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesAlgo.h"

//----------------------------------------------------------------------
// Constructor
TPCLinesAlgo::TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset):
    fTPCLinesPset(tpcLinesAlgoPset),
    fNTotalHits(0),
    fHitList({}),
    fVertex(SPoint(-1, -1), ""),
    fHoughAlgo(tpcLinesAlgoPset.HoughAlgorithmPset),
    fTrackFinder(tpcLinesAlgoPset.TrackFinderAlgorithmPset),
    fVertexFinder(tpcLinesAlgoPset.VertexFinderAlgorithmPset),
    fDisplay(TPCLinesDisplay(tpcLinesAlgoPset.Verbose>0, tpcLinesAlgoPset.OutputPath))
{
}


//----------------------------------------------------------------------
// sy function
void TPCLinesAlgo::Display(std::string name){

    fDisplay.Show(name, fHitList, LineEquation(0, 0), {}, fFinalTrackList, fMainDirection, fAngleList, fVertex, fVertexTrue, fOrigins);

    return;
}


//----------------------------------------------------------------------
// Get the view with the less Chi2
int TPCLinesAlgo::GetBestView(std::vector<int> *_View, std::vector<double> *_Chi2){

        std::map<int, double> meanChi2Map;

        std::map<int, double> sumChi2Map;
        std::map<int, int> countMap;

        for (size_t i = 0; i < _View->size(); ++i) {
            int viewValue = (*_View)[i];
            double chi2Value = (*_Chi2)[i];

            sumChi2Map[viewValue] += chi2Value;
            countMap[viewValue]++;
        }

        int minChi2View = -1;
        double minChi2=1e4;
        for (const auto& entry : sumChi2Map) {
            int viewValue = entry.first;
            double sumChi2 = entry.second;
            int count = countMap[viewValue];

            
            if (count != 0) {
                double meanChi2 = sumChi2/count; 
                if(meanChi2<minChi2){
                    minChi2=meanChi2;
                    minChi2View=viewValue;
                }
            }

        }

        return minChi2View;
}


//----------------------------------------------------------------------
// Set the input hits
bool TPCLinesAlgo::SetHitList(int view,
                            std::vector<int>& vertex,
                            std::vector<int>& vertexTrue,
                            std::vector<SHit> hits)
                            /*std::vector<int> *_X,
                            std::vector<double> *_Y,
                            std::vector<double> *_Int,
                            std::vector<double> *_Wi,
                            std::vector<double> *_ST,
                            std::vector<double> *_ET,
                            std::vector<int> *_View,
                            std::vector<double> *_Chi2)*/
{   
    // reset the class variables
    fHitList.clear();
    fNTotalHits=0;

    // Define the reco vertex
    double vertexX = vertex[2];
    if (view == 0) vertexX = vertex[0];
    if (view == 1) vertexX = vertex[1];
    double vertexY = vertex[3] / fTPCLinesPset.DriftConversion;
    
    // Define the true vertex
    double vertexXTrue = vertexTrue[2];
    if (view == 0) vertexXTrue = vertexTrue[0];
    if (view == 1) vertexXTrue = vertexTrue[1];
    double vertexYTrue = vertexTrue[3] / fTPCLinesPset.DriftConversion;
    if(fTPCLinesPset.Verbose>=2) std::cout << "  ** Vertex XY: " << vertexX << " " << vertexY << std::endl;
    
    if (vertexX != -1){

        // Get hits in the ROI circle around the vertex
        // Hits are already filtered by plane type
        for (size_t i = 0; i < hits.size(); i++) { 
            
            double x = hits[i].X();
            double y = hits[i].Y();
            
            double d = std::sqrt(std::pow(x - vertexX, 2) + std::pow(y - vertexY, 2));
            
            if (d < fTPCLinesPset.MaxRadius) {
                fHitList.push_back(hits[i]);
                /*if(hits[i].Width()>3.5){
                    SHit newHit1(hits.size()+i, hits[i].X(), hits[i].Y()+hits[i].Width(), hits[i].Width()/2, hits[i].StartT(), hits[i].EndT(), hits[i].Integral(), hits[i].Chi2());
                    SHit newHit2(2*hits.size()+i, hits[i].X(), hits[i].Y()-hits[i].Width(), hits[i].Width()/2, hits[i].StartT(), hits[i].EndT(), hits[i].Integral(), hits[i].Chi2());
                    fHitList.push_back(newHit1);
                    fHitList.push_back(newHit2);
                }*/
            }
        }

        if ( (int)fHitList.size() > fTPCLinesPset.MinTrackHits) {
            
            // Shift hits to have origin in (0,0), leave 3 ticks under/overflow
            double minX = std::min_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.X()<h2.X();})->X() - 3;
            double minY = std::min_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.Y()<h2.Y();})->Y() - 3;
            double maxX = std::max_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.X()>h2.X();})->X() - 3;        

            for (size_t i = 0; i < fHitList.size(); i++) {
                fHitList[i] = SHit(
                                i,
                                fHitList[i].X() - minX,
                                fHitList[i].Y() - minY,
                                fHitList[i].Width(),
                                fHitList[i].Integral(),
                                fHitList[i].StartT() - minY,
                                fHitList[i].EndT() - minY,
                                fHitList[i].Chi2()
                                );

            }
            fNTotalHits = fHitList.size();

            // create vertex with common origin
            vertexX = vertexX - minX;
            vertexY = vertexY - minY;
            vertexXTrue = vertexXTrue - minX;
            vertexYTrue = vertexYTrue - minY;
            fVertex = SVertex(SPoint(vertexX, vertexY), std::to_string(view));
            fVertexTrue = SVertex(SPoint(vertexXTrue, vertexYTrue), std::to_string(view));

            // set edge variables
            fMaxX = maxX - minX;
            fMinX = minX;
            fMinY = minY;
            std::cout << "  ** Origin vertex: " << fVertex;
            std::cout << "  ** NInputHits:"<<fNTotalHits<<std::endl;

            return true;
        }
        else {
            std::cout << "   SKIPPED NHits selected near the vertex " << fHitList.size() << std::endl;
            return false;
        }

    }

    return false;
}


//----------------------------------------------------------------------
// Function to get the average hit density
double TPCLinesAlgo::GetAverageHitDensity(){
    // initialize hit map
    std::map<int, int> fHitMap;
    for(int i=0; i<=fMaxX; i++){
        fHitMap[i]=0;
    }

    for(SHit & h:fHitList){
        fHitMap[h.X()]+=1;
    }

    double sum=0;
    int nwires=0;
    for(int i=0; i<=fMaxX; i++){
        if(fHitMap[i]>0 and 0.3*std::abs(i-fVertex.X())<40 ){
            sum+=fHitMap[i];
            nwires++;
        }
        
    }

    return sum/nwires;
}

//----------------------------------------------------------------------
// Remove isolated hits
std::vector<SHit> TPCLinesAlgo::RemoveIsolatedHits(std::vector<SHit> hitListForHough, std::vector<SHit>& discardedHits, double maxD, int minNeighbours) {

    if(fTPCLinesPset.Verbose>=2) std::cout << "\n\n +-+-+-+-+-+-+- Removing isolated hits +-+-+-+-+-+-+-" << std::endl;

    std::vector<SHit> hitListForHoughWithoutIsolatedHits;

    for (SHit& hit : hitListForHough) {
        int NNeighboursHits = 0;
        for (SHit& hit2 : hitListForHough) {
            if (hit.Id() == hit2.Id()) {
                continue;
            }
            //double d = GetHitDistance1D(hit, hit2);
            double d = std::abs(hit.X()-hit2.X());
            if (d <= maxD) {
                NNeighboursHits++;
            }

        }
        
        if (NNeighboursHits >= minNeighbours) {
            hitListForHoughWithoutIsolatedHits.push_back(hit);
        }
        else{
            discardedHits.push_back(hit);
        }
    }

    return hitListForHoughWithoutIsolatedHits;
}


//----------------------------------------------------------------------
// Merge isolated hits
std::vector<SLinearCluster> TPCLinesAlgo::MergeIsolatedHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, double dCleaning1D, double dTh){
    if(fTPCLinesPset.Verbose>=2) std::cout << "\n\n +-+-+-+-+-+-+- Isolated hit merger +-+-+-+-+-+-+-" << std::endl;
    std::vector<SLinearCluster> newRecoTrackList;

    // initializations
    std::vector<std::vector<int>> hitToTrackDict(hitList.size());
    std::vector<bool> mergedHitsCounter(hitList.size(), false);
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        hitToTrackDict[hitix] = {};
    }

    // first collect potential hits to merge
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        SHit hit = hitList[hitix];
        double minD = 1e3;
        for (size_t trkix = 0; trkix < recoTrackList.size(); ++trkix) {
            SLinearCluster track = recoTrackList[trkix];
            double trackComp = track.GetCompactness();

            double hitTrackDist = track.GetHitCluster().GetMinDistanceToCluster1D(hit);
            if (hitTrackDist < dCleaning1D) {
                double d = track.GetTrackEquation().GetDistance(SPoint(hit.X(), hit.Y()));
                //double hypoY = track.GetTrackEquation().Slope() * hit.X() + track.GetTrackEquation().Intercept();
                if (d < dTh * trackComp && d < minD) {
                    minD = d;
                    hitToTrackDict[hitix].push_back(trkix);
                }
            }
        }
    }

    for (size_t trkix = 0; trkix < recoTrackList.size(); ++trkix) {
        SLinearCluster track = recoTrackList[trkix];
        double trackComp = track.GetCompactness();
        double trackConn = track.GetConnectedness();

        // Get average residual
        double trackAvResidual=0.;
        LineEquation trk_eq = track.GetTrackEquation();
        for (SHit hit : track.GetHits()) {
            trackAvResidual = trackAvResidual +  trk_eq.GetDistance(SPoint(hit.X(), hit.Y()));
        }
        trackAvResidual/=track.NHits();
        if(fTPCLinesPset.Verbose>=2) std::cout << "\n Merging track " << trkix << " " <<track.GetMinX()<<" "<<track.GetMaxX()<< " COMP " << trackComp << " CONN " << trackConn << " AvResidual: "<<trackAvResidual<<std::endl;

        std::vector<SHit> candidateHits;
        std::vector<int> candidateHitsIx;
        for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
            if (!mergedHitsCounter[hitix]) {
                if (std::find(hitToTrackDict[hitix].begin(), hitToTrackDict[hitix].end(), trkix) != hitToTrackDict[hitix].end()) {
                    candidateHits.push_back(hitList[hitix]);
                    candidateHitsIx.push_back(hitix);
                }
            }
        }

        std::vector<double> candidateHitsDist;
        for (SHit hit : candidateHits) {
            double distToCluster = track.GetHitCluster().GetMinDistanceToCluster(hit);
            candidateHitsDist.push_back(distToCluster);
        }

        // sort hits by distance to the cluster
        std::vector<int> sortedCandidateHits(candidateHits.size());
        std::iota(sortedCandidateHits.begin(), sortedCandidateHits.end(), 0);
        std::sort(sortedCandidateHits.begin(), sortedCandidateHits.end(), [&](int i, int j)->bool {
            return candidateHitsDist[i] < candidateHitsDist[j];
        });

        std::vector<int> sortedCandidateHitsIx(candidateHitsIx.size());
        std::iota(sortedCandidateHitsIx.begin(), sortedCandidateHitsIx.end(), 0);
        std::sort(sortedCandidateHitsIx.begin(), sortedCandidateHitsIx.end(), [&
        ](int i, int j) {
            return candidateHitsDist[i] < candidateHitsDist[j];
        });

        std::vector<SHit> newHitList = track.GetHits();
        SCluster currentCluster = track.GetHitCluster();

        if(fTPCLinesPset.Verbose>=2) std::cout << "  candidates hits: " << std::endl;
        for (size_t ix = 0; ix < candidateHits.size(); ++ix) {
            SHit hit = candidateHits[sortedCandidateHits[ix]];

            LineEquation trk_eqClosest = (std::abs(hit.X()-track.GetMinX())<std::abs(hit.X()-track.GetMaxX()))? track.GetTrackEquationStart():track.GetTrackEquationEnd();

            double minD = currentCluster.GetMinDistanceToCluster(hit);
            double minDConn = currentCluster.GetMinDistanceToClusterOverlap(hit);

            double hypoY = track.GetTrackEquation().Slope() * hit.X() + track.GetTrackEquation().Intercept();
            double minDHypo = std::abs(hypoY-hit.Y());
            
            double residual =  trk_eqClosest.GetDistance(SPoint(hit.X(), hit.Y()));

            if(fTPCLinesPset.Verbose>=2) std::cout << "       " << hit.Id() << "  X=" << hit.X() << " Y=" << hit.Y() << " " << " d " << minD << " dConn " << minDConn << " dCompHypo" << minDHypo << " residual: "<<residual << std::endl;
            if (minDConn < dTh * trackConn && minDHypo < dTh * trackComp && minD < dTh * trackComp){// && residual < dTh*trackAvResidual) {
                newHitList.push_back(hit);
                mergedHitsCounter[candidateHitsIx[sortedCandidateHits[ix]]] = true;
                currentCluster = SCluster(newHitList);
            }
        }

        newRecoTrackList.push_back(SLinearCluster(newHitList));
    
    }

    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Main function
void TPCLinesAlgo::AnaView(std::string eventLabel)
{
    // --- Objects for the Hough tracks
    int trkIter = 0;
    std::vector<SHit> hitListForHough = fHitList;
    std::vector<SHit> discardedHits;
    std::vector<SLinearCluster> finalLinearClusterV;

    trkIter = (int)hitListForHough.size()>fTPCLinesPset.MinTrackHits ? 0 : fTPCLinesPset.MaxHoughTracks;
    std::cout<<trkIter<<" "<<fHitList.size()<<" "<<hitListForHough.size()<<std::endl;
    // --- Loop over the hough tracks
    while(trkIter<fTPCLinesPset.MaxHoughTracks){

        if(fTPCLinesPset.Verbose>=2) std::cout<<" **** Track finder iteration: "<<trkIter<< " # hough hits:"<<hitListForHough.size()<<std::endl;

        // -- Get the best Hough line       
        HoughLine houghLine = fHoughAlgo.GetBestHoughLine(hitListForHough, fVertex);
        if(fTPCLinesPset.Verbose>=3)
            if(fTPCLinesPset.Verbose>=2) std::cout<<"    Hough line results: "<<houghLine.NHits()<<" Score: "<<houghLine.Score()<<std::endl;

        // -- Skip if not enough hits
        if(houghLine.NHits() < fTPCLinesPset.MinTrackHits){
            if(fTPCLinesPset.Verbose>=3)
                if(fTPCLinesPset.Verbose>=2) std::cout<<"   Hough lines has <"<<fTPCLinesPset.MinTrackHits<<", ending the search\n";
            trkIter = fTPCLinesPset.MaxHoughTracks;
        }

        // -- Call the track finfder otherwise
        else{
            std::vector<SLinearCluster> linearClusterV = fTrackFinder.ReconstructTracksFromHoughDirection(hitListForHough, houghLine.GetLineEquation(), trkIter);

            // -- check the found tracks has enough hits
            if(linearClusterV.size()!=0){
                
                finalLinearClusterV.insert(finalLinearClusterV.begin(), linearClusterV.begin(), linearClusterV.end());
                
                if(fTPCLinesPset.Verbose>=2){
                    std::cout<<"   Found "<<linearClusterV.size()<<" new tracks! Current track list: \n";
                    for(SLinearCluster &trk:finalLinearClusterV){
                        std::cout<<"   "<<trk.GetId()<<" "<<trk.GetMinX()<<" "<<trk.GetMaxX()<<std::endl;
                    }
                }

                // -- check there's enough hits for the next Hough iteration
                if((int)hitListForHough.size() < fTPCLinesPset.MinTrackHits){
                    if(fTPCLinesPset.Verbose>=3)
                        if(fTPCLinesPset.Verbose>=2) std::cout<<"   Remaining hits are "<<hitListForHough.size()<<", ending the search\n";
                    trkIter = fTPCLinesPset.MaxHoughTracks;
                }
            }
            // end Hough loop otherwise
            else{
                if(fTPCLinesPset.Verbose>=3)
                    if(fTPCLinesPset.Verbose>=2) std::cout<<"   Not found new tracks\n";
                trkIter = fTPCLinesPset.MaxHoughTracks;
            }
            
            // -- Remove isolated hits
            if (fTPCLinesPset.RemoveIsolatedHits) {
                hitListForHough = RemoveIsolatedHits(hitListForHough, discardedHits, fTPCLinesPset.MaxNeighbourDistance, fTPCLinesPset.MinNeighboursHits);
            }
            
        }

        trkIter++;
    } // -- end Hough loop


    for(size_t ix = 0; ix<finalLinearClusterV.size(); ix++){
        std::cout<<" ix:"<<ix<<" Min/MaxX:"<<finalLinearClusterV[ix].GetMinX()<<" "<<finalLinearClusterV[ix].GetMaxX()<<std::endl;
    }
    
    
    // Slope track merger
    // sort by minX
    std::sort(finalLinearClusterV.begin(), finalLinearClusterV.end(), [&](SLinearCluster& l1, SLinearCluster& l2) {return l1.GetMinX() < l2.GetMinX();} );    
    finalLinearClusterV = TPCLinesDirectionUtils::SlopeTrackMerger(finalLinearClusterV, 2, 15, fTPCLinesPset.Verbose); 



    // Isolated hit merger
    std::vector<SHit> remainingHits = hitListForHough;
    remainingHits.insert(remainingHits.end(), discardedHits.begin(), discardedHits.end());
    //fsy.Show(eventLabel+"_BeforeIsolatedMerger_TPCLines_", fHitList, LineEquation(0, 0), {}, finalLinearClusterV);
    finalLinearClusterV = MergeIsolatedHits(finalLinearClusterV, remainingHits, 10);
    //fsy.Show(eventLabel+"_AfterIsolatedMerger_TPCLines_", fHitList, LineEquation(0, 0), {}, finalLinearClusterV);



    
    // Characterize the tracks
    for(size_t ix = 0; ix<finalLinearClusterV.size(); ix++){
        finalLinearClusterV[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
        finalLinearClusterV[ix].AssignId(ix);
    }

    // ------ Get the short/vertex tracks and remove them
    std::map<int, std::vector<int>> shortToLongTrackDict;
    std::map<int, int> shortConnectionTrackDict;
    std::vector<SLinearCluster> shortTrackList = TPCLinesDirectionUtils::GetVertexTracks(finalLinearClusterV, shortToLongTrackDict, shortConnectionTrackDict, 6, 3, 5, fTPCLinesPset.Verbose); 
    // remove the short tracks
    std::vector<SLinearCluster> newTrackList;
    for (SLinearCluster& track : finalLinearClusterV) {
    
        std::cout<<"TRACK ID: "<<track.GetId()<<std::endl;
        if (shortToLongTrackDict.find(track.GetId()) == shortToLongTrackDict.end()) {
            std::cout<<"   Merging Short for : "<<track.GetId()<<std::endl;

            SLinearCluster newTrack = track;
            
            for(auto & sTrackId:shortToLongTrackDict){
                if( sTrackId.second[0]==track.GetId() || sTrackId.second[1]==track.GetId() ) {
                    std::vector<SHit> auxHits = track.GetHits();
                    std::vector<SHit> auxHitsShort;
                    for(SLinearCluster sTrack:shortTrackList){
                        
                        if(sTrack.GetId()==sTrackId.first){
                            auxHitsShort = sTrack.GetHits(); 
                        }
                    }
                    
                    auxHits.insert(auxHits.begin(), auxHitsShort.begin(), auxHitsShort.end());
                    newTrack = SLinearCluster(auxHits);
                }
                
            }
            
            newTrack.FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            newTrack.AssignId(track.GetId());
            newTrackList.push_back(newTrack);
        }
    }

    finalLinearClusterV.clear();
    for(SLinearCluster & trk: newTrackList){
        finalLinearClusterV.push_back(trk);
    }

    std::vector<SLinearCluster> NewTrackList;
    NewTrackList.clear();
    SLinearCluster mainDirection;
    std::vector<SOrigin> intersectionsInBall;
    intersectionsInBall.clear();
    std::vector<STriangle> vertexList;
    vertexList.clear();
    std::vector<SOrigin> associatedOrigins;
    associatedOrigins.clear();

    if(finalLinearClusterV.size()>0){
        
        // ------- Get the parallel tracks
        std::vector<std::vector<SLinearCluster>> parallelTracks = TPCLinesDirectionUtils::GetParallelTracks(finalLinearClusterV, -2, 15, 15, fTPCLinesPset.Verbose);
        NewTrackList.clear();

        std::cout<<" Final tracks for origin finder: \n";
        for(size_t ix = 0; ix<parallelTracks.size(); ix++){ 
            std::cout<<"\n "<<parallelTracks[ix][0].GetId();
            std::vector<SHit> hitList = parallelTracks[ix][0].GetHits();
            for(size_t jx = 1; jx<parallelTracks[ix].size(); jx++){
                std::cout<<" "<<parallelTracks[ix][jx].GetId();
                std::vector<SHit> moreHits = parallelTracks[ix][jx].GetHits();
                hitList.insert(hitList.end(), moreHits.begin(), moreHits.end());
            }
            SLinearCluster newTrack(hitList);
            newTrack.FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            newTrack.AssignId(ix);
            NewTrackList.push_back( newTrack );    
        }

        //Find secondary vertexes
        std::cout<<" We have "<<NewTrackList.size()<<" tracks\n";
        if(NewTrackList.size()>0){
            intersectionsInBall = fVertexFinder.GetAngleVertices(NewTrackList, fVertex.Point(), vertexList, associatedOrigins, mainDirection);
        }

        // Set the main vertex
        fMainVertex = mainDirection.GetStartPoint();
        float d1 = std::hypot( 0.3*(mainDirection.GetStartPoint().X() - fVertex.Point().X()), 0.075*(mainDirection.GetStartPoint().Y() - fVertex.Point().Y()) );
        float d2 = std::hypot( 0.3*(mainDirection.GetEndPoint().X() - fVertex.Point().X()), 0.075*(mainDirection.GetEndPoint().Y() - fVertex.Point().Y()) );
        fMainVertex = (d1<d2)? mainDirection.GetStartPoint() : mainDirection.GetEndPoint();

        // Save final results
        fFinalTrackList = NewTrackList;
        fMainDirection = mainDirection;
        fAngleList = vertexList;
        fOrigins = intersectionsInBall;

        for(SOrigin & ori: intersectionsInBall){
            std::cout<<ori;
        }

    }



    
    
    double hitDensity = GetAverageHitDensity();
    std::cout<<"Hit density: "<<hitDensity<<std::endl;

    SEvent recoEvent(NewTrackList, intersectionsInBall, vertexList, associatedOrigins, hitDensity);
    fRecoEvent = recoEvent;

    return;
}