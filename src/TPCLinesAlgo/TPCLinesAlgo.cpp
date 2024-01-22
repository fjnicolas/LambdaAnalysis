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
TPCLinesAlgo::TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, FRANSPsetType fransPset):
    fTPCLinesPset(tpcLinesAlgoPset),
    fNTotalHits(0),
    fHitList({}),
    fHitListOutROI({}),
    fUnmatchedHits({}),
    fVertex(SPoint(-1, -1), ""),
    fHoughAlgo(tpcLinesAlgoPset.HoughAlgorithmPset),
    fTrackFinder(tpcLinesAlgoPset.TrackFinderAlgorithmPset),
    fVertexFinder(tpcLinesAlgoPset.VertexFinderAlgorithmPset),
    fFRANSAlgo(fransPset, tpcLinesAlgoPset.View),
    fTriangleCalo(tpcLinesAlgoPset.CaloAlgorithmPset),
    fHasTriangle(false),
    fDisplay(TPCLinesDisplay(tpcLinesAlgoPset.Verbose>0, tpcLinesAlgoPset.OutputPath))
{
}


//----------------------------------------------------------------------
// Display function
void TPCLinesAlgo::Display(std::string name, TCanvas *canvas, TCanvas *canvasCalo, TCanvas *canvasFRANS){

    std::vector<SHit> hitList = fHitList;
    hitList.insert(hitList.end(), fHitListOutROI.begin(), fHitListOutROI.end());

    fDisplay.Show(false, name, hitList, LineEquation(0, 0), fUnmatchedHits, fFinalTrackList, fMainDirection, fAngleList, fVertex, fVertexTrue, fOrigins, canvas);

    if(fHasTriangle){
        fTriangleCalo.Display(canvasCalo);
        fFRANSAlgo.Display(canvasFRANS);
    }

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
std::vector<SLinearCluster> TPCLinesAlgo::MergeIsolatedHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, double dCleaning1D, std::vector<SHit> & discardedHits, double dTh)
{
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
            if (minDConn < dTh * trackConn && minDHypo < dTh * trackComp && minD < dTh * trackComp){
                newHitList.push_back(hit);
                mergedHitsCounter[candidateHitsIx[sortedCandidateHits[ix]]] = true;
                currentCluster = SCluster(newHitList);
            }
        }

        newRecoTrackList.push_back(SLinearCluster(newHitList));
    
    }

    // add non merged hits to the discarded vector
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        if (!mergedHitsCounter[hitix]) {
            discardedHits.push_back(hitList[hitix]);
        }
    }

    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Merge out of ROI hits
std::vector<SLinearCluster> TPCLinesAlgo::MergeOutOfROIHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList)
{

    // map of vector of hits 
    std::map<int, std::vector<SHit>> newHitsMap;
    for(size_t k=0; k<recoTrackList.size(); k++){
        newHitsMap[k] = {};
    }

    for(SHit & h:fHitListOutROI){
        for(size_t k=0; k<recoTrackList.size(); k++){
            std::cout<<" MAIN CLUSTER ID "<<recoTrackList[k].GetHitCluster().GetMainClusterID()<<std::endl;
            if( h.ClusterId()==recoTrackList[k].GetHitCluster().GetMainClusterID() ){
                newHitsMap[k].push_back(h);
            }
        }
    }


    std::vector<SLinearCluster> newRecoTrackList;
    for(size_t k=0; k<recoTrackList.size(); k++){
        std::vector<SHit> newHitList = recoTrackList[k].GetHits();
        newHitList.insert(newHitList.end(), newHitsMap[k].begin(), newHitsMap[k].end());
        SLinearCluster newCluster(newHitList);
        newRecoTrackList.push_back(newCluster);
    }



    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Create out of ROI clusters
std::vector<SLinearCluster> TPCLinesAlgo::CreateOutOfROIClusters(std::vector<SLinearCluster> recoTrackList)
{

    std::cout<<"UNUSED CLUSTERS:\n";

    // get used clusters IDs
    std::unordered_set<int> usedClusterID;
    for(SLinearCluster &cluster:recoTrackList){
        std::vector<int> clusterIDs = cluster.GetHitCluster().GetClusterIDs();
        for(int &id:clusterIDs){
            std::cout<<cluster.NHits()<<" "<<id<<std::endl;
            usedClusterID.insert(id);
        }
    }

    // loop over the out of ROI hits
    std::vector<SHit> hitsToAdd;
    std::unordered_set<int> freeClusterIDs;
    for(SHit &h:fHitListOutROI){
        // check the hit cluster ID is not in the unordered set
        if( usedClusterID.find(h.ClusterId()) == usedClusterID.end() ){
            hitsToAdd.push_back(h);
            freeClusterIDs.insert(h.ClusterId());
        }
    }


    std::map<int, std::vector<SHit>> newHitsMap;
    for (const auto& id:freeClusterIDs){
        newHitsMap[id] = {};
        std::cout<<" Unused cluster ID "<<id<<std::endl;
    }

    for(SHit &h:hitsToAdd){
        newHitsMap[h.ClusterId()].push_back(h);
    }

    std::vector<SLinearCluster> newRecoTrackList = recoTrackList;
    for(auto& clusterPair:newHitsMap){
        SLinearCluster newCluster(clusterPair.second);
        newRecoTrackList.push_back(newCluster);
    }


    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Set the input hits
bool TPCLinesAlgo::SetHitList(int view,
                            std::vector<int>& vertex,
                            std::vector<int>& vertexTrue,
                            std::vector<SHit> hits)
{   
    // reset the class variables
    fHitList.clear();
    fHitListOutROI.clear();
    fNTotalHits=0;
    fClusterIdCounter.clear();

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
            }
            else if(std::abs(x-vertexX)<1500){ // make sure is the same TPC
                fHitListOutROI.push_back(hits[i]);
            }
        }

        if ( (int)fHitList.size() >= fTPCLinesPset.MinTrackHits) {
            
            // Shift hits to have origin in (0,0), leave 3 ticks under/overflow
            double minX = std::min_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.X()<h2.X();})->X() - 3;
            double minY = std::min_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.Y()<h2.Y();})->Y() - 3;
            double maxX = std::max_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.X()>h2.X();})->X() - 3;        
            // Fill cluster counter map
            int maxClusterId = std::max_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.ClusterId()<h2.ClusterId();} )->ClusterId();        

            // Set cluster id counter
            for (int cID = 0; cID <= maxClusterId; cID++) {
                fClusterIdCounter[cID] = 0;
            }

            for (size_t i = 0; i < fHitList.size(); i++) {
                SHit h = fHitList[i];
                fHitList[i] = SHit(i, h.X() - minX, h.Y() - minY,
                                    h.Width(), h.Integral(),
                                    h.StartT() - minY, h.EndT() - minY,
                                    h.Chi2(), h.ClusterId()
                                    );
                fClusterIdCounter[h.ClusterId()]++;

            }

            for (size_t i = 0; i < fHitListOutROI.size(); i++) {
                SHit h = fHitListOutROI[i];
                fHitListOutROI[i] = SHit(fHitList.size()+i, h.X() - minX, h.Y() - minY,
                                        h.Width(), h.Integral(),
                                        h.StartT() - minY, h.EndT() - minY,
                                        h.Chi2(), h.ClusterId()
                                    );
            }


            fNTotalHits = fHitList.size();
            // cout cluster ID counter
            std::cout<<" CLM  "<<maxClusterId<<std::endl;
            for(auto & co:fClusterIdCounter){
                std::cout<<co.first<<":"<<co.second<<"  ";
            }


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
// Get the hits in a given cluster
std::vector<SHit> TPCLinesAlgo::GetHitsInCluster(int clusterId){

    std::vector<SHit> hits;
    for(SHit &h:fHitList){
        if(h.ClusterId()==clusterId) hits.push_back(h);
    }

    return hits;
}


//----------------------------------------------------------------------
// Main function
void TPCLinesAlgo::AnaView(std::string eventLabel)
{
    // --- Reconstructed objects
    std::vector<SHit> discardedHits;
    std::vector<SLinearCluster> finalLinearClusterV;

    for(std::pair<int, int> clusterPair: fClusterIdCounter){
            
        //std::vector<SHit> hitListForHough = fHitList;
        std::vector<SHit> hitListForHough = GetHitsInCluster(clusterPair.first);

        std::cout<<" Ana cluster "<<clusterPair.first<<" "<<clusterPair.second<<std::endl;
    
        int trkIterCluster = (int)hitListForHough.size()>=fTPCLinesPset.MinTrackHits ? 0 : fTPCLinesPset.MaxHoughTracks;
        std::cout<<trkIterCluster<<" "<<fHitList.size()<<" "<<hitListForHough.size()<<std::endl;
        // --- Loop over the hough tracks
        while(trkIterCluster<fTPCLinesPset.MaxHoughTracks){

            if(fTPCLinesPset.Verbose>=2) std::cout<<" **** Track finder iteration: "<<trkIterCluster<< " # hough hits:"<<hitListForHough.size()<<std::endl;

            // -- Get the best Hough line       
            HoughLine houghLine = fHoughAlgo.GetBestHoughLine(hitListForHough, fVertex);
            if(fTPCLinesPset.Verbose>=3)
                if(fTPCLinesPset.Verbose>=2) std::cout<<"    Hough line results: "<<houghLine.NHits()<<" Score: "<<houghLine.Score()<<std::endl;

            // -- Skip if not enough hits
            if(houghLine.NHits() < fTPCLinesPset.MinTrackHits){
                if(fTPCLinesPset.Verbose>=3)
                    if(fTPCLinesPset.Verbose>=2) std::cout<<"   Hough lines has <"<<fTPCLinesPset.MinTrackHits<<", ending the search\n";
                trkIterCluster = fTPCLinesPset.MaxHoughTracks;
            }

            // -- Call the track finfder otherwise
            else{
                std::vector<SLinearCluster> linearClusterV = fTrackFinder.ReconstructTracksFromHoughDirection(hitListForHough, houghLine.GetLineEquation(), trkIterCluster);

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
                        trkIterCluster = fTPCLinesPset.MaxHoughTracks;
                    }
                }
                // end Hough loop otherwise
                else{
                    if(fTPCLinesPset.Verbose>=3)
                        if(fTPCLinesPset.Verbose>=2) std::cout<<"   Not found new tracks\n";
                    trkIterCluster = fTPCLinesPset.MaxHoughTracks;
                }
                
                // -- Remove isolated hits
                if (fTPCLinesPset.RemoveIsolatedHits) {
                    hitListForHough = RemoveIsolatedHits(hitListForHough, discardedHits, fTPCLinesPset.MaxNeighbourDistance, fTPCLinesPset.MinNeighboursHits);
                }
                
            }

            trkIterCluster++;
        } // -- end Hough loop

        discardedHits.insert(discardedHits.end(), hitListForHough.begin(), hitListForHough.end());

    }

    for(size_t ix = 0; ix<finalLinearClusterV.size(); ix++){
        std::cout<<" ix:"<<ix<<" Min/MaxX:"<<finalLinearClusterV[ix].GetMinX()<<" "<<finalLinearClusterV[ix].GetMaxX()<<std::endl;
    }
    
    
    // Slope track merger
    // sort by minX
    std::sort(finalLinearClusterV.begin(), finalLinearClusterV.end(), [&](SLinearCluster& l1, SLinearCluster& l2) {return l1.GetMinX() < l2.GetMinX();} );    
    finalLinearClusterV = TPCLinesDirectionUtils::SlopeTrackMerger(finalLinearClusterV, 2, 15, fTPCLinesPset.Verbose); 

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
    
        if (shortToLongTrackDict.find(track.GetId()) == shortToLongTrackDict.end()) {

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



        // Isolated hit merger
        std::vector<SHit> remainingHits;// = hitListForHough;
        remainingHits.insert(remainingHits.end(), discardedHits.begin(), discardedHits.end());
        discardedHits.clear();
        NewTrackList = MergeIsolatedHits(NewTrackList, remainingHits, 10, discardedHits);
        // Characterize the tracks
        for(size_t ix = 0; ix<NewTrackList.size(); ix++){
            NewTrackList[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            NewTrackList[ix].AssignId(ix);
        }

        // Out of ROI hits merger
        /*NewTrackList = MergeOutOfROIHits(NewTrackList, {});
        // Characterize the tracks
        for(size_t ix = 0; ix<NewTrackList.size(); ix++){
            NewTrackList[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            NewTrackList[ix].AssignId(ix);
        }*/


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

    NewTrackList = CreateOutOfROIClusters(NewTrackList);
    // Characterize the tracks
    for(size_t ix = 0; ix<NewTrackList.size(); ix++){
        NewTrackList[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
        NewTrackList[ix].AssignId(ix);
    }

    double hitDensity = GetAverageHitDensity();
    std::cout<<"Hit density: "<<hitDensity<<std::endl;

    discardedHits.insert(discardedHits.end(), fHitListOutROI.begin(), fHitListOutROI.end());
    fUnmatchedHits = discardedHits;
    SEvent recoEvent(NewTrackList, intersectionsInBall, vertexList, associatedOrigins, hitDensity, discardedHits);
    fRecoEvent = recoEvent;

    return;
}



void TPCLinesAlgo::FillLambdaAnaTree(LambdaAnaTree &lambdaAnaTree){



    // --- Hits for FRANS ---
    std::vector<SHit> allHits = fHitList;
    if(fHitListOutROI.size()>0)
        allHits.insert(allHits.end(), fHitListOutROI.begin(), fHitListOutROI.end());
    
    // --- FRANS with PANDORA vertex ---
    std::cout<<" FRANS PANDORA with "<<allHits.size()<<" hits\n";
    fFRANSAlgo.Fill(allHits, fVertex);
    double FRANSScorePANDORA = fFRANSAlgo.Score();
    lambdaAnaTree.fFRANSScorePANDORA = FRANSScorePANDORA;
    
    // --- FRANS with custom origin: get the best triangle ---
    std::vector<STriangle> angleList = fRecoEvent.GetAngleList();
    std::vector<SOrigin> associatedOrigins = fRecoEvent.GetAssociatedOrigins();
    int bestTriangleIx = -1;
    double bestFRANSScore = -1000;
    fHasTriangle = false;
    for(size_t orix=0; orix<angleList.size(); orix++){
        SVertex fVertexMine(associatedOrigins[orix].GetPoint(), "");
        
        std::cout<<" FRANS with "<<allHits.size()<<" hits\n";
        fFRANSAlgo.Fill(allHits, fVertexMine);
        double score = fFRANSAlgo.Score();

        if(score>bestFRANSScore){
            bestFRANSScore = score;
            bestTriangleIx = orix;
        }
    }

    // --- Number of origins variables ---
    lambdaAnaTree.fNOrigins = fRecoEvent.GetNOrigins();
    lambdaAnaTree.fNOriginsMult1 = fRecoEvent.GetNOriginsMult(1);
    lambdaAnaTree.fNOriginsMult2 = fRecoEvent.GetNOriginsMult(2);
    lambdaAnaTree.fNOriginsMultGT3 = fRecoEvent.GetNOriginsMultGt(3);
    lambdaAnaTree.fNOriginsPairOneTwo = fRecoEvent.GetNOriginsMult(2) * fRecoEvent.GetNOriginsMult(1);
    lambdaAnaTree.fNAngles = fRecoEvent.GetNAngles(); 

    // --- Best triangle variables ---
    if(bestTriangleIx!=-1){
        fHasTriangle = true;
        STriangle bestTriangle = angleList[bestTriangleIx];

        // --- Triangle variables ---
        lambdaAnaTree.fAngleFRANSScore = bestFRANSScore;
        lambdaAnaTree.fAngleGap = bestTriangle.GetGap();
        lambdaAnaTree.fAngleDecayContainedDiff = bestTriangle.GetDecayAngleDifference();
        lambdaAnaTree.fAngleNHits = bestTriangle.GetNHitsTriangle();
        lambdaAnaTree.fAngleNHitsTrack1 = bestTriangle.GetNHitsTrack1();
        lambdaAnaTree.fAngleNHitsTrack2 = bestTriangle.GetNHitsTrack2();
        lambdaAnaTree.fAngleNHitsMainTrack = bestTriangle.GetNHitsMainTrack();
        lambdaAnaTree.fAngleLengthTrack1 = bestTriangle.GetLengthTrack1();
        lambdaAnaTree.fAngleLengthTrack2 = bestTriangle.GetLengthTrack2();
        lambdaAnaTree.fAngleLengthMainTrack = bestTriangle.GetLengthMainTrack();
        lambdaAnaTree.fAngleMinNHits = std::min(bestTriangle.GetNHitsTrack1(), bestTriangle.GetNHitsTrack2());
        lambdaAnaTree.fAngleOpeningAngle = bestTriangle.GetOpeningAngle();
        
        // --- Unassociated origins ---
        SEvent notAssociatedRecoEvent = fRecoEvent.UnassociatedOrigins(bestTriangle);
        lambdaAnaTree.fNUnOrigins = notAssociatedRecoEvent.GetNOrigins();
        lambdaAnaTree.fNUnOriginsMult1 = notAssociatedRecoEvent.GetNOriginsMult(1);
        lambdaAnaTree.fNUnOriginsMult2 = notAssociatedRecoEvent.GetNOriginsMult(2);
        lambdaAnaTree.fNUnOriginsMultGT3 = notAssociatedRecoEvent.GetNOriginsMultGt(3);
        if( bestTriangle.GetLengthMainTrack()> bestTriangle.GetLengthTrack1() && bestTriangle.GetLengthMainTrack()> bestTriangle.GetLengthTrack2() ){
            lambdaAnaTree.fAngleLongestIsMain = true;
        }
        else{
            lambdaAnaTree.fAngleLongestIsMain = false;
        }
        lambdaAnaTree.fAngleCoveredArea = bestTriangle.ComputeCoveredArea();

        // --- Triangle cleaness ---
        int nDirtHitsInTriangle = 0;
        double nFractionDirtHitsInTriangle = 0;
        int nDirtHitsInTriangleWires = 0;
        double nFractionDirtHitsInTriangleWires = 0;
        fRecoEvent.FreeHitsAroundTriangle(bestTriangle,
                                        nDirtHitsInTriangle,
                                        nFractionDirtHitsInTriangle,
                                        nDirtHitsInTriangleWires,
                                        nFractionDirtHitsInTriangleWires);
          
        lambdaAnaTree.fAngleDirtHits = nDirtHitsInTriangle;
        lambdaAnaTree.fAngleDirtHitsWires = nDirtHitsInTriangleWires;
        lambdaAnaTree.fAngleDirtHitsRatio = nFractionDirtHitsInTriangle;
        lambdaAnaTree.fAngleDirtHitsWiresRatio = nFractionDirtHitsInTriangleWires;

        // --- Unassociated hits ---
        int nFreeHits = 0;
        int nUnassociatedHits = 0;
        fRecoEvent.GetUnassociatedHits(bestTriangle, nFreeHits, nUnassociatedHits);               
        lambdaAnaTree.fNUnassociatedHits = nUnassociatedHits;

                
        // --- Triangle calorimetry ---
        fTriangleCalo.SetTriangle(bestTriangle);
        fTriangleCalo.JointFitAnalysis();

        lambdaAnaTree.fAnglePassFit = fTriangleCalo.PassFit();
        lambdaAnaTree.fAngleTwoLinesChi2 = fTriangleCalo.TwoLinesChi2();
        lambdaAnaTree.fAngleNVertexHits = fTriangleCalo.NVertexHits();
        lambdaAnaTree.fAngleNBulkHits = fTriangleCalo.NBulkHits();
        lambdaAnaTree.fAngleVertexHitIntegralRatio = fTriangleCalo.VertexHitIntegralRatio();
        lambdaAnaTree.fAngleVertexHitIntegralDifference = fTriangleCalo.VertexHitIntegralDifference();
        lambdaAnaTree.fAngleVertexHitIntegralRelativeDifference = fTriangleCalo.VertexHitIntegralRelativeDifference();
        lambdaAnaTree.fAngleTrackLength1 = fTriangleCalo.TrackLength1();
        lambdaAnaTree.fAngleTrackLength2 = fTriangleCalo.TrackLength2();
        lambdaAnaTree.fAngleTrackLengthRatio = fTriangleCalo.TrackLengthRatio();
        lambdaAnaTree.fAngleResidualRange1RMS = fTriangleCalo.ResidualRange1RMS();
        lambdaAnaTree.fAngleResidualRange2RMS = fTriangleCalo.ResidualRange2RMS();
        lambdaAnaTree.fAngleResidualRangeMinRMS = fTriangleCalo.ResidualRangeMinRMS();
        lambdaAnaTree.fAngleChargeRatioAverage = fTriangleCalo.ChargeRatioAverage();
        lambdaAnaTree.fAngleChargeDifferenceAverage = fTriangleCalo.ChargeDifferenceAverage();
        lambdaAnaTree.fAngleChargeRelativeDifferenceAverage = fTriangleCalo.ChargeRelativeDifferenceAverage();
        lambdaAnaTree.fAnglePassChargeFit = fTriangleCalo.PassChargeFit();
        lambdaAnaTree.fAngleBandOverlap = fTriangleCalo.BandOverlap();
        lambdaAnaTree.fAngleBandCrossHits = fTriangleCalo.BandCrossHits();
        lambdaAnaTree.fAngleChargeRatioFit = fTriangleCalo.ChargeRatioFit();
        lambdaAnaTree.fAngleChargeDifferenceFit = fTriangleCalo.ChargeDifferenceFit();
        lambdaAnaTree.fAngleChargeRelativeDifferenceFit = fTriangleCalo.ChargeRelativeDifferenceFit();
        lambdaAnaTree.fAngleChargeRatioIntegral = fTriangleCalo.ChargeRatioIntegral();
        lambdaAnaTree.fAngleChargeDifferenceIntegral = fTriangleCalo.ChargeDifferenceIntegral();
        lambdaAnaTree.fAngleChargeRelativeDifferenceIntegral = fTriangleCalo.ChargeRelativeDifferenceIntegral();   
    }

}