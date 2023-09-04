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
#include "TPCLinesDirectionRecoUtils.cpp"

//----------------------------------------------------------------------
// Constructor
TPCLinesAlgo::TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, std::string displayPath):
    fTPCLinesPset(tpcLinesAlgoPset),
    fNTotalHits(0),
    fHitList({}),
    fVertex(SPoint(-1, -1), ""),
    fHoughAlgo(tpcLinesAlgoPset.HoughAlgorithmPset),
    fTrackFinder(tpcLinesAlgoPset.TrackFinderAlgorithmPset),
    fVertexFinder(tpcLinesAlgoPset.VertexFinderAlgorithmPset),
    fDisplay(TPCLinesDisplay(tpcLinesAlgoPset.Verbose>0, displayPath))
{}


//----------------------------------------------------------------------
// Display function
void TPCLinesAlgo::Display(){

    fDisplay.Show("Final reco", fHitList,LineEquation(0, 0),fHitList);

    return;
}


//----------------------------------------------------------------------
// Set the input hits
void TPCLinesAlgo::SetHitList(std::string view,
                            std::vector<int>& vertex,
                            std::vector<int> *_X,
                            std::vector<double> *_Y,
                            std::vector<double> *_Int,
                            std::vector<double> *_Wi,
                            std::vector<double> *_ST,
                            std::vector<double> *_ET,
                            std::string eventLabel)
{   
    // reset variables
    fHitList.clear();
    fNTotalHits=0;

    size_t nTotalHits = _X->size();

    // Define the vertex
    double vertexX = vertex[2];
    if (view == "U0" || view == "U1") vertexX = vertex[0];
    if (view == "V0" || view == "V1") vertexX = vertex[1];
    double vertexY = vertex[3] / fTPCLinesPset.DriftConversion;
    std::cout << "  ** Vertex XY: " << vertexX << " " << vertexY << std::endl;
    if (vertexX != -1){

        // Get hits in the selected plane
        std::vector<double> filteredX, filteredY, filteredInt, filteredST, filteredET, filteredWi;
        for (int i = 0; i < nTotalHits; i++) {
            int x = _X->at(i);
            
            // filter channels for the view
            if ( x > fChB[view][0] && x <= fChB[view][1]) {
                double y = _Y->at(i)/fTPCLinesPset.DriftConversion;
                
                // filter distance to vertex
                double d = std::sqrt(std::pow(x - vertexX, 2) + std::pow(y - vertexY, 2));
                if (d < fTPCLinesPset.MaxRadius) {
                    filteredX.push_back( x );
                    filteredY.push_back( y );
                    filteredInt.push_back( _Int->at(i) );
                    filteredWi.push_back( _Wi->at(i) / fTPCLinesPset.DriftConversion );
                    filteredST.push_back( _ST->at(i) / fTPCLinesPset.DriftConversion );
                    filteredET.push_back( _ET->at(i) / fTPCLinesPset.DriftConversion );
                }
            }
        }

        if (filteredX.size() > 3) {
            // Shift hits to have origin in (0,0)
            double minX = *std::min_element(filteredX.begin(), filteredX.end()) - 3;
            double minY = *std::min_element(filteredY.begin(), filteredY.end()) - 3;

            
            for (int i = 0; i < filteredX.size(); i++) {
                SHit hit(i, filteredX[i] - minX, filteredY[i] - minY, filteredWi[i], filteredInt[i], filteredST[i] - minY, filteredET[i] - minY);
                fHitList.push_back(hit);
            }
            fNTotalHits = fHitList.size();

            // create vertex with common origin
            vertexX = vertexX - minX;
            vertexY = vertexY - minY;

            fVertex = SVertex(SPoint(vertexX, vertexY), view);
            std::cout << "  ** Origin vertex: " << fVertex;
            std::cout << "  ** NInputHits:"<<fNTotalHits<<std::endl;

        }
        else {
            std::cout << "   SKIPPED NHits selected near the vertex " << filteredX.size() << std::endl;
        }

    }
}


//----------------------------------------------------------------------
// Remove isolated hits
std::vector<SHit> TPCLinesAlgo::RemoveIsolatedHits(std::vector<SHit> hitListForHough, std::vector<SHit>& discardedHits, double maxD, int minNeighbours) {

    std::cout << "\n\n +-+-+-+-+-+-+- Removing isolated hits +-+-+-+-+-+-+-" << std::endl;

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
    std::cout << "\n\n +-+-+-+-+-+-+- Isolated hit merger +-+-+-+-+-+-+-" << std::endl;
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
                double hypoY = track.GetTrackEquation().Slope() * hit.X() + track.GetTrackEquation().Intercept();
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
        std::cout << "\n Merging track " << trkix << " COMP " << trackComp << " CONN " << trackConn << std::endl;

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

        std::cout << "  candidates hits: " << std::endl;
        for (size_t ix = 0; ix < candidateHits.size(); ++ix) {
            SHit hit = candidateHits[sortedCandidateHits[ix]];

            double minD = currentCluster.GetMinDistanceToCluster(hit);
            double minDConn = currentCluster.GetMinDistanceToClusterOverlap(hit);
            std::cout << hit.X() << " " << hit.Y() << " " << minD << " " << minDConn << std::endl;
            std::cout << "       " << hit.Id() << " " << hit.X() << " " << hit.Y() << " d " << minD << " dConn " << minDConn << std::endl;
            if (minDConn < dTh * trackConn) {
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
int TPCLinesAlgo::AnaView(std::string eventLabel)
{

    // --- Objects for the Hough tracks
    int trkIter = 0;
    std::vector<SHit> hitListForHough = fHitList;
    std::vector<SHit> discardedHits;
    std::vector<SLinearCluster> finalLinearClusterV;

    // --- Loop over the hough tracks
    while(trkIter<fTPCLinesPset.MaxHoughTracks){

        std::cout<<" **** Track finder iteration: "<<trkIter<< " # hough hits:"<<hitListForHough.size()<<std::endl;

        // -- Get the best Hough line       
        HoughLine houghLine = fHoughAlgo.GetBestHoughLine(hitListForHough, fVertex);
        if(fTPCLinesPset.Verbose>=3)
            std::cout<<"    Hough line results: "<<houghLine.NHits()<<" Score: "<<houghLine.Score()<<std::endl;

        // -- Skip if not enough hits
        if(houghLine.NHits() < fTPCLinesPset.MinTrackHits){
            if(fTPCLinesPset.Verbose>=3)
                std::cout<<"   Hough lines has <"<<fTPCLinesPset.MinTrackHits<<", ending the search\n";
            trkIter = fTPCLinesPset.MaxHoughTracks;
        }
        // -- Call the track finfder otherwise
        else{
            std::vector<SLinearCluster> linearClusterV = fTrackFinder.ReconstructTracksFromHoughDirection(hitListForHough, houghLine.GetLineEquation(), trkIter);

            // -- check the found tracks has enough hits
            if(linearClusterV.size()!=0){
                if(fTPCLinesPset.Verbose>=2)
                    std::cout<<"   Found "<<linearClusterV.size()<<" new tracks!\n";
                finalLinearClusterV.insert(finalLinearClusterV.begin(), linearClusterV.begin(), linearClusterV.end());
                
                // -- check there's enough hits for the next Hough iteration
                if(hitListForHough.size() < fTPCLinesPset.MinTrackHits){
                    if(fTPCLinesPset.Verbose>=3)
                        std::cout<<"   Remaining hits are "<<hitListForHough.size()<<", ending the search\n";
                    trkIter = fTPCLinesPset.MaxHoughTracks;
                }
            }
            // end Hough loop otherwise
            else{
                if(fTPCLinesPset.Verbose>=3)
                    std::cout<<"   Not found new tracks\n";
                trkIter = fTPCLinesPset.MaxHoughTracks;
            }
            
            // -- Remove isolated hits
            if (fTPCLinesPset.RemoveIsolatedHits) {
                hitListForHough = RemoveIsolatedHits(hitListForHough, discardedHits, fTPCLinesPset.MaxNeighbourDistance, fTPCLinesPset.MinNeighboursHits);
            }
            
        }

        trkIter++;
    } // -- end Hough loop
    
    // Slope track merger
    // sort by minX
    std::sort(finalLinearClusterV.begin(), finalLinearClusterV.end(), [&](SLinearCluster& l1, SLinearCluster& l2) {return l1.GetMinX() < l2.GetMinX();} );    
    finalLinearClusterV = TPCLinesDirectionUtils::SlopeTrackMerger(finalLinearClusterV, 2, 15); 


    // Isolated hit merger
    std::vector<SHit> remainingHits = hitListForHough; 
    remainingHits.insert(remainingHits.end(), discardedHits.begin(), discardedHits.end());   
    finalLinearClusterV = MergeIsolatedHits(finalLinearClusterV, remainingHits, 10);

   
    std::cout<<"\n\n +-+-+-+-+-+-+- Characterizing tracks +-+-+-+-+-+-+-\n";
    for(size_t ix = 0; ix<finalLinearClusterV.size(); ix++){
        //// !!!! WARNING
        finalLinearClusterV[ix].FillResidualHits();
        finalLinearClusterV[ix].AssignId(ix);
    }
        
    
    //Find secondary vertexes
    std::vector<STriangle> vertexList;
    std::vector<SPoint> originList;
    SLinearCluster mainDirection;
    
    fVertexFinder.GetOrigins(finalLinearClusterV, vertexList, originList, mainDirection);

    bool accepted = vertexList.size()>0;
    std::string outNamePreffix = accepted? "Accepted Final Reco":"Rejected Final Reco";
    
    fDisplay.Show(outNamePreffix+eventLabel, fHitList, LineEquation(0, 0), {}, finalLinearClusterV, mainDirection, vertexList, fVertex);
    

    return vertexList.size();
}