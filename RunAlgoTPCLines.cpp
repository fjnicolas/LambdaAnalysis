// Objects
#include "src/SObjects/TPCLinesParameters.cpp"
#include "src/SObjects/TPCSimpleHits.cpp"
#include "src/SObjects/TPCSimpleLines.cpp"
#include "src/SObjects/TPCSimpleClusters.cpp"
#include "src/SObjects/TPCSimpleTriangles.cpp"

// Algorithms
#include "src/TPCLinesHough.cpp"
#include "src/TPCLinesTrackFinder.cpp"
#include "src/TPCLinesVertexFinder.cpp"
#include "src/TPCLinesAlgo.cpp"
/*#include "src/SObjects/TPCLinesDistanceUtils.cpp"
#include "src/SObjects/TPCLinesPCA.cpp"
#include "src/TPCLinesDBSCAN.cpp"
#include "src/TPCLinesDirectionRecoUtils.cpp"*/


// Display
#include "src/TPCLinesDisplay.cpp"


void RunAlgoTPCLines(int Debug=0, int DebugMode=-1, int n=1e6, int nskip=-1, int event=-1, int sr=-1, std::string file_name="", const char *directory_path=".", const char *ext=".root")
{

    //file_name = "analyzeItOutput_CCQE_R1-1_SR1-103.root";
    //file_name = "analyzeItOutput_V0Lambda_R1-1_SR1-195.root";
    std::cout<<"DEBUG "<<Debug<<std::endl;
    // program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;

    std::string fView="V";
    std::string fAppDisplayPath = "plots/";
    if(DebugMode==0) fAppDisplayPath = "plotsbg";
    else if(DebugMode==1) fAppDisplayPath = "plotssignal";

    std::vector<double> fReadoutWindow = {-200, 1500};
    int fStampTime = -200;
    double fSamplingFrequency = 500;
    double fSamplingTime = 0.5;
    int fReadoutWindowSize = 3400;


    // Algorithm paramters
    // Track finder parameters
    double fMaxDTube = 10;
    double fMaxDCluster = fMaxDTube/2;
    bool fSingleWireMode = false;
    int fMinClusterHits = 3;
    double fDCleaning=2.;
    double fClusterCompletenessCut=0.8;
    double fClusterAngleCut = 5;
    bool fCaptureMissingHits = true;
    int fMinTrackHits = 3;
    int fVerboseTrack = Debug;
    TrackFinderAlgorithmPsetType fPsetTrackFinder(fMaxDTube, fMaxDCluster, fSingleWireMode, fMinClusterHits,
                                    fDCleaning, fClusterCompletenessCut, fClusterAngleCut,
                                    fCaptureMissingHits, fMinTrackHits, fVerboseTrack);
    // Hough algorithm parameters
    double fMaxRadiusLineHypothesis = 25; //in cm
    double fThetaRes=25; //degrees
    double fMaxDistanceTube = 10;
    int fMinHoughHits = 3;
    int fVerboseHough = Debug;
    HoughAlgorithmPsetType fPsetHough(fMaxRadiusLineHypothesis,
                                fThetaRes,
                                fMaxDistanceTube,
                                fMinHoughHits,
                                fVerboseHough);
    //Ana view parameters
    double fMaxRadius = 250;
    double fDriftConversion = 1.;
    int fMaxHoughTracks = 15;
    bool fRemoveIsolatedHits = true;
    double fMaxNeighbourDistance = 2;
    double fMinNeighboursHits = 2;
    int fVerbose = Debug;
    int fDebugMode = DebugMode;
    TPCLinesAlgoPsetType fPsetAnaView(fMaxRadius, fDriftConversion, fMaxHoughTracks, fMinTrackHits,
    fRemoveIsolatedHits, fMaxNeighbourDistance, fMinNeighboursHits,
    fVerbose, fDebugMode, fPsetHough, fPsetTrackFinder);
        
    // Get the candidate Files
    std::vector<TString> fFilePaths;
    TSystemDirectory dir(".", ".");
    TList *files = dir.GetListOfFiles();
    TString targetFileName(file_name);
    std::cout<<" Target file name "<<targetFileName<<std::endl;
    gSystem->Exec("ls");
    if (files){
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            
            if (!file->IsDirectory() && fname.EndsWith(ext)){
                if(fname.Contains(targetFileName)){
                    fFilePaths.push_back(fname);
                }
                std::cout << fname << " Target:" << targetFileName <<std::endl;
                std::cout<<" Contains: "<<fname.Contains(targetFileName)<<std::endl;
                
                
            }
        }
    }


    // TPC LINES ALGORITHM
    TPCLinesAlgo _TPCLinesAlgo(fPsetAnaView, fAppDisplayPath);

    // TTree loop
    int nEvents = 0;
    int nProcessedEvents = 0;
    int nEventsSkipped = 0;
    int nEventsSelected = 0;
    std::vector<int> fNVertex;
    int nEntries = 0;
    for (const auto& filepath : fFilePaths) {

        std::cout<<" ANALYZING THE FILE: "<<filepath<<std::endl;

        TFile* f= new TFile(filepath);
	    TTree* tree = (TTree*)f->Get("ana/AnaTPCTree");
        
        // event ID
        int eventID;
        int subrunID;
        int runID;
        tree->SetBranchAddress("EventID",&eventID);
        tree->SetBranchAddress("RunID",&runID);
        tree->SetBranchAddress("SubRunID",&subrunID);
        // true Vertex
        double nuvE, nuvT, nuvX, nuvY, nuvZ;
        int nuvU, nuvV, nuvC, nuvTimeTick;
        int TPC;
        double nuvDriftTime;
        tree->SetBranchAddress("TrueVEnergy",&nuvE);
        tree->SetBranchAddress("TrueVt",&nuvT);
        tree->SetBranchAddress("TrueVx",&nuvX);
        tree->SetBranchAddress("TrueVy",&nuvY);
        tree->SetBranchAddress("TrueVz",&nuvZ);
        tree->SetBranchAddress("TrueVU",&nuvU);
        tree->SetBranchAddress("TrueVV",&nuvV);
        tree->SetBranchAddress("TrueVC",&nuvC);
        tree->SetBranchAddress("TrueVTimeTick",&nuvTimeTick);
        // reco Vertex
        double recnuvX = -1, recnuvY = -1, recnuvZ = -1;
        int recnuvU = -1, recnuvV = -1, recnuvC = -1, recnuvTimeTick = -1;
        tree->SetBranchAddress("RecoVx",&recnuvX);
        tree->SetBranchAddress("RecoVy",&recnuvY);
        tree->SetBranchAddress("RecoVz",&recnuvZ);
        tree->SetBranchAddress("RecoVU",&recnuvU);
        tree->SetBranchAddress("RecoVV",&recnuvV);
        tree->SetBranchAddress("RecoVC",&recnuvC);
        tree->SetBranchAddress("RecoVTimeTick",&recnuvTimeTick);

        std::vector<double> * hitsIntegral=new std::vector<double>;
        std::vector<double> * hitsPeakTime=new std::vector<double>;
        std::vector<int> * hitsChannel=new std::vector<int>;
        std::vector<double> * hitsRMS=new std::vector<double>;
        //std::vector<double> * hitsWidth=new std::vector<double>;
        std::vector<double> * hitsStartT=new std::vector<double>;
        std::vector<double> * hitsEndT=new std::vector<double>;
        std::vector<double> * hitsChi2=new std::vector<double>;
        std::vector<double> * hitsNDF=new std::vector<double>;
        std::vector<int> * hitsClusterID=new std::vector<int>;        
        tree->SetBranchAddress("HitsIntegral",&hitsIntegral);
        tree->SetBranchAddress("HitsPeakTime",&hitsPeakTime);
        tree->SetBranchAddress("HitsChannel",&hitsChannel);
        tree->SetBranchAddress("HitsRMS",&hitsRMS);
        //tree->SetBranchAddress("HitsWidth",&hitsWidth);
        tree->SetBranchAddress("HitsStartT",&hitsStartT);
        tree->SetBranchAddress("HitsEndT",&hitsEndT);
        tree->SetBranchAddress("HitsChi2",&hitsChi2);
        tree->SetBranchAddress("HitsNDF",&hitsNDF);
        tree->SetBranchAddress("HitsClusterID",&hitsClusterID);


        for (int entry = 0; entry < tree->GetEntries(); entry++) {
            if (fNEv > 0 && nEvents >= fNEv) continue;
            tree->GetEntry(entry);
            if (eventID != fEv && fEv != -1) continue;
            if (subrunID != fSubRun && fSubRun != -1) continue;
            nEntries++;
            std::cout << "Analyzing: " << runID << " " << subrunID << " " << eventID << std::endl;
            if (nEntries <= fNEvSkip) continue;
            std::string eventLabel = "R="+std::to_string(runID)+" SR="+std::to_string(subrunID)+" E="+std::to_string(eventID);
            nEvents++;

            // True vertex
            TPC = (nuvX > 0) ? 1 : 0;
            nuvDriftTime = nuvTimeTick * fSamplingTime + fStampTime;
            std::vector<double> VertexXYZ = {nuvX, nuvY, nuvZ};
            std::vector<int> VertexUVYT = {nuvU, nuvV, nuvC, nuvTimeTick};
            std::cout << "   - True vertex (X, Y, Z, T) " << VertexXYZ[0] << " " << VertexXYZ[1] << " " << VertexXYZ[2] << " " << nuvT << " (U, V, C, TT): " << VertexUVYT[0] << " " << VertexUVYT[1] << " " << VertexUVYT[2] << " " << VertexUVYT[3] << std::endl;

            // Reco vertex
            double recnuvDriftTime = recnuvTimeTick * fSamplingTime + fStampTime;
            std::vector<double> RecoVertexXYZ = {recnuvX, recnuvY, recnuvZ};
            std::vector<int> RecoVertexUVYT = {recnuvU, recnuvV, recnuvC, recnuvTimeTick};
            std::cout << "  - Reco vertex (X, Y, Z) " << RecoVertexXYZ[0] << " " << RecoVertexXYZ[1] << " " << RecoVertexXYZ[2] << " (U, V, C, TT): " << RecoVertexUVYT[0] << " " << RecoVertexUVYT[1] << " " << RecoVertexUVYT[2] << " "<< RecoVertexUVYT[3]  << std::endl;

            size_t nhits = hitsChannel->size();
            std::cout << "  - NHits: " << nhits << std::endl;
            if(nhits<=3){
                std::cout<<"   SKIPPED NHits\n";
                nEventsSkipped++;
                continue;
            }
            if(recnuvU==-1){
               std::cout<<"    SKIPPED RECOVERTEX\n"; 
               continue;
               nEventsSkipped++;
            }

            std::string view = fView+std::to_string(TPC);

            // Set the hits
            // note! we use hitRMS for the width
            _TPCLinesAlgo.SetHitList(view, RecoVertexUVYT, 
                                    hitsChannel,
                                    hitsPeakTime,
                                    hitsIntegral, 
                                    hitsRMS,
                                    hitsStartT, 
                                    hitsEndT,
                                    "");

            // Analyze
            int nOrigins = _TPCLinesAlgo.AnaView(eventLabel);
            nProcessedEvents++;
            if(nOrigins>0)
                nEventsSelected+=1;

            std::cout << " ********** Final status-" << std::endl;
            std::cout << " ... NTotalEvents=" << nEvents << " NSkipped=" << nEventsSkipped << std::endl;
            std::cout << " ... NProcessed=" << nProcessedEvents << " NSelected=" << nEventsSelected;
            std::cout << " efficiency=" << static_cast<double>(nEventsSelected) / nProcessedEvents;
            std::cout << " efficiency all=" << static_cast<double>(nEventsSelected) / nEvents << std::endl;
        }
    }


    #include <iostream>
    
    std::cout << " ********** Final status-" << std::endl;
    std::cout << " ... NTotalEvents=" << nEvents << " NSkipped=" << nEventsSkipped << std::endl;
    std::cout << " ... NProcessed=" << nProcessedEvents << " NSelected=" << nEventsSelected;
    std::cout << " efficiency=" << static_cast<double>(nEventsSelected) / nProcessedEvents;
    std::cout << " efficiency all=" << static_cast<double>(nEventsSelected) / nEvents << std::endl;


    return 0;

}


